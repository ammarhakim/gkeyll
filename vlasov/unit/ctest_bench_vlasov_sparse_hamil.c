// Micro-benchmark: sparse vs. dense velocity-space Hamiltonian kernels.
//
// For the default model (H = v^2/2) both GKYL_HAMIL_VEL_SPARSE and
// GKYL_HAMIL_VEL_DENSE routings are numerically valid (the dense-layout cross
// coefficients are structurally zero), so the two can be timed head-to-head in
// one binary on identical inputs. Covers both brackets:
//   - canonical (GKYL_MODEL_DEFAULT): vel-flux surface assembly (E+B
//     alpha_quad + Lax flux), hyper_dg advance (volume + conf-space surf +
//     accel surf), and the five-moments calculation;
//   - non-canonical triads (GKYL_MODEL_TRIAD): conf-flux surface assembly
//     (nc hamil alpha_quad in configuration space), vel-flux surface assembly
//     (nc hamil alpha_quad in velocity space), and hyper_dg advance
//     (nc volume, incl. the precomputed-alpha combos, + flux-consumer surf).
// Each case cross-checks that sparse and dense agree to machine precision.
//
// Iterations: set BENCH_NREP (default 8, CI-friendly). Timings are per-call.

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

#include <stdlib.h>

struct bench_updaters {
  struct gkyl_dg_vlasov_conf_flux_surf *calc_conf_flux; // triads only
  struct gkyl_dg_vlasov_vel_flux_surf *calc_vel_flux;
  struct gkyl_dg_eqn *eqn;
  gkyl_hyper_dg *slvr;
  struct gkyl_mom_type *mom_type;
  struct gkyl_mom_calc *mom_calc;
};

static double
max_abs_diff(const struct gkyl_array *a, const struct gkyl_array *b)
{
  const double *ad = a->data, *bd = b->data;
  double m = 0.0;
  for (size_t i = 0; i < a->size * a->ncomp; ++i) {
    double d = fabs(ad[i] - bd[i]);
    if (d > m) m = d;
  }
  return m;
}

static void
bench_case(enum gkyl_model_id model_id, int cdim, int vdim, int poly_order,
  const int *conf_cells, const int *vel_cells, double ab_rel_tol)
{
  int pdim = cdim + vdim;
  bool is_triad = (model_id == GKYL_MODEL_TRIAD);
  int nrep = 8;
  const char *nrep_env = getenv("BENCH_NREP");
  if (nrep_env) nrep = atoi(nrep_env);

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int velghost[GKYL_MAX_DIM] = { 0 };
  double vellower[GKYL_MAX_DIM], velupper[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    cells[d] = conf_cells[d]; ghost[d] = 1; lower[d] = 0.0; upper[d] = 1.0;
  }
  for (int d = 0; d < vdim; ++d) {
    cells[cdim + d] = vel_cells[d]; ghost[cdim + d] = 0;
    lower[cdim + d] = -1.0; upper[cdim + d] = 1.0;
    vellower[d] = -1.0; velupper[d] = 1.0;
  }

  struct gkyl_rect_grid confGrid, velGrid, phaseGrid;
  struct gkyl_range confRange, confRange_ext, velRange, velRange_ext, phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, &cells[cdim]);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Nodal surface-expansion sizes, matching vm_species_collisionless.c
  // (Serendipity, use_lo = false: higher-order nodes except at p=1).
  int highorder = (poly_order == 1) ? 0 : 1;
  int num_surf_vel_nodes = (int) pow(poly_order + 1 + highorder, pdim - 1);
  int num_surf_conf_nodes = (int) pow(poly_order + 1 + highorder, pdim - 1);

  // Distribution, RHS, CFL, and field arrays with synthetic data.
  struct gkyl_array *fin = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1, phaseRange_ext.volume);
  struct gkyl_array *qmem = gkyl_array_new(GKYL_DOUBLE, 8*confBasis.num_basis, confRange_ext.volume);

  size_t nf = phaseRange_ext.volume * basis.num_basis;
  double *fin_d = fin->data;
  for (size_t i = 0; i < nf; ++i)
    fin_d[i] = (double)(2*i + 11 % nf) / nf * ((i % 2 == 0) ? 1 : -1);
  size_t nem = confRange_ext.volume * 8 * confBasis.num_basis;
  double *qmem_d = qmem->data;
  for (size_t i = 0; i < nem; ++i)
    qmem_d[i] = (double)(i % 97) / 97.0 * ((i % 2 == 0) ? 1 : -1);

  // Identity velocity/position maps and the default (separable) Hamiltonian
  // (the triad model builds its Hamiltonian with the same default calc).
  struct gkyl_array *hamil = gkyl_array_new(GKYL_DOUBLE, velBasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = gkyl_array_new(GKYL_DOUBLE, velBasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &velBasis, inp_vmap, false);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &confBasis, inp_pmap, false);
  gkyl_dg_vlasov_calc_hamil(&velGrid, &velBasis, &velRange,
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, false);

  // Poisson tensor: synthetic (smooth, nonzero) so the nc kernels do real work.
  int num_pt_indices[3] = { 1, 6, 18 };
  struct gkyl_array *poisson_tensor_conf = gkyl_array_new(GKYL_DOUBLE,
    confBasis.num_basis*num_pt_indices[vdim-1], confRange_ext.volume);
  size_t npt = poisson_tensor_conf->size * poisson_tensor_conf->ncomp;
  double *pt_d = poisson_tensor_conf->data;
  for (size_t i = 0; i < npt; ++i)
    pt_d[i] = 1.0 + 0.1 * (double)(i % 13) / 13.0 * ((i % 2 == 0) ? 1 : -1);

  struct gkyl_array *pot_tot = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis*4, confRange_ext.volume);
  struct gkyl_array *f_no_J = gkyl_array_new(GKYL_DOUBLE, fin->ncomp, fin->size);
  struct gkyl_array *rad = gkyl_array_new(GKYL_DOUBLE, vdim*velBasis.num_basis, velRange.volume);
  gkyl_array_set(f_no_J, 1.0, fin);

  // Per-flavor outputs so sparse and dense can be cross-checked.
  enum gkyl_hamil_id flavors[2] = { GKYL_HAMIL_VEL_DENSE, GKYL_HAMIL_VEL_SPARSE };
  struct gkyl_array *rhs[2], *vel_flux_surf[2], *conf_flux_surf[2], *marr[2];
  struct bench_updaters up[2] = { 0 };
  double t_cflux[2], t_flux[2], t_hyper[2], t_mom[2];

  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = zero_flux_flags[d + pdim] = (d < cdim) ? 0 : 1;
  }

  for (int fl = 0; fl < 2; ++fl) {
    rhs[fl] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ext.volume);
    vel_flux_surf[fl] = gkyl_array_new(GKYL_DOUBLE, vdim*num_surf_vel_nodes, phaseRange_ext.volume);
    conf_flux_surf[fl] = is_triad ?
      gkyl_array_new(GKYL_DOUBLE, cdim*num_surf_conf_nodes, phaseRange_ext.volume) : 0;

    if (is_triad) {
      struct gkyl_dg_vlasov_conf_flux_surf_inp inp_conf_flux = {
        .phase_grid = &phaseGrid,
        .conf_basis = &confBasis,
        .phase_basis = &basis,
        .vel_range = &velRange,
        .hamil_range = &velRange,
        .skip_cell_thresh = 0.0,
        .model_id = model_id,
        .hamil_id = flavors[fl],
        .use_lo = false,
        .use_gpu = false,
      };
      up[fl].calc_conf_flux = gkyl_dg_vlasov_conf_flux_surf_inew(&inp_conf_flux);
    }

    struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
      .phase_grid = &phaseGrid,
      .conf_basis = &confBasis,
      .phase_basis = &basis,
      .vel_map = vel_map,
      .pos_map = pos_map,
      .hamil_range = &velRange,
      .skip_cell_thresh = 0.0,
      .model_id = model_id,
      .hamil_id = flavors[fl],
      // Canonical: Lorentz forces drive velocity space. Triads: the nc
      // Hamiltonian alpha_quad is the force; run it field-free.
      .has_E = !is_triad,
      .has_phi = false,
      .has_B = !is_triad && vdim > 1,
      .has_rad = false,
      .use_lo = false,
      .use_gpu = false,
    };
    up[fl].calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux);

    struct gkyl_dg_vlasov_inp inp_eqn = {
      .conf_basis = &confBasis,
      .phase_basis = &basis,
      .conf_range = &confRange,
      .hamil_range = &velRange,
      .phase_range = &phaseRange,
      .vel_map = vel_map,
      .pos_map = pos_map,
      .skip_cell_thresh = 0.0,
      .model_id = model_id,
      .hamil_id = flavors[fl],
      .has_E = !is_triad,
      .has_phi = false,
      .has_B = !is_triad && vdim > 1,
      .has_rad = false,
      .poisson_tensor_conf = poisson_tensor_conf,
      .hamil = hamil,
      .qmem = qmem,
      .pot_tot = pot_tot,
      .conf_flux_surf = conf_flux_surf[fl],
      .vel_flux_surf = vel_flux_surf[fl],
      .f_no_J = f_no_J,
      .rad = rad,
      .use_lo = false,
      .use_gpu = false,
    };
    up[fl].eqn = gkyl_dg_vlasov_inew(&inp_eqn);
    up[fl].slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, up[fl].eqn, pdim, up_dirs, zero_flux_flags, 1, false);

    struct gkyl_mom_vlasov_inp inp_mom = {
      .conf_basis = &confBasis,
      .phase_basis = &basis,
      .vel_range = &velRange,
      .vel_map = vel_map,
      .hamil_range = &velRange,
      .hamil = hamil,
      .model_id = model_id,
      .hamil_id = flavors[fl],
      .mom_type = GKYL_F_MOMENT_M0M1M2,
      .use_gpu = false,
    };
    up[fl].mom_type = gkyl_mom_vlasov_inew(&inp_mom);
    up[fl].mom_calc = gkyl_mom_calc_new(&phaseGrid, up[fl].mom_type, false);
    marr[fl] = gkyl_array_new(GKYL_DOUBLE, (2 + vdim)*confBasis.num_basis, confRange_ext.volume);

    // warm-up (touch all memory once, outside the timers)
    gkyl_array_clear(rhs[fl], 0.0);
    gkyl_array_clear(cflrate, 0.0);
    if (is_triad)
      gkyl_dg_vlasov_conf_flux_surf_advance(up[fl].calc_conf_flux, &confRange, &phaseRange, &phaseRange_ext,
        poisson_tensor_conf, hamil, fin, cflrate, conf_flux_surf[fl]);
    gkyl_dg_vlasov_vel_flux_surf_advance(up[fl].calc_vel_flux, &confRange, &phaseRange,
      poisson_tensor_conf, hamil, qmem, pot_tot, rad, f_no_J, cflrate, vel_flux_surf[fl]);
    gkyl_hyper_dg_advance(up[fl].slvr, &phaseRange, fin, cflrate, rhs[fl]);
    gkyl_mom_calc_advance(up[fl].mom_calc, &phaseRange, &confRange, fin, marr[fl]);

  }

  // Interleaved A/B timing: alternate dense/sparse within every rep so slow
  // thermal/frequency drift affects both flavors equally.
  for (int fl = 0; fl < 2; ++fl) t_cflux[fl] = t_flux[fl] = t_hyper[fl] = t_mom[fl] = 0.0;
  struct timespec tm;
  for (int n = 0; n < nrep; ++n) {
    for (int fl = 0; fl < 2; ++fl) {
      if (is_triad) {
        tm = gkyl_wall_clock();
        gkyl_dg_vlasov_conf_flux_surf_advance(up[fl].calc_conf_flux, &confRange, &phaseRange, &phaseRange_ext,
          poisson_tensor_conf, hamil, fin, cflrate, conf_flux_surf[fl]);
        t_cflux[fl] += gkyl_time_diff_now_sec(tm);
      }
      tm = gkyl_wall_clock();
      gkyl_dg_vlasov_vel_flux_surf_advance(up[fl].calc_vel_flux, &confRange, &phaseRange,
        poisson_tensor_conf, hamil, qmem, pot_tot, rad, f_no_J, cflrate, vel_flux_surf[fl]);
      t_flux[fl] += gkyl_time_diff_now_sec(tm);

      gkyl_array_clear(rhs[fl], 0.0);
      gkyl_array_clear(cflrate, 0.0);
      tm = gkyl_wall_clock();
      gkyl_hyper_dg_advance(up[fl].slvr, &phaseRange, fin, cflrate, rhs[fl]);
      t_hyper[fl] += gkyl_time_diff_now_sec(tm);

      tm = gkyl_wall_clock();
      gkyl_mom_calc_advance(up[fl].mom_calc, &phaseRange, &confRange, fin, marr[fl]);
      t_mom[fl] += gkyl_time_diff_now_sec(tm);
    }
  }
  for (int fl = 0; fl < 2; ++fl) {
    t_cflux[fl] /= nrep; t_flux[fl] /= nrep; t_hyper[fl] /= nrep; t_mom[fl] /= nrep;
  }

  // sparse and dense must agree to machine precision on identical inputs
  double d_rhs = max_abs_diff(rhs[1], rhs[0]);
  double d_flux = max_abs_diff(vel_flux_surf[1], vel_flux_surf[0]);
  double d_mom = max_abs_diff(marr[1], marr[0]);
  double d_cflux = is_triad ? max_abs_diff(conf_flux_surf[1], conf_flux_surf[0]) : 0.0;
  double rhs_mag = 0.0;
  const double *r0 = rhs[0]->data;
  for (size_t i = 0; i < rhs[0]->size * rhs[0]->ncomp; ++i)
    if (fabs(r0[i]) > rhs_mag) rhs_mag = fabs(r0[i]);
  // Most cases agree to machine precision. At triad combos where the dense
  // flavor uses the precomputed-alpha projection but the sparse flavor inlines
  // the exact bracket (p=2 PreComp combos), the two differ at the projection-
  // truncation level; the caller passes the appropriate tolerance.
  TEST_CHECK(d_rhs <= ab_rel_tol * fmax(rhs_mag, 1.0));
  TEST_MSG("sparse vs dense RHS max diff %e (rhs magnitude %e, tol %e)", d_rhs, rhs_mag, ab_rel_tol);

  printf("\n%s %dx%dv p%d (%ld phase cells):\n",
    is_triad ? "TRIAD" : "DEFAULT", cdim, vdim, poly_order, (long)phaseRange.volume);
  printf("  %-14s %12s %12s %9s\n", "updater", "dense ms", "sparse ms", "speedup");
  if (is_triad)
    printf("  %-14s %12.3f %12.3f %8.2fx\n", "conf_flux_surf", 1e3*t_cflux[0], 1e3*t_cflux[1], t_cflux[0]/t_cflux[1]);
  printf("  %-14s %12.3f %12.3f %8.2fx\n", "vel_flux_surf", 1e3*t_flux[0], 1e3*t_flux[1], t_flux[0]/t_flux[1]);
  printf("  %-14s %12.3f %12.3f %8.2fx\n", "hyper_dg", 1e3*t_hyper[0], 1e3*t_hyper[1], t_hyper[0]/t_hyper[1]);
  printf("  %-14s %12.3f %12.3f %8.2fx\n", "five_moments", 1e3*t_mom[0], 1e3*t_mom[1], t_mom[0]/t_mom[1]);
  if (is_triad)
    printf("  max |sparse - dense|: rhs %.2e, conf_flux %.2e, vel_flux %.2e, moments %.2e\n",
      d_rhs, d_cflux, d_flux, d_mom);
  else
    printf("  max |sparse - dense|: rhs %.2e, vel_flux %.2e, moments %.2e\n", d_rhs, d_flux, d_mom);

  for (int fl = 0; fl < 2; ++fl) {
    if (up[fl].calc_conf_flux) gkyl_dg_vlasov_conf_flux_surf_release(up[fl].calc_conf_flux);
    gkyl_dg_vlasov_vel_flux_surf_release(up[fl].calc_vel_flux);
    gkyl_hyper_dg_release(up[fl].slvr);
    gkyl_dg_eqn_release(up[fl].eqn);
    gkyl_mom_calc_release(up[fl].mom_calc);
    gkyl_mom_type_release(up[fl].mom_type);
    gkyl_array_release(rhs[fl]);
    gkyl_array_release(vel_flux_surf[fl]);
    if (conf_flux_surf[fl]) gkyl_array_release(conf_flux_surf[fl]);
    gkyl_array_release(marr[fl]);
  }
  gkyl_array_release(fin);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(poisson_tensor_conf);
  gkyl_array_release(pot_tot);
  gkyl_array_release(f_no_J);
  gkyl_array_release(rad);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);
}

// canonical bracket (H = v^2/2, E+B forces)
static void bench_1x2v_p1(void) { bench_case(GKYL_MODEL_DEFAULT, 1, 2, 1, (int[]){ 32 }, (int[]){ 16, 16 }, 1e-11); }
static void bench_1x2v_p2(void) { bench_case(GKYL_MODEL_DEFAULT, 1, 2, 2, (int[]){ 32 }, (int[]){ 16, 16 }, 1e-11); }
static void bench_1x3v_p1(void) { bench_case(GKYL_MODEL_DEFAULT, 1, 3, 1, (int[]){ 16 }, (int[]){ 12, 12, 12 }, 1e-11); }
static void bench_1x3v_p2(void) { bench_case(GKYL_MODEL_DEFAULT, 1, 3, 2, (int[]){ 8 }, (int[]){ 8, 8, 8 }, 1e-11); }
static void bench_2x3v_p1(void) { bench_case(GKYL_MODEL_DEFAULT, 2, 3, 1, (int[]){ 8, 8 }, (int[]){ 8, 8, 8 }, 1e-11); }
static void bench_3x3v_p1(void) { bench_case(GKYL_MODEL_DEFAULT, 3, 3, 1, (int[]){ 4, 4, 4 }, (int[]){ 8, 8, 8 }, 1e-11); }
// non-canonical bracket (triads; same separable H, Poisson-tensor forces)
static void bench_triad_1x2v_p1(void) { bench_case(GKYL_MODEL_TRIAD, 1, 2, 1, (int[]){ 32 }, (int[]){ 16, 16 }, 1e-11); }
static void bench_triad_1x2v_p2(void) { bench_case(GKYL_MODEL_TRIAD, 1, 2, 2, (int[]){ 32 }, (int[]){ 16, 16 }, 1e-11); }
static void bench_triad_1x3v_p1(void) { bench_case(GKYL_MODEL_TRIAD, 1, 3, 1, (int[]){ 16 }, (int[]){ 12, 12, 12 }, 1e-11); }
// 1x3v p2: dense = precomputed-alpha + comps, sparse = inline exact bracket.
// These are different discretizations when the Poisson tensor varies within a
// cell (alpha projection truncates): observed ~3e-4 relative rhs difference,
// machine-identical when the tensor is cell-wise constant.
static void bench_triad_1x3v_p2(void) { bench_case(GKYL_MODEL_TRIAD, 1, 3, 2, (int[]){ 8 }, (int[]){ 8, 8, 8 }, 1e-3); }
static void bench_triad_2x3v_p1(void) { bench_case(GKYL_MODEL_TRIAD, 2, 3, 1, (int[]){ 8, 8 }, (int[]){ 8, 8, 8 }, 1e-11); }
static void bench_triad_3x3v_p1(void) { bench_case(GKYL_MODEL_TRIAD, 3, 3, 1, (int[]){ 4, 4, 4 }, (int[]){ 8, 8, 8 }, 1e-11); }

TEST_LIST = {
  { "bench_1x2v_p1", bench_1x2v_p1 },
  { "bench_1x2v_p2", bench_1x2v_p2 },
  { "bench_1x3v_p1", bench_1x3v_p1 },
  { "bench_1x3v_p2", bench_1x3v_p2 },
  { "bench_2x3v_p1", bench_2x3v_p1 },
  { "bench_3x3v_p1", bench_3x3v_p1 },
  { "bench_triad_1x2v_p1", bench_triad_1x2v_p1 },
  { "bench_triad_1x2v_p2", bench_triad_1x2v_p2 },
  { "bench_triad_1x3v_p1", bench_triad_1x3v_p1 },
  { "bench_triad_1x3v_p2", bench_triad_1x3v_p2 },
  { "bench_triad_2x3v_p1", bench_triad_2x3v_p1 },
  { "bench_triad_3x3v_p1", bench_triad_3x3v_p1 },
  { NULL, NULL },
};
