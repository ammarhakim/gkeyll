#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_velocity_map.h>

#include <float.h>
#include <math.h>
#include <stdbool.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct lte_moment_ctx {
  int cdim;
  int vdim;
  bool fluctuate;
  double lower[GKYL_MAX_DIM];
  double upper[GKYL_MAX_DIM];
  double n;
  double u[3];
  double temp;
};

static double
unit_coord(const struct lte_moment_ctx *m, const double *xn, int d)
{
  return (xn[d]-m->lower[d])/(m->upper[d]-m->lower[d]);
}

static double
fluct_x(const struct lte_moment_ctx *m, const double *xn)
{
  if (!m->fluctuate)
    return 0.0;

  double x0 = unit_coord(m, xn, 0);
  double fluct = sin(2.0*GKYL_PI*x0);
  if (m->cdim > 1) {
    double x1 = unit_coord(m, xn, 1);
    fluct += 0.5*cos(2.0*GKYL_PI*x1);
    fluct += 0.25*sin(2.0*GKYL_PI*(x0+x1));
  }
  if (m->cdim > 2) {
    double x2 = unit_coord(m, xn, 2);
    fluct += 0.25*cos(2.0*GKYL_PI*x2);
  }

  return fluct;
}

static void
eval_lte_n(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  fout[0] = m->n*(1.0 + 0.12*fluct_x(m, xn));
}

static void
eval_lte_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  double x0 = m->fluctuate ? unit_coord(m, xn, 0) : 0.0;
  double x1 = m->fluctuate && m->cdim > 1 ? unit_coord(m, xn, 1) : 0.0;

  for (int d=0; d<m->vdim; ++d) {
    double du = 0.0;
    if (m->fluctuate) {
      du = 0.06*cos(2.0*GKYL_PI*(d+1)*x0)/(d+1);
      if (m->cdim > 1)
        du += 0.035*sin(2.0*GKYL_PI*(x1 + 0.25*d));
    }
    fout[d] = m->u[d] + du;
  }
}

static void
eval_lte_temp(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  double x0 = m->fluctuate ? unit_coord(m, xn, 0) : 0.0;
  double fluct = 0.0;
  if (m->fluctuate) {
    fluct = cos(2.0*GKYL_PI*x0);
    if (m->cdim > 1) {
      double x1 = unit_coord(m, xn, 1);
      fluct += 0.5*sin(2.0*GKYL_PI*x1);
    }
  }
  fout[0] = m->temp*(1.0 + 0.08*fluct);
}

static struct lte_moment_ctx
make_moment_ctx(enum gkyl_model_id model_id, int cdim, int vdim, bool fluctuate,
  const double *lower, const double *upper)
{
  struct lte_moment_ctx m = { 0 };
  m.cdim = cdim;
  m.vdim = vdim;
  m.fluctuate = fluctuate;
  for (int d=0; d<cdim; ++d) {
    m.lower[d] = lower[d];
    m.upper[d] = upper[d];
  }

  if (model_id == GKYL_MODEL_SR) {
    // Same constant Maxwell-Juttner moments used in ctest_correct_mj_integrated.
    m.n = 1.0;
    m.temp = 1.0;
    if (vdim == 1) {
      m.u[0] = 1.0;
    }
    else if (vdim == 2) {
      m.u[0] = 1.0;
      m.u[1] = 0.5;
    }
    else {
      m.u[0] = 1.5;
      m.u[1] = 1.0;
      m.u[2] = 0.5;
    }
  }
  else {
    // Same non-relativistic Maxwellian density/drifts as ctest_correct_maxwellian,
    // with the 3v drift extrapolated from the 1v/2v cases.
    m.n = 1.25;
    m.temp = 1.0;
    m.u[0] = 0.5;
    m.u[1] = vdim > 1 ? 0.25 : 0.0;
    m.u[2] = vdim > 2 ? 0.125 : 0.0;
  }

  return m;
}

static void
component_label(int comp, int vdim, char label[16])
{
  if (comp == 0)
    snprintf(label, 16, "n");
  else if (comp == vdim+1)
    snprintf(label, 16, "T/m");
  else
    snprintf(label, 16, "u%d", comp-1);
}

static void
calc_modal_scaled_errors(const struct gkyl_array *moms, const struct gkyl_array *moms_target,
  const struct gkyl_range *conf_local, int num_comp, int num_conf_basis, double *err)
{
  for (int c=0; c<num_comp*num_conf_basis; ++c)
    err[c] = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_local);
  while (gkyl_range_iter_next(&iter)) {
    long midx = gkyl_range_idx(conf_local, iter.idx);
    const double *moms_d = gkyl_array_cfetch(moms, midx);
    const double *target_d = gkyl_array_cfetch(moms_target, midx);

    for (int d=0; d<num_comp; ++d) {
      double scale = fabs(target_d[d*num_conf_basis]);
      if (d == 0 || d == num_comp-1)
        scale = scale > DBL_MIN ? scale : 1.0;
      else
        scale = scale < 1.0 ? 1.0 : scale;

      for (int k=0; k<num_conf_basis; ++k) {
        double scaled_err = fabs(moms_d[d*num_conf_basis+k] - target_d[d*num_conf_basis+k])/scale;
        err[d*num_conf_basis+k] = fmax(err[d*num_conf_basis+k], scaled_err);
      }
    }
  }
}

static double
avg_convergence_rate(double initial_err, double final_err, int niter)
{
  if (niter <= 0 || initial_err <= DBL_MIN)
    return -1.0;
  if (final_err <= DBL_MIN)
    return 0.0;

  return pow(final_err/initial_err, 1.0/niter);
}

static void
format_rate(double initial_err, double final_err, int niter, char out[16])
{
  double rate = avg_convergence_rate(initial_err, final_err, niter);
  if (rate < 0.0 || !isfinite(rate))
    snprintf(out, 16, "%11s", "n/a");
  else
    snprintf(out, 16, "%11.3e", rate);
}

static double
status_max_error(const struct gkyl_vlasov_lte_correct_status *status, int num_comp)
{
  double max_error = 0.0;
  for (int i=0; i<num_comp; ++i)
    max_error = fmax(max_error, status->error[i]);
  return max_error;
}

static void
print_status_row(const char *case_name, const char *model_name, int cdim, int vdim, int poly_order,
  const struct gkyl_vlasov_lte_correct_status *picard_status,
  const struct gkyl_vlasov_lte_correct_status *anderson_status,
  double picard_time, double anderson_time)
{
  int num_comp = vdim+2;
  double picard_error = status_max_error(picard_status, num_comp);
  double anderson_error = status_max_error(anderson_status, num_comp);
  double picard_time_per_iter = picard_time/picard_status->num_iter;
  double anderson_time_per_iter = anderson_time/anderson_status->num_iter;

  char dim_label[16];
  snprintf(dim_label, sizeof dim_label, "%dx%dv", cdim, vdim);

  printf("%-8s %-10s %5s %2d %7d %7d %10.3e %10.3e %10.3e %10.3e %7s %7s %11.4e %11.4e\n",
    case_name, model_name, dim_label, poly_order,
    picard_status->num_iter, anderson_status->num_iter,
    picard_time, anderson_time, 1.0e3*picard_time_per_iter, 1.0e3*anderson_time_per_iter,
    picard_status->iter_converged ? "fail" : "ok",
    anderson_status->iter_converged ? "fail" : "ok",
    picard_error, anderson_error);

  TEST_CHECK(isfinite(picard_error));
  TEST_CHECK(isfinite(anderson_error));
  TEST_CHECK(isfinite(picard_time));
  TEST_CHECK(isfinite(anderson_time));
  TEST_CHECK(picard_status->num_iter > 0);
  TEST_CHECK(anderson_status->num_iter > 0);
}

static void
print_modal_convergence_table(const char *case_name, const char *model_name,
  int cdim, int vdim, int poly_order, int num_conf_basis,
  const struct gkyl_vlasov_lte_correct_status *picard_status,
  const struct gkyl_vlasov_lte_correct_status *anderson_status,
  const double *initial_err, const double *picard_err, const double *anderson_err)
{
  char dim_label[16];
  snprintf(dim_label, sizeof dim_label, "%dx%dv", cdim, vdim);

  printf("  modal convergence for %-8s %-10s %5s p%d "
    "(scaled max over cells, rate=(final/initial)^(1/iter))\n",
    case_name, model_name, dim_label, poly_order);
  printf("    %-5s %4s %11s %11s %11s %11s %11s %11s\n",
    "comp", "mode", "init", "Perr", "Prate", "Aerr", "Arate", "A/Perr");

  int num_comp = vdim+2;
  for (int d=0; d<num_comp; ++d) {
    char comp_label[16];
    component_label(d, vdim, comp_label);

    for (int k=0; k<num_conf_basis; ++k) {
      int idx = d*num_conf_basis+k;
      char prate[16], arate[16];
      format_rate(initial_err[idx], picard_err[idx], picard_status->num_iter, prate);
      format_rate(initial_err[idx], anderson_err[idx], anderson_status->num_iter, arate);

      double ap_ratio = picard_err[idx] > DBL_MIN ? anderson_err[idx]/picard_err[idx] : 0.0;
      printf("    %-5s %4d %11.3e %11.3e %11s %11.3e %11s %11.3e\n",
        comp_label, k, initial_err[idx], picard_err[idx], prate,
        anderson_err[idx], arate, ap_ratio);

      TEST_CHECK(isfinite(initial_err[idx]));
      TEST_CHECK(isfinite(picard_err[idx]));
      TEST_CHECK(isfinite(anderson_err[idx]));
      TEST_CHECK(isfinite(ap_ratio));
    }
  }
}

static void
run_iter_count_case(enum gkyl_model_id model_id, int cdim, int vdim, int poly_order,
  bool fluctuate)
{
  const int ndim = cdim+vdim;
  const int max_iter = 400;
  const double eps = 1.0e-12;
  const double vmax = model_id == GKYL_MODEL_SR ? 10.0 : 6.0;

  double lower[GKYL_MAX_DIM] = { 0.0 };
  double upper[GKYL_MAX_DIM] = { 0.0 };
  int cells[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) {
    lower[d] = 0.0;
    upper[d] = 1.0;
    cells[d] = fluctuate ? (cdim == 1 ? 4 : 3) : 2;
  }
  for (int d=0; d<vdim; ++d) {
    lower[cdim+d] = -vmax;
    upper[cdim+d] = vmax;
    cells[cdim+d] = vdim == 3 ? 12 : (fluctuate && vdim > 1 ? 12 : 24);
  }

  double confLower[GKYL_MAX_DIM] = { 0.0 };
  double confUpper[GKYL_MAX_DIM] = { 0.0 };
  int confCells[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }

  double velLower[3] = { 0.0 };
  double velUpper[3] = { 0.0 };
  int velCells[3] = { 0 };
  for (int d=0; d<vdim; ++d) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  struct gkyl_rect_grid grid, confGrid, vel_grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d)
    confGhost[d] = 1;
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d)
    ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct lte_moment_ctx moment_ctx = make_moment_ctx(model_id, cdim, vdim,
    fluctuate, confLower, confUpper);

  struct gkyl_array *m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *u = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *temp = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_lte_n, &moment_ctx);
  gkyl_proj_on_basis *proj_u = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_lte_u, &moment_ctx);
  gkyl_proj_on_basis *proj_temp = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_lte_temp, &moment_ctx);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_u, 0.0, &confLocal, u);
  gkyl_proj_on_basis_advance(proj_temp, 0.0, &confLocal, temp);

  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, u, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, temp, (vdim+1)*confBasis.num_basis, &confLocal);

  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&vel_grid,
    &velLocal, &velBasis, inp_vmap, false);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal,
    model_id, vel_map, hamil, gamma_inv, false);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid,
    .conf_basis = &confBasis,
    .vel_basis = &velBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .hamil = hamil,
    .hamil_range = &velLocal,
    .gamma_inv = gamma_inv,
    .model_id = model_id,
    .use_extended_hamil_def = false,
    .use_gpu = false,
    .vel_map = vel_map,
  };
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);

  struct gkyl_array *distf0 = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_picard = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_anderson = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *moms_initial = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_picard = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_anderson = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf0);
  gkyl_array_copy(distf_picard, distf0);
  gkyl_array_copy(distf_anderson, distf0);

  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid,
    .conf_basis = &confBasis,
    .vel_basis = &velBasis,
    .phase_basis = &basis,
    .conf_range = &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .hamil = hamil,
    .hamil_range = &velLocal,
    .gamma_inv = gamma_inv,
    .model_id = model_id,
    .use_extended_hamil_def = false,
    .use_gpu = false,
    .vel_map = vel_map,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf0, moms_initial);

  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid,
    .conf_basis = &confBasis,
    .vel_basis = &velBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .hamil = hamil,
    .hamil_range = &velLocal,
    .gamma_inv = gamma_inv,
    .model_id = model_id,
    .use_picard = true,
    .use_last_converged = true,
    .use_extended_hamil_def = false,
    .use_gpu = false,
    .max_iter = max_iter,
    .eps = eps,
    .vel_map = vel_map,
  };
  gkyl_vlasov_lte_correct *picard_corr = gkyl_vlasov_lte_correct_inew(&inp_corr);
  struct timespec picard_tm = gkyl_wall_clock();
  struct gkyl_vlasov_lte_correct_status picard_status =
    gkyl_vlasov_lte_correct_all_moments(picard_corr, distf_picard, moms, &local, &confLocal);
  double picard_time = gkyl_time_diff_now_sec(picard_tm);
  gkyl_vlasov_lte_correct_release(picard_corr);
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf_picard, moms_picard);

  inp_corr.use_picard = false;
  gkyl_vlasov_lte_correct *anderson_corr = gkyl_vlasov_lte_correct_inew(&inp_corr);
  struct timespec anderson_tm = gkyl_wall_clock();
  struct gkyl_vlasov_lte_correct_status anderson_status =
    gkyl_vlasov_lte_correct_all_moments(anderson_corr, distf_anderson, moms, &local, &confLocal);
  double anderson_time = gkyl_time_diff_now_sec(anderson_tm);
  gkyl_vlasov_lte_correct_release(anderson_corr);
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf_anderson, moms_anderson);

  int num_comp = vdim+2;
  int modal_err_size = num_comp*confBasis.num_basis;
  double *initial_err = gkyl_malloc(sizeof(double[modal_err_size]));
  double *picard_err = gkyl_malloc(sizeof(double[modal_err_size]));
  double *anderson_err = gkyl_malloc(sizeof(double[modal_err_size]));

  calc_modal_scaled_errors(moms_initial, moms, &confLocal, num_comp, confBasis.num_basis,
    initial_err);
  calc_modal_scaled_errors(moms_picard, moms, &confLocal, num_comp, confBasis.num_basis,
    picard_err);
  calc_modal_scaled_errors(moms_anderson, moms, &confLocal, num_comp, confBasis.num_basis,
    anderson_err);

  const char *case_name = fluctuate ? "fluct" : "const";
  const char *model_name = model_id == GKYL_MODEL_SR ? "MJ" : "Maxwellian";
  print_status_row(case_name, model_name, cdim, vdim, poly_order,
    &picard_status, &anderson_status, picard_time, anderson_time);
  print_modal_convergence_table(case_name, model_name, cdim, vdim, poly_order,
    confBasis.num_basis, &picard_status, &anderson_status, initial_err,
    picard_err, anderson_err);

  gkyl_free(anderson_err);
  gkyl_free(picard_err);
  gkyl_free(initial_err);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_array_release(moms_anderson);
  gkyl_array_release(moms_picard);
  gkyl_array_release(moms_initial);
  gkyl_array_release(distf_anderson);
  gkyl_array_release(distf_picard);
  gkyl_array_release(distf0);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_proj_on_basis_release(proj_temp);
  gkyl_proj_on_basis_release(proj_u);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_array_release(moms);
  gkyl_array_release(temp);
  gkyl_array_release(u);
  gkyl_array_release(m0);
}

static void
test_lte_correct_iteration_counts(void)
{
  printf("\nLTE correction iteration counts (eps=1e-12, max_iter=400)\n");
  printf("%-8s %-10s %5s %2s %7s %7s %10s %10s %10s %10s %7s %7s %11s %11s\n",
    "case", "model", "dims", "p", "Pit", "Ait", "Ptime(s)", "Atime(s)",
    "Pms/it", "Ams/it", "Pstat", "Astat", "Perr0", "Aerr0");

  for (int vdim=1; vdim<=3; ++vdim)
    run_iter_count_case(GKYL_MODEL_DEFAULT, 1, vdim, 2, false);
  for (int vdim=1; vdim<=3; ++vdim)
    run_iter_count_case(GKYL_MODEL_SR, 1, vdim, 2, false);
}

static void
test_lte_correct_fluctuating_iteration_counts(void)
{
  printf("\nLTE correction fluctuating moment cases (eps=1e-12, max_iter=400)\n");
  printf("%-8s %-10s %5s %2s %7s %7s %10s %10s %10s %10s %7s %7s %11s %11s\n",
    "case", "model", "dims", "p", "Pit", "Ait", "Ptime(s)", "Atime(s)",
    "Pms/it", "Ams/it", "Pstat", "Astat", "Perr0", "Aerr0");

  for (int poly_order=1; poly_order<=2; ++poly_order) {
    run_iter_count_case(GKYL_MODEL_DEFAULT, 1, 1, poly_order, true);
    run_iter_count_case(GKYL_MODEL_DEFAULT, 2, 2, poly_order, true);
    run_iter_count_case(GKYL_MODEL_DEFAULT, 1, 3, poly_order, true);
    run_iter_count_case(GKYL_MODEL_SR, 1, 1, poly_order, true);
    run_iter_count_case(GKYL_MODEL_SR, 2, 2, poly_order, true);
    run_iter_count_case(GKYL_MODEL_SR, 1, 3, poly_order, true);
  }
}

TEST_LIST = {
  { "lte_correct_iteration_counts", test_lte_correct_iteration_counts },
  { "lte_correct_fluctuating_iteration_counts", test_lte_correct_fluctuating_iteration_counts },
  { NULL, NULL },
};
