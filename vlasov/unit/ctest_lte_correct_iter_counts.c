#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_velocity_map.h>

#include <math.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct lte_moment_ctx {
  int vdim;
  double n;
  double u[3];
  double temp;
};

static void
eval_lte_n(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  fout[0] = m->n;
}

static void
eval_lte_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  for (int d=0; d<m->vdim; ++d)
    fout[d] = m->u[d];
}

static void
eval_lte_temp(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct lte_moment_ctx *m = ctx;
  fout[0] = m->temp;
}

static struct lte_moment_ctx
make_moment_ctx(enum gkyl_model_id model_id, int vdim)
{
  struct lte_moment_ctx m = { 0 };
  m.vdim = vdim;

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

static double
status_max_error(const struct gkyl_vlasov_lte_correct_status *status, int num_comp)
{
  double max_error = 0.0;
  for (int i=0; i<num_comp; ++i)
    max_error = fmax(max_error, status->error[i]);
  return max_error;
}

static void
print_status_row(const char *model_name, int vdim,
  const struct gkyl_vlasov_lte_correct_status *picard_status,
  const struct gkyl_vlasov_lte_correct_status *anderson_status,
  double picard_time, double anderson_time)
{
  int num_comp = vdim+2;
  double picard_error = status_max_error(picard_status, num_comp);
  double anderson_error = status_max_error(anderson_status, num_comp);
  double picard_time_per_iter = picard_time/picard_status->num_iter;
  double anderson_time_per_iter = anderson_time/anderson_status->num_iter;

  printf("%-10s %2dv %7d %7d %10.3e %10.3e %10.3e %10.3e %7s %7s %11.4e %11.4e\n",
    model_name, vdim, picard_status->num_iter, anderson_status->num_iter,
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
run_iter_count_case(enum gkyl_model_id model_id, int vdim)
{
  const int cdim = 1;
  const int ndim = cdim+vdim;
  const int poly_order = 2;
  const int max_iter = 200;
  const double eps = 1.0e-12;
  const double vmax = model_id == GKYL_MODEL_SR ? 10.0 : 6.0;

  double lower[GKYL_MAX_DIM] = { 0.0 };
  double upper[GKYL_MAX_DIM] = { 0.0 };
  int cells[GKYL_MAX_DIM] = { 0 };
  lower[0] = 0.1;
  upper[0] = 1.0;
  cells[0] = 2;
  for (int d=0; d<vdim; ++d) {
    lower[1+d] = -vmax;
    upper[1+d] = vmax;
    cells[1+d] = vdim == 3 ? 12 : 24;
  }

  double confLower[] = { lower[0] };
  double confUpper[] = { upper[0] };
  int confCells[] = { cells[0] };

  double velLower[3] = { 0.0 };
  double velUpper[3] = { 0.0 };
  int velCells[3] = { 0 };
  for (int d=0; d<vdim; ++d) {
    velLower[d] = lower[1+d];
    velUpper[d] = upper[1+d];
    velCells[d] = cells[1+d];
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

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[GKYL_MAX_DIM] = { 0 };
  ghost[0] = confGhost[0];
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct lte_moment_ctx moment_ctx = make_moment_ctx(model_id, vdim);

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

  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf0);
  gkyl_array_copy(distf_picard, distf0);
  gkyl_array_copy(distf_anderson, distf0);

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

  inp_corr.use_picard = false;
  gkyl_vlasov_lte_correct *anderson_corr = gkyl_vlasov_lte_correct_inew(&inp_corr);
  struct timespec anderson_tm = gkyl_wall_clock();
  struct gkyl_vlasov_lte_correct_status anderson_status =
    gkyl_vlasov_lte_correct_all_moments(anderson_corr, distf_anderson, moms, &local, &confLocal);
  double anderson_time = gkyl_time_diff_now_sec(anderson_tm);
  gkyl_vlasov_lte_correct_release(anderson_corr);

  print_status_row(model_id == GKYL_MODEL_SR ? "MJ" : "Maxwellian",
    vdim, &picard_status, &anderson_status, picard_time, anderson_time);

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
  printf("\nLTE correction iteration counts (p2, eps=1e-12, max_iter=200)\n");
  printf("%-10s %3s %7s %7s %10s %10s %10s %10s %7s %7s %11s %11s\n",
    "model", "vdim", "Pit", "Ait", "Ptime(s)", "Atime(s)",
    "Pms/it", "Ams/it", "Pstat", "Astat", "Perr", "Aerr");

  for (int vdim=1; vdim<=3; ++vdim)
    run_iter_count_case(GKYL_MODEL_DEFAULT, vdim);
  for (int vdim=1; vdim<=3; ++vdim)
    run_iter_count_case(GKYL_MODEL_SR, vdim);
}

TEST_LIST = {
  { "lte_correct_iteration_counts", test_lte_correct_iteration_counts },
  { NULL, NULL },
};
