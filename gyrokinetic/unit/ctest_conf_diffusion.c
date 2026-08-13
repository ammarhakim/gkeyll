#include <math.h>

#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_conf_diffusion.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static void
copy_periodic_ghosts(const struct gkyl_range *local,
  const struct gkyl_range *local_ext, struct gkyl_array *field)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local_ext);
  while (gkyl_range_iter_next(&iter)) {
    int src_idx[GKYL_MAX_CDIM];
    bool is_ghost = false;
    for (int d=0; d<local->ndim; ++d) {
      src_idx[d] = iter.idx[d];
      if (src_idx[d] < local->lower[d]) {
        src_idx[d] = local->upper[d];
        is_ghost = true;
      }
      else if (src_idx[d] > local->upper[d]) {
        src_idx[d] = local->lower[d];
        is_ghost = true;
      }
    }
    if (!is_ghost)
      continue;

    double *dst = gkyl_array_fetch(field,
      gkyl_range_idx(local_ext, iter.idx));
    const double *src = gkyl_array_cfetch(field,
      gkyl_range_idx(local_ext, src_idx));
    for (int k=0; k<field->ncomp; ++k)
      dst[k] = src[k];
  }
}

static void
set_cartesian_geometry(int cdim, const struct gkyl_basis *basis,
  const struct gkyl_range *range, const struct gkyl_array *D,
  struct gkyl_array *tensor, struct gkyl_array *jacobgeo_inv)
{
  const double modal_norm = pow(sqrt(2.0), cdim);
  gkyl_array_clear(tensor, 0.0);
  gkyl_array_clear(jacobgeo_inv, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *D_i = gkyl_array_cfetch(D, loc);
    double *K_i = gkyl_array_fetch(tensor, loc);
    double *Jinv_i = gkyl_array_fetch(jacobgeo_inv, loc);
    Jinv_i[0] = modal_norm;
    for (int d=0; d<cdim; ++d)
      for (int k=0; k<basis->num_basis; ++k)
        K_i[(d*cdim+d)*basis->num_basis+k] = D_i[k];
  }
}

static void
test_conf_diffusion(int cdim)
{
  const int cells[] = { 6, 5, 4 };
  const double lower[] = { 0.0, 0.0, 0.0 };
  const double upper[] = { 1.0, 1.25, 1.5 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);

  const int ghost[] = { 1, 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *diffusion_coeff = gkyl_array_new(GKYL_DOUBLE,
    basis.num_basis, local_ext.volume);
  struct gkyl_array *diffusion_tensor = gkyl_array_new(GKYL_DOUBLE,
    cdim*cdim*basis.num_basis, local_ext.volume);
  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE,
    basis.num_basis, local_ext.volume);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE,
    basis.num_basis, local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE,
    basis.num_basis, local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1,
    local_ext.volume);
  const int periodic_flags[2*GKYL_MAX_CDIM] = { 0 };

  struct gkyl_dg_updater_conf_diffusion *up =
    gkyl_dg_updater_conf_diffusion_new(&grid, &basis, &local_ext,
      diffusion_tensor, jacobgeo_inv, periodic_flags, false);

  double modal_norm = pow(sqrt(2.0),cdim);

  // A spatially uniform field must remain in the null space even when D is a
  // discontinuous, spatially nonuniform DG field.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    double *D = gkyl_array_fetch(diffusion_coeff, loc);
    double *fin = gkyl_array_fetch(f, loc);
    int idx_sum = 0;
    for (int d=0; d<cdim; ++d)
      idx_sum += iter.idx[d];
    D[0] = modal_norm*(1.0+0.01*idx_sum);
    fin[0] = modal_norm*2.5;
    for (int k=1; k<basis.num_basis; ++k) {
      D[k] = 0.01*cos(k+idx_sum);
      fin[k] = 0.0;
    }
  }
  set_cartesian_geometry(cdim, &basis, &local_ext, diffusion_coeff,
    diffusion_tensor, jacobgeo_inv);

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double rdx2_sum = 0.0;
  for (int d=0; d<cdim; ++d)
    rdx2_sum += pow(2.0/grid.dx[d],2.0);

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *rhs_i = gkyl_array_cfetch(rhs, loc);
    const double *D = gkyl_array_cfetch(diffusion_coeff, loc);
    const double *cfl = gkyl_array_cfetch(cflrate, loc);
    double expected_cfl = 4.0*(D[0]/modal_norm)*rdx2_sum;

    for (int k=0; k<basis.num_basis; ++k)
      TEST_CHECK(gkyl_compare(rhs_i[k], 0.0, 1e-11));
    TEST_CHECK(gkyl_compare(cfl[0], expected_cfl, 1e-12));
  }

  // Exercise every modal input in the full div(D grad(f)) update. With
  // periodic ghosts, the sum of the cell-average RHS coefficients must vanish.
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    double phase = 0.0;
    for (int d=0; d<cdim; ++d) {
      double x = grid.lower[d]
        +(iter.idx[d]-local.lower[d]+0.5)*grid.dx[d];
      phase += 2.0*M_PI*x/(grid.upper[d]-grid.lower[d]);
    }
    double *D = gkyl_array_fetch(diffusion_coeff, loc);
    double *fin = gkyl_array_fetch(f, loc);
    D[0] = modal_norm*(1.0+0.1*cos(phase));
    fin[0] = modal_norm*(1.0+0.2*sin(phase));
    for (int k=1; k<basis.num_basis; ++k) {
      D[k] = 0.01*cos((k+1)*phase);
      fin[k] = 0.03*sin((k+1)*phase);
    }
  }
  copy_periodic_ghosts(&local, &local_ext, diffusion_coeff);
  copy_periodic_ghosts(&local, &local_ext, f);
  set_cartesian_geometry(cdim, &basis, &local_ext, diffusion_coeff,
    diffusion_tensor, jacobgeo_inv);

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double rhs_avg_sum = 0.0, f_dot_rhs = 0.0;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *rhs_i = gkyl_array_cfetch(rhs, loc);
    const double *fin = gkyl_array_cfetch(f, loc);
    rhs_avg_sum += rhs_i[0];
    for (int k=0; k<basis.num_basis; ++k)
      f_dot_rhs += fin[k]*rhs_i[k];
  }
  TEST_CHECK(gkyl_compare(rhs_avg_sum, 0.0, 1e-10));
  TEST_CHECK(f_dot_rhs < 0.0);

  struct gkyl_dg_updater_conf_diffusion_tm tm =
    gkyl_dg_updater_conf_diffusion_get_tm(up);
  TEST_CHECK(tm.diffusion_tm >= 0.0);

  gkyl_dg_updater_conf_diffusion_release(up);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs);
  gkyl_array_release(f);
  gkyl_array_release(jacobgeo_inv);
  gkyl_array_release(diffusion_tensor);
  gkyl_array_release(diffusion_coeff);
}

static void test_conf_diffusion_1x_p1(void) { test_conf_diffusion(1); }
static void test_conf_diffusion_2x_p1(void) { test_conf_diffusion(2); }
static void test_conf_diffusion_3x_p1(void) { test_conf_diffusion(3); }

struct manufactured_ctx {
  int cdim;
  double density_amp;
  double diffusion_amp;
};

static void
eval_manufactured_density(double t, const double *xn, double *fout, void *vctx)
{
  struct manufactured_ctx *ctx = vctx;
  double cos_prod = 1.0;
  for (int d=0; d<ctx->cdim; ++d)
    cos_prod *= cos(xn[d]);
  fout[0] = 1.0+ctx->density_amp*cos_prod;
}

static void
eval_manufactured_diffusion(double t, const double *xn, double *fout,
  void *vctx)
{
  struct manufactured_ctx *ctx = vctx;
  double sin_sum = 0.0;
  for (int d=0; d<ctx->cdim; ++d)
    sin_sum += sin(xn[d]);
  fout[0] = 1.0+ctx->diffusion_amp*sin_sum/ctx->cdim;
}

static void
eval_manufactured_rhs(double t, const double *xn, double *fout, void *vctx)
{
  struct manufactured_ctx *ctx = vctx;
  double cos_prod = 1.0, sin_sum = 0.0, grad_D_dot_grad_n = 0.0;
  for (int d=0; d<ctx->cdim; ++d) {
    cos_prod *= cos(xn[d]);
    sin_sum += sin(xn[d]);
  }
  const double D = 1.0+ctx->diffusion_amp*sin_sum/ctx->cdim;
  for (int d=0; d<ctx->cdim; ++d) {
    double other_cos_prod = 1.0;
    for (int j=0; j<ctx->cdim; ++j)
      if (j != d) other_cos_prod *= cos(xn[j]);
    const double dD = ctx->diffusion_amp*cos(xn[d])/ctx->cdim;
    const double dn = -ctx->density_amp*sin(xn[d])*other_cos_prod;
    grad_D_dot_grad_n += dD*dn;
  }
  const double lap_n = -ctx->cdim*ctx->density_amp*cos_prod;
  fout[0] = grad_D_dot_grad_n+D*lap_n;
}

static double
manufactured_error(int cdim, int ncells)
{
  int cells[] = { ncells, ncells, ncells };
  const double lower[] = { 0.0, 0.0, 0.0 };
  const double upper[] = { 2.0*M_PI, 2.0*M_PI, 2.0*M_PI };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  const int ghost[] = { 1, 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *D = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *tensor = gkyl_array_new(GKYL_DOUBLE,
    cdim*cdim*basis.num_basis, local_ext.volume);
  struct gkyl_array *jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE,
    basis.num_basis, local_ext.volume);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs_exact = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1,
    local_ext.volume);
  const int periodic_flags[2*GKYL_MAX_CDIM] = { 0 };

  struct manufactured_ctx ctx = {
    .cdim = cdim, .density_amp = 0.2, .diffusion_amp = 0.3,
  };
  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis, 2,
    1, eval_manufactured_density, &ctx);
  gkyl_proj_on_basis *proj_D = gkyl_proj_on_basis_new(&grid, &basis, 2,
    1, eval_manufactured_diffusion, &ctx);
  gkyl_proj_on_basis *proj_rhs = gkyl_proj_on_basis_new(&grid, &basis, 3,
    1, eval_manufactured_rhs, &ctx);
  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, f);
  gkyl_proj_on_basis_advance(proj_D, 0.0, &local, D);
  gkyl_proj_on_basis_advance(proj_rhs, 0.0, &local, rhs_exact);
  copy_periodic_ghosts(&local, &local_ext, f);
  copy_periodic_ghosts(&local, &local_ext, D);
  set_cartesian_geometry(cdim, &basis, &local_ext, D, tensor,
    jacobgeo_inv);

  struct gkyl_dg_updater_conf_diffusion *up =
    gkyl_dg_updater_conf_diffusion_new(&grid, &basis, &local_ext, tensor,
      jacobgeo_inv, periodic_flags, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double err_sq = 0.0, exact_sq = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *rhs_i = gkyl_array_cfetch(rhs, loc);
    const double *exact_i = gkyl_array_cfetch(rhs_exact, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      err_sq += pow(rhs_i[k]-exact_i[k], 2);
      exact_sq += pow(exact_i[k], 2);
    }
  }

  gkyl_dg_updater_conf_diffusion_release(up);
  gkyl_proj_on_basis_release(proj_rhs);
  gkyl_proj_on_basis_release(proj_D);
  gkyl_proj_on_basis_release(proj_f);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs_exact);
  gkyl_array_release(rhs);
  gkyl_array_release(f);
  gkyl_array_release(jacobgeo_inv);
  gkyl_array_release(tensor);
  gkyl_array_release(D);
  return sqrt(err_sq/exact_sq);
}

static void
test_conf_diffusion_manufactured(int cdim)
{
  const double coarse = manufactured_error(cdim, 8);
  const double fine = manufactured_error(cdim, 16);
  const double order = log(coarse/fine)/log(2.0);
  TEST_MSG("cdim=%d: coarse error=%g fine error=%g order=%g",
    cdim, coarse, fine, order);
  TEST_CHECK(fine < coarse);
  TEST_CHECK(fine < 0.05);
  // This compares the second-derivative DG residual itself. For p1 its L2
  // consistency is first order; the evolved solution may converge faster.
  TEST_CHECK(order > 0.9);
}

static void test_conf_diffusion_manufactured_2x_p1(void)
{ test_conf_diffusion_manufactured(2); }
static void test_conf_diffusion_manufactured_3x_p1(void)
{ test_conf_diffusion_manufactured(3); }

struct dirichlet_1x_ctx {
  double lower, upper, length, loss, diffusion;
};

static void
eval_dirichlet_steady_1x(double t, const double *z, double *fout, void *vctx)
{
  struct dirichlet_1x_ctx *ctx = vctx;
  if (ctx->loss == 0.0) {
    fout[0] = ctx->lower+(ctx->upper-ctx->lower)*z[0]/ctx->length;
    return;
  }

  const double k = sqrt(ctx->loss/ctx->diffusion);
  const double sinh_kL = sinh(k*ctx->length);
  fout[0] = ctx->lower*sinh(k*(ctx->length-z[0]))/sinh_kL
    +ctx->upper*sinh(k*z[0])/sinh_kL;
}

static void
set_dirichlet_ghost_1x(const struct gkyl_basis *basis,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  double lower, double upper, struct gkyl_array *f,
  struct gkyl_array *tensor)
{
  const int nb = basis->num_basis;
  const double values[2] = { lower, upper };
  for (int edge=0; edge<2; ++edge) {
    int skin_idx[] = { edge == 0 ? local->lower[0] : local->upper[0] };
    int ghost_idx[] = { edge == 0 ? local->lower[0]-1 : local->upper[0]+1 };
    const double *fs = gkyl_array_cfetch(f,
      gkyl_range_idx(local_ext, skin_idx));
    double *fg = gkyl_array_fetch(f, gkyl_range_idx(local_ext, ghost_idx));
    basis->flip_odd_sign(0, fs, fg);
    for (int k=0; k<nb; ++k) fg[k] = -fg[k];
    fg[0] += 2.0*sqrt(2.0)*values[edge];

    const double *Ks = gkyl_array_cfetch(tensor,
      gkyl_range_idx(local_ext, skin_idx));
    double *Kg = gkyl_array_fetch(tensor,
      gkyl_range_idx(local_ext, ghost_idx));
    basis->flip_odd_sign(0, Ks, Kg);
  }
}

static double
dirichlet_steady_residual_1x(int ncells, double loss)
{
  const int cells[] = { ncells }, ghost[] = { 1 };
  const double lower[] = { 0.0 }, upper[] = { 1.0 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *D = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *K = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *Jinv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1,
    local_ext.volume);

  struct dirichlet_1x_ctx ctx = { .lower = 1.0, .upper = 1.7,
    .length = 1.0, .loss = loss, .diffusion = 0.4 };
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 3, 1,
    eval_dirichlet_steady_1x, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, f);
  gkyl_array_shiftc(D, sqrt(2.0)*ctx.diffusion, 0);
  set_cartesian_geometry(1, &basis, &local_ext, D, K, Jinv);
  set_dirichlet_ghost_1x(&basis, &local, &local_ext, ctx.lower, ctx.upper,
    f, K);

  const int ghost_surface_flags[2*GKYL_MAX_CDIM] = { 0 };
  struct gkyl_dg_updater_conf_diffusion *up =
    gkyl_dg_updater_conf_diffusion_new(&grid, &basis, &local_ext, K, Jinv,
      ghost_surface_flags, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double residual_sq = 0.0, scale_sq = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *fi = gkyl_array_cfetch(f, loc);
    const double *ri = gkyl_array_cfetch(rhs, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      const double residual = ri[k]-loss*fi[k];
      residual_sq += residual*residual;
      scale_sq += fi[k]*fi[k];
    }
  }

  gkyl_dg_updater_conf_diffusion_release(up);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs);
  gkyl_array_release(f);
  gkyl_array_release(Jinv);
  gkyl_array_release(K);
  gkyl_array_release(D);
  return sqrt(residual_sq/scale_sq);
}

static void
test_conf_diffusion_dirichlet_steady_1x_p1(void)
{
  const double linear = dirichlet_steady_residual_1x(16, 0.0);
  const double coarse = dirichlet_steady_residual_1x(16, 0.3);
  const double fine = dirichlet_steady_residual_1x(32, 0.3);
  TEST_MSG("Dirichlet steady residuals: linear=%g cosh coarse=%g fine=%g",
    linear, coarse, fine);
  TEST_CHECK(linear < 1e-10);
  TEST_CHECK(fine < coarse);
}

static void
test_conf_diffusion_zero_flux_2x_p1(void)
{
  const int cdim = 2, cells[] = { 7, 6 };
  const double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.25 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  const int ghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *D = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *tensor = gkyl_array_new(GKYL_DOUBLE,
    cdim*cdim*basis.num_basis, local_ext.volume);
  struct gkyl_array *Jinv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1,
    local_ext.volume);
  const int zero_flux_flags[2*GKYL_MAX_CDIM] = { 1, 1, 1, 1 };
  const double modal_norm = pow(sqrt(2.0), cdim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    double x = lower[0]+(iter.idx[0]-local.lower[0]+0.5)*grid.dx[0];
    double y = lower[1]+(iter.idx[1]-local.lower[1]+0.5)*grid.dx[1];
    double *D_i = gkyl_array_fetch(D, loc);
    double *f_i = gkyl_array_fetch(f, loc);
    D_i[0] = modal_norm*(0.7+0.1*x+0.05*y);
    f_i[0] = modal_norm*(1.0+0.2*x*x-0.1*y);
    for (int k=1; k<basis.num_basis; ++k) {
      D_i[k] = 0.01*cos(k+x-y);
      f_i[k] = 0.02*sin(k+x+y);
    }
  }
  set_cartesian_geometry(cdim, &basis, &local_ext, D, tensor, Jinv);
  // Add symmetric cross terms so the natural boundary condition exercises
  // the full normal flux K^{ij} d_j f, as required by mapped geometry.
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *D_i = gkyl_array_cfetch(D, loc);
    double *K_i = gkyl_array_fetch(tensor, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      K_i[(0*cdim+1)*basis.num_basis+k] = 0.12*D_i[k];
      K_i[(1*cdim+0)*basis.num_basis+k] = 0.12*D_i[k];
    }
  }

  struct gkyl_dg_updater_conf_diffusion *up =
    gkyl_dg_updater_conf_diffusion_new(&grid, &basis, &local_ext, tensor,
      Jinv, zero_flux_flags, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double rhs_avg_sum = 0.0;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *rhs_i = gkyl_array_cfetch(rhs, loc);
    rhs_avg_sum += rhs_i[0];
  }
  TEST_MSG("zero-flux cell-average RHS sum=%g", rhs_avg_sum);
  TEST_CHECK(fabs(rhs_avg_sum) < 1e-11);

  gkyl_dg_updater_conf_diffusion_release(up);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs);
  gkyl_array_release(f);
  gkyl_array_release(Jinv);
  gkyl_array_release(tensor);
  gkyl_array_release(D);
}

struct mapped_ctx {
  int cdim;
  double map_amp, shear, diffusion, density_amp;
};

static void
eval_mapped_density(double t, const double *z, double *fout, void *vctx)
{
  struct mapped_ctx *ctx = vctx;
  const double X = z[0]+ctx->map_amp*sin(z[0])+ctx->shear*sin(z[1]);
  double prod = cos(X)*cos(z[1]);
  if (ctx->cdim == 3) prod *= cos(z[2]);
  fout[0] = 1.0+ctx->density_amp*prod;
}

static void
eval_mapped_jacobgeo_inv(double t, const double *z, double *fout,
  void *vctx)
{
  struct mapped_ctx *ctx = vctx;
  fout[0] = 1.0/(1.0+ctx->map_amp*cos(z[0]));
}

static void
eval_mapped_tensor(double t, const double *z, double *fout, void *vctx)
{
  struct mapped_ctx *ctx = vctx;
  const double A = 1.0+ctx->map_amp*cos(z[0]);
  const double B = ctx->shear*cos(z[1]);
  for (int i=0; i<ctx->cdim*ctx->cdim; ++i) fout[i] = 0.0;
  fout[0] = ctx->diffusion*(1.0+B*B)/A;
  fout[1] = fout[ctx->cdim] = -ctx->diffusion*B;
  fout[ctx->cdim+1] = ctx->diffusion*A;
  if (ctx->cdim == 3)
    fout[2*ctx->cdim+2] = ctx->diffusion*A;
}

static void
eval_mapped_rhs(double t, const double *z, double *fout, void *vctx)
{
  struct mapped_ctx *ctx = vctx;
  const double X = z[0]+ctx->map_amp*sin(z[0])+ctx->shear*sin(z[1]);
  const double J = 1.0+ctx->map_amp*cos(z[0]);
  double prod = cos(X)*cos(z[1]);
  if (ctx->cdim == 3) prod *= cos(z[2]);
  fout[0] = -ctx->cdim*J*ctx->diffusion*ctx->density_amp*prod;
}

static double
mapped_manufactured_error(int cdim, int ncells, double shear,
  double *conservation)
{
  int cells[] = { ncells, ncells, ncells };
  const double lower[] = { 0.0, 0.0, 0.0 };
  const double upper[] = { 2.0*M_PI, 2.0*M_PI, 2.0*M_PI };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  const int ghost[] = { 1, 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *K = gkyl_array_new(GKYL_DOUBLE,
    cdim*cdim*basis.num_basis, local_ext.volume);
  struct gkyl_array *Jinv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *rhs_exact = gkyl_array_new(GKYL_DOUBLE, basis.num_basis,
    local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1,
    local_ext.volume);
  const int periodic_flags[2*GKYL_MAX_CDIM] = { 0 };

  struct mapped_ctx ctx = { .cdim = cdim, .map_amp = 0.2,
    .shear = shear, .diffusion = 0.4, .density_amp = 0.1 };
  gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis, 3,
    1, eval_mapped_density, &ctx);
  gkyl_proj_on_basis *proj_K = gkyl_proj_on_basis_new(&grid, &basis, 3,
    cdim*cdim, eval_mapped_tensor, &ctx);
  gkyl_proj_on_basis *proj_Jinv = gkyl_proj_on_basis_new(&grid, &basis, 3,
    1, eval_mapped_jacobgeo_inv, &ctx);
  gkyl_proj_on_basis *proj_rhs = gkyl_proj_on_basis_new(&grid, &basis, 3,
    1, eval_mapped_rhs, &ctx);
  gkyl_proj_on_basis_advance(proj_f, 0.0, &local, f);
  gkyl_proj_on_basis_advance(proj_K, 0.0, &local, K);
  gkyl_proj_on_basis_advance(proj_Jinv, 0.0, &local, Jinv);
  gkyl_proj_on_basis_advance(proj_rhs, 0.0, &local, rhs_exact);
  copy_periodic_ghosts(&local, &local_ext, f);
  copy_periodic_ghosts(&local, &local_ext, K);
  copy_periodic_ghosts(&local, &local_ext, Jinv);

  struct gkyl_dg_updater_conf_diffusion *up =
    gkyl_dg_updater_conf_diffusion_new(&grid, &basis, &local_ext, K, Jinv,
      periodic_flags, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(up, &local, f, cflrate, rhs);

  double err_sq = 0.0, exact_sq = 0.0, rhs_sum = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *rhs_i = gkyl_array_cfetch(rhs, loc);
    const double *exact_i = gkyl_array_cfetch(rhs_exact, loc);
    rhs_sum += rhs_i[0];
    for (int k=0; k<basis.num_basis; ++k) {
      err_sq += pow(rhs_i[k]-exact_i[k], 2);
      exact_sq += pow(exact_i[k], 2);
    }
  }
  *conservation = fabs(rhs_sum);

  gkyl_dg_updater_conf_diffusion_release(up);
  gkyl_proj_on_basis_release(proj_rhs);
  gkyl_proj_on_basis_release(proj_Jinv);
  gkyl_proj_on_basis_release(proj_K);
  gkyl_proj_on_basis_release(proj_f);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs_exact);
  gkyl_array_release(rhs);
  gkyl_array_release(Jinv);
  gkyl_array_release(K);
  gkyl_array_release(f);
  return sqrt(err_sq/exact_sq);
}

static void
test_conf_diffusion_mapped(int cdim, double shear)
{
  double cons_coarse = 0.0, cons_fine = 0.0;
  const double coarse = mapped_manufactured_error(cdim, 8, shear,
    &cons_coarse);
  const double fine = mapped_manufactured_error(cdim, 16, shear, &cons_fine);
  const double order = log(coarse/fine)/log(2.0);
  TEST_MSG("cdim=%d shear=%g: coarse=%g fine=%g order=%g cons=%g",
    cdim, shear, coarse, fine, order, cons_fine);
  TEST_CHECK(fine < coarse);
  TEST_CHECK(fine < 0.06);
  TEST_CHECK(order > 0.8);
  TEST_CHECK(cons_coarse < 1e-11);
  TEST_CHECK(cons_fine < 1e-10);
}

static void test_conf_diffusion_mapped_orthogonal_2x_p1(void)
{ test_conf_diffusion_mapped(2, 0.0); }
static void test_conf_diffusion_mapped_sheared_2x_p1(void)
{ test_conf_diffusion_mapped(2, 0.25); }
static void test_conf_diffusion_mapped_sheared_3x_p1(void)
{ test_conf_diffusion_mapped(3, 0.25); }

TEST_LIST = {
  { "conf_diffusion_1x_p1", test_conf_diffusion_1x_p1 },
  { "conf_diffusion_2x_p1", test_conf_diffusion_2x_p1 },
  { "conf_diffusion_3x_p1", test_conf_diffusion_3x_p1 },
  { "conf_diffusion_manufactured_2x_p1",
    test_conf_diffusion_manufactured_2x_p1 },
  { "conf_diffusion_manufactured_3x_p1",
    test_conf_diffusion_manufactured_3x_p1 },
  { "conf_diffusion_dirichlet_steady_1x_p1",
    test_conf_diffusion_dirichlet_steady_1x_p1 },
  { "conf_diffusion_zero_flux_2x_p1",
    test_conf_diffusion_zero_flux_2x_p1 },
  { "conf_diffusion_mapped_orthogonal_2x_p1",
    test_conf_diffusion_mapped_orthogonal_2x_p1 },
  { "conf_diffusion_mapped_sheared_2x_p1",
    test_conf_diffusion_mapped_sheared_2x_p1 },
  { "conf_diffusion_mapped_sheared_3x_p1",
    test_conf_diffusion_mapped_sheared_3x_p1 },
  { NULL, NULL },
};
