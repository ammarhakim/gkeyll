// Tests for the Vlasov FPO routines to compute moments and boundary corrections
// and resulting corrections to drag and diffusion coefficients.
//
// Included are tests for p=1 hybrid and p=2 serendipity.
// As shown in Rodman 2025 PhD Thesis, the p=2 computation enforces
// conservation of momentum and energy, so this result is taken as the
// "correct" result for the hybrid basis.

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>

#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_fpo_vlasov_coeff_correct.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_fpo_vlasov.h>
#include <gkyl_dg_updater_moment.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

// Allocate array
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  if (use_gpu) {
    struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
    return a;
  }
  else {
    struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
    return a;
  }
}

static inline double
bump_maxwellian(double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA, double bUx, double bUy, double bUz, double bS, double bVt)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  double bv2 = (vx - bUx)*(vx - bUx) + (vy - bUy)*(vy - bUy) + (vz - bUz)*(vz - bUz);
  return n/pow(sqrt(2*M_PI*vt*vt), 3)*exp(-v2/(2*vt*vt)) + n/pow(sqrt(2*M_PI*bVt*bVt), 3)*exp(-bv2/(2*bVt*bVt))*(bA*bA)/(bv2 + bS*bS);
}

void
eval_distf_square(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  // double v2 = vx*vx+vy*vy+vz*vz;
  // fout[0] = 2.5/pow(sqrt(2*M_PI*1.0*1.0), 3)*exp(-v2/(2.0*1.0*1.0));
  double width = 2.0; // corresponds to a final vth of 2/3
  if(vx>-width && vx<width && vy>-width && vy<width && vz>-width && vz<width) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_distf_bump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = bump_maxwellian(1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vth0, 
    sqrt(0.15), 4.0*app->vth0, 0.0, 0.0, 0.14, 3.0*app->vth0);
}

void
eval_drag_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vy*vz;
  fout[1] = vx*vz;
  fout[2] = vy*vz;
}

void
eval_diff_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vx*vx;
  fout[1] = vx*vy;
  fout[2] = vx*vz;
  fout[3] = vy*vx;
  fout[4] = vy*vy;
  fout[5] = vz*vz;
  fout[6] = vz*vx;
  fout[7] = vz*vy;
  fout[8] = vz*vz;
}

struct fpo_ctx
create_ctx() {
  struct fpo_ctx ctx = {
    .n0 = 1.0,
    .ux0 = 0.0,
    .uy0 = 0.0,
    .uz0 = 0.0,
    .vth0 = 1.0,
    .gamma0 = 1.0,
  };
return ctx;
}

void test_1x3v_square(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_square = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_square, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 5.684341886080801e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -6.239229700423411e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 4.263256414560601e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 1.032840991473699e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 2.842170943040401e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -7.363988580231177e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 7.542472332656512e+02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -2.076915513147054e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -1.776356839400249e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.949759281382314e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.332267629550187e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -3.227628098355307e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], -8.881784197001246e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 2.301246431322241e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -7.856742013183861e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 1.583657215757345e-16, 1e-14) );
  }

  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -4.074735897754556e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 1.421085471520200e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -3.508797065679343e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 2.486899575160351e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -7.874465007375168e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 8.221294842595601e+02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -3.027827215005209e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 9.396855349118749e-34, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.273354968048298e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -4.440892098500623e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 1.096499083024794e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], -7.771561172376091e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 2.460770314804738e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -8.563848794370413e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 2.522008558586881e-16, 1e-14) );
  }

  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 2.842170943040401e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -5.932662803013587e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 4.102189474635924e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -8.526512829121202e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.283714983449953e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 1.266898361936780e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -2.288509781858262e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[8], -3.487259985304282e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 8.447569012575273e+02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], 5.122523680855253e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[11], -1.926714405166172e-13, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[17], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[20], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[23], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -8.881784197001260e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.853957125941748e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.281934210823729e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.664535259100378e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 4.011609323281108e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -3.959057381051981e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 2.388673753480699e-33, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 7.151593068307076e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[8], 1.089768745407589e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -8.799551054765915e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], -1.167851871009217e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[11], 4.944422314733905e-16, 1e-14) );
  }

  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -5.329070518200751e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 1.418697319828291e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -8.475692026015594e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.131628207280301e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.238225742785780e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -8.464166370424538e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -7.105427357601002e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -1.089352990289268e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[8], -1.290542022163894e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 8.447569012575279e+02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], -3.574851739306521e-14, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[11], -2.081151197511259e-13, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[17], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[20], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[23], 0.000000000000000e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 1.665334536937738e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -4.433429124463419e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 2.648653758129907e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 6.661338147750951e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.869455446205570e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 2.645051990757685e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 2.220446049250317e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 3.404228094653969e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[8], 4.032943819262180e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -8.799551054765931e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], 3.089548240857505e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[11], 6.553138901661876e-16, 1e-14) );
  }
}

void test_1x3v_bump(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_bump = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_bump, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.117608256804482e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -9.043347143054705e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 2.164934898019055e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 7.959408156657205e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 1.613292832658431e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -5.629826622284888e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.052992260884340e+01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -2.797617919019495e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 8.071672647588428e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 6.439424529985747e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 2.058290061827950e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 3.493353320748819e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -2.995108501491206e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 1.200830000828715e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.228669059035373e-02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -4.453443140625313e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.272705566928377e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 1.475248667780351e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 4.163336342344337e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 7.422648601811846e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 4.163336342344337e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 8.628671013942246e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -9.754905928558748e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -9.432063960551569e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -2.469249206173595e-02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.013569553329405e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -2.120423238333841e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -7.346174974803201e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], -1.040750456965584e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 5.621165337886929e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -7.070579733856278e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -3.673393749367210e-17, 1e-14) );
  }

  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.144265799633754e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -1.535102825108805e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 1.745565497701662e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 2.693554579729008e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 2.864462139706703e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -8.025342156315925e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.043430781841897e+01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -5.349764606312907e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.289521729888946e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -1.944063469141467e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -4.294457042579303e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 4.305397722447073e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -7.953639374717880e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 5.403939358927745e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.715808691955579e-02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 3.892519044015941e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.529011550101355e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 2.596214490289867e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 5.664956351236761e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 3.685592689916960e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 5.664956351236761e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 1.280540214833638e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -5.214826735047042e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -2.780865175746028e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -6.646442379735286e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.866426810682822e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.656301161573374e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -2.329426572994407e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], -2.749660043071186e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 8.069949608969555e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -6.876279199501303e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 8.832320457228894e-18, 1e-14) );
  }

  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.168537534358513e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 4.176043644886353e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 2.384301307574831e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -3.179748131465487e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -3.052165688974618e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 1.872028733900613e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -5.603156827405087e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 6.020136031685679e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[8], 9.828299788862179e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 2.044372332849935e+01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], -1.279929899127277e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[11], -3.657199000944288e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.718760928113808e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -5.492072533415933e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -4.764230807601693e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -2.934122129288896e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 2.736664154823537e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], -5.204406542692241e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 5.421010862427522e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 3.682366711745184e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 9.424226929452060e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.887504371245527e-02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 9.548686829671502e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 3.610485862939490e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.622842924371904e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 9.962005938584539e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -1.234418458497475e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -8.651933336434325e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 2.357315731195748e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[17], -4.309219854823383e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], 8.890457814381136e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 1.367616714833900e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[20], -6.823357878312430e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], -2.642126207483761e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], 1.428801915160532e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[23], 1.528891786927989e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -1.234831674791840e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -2.886611406995979e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -2.471545446220997e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 3.085214258020157e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.134014061579849e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -1.876680387611842e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 5.625718686789150e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -5.859078461986218e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[8], -9.871066699497807e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -6.791741558271588e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], 1.364776749417889e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[11], -5.787489328949392e-17, 1e-14) );
  }

  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.150017324508912e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -2.546528338451006e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -6.390355338505458e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 6.071532165918825e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -2.267473800546111e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -3.786915905595604e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -1.354168513434395e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 5.442909171135572e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[8], -7.840257220083711e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 2.043685901847408e+01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], 5.644556045385250e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[11], -2.132086621473017e-15, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.972046254947966e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 1.916788010228144e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -5.594736960135818e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -3.976819687358940e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -7.820544416664121e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], -5.314838274467147e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -1.266314256145179e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -1.100709950012486e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 3.153654747183855e-20, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.988818501979186e-02, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], -6.914046284241020e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 2.324342125145657e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.658940327502584e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 7.638194058889744e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -2.077400575207255e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -1.307818870560640e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 7.669292658022297e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[17], 2.054275056678685e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], -5.488773498207866e-18, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], -8.562274170981868e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[20], 3.564701589537615e-19, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], -5.186055989012229e-01, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], -2.942385704949726e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[23], 8.798763347754168e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 2.417204829760448e-03, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 3.161475558640719e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 4.965021343865173e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -7.135905991232719e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 2.263024721756329e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 3.806509604616858e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 1.296315688675056e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -5.436800278378912e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[8], 7.805672338358357e-17, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -6.872219553554434e+00, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], -4.857415646530594e-16, 1e-14) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[11], 5.814896626746330e-16, 1e-14) );
  }}

void test_1x3v_square_p1_4() { test_1x3v_square(1, 4, false); }
void test_1x3v_square_p1_8() { test_1x3v_square(1, 8, false); }
void test_1x3v_bump_p1_4() { test_1x3v_bump(1, 4, false); }
void test_1x3v_bump_p1_8() { test_1x3v_bump(1, 8, false); }

void test_1x3v_square_p2_4() { test_1x3v_square(2, 4, false); }
void test_1x3v_square_p2_8() { test_1x3v_square(2, 8, false); }
void test_1x3v_bump_p2_4() { test_1x3v_bump(2, 4, false); }
void test_1x3v_bump_p2_8() { test_1x3v_bump(2, 8, false); }

void test_1x3v_square_p1_4_cu() { test_1x3v_square(1, 4, true); }
void test_1x3v_square_p1_8_cu() { test_1x3v_square(1, 8, true); }
void test_1x3v_bump_p1_4_cu() { test_1x3v_bump(1, 4, true); }
void test_1x3v_bump_p1_8_cu() { test_1x3v_bump(1, 8, true); }

void test_1x3v_square_p2_4_cu() { test_1x3v_square(2, 4, true); }
void test_1x3v_square_p2_8_cu() { test_1x3v_square(2, 8, true); }
void test_1x3v_bump_p2_4_cu() { test_1x3v_bump(2, 4, true); }
void test_1x3v_bump_p2_8_cu() { test_1x3v_bump(2, 8, true); }

TEST_LIST = {
  { "test_1x3v_square_p1_4", test_1x3v_square_p1_4 },
  { "test_1x3v_square_p1_8", test_1x3v_square_p1_8 },
  { "test_1x3v_bump_p1_4", test_1x3v_bump_p1_4 },
  { "test_1x3v_bump_p1_8", test_1x3v_bump_p1_8 },
  { "test_1x3v_square_p2_4", test_1x3v_square_p2_4 },
  { "test_1x3v_square_p2_8", test_1x3v_square_p2_8 },
  { "test_1x3v_bump_p2_4", test_1x3v_bump_p2_4 },
  { "test_1x3v_bump_p2_8", test_1x3v_bump_p2_8 },

  #ifdef GKYL_HAVE_CUDA
  { "test_1x3v_square_p1_4_cu", test_1x3v_square_p1_4_cu },
  { "test_1x3v_square_p1_8_cu", test_1x3v_square_p1_8_cu },
  { "test_1x3v_bump_p1_4_cu", test_1x3v_bump_p1_4_cu },
  { "test_1x3v_bump_p1_8_cu", test_1x3v_bump_p1_8_cu },
  { "test_1x3v_square_p2_4_cu", test_1x3v_square_p2_4_cu },
  { "test_1x3v_square_p2_8_cu", test_1x3v_square_p2_8_cu },
  { "test_1x3v_bump_p2_4_cu", test_1x3v_bump_p2_4_cu },
  { "test_1x3v_bump_p2_8_cu", test_1x3v_bump_p2_8_cu },
  #endif
  { NULL, NULL }
};
