// Tests for the Vlasov FPO Maxwellian Rosenbluth potential projection and
// drag and diffusion coefficient computation

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
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_gamma(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->gamma0;
}

void eval_lte_moms(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->n0;
  fout[1] = app->ux0;
  fout[2] = app->uy0;
  fout[3] = app->uz0;
  fout[4] = app->vth0*app->vth0;
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

void test_1x3v_p2(int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {1, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 5.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  gkyl_proj_on_basis *proj_gamma = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 1, eval_gamma, &ctx);
  gkyl_proj_on_basis *proj_lte = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 5, eval_lte_moms, &ctx);

  struct gkyl_array *lte_moms, *gamma, *h, *g, *h_surf, *g_surf, *dhdv_surf, *dgdv_surf, *d2gdv2_surf;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  lte_moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  gamma = mkarr(conf_basis.num_basis, conf_range_ext.volume);
  h = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  g = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  h_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  g_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dhdv_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dgdv_surf = mkarr(2*vdim*surf_basis.num_basis, phase_range_ext.volume);
  d2gdv2_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);

  // Array of relative offsets for recovery stencils
  long offsets[36];

  // Initialize updaters for potentials and coeffs
  gkyl_proj_maxwellian_pots_on_basis *pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
    &phase_grid, &conf_range, &phase_range, &conf_basis, &phase_basis, poly_order+1, 0);

  gkyl_fpo_vlasov_coeff_recovery *coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(&phase_grid,
    &phase_basis, &phase_range_ext, offsets, 0);

  // Project moments and compute potentials
  gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma);
  gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms);
  gkyl_proj_maxwellian_pots_on_basis_advance(pot_slvr, &phase_range, &conf_range, 
    lte_moms, h, g, h_surf, g_surf,
    dhdv_surf, dgdv_surf, d2gdv2_surf);

  // Compute drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, &phase_range,
    &conf_range, gamma, h, dhdv_surf, drag_coeff, drag_coeff_surf, 0);

  gkyl_calc_fpo_diff_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, 
    &phase_range, &conf_range, gamma,
    g, g_surf, dgdv_surf, d2gdv2_surf, 
    diff_coeff, diff_coeff_surf, 0);

  const char *fmt_drag = "ctest_fpo_drag_coeff_%d.gkyl";
  const char *fmt_diff = "ctest_fpo_diff_coeff_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt_drag, NV);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof(fileNm), fmt_drag, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, drag_coeff, fileNm);
  snprintf(fileNm, sizeof(fileNm), fmt_diff, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, diff_coeff, fileNm);

  // Compare values in corner cell and interior cell (of velocity space)
  // Checking first five components of a_x, a_y, D_xx, D_xy, and D_yx
  int idx_corner[] = {1, 1, 1, 1};
  int idx_int[] = {1, 2, 2, 2};
  long lin_corner = gkyl_range_idx(&phase_range, idx_corner);
  long lin_int = gkyl_range_idx(&phase_range, idx_int);

  double *drag_corner, *drag_int;
  double *diff_corner, *diff_int;
  drag_corner = gkyl_array_fetch(drag_coeff, lin_corner);
  drag_int = gkyl_array_fetch(drag_coeff, lin_int);
  diff_corner = gkyl_array_fetch(diff_coeff, lin_corner);
  diff_int = gkyl_array_fetch(diff_coeff, lin_int);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.1095780224753596, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.0006267138106831, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 0.1095780224753595, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 0.0006267138106833, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 0.5997033462078907, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 0.0000000000000003, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -0.1738666273128693, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 0.1892738907386429, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 0.1892738907386429, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.5997033462078906, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.1892738907386430, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.1738666273128693, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.1892738907386430, 1e-12) );
    
    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.4100999889255052, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.0000000000000009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.0752757190222467, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.0027322422345083, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.0027322422345056, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.1883902242384373, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -0.0000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.0015462754276331, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.0069468342856277, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.0334615672390714, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.1883902242384379, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.0069468342856311, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.0015462754276336, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.0334615672390713, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450405, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.0000000000000006, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.3227918552176465, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.0770701235291489, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.0770701235291491, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.1888164518268585, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.0735188324952567, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.0735188324952577, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.0414524066212169, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.1888164518268585, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.0735188324952561, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.0735188324952565, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.0414524066212169, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.0804379082040003, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.0000347913213152, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 0.0804379082040009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 0.0000347913213154, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 0.1576774329802843, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 0.0001724567096385, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.1576774329802842, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 0.0001724567096392, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.0182634189498116, 1e-12) );
    
    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.3518364524562119, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.0000000000000015, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.0280494232833636, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.0005618962338415, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.0005618962338392, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.1673341378406684, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.0000000000000004, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.0005142107433202, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.0025349166383694, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.0132056938305365, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.1673341378406678, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], -0.0000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.0025349166383768, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.0005142107433286, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.0132056938305363, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.4927165215706708, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.0000000000000075, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.0527284542989293, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.0022095876646917, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.0022095876646904, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.2178235155393529, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.0019949327795594, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.0019949327795655, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.0230457152585550, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.2178235155393529, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.0019949327795587, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.0019949327795674, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.0230457152585550, 1e-12) );
  }

  // Release memory
  gkyl_array_release(lte_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(h);
  gkyl_array_release(g);
  gkyl_array_release(h_surf);
  gkyl_array_release(g_surf);
  gkyl_array_release(dhdv_surf);
  gkyl_array_release(dgdv_surf);
  gkyl_array_release(d2gdv2_surf);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);

  gkyl_proj_maxwellian_pots_on_basis_release(pot_slvr);
  gkyl_fpo_vlasov_coeff_recovery_release(coeff_recovery);
}

void test_1x3v_p2_cu(int NV)
{
  bool use_gpu = true;
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {1, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 5.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  gkyl_proj_on_basis *proj_gamma = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 1, eval_gamma, &ctx);
  gkyl_proj_on_basis *proj_lte = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 5, eval_lte_moms, &ctx);

  struct gkyl_array *lte_moms, *gamma, *h, *g, *h_surf, *g_surf, *dhdv_surf, *dgdv_surf, *d2gdv2_surf;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  lte_moms = mkarr_cu(5*conf_basis.num_basis, conf_range_ext.volume);
  gamma = mkarr_cu(conf_basis.num_basis, conf_range_ext.volume);
  h = mkarr_cu(phase_basis.num_basis, phase_range_ext.volume);
  g = mkarr_cu(phase_basis.num_basis, phase_range_ext.volume);
  h_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  g_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dhdv_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dgdv_surf = mkarr_cu(2*vdim*surf_basis.num_basis, phase_range_ext.volume);
  d2gdv2_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  drag_coeff = mkarr_cu(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr_cu(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr_cu(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);

  // Array of relative offsets for recovery stencils
  long offsets[36];

  // Initialize updaters for potentials and coeffs
  gkyl_proj_maxwellian_pots_on_basis *pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
    &phase_grid, &conf_range, &phase_range, &conf_basis, &phase_basis, poly_order+1, use_gpu);

  gkyl_fpo_vlasov_coeff_recovery *coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(&phase_grid,
    &phase_basis, &phase_range_ext, offsets, use_gpu);

  // Project moments and compute potentials
  gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma);
  gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms);
  gkyl_proj_maxwellian_pots_on_basis_advance(pot_slvr, &phase_range, &conf_range, 
    lte_moms, h, g, h_surf, g_surf,
    dhdv_surf, dgdv_surf, d2gdv2_surf);

  // Compute drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, &phase_range,
    &conf_range, gamma, h, dhdv_surf, drag_coeff, drag_coeff_surf, use_gpu);

  gkyl_calc_fpo_diff_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, 
    &phase_range, &conf_range, gamma,
    g, g_surf, dgdv_surf, d2gdv2_surf, 
    diff_coeff, diff_coeff_surf, use_gpu);

  // Compare values in corner cell and interior cell (of velocity space)
  // Checking first five components of a_x, a_y, D_xx, D_xy, and D_yx
  int idx_corner[] = {1, 1, 1, 1};
  int idx_int[] = {1, 2, 2, 2};
  long lin_corner = gkyl_range_idx(&phase_range, idx_corner);
  long lin_int = gkyl_range_idx(&phase_range, idx_int);

  double *drag_corner, *drag_int;
  double *diff_corner, *diff_int;
  drag_corner = gkyl_array_fetch(drag_coeff, lin_corner);
  drag_int = gkyl_array_fetch(drag_coeff, lin_int);
  diff_corner = gkyl_array_fetch(diff_coeff, lin_corner);
  diff_int = gkyl_array_fetch(diff_coeff, lin_int);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.1095780224753596, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.0006267138106831, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 0.1095780224753595, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 0.0006267138106833, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 0.0212899929067877, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 0.5997033462078907, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 0.0000000000000003, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -0.1738666273128693, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 0.1892738907386429, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 0.1892738907386429, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.5997033462078906, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.1892738907386430, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.1738666273128693, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.1892738907386430, 1e-12) );
    
    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.4100999889255052, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.0000000000000009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.0752757190222467, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.0027322422345083, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.0027322422345056, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.1883902242384373, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -0.0000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.0015462754276331, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.0069468342856277, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.0334615672390714, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.1883902242384379, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.0069468342856311, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.0015462754276336, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.0334615672390713, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450405, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.0000000000000006, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.3227918552176465, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.0770701235291489, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.0770701235291491, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.1888164518268585, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.0735188324952567, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.0735188324952577, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.0414524066212169, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.1888164518268585, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.0735188324952561, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.0735188324952565, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.0414524066212169, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.0804379082040003, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.0000347913213152, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 0.0804379082040009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 0.0000347913213154, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 0.0066446451201719, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 0.1576774329802843, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 0.0001724567096385, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.1576774329802842, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.0000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.0182634189498116, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 0.0001724567096392, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.0182634189498116, 1e-12) );
    
    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.3518364524562119, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.0000000000000015, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.0280494232833636, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.0005618962338415, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.0005618962338392, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.1673341378406684, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.0000000000000004, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.0005142107433202, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.0025349166383694, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.0132056938305365, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.1673341378406678, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], -0.0000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.0025349166383768, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.0005142107433286, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.0132056938305363, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.4927165215706708, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.0000000000000075, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.0527284542989293, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.0022095876646917, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.0022095876646904, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.2178235155393529, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.0019949327795594, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.0019949327795655, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.0230457152585550, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.2178235155393529, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], 0.0000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.0019949327795587, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.0019949327795674, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.0230457152585550, 1e-12) );
  }

  // Release memory
  gkyl_array_release(lte_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(h);
  gkyl_array_release(g);
  gkyl_array_release(h_surf);
  gkyl_array_release(g_surf);
  gkyl_array_release(dhdv_surf);
  gkyl_array_release(dgdv_surf);
  gkyl_array_release(d2gdv2_surf);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);

  gkyl_proj_maxwellian_pots_on_basis_release(pot_slvr);
  gkyl_fpo_vlasov_coeff_recovery_release(coeff_recovery);
}

void test_1x3v_p2_4() { test_1x3v_p2(4); }
void test_1x3v_p2_8() { test_1x3v_p2(8); }
void test_1x3v_p2_4_cu() { test_1x3v_p2_cu(4); }
void test_1x3v_p2_8_cu() { test_1x3v_p2_cu(8); }

TEST_LIST = {
  { "test_1x3v_p2_4", test_1x3v_p2_4 },
  { "test_1x3v_p2_8", test_1x3v_p2_8 },
  #ifdef GKYL_HAVE_CUDA
  { "test_1x3v_p2_4_cu", test_1x3v_p2_4_cu }
  { "test_1x3v_p2_8_cu", test_1x3v_p2_8_cu }
  #endif
  { NULL, NULL }
};
