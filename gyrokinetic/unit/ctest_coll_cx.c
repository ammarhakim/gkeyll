#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_cx.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double d_ion_mass = GKYL_PROTON_MASS*2.01410177811;
double B0 = 0.5;
double check_fac = 1.0e10;

void eval_n(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_T_over_m_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/d_ion_mass;  
}
void eval_T_over_m_neut(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/d_ion_mass;  
}
void eval_bfield(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_jac(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

// test 2x2v / 2x3v
void
test_coll_cx_d(bool use_gpu)
{
  int charge_state = 0;
  // use vt = 40 eV for all grids
  double vmax_ion = 4.*sqrt(40.*echarge/d_ion_mass);
  double vmin_ion = -vmax_ion;
  double mumax_ion = 12.*40.*echarge/(2.*B0);
  int poly_order = 1;
  int cdim = 2, vdim_gk = 2, vdim_vl = 3;
  int pdim_gk = cdim + vdim_gk, pdim_vl = cdim + vdim_vl;

  // for gk grids 
  double lower_ion[] = {-2.0,-2.0,vmin_ion,0.0}, upper_ion[] = {2.0,2.0,vmax_ion,mumax_ion};
  int ghost_gk[] = {0, 0, 0, 0};
  int cells_gk[] = {2, 2, 16, 8};

  // for vlasov grid
  double lower_vl[] = {-2.0,-2.0,vmin_ion,vmin_ion,vmin_ion}, upper_vl[] = {2.0,2.0,vmax_ion,vmax_ion,vmax_ion};
  int ghost_vl[] = {0, 0, 0, 0, 0};
  int cells_vl[] = {2, 2, 16, 16, 16};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_ion, upper_ion, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange_ext, &confRange);

  // ion phase grid
  struct gkyl_rect_grid phaseGrid_ion;
  struct gkyl_range phaseRange_ion, phaseRange_ext_ion;
  gkyl_rect_grid_init(&phaseGrid_ion, pdim_gk, lower_ion, upper_ion, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_ion, ghost_gk, &phaseRange_ext_ion, &phaseRange_ion);

  // vlasov phase grid
  struct gkyl_rect_grid phaseGrid_vl;
  struct gkyl_range phaseRange_vl, phaseRange_ext_vl;
  gkyl_rect_grid_init(&phaseGrid_vl, pdim_vl, lower_vl, upper_vl, cells_vl);
  gkyl_create_grid_ranges(&phaseGrid_vl, ghost_vl, &phaseRange_ext_vl, &phaseRange_vl);

  // basis functions
  struct gkyl_basis phaseBasis_vl, phaseBasis_gk, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis_gk, cdim, vdim_gk);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // vlasov basis func
  gkyl_cart_modal_hybrid(&phaseBasis_vl, cdim, vdim_vl);

  // projection updater for moments
  gkyl_proj_on_basis *proj_n = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_n, NULL);
  gkyl_proj_on_basis *proj_T_over_m_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_ion, NULL);
  gkyl_proj_on_basis *proj_T_over_m_neut = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_neut, NULL);

  double vt_sq_ion_min = 1*echarge/d_ion_mass;
  double vt_sq_neut_min = 1*echarge/d_ion_mass;
  struct gkyl_dg_cx_inp cx_inp = {
    .grid = &phaseGrid_ion,
    .cbasis = &basis,
    .pbasis_gk = &phaseBasis_gk,
    .pbasis_vl = &phaseBasis_vl,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_ion,
    .mass_ion = d_ion_mass,
    .mass_neut = d_ion_mass,
    .vt_sq_ion_min = vt_sq_ion_min, 
    .vt_sq_neut_min = vt_sq_neut_min, 
    .type_ion = GKYL_ION_D,
  };
  
  // coll struct.
  struct gkyl_dg_cx *coll_cx_up = gkyl_dg_cx_new(&cx_inp, use_gpu);;

  struct gkyl_array *n_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_ion = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *n_neut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_neut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 5*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_cx = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // arrays necessary for prim_vars
  struct gkyl_array *u_par_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(proj_n, 0.0, &confRange, n_ion);
  gkyl_proj_on_basis_advance(proj_T_over_m_ion, 0.0, &confRange, T_over_m_ion);
  gkyl_array_set_offset(moms_ion, 1.0, n_ion, 0);
  gkyl_array_set_offset(moms_ion, 1.0, T_over_m_ion, 2*basis.num_basis);

  gkyl_proj_on_basis_advance(proj_n, 0.0, &confRange, n_neut);
  gkyl_proj_on_basis_advance(proj_T_over_m_neut, 0.0, &confRange, T_over_m_neut);
  gkyl_array_set_offset(moms_neut, 1.0, n_neut, 0);
  gkyl_array_set_offset(moms_neut, 1.0, T_over_m_neut, (1+vdim_vl)*basis.num_basis);

  // project b_i
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_neut_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, (2+vdim_vl)*basis.num_basis, confRange.volume);
    struct gkyl_array *moms_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *u_par_ion_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_cx_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_ion_cu, moms_ion);
    gkyl_array_copy(moms_neut_cu, moms_neut);
    gkyl_array_copy(u_par_ion_cu, u_par_ion);
    gkyl_array_copy(coef_cx_cu, coef_cx);
	
    gkyl_dg_cx_coll(coll_cx_up, moms_ion_cu, moms_neut_cu, 
      u_par_ion_cu, coef_cx_cu, 0);
    gkyl_array_copy(coef_cx, coef_cx_cu);

    gkyl_array_release(moms_neut_cu); 
    gkyl_array_release(moms_ion_cu);
    gkyl_array_release(u_par_ion_cu); 
    gkyl_array_release(coef_cx_cu);
  }
  else {
    gkyl_dg_cx_coll(coll_cx_up, moms_ion, moms_neut, 
      u_par_ion, coef_cx, 0);    
  }
  const double *cv_cx = gkyl_array_cfetch(coef_cx, gkyl_range_idx(&confRange, (int[2]) { 1, 1}));

  // test against predicted value
  double p1_vals[] = {3.242709205939892e-14, 0.000000000000000e+00,
		      0.000000000000000e+00, 0.000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_cx[i]*check_fac, 1e-12) );
  }

  gkyl_array_release(n_ion); 
  gkyl_array_release(T_over_m_ion); 
  gkyl_array_release(n_neut); 
  gkyl_array_release(T_over_m_neut);
  gkyl_array_release(moms_neut); 
  gkyl_array_release(moms_ion);
  gkyl_array_release(u_par_ion);
  gkyl_array_release(b_i); 
  gkyl_array_release(b_x); 
  gkyl_array_release(b_y); 
  gkyl_array_release(b_z);
  gkyl_array_release(coef_cx);
  gkyl_proj_on_basis_release(proj_n); 
  gkyl_proj_on_basis_release(proj_T_over_m_ion); 
  gkyl_proj_on_basis_release(proj_T_over_m_neut); 
  gkyl_dg_cx_release(coll_cx_up); 
}


void coll_cx_d() { test_coll_cx_d(false); }

#ifdef GKYL_HAVE_CUDA
void coll_cx_d_gpu() { test_coll_cx_d(true); }
#endif

TEST_LIST = {
  { "coll_cx_d", coll_cx_d },
#ifdef GKYL_HAVE_CUDA
  { "coll_cx_d_gpu", coll_cx_d_gpu },
#endif
  { NULL, NULL },
};
