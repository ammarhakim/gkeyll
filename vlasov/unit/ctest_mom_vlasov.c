// Test calculation of Vlasov moments of a distribution function.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

void
test_mom_vlasov()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis, vbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&vbasis, 3, poly_order); // 3V
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

  struct gkyl_array hamil;
  struct gkyl_range velLocal;

  struct gkyl_mom_vlasov_inp inp_mom_m2ij = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = &hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M2IJ, 
    .use_gpu = false,
  };

 struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_inew(&inp_mom_m2ij);

  TEST_CHECK( m2ij->cdim == 1 );
  TEST_CHECK( m2ij->pdim == 4 );
  TEST_CHECK( m2ij->poly_order == 2 );
  TEST_CHECK( m2ij->num_config == cbasis.num_basis );
  TEST_CHECK( m2ij->num_phase == pbasis.num_basis );
  TEST_CHECK( m2ij->num_mom == 6 );

  struct gkyl_mom_vlasov_inp inp_mom_m3ijk = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = &hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M3IJK, 
    .use_gpu = false,
  };

 struct gkyl_mom_type *m3ijk = gkyl_mom_vlasov_inew(&inp_mom_m3ijk);
  TEST_CHECK( m3ijk->num_mom == 10 );
 
}

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = (x*x)*(vx-0.5)*(vx-0.5);
}

void evalFuncSimple(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = 1./(2.0-(-2.0));
}

void evalFunc4d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double vx = xn[2], vy = xn[3];
  fout[0] = (vx-1)*(vy-2);
}

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

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void
test_1x1v_p1()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {4, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M0M1M2;
  struct gkyl_mom_type *int_t = gkyl_int_mom_vlasov_inew(&inp_mom);
  
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  gkyl_mom_calc *intcalc = gkyl_mom_calc_new(&grid, int_t, false);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *int_mom;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  int_mom = mkarr(vdim+2, confLocal_ext.volume);
  
  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
  gkyl_mom_calc_advance(intcalc, &local, &confLocal, distf, int_mom);
 
  // Check M0.
  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  2.089893375506907e+01, m00[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -7.756717518813397e+00, m00[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009867e+00, m01[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -2.585572506271132e+00, m01[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009867e+00, m02[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.585572506271132e+00, m02[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.089893375506907e+01, m03[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  7.756717518813396e+00, m03[1], 1e-12) );

  // Check M1i.
  double *m1i0 = gkyl_array_fetch(m1i, 0+confGhost[0]); double *m1i1 = gkyl_array_fetch(m1i, 1+confGhost[0]);
  double *m1i2 = gkyl_array_fetch(m1i, 2+confGhost[0]); double *m1i3 = gkyl_array_fetch(m1i, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare( -1.319932658214889e+01, m1i0[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  4.898979485566356e+00, m1i0[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.885618083164127e+00, m1i1[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.632993161855452e+00, m1i1[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.885618083164127e+00, m1i2[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.632993161855452e+00, m1i2[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.319932658214889e+01, m1i3[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -4.898979485566354e+00, m1i3[1], 1e-12) );

  // Check M2.
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  2.969848480983499e+01, m20[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.102270384252430e+01, m20[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  4.242640687119285e+00, m21[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -3.674234614174766e+00, m21[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  4.242640687119285e+00, m22[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  3.674234614174766e+00, m22[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.969848480983499e+01, m23[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.102270384252430e+01, m23[1], 1e-12) );

  // sum the integrated moments
  double red_mom[vdim+2];
  gkyl_array_reduce_range(red_mom, int_mom, GKYL_SUM, &confLocal);

  TEST_CHECK( gkyl_compare(  33.7777777777777786, red_mom[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -21.3333333333333357, red_mom[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  48.0000000000000000, red_mom[2], 1e-12) );

  //// Write the moment array to file.
  //const char *fmt = "%s-%s-%s_%d.gkyl";
  //char name[] = "mom_calc";
  //char speciesName[] = "distf";
  //char momName[] = gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0];
  //int frame = 0;
  //int sz = snprintf(0, 0, fmt, name, speciesName, momName, frame);
  //char fileNm[sz+1]; // ensures no buffer overflow
  //snprintf(fileNm, sizeof fileNm, fmt, name, speciesName, momName, frame);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, "ctest_mom_vlasov_distf_0.gkyl");

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(int_mom);
  
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_calc_release(intcalc);
  
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);
  gkyl_mom_type_release(int_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x1v_p2()
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {4, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M0M1M2;
  struct gkyl_mom_type *int_t = gkyl_int_mom_vlasov_inew(&inp_mom);
  
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  gkyl_mom_calc *intcalc = gkyl_mom_calc_new(&grid, int_t, false);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *int_mom;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  int_mom = mkarr(vdim+2, confLocal_ext.volume);
  
  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
  gkyl_mom_calc_advance(intcalc, &local, &confLocal, distf, int_mom);
 
  // Check M0.
  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  2.089893375506907e+01, m00[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -7.756717518813399e+00, m00[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.675919504799956e-01, m00[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009868e+00, m01[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -2.585572506271133e+00, m01[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.675919504799922e-01, m01[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009868e+00, m02[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.585572506271133e+00, m02[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.675919504799922e-01, m02[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.089893375506907e+01, m03[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  7.756717518813399e+00, m03[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.675919504799960e-01, m03[2], 1e-12) );

  // Check M1i.
  double *m1i0 = gkyl_array_fetch(m1i, 0+confGhost[0]); double *m1i1 = gkyl_array_fetch(m1i, 1+confGhost[0]);
  double *m1i2 = gkyl_array_fetch(m1i, 2+confGhost[0]); double *m1i3 = gkyl_array_fetch(m1i, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  -1.759910210953185e+01, m1i0[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   6.531972647421806e+00, m1i0[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -5.621826951410490e-01, m1i0[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.514157444218835e+00, m1i1[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   2.177324215807270e+00, m1i1[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -5.621826951410490e-01, m1i1[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.514157444218835e+00, m1i2[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.177324215807270e+00, m1i2[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -5.621826951410460e-01, m1i2[2], 1e-12) );
  TEST_CHECK( gkyl_compare(  -1.759910210953185e+01, m1i3[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -6.531972647421807e+00, m1i3[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -5.621826951410495e-01, m1i3[2], 1e-12) );

  // Check M2.
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(   2.331881029512970e+01, m20[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -8.654863757833898e+00, m20[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   7.261526478905216e-01, m20[2], 1e-12) );
  TEST_CHECK( gkyl_compare(   3.331258613589957e+00, m21[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.884954585944633e+00, m21[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   7.261526478905216e-01, m21[2], 1e-12) );
  TEST_CHECK( gkyl_compare(   3.331258613589958e+00, m22[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   2.884954585944633e+00, m22[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   7.261526478905216e-01, m22[2], 1e-12) );
  TEST_CHECK( gkyl_compare(   2.331881029512970e+01, m23[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   8.654863757833896e+00, m23[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   7.261526478905216e-01, m23[2], 1e-12) );

  // sum the integrated moments
  double red_mom[vdim+2];
  gkyl_array_reduce_range(red_mom, int_mom, GKYL_SUM, &confLocal);

  TEST_CHECK( gkyl_compare(  33.7777777777778, red_mom[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -28.4444444444444, red_mom[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  37.6888888888888971, red_mom[2], 1e-12) );

  //// Write the moment array to file.
  //const char *fmt = "%s-%s-%s_%d.gkyl";
  //char name[] = "mom_calc";
  //char speciesName[] = "distf";
  //char momName[] = gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0];
  //int frame = 0;
  //int sz = snprintf(0, 0, fmt, name, speciesName, momName, frame);
  //char fileNm[sz+1]; // ensures no buffer overflow
  //snprintf(fileNm, sizeof fileNm, fmt, name, speciesName, momName, frame);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, "ctest_mom_vlasov_distf_0.gkyl");

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(int_mom);
  
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_calc_release(intcalc);
  
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);
  gkyl_mom_type_release(int_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x2v_p1()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {4, 64, 64};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  
  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
 
  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m0ptr = gkyl_array_fetch(m0, linc);
    TEST_CHECK( gkyl_compare( 5.656854249492378, m0ptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m0ptr[k]);
  }}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m1iptr = gkyl_array_fetch(m1i, linc);
    for (unsigned int k=0; k<m1i->ncomp; ++k) {
      TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m1iptr[k]);
  }}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m2ptr = gkyl_array_fetch(m2, linc);
    TEST_CHECK( gkyl_compare( 7.5461551804752e+00, m2ptr[0], 1e-12) );
    TEST_MSG("Expected: %.13e in cell (%d)", 7.5461551804752e+00, i);
    TEST_MSG("Produced: %.13e", m2ptr[0]);
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m2ptr[k]);
  }}

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x2v_p2()
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {4, 64, 64};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0 , 0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  
  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
 
  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m0ptr = gkyl_array_fetch(m0, linc);
    TEST_CHECK( gkyl_compare( 5.656854249492378, m0ptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m0ptr[k]);
  }}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m1iptr = gkyl_array_fetch(m1i, linc);
    for (unsigned int k=0; k<m1i->ncomp; ++k) {
      TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m1iptr[k]);
  }}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m2ptr = gkyl_array_fetch(m2, linc);
    TEST_CHECK( gkyl_compare( 7.5424723326565e+00, m2ptr[0], 1e-12) );
    TEST_MSG("Expected: %.13e in cell (%d)", 7.5424723326565e+00, i);
    TEST_MSG("Produced: %.13e", m2ptr[0]);
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
      TEST_MSG("Produced: %.13e", m2ptr[k]);
  }}

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}


void
test_2x2v_p1()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {4, 4, 64, 64};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 2;
  
  double velLower[] = {lower[2], lower[3]}, velUpper[] = {upper[2], upper[3]};
  int velCells[] = {cells[2], cells[3]};
  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  
  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
 
  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 8.0, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m0ptr[k]);
  }}}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);
      for (unsigned int k=0; k<m1i->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m1iptr[k]);
  }}}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 1.0671875000000e+01, m2ptr[0], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 1.0671875000000e+01, i);
      TEST_MSG("Produced: %.13e", m2ptr[0]);
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m2ptr[k]);
  }}}

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_big_2x2v_p2()
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {22, 22, 32, 32};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 2;

  double velLower[] = {lower[2], lower[3]}, velUpper[] = {upper[2], upper[3]};
  int velCells[] = {cells[2], cells[3]};
  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc4d, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  
  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
 
  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 64.0, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m0ptr[k]);
  }}}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);

      TEST_CHECK( gkyl_compare( -85.33333333333333333, m1iptr[0], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", -85.33333333333333333, i);
      TEST_MSG("Produced: %.13e", m1iptr[0]);

      TEST_CHECK( gkyl_compare( -42.6666666666666667, m1iptr[8], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", -42.6666666666666667, i);
      TEST_MSG("Produced: %.13e", m1iptr[8]);
    }}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 8.5333333333333e+01, m2ptr[0], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 8.5333333333333e+01, i);
      TEST_MSG("Produced: %.13e", m2ptr[0]);
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m2ptr[k]);
  }}}

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_2x3v_p1()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0, 2.0};
  int cells[] = {4, 4, 64, 64, 64};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 3, cdim = 2;

  double velLower[] = {lower[2], lower[3], lower[4]}, velUpper[] = {upper[2], upper[3], upper[4]};
  int velCells[] = {cells[2], cells[3], cells[4]};
  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0, 0, 0 };
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &velBasis, &velLocal, 
    GKYL_MODEL_DEFAULT, 0, hamil, gamma_inv, false); 

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_range = &velLocal,
    .use_vmap = 0, 
    .vmap = 0, 
    .hamil_range = &velLocal,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .mom_type = GKYL_F_MOMENT_M0, 
    .use_gpu = false,
  };

  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M1;
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_inew(&inp_mom);
  inp_mom.mom_type = GKYL_F_MOMENT_M2;
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t, false);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t, false);
  
  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, &local, &confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);
 
  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 31.999999999999996, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m0ptr[k]);
  }}}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);
      for (unsigned int k=0; k<m1i->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m1iptr[k]);
  }}}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 6.4031250000000e+01, m2ptr[0], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", 6.4031250000000e+01, i);
      TEST_MSG("Produced: %.13e", m2ptr[0]);
      for (unsigned int k=1; k<confBasis.num_basis; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
        TEST_MSG("Expected: %.13e in cell (%d)", 0., i);
        TEST_MSG("Produced: %.13e", m2ptr[k]);
  }}}

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

#ifdef GKYL_HAVE_CUDA
int cu_mom_vlasov_test(const struct gkyl_mom_type *mom);

void
test_cu_mom_vlasov()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M2IJ, true);

  // Removing these lines because of the issues with linking *.cu files at the /vlasov layer and above.
  // -- JG (05/13/2025)
  //int nfail = cu_mom_vlasov_test(m2ij->on_dev);
  //TEST_CHECK( nfail == 0 );

  gkyl_cu_free(m2ij);
}

void
test_1x1v_p1_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {4, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create distribution function and moment arrays
  struct gkyl_array *distf, *m0, *m1i, *m2;
  distf = mkarr(basis.num_basis, local_ext.volume);
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // initialize device arrays
  struct gkyl_array *distf_cu, *m0_cu, *m1i_cu, *m2_cu;
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_copy(distf_cu, distf);

  // Create vlasov moment and mom_calc objects on device.
  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0, true);
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M1, true);
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M2, true);
  gkyl_mom_calc *m0Calc = gkyl_mom_calc_new(&grid, vmM0_t, true);
  gkyl_mom_calc *m1iCalc = gkyl_mom_calc_new(&grid, vmM1i_t, true);
  gkyl_mom_calc *m2Calc = gkyl_mom_calc_new(&grid, vmM2_t, true);

  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0Calc, &local, &confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1iCalc, &local, &confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2Calc, &local, &confLocal, distf_cu, m2_cu);
  gkyl_array_copy(m0, m0_cu);
  gkyl_array_copy(m1i, m1i_cu);
  gkyl_array_copy(m2, m2_cu);

  // Check M0.
  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare( 20.898933755069070, m00[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -7.756717518813397, m00[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009867, m01[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -2.585572506271132, m01[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.985561965009867, m02[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  2.585572506271132, m02[1], 1e-12) );
  TEST_CHECK( gkyl_compare( 20.898933755069070, m03[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  7.756717518813396, m03[1], 1e-12) );

  // Check M1i.
  double *m1i0 = gkyl_array_fetch(m1i, 0+confGhost[0]); double *m1i1 = gkyl_array_fetch(m1i, 1+confGhost[0]);
  double *m1i2 = gkyl_array_fetch(m1i, 2+confGhost[0]); double *m1i3 = gkyl_array_fetch(m1i, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare( -17.599102109531849, m1i0[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   6.531972647421810, m1i0[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.514157444218835, m1i1[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   2.177324215807269, m1i1[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.514157444218835, m1i2[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -2.177324215807269, m1i2[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -17.599102109531849, m1i3[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -6.531972647421808, m1i3[1], 1e-12) );

  // Check M2.
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  45.464347116290611, m20[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -16.874262672506337, m20[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   6.494906730898657, m21[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -5.624754224168777, m21[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   6.494906730898657, m22[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   5.624754224168777, m22[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  45.464347116290611, m23[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  16.874262672506337, m23[1], 1e-12) );

  // free allocated memory.
  gkyl_array_release(distf); gkyl_array_release(distf_cu);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);

  gkyl_mom_calc_release(m0Calc); gkyl_mom_calc_release(m1iCalc); gkyl_mom_calc_release(m2Calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);  

  gkyl_proj_on_basis_release(projDistf);
}

void
test_1x2v_p1_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {4, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis and copy to device
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);

  // Create vlasov moment and mom_calc objects on device.
  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0, true);
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M1, true);
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M2, true);
  gkyl_mom_calc *m0Calc = gkyl_mom_calc_new(&grid, vmM0_t, true);
  gkyl_mom_calc *m1iCalc = gkyl_mom_calc_new(&grid, vmM1i_t, true);
  gkyl_mom_calc *m2Calc = gkyl_mom_calc_new(&grid, vmM2_t, true);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *m0_cu, *m1i_cu, *m2_cu;;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0Calc, &local, &confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1iCalc, &local, &confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2Calc, &local, &confLocal, distf_cu, m2_cu);
  gkyl_array_copy(m0, m0_cu);
  gkyl_array_copy(m1i, m1i_cu);
  gkyl_array_copy(m2, m2_cu);

  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m0ptr = gkyl_array_fetch(m0, linc);
    TEST_CHECK( gkyl_compare( 5.656854249492378, m0ptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
  }}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m1iptr = gkyl_array_fetch(m1i, linc);
    for (unsigned int k=0; k<m1i->ncomp; ++k) {
      TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
  }}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[2] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *m2ptr = gkyl_array_fetch(m2, linc);
    TEST_CHECK( gkyl_compare( 15.084944665313010, m2ptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
  }}

  // free allocated memory.
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);

  gkyl_mom_calc_release(m0Calc); gkyl_mom_calc_release(m1iCalc); gkyl_mom_calc_release(m2Calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);  

  gkyl_proj_on_basis_release(projDistf);
}

void
test_2x2v_p1_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {4, 4, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis and copy to device
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);

  // Create vlasov moment and mom_calc objects on device.
  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0, true);
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M1, true);
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M2, true);
  gkyl_mom_calc *m0Calc = gkyl_mom_calc_new(&grid, vmM0_t, true);
  gkyl_mom_calc *m1iCalc = gkyl_mom_calc_new(&grid, vmM1i_t, true);
  gkyl_mom_calc *m2Calc = gkyl_mom_calc_new(&grid, vmM2_t, true);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *m0_cu, *m1i_cu, *m2_cu;;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0Calc, &local, &confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1iCalc, &local, &confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2Calc, &local, &confLocal, distf_cu, m2_cu);
  gkyl_array_copy(m0, m0_cu);
  gkyl_array_copy(m1i, m1i_cu);
  gkyl_array_copy(m2, m2_cu);

  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 8.0, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<m0->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
  }}}

  //// Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);
      for (unsigned int k=0; k<m1i->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
  }}}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 21.333333333333336, m2ptr[0], 1e-12) );
      for (unsigned int k=1; k<m2->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
  }}}

  // free allocated memory.
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);

  gkyl_mom_calc_release(m0Calc); gkyl_mom_calc_release(m1iCalc); gkyl_mom_calc_release(m2Calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);  

  gkyl_proj_on_basis_release(projDistf);
}

void
test_2x3v_p1_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0, 2.0};
  int cells[] = {4, 4, 2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 3, cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFuncSimple, NULL);

  // create distribution function
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis and copy to device
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);

  // Create vlasov moment and mom_calc objects on device.
  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0, true);
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M1, true);
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M2, true);
  gkyl_mom_calc *m0Calc = gkyl_mom_calc_new(&grid, vmM0_t, true);
  gkyl_mom_calc *m1iCalc = gkyl_mom_calc_new(&grid, vmM1i_t, true);
  gkyl_mom_calc *m2Calc = gkyl_mom_calc_new(&grid, vmM2_t, true);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *m0_cu, *m1i_cu, *m2_cu;;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0Calc, &local, &confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1iCalc, &local, &confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2Calc, &local, &confLocal, distf_cu, m2_cu);
  gkyl_array_copy(m0, m0_cu);
  gkyl_array_copy(m1i, m1i_cu);
  gkyl_array_copy(m2, m2_cu);

  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 31.999999999999996, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<m0->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
  }}}

  // Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);
      for (unsigned int k=0; k<m1i->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m1iptr[k], 1e-12) );
  }}}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 127.999999999999972, m2ptr[0], 1e-12) );
      for (unsigned int k=1; k<m2->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
  }}}

  // free allocated memory.
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);

  gkyl_mom_calc_release(m0Calc); gkyl_mom_calc_release(m1iCalc); gkyl_mom_calc_release(m2Calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_proj_on_basis_release(projDistf);
}

void
test_big_2x2v_p2_cu()
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {44, 44, 32, 32};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc4d, NULL);

  // create distribution function
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis and copy to device
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);

  // Create vlasov moment and mom_calc objects on device.
  struct gkyl_mom_type *vmM0_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0, true);
  struct gkyl_mom_type *vmM1i_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M1, true);
  struct gkyl_mom_type *vmM2_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M2, true);
  gkyl_mom_calc *m0Calc = gkyl_mom_calc_new(&grid, vmM0_t, true);
  gkyl_mom_calc *m1iCalc = gkyl_mom_calc_new(&grid, vmM1i_t, true);
  gkyl_mom_calc *m2Calc = gkyl_mom_calc_new(&grid, vmM2_t, true);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *m0_cu, *m1i_cu, *m2_cu;;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0Calc, &local, &confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1iCalc, &local, &confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2Calc, &local, &confLocal, distf_cu, m2_cu);
  gkyl_array_copy(m0, m0_cu);
  gkyl_array_copy(m1i, m1i_cu);
  gkyl_array_copy(m2, m2_cu);

  // Check M0.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m0ptr = gkyl_array_fetch(m0, linc);
      TEST_CHECK( gkyl_compare( 64.0, m0ptr[0], 1e-12) );
      for (unsigned int k=1; k<m0->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m0ptr[k], 1e-12) );
  }}}

  //// Check M1i.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m1iptr = gkyl_array_fetch(m1i, linc);

      TEST_CHECK( gkyl_compare( -85.33333333333333333, m1iptr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( -42.6666666666666667, m1iptr[8], 1e-12) );
  }}

  // Check M2.
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    for (unsigned int j=1; j<cells[1]+1; ++j) {
      int cidx[2] = {i,j};
      long linc = gkyl_range_idx(&confLocal, cidx);
      double *m2ptr = gkyl_array_fetch(m2, linc);
      TEST_CHECK( gkyl_compare( 170.666666666666667, m2ptr[0], 1e-12) );
      for (unsigned int k=1; k<m2->ncomp; ++k) {
        TEST_CHECK( gkyl_compare( 0., m2ptr[k], 1e-12) );
  }}}

  // free allocated memory.
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);

  gkyl_mom_calc_release(m0Calc); gkyl_mom_calc_release(m1iCalc); gkyl_mom_calc_release(m2Calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);  

  gkyl_proj_on_basis_release(projDistf);
}

#endif

TEST_LIST = {
  { "mom_vlasov", test_mom_vlasov },
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x2v_p1", test_1x2v_p1 },
  { "test_1x2v_p2", test_1x2v_p2 },
  { "test_2x2v_p1", test_2x2v_p1 },
  { "test_big_2x2v_p2", test_big_2x2v_p2 },  
  { "test_2x3v_p1", test_2x3v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "cu_mom_vlasov", test_cu_mom_vlasov },
  { "test_1x1v_p1_cu", test_1x1v_p1_cu },
  { "test_1x2v_p1_cu", test_1x2v_p1_cu },
  { "test_2x2v_p1_cu", test_2x2v_p1_cu },
  { "test_2x3v_p1_cu", test_2x3v_p1_cu },
//  { "test_big_2x2v_p2_cu", test_big_2x2v_p2_cu },
#endif
  { NULL, NULL },
};
