#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_triad_math_1v()
{
  // Test: Can the correct bracket be constructed for:
  // 1D streched coordinates: x(z^1) = (z^1)^2 \sigma_x
  int vdim = 1;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[1];

  // Randomly choosen coordinates
  double z1 = 3.2;

  // Assign covaraint tangent basis 
  vierbein[0] = 2*z1;

  // 1. Test kernel_vierbein_inv_1v()
  double vierbein_inv[1];
  kernel_vierbein_inv_1v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0 / (2.0 * z1), eps));

  // 2. Test kernel_conf_poisson_tensor_1v()
  double conf_poisson_tensor[1];

  // Assign the triad basis gradient
  double vierbein_gradient[1];

  // gradient in z1
  vierbein_gradient[0] = 2.0;

  kernel_conf_poisson_tensor_vierbein_1v(vierbein_inv,
    vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0 / (2*z1), eps));

}

void eval_vierbein_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void eval_vierbein_gradient_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
test_triad_1x1v_flat_conf(int poly_order)
{

  double lower[] = {0.1, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom;
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_1v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_1v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[2], sqrt(0.0), 1e-12) );
  }


  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}




void eval_vierbein_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;
}

void eval_vierbein_gradient_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;

  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
test_triad_1x2v_flat_conf(int poly_order)
{

  double lower[] = {0.1, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom;
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_2v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    /* Shape of conf_poisson_tensor non-zero elements.
    // xx - components are not used by 1x2v, but might still be saved/ non-zero
    // -- - compoentns are zero/implied by symmetry
    // Used - components actualled used by the kernels 
    // 
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | 00 | 01 |     | -- | xx |Used|Used|     
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | 02 | 03 |     | xx | xx | xx | xx |     
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | -- |04+j|     | -- | xx | -- |Used|      
    +----+----+----+----+     +----+----+----+----+     
    | -- | -- | -- | -- |     | -- | xx | -- | -- |     
    +----+----+----+----+     +----+----+----+----+     
    */

    // pt at the first and second points of conflocal
    double conf_pt_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0,    // PT-Comp 00
      0.0, 0.0, 0.0, // PT-Comp 01
      0.0, 0.0, 0.0, // PT-Comp 02
      1.4142135623730951e+00, 0.0, 0.0, // PT-Comp 03
      0.0, 0.0, 0.0, // PT-Comp 04 + 0
      0.0, 0.0, 0.0  }; // PT-Comp 05

    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = conf_pt_vals[test_idx];
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }


  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void eval_vierbein_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = q_r;
}

void eval_vierbein_gradient_annulus_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0;

  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
test_triad_1x2v_annulus_conf(int poly_order)
{

  double lower[] = {0.1, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

    struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Construct the input map
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom;
  inp_triad_geom.use_vierbein = true;
  inp_triad_geom.use_preset_geom = false;
  inp_triad_geom.eval_vierbein = eval_vierbein_annulus_2v; 
  inp_triad_geom.eval_vierbein_gradient = eval_vierbein_gradient_annulus_2v; 
  inp_triad_geom.eval_vierbein_ctx = 0; 
  inp_triad_geom.eval_vierbein_gradient_ctx = 0;  

  // Make the memory for arrays (modal)
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_triad_geom, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric vierbein
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      5.6865230654862211e+00, -3.3402132856134248e+00, 1.1941468087349127e+00, 0.0, 0.0, 0.0,
      5.6865230654862211e+00, -3.3402132856134248e+00, 1.1941468087349127e+00  };
    double conf_pt_pnt2_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.8807796203407585e+00, -3.3402132856134242e-01, 5.0077124237270627e-02, 0.0, 0.0, 0.0,
      1.8807796203407585e+00, -3.3402132856134242e-01, 5.0077124237270627e-02  };

    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt2_vals[test_idx];
        //printf("conf_poisson_tensor_d[%d]: %1.16e\n",test_idx,conf_poisson_tensor_d[test_idx]);
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }

  // Release the memory
  gkyl_array_release(conf_poisson_tensor);
}

void
test_triad_math_2v()
{
  // Test: Can the correct bracket be constructed for:
  // 2D Annular Disk Coordiantes
  int vdim = 2;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[4];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;

  // Assign covaraint tangent basis 
  vierbein[0] = 1.0;
  vierbein[1] = 0.0;
  vierbein[2] = 0.0;
  vierbein[3] = r;

  // 1. Test kernel_vierbein_inv_2v()
  double vierbein_inv[4];
  kernel_vierbein_inv_2v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[3], 1.0/r, eps));

  // 2. Test kernel_conf_poisson_tensor_2v()
  double conf_poisson_tensor[6];

  // Assign the triad basis gradient
  double vierbein_gradient[8];

  // gradient in r
  vierbein_gradient[0] = 0.0;
  vierbein_gradient[1] = 0.0;

  vierbein_gradient[2] = 0.0;
  vierbein_gradient[3] = 1.0;

  // gradient in theta
  vierbein_gradient[4] = 0.0;
  vierbein_gradient[5] = 0.0;

  vierbein_gradient[6] = 0.0;
  vierbein_gradient[7] = 0.0;

  kernel_conf_poisson_tensor_vierbein_2v(vierbein_inv, vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[3], 1/r, eps));

  // Pi_{pp} block, p_r coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[4], 0.0, eps));

  // Pi_{pp} block, p_theta coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[5], 1.0/r, eps));
}

void
test_triad_math_3v()
{

  // Test: Can the correct bracket be constructed for:
  // 3D Spherical Coordiantes
  int vdim = 3;

  // Accuracy of the tests
  double eps = 1e-12;
  double vierbein[9];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;
  double phi = 0.8;

  // Assign covaraint tangent basis 
  vierbein[0] = 1.0;
  vierbein[1] = 0.0;
  vierbein[2] = 0.0;

  vierbein[3] = 0.0;
  vierbein[4] = r;
  vierbein[5] = 0.0;

  vierbein[6] = 0.0;
  vierbein[7] = 0.0;
  vierbein[8] = r * sin(theta);

  // 1. Test kernel_vierbein_inv_3v()
  double vierbein_inv[9];
  kernel_vierbein_inv_3v(vierbein, vierbein_inv);

  // vierbein_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(vierbein_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[4], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(vierbein_inv[8], 1.0/(r*sin(theta)), eps));

  // 5. Test kernel_conf_poisson_tensor_3v()
  double conf_poisson_tensor[18];

  // Assign the triad basis gradient
  double vierbein_gradient[27];

  // gradient in r
  vierbein_gradient[0] = 0.0;
  vierbein_gradient[1] = 0.0;
  vierbein_gradient[2] = 0.0;

  vierbein_gradient[3] = 0.0;
  vierbein_gradient[4] = 1.0;
  vierbein_gradient[5] = 0.0;

  vierbein_gradient[6] = 0.0;
  vierbein_gradient[7] = 0.0;
  vierbein_gradient[8] = sin(theta);

  // gradient in theta
  vierbein_gradient[9] = 0.0;
  vierbein_gradient[10] = 0.0;
  vierbein_gradient[11] = 0.0;

  vierbein_gradient[12] = 0.0;
  vierbein_gradient[13] = 0.0;
  vierbein_gradient[14] = 0.0;

  vierbein_gradient[15] = 0.0;
  vierbein_gradient[16] = 0.0;
  vierbein_gradient[17] = r*cos(theta);

  // gradient in phi
  vierbein_gradient[18] = 0.0;
  vierbein_gradient[19] = 0.0;
  vierbein_gradient[20] = 0.0;

  vierbein_gradient[21] = 0.0;
  vierbein_gradient[22] = 0.0;
  vierbein_gradient[23] = 0.0;

  vierbein_gradient[24] = 0.0;
  vierbein_gradient[25] = 0.0;
  vierbein_gradient[26] = 0.0;

  kernel_conf_poisson_tensor_vierbein_3v(vierbein_inv, vierbein_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[4], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[8], 1.0/(r*sin(theta)), eps));

  // Pi_{pp} block, p_r coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[9], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[10], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[11], 0.0, eps));

  // Pi_{pp} block, p_theta coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[12], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[13], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[14], 0.0, eps));

  // Pi_{pp} block, p_phi coeffs 
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[15], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[16], 1.0/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[17], cos(theta)/(r*sin(theta)), eps));

}

void test_triad_1v() { test_triad_math_1v(); }
void test_triad_1x1v_flat() { test_triad_1x1v_flat_conf(2); }
void test_triad_1x2v_flat() { test_triad_1x2v_flat_conf(2); }
void test_triad_1x2v_annulus() { test_triad_1x2v_annulus_conf(2); }
void test_triad_2v() { test_triad_math_2v(); }
void test_triad_3v() { test_triad_math_3v(); }

TEST_LIST = {
  { "test_triad_1v", test_triad_1v}, 
  { "test_triad_1x1v_flat", test_triad_1x1v_flat}, 
  { "test_triad_1x2v_flat", test_triad_1x2v_flat}, 
  { "test_triad_1x2v_annulus", test_triad_1x2v_annulus}, 
  { "test_triad_2v", test_triad_2v},
  { "test_triad_3v", test_triad_3v},
  {NULL, NULL}
};
