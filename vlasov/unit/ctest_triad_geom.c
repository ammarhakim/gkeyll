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
test_triad_math_1v(int vdim)
{
  // Test: Can the correct bracket be constructed for:
  // 1D streched coordinates: x(z^1) = (z^1)^2 \sigma_x

  // Accuracy of the tests
  double eps = 1e-12;

  // 1. Test kernel_metric_1v()
  double cov_tangent_basis[1];
  double h_ij[1];

  // Randomly choosen coordinates
  double z1 = 3.2;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = 2*z1;

  kernel_metric_1v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 4.0 * z1 * z1, eps));

  // 2. Test kernel_metric_inv_1v()
  double h_ij_inv[1];
  kernel_metric_inv_1v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0 / (4.0 * z1 * z1), eps));

  // 3. Test kernel_metric_det_1v()
  double det_h[1];
  kernel_metric_det_1v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], 2.0 * z1, eps));

  // 4. Test compute_nu_inv_1v()

  double triad_basis[1];
  double nu_inv[1];

  // Assign the triad basis 
  triad_basis[0] = 1;

  compute_nu_inv_1v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0 / (2.0 * z1), eps));

  // 5. Test kernel_conf_poisson_tensor_1v()
  double conf_poisson_tensor[1];

  // Assign the triad basis gradient
  double triad_basis_gradient[1];

  // gradient in r
  triad_basis_gradient[0] = 0.0;

  kernel_conf_poisson_tensor_1v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0 / (2*z1), eps));

}

void eval_triad_basis_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1;
}

void eval_cov_tangent_basis_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1;
}

void eval_triad_basis_gradient_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0;
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
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors;
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis_1v; 
  inp_basis_vectors.eval_triad_basis = eval_triad_basis_1v; 
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient_1v; 
  inp_basis_vectors.ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *h_ij; // Specified metric inverse for canonical poisson bracket
  struct gkyl_array *h_ij_inv; // Specified metric inverse for canonical poisson bracket
  struct gkyl_array *det_h; // Specified metric determinant
  struct gkyl_array *cov_tangent_basis; // Covariant tangent basis
  struct gkyl_array *triad_basis; // Triad basis
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  cov_tangent_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  triad_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  h_ij = mkarr(basis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  h_ij_inv = mkarr(basis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(basis.num_basis, local_ext.volume);
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_basis_vectors, cov_tangent_basis,
    triad_basis, h_ij, h_ij_inv, det_h, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *h_ij_d = gkyl_array_cfetch(h_ij, gkyl_range_idx(&confLocal, iter.idx));
    const double *h_ij_inv_d = gkyl_array_cfetch(h_ij_inv, gkyl_range_idx(&confLocal, iter.idx));
    const double *det_h_d = gkyl_array_cfetch(det_h, gkyl_range_idx(&confLocal, iter.idx));
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;
    TEST_CHECK( gkyl_compare_double(h_ij_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(h_ij_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(h_ij_d[2], sqrt(0.0), 1e-12) );

    TEST_CHECK( gkyl_compare_double(h_ij_inv_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(h_ij_inv_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(h_ij_inv_d[2], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(det_h_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(det_h_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(det_h_d[2], sqrt(0.0), 1e-12) );

    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[0], sqrt(2.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[1], sqrt(0.0), 1e-12) );
    TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[2], sqrt(0.0), 1e-12) );
  }


  // Release the memory
  gkyl_array_release(h_ij);
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(cov_tangent_basis);
  gkyl_array_release(triad_basis);
  gkyl_array_release(conf_poisson_tensor);
}


void eval_triad_basis_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = cos(q_theta);
  fout[1] = sin(q_theta);
  fout[2] = - q_r * sin(q_theta);
  fout[3] = q_r * cos(q_theta);
}

void eval_cov_tangent_basis_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = cos(q_theta);
  fout[1] = sin(q_theta);
  fout[2] = - sin(q_theta);
  fout[3] = cos(q_theta);
}

void eval_triad_basis_gradient_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double q_r = xn[0];
  double q_theta = 0.0;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;

  fout[4] = sin(q_theta);
  fout[5] = cos(q_theta);
  fout[6] = - cos(q_theta);
  fout[7] = sin(q_theta);
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
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors;
  inp_basis_vectors.eval_cov_tangent_basis = eval_cov_tangent_basis_2v; 
  inp_basis_vectors.eval_triad_basis = eval_triad_basis_2v; 
  inp_basis_vectors.eval_triad_basis_gradient = eval_triad_basis_gradient_2v; 
  inp_basis_vectors.ctx = 0; 

  // Make the memory for arrays (modal)
  struct gkyl_array *h_ij; // Specified metric inverse for canonical poisson bracket
  struct gkyl_array *h_ij_inv; // Specified metric inverse for canonical poisson bracket
  struct gkyl_array *det_h; // Specified metric determinant
  struct gkyl_array *cov_tangent_basis; // Covariant tangent basis
  struct gkyl_array *triad_basis; // Triad basis
  struct gkyl_array *conf_poisson_tensor; // Configuration space Poisson tensor representation

  // Size of the PT
  int num_pt_indices[3] = { 1 , 6, 18 }; 

  // Allocate arrays for covariant tangent basis 
  cov_tangent_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  triad_basis = mkarr(basis.num_basis*vdim*vdim, local_ext.volume);
  h_ij = mkarr(basis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  h_ij_inv = mkarr(basis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(basis.num_basis, local_ext.volume);
  conf_poisson_tensor = mkarr(basis.num_basis*num_pt_indices[vdim-1], local_ext.volume);

  // Construct the Geometry for this configuration
  gkyl_vlasov_triad_geom_new(&confGrid, &confLocal, confBasis, 
    &grid, &local, basis, inp_basis_vectors, cov_tangent_basis,
    triad_basis, h_ij, h_ij_inv, det_h, conf_poisson_tensor);

  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *h_ij_d = gkyl_array_cfetch(h_ij, gkyl_range_idx(&confLocal, iter.idx));
    const double *h_ij_inv_d = gkyl_array_cfetch(h_ij_inv, gkyl_range_idx(&confLocal, iter.idx));
    const double *det_h_d = gkyl_array_cfetch(det_h, gkyl_range_idx(&confLocal, iter.idx));
    const double *conf_poisson_tensor_d = gkyl_array_cfetch(conf_poisson_tensor, gkyl_range_idx(&confLocal, iter.idx));
    
    int NC = confBasis.num_basis;

    // index size of the symmetric h_ij
    int n_sym = vdim*(vdim+1)/2; 

    // Double precision values of the expected output
    double h_ij_vals[9] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4142135623730951e+00, 0.0, 0.0 };
    double h_ij_inv_vals[9] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4142135623730951e+00, 0.0, 0.0 };
    double det_h_vals[3] = { 1.4142135623730951e+00, 0.0, 0.0 };

    // pt at the first and second points of conflocal
    double conf_pt_pnt1_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      4.5961940777125598e-01, 1.8371173070873839e-01, -1.3360012858723191e-17, 0.0, 0.0, 0.0,
      4.5961940777125598e-01, 1.8371173070873839e-01, -1.3360012858723191e-17  };
    double conf_pt_pnt2_vals[18] = { 1.4142135623730951e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0960155108391487e+00, 1.8371173070873839e-01, -1.3360012858723191e-17, 0.0, 0.0, 0.0,
      1.0960155108391487e+00, 1.8371173070873839e-01, -1.3360012858723191e-17  };


    // Test checks
    for (int k = 0; k<n_sym; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        TEST_CHECK( gkyl_compare_double(h_ij_d[test_idx],     h_ij_vals[test_idx],     1e-12) );
        TEST_CHECK( gkyl_compare_double(h_ij_inv_d[test_idx], h_ij_inv_vals[test_idx],  1e-12) );
      }
    }

    for (int c = 0; c < NC; ++c) {
      TEST_CHECK( gkyl_compare_double(det_h_d[c], det_h_vals[c], 1e-12) );
    }

    for (int k = 0; k<num_pt_indices[vdim-1]; ++k) {
      for (int m = 0; m<NC; ++m) {
        int test_idx = k*NC + m;
        double expected = 0;
        if (iter.idx[0] == 1) expected = conf_pt_pnt1_vals[test_idx];
        if (iter.idx[0] == 2) expected = conf_pt_pnt2_vals[test_idx];
        TEST_CHECK( gkyl_compare_double(conf_poisson_tensor_d[test_idx], expected, 1e-12) );
      }
    }
  }


  // Release the memory
  gkyl_array_release(h_ij);
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(cov_tangent_basis);
  gkyl_array_release(triad_basis);
  gkyl_array_release(conf_poisson_tensor);
}

void
test_triad_math_2v(int vdim)
{
  // Test: Can the correct bracket be constructed for:
  // 2D Annular Disk Coordiantes

  // Accuracy of the tests
  double eps = 1e-12;

  // 1. Test kernel_metric_2v()
  double cov_tangent_basis[4];
  double h_ij[3];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = cos(theta);
  cov_tangent_basis[1] = sin(theta);
  cov_tangent_basis[2] = - r * sin(theta);
  cov_tangent_basis[3] = r * cos(theta);

  kernel_metric_2v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[2], r*r, eps));

  // 2. Test kernel_metric_inv_2v()
  double h_ij_inv[3];
  kernel_metric_inv_2v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[2], 1.0/(r*r), eps));

  // 3. Test kernel_metric_det_2v()
  double det_h[1];
  kernel_metric_det_2v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], r, eps));

  // 4. Test compute_nu_inv_2v()

  double triad_basis[4];
  double nu_inv[4];

  // Assign the triad basis 
  triad_basis[0] = cos(theta);
  triad_basis[1] = sin(theta);
  triad_basis[2] = - sin(theta);
  triad_basis[3] = cos(theta);

  compute_nu_inv_2v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[3], 1/r, eps));

  // 5. Test kernel_conf_poisson_tensor_2v()
  double conf_poisson_tensor[6];

  // Assign the triad basis gradient
  double triad_basis_gradient[8];

  // gradient in r
  triad_basis_gradient[0] = 0.0;
  triad_basis_gradient[1] = 0.0;

  triad_basis_gradient[2] = 0.0;
  triad_basis_gradient[3] = 0.0;

  // gradient in theta
  triad_basis_gradient[4] = -sin(theta);
  triad_basis_gradient[5] = cos(theta);

  triad_basis_gradient[6] = -cos(theta);
  triad_basis_gradient[7] = -sin(theta);

  kernel_conf_poisson_tensor_2v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

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
test_triad_math_3v(int vdim)
{
  // Test: Can the correct bracket be constructed for:
  // 3D Spherical Coordiantes

  // Accuracy of the tests
  double eps = 1e-12;

  // 1. Test kernel_metric_3v()
  double cov_tangent_basis[9];
  double h_ij[6];

  // Randomly choosen coordinates
  double r = 3.2;
  double theta = 1.4;
  double phi = 0.8;

  // Assign covaraint tangent basis 
  cov_tangent_basis[0] = sin(theta) * cos(phi);
  cov_tangent_basis[1] = sin(theta) * sin(phi);
  cov_tangent_basis[2] = cos(theta);

  cov_tangent_basis[3] = r * cos(theta) * cos(phi);
  cov_tangent_basis[4] = r * cos(theta) * sin(phi);
  cov_tangent_basis[5] = -r * sin(theta);

  cov_tangent_basis[6] = - r * sin(theta) * sin(phi);
  cov_tangent_basis[7] = r * sin(theta) * cos(phi);
  cov_tangent_basis[8] = 0.0;

  kernel_metric_3v(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[3], r*r, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[5], r*r*sin(theta)*sin(theta), eps));

  // 2. Test kernel_metric_inv_3v()
  double h_ij_inv[6];
  kernel_metric_inv_3v(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[3], 1.0/(r*r), eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[5], 1.0/(r*r*sin(theta)*sin(theta)), eps));

  // 3. Test kernel_metric_det_3v()
  double det_h[1];
  kernel_metric_det_3v(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], r*r*sin(theta), eps));

  // 4. Test compute_nu_inv_3v()

  double triad_basis[9];
  double nu_inv[9];

  // Assign the triad basis 
  triad_basis[0] = sin(theta) * cos(phi);
  triad_basis[1] = sin(theta) * sin(phi);
  triad_basis[2] = cos(theta);

  triad_basis[3] = cos(theta) * cos(phi);
  triad_basis[4] = cos(theta) * sin(phi);
  triad_basis[5] = - sin(theta);

  triad_basis[6] = - sin(phi);
  triad_basis[7] = cos(phi);
  triad_basis[8] = 0.0;

  compute_nu_inv_3v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[4], 1/r, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[8], 1/(r*sin(theta)), eps));

  // 5. Test kernel_conf_poisson_tensor_3v()
  double conf_poisson_tensor[18];

  // Assign the triad basis gradient
  double triad_basis_gradient[27];

  // gradient in r
  triad_basis_gradient[0] = 0.0;
  triad_basis_gradient[1] = 0.0;
  triad_basis_gradient[2] = 0.0;

  triad_basis_gradient[3] = 0.0;
  triad_basis_gradient[4] = 0.0;
  triad_basis_gradient[5] = 0.0;

  triad_basis_gradient[6] = 0.0;
  triad_basis_gradient[7] = 0.0;
  triad_basis_gradient[8] = 0.0;

  // gradient in theta
  triad_basis_gradient[9] = cos(phi)*cos(theta);
  triad_basis_gradient[10] = cos(theta)*sin(phi);
  triad_basis_gradient[11] = -sin(theta);

  triad_basis_gradient[12] = -cos(phi)*sin(theta);
  triad_basis_gradient[13] = -sin(phi)*sin(theta);
  triad_basis_gradient[14] = -cos(theta);

  triad_basis_gradient[15] = 0.0;
  triad_basis_gradient[16] = 0.0;
  triad_basis_gradient[17] = 0.0;

  // gradient in phi
  triad_basis_gradient[18] = -sin(phi)*sin(theta);
  triad_basis_gradient[19] = cos(phi)*sin(theta);
  triad_basis_gradient[20] = 0.0;

  triad_basis_gradient[21] = -cos(theta)*sin(phi);
  triad_basis_gradient[22] = cos(phi)*cos(theta);
  triad_basis_gradient[23] = 0.0;

  triad_basis_gradient[24] = -cos(phi);
  triad_basis_gradient[25] = -sin(phi);
  triad_basis_gradient[26] = 0.0;

  kernel_conf_poisson_tensor_3v(h_ij_inv, triad_basis, cov_tangent_basis,
    triad_basis_gradient, conf_poisson_tensor);

  // Pi_{xx} block
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[4], 1/r, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(conf_poisson_tensor[8], 1/(r*sin(theta)), eps));

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

void test_triad_1v() { test_triad_math_1v(1); }
void test_triad_1x1v_flat() { test_triad_1x1v_flat_conf(2); }
void test_triad_1x2v_annulus() { test_triad_1x2v_annulus_conf(2); }
void test_triad_2v() { test_triad_math_2v(2); }
void test_triad_3v() { test_triad_math_3v(3); }

TEST_LIST = {
  { "test_triad_1v", test_triad_1v}, 
  { "test_triad_1x1v_flat", test_triad_1x1v_flat}, 
  { "test_triad_1x2v_annulus", test_triad_1x2v_annulus}, 
  { "test_triad_2v", test_triad_2v},
  { "test_triad_3v", test_triad_3v},
  {NULL, NULL}
};
