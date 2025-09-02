#include <acutest.h>

#include <math.h>

#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>

void
test_triad_math(int vdim, int cdim)
{
  // Test: Can the correct bracket be constructed for:
  // 3D Spherical Coordiantes

  // Accuracy of the tests
  double eps = 1e-12;

  // 1. Test kernel_metric_3x()
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

  kernel_metric_3x(cov_tangent_basis, h_ij);

  // h_ij only has symmetric entries
  TEST_CHECK(gkyl_compare_double(h_ij[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[3], r*r, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij[5], r*r*sin(theta)*sin(theta), eps));

  // 2. Test kernel_metric_inv_3x()
  double h_ij_inv[6];
  kernel_metric_inv_3x(h_ij, h_ij_inv);

  // h_ij_inv only has symmetric entries as well
  TEST_CHECK(gkyl_compare_double(h_ij_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[3], 1.0/(r*r), eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[4], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(h_ij_inv[5], 1.0/(r*r*sin(theta)*sin(theta)), eps));

  // 3. Test kernel_metric_det_3x()
  double det_h[1];
  kernel_metric_det_3x(h_ij, det_h);

  TEST_CHECK(gkyl_compare_double(det_h[0], r*r*sin(theta), eps));

  // 4. Test compute_nu_inv_3x()

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

  compute_nu_inv_3x(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  TEST_CHECK(gkyl_compare_double(nu_inv[0], 1.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[1], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[2], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[3], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[4], 1/r, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[5], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[6], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[7], 0.0, eps));
  TEST_CHECK(gkyl_compare_double(nu_inv[8], 1/(r*sin(theta)), eps));

  // 5. Test kernel_conf_poisson_tensor_3x3v()
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

  kernel_conf_poisson_tensor_3x3v(h_ij_inv, triad_basis, cov_tangent_basis,
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

void test_triad_3x3v() { test_triad_math(3, 3); }

TEST_LIST = {
  { "test_triad_3x3v", test_triad_3x3v},
  {NULL, NULL}
};
