// Unit tests for the characteristic-based outflow flux kernel
// outflow_flux_x_2x_ser_p1 (and y variant). The kernel's job is to compute
//   F_outflow = A^outgoing . U_skin
// where A^outgoing is the GR Maxwell flux Jacobian projected onto the
// outgoing characteristics at a domain boundary. We test two consistency
// properties of the kernel that don't depend on getting the eigenvector math
// "right" externally — they just check that the kernel constructs the right
// outgoing-eigenvalue vector and dispatches to A_dQ_<dir>_calc as intended:
//
//   (A) When the metric is set so all three eigenvalues l1, l2, l3 at a face
//       node are POSITIVE (e.g., a strongly negative shift), edge=upper
//       should keep all three (A^outgoing = A_full); the output should
//       equal a direct call to A_dQ_x_calc with the full eigenvalues.
//       Conversely edge=lower should produce a ZERO flux at that node
//       (no negative eigenvalues to keep).
//
//   (B) Mirror: all eigenvalues negative -> edge=lower keeps all,
//       edge=upper produces zero.
//
// These together exercise the full kernel pipeline (modal->nodal eval,
// eigenvalue construction, A_dQ dispatch, output layout, cleaning-slot
// zeroing) against the kernel's own intended behavior, and they catch any
// dispatch-table or signature mistakes.

#include <acutest.h>
#include <math.h>

#include <gkyl_dg_gr_maxwell_kernels.h>

#define NUM_QUAD_FACE 2     // (p+1) Gauss nodes per surface in 2D p=1
#define FACE_FLUX_SIZE 16   // 8 components * 2 quad nodes for a 2x_ser_p1 face slot

// Direct evaluation of a 2x ser_p1 modal expansion at a single (xi, eta) node.
static inline double
mode_at(const double *c, double xi, double eta)
{
  // Basis: {1/2, sqrt(3)/2 * xi, sqrt(3)/2 * eta, 3/2 * xi*eta}.
  return 0.5*c[0]
       + 0.8660254037844386*xi*c[1]
       + 0.8660254037844386*eta*c[2]
       + 1.5*xi*eta*c[3];
}

// Build a manual reference: reproduce what outflow_flux_x_2x_ser_p1 should
// compute. This is intentionally an independent reimplementation so that a
// kernel bug that happens to also exist in the test is unlikely.
static void
reference_outflow_flux_x(int edge,
  const double *lapse_n, const double *shift_n, const double *h_ij_n,
  const double *det_h_n,
  const double *field_skin,
  double out_ref[FACE_FLUX_SIZE])
{
  for (int k = 0; k < FACE_FLUX_SIZE; ++k) out_ref[k] = 0.0;

  static const double xi_face[2]    = { -1.0, +1.0 };
  static const double eta_quad[2]   = { -0.5773502691896258, 0.5773502691896258 };

  double xi = xi_face[edge];
  for (int i = 0; i < NUM_QUAD_FACE; ++i) {
    double eta = eta_quad[i];

    double U_n[6];
    for (int comp = 0; comp < 6; ++comp) {
      U_n[comp] = mode_at(&field_skin[comp*4], xi, eta);
    }

    double sx = shift_n[0*2 + i], sy = shift_n[1*2 + i], sz = shift_n[2*2 + i];
    double hxx = h_ij_n[0*2 + i], hxy = h_ij_n[1*2 + i], hxz = h_ij_n[2*2 + i];
    double hyy = h_ij_n[3*2 + i], hyz = h_ij_n[4*2 + i], hzz = h_ij_n[5*2 + i];
    double alpha = lapse_n[i];
    double Jc = det_h_n[i];

    double sqrt_term = sqrt(hyy*hzz - hyz*hyz)/Jc;
    double l1 = -sx;
    double l2 = -sx + alpha*sqrt_term;
    double l3 = -sx - alpha*sqrt_term;

    double lambda[6];
    if (edge == 0) {
      lambda[0] = (l1 < 0.0) ? l1 : 0.0;
      lambda[1] = (l1 < 0.0) ? l1 : 0.0;
      lambda[2] = (l2 < 0.0) ? l2 : 0.0;
      lambda[3] = (l2 < 0.0) ? l2 : 0.0;
      lambda[4] = (l3 < 0.0) ? l3 : 0.0;
      lambda[5] = (l3 < 0.0) ? l3 : 0.0;
    } else {
      lambda[0] = (l1 > 0.0) ? l1 : 0.0;
      lambda[1] = (l1 > 0.0) ? l1 : 0.0;
      lambda[2] = (l2 > 0.0) ? l2 : 0.0;
      lambda[3] = (l2 > 0.0) ? l2 : 0.0;
      lambda[4] = (l3 > 0.0) ? l3 : 0.0;
      lambda[5] = (l3 > 0.0) ? l3 : 0.0;
    }

    double F_n[6] = {0.0};
    A_dQ_x_calc(alpha, sx, sy, sz, hxx, hxy, hxz, hyy, hyz, hzz, Jc,
      lambda, U_n, F_n);

    out_ref[0  + i] = F_n[0];
    out_ref[2  + i] = F_n[1];
    out_ref[4  + i] = F_n[2];
    out_ref[6  + i] = F_n[3];
    out_ref[8  + i] = F_n[4];
    out_ref[10 + i] = F_n[5];
  }
}

// Test driver: pick a metric and a field, call the kernel + the reference,
// compare term-by-term.
static void
run_xface_consistency(const char *name, int edge,
  const double *lapse_n, const double *shift_n, const double *h_ij_n,
  const double *det_h_n,
  const double *field_skin)
{
  double dxv[2] = { 0.5, 0.5 }; // Cell sizes; only enters the cflrate scaling.

  double out_kernel[FACE_FLUX_SIZE] = {0.0};
  // Pre-fill cleaning slots with garbage to verify the kernel zeroes them.
  out_kernel[12] = out_kernel[13] = out_kernel[14] = out_kernel[15] = 1.0e9;

  double cfl = outflow_flux_x_2x_ser_p1(dxv, edge, lapse_n, shift_n,
    h_ij_n, det_h_n, field_skin, out_kernel);

  double out_ref[FACE_FLUX_SIZE];
  reference_outflow_flux_x(edge, lapse_n, shift_n, h_ij_n, det_h_n,
    field_skin, out_ref);

  // EM components 0..11 must match.
  for (int k = 0; k < 12; ++k) {
    bool ok = fabs(out_kernel[k] - out_ref[k]) < 1.0e-12 * (1.0 + fabs(out_ref[k]));
    if (!ok) {
      printf("[%s] mismatch at out[%d]: kernel=%1.12e ref=%1.12e\n",
        name, k, out_kernel[k], out_ref[k]);
    }
    TEST_CHECK(ok);
  }
  // Cleaning slots must be zeroed.
  for (int k = 12; k < 16; ++k) {
    TEST_CHECK(out_kernel[k] == 0.0);
  }
  TEST_CHECK(cfl >= 0.0);
}

// All-positive eigenvalues: shift_x = -1, lapse and metric small. l1, l2, l3
// all > 0. Outflow at upper edge keeps everything; lower edge produces 0.
static void
test_xface_all_positive(void)
{
  double lapse_n[2] = { 0.1, 0.1 };
  double shift_n[6] = { -1.0, -1.0,  0.0, 0.0,  0.0, 0.0 };
  double h_ij_n[12] = {
    1.0, 1.0,  // h_xx
    0.0, 0.0,  // h_xy
    0.0, 0.0,  // h_xz
    1.0, 1.0,  // h_yy
    0.0, 0.0,  // h_yz
    1.0, 1.0,  // h_zz
  };
  double det_h_n[2] = { 1.0, 1.0 };
  // Field: nonzero in JD^x mode 0 only.
  double field_skin[8*4] = {0.0};
  field_skin[0*4 + 0] = 1.0;
  field_skin[3*4 + 0] = 0.5; // also some JB^x

  run_xface_consistency("all-positive, edge=upper", 1,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
  run_xface_consistency("all-positive, edge=lower", 0,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
}

// All-negative eigenvalues: shift_x = +1, same metric. l1 = -1, l2 ≈ -1+0.1, l3 ≈ -1-0.1.
// Outflow at lower edge keeps everything; upper edge produces 0.
static void
test_xface_all_negative(void)
{
  double lapse_n[2] = { 0.1, 0.1 };
  double shift_n[6] = { 1.0, 1.0,  0.0, 0.0,  0.0, 0.0 };
  double h_ij_n[12] = {
    1.0, 1.0,
    0.0, 0.0,
    0.0, 0.0,
    1.0, 1.0,
    0.0, 0.0,
    1.0, 1.0,
  };
  double det_h_n[2] = { 1.0, 1.0 };
  double field_skin[8*4] = {0.0};
  field_skin[2*4 + 0] = 1.0; // JD^z mode 0
  field_skin[4*4 + 0] = -0.5; // JB^y mode 0

  run_xface_consistency("all-negative, edge=lower", 0,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
  run_xface_consistency("all-negative, edge=upper", 1,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
}

// Mixed-sign eigenvalues: shift_x = +0.3, lapse = 1, l1 = -0.3, l2 = +0.7, l3 = -1.3.
// Both edge=lower and edge=upper produce nontrivial nonzero output (different
// projections of A onto outgoing characteristics). This exercises the
// general case the simulation will see in the ergoregion.
static void
test_xface_mixed_sign(void)
{
  double lapse_n[2] = { 1.0, 1.0 };
  double shift_n[6] = { 0.3, 0.3,  0.0, 0.0,  0.0, 0.0 };
  double h_ij_n[12] = {
    1.0, 1.0,
    0.0, 0.0,
    0.0, 0.0,
    1.0, 1.0,
    0.0, 0.0,
    1.0, 1.0,
  };
  double det_h_n[2] = { 1.0, 1.0 };
  double field_skin[8*4] = {0.0};
  // Excite multiple modes and components for a more stringent test.
  field_skin[0*4 + 0] = 0.7;  field_skin[0*4 + 1] = 0.1;
  field_skin[1*4 + 0] = -0.4; field_skin[1*4 + 2] = 0.05;
  field_skin[2*4 + 0] = 0.2;  field_skin[2*4 + 3] = 0.02;
  field_skin[3*4 + 0] = 0.6;  field_skin[3*4 + 1] = -0.03;
  field_skin[4*4 + 0] = -0.5; field_skin[4*4 + 2] = 0.04;
  field_skin[5*4 + 0] = 0.9;  field_skin[5*4 + 3] = -0.06;

  run_xface_consistency("mixed-sign, edge=upper", 1,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
  run_xface_consistency("mixed-sign, edge=lower", 0,
    lapse_n, shift_n, h_ij_n, det_h_n, field_skin);
}

// Sanity: with U_skin == 0 the outflow flux must be identically zero.
static void
test_xface_zero_field(void)
{
  double lapse_n[2] = { 1.0, 1.0 };
  double shift_n[6] = { 0.3, 0.3,  0.0, 0.0,  0.0, 0.0 };
  double h_ij_n[12] = {
    1.5, 1.5,  0.0, 0.0,  0.0, 0.0,
    2.0, 2.0,  0.0, 0.0,  3.0, 3.0,
  };
  double det_h_n[2] = { sqrt(1.5*2.0*3.0), sqrt(1.5*2.0*3.0) };
  double field_skin[8*4] = {0.0};

  double out[FACE_FLUX_SIZE] = {0.0};
  for (int k = 12; k < 16; ++k) out[k] = 1.0e9;
  double dxv[2] = { 1.0, 1.0 };
  double cfl = outflow_flux_x_2x_ser_p1(dxv, 0, lapse_n, shift_n,
    h_ij_n, det_h_n, field_skin, out);
  for (int k = 0; k < 12; ++k) TEST_CHECK(out[k] == 0.0);
  for (int k = 12; k < 16; ++k) TEST_CHECK(out[k] == 0.0);
  TEST_CHECK(cfl >= 0.0);
}

TEST_LIST = {
  { "xface_all_positive", test_xface_all_positive },
  { "xface_all_negative", test_xface_all_negative },
  { "xface_mixed_sign",   test_xface_mixed_sign },
  { "xface_zero_field",   test_xface_zero_field },
  { NULL, NULL },
};
