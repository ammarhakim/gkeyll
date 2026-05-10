// Per-quadrature-point correctness tests for the DG GR Maxwell tetrad-Roe
// surface-flux kernels.
//
// Mirrors the pattern of moments/unit/ctest_wv_gr_maxwell_tetrad.c (which
// verifies the finite-volume tetrad-Roe scheme's algebraic identities) but
// at the DG-quadrature-point level: the math we test here is what runs at
// each face GL node inside lax_flux_x_2x_ser_p1_tetrad_roe.c and the y-face
// analog.
//
// What we verify:
//   1. M^T M = h_ij at the face GL node (correct tetrad transform).
//   2. M M^{-1} = I (correct closed-form inverse formula).
//   3. Paper eq. 57 right eigenvectors mutually orthogonal in flat tetrad
//      (with the expected norms), for both x- and y-direction Maxwell.
//   4. Sum_p alpha_p r_p = dU_t (wave decomposition is complete, i.e., the 6
//      eigenvectors span the conserved-variable space).
//   5. Tetrad-Roe kernel produces zero dissipation when U_L = U_R (zero-jump
//      sanity check).
//   6. Tetrad-Roe kernel reduces to the expected flat-Maxwell Roe flux in the
//      Cartesian flat limit (h_ij = delta, alpha = 1, beta = 0).
//
// Test points span Cartesian flat, diagonal stretched, Kerr-Schild interior,
// and Kerr-Schild near-pole. The near-pole test deliberately stays slightly
// off the pole (sin theta = 0.1) so the y-face M^{-1} is still finite; the
// at-the-pole behavior is handled separately by the kernel's theta_pole=1
// reflective fallback (covered by test 5/6 with theta_pole=1).

#include <acutest.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_util.h>

// ---- Tetrad-matrix helpers (same closed forms as the production kernels) ----

// Radial (x) face: h_rtheta = h_thetaphi = 0 case from
// lax_flux_x_2x_ser_p1_tetrad_roe.c.
static void
build_M_x_face(double h_rp, double h_yy, double h_zz, double Jc,
  double M[3][3], double Minv[3][3])
{
  // h_yz = 0 in our metric so the sqrt simplifies.
  double sqrt_gamma_xx = sqrt(h_yy*h_zz) / Jc;
  double inv_sqrt_gamma_xx = 1.0 / sqrt_gamma_xx;
  double sqrt_h_yy = sqrt(h_yy);
  double sqrt_h_zz = sqrt(h_zz);

  // M^a_i (rows = tetrad index a, cols = coord index i).
  M[0][0] = inv_sqrt_gamma_xx; M[0][1] = 0.0;       M[0][2] = 0.0;
  M[1][0] = 0.0;               M[1][1] = sqrt_h_yy; M[1][2] = 0.0;
  M[2][0] = h_rp / sqrt_h_zz;  M[2][1] = 0.0;       M[2][2] = sqrt_h_zz;

  // (M^{-1})^i_a (rows = i, cols = a).
  Minv[0][0] = sqrt_gamma_xx;                Minv[0][1] = 0.0;            Minv[0][2] = 0.0;
  Minv[1][0] = 0.0;                          Minv[1][1] = 1.0/sqrt_h_yy;  Minv[1][2] = 0.0;
  Minv[2][0] = -h_rp * sqrt_gamma_xx / h_zz; Minv[2][1] = 0.0;            Minv[2][2] = 1.0/sqrt_h_zz;
}

// Theta (y) face: same metric assumption, from
// lax_flux_y_2x_ser_p1_tetrad_roe.c.
static void
build_M_y_face(double h_rr, double h_rp, double h_tt, double h_pp,
  double M[3][3], double Minv[3][3])
{
  double D_theta = sqrt(h_rr*h_pp - h_rp*h_rp);
  double sqrt_h_rr = sqrt(h_rr);
  double sqrt_h_tt = sqrt(h_tt);

  // tetrad x = e_theta (face normal); tetrad y = e_r (in-face); tetrad z =
  // e_phi orthogonalized within Sigma_theta.
  M[0][0] = 0.0;               M[0][1] = sqrt_h_tt; M[0][2] = 0.0;
  M[1][0] = sqrt_h_rr;         M[1][1] = 0.0;       M[1][2] = h_rp / sqrt_h_rr;
  M[2][0] = 0.0;               M[2][1] = 0.0;       M[2][2] = D_theta / sqrt_h_rr;

  Minv[0][0] = 0.0;          Minv[0][1] = 1.0 / sqrt_h_rr; Minv[0][2] = -h_rp / (sqrt_h_rr * D_theta);
  Minv[1][0] = 1.0/sqrt_h_tt; Minv[1][1] = 0.0;             Minv[1][2] = 0.0;
  Minv[2][0] = 0.0;          Minv[2][1] = 0.0;             Minv[2][2] = sqrt_h_rr / D_theta;
}

// ---- Verification helpers ----

static void
check_M_satisfies_metric(double M[3][3], double h_ij[3][3], double tol,
  const char *label)
{
  double MtM[3][3] = {{0}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int a = 0; a < 3; ++a)
        MtM[i][j] += M[a][i] * M[a][j];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double diff = fabs(MtM[i][j] - h_ij[i][j]);
      if (diff > tol) {
        printf("  FAIL [%s]: M^T M[%d][%d] = %e, h_ij = %e, diff = %e\n",
          label, i, j, MtM[i][j], h_ij[i][j], diff);
      }
      TEST_CHECK(diff < tol);
    }
  }
}

static void
check_M_inverse_correct(double M[3][3], double Minv[3][3], double tol,
  const char *label)
{
  double MMinv[3][3] = {{0}};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      for (int i = 0; i < 3; ++i)
        MMinv[a][b] += M[a][i] * Minv[i][b];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      double expected = (a == b) ? 1.0 : 0.0;
      double diff = fabs(MMinv[a][b] - expected);
      if (diff > tol) {
        printf("  FAIL [%s]: (M M^-1)[%d][%d] = %e, expected = %e\n",
          label, a, b, MMinv[a][b], expected);
      }
      TEST_CHECK(diff < tol);
    }
  }
}

static void
check_eigenvectors_orthogonal(double r[6][6], double expected_norm2[6],
  double tol, const char *label)
{
  for (int p = 0; p < 6; ++p) {
    for (int q = 0; q < 6; ++q) {
      double dot = 0.0;
      for (int i = 0; i < 6; ++i) dot += r[p][i] * r[q][i];
      double expected = (p == q) ? expected_norm2[p] : 0.0;
      double diff = fabs(dot - expected);
      if (diff > tol) {
        printf("  FAIL [%s]: r_%d . r_%d = %e, expected = %e\n",
          label, p, q, dot, expected);
      }
      TEST_CHECK(diff < tol);
    }
  }
}

// X-direction wave-decomposition reconstruction (paper eq. 57 conventions).
static void
check_x_reconstruction(double dU_t[6], double tol, const char *label)
{
  double a_m1 = (-dU_t[1] + dU_t[5]) * 0.5;
  double a_m2 = ( dU_t[2] + dU_t[4]) * 0.5;
  double a_01 =   dU_t[3];
  double a_02 =   dU_t[0];
  double a_p1 = ( dU_t[1] + dU_t[5]) * 0.5;
  double a_p2 = (-dU_t[2] + dU_t[4]) * 0.5;

  double rec[6] = {0};
  rec[1] += a_m1 * (-1); rec[5] += a_m1 * (1);    // r_minus_1 = (0,-1, 0,0,0,1)
  rec[2] += a_m2 * (1);  rec[4] += a_m2 * (1);    // r_minus_2 = (0, 0, 1,0,1,0)
  rec[3] += a_01 * (1);                           // r_zero_1  = (0, 0, 0,1,0,0)
  rec[0] += a_02 * (1);                           // r_zero_2  = (1, 0, 0,0,0,0)
  rec[1] += a_p1 * (1);  rec[5] += a_p1 * (1);    // r_plus_1  = (0, 1, 0,0,0,1)
  rec[2] += a_p2 * (-1); rec[4] += a_p2 * (1);    // r_plus_2  = (0, 0,-1,0,1,0)

  for (int i = 0; i < 6; ++i) {
    double diff = fabs(rec[i] - dU_t[i]);
    if (diff > tol) {
      printf("  FAIL [%s]: rec[%d] = %e, dU_t[%d] = %e, diff = %e\n",
        label, i, rec[i], i, dU_t[i], diff);
    }
    TEST_CHECK(diff < tol);
  }
}

// (No y-direction-specific reconstruction needed: the wave decomposition in
// the tetrad framework is universal -- the eigenvectors of the flat-Maxwell
// flux Jacobian along tetrad-x = face-normal are paper eq. 57 regardless of
// which physical face the kernel is operating on. Both x- and y-face kernels
// now use those same eigenvectors. check_x_reconstruction is the universal
// test.)

// ---- Tests ----

// Test 1: M^T M = h_ij and M M^{-1} = I across representative metric points.
static void
test_tetrad_M_satisfies_metric_and_inverse(void)
{
  double tol = 1.0e-12;

  // 1a: Cartesian flat.
  printf("\n  [point] Cartesian flat\n");
  {
    double h_ij[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double M[3][3], Minv[3][3];

    build_M_x_face(/*h_rp=*/0.0, /*h_yy=*/1.0, /*h_zz=*/1.0, /*Jc=*/1.0, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "x-face flat");
    check_M_inverse_correct(M, Minv, tol, "x-face flat");

    build_M_y_face(/*h_rr=*/1.0, /*h_rp=*/0.0, /*h_tt=*/1.0, /*h_pp=*/1.0, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "y-face flat");
    check_M_inverse_correct(M, Minv, tol, "y-face flat");
  }

  // 1b: Diagonal stretched.
  printf("  [point] diagonal stretched (a=2, b=3, c=4)\n");
  {
    double a = 2.0, b = 3.0, c = 4.0;
    double h_ij[3][3] = {{a*a,0,0},{0,b*b,0},{0,0,c*c}};
    double M[3][3], Minv[3][3];

    build_M_x_face(/*h_rp=*/0.0, b*b, c*c, a*b*c, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "x-face diagonal");
    check_M_inverse_correct(M, Minv, tol, "x-face diagonal");

    build_M_y_face(a*a, /*h_rp=*/0.0, b*b, c*c, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "y-face diagonal");
    check_M_inverse_correct(M, Minv, tol, "y-face diagonal");
  }

  // 1c: Spherical Kerr-Schild interior (r=1.5, sin^2 theta = 0.5, M_BH=1, a=0.95).
  printf("  [point] Kerr-Schild interior (r=1.5, sin^2 theta=0.5)\n");
  {
    double rho_sq = 1.5*1.5 + 0.95*0.95 * 0.5;
    double f = 2.0 * 1.5 / rho_sq;
    double h_rr = 1.0 + f;
    double h_tt = rho_sq;
    double h_pp = 0.5 * (rho_sq + 0.95*0.95 * (1.0 + f) * 0.5);
    double h_rp = -0.95 * (1.0 + f) * 0.5;
    double det_h = h_tt * (h_rr * h_pp - h_rp * h_rp);
    double Jc = sqrt(det_h);

    double h_ij[3][3] = {{h_rr, 0, h_rp}, {0, h_tt, 0}, {h_rp, 0, h_pp}};
    double M[3][3], Minv[3][3];

    build_M_x_face(h_rp, h_tt, h_pp, Jc, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "x-face Kerr-Schild");
    check_M_inverse_correct(M, Minv, tol, "x-face Kerr-Schild");

    build_M_y_face(h_rr, h_rp, h_tt, h_pp, M, Minv);
    check_M_satisfies_metric(M, h_ij, tol, "y-face Kerr-Schild");
    check_M_inverse_correct(M, Minv, tol, "y-face Kerr-Schild");
  }

  // 1d: Near-pole Kerr-Schild (sin theta = 0.1).
  printf("  [point] Kerr-Schild near-pole (r=1.5, sin theta = 0.1)\n");
  {
    double sin_th = 0.1;
    double cos_sq = 1.0 - sin_th*sin_th;
    double rho_sq = 1.5*1.5 + 0.95*0.95 * cos_sq;
    double f = 2.0 * 1.5 / rho_sq;
    double h_rr = 1.0 + f;
    double h_tt = rho_sq;
    double h_pp = sin_th*sin_th * (rho_sq + 0.95*0.95 * (1.0 + f) * sin_th*sin_th);
    double h_rp = -0.95 * (1.0 + f) * sin_th*sin_th;
    double det_h = h_tt * (h_rr * h_pp - h_rp * h_rp);
    double Jc = sqrt(det_h);

    double h_ij[3][3] = {{h_rr, 0, h_rp}, {0, h_tt, 0}, {h_rp, 0, h_pp}};
    double M[3][3], Minv[3][3];
    // h_pp ~ 0.01 here, so use a slightly looser tolerance.
    double pole_tol = 1.0e-11;

    build_M_x_face(h_rp, h_tt, h_pp, Jc, M, Minv);
    check_M_satisfies_metric(M, h_ij, pole_tol, "x-face near-pole");
    check_M_inverse_correct(M, Minv, pole_tol, "x-face near-pole");

    build_M_y_face(h_rr, h_rp, h_tt, h_pp, M, Minv);
    check_M_satisfies_metric(M, h_ij, pole_tol, "y-face near-pole");
    check_M_inverse_correct(M, Minv, pole_tol, "y-face near-pole");
  }
}

// Test 2: Eigenvectors orthogonal in flat tetrad.
static void
test_tetrad_eigenvectors_orthogonal(void)
{
  double tol = 1.0e-14;

  // Paper eq. 57 right eigenvectors for the flat-Maxwell flux Jacobian along
  // tetrad-x = face-normal direction. UNIVERSAL in tetrad framework -- both
  // x- and y-face kernels use exactly these eigenvectors (only the M used to
  // produce the tetrad-frame jumps differs between the two kernels).
  double r_universal[6][6] = {
    { 0, -1,  0, 0,  0, 1},  // r_minus_1
    { 0,  0,  1, 0,  1, 0},  // r_minus_2
    { 0,  0,  0, 1,  0, 0},  // r_zero_1
    { 1,  0,  0, 0,  0, 0},  // r_zero_2
    { 0,  1,  0, 0,  0, 1},  // r_plus_1
    { 0,  0, -1, 0,  1, 0},  // r_plus_2
  };
  double norm2[6] = { 2, 2, 1, 1, 2, 2 };
  check_eigenvectors_orthogonal(r_universal, norm2, tol, "tetrad eigenvectors");
}

// Test 3: Wave-decomposition reconstruction is complete.
static void
test_tetrad_wave_decomposition_complete(void)
{
  double tol = 1.0e-14;

  // Paper eq. 57 wave decomposition is universal in the tetrad framework:
  // the eigenvectors for the face-normal flat-Maxwell flux Jacobian apply to
  // both x- and y-face kernels.
  double dU1[6] = { 0.7, -1.3, 0.4, 1.9, -0.2, 0.5 };
  check_x_reconstruction(dU1, tol, "tetrad-mixed");
  double dU2[6] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  check_x_reconstruction(dU2, tol, "tetrad-Dx-only");
  double dU3[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  check_x_reconstruction(dU3, tol, "tetrad-Bz-only");
}

// Test 4: Tetrad-Roe kernel produces zero dissipation when U_L = U_R (the
// jump vanishes, so each wave strength is zero, so all dissipation is zero).
// Verifies by calling the kernel directly with crafted inputs.
static void
test_tetrad_kernel_zero_jump(void)
{
  double tol = 1.0e-14;

  // Cartesian flat metric at GL nodes (both nodes the same for this test).
  // Layout per face: 6 h_ij components, 2 nodes each (component-major).
  double h_ij[12]    = { 1, 1,  0, 0,  0, 0,   1, 1,  0, 0,   1, 1};   // h_xx, h_xy, h_xz, h_yy, h_yz, h_zz
  double det_h[2]    = { 1.0, 1.0 };
  double lapse[2]    = { 1.0, 1.0 };
  double shift[6]    = { 0, 0,  0, 0,  0, 0 };

  // To produce U_L = U_R AT THE FACE GL NODE, we must populate only the
  // cell-average (basis-0) mode. Non-zero slope modes (basis 1,2,3) evaluate
  // with opposite signs on the L cell (face at xi=+1) and R cell (face at
  // xi=-1), so identical modal coefficients give DIFFERENT node values.
  // With only basis-0 populated, the node value is just 0.5*q[0] on both
  // sides (matching coefficients), so the jump vanishes.
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  for (int comp = 0; comp < 6; ++comp) {
    double avg = 0.31415 * (double)(comp+1);
    field_con_l[comp*4 + 0] = avg;  // basis 0 (cell-average) only
    field_con_r[comp*4 + 0] = avg;
  }

  // Centered-flux input (precomputed; we only test that the kernel doesn't
  // add dissipation, so the centered piece is whatever flux_l == flux_r).
  double flux_l[12], flux_r[12];
  for (int k = 0; k < 12; ++k) flux_l[k] = flux_r[k] = 0.123 * (double)(k+1);

  double dxv[2] = { 1.0, 1.0 };
  double conf_flux_surf[32] = {0};

  // ---- Radial face ----
  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  // For equal L and R, jump = 0, so dissipation = 0 and the output equals
  // the centered flux (1/2 (F_L + F_R) = F_L). Check per component, per node.
  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = flux_l[comp*2 + n];   // == flux_r[comp*2+n]
      double got = conf_flux_surf[n + comp*2];
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [x-zero-jump]: comp %d node %d: got %e expected %e diff %e\n",
          comp, n, got, expected, diff);
      }
      TEST_CHECK(diff < tol);
    }
  }
  // Auxiliary slots (12..15) zeroed by the kernel.
  for (int k = 12; k < 16; ++k) {
    TEST_CHECK(conf_flux_surf[k] == 0.0);
  }

  // ---- Theta face ----
  for (int k = 0; k < 32; ++k) conf_flux_surf[k] = 0.0;
  lax_flux_y_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = flux_l[comp*2 + n];
      double got = conf_flux_surf[16 + n + comp*2];   // dir=1 slice offset
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [y-zero-jump]: comp %d node %d: got %e expected %e diff %e\n",
          comp, n, got, expected, diff);
      }
      TEST_CHECK(diff < tol);
    }
  }

  // ---- Theta face at theta_pole ----
  // theta_pole=1 should also produce centered-flux only (zero dissipation).
  for (int k = 0; k < 32; ++k) conf_flux_surf[k] = 0.0;
  lax_flux_y_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/1,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = flux_l[comp*2 + n];
      double got = conf_flux_surf[16 + n + comp*2];
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [y-pole-zero-jump]: comp %d node %d: got %e expected %e diff %e\n",
          comp, n, got, expected, diff);
      }
      TEST_CHECK(diff < tol);
    }
  }
  // The theta_pole=1 branch should also produce centered-flux output even
  // when U_L != U_R. Sanity-check this with a non-zero jump at the pole:
  // the output should not depend on the jump (only on the centered flux).
  // Use a basis-0 jump (so the node-value jump is well-defined).
  for (int comp = 0; comp < 6; ++comp) {
    field_con_r[comp*4 + 0] = field_con_l[comp*4 + 0] + 0.7;
  }
  for (int k = 0; k < 32; ++k) conf_flux_surf[k] = 0.0;
  lax_flux_y_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/1,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);
  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = 0.5*(flux_l[comp*2+n] + flux_r[comp*2+n]);
      double got = conf_flux_surf[16 + n + comp*2];
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [y-pole-nonzero-jump]: comp %d node %d: got %e expected %e\n",
          comp, n, got, expected);
      }
      TEST_CHECK(diff < tol);
    }
  }
}

// ---- Additional tests for the co-moving frame and back-transformation ----

// Apply the flat-Maxwell flux Jacobian for the x-direction to a 6-vector v.
// J^x is the Jacobian of F^x[D, B] = (0, B_z, -B_y, 0, -D_z, D_y), so
//   J^x[1][5] = +1   (d F_Dy / d B_z)
//   J^x[2][4] = -1   (d F_Dz / d B_y)
//   J^x[4][2] = -1   (d F_By / d D_z)
//   J^x[5][1] = +1   (d F_Bz / d D_y)
// All other entries are zero.
static void
apply_J_x_flat(double v[6], double Jv[6])
{
  Jv[0] =  0.0;
  Jv[1] = +v[5];
  Jv[2] = -v[4];
  Jv[3] =  0.0;
  Jv[4] = -v[2];
  Jv[5] = +v[1];
}

// (apply_J_y_flat removed: the kernel uses paper eq. 57 universally in tetrad
// basis, which is for the FACE-NORMAL flat-Maxwell flux Jacobian. That
// Jacobian is the same form regardless of which physical face we are on.
// apply_J_x_flat above covers all the kernel's eigenstructure.)

// Test 5: each r_p is a proper eigenvector with the claimed eigenvalue.
// Verifies J · r_p = lambda_p · r_p directly for all 6 eigenvectors of each
// direction. This catches any inconsistency between the eigenvectors and the
// eigenvalues we use for dissipation.
static void
test_tetrad_eigenvalue_eigenvector_consistency(void)
{
  double tol = 1.0e-14;

  printf("\n  [direction] x (radial face)\n");
  // Right eigenvectors (paper eq. 57 for x-direction) and claimed eigenvalues.
  double r_x[6][6] = {
    { 0, -1,  0, 0,  0, 1},  // lambda = -1
    { 0,  0,  1, 0,  1, 0},  // lambda = -1
    { 0,  0,  0, 1,  0, 0},  // lambda =  0
    { 1,  0,  0, 0,  0, 0},  // lambda =  0
    { 0,  1,  0, 0,  0, 1},  // lambda = +1
    { 0,  0, -1, 0,  1, 0},  // lambda = +1
  };
  double lambda_x[6] = { -1, -1, 0, 0, +1, +1 };
  for (int p = 0; p < 6; ++p) {
    double Jr[6];
    apply_J_x_flat(r_x[p], Jr);
    for (int i = 0; i < 6; ++i) {
      double expected = lambda_x[p] * r_x[p][i];
      double diff = fabs(Jr[i] - expected);
      if (diff > tol) {
        printf("  FAIL [x]: r_%d, J*r[%d] = %e, lambda*r[%d] = %e\n",
          p, i, Jr[i], i, expected);
      }
      TEST_CHECK(diff < tol);
    }
  }

  // (y-face tests use the SAME paper eq. 57 eigenstructure -- the wave
  // decomposition is for the face-normal direction, which is "tetrad x" in
  // both kernels, and the flat-Maxwell flux Jacobian along that axis is
  // identical regardless of which coord direction is the face normal.)
}

// Test 6: Roe flux differencing identity: sum_p lambda_p alpha_p r_p = J * dU_t.
// This is the discrete analog of the moments-test "f_R - f_L = amdq + apdq"
// check, expressed in tetrad-basis quantities. Catches sign errors,
// eigenpair mismatches, and decomposition errors in one identity.
static void
test_tetrad_flux_differencing(void)
{
  double tol = 1.0e-13;

  printf("\n  [direction] x (radial face)\n");
  {
    // Pick an arbitrary jump (any non-trivial vector).
    double dU_t[6] = { 0.31, -0.7, 1.2, -0.5, 0.95, 0.4 };

    // Wave strengths (paper eq. 57).
    double a_m1 = (-dU_t[1] + dU_t[5]) * 0.5;
    double a_m2 = ( dU_t[2] + dU_t[4]) * 0.5;
    double a_01 =   dU_t[3];
    double a_02 =   dU_t[0];
    double a_p1 = ( dU_t[1] + dU_t[5]) * 0.5;
    double a_p2 = (-dU_t[2] + dU_t[4]) * 0.5;
    double lambda[6] = { -1, -1, 0, 0, +1, +1 };
    double alpha[6]  = { a_m1, a_m2, a_01, a_02, a_p1, a_p2 };
    double r[6][6] = {
      { 0, -1,  0, 0,  0, 1},
      { 0,  0,  1, 0,  1, 0},
      { 0,  0,  0, 1,  0, 0},
      { 1,  0,  0, 0,  0, 0},
      { 0,  1,  0, 0,  0, 1},
      { 0,  0, -1, 0,  1, 0},
    };

    // Sum lambda_p alpha_p r_p.
    double sum[6] = {0};
    for (int p = 0; p < 6; ++p)
      for (int i = 0; i < 6; ++i)
        sum[i] += lambda[p] * alpha[p] * r[p][i];

    // Compare to J^x dU_t.
    double JdU[6];
    apply_J_x_flat(dU_t, JdU);
    for (int i = 0; i < 6; ++i) {
      double diff = fabs(sum[i] - JdU[i]);
      if (diff > tol) {
        printf("  FAIL [x]: i=%d sum=%e J*dU=%e diff=%e\n",
          i, sum[i], JdU[i], diff);
      }
      TEST_CHECK(diff < tol);
    }
  }

  // (y-face flux differencing uses the SAME paper eq. 57 wave decomposition
  // -- the identity sum_p lambda_p alpha_p r_p = J * dU_t holds in tetrad
  // basis for the face-normal flux Jacobian, which is universal across both
  // kernels. The x-direction test above covers it.)
}

// Test 7: norm preservation under back-transform.
// The forward transform V^a = M^a_i V^i takes a coord vector to its tetrad
// components. The inverse V^i = M^{-1}^i_a V^a takes tetrad components back
// to coord. The fundamental property of the orthonormal tetrad is that the
// magnitude is preserved: delta_ab V^a V^b = h_ij V^i V^j.
//
// Equivalently, given any tetrad vector V^a, the back-transformed coord
// vector V^i must satisfy h_ij V^i V^j = delta_ab V^a V^b. This test
// catches errors in M^{-1} that don't show up in M^T M = h_ij or M M^{-1} = I.
static void
test_tetrad_norm_preservation(void)
{
  double tol = 1.0e-12;

  // Use the Kerr-Schild interior point (non-trivial off-diagonal h_rphi).
  double rho_sq = 1.5*1.5 + 0.95*0.95 * 0.5;
  double f = 2.0 * 1.5 / rho_sq;
  double h_rr = 1.0 + f;
  double h_tt = rho_sq;
  double h_pp = 0.5 * (rho_sq + 0.95*0.95 * (1.0 + f) * 0.5);
  double h_rp = -0.95 * (1.0 + f) * 0.5;
  double det_h = h_tt * (h_rr * h_pp - h_rp * h_rp);
  double Jc = sqrt(det_h);

  double h_ij[3][3] = {{h_rr, 0, h_rp}, {0, h_tt, 0}, {h_rp, 0, h_pp}};
  double M[3][3], Minv[3][3];

  // Several test tetrad vectors.
  double V_t[5][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {0.7, -1.3, 0.9},
    {-0.5, 0.0, 1.1},
  };

  printf("\n  [direction] x (radial face)\n");
  build_M_x_face(h_rp, h_tt, h_pp, Jc, M, Minv);
  for (int t = 0; t < 5; ++t) {
    // V^i = (M^{-1})^i_a V^a
    double V[3] = {0};
    for (int i = 0; i < 3; ++i)
      for (int a = 0; a < 3; ++a)
        V[i] += Minv[i][a] * V_t[t][a];

    // Tetrad norm: delta_ab V^a V^b
    double norm_t = 0.0;
    for (int a = 0; a < 3; ++a) norm_t += V_t[t][a] * V_t[t][a];

    // Coord norm: h_ij V^i V^j
    double norm_c = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        norm_c += h_ij[i][j] * V[i] * V[j];

    double diff = fabs(norm_t - norm_c);
    if (diff > tol) {
      printf("  FAIL [x t=%d]: norm_tetrad = %e, norm_coord = %e, diff = %e\n",
        t, norm_t, norm_c, diff);
    }
    TEST_CHECK(diff < tol);
  }

  printf("  [direction] y (theta face)\n");
  build_M_y_face(h_rr, h_rp, h_tt, h_pp, M, Minv);
  for (int t = 0; t < 5; ++t) {
    double V[3] = {0};
    for (int i = 0; i < 3; ++i)
      for (int a = 0; a < 3; ++a)
        V[i] += Minv[i][a] * V_t[t][a];

    double norm_t = 0.0;
    for (int a = 0; a < 3; ++a) norm_t += V_t[t][a] * V_t[t][a];

    double norm_c = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        norm_c += h_ij[i][j] * V[i] * V[j];

    double diff = fabs(norm_t - norm_c);
    if (diff > tol) {
      printf("  FAIL [y t=%d]: norm_tetrad = %e, norm_coord = %e, diff = %e\n",
        t, norm_t, norm_c, diff);
    }
    TEST_CHECK(diff < tol);
  }
}

// Test 8: kernel reduces to standard flat-Maxwell Roe in the Cartesian flat
// limit (h = delta, alpha = 1, beta = 0). This is the most direct end-to-end
// check that the entire pipeline (modal-to-nodal, M, wave decomposition,
// dissipation, back-transform M^{-1}, output writes) is plumbed correctly.
//
// In flat space with no shift, the analytic flat-Maxwell Roe dissipation per
// component is:
//   diss_Dx = 0,                     diss_Dy = (dU_y)/2 * |1| evaluated only on D_y itself,
// or more precisely:
//   F* = (1/2)(F_L + F_R) - (1/2) sum_p |lambda_p| alpha_p r_p
// with lambda_p = +-1, 0 and alpha_p the wave strengths from the basis-0 jump.
// For our basis-0-only inputs, the node value = 0.5*q[0], so dU_node = 0.5
// times the modal jump.
static void
test_tetrad_kernel_flat_limit(void)
{
  double tol = 1.0e-13;

  // Cartesian flat geometry at both face nodes.
  double h_ij[12]    = { 1, 1,  0, 0,  0, 0,   1, 1,  0, 0,   1, 1};
  double det_h[2]    = { 1.0, 1.0 };
  double lapse[2]    = { 1.0, 1.0 };
  double shift[6]    = { 0, 0,  0, 0,  0, 0 };
  double dxv[2] = { 1.0, 1.0 };

  // Pick basis-0-only L and R states with a known jump in each component.
  // Field layout: 6 components x 4 modes per cell.
  double dU0[6] = { 0.3, -0.7, 1.2, -0.5, 0.95, 0.4 };  // jump in cell-avg per component
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  for (int comp = 0; comp < 6; ++comp) {
    field_con_l[comp*4 + 0] = 0.0;
    field_con_r[comp*4 + 0] = dU0[comp];   // basis-0 modes give jump dU0 directly at node
  }
  // At a face GL node, node value = 0.5 * q[0] (basis-0 normalization).
  // So the L node value = 0, R node value = 0.5*dU0[comp]. Jump = 0.5*dU0[comp].

  // Centered flux: arbitrary (kernel passes it through unchanged).
  double flux_l[12], flux_r[12];
  for (int k = 0; k < 12; ++k) flux_l[k] = flux_r[k] = 0.0;

  double conf_flux_surf[32] = {0};

  // ---- Radial face (x-direction Roe) ----
  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  // Expected dissipation per component (flat Maxwell, x-direction):
  //   diss_Dx = 0  (D_x rest mode)
  //   diss_Dy = (1/2) * dD_y_node    (= (1/2) * (1/2)*dU0[1] = dU0[1]/4)
  //   diss_Dz = (1/2) * dD_z_node
  //   diss_Bx = 0
  //   diss_By = (1/2) * dB_y_node
  //   diss_Bz = (1/2) * dB_z_node
  // dU_node[comp] = 0.5 * dU0[comp]. Output = -dissipation (since centered flux is 0).
  double expected_diss[6] = {
    0.0,
    0.5 * 0.5 * dU0[1],
    0.5 * 0.5 * dU0[2],
    0.0,
    0.5 * 0.5 * dU0[4],
    0.5 * 0.5 * dU0[5],
  };
  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = -expected_diss[comp];  // out = 0 (centered) - diss
      double got = conf_flux_surf[n + comp*2];
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [x flat]: comp %d node %d: got %e expected %e diff %e\n",
          comp, n, got, expected, diff);
      }
      TEST_CHECK(diff < tol);
    }
  }

  // ---- Theta face (y-direction Roe) ----
  for (int k = 0; k < 32; ++k) conf_flux_surf[k] = 0.0;
  lax_flux_y_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  // y-face tetrad-Roe in flat space: M is a permutation (tetrad x = coord
  // theta, tetrad y = coord r, tetrad z = coord phi). The kernel uses paper
  // eq. 57 for the face-normal (tetrad x) flat-Maxwell Jacobian, so the rest
  // modes in TETRAD basis are D^x_t = D^theta and B^x_t = B^theta. After
  // back-transforming with M^{-1} = M^T (permutation) those rest modes stay
  // in the theta direction in coord. So expected coord-basis dissipation:
  //   diss_Dr = (1/2)*dD^r_node,  diss_Dtheta = 0,  diss_Dphi = (1/2)*dD^phi_node,
  //   diss_Br = (1/2)*dB^r_node,  diss_Btheta = 0,  diss_Bphi = (1/2)*dB^phi_node.
  // (Coord components 0..5 = D^r, D^theta, D^phi, B^r, B^theta, B^phi.)
  double expected_diss_y[6] = {
    0.5 * 0.5 * dU0[0],   // D^r dissipation
    0.0,                   // D^theta = rest mode in tetrad x = face normal
    0.5 * 0.5 * dU0[2],   // D^phi dissipation
    0.5 * 0.5 * dU0[3],   // B^r dissipation
    0.0,                   // B^theta = rest mode
    0.5 * 0.5 * dU0[5],   // B^phi dissipation
  };
  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double expected = -expected_diss_y[comp];
      double got = conf_flux_surf[16 + n + comp*2];
      double diff = fabs(got - expected);
      if (diff > tol) {
        printf("  FAIL [y flat]: comp %d node %d: got %e expected %e diff %e\n",
          comp, n, got, expected, diff);
      }
      TEST_CHECK(diff < tol);
    }
  }
}

// Test 9: moving-interface eigenvalue shift.
// In the co-moving frame, the modified flux Jacobian J_total = J_flat -
// (beta_tilde^x/alpha) * I has eigenvalues shifted by -beta_tilde^x/alpha
// while eigenvectors are unchanged. Set up Cartesian metric (h = delta) but
// non-zero shift, so beta_tilde^x = beta^x in flat. Verify the kernel
// computes dissipation magnitudes |lambda_p - bxa| as expected.
//
// We test this by setting dU jumps that select individual waves and
// comparing the resulting dissipation magnitude.
static void
test_tetrad_kernel_moving_interface(void)
{
  double tol = 1.0e-12;

  // Cartesian flat metric, non-zero radial shift.
  double h_ij[12]    = { 1, 1,  0, 0,  0, 0,   1, 1,  0, 0,   1, 1};
  double det_h[2]    = { 1.0, 1.0 };
  double lapse[2]    = { 1.0, 1.0 };
  double beta_x      = 0.4;
  double shift[6]    = { beta_x, beta_x,  0, 0,  0, 0 };
  double dxv[2]      = { 1.0, 1.0 };

  // Select a pure r_minus_1 eigenmode (paper eq. 57 x-direction):
  // r_minus_1 = (0, -1, 0, 0, 0, 1). Set basis-0 jump dU = r_minus_1.
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  field_con_r[1*4 + 0] = -1.0;  // D_y basis-0 = -1
  field_con_r[5*4 + 0] = +1.0;  // B_z basis-0 = +1

  double flux_l[12], flux_r[12];
  for (int k = 0; k < 12; ++k) flux_l[k] = flux_r[k] = 0.0;

  double conf_flux_surf[32] = {0};
  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  // For a pure r_minus_1 jump, only the lambda = -1 wave is excited.
  // Its modified eigenvalue is -1 - bxa with bxa = beta_x/(alpha*sqrt(gamma_xx))
  // = beta_x (in flat space, alpha=1, gamma_xx=1).
  // |lambda_modified| = |-1 - bxa| = 1 + bxa (positive since bxa = 0.4).
  // Coord-time scale factor alpha*sqrt(gamma_xx) = 1, so |lambda|_coord = 1 + bxa.
  // Wave strength alpha_minus_1 = (-dD_y_node + dB_z_node)/2 = (-(-0.5) + 0.5)/2
  // = 0.5. (basis-0 mode -> 0.5*coefficient at node)
  // Dissipation in tetrad: 0.5 * |lambda| * alpha_p * r_p
  //   diss_t[D_y] = 0.5 * (1 + bxa) * 0.5 * (-1) = -(1+bxa)/4
  //   diss_t[B_z] = 0.5 * (1 + bxa) * 0.5 * (+1) = (1+bxa)/4
  // Back-transform with M = I (flat space): diss_coord = diss_tetrad.
  // Output = 0 (centered) - diss_coord:
  //   out[D_y] = (1+bxa)/4
  //   out[B_z] = -(1+bxa)/4
  // All other components: 0.
  double expected_Dy = +(1.0 + beta_x) * 0.25;
  double expected_Bz = -(1.0 + beta_x) * 0.25;
  for (int n = 0; n < 2; ++n) {
    double diff_Dx = fabs(conf_flux_surf[n + 0*2] - 0.0);
    double diff_Dy = fabs(conf_flux_surf[n + 1*2] - expected_Dy);
    double diff_Dz = fabs(conf_flux_surf[n + 2*2] - 0.0);
    double diff_Bx = fabs(conf_flux_surf[n + 3*2] - 0.0);
    double diff_By = fabs(conf_flux_surf[n + 4*2] - 0.0);
    double diff_Bz = fabs(conf_flux_surf[n + 5*2] - expected_Bz);
    if (diff_Dx > tol || diff_Dy > tol || diff_Dz > tol ||
        diff_Bx > tol || diff_By > tol || diff_Bz > tol) {
      printf("  FAIL [moving-interface, n=%d]:\n", n);
      printf("    Dx: got %e (diff %e)\n", conf_flux_surf[n+0*2], diff_Dx);
      printf("    Dy: got %e expected %e (diff %e)\n", conf_flux_surf[n+1*2], expected_Dy, diff_Dy);
      printf("    Dz: got %e (diff %e)\n", conf_flux_surf[n+2*2], diff_Dz);
      printf("    Bx: got %e (diff %e)\n", conf_flux_surf[n+3*2], diff_Bx);
      printf("    By: got %e (diff %e)\n", conf_flux_surf[n+4*2], diff_By);
      printf("    Bz: got %e expected %e (diff %e)\n", conf_flux_surf[n+5*2], expected_Bz, diff_Bz);
    }
    TEST_CHECK(diff_Dx < tol);
    TEST_CHECK(diff_Dy < tol);
    TEST_CHECK(diff_Dz < tol);
    TEST_CHECK(diff_Bx < tol);
    TEST_CHECK(diff_By < tol);
    TEST_CHECK(diff_Bz < tol);
  }

  // Now try a pure r_plus_1 jump (eigenvector (0, 1, 0, 0, 0, 1)). The wave
  // strength alpha_plus_1 = (dD_y + dB_z)/2 = (0.5 + 0.5)/2 = 0.5.
  // Modified eigenvalue: +1 - bxa = 1 - 0.4 = 0.6, |lambda| = 0.6.
  // Dissipation in tetrad: 0.5 * 0.6 * 0.5 * (eigvec components)
  //   diss_t[D_y] = 0.5 * 0.6 * 0.5 * (+1) = 0.15
  //   diss_t[B_z] = 0.5 * 0.6 * 0.5 * (+1) = 0.15
  // Output: -diss_coord = -0.15 each.
  for (int k = 0; k < 24; ++k) field_con_r[k] = 0.0;
  field_con_r[1*4 + 0] = +1.0;  // D_y basis-0 = +1
  field_con_r[5*4 + 0] = +1.0;  // B_z basis-0 = +1
  for (int k = 0; k < 32; ++k) conf_flux_surf[k] = 0.0;

  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  double expected_Dy_plus = -(1.0 - beta_x) * 0.25;
  double expected_Bz_plus = -(1.0 - beta_x) * 0.25;
  for (int n = 0; n < 2; ++n) {
    double diff_Dy = fabs(conf_flux_surf[n + 1*2] - expected_Dy_plus);
    double diff_Bz = fabs(conf_flux_surf[n + 5*2] - expected_Bz_plus);
    if (diff_Dy > tol || diff_Bz > tol) {
      printf("  FAIL [moving-interface plus, n=%d]:\n", n);
      printf("    Dy: got %e expected %e (diff %e)\n", conf_flux_surf[n+1*2], expected_Dy_plus, diff_Dy);
      printf("    Bz: got %e expected %e (diff %e)\n", conf_flux_surf[n+5*2], expected_Bz_plus, diff_Bz);
    }
    TEST_CHECK(diff_Dy < tol);
    TEST_CHECK(diff_Bz < tol);
  }
}

// ---- End-to-end coverage tests for the curved-metric pipeline ----
//
// Most of the kernel logic that is currently only IMPLICITLY tested by
// flat_limit and moving_interface (Cartesian metric where M is a permutation):
//
//   step 5.  forward transform of jumps:  dU_t = M . dU
//   step 7.  eigenvalue scaling: |lambda|_coord = |lambda_tetrad| * alpha * sqrt(gamma^xx)
//   step 10. back-transform of dissipation: diss = M^{-1} . diss_t
//
// The "decomposition consistency" tests below cover all three end-to-end in
// the curved metric (h_rphi != 0). Strategy:
//
//   1. Set up basis-0-only L and R fields with a known coord-basis jump dU0.
//   2. Run the kernel.
//   3. Decompose the INPUT jump dU0 (via M, then paper eq. 57) into wave
//      strengths alpha_p_in.
//   4. Decompose the OUTPUT (via M, then paper eq. 57) into wave strengths
//      alpha_p_out.
//   5. Verify alpha_p_out = -0.5 * |lambda_p|_coord * alpha_p_in for each p.
//
// If any of M, M^{-1}, the wave decomposition, the eigenvalue scaling, or
// the moving-interface correction is implemented wrong in the kernel, this
// relationship breaks. Tests K and L use the radial face (with and without
// shift); test M uses the theta face.

// Helper: paper eq. 57 wave-strength projection from a tetrad-frame jump.
static void
compute_wave_strengths(double dU_t[6], double a[6])
{
  a[0] = (-dU_t[1] + dU_t[5]) * 0.5;  // alpha_minus_1
  a[1] = ( dU_t[2] + dU_t[4]) * 0.5;  // alpha_minus_2
  a[2] =   dU_t[3];                    // alpha_zero_1
  a[3] =   dU_t[0];                    // alpha_zero_2
  a[4] = ( dU_t[1] + dU_t[5]) * 0.5;  // alpha_plus_1
  a[5] = (-dU_t[2] + dU_t[4]) * 0.5;  // alpha_plus_2
}

// Helper: forward-transform a (D, B) 6-vector to tetrad components using M.
static void
forward_transform_DB(double M[3][3], double dU[6], double dU_t[6])
{
  // D-block (components 0..2), then B-block (components 3..5).
  for (int a = 0; a < 3; ++a) {
    dU_t[a] = 0.0;
    dU_t[a+3] = 0.0;
    for (int i = 0; i < 3; ++i) {
      dU_t[a]   += M[a][i] * dU[i];     // M . D
      dU_t[a+3] += M[a][i] * dU[i+3];   // M . B
    }
  }
}

// Test J: modal-to-nodal coefficient sanity.
//
// Verifies the kernel-internal modal-to-nodal evaluation correctly produces
// the analytic node value from the basis-function expansion at the GL node,
// AND correctly distinguishes the L-side (xi=+1) from R-side (xi=-1) of the
// shared face.
//
// Setup: Cartesian flat metric (M = I). L AND R both populate the same
// basis-1 mode of D_y: q_l = q_r = (0, q1, 0, 0). At the SHARED face GL nodes:
//   L value (its own xi=+1, eta=+-1/sqrt(3)): 0.5*0 + (sqrt(3)/2)*q1 +- 0.5*0 +- (sqrt(3)/2)*0
//                                              = (sqrt(3)/2) * q1   [= +0.866*q1]
//   R value (its own xi=-1, eta=+-1/sqrt(3)): 0.5*0 - (sqrt(3)/2)*q1 +- 0.5*0 -+ (sqrt(3)/2)*0
//                                              = -(sqrt(3)/2) * q1  [= -0.866*q1]
//   Jump R - L = -sqrt(3) * q1
//
// If the kernel's modal-to-nodal coefficients are correct, the L cell's
// q1=0 and R cell's q1=0 give jump=0, but identical q1 in both gives a
// nontrivial jump because the cells evaluate at OPPOSITE local xi.
//
// Expected output (flat-Maxwell x-direction Roe with paper eq. 57):
//   alpha_minus_1 = (-dD_y_node + dB_z_node)/2 = (sqrt(3)*q1)/2
//   alpha_plus_1  = ( dD_y_node + dB_z_node)/2 = (-sqrt(3)*q1)/2
//   diss_Dy_t = 0.5 * (-|lambda_-|*alpha_m1 + |lambda_+|*alpha_p1)
//             = 0.5 * (-(sqrt(3)*q1)/2 + (-(sqrt(3)*q1)/2))
//             = -sqrt(3)*q1 / 2
//   out[D_y] = -diss_Dy_t = sqrt(3)*q1/2
//   diss_Bz_t = 0.5 * ( |lambda_-|*alpha_m1 + |lambda_+|*alpha_p1)
//             = 0.5 * ((sqrt(3)*q1)/2 - (sqrt(3)*q1)/2) = 0
//   out[B_z] = 0
static void
test_tetrad_modal_to_nodal_sanity(void)
{
  double tol = 1.0e-13;

  double h_ij[12]    = { 1, 1,  0, 0,  0, 0,   1, 1,  0, 0,   1, 1};
  double det_h[2]    = { 1.0, 1.0 };
  double lapse[2]    = { 1.0, 1.0 };
  double shift[6]    = { 0, 0,  0, 0,  0, 0 };
  double dxv[2]      = { 1.0, 1.0 };

  // Both cells: same basis-1 (linear-in-xi) mode of D_y. Shared face has
  // L_node = +(sqrt(3)/2)*q1, R_node = -(sqrt(3)/2)*q1, so jump = -sqrt(3)*q1.
  double q1 = 0.7;
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  field_con_l[1*4 + 1] = q1;   // D_y, basis 1 = +q1
  field_con_r[1*4 + 1] = q1;   // R cell same modal coefficient (different physical value at face)

  double flux_l[12] = {0};
  double flux_r[12] = {0};
  double conf_flux_surf[32] = {0};

  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  double sqrt3 = sqrt(3.0);
  double expected_out_Dy = sqrt3 * q1 * 0.5;

  for (int n = 0; n < 2; ++n) {
    double got_Dy = conf_flux_surf[n + 1*2];
    double diff = fabs(got_Dy - expected_out_Dy);
    if (diff > tol) {
      printf("  FAIL [modal-to-nodal n=%d]: out[D_y] = %e, expected %e (sqrt(3)*q1/2), diff %e\n",
        n, got_Dy, expected_out_Dy, diff);
    }
    TEST_CHECK(diff < tol);

    // Other components: should all be 0 (rest modes or no-jump components).
    int comps_other[5] = { 0, 2, 3, 4, 5 };
    for (int k = 0; k < 5; ++k) {
      int comp = comps_other[k];
      double got = conf_flux_surf[n + comp*2];
      if (fabs(got) > tol) {
        printf("  FAIL [modal-to-nodal n=%d comp=%d]: got %e, expected 0\n",
          n, comp, got);
      }
      TEST_CHECK(fabs(got) < tol);
    }
  }
  printf("\n  Modal-to-nodal coefficients consistent with basis-function values at GL nodes.\n");
}

// Test K: decomposition consistency in CURVED metric (radial face, no shift).
// Stresses M (with non-trivial h_rphi off-diagonal) and the back-transform
// M^{-1} in a setting where flat_limit can't catch errors.
static void
test_tetrad_decomp_consistency_x_curved_no_shift(void)
{
  double tol = 1.0e-12;

  // Kerr-Schild interior at r=1.5, sin^2 theta=0.5.
  double rho_sq = 1.5*1.5 + 0.95*0.95 * 0.5;
  double f = 2.0 * 1.5 / rho_sq;
  double h_rr = 1.0 + f;
  double h_tt = rho_sq;
  double h_pp = 0.5 * (rho_sq + 0.95*0.95 * (1.0 + f) * 0.5);
  double h_rp = -0.95 * (1.0 + f) * 0.5;
  double det_h_val = h_tt * (h_rr * h_pp - h_rp * h_rp);
  double Jc = sqrt(det_h_val);
  double alpha = 1.0 / sqrt(1.0 + f);
  // Force shift = 0 so we isolate the curved-metric back-transform from the
  // moving-interface piece (test L below adds the shift back in).

  double h_ij[12] = {
    h_rr, h_rr,  0,  0,    h_rp, h_rp,
    h_tt, h_tt,  0,  0,    h_pp, h_pp
  };
  double det_h_arr[2] = { Jc, Jc };
  double lapse[2]    = { alpha, alpha };
  double shift[6]    = { 0, 0,  0, 0,  0, 0 };
  double dxv[2]      = { 1.0, 1.0 };

  double dU0[6] = { 0.31, -0.7, 1.2, -0.5, 0.95, 0.4 };
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  for (int comp = 0; comp < 6; ++comp) field_con_r[comp*4 + 0] = dU0[comp];

  double flux_l[12] = {0};
  double flux_r[12] = {0};
  double conf_flux_surf[32] = {0};

  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h_arr, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  // Expected node jump = 0.5 * dU0 (basis-0 coefficient -> 0.5 * value at GL node).
  double dU_node[6];
  for (int i = 0; i < 6; ++i) dU_node[i] = 0.5 * dU0[i];

  // Build M, M^{-1} and project to tetrad.
  double M[3][3], Minv[3][3];
  build_M_x_face(h_rp, h_tt, h_pp, Jc, M, Minv);

  double dU_t[6];
  forward_transform_DB(M, dU_node, dU_t);

  double alpha_in[6];
  compute_wave_strengths(dU_t, alpha_in);

  // |lambda_p|_coord = |lambda_p_tetrad| * alpha * sqrt(gamma^xx); shift = 0,
  // so no moving-interface modification.
  double sqrt_gamma_xx = sqrt(h_tt * h_pp) / Jc;  // h_yz = 0 in our metric
  double scale_t2c = alpha * sqrt_gamma_xx;
  double abs_lambda_coord[6] = { scale_t2c, scale_t2c, 0, 0, scale_t2c, scale_t2c };

  // Decompose output and verify consistency.
  for (int n = 0; n < 2; ++n) {
    double output[6];
    for (int comp = 0; comp < 6; ++comp) output[comp] = conf_flux_surf[n + comp*2];

    double output_t[6];
    forward_transform_DB(M, output, output_t);

    double alpha_out[6];
    compute_wave_strengths(output_t, alpha_out);

    for (int p = 0; p < 6; ++p) {
      double expected = -0.5 * abs_lambda_coord[p] * alpha_in[p];
      double diff = fabs(alpha_out[p] - expected);
      if (diff > tol) {
        printf("  FAIL [n=%d wave p=%d]: alpha_out=%e expected=%e diff=%e\n",
          n, p, alpha_out[p], expected, diff);
        printf("    alpha_in[p]=%e, |lambda|_coord=%e\n", alpha_in[p], abs_lambda_coord[p]);
      }
      TEST_CHECK(diff < tol);
    }
  }
  printf("\n  Decomposition consistent in curved metric (radial, no shift).\n");
}

// Test L: decomposition consistency in CURVED metric (radial face, with shift).
// Adds the moving-interface piece on top of test K. Stresses the eigenvalue
// modification beta_tilde^x/alpha and the eigenvalue scaling alpha*sqrt(gamma).
static void
test_tetrad_decomp_consistency_x_curved_with_shift(void)
{
  double tol = 1.0e-12;

  double rho_sq = 1.5*1.5 + 0.95*0.95 * 0.5;
  double f = 2.0 * 1.5 / rho_sq;
  double h_rr = 1.0 + f;
  double h_tt = rho_sq;
  double h_pp = 0.5 * (rho_sq + 0.95*0.95 * (1.0 + f) * 0.5);
  double h_rp = -0.95 * (1.0 + f) * 0.5;
  double det_h_val = h_tt * (h_rr * h_pp - h_rp * h_rp);
  double Jc = sqrt(det_h_val);
  double alpha = 1.0 / sqrt(1.0 + f);
  // Real Kerr-Schild radial shift: beta^r = 2Mr/(rho^2 + 2Mr).
  double beta_r = 2.0 * 1.5 / (rho_sq + 2.0 * 1.5);

  double h_ij[12] = {
    h_rr, h_rr,  0,  0,    h_rp, h_rp,
    h_tt, h_tt,  0,  0,    h_pp, h_pp
  };
  double det_h_arr[2] = { Jc, Jc };
  double lapse[2]    = { alpha, alpha };
  double shift[6]    = { beta_r, beta_r,  0, 0,  0, 0 };
  double dxv[2]      = { 1.0, 1.0 };

  double dU0[6] = { 0.31, -0.7, 1.2, -0.5, 0.95, 0.4 };
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  for (int comp = 0; comp < 6; ++comp) field_con_r[comp*4 + 0] = dU0[comp];

  double flux_l[12] = {0};
  double flux_r[12] = {0};
  double conf_flux_surf[32] = {0};

  lax_flux_x_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h_arr, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  double dU_node[6];
  for (int i = 0; i < 6; ++i) dU_node[i] = 0.5 * dU0[i];

  double M[3][3], Minv[3][3];
  build_M_x_face(h_rp, h_tt, h_pp, Jc, M, Minv);

  double dU_t[6];
  forward_transform_DB(M, dU_node, dU_t);

  double alpha_in[6];
  compute_wave_strengths(dU_t, alpha_in);

  // beta_tilde^x/alpha = beta^x_coord / (alpha * sqrt(gamma^xx)).
  double sqrt_gamma_xx = sqrt(h_tt * h_pp) / Jc;
  double bxa = beta_r / (alpha * sqrt_gamma_xx);
  double scale_t2c = alpha * sqrt_gamma_xx;
  double abs_lambda_coord[6] = {
    fabs(-1.0 - bxa) * scale_t2c,    // r_minus_1
    fabs(-1.0 - bxa) * scale_t2c,    // r_minus_2
    fabs(    -bxa)   * scale_t2c,    // r_zero_1
    fabs(    -bxa)   * scale_t2c,    // r_zero_2
    fabs( 1.0 - bxa) * scale_t2c,    // r_plus_1
    fabs( 1.0 - bxa) * scale_t2c,    // r_plus_2
  };

  for (int n = 0; n < 2; ++n) {
    double output[6];
    for (int comp = 0; comp < 6; ++comp) output[comp] = conf_flux_surf[n + comp*2];

    double output_t[6];
    forward_transform_DB(M, output, output_t);

    double alpha_out[6];
    compute_wave_strengths(output_t, alpha_out);

    for (int p = 0; p < 6; ++p) {
      double expected = -0.5 * abs_lambda_coord[p] * alpha_in[p];
      double diff = fabs(alpha_out[p] - expected);
      if (diff > tol) {
        printf("  FAIL [n=%d wave p=%d]: alpha_out=%e expected=%e diff=%e\n",
          n, p, alpha_out[p], expected, diff);
        printf("    alpha_in[p]=%e, |lambda|_coord=%e, bxa=%e\n",
          alpha_in[p], abs_lambda_coord[p], bxa);
      }
      TEST_CHECK(diff < tol);
    }
  }
  printf("\n  Decomposition consistent in curved metric (radial, with shift bxa = %.4f).\n", bxa);
}

// Test M: decomposition consistency in CURVED metric (theta face).
// Stresses the y-face M (h_rphi mixes differently), and verifies the kernel
// uses paper eq. 57 universally (after the x-kernel-bug fix).
static void
test_tetrad_decomp_consistency_y_curved(void)
{
  double tol = 1.0e-12;

  double rho_sq = 1.5*1.5 + 0.95*0.95 * 0.5;
  double f = 2.0 * 1.5 / rho_sq;
  double h_rr = 1.0 + f;
  double h_tt = rho_sq;
  double h_pp = 0.5 * (rho_sq + 0.95*0.95 * (1.0 + f) * 0.5);
  double h_rp = -0.95 * (1.0 + f) * 0.5;
  double det_h_val = h_tt * (h_rr * h_pp - h_rp * h_rp);
  double Jc = sqrt(det_h_val);
  double alpha = 1.0 / sqrt(1.0 + f);
  // Kerr-Schild has beta^theta = 0 (the y-face moving interface velocity is
  // automatically zero in our metric).

  double h_ij[12] = {
    h_rr, h_rr,  0,  0,    h_rp, h_rp,
    h_tt, h_tt,  0,  0,    h_pp, h_pp
  };
  double det_h_arr[2] = { Jc, Jc };
  double lapse[2]    = { alpha, alpha };
  double shift[6]    = { 0, 0,  0, 0,  0, 0 };
  double dxv[2]      = { 1.0, 1.0 };

  double dU0[6] = { 0.31, -0.7, 1.2, -0.5, 0.95, 0.4 };
  double field_con_l[24] = {0};
  double field_con_r[24] = {0};
  for (int comp = 0; comp < 6; ++comp) field_con_r[comp*4 + 0] = dU0[comp];

  double flux_l[12] = {0};
  double flux_r[12] = {0};
  double conf_flux_surf[32] = {0};

  lax_flux_y_2x_ser_p1_tetrad_roe(dxv, /*theta_pole=*/0,
    lapse, shift, h_ij, det_h_arr, flux_l, flux_r,
    field_con_l, field_con_r, conf_flux_surf);

  double dU_node[6];
  for (int i = 0; i < 6; ++i) dU_node[i] = 0.5 * dU0[i];

  // Y-face M (different from x-face M).
  double M[3][3], Minv[3][3];
  build_M_y_face(h_rr, h_rp, h_tt, h_pp, M, Minv);

  double dU_t[6];
  forward_transform_DB(M, dU_node, dU_t);

  double alpha_in[6];
  compute_wave_strengths(dU_t, alpha_in);

  // y-face: |lambda|_coord = |lambda_tetrad| * alpha * D_theta / Jc, with
  // bxa = 0 since beta^theta = 0.
  double D_theta = sqrt(h_rr * h_pp - h_rp * h_rp);
  double scale_t2c = alpha * D_theta / Jc;
  double abs_lambda_coord[6] = { scale_t2c, scale_t2c, 0, 0, scale_t2c, scale_t2c };

  // Output for y-face goes to slot [16 + comp*2 + n].
  for (int n = 0; n < 2; ++n) {
    double output[6];
    for (int comp = 0; comp < 6; ++comp) output[comp] = conf_flux_surf[16 + n + comp*2];

    double output_t[6];
    forward_transform_DB(M, output, output_t);

    double alpha_out[6];
    compute_wave_strengths(output_t, alpha_out);

    for (int p = 0; p < 6; ++p) {
      double expected = -0.5 * abs_lambda_coord[p] * alpha_in[p];
      double diff = fabs(alpha_out[p] - expected);
      if (diff > tol) {
        printf("  FAIL [y n=%d wave p=%d]: alpha_out=%e expected=%e diff=%e\n",
          n, p, alpha_out[p], expected, diff);
        printf("    alpha_in[p]=%e, |lambda|_coord=%e\n", alpha_in[p], abs_lambda_coord[p]);
      }
      TEST_CHECK(diff < tol);
    }
  }
  printf("\n  Decomposition consistent in curved metric (theta face).\n");
}

// ---- Tetrad-Roe characteristic outflow kernel tests ----

// Test N: outflow kernel zero-output for a state with NO outgoing
// characteristic content. For lower-edge outflow at Cartesian flat with no
// shift, eigenvalues are -1, 0, +1. Outgoing = negative. The wave r_plus_1
// (eigenvalue +1) is incoming, so its amplitude should be zeroed.
// Construct: skin state in basis-0 such that ONLY r_plus_1 is excited.
// Expected: outflow kernel output = 0 (only incoming wave -> nothing keeps).
static void
test_outflow_kernel_only_incoming(void)
{
  double tol = 1.0e-14;

  // Cartesian flat, no shift.
  double h_ij[12] = { 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1 };
  double det_h[2] = { 1.0, 1.0 };
  double lapse[2] = { 1.0, 1.0 };
  double shift[6] = { 0, 0, 0, 0, 0, 0 };
  double dxv[2] = { 1.0, 1.0 };

  // r_plus_1 in tetrad = (0, 1, 0, 0, 0, 1) means D^y = 1, B^z = 1.
  // In flat M = I, this corresponds to D^y_coord = 1, B^z_coord = 1.
  // Skin field in basis 0 only.
  double field_con_skin[24] = {0};
  field_con_skin[1*4 + 0] = +1.0;  // D_y basis 0 = +1
  field_con_skin[5*4 + 0] = +1.0;  // B_z basis 0 = +1
  // At face GL node, node value = 0.5 * basis_0_coeff = 0.5 each. So
  // U_node = (0, 0.5, 0, 0, 0, 0.5). This excites r_plus_1 only.

  double conf_flux_surf[32] = {0};
  outflow_flux_x_2x_ser_p1_tetrad_roe(dxv, /*edge=*/0, lapse, shift, h_ij,
    det_h, field_con_skin, conf_flux_surf);

  // For edge=0 lower outflow: outgoing = negative. r_plus_1 has lambda_+ = +1
  // (positive) -> incoming -> zeroed. So no flux should be written.
  for (int comp = 0; comp < 6; ++comp) {
    for (int n = 0; n < 2; ++n) {
      double got = conf_flux_surf[n + comp*2];
      if (fabs(got) > tol) {
        printf("  FAIL [edge=0 r_plus_1]: comp %d node %d: got %e (should be 0)\n",
          comp, n, got);
      }
      TEST_CHECK(fabs(got) < tol);
    }
  }
  printf("\n  Lower-edge outflow correctly zeros pure-incoming (r_plus_1) flux.\n");
}

// Test O: outflow kernel non-zero output for a state with PURE outgoing
// content. Lower edge, r_minus_1 has eigenvalue -1 (outgoing). Skin state
// excites only r_minus_1; outflow kernel should give the full A^outgoing . U
// contribution.
static void
test_outflow_kernel_only_outgoing(void)
{
  double tol = 1.0e-13;

  double h_ij[12] = { 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1 };
  double det_h[2] = { 1.0, 1.0 };
  double lapse[2] = { 1.0, 1.0 };
  double shift[6] = { 0, 0, 0, 0, 0, 0 };
  double dxv[2] = { 1.0, 1.0 };

  // r_minus_1 in tetrad = (0, -1, 0, 0, 0, 1) means D^y = -1, B^z = +1.
  double field_con_skin[24] = {0};
  field_con_skin[1*4 + 0] = -1.0;  // D_y basis 0 = -1
  field_con_skin[5*4 + 0] = +1.0;  // B_z basis 0 = +1
  // U_node = (0, -0.5, 0, 0, 0, 0.5). alpha_minus_1 = (-(-0.5)+0.5)/2 = 0.5.

  double conf_flux_surf[32] = {0};
  outflow_flux_x_2x_ser_p1_tetrad_roe(dxv, /*edge=*/0, lapse, shift, h_ij,
    det_h, field_con_skin, conf_flux_surf);

  // Expected output:
  //   F = lambda_- * alpha_- * r_minus_1 (with lambda_-_coord = -1 in flat)
  //   F^a_tetrad = -1 * 0.5 * (0, -1, 0, 0, 0, 1) = (0, +0.5, 0, 0, 0, -0.5)
  // In flat M=I, F^i_coord = F^a_tetrad.
  double expected_Dy = +0.5;
  double expected_Bz = -0.5;

  for (int n = 0; n < 2; ++n) {
    double diff_Dy = fabs(conf_flux_surf[n + 1*2] - expected_Dy);
    double diff_Bz = fabs(conf_flux_surf[n + 5*2] - expected_Bz);
    if (diff_Dy > tol || diff_Bz > tol) {
      printf("  FAIL [r_minus_1 n=%d]: D_y got %e exp %e, B_z got %e exp %e\n",
        n, conf_flux_surf[n+1*2], expected_Dy, conf_flux_surf[n+5*2], expected_Bz);
    }
    TEST_CHECK(diff_Dy < tol);
    TEST_CHECK(diff_Bz < tol);
    // Other components should be 0.
    for (int comp = 0; comp < 6; ++comp) {
      if (comp == 1 || comp == 5) continue;
      double got = conf_flux_surf[n + comp*2];
      if (fabs(got) > tol) {
        printf("  FAIL [r_minus_1 n=%d comp=%d]: got %e (should be 0)\n",
          n, comp, got);
      }
      TEST_CHECK(fabs(got) < tol);
    }
  }
  printf("\n  Lower-edge outflow correctly applies A^outgoing.U_skin for r_minus_1.\n");
}

// Test P: upper-edge analog. For edge=1, outgoing = positive eigenvalue.
// Use r_plus_1 as the pure outgoing wave there.
static void
test_outflow_kernel_upper_edge(void)
{
  double tol = 1.0e-13;

  double h_ij[12] = { 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1 };
  double det_h[2] = { 1.0, 1.0 };
  double lapse[2] = { 1.0, 1.0 };
  double shift[6] = { 0, 0, 0, 0, 0, 0 };
  double dxv[2] = { 1.0, 1.0 };

  // r_plus_1 = (0, 1, 0, 0, 0, 1). At upper edge, this is outgoing.
  double field_con_skin[24] = {0};
  field_con_skin[1*4 + 0] = +1.0;  // D_y basis 0
  field_con_skin[5*4 + 0] = +1.0;  // B_z basis 0

  double conf_flux_surf[32] = {0};
  outflow_flux_x_2x_ser_p1_tetrad_roe(dxv, /*edge=*/1, lapse, shift, h_ij,
    det_h, field_con_skin, conf_flux_surf);

  // alpha_plus_1 = ( 0.5 + 0.5)/2 = 0.5
  // F^a = +1 * 0.5 * (0, 1, 0, 0, 0, 1) = (0, 0.5, 0, 0, 0, 0.5)
  double expected_Dy = +0.5;
  double expected_Bz = +0.5;

  for (int n = 0; n < 2; ++n) {
    double diff_Dy = fabs(conf_flux_surf[n + 1*2] - expected_Dy);
    double diff_Bz = fabs(conf_flux_surf[n + 5*2] - expected_Bz);
    if (diff_Dy > tol || diff_Bz > tol) {
      printf("  FAIL [upper r_plus_1 n=%d]: D_y got %e exp %e, B_z got %e exp %e\n",
        n, conf_flux_surf[n+1*2], expected_Dy, conf_flux_surf[n+5*2], expected_Bz);
    }
    TEST_CHECK(diff_Dy < tol);
    TEST_CHECK(diff_Bz < tol);
  }
  printf("\n  Upper-edge outflow correctly applies A^outgoing.U_skin for r_plus_1.\n");
}

TEST_LIST = {
  { "M_satisfies_metric_and_inverse", test_tetrad_M_satisfies_metric_and_inverse },
  { "eigenvectors_orthogonal",        test_tetrad_eigenvectors_orthogonal       },
  { "wave_decomposition_complete",    test_tetrad_wave_decomposition_complete   },
  { "kernel_zero_jump",               test_tetrad_kernel_zero_jump              },
  { "eigenvalue_consistency",         test_tetrad_eigenvalue_eigenvector_consistency },
  { "flux_differencing",              test_tetrad_flux_differencing             },
  { "norm_preservation",              test_tetrad_norm_preservation             },
  { "kernel_flat_limit",              test_tetrad_kernel_flat_limit             },
  { "kernel_moving_interface",        test_tetrad_kernel_moving_interface       },
  { "modal_to_nodal_sanity",          test_tetrad_modal_to_nodal_sanity         },
  { "decomp_x_curved_no_shift",       test_tetrad_decomp_consistency_x_curved_no_shift },
  { "decomp_x_curved_with_shift",     test_tetrad_decomp_consistency_x_curved_with_shift },
  { "decomp_y_curved",                test_tetrad_decomp_consistency_y_curved   },
  { "outflow_only_incoming",          test_outflow_kernel_only_incoming         },
  { "outflow_only_outgoing",          test_outflow_kernel_only_outgoing         },
  { "outflow_upper_edge",             test_outflow_kernel_upper_edge            },
  { NULL, NULL },
};
