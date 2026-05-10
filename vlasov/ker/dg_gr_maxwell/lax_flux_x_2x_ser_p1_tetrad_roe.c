#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Tetrad-frame Roe flux on the radial (x-direction) face, 2D ser p=1.
//
// Implements the per-quadrature-point version of the tetrad-first algorithm
// of Gorard et al. (2025). At each face Gauss-Legendre node we:
//
//   1. Build M^a_i, the orthonormal-tetrad change-of-basis matrix at that
//      physical point. For our spherical Kerr-Schild metric (h_rtheta =
//      h_thetaphi = 0) the closed form is
//          M = [[ 1/sqrt(gamma^xx), 0,            0           ],
//               [ 0,                sqrt(h_yy),   0           ],
//               [ h_xz/sqrt(h_zz),  0,            sqrt(h_zz)  ]],
//      where gamma^xx = (h_yy h_zz - h_yz^2)/det h. Verified to satisfy
//      M^T M = h_ij (analytical check; this carries flat tetrad metric back
//      to curved spatial metric in the inverse direction).
//   2. Transform jumps to tetrad: Delta tilde U^a = M^a_i Delta U^i, applied
//      independently to the D and B 3-vectors.
//   3. Decompose into 6 flat-Maxwell eigenwaves using paper eq. 57 right
//      eigenvectors. Eigenvectors are mutually orthogonal, so the wave
//      strengths are simple inner products (no matrix inversion).
//      Eigenvectors (in tetrad components ordered (D^x, D^y, D^z, B^x, B^y, B^z)):
//          r_minus_1 = (0, -1, 0, 0,  0, 1),     lambda = -1
//          r_minus_2 = (0,  0, 1, 0,  1, 0),     lambda = -1
//          r_zero_1  = (0,  0, 0, 1,  0, 0),     lambda =  0
//          r_zero_2  = (1,  0, 0, 0,  0, 0),     lambda =  0
//          r_plus_1  = (0,  1, 0, 0,  0, 1),     lambda = +1
//          r_plus_2  = (0,  0,-1, 0,  1, 0),     lambda = +1
//      Wave strengths (using r_p / |r_p|^2 as left eigenvectors):
//          alpha_minus_1 = (-dDy_t + dBz_t)/2
//          alpha_minus_2 = ( dDz_t + dBy_t)/2
//          alpha_zero_1  =   dBx_t
//          alpha_zero_2  =   dDx_t
//          alpha_plus_1  = ( dDy_t + dBz_t)/2
//          alpha_plus_2  = (-dDz_t + dBy_t)/2
//   4. The face hypersurface drifts in tetrad coords with velocity
//      tilde beta^x / alpha = beta^x / (alpha sqrt(gamma^xx)). The modified
//      eigenvalues for the moving-interface Riemann problem are
//      (eigenvalues of flat Maxwell) - tilde beta^x/alpha. Eigenvectors are
//      unchanged. So the dissipation magnitudes are
//          |lambda_minus^mod| = |-1 - tilde beta^x/alpha|
//          |lambda_zero^mod|  = |    -tilde beta^x/alpha|
//          |lambda_plus^mod|  = | 1 - tilde beta^x/alpha|.
//   5. Roe dissipation in tetrad: tilde diss^a = (1/2) sum_p |lambda_p^mod|
//      alpha_p r_p^a. Back-transform with M^{-1} to coord basis:
//          M^{-1} = [[ sqrt(gamma^xx),                 0,             0          ],
//                    [ 0,                              1/sqrt(h_yy),  0          ],
//                    [-h_xz sqrt(gamma^xx)/h_zz,       0,             1/sqrt(h_zz)]].
//   6. Output: F^i* = (1/2)(F^i_L + F^i_R) - diss^i, with the centered flux
//      taken from the existing precomputed flux_l, flux_r (which are in coord
//      basis, already including the curved Maxwell shift-advection terms).
//      M F^coord exactly equals tilde F^x_total = tilde F^x_flat - (tilde beta^x/alpha)
//      tilde U, so M^{-1} M (1/2)(F_L + F_R) recovers (1/2)(F^i_L + F^i_R)
//      identically -- no correction to the centered flux is needed.
//
// Aux components (phi, psi divergence cleaning) are decoupled from the EM
// system and zeroed in the dir=0 slot, matching the existing convention.
//
// CFL: returns 1.5 * (2/dx) * max wave speed across the 2 face nodes,
// matching the existing curved-LLF kernel's return convention.
GKYL_CU_DH double
lax_flux_x_2x_ser_p1_tetrad_roe(const double *dxv, const int theta_pole,
  const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *det_h,
  const double *flux_l, const double *flux_r,
  const double *field_con_l, const double *field_con_r,
  double* GKYL_RESTRICT conf_flux_surf)
{
  double dx10 = 2.0 / dxv[0];

  // Conserved-variable pointers (J*D and J*B; layout: 6 components x 4 modes).
  const double *JDx_l = &field_con_l[0];
  const double *JDy_l = &field_con_l[4];
  const double *JDz_l = &field_con_l[8];
  const double *JBx_l = &field_con_l[12];
  const double *JBy_l = &field_con_l[16];
  const double *JBz_l = &field_con_l[20];

  const double *JDx_r = &field_con_r[0];
  const double *JDy_r = &field_con_r[4];
  const double *JDz_r = &field_con_r[8];
  const double *JBx_r = &field_con_r[12];
  const double *JBy_r = &field_con_r[16];
  const double *JBz_r = &field_con_r[20];

  // h_ij and shift slicing per face GL node (component-major, 2 nodes per
  // component). Same layout as the existing curved-LLF kernel.
  const double *h_xx_n = &h_ij_nodal[0];
  const double *h_xy_n = &h_ij_nodal[2];
  const double *h_xz_n = &h_ij_nodal[4];
  const double *h_yy_n = &h_ij_nodal[6];
  const double *h_yz_n = &h_ij_nodal[8];
  const double *h_zz_n = &h_ij_nodal[10];
  const double *shift_x_n = &shift_nodal[0];

  double *out = conf_flux_surf;
  double cflrate = 0.0;

  for (int n = 0; n < 2; ++n) {
    // ---- Modal-to-nodal eval at this face GL node ----
    // Coefficients copied verbatim from dg_gr_maxwell_alpha_quad_x_2x_ser_p1.c
    // so the modal-to-nodal evaluation matches what the rest of the pipeline uses.
    double JDx_L, JDy_L, JDz_L, JBx_L, JBy_L, JBz_L;
    double JDx_R, JDy_R, JDz_R, JBx_R, JBy_R, JBz_R;
    if (n == 0) {
      JDx_L = -(0.8660254037844386)*JDx_l[3] - 0.5*JDx_l[2] + 0.8660254037844386*JDx_l[1] + 0.5*JDx_l[0];
      JDy_L = -(0.8660254037844386)*JDy_l[3] - 0.5*JDy_l[2] + 0.8660254037844386*JDy_l[1] + 0.5*JDy_l[0];
      JDz_L = -(0.8660254037844386)*JDz_l[3] - 0.5*JDz_l[2] + 0.8660254037844386*JDz_l[1] + 0.5*JDz_l[0];
      JBx_L = -(0.8660254037844386)*JBx_l[3] - 0.5*JBx_l[2] + 0.8660254037844386*JBx_l[1] + 0.5*JBx_l[0];
      JBy_L = -(0.8660254037844386)*JBy_l[3] - 0.5*JBy_l[2] + 0.8660254037844386*JBy_l[1] + 0.5*JBy_l[0];
      JBz_L = -(0.8660254037844386)*JBz_l[3] - 0.5*JBz_l[2] + 0.8660254037844386*JBz_l[1] + 0.5*JBz_l[0];

      JDx_R =   0.8660254037844386*JDx_r[3] - 0.5*JDx_r[2] - 0.8660254037844386*JDx_r[1] + 0.5*JDx_r[0];
      JDy_R =   0.8660254037844386*JDy_r[3] - 0.5*JDy_r[2] - 0.8660254037844386*JDy_r[1] + 0.5*JDy_r[0];
      JDz_R =   0.8660254037844386*JDz_r[3] - 0.5*JDz_r[2] - 0.8660254037844386*JDz_r[1] + 0.5*JDz_r[0];
      JBx_R =   0.8660254037844386*JBx_r[3] - 0.5*JBx_r[2] - 0.8660254037844386*JBx_r[1] + 0.5*JBx_r[0];
      JBy_R =   0.8660254037844386*JBy_r[3] - 0.5*JBy_r[2] - 0.8660254037844386*JBy_r[1] + 0.5*JBy_r[0];
      JBz_R =   0.8660254037844386*JBz_r[3] - 0.5*JBz_r[2] - 0.8660254037844386*JBz_r[1] + 0.5*JBz_r[0];
    } else {
      JDx_L =   0.8660254037844386*JDx_l[3] + 0.5*JDx_l[2] + 0.8660254037844386*JDx_l[1] + 0.5*JDx_l[0];
      JDy_L =   0.8660254037844386*JDy_l[3] + 0.5*JDy_l[2] + 0.8660254037844386*JDy_l[1] + 0.5*JDy_l[0];
      JDz_L =   0.8660254037844386*JDz_l[3] + 0.5*JDz_l[2] + 0.8660254037844386*JDz_l[1] + 0.5*JDz_l[0];
      JBx_L =   0.8660254037844386*JBx_l[3] + 0.5*JBx_l[2] + 0.8660254037844386*JBx_l[1] + 0.5*JBx_l[0];
      JBy_L =   0.8660254037844386*JBy_l[3] + 0.5*JBy_l[2] + 0.8660254037844386*JBy_l[1] + 0.5*JBy_l[0];
      JBz_L =   0.8660254037844386*JBz_l[3] + 0.5*JBz_l[2] + 0.8660254037844386*JBz_l[1] + 0.5*JBz_l[0];

      JDx_R = -(0.8660254037844386)*JDx_r[3] + 0.5*JDx_r[2] - 0.8660254037844386*JDx_r[1] + 0.5*JDx_r[0];
      JDy_R = -(0.8660254037844386)*JDy_r[3] + 0.5*JDy_r[2] - 0.8660254037844386*JDy_r[1] + 0.5*JDy_r[0];
      JDz_R = -(0.8660254037844386)*JDz_r[3] + 0.5*JDz_r[2] - 0.8660254037844386*JDz_r[1] + 0.5*JDz_r[0];
      JBx_R = -(0.8660254037844386)*JBx_r[3] + 0.5*JBx_r[2] - 0.8660254037844386*JBx_r[1] + 0.5*JBx_r[0];
      JBy_R = -(0.8660254037844386)*JBy_r[3] + 0.5*JBy_r[2] - 0.8660254037844386*JBy_r[1] + 0.5*JBy_r[0];
      JBz_R = -(0.8660254037844386)*JBz_r[3] + 0.5*JBz_r[2] - 0.8660254037844386*JBz_r[1] + 0.5*JBz_r[0];
    }

    // Coord-basis jumps in (J D, J B). The J factor is a positive scalar at
    // this node and rides through the M-transformation transparently.
    double dDx = JDx_R - JDx_L;
    double dDy = JDy_R - JDy_L;
    double dDz = JDz_R - JDz_L;
    double dBx = JBx_R - JBx_L;
    double dBy = JBy_R - JBy_L;
    double dBz = JBz_R - JBz_L;

    // ---- Geometry at this node ----
    double h_yy = h_yy_n[n], h_yz = h_yz_n[n], h_zz = h_zz_n[n];
    double h_xz = h_xz_n[n];
    double Jc   = det_h[n];
    double alpha   = lapse_nodal[n];
    double beta_x  = shift_x_n[n];

    // Suppress -Wunused-variable; h_xx, h_xy not needed for our (x-face,
    // h_xy = h_yz = 0) tetrad construction.
    (void)h_xx_n;
    (void)h_xy_n;

    // ---- M^a_i and M^{-1} (closed form for h_xy = h_yz = 0) ----
    // sqrt(gamma^xx) = sqrt(h_yy h_zz - h_yz^2)/Jc -- same expression that
    // appears as the curved-Maxwell sqrt-term in alpha_quad_x.
    double sqrt_gamma_xx     = sqrt(h_yy*h_zz - h_yz*h_yz) / Jc;
    double inv_sqrt_gamma_xx = 1.0 / sqrt_gamma_xx;
    double sqrt_h_yy = sqrt(h_yy);
    double sqrt_h_zz = sqrt(h_zz);
    double inv_sqrt_h_yy = 1.0 / sqrt_h_yy;
    double inv_sqrt_h_zz = 1.0 / sqrt_h_zz;

    // M^a_i:
    double M_xx = inv_sqrt_gamma_xx;
    double M_yy = sqrt_h_yy;
    double M_zz = sqrt_h_zz;
    double M_zx = h_xz / sqrt_h_zz;

    // (M^{-1})^i_a:
    double Minv_xx = sqrt_gamma_xx;
    double Minv_yy = inv_sqrt_h_yy;
    double Minv_zx = -h_xz * sqrt_gamma_xx / h_zz;
    double Minv_zz = inv_sqrt_h_zz;

    // ---- Tetrad jumps: dU_t^a = M^a_i dU^i ----
    // Rows of M:
    //   row x: (M_xx, 0,    0   )
    //   row y: (0,    M_yy, 0   )
    //   row z: (M_zx, 0,    M_zz)
    double dDx_t = M_xx * dDx;
    double dDy_t = M_yy * dDy;
    double dDz_t = M_zx * dDx + M_zz * dDz;
    double dBx_t = M_xx * dBx;
    double dBy_t = M_yy * dBy;
    double dBz_t = M_zx * dBx + M_zz * dBz;

    // ---- Wave strengths (paper eq. 57; eigenvectors orthogonal so left
    // eigenvectors are r_p / |r_p|^2; the |r_p|^2 = 2 for r_minus, r_plus
    // and 1 for r_zero) ----
    double alpha_m1 = ( -dDy_t + dBz_t) * 0.5;  // r_minus_1
    double alpha_m2 = (  dDz_t + dBy_t) * 0.5;  // r_minus_2
    double alpha_01 =    dBx_t;                  // r_zero_1
    double alpha_02 =    dDx_t;                  // r_zero_2
    double alpha_p1 = (  dDy_t + dBz_t) * 0.5;  // r_plus_1
    double alpha_p2 = ( -dDz_t + dBy_t) * 0.5;  // r_plus_2

    // ---- Modified eigenvalue magnitudes with moving-interface correction ----
    // In tetrad time/distance, eigenvalues of the moving-interface modified
    // flat-Maxwell Jacobian are 0 - tilde beta^x/alpha and +-1 - tilde beta^x/alpha,
    // where tilde beta^x/alpha = beta^x / (alpha sqrt(gamma^xx)). To use these
    // as dissipation magnitudes against the COORD-basis centered flux, we
    // convert from tetrad to coord time units. A wave at unit tetrad speed
    // covers (alpha * sqrt(gamma^xx)) in coord time per coord distance, so
    // multiply each tetrad |lambda| by alpha * sqrt(gamma^xx). Verified
    // analytically: max |lambda|_coord = alpha * sqrt(gamma^xx) * (1 + |bxa|)
    // = alpha * sqrt(gamma^xx) + |beta^x|, which matches the coord-basis
    // curved-Maxwell eigenvalue formula used by alpha_quad_x.
    double bxa = beta_x * inv_sqrt_gamma_xx / alpha;
    double scale_t2c = alpha * sqrt_gamma_xx;
    double absL_m = fabs(-1.0 - bxa) * scale_t2c;
    double absL_0 = fabs( -bxa)      * scale_t2c;
    double absL_p = fabs( 1.0 - bxa) * scale_t2c;

    // CFL contribution: max wave speed in coord-time/coord-distance units.
    // The wrapper applies (2/dx) downstream.
    if (absL_m > cflrate) cflrate = absL_m;
    if (absL_p > cflrate) cflrate = absL_p;
    if (absL_0 > cflrate) cflrate = absL_0;

    // ---- Roe dissipation in tetrad: (1/2) sum_p |lambda_p^mod| alpha_p r_p^a ----
    // For each component (D^x, D^y, D^z, B^x, B^y, B^z in tetrad), the
    // dissipation pulls only the eigenvectors that have a non-zero entry in
    // that component slot. From paper eq. 57:
    //   D_x slot: only r_zero_2 is non-zero (entry = +1).
    //   D_y slot: r_minus_1 (-1), r_plus_1 (+1).
    //   D_z slot: r_minus_2 (+1), r_plus_2 (-1).
    //   B_x slot: only r_zero_1 (+1).
    //   B_y slot: r_minus_2 (+1), r_plus_2 (+1).
    //   B_z slot: r_minus_1 (+1), r_plus_1 (+1).
    double diss_Dx_t = 0.5 * absL_0 * alpha_02;
    double diss_Dy_t = 0.5 * (-absL_m * alpha_m1 + absL_p * alpha_p1);
    double diss_Dz_t = 0.5 * ( absL_m * alpha_m2 - absL_p * alpha_p2);
    double diss_Bx_t = 0.5 * absL_0 * alpha_01;
    double diss_By_t = 0.5 * ( absL_m * alpha_m2 + absL_p * alpha_p2);
    double diss_Bz_t = 0.5 * ( absL_m * alpha_m1 + absL_p * alpha_p1);

    // ---- Back-transform to coord basis: diss^i = (M^{-1})^i_a diss^a_t ----
    // Rows of M^{-1}:
    //   row x: (Minv_xx, 0,        0      )
    //   row y: (0,       Minv_yy,  0      )
    //   row z: (Minv_zx, 0,        Minv_zz)
    double diss_Dx = Minv_xx * diss_Dx_t;
    double diss_Dy = Minv_yy * diss_Dy_t;
    double diss_Dz = Minv_zx * diss_Dx_t + Minv_zz * diss_Dz_t;
    double diss_Bx = Minv_xx * diss_Bx_t;
    double diss_By = Minv_yy * diss_By_t;
    double diss_Bz = Minv_zx * diss_Bx_t + Minv_zz * diss_Bz_t;

    // ---- Output: F* = (1/2)(F_L + F_R) - dissipation ----
    // flux_l/r layout: 2 nodes x 6 components, indexed as flux[comp*2 + n].
    // Output layout: same (out[comp*2 + n]).
    out[n + 0*2]  = 0.5*(flux_l[0*2+n] + flux_r[0*2+n]) - diss_Dx;  // D^x equation (flat F^x = 0; centered = 0 + shift terms)
    out[n + 1*2]  = 0.5*(flux_l[1*2+n] + flux_r[1*2+n]) - diss_Dy;  // D^y equation
    out[n + 2*2]  = 0.5*(flux_l[2*2+n] + flux_r[2*2+n]) - diss_Dz;  // D^z equation
    out[n + 3*2]  = 0.5*(flux_l[3*2+n] + flux_r[3*2+n]) - diss_Bx;  // B^x equation
    out[n + 4*2]  = 0.5*(flux_l[4*2+n] + flux_r[4*2+n]) - diss_By;  // B^y equation
    out[n + 5*2]  = 0.5*(flux_l[5*2+n] + flux_r[5*2+n]) - diss_Bz;  // B^z equation
  }

  // Aux components (phi, psi divergence cleaning) decoupled; zero them in
  // the dir=0 slot, matching the existing curved-LLF kernel's convention.
  out[12] = 0.0;
  out[13] = 0.0;
  out[14] = 0.0;
  out[15] = 0.0;

  // Suppress -Wunused-parameter for theta_pole; radial faces do not encounter
  // a theta-pole (GL nodes sit at eta_y = +-1/sqrt(3), never on the pole).
  (void)theta_pole;

  return 1.5 * dx10 * cflrate;
}
