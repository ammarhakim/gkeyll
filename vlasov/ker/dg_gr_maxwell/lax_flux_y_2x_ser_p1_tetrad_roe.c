#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Tetrad-frame Roe flux on the theta-direction face, 2D ser p=1.
//
// The y-face Gram-Schmidt tetrad (with face normal e_y_tetrad = unit normal to
// Sigma_theta in the spatial slice, e_x_tetrad = e_r in-face, e_z_tetrad =
// orthogonalized e_phi within Sigma_theta) gives, for our spherical Kerr-
// Schild metric (h_rtheta = h_thetaphi = 0):
//
//   M^a_i = [[ 0,                sqrt(h_theta_theta), 0                          ],
//            [ sqrt(h_rr),       0,                    h_rphi/sqrt(h_rr)         ],
//            [ 0,                0,                    D_theta/sqrt(h_rr)        ]]
//
// where D_theta = sqrt(h_rr h_phi_phi - h_rphi^2). Analytically verified
// M^T M = h_ij. The inverse is
//
//   M^{-1}^i_a = [[ 0,                       1/sqrt(h_rr),     -h_rphi/(sqrt(h_rr) D_theta) ],
//                 [ 1/sqrt(h_theta_theta),   0,                 0                            ],
//                 [ 0,                       0,                 sqrt(h_rr)/D_theta           ]].
//
// At the theta pole D_theta proportional to |sin theta| -> 0 and the inverse
// transform diverges. We handle this by zeroing the dissipation at theta-pole
// faces (theta_pole=1), matching the existing curved-LLF kernel convention.
//
// For our spherical Kerr-Schild, beta^theta = 0 (the shift is purely radial),
// so the moving-interface velocity tilde beta^y_face/alpha = M^x_tet_i beta^i
// / alpha = sqrt(h_theta_theta) beta^theta / alpha = 0. The flat-Maxwell
// eigenvalues are therefore exactly +-1, 0 (no shift modification).
//
// IMPORTANT: with the tetrad oriented so that "tetrad x" = face normal (=
// e_theta for the y-face), the flat-Maxwell flux Jacobian whose eigenvectors
// we need is the one along tetrad x (the face normal direction). Paper eq. 57
// gives those eigenvectors; they are the SAME as for the x-face kernel, just
// applied to a different M-transformed jump. Using the y-direction coord-
// frame flat-Maxwell eigenvectors here would be incorrect (those decompose
// the flux Jacobian along the COORD y axis, not the tetrad-x = face-normal
// axis). The flat-limit test in ctest_dg_gr_maxwell_tetrad.c originally
// caught this; we now use paper eq. 57 universally.
//
// Paper eq. 57 right eigenvectors of the flat-Maxwell flux Jacobian along
// the face-normal direction, in tetrad components (D^x, D^y, D^z, B^x, B^y, B^z):
//   r_minus_1 = ( 0,-1, 0, 0, 0, 1),    lambda = -1
//   r_minus_2 = ( 0, 0, 1, 0, 1, 0),    lambda = -1
//   r_zero_1  = ( 0, 0, 0, 1, 0, 0),    lambda =  0   (B^x_tet rest mode)
//   r_zero_2  = ( 1, 0, 0, 0, 0, 0),    lambda =  0   (D^x_tet rest mode)
//   r_plus_1  = ( 0, 1, 0, 0, 0, 1),    lambda = +1
//   r_plus_2  = ( 0, 0,-1, 0, 1, 0),    lambda = +1
//
// Wave strengths (left eigenvectors r_p / |r_p|^2):
//   alpha_minus_1 = (-dD^y_t + dB^z_t)/2
//   alpha_minus_2 = ( dD^z_t + dB^y_t)/2
//   alpha_zero_1  =   dB^x_t
//   alpha_zero_2  =   dD^x_t
//   alpha_plus_1  = ( dD^y_t + dB^z_t)/2
//   alpha_plus_2  = (-dD^z_t + dB^y_t)/2
//
// Roe dissipation in tetrad: tilde diss^a = (1/2) sum_p |lambda_p| alpha_p r_p^a:
//   diss_Dx_t = (1/2) |lambda_0| alpha_zero_2
//   diss_Dy_t = (1/2) (-|lambda_-| alpha_minus_1 + |lambda_+| alpha_plus_1)
//   diss_Dz_t = (1/2) ( |lambda_-| alpha_minus_2 - |lambda_+| alpha_plus_2)
//   diss_Bx_t = (1/2) |lambda_0| alpha_zero_1
//   diss_By_t = (1/2) ( |lambda_-| alpha_minus_2 + |lambda_+| alpha_plus_2)
//   diss_Bz_t = (1/2) ( |lambda_-| alpha_minus_1 + |lambda_+| alpha_plus_1)
//
// Convert tetrad eigenvalues to coord time/distance: multiply by alpha *
// sqrt(gamma^yy). gamma^yy = (h_rr h_phi_phi - h_rphi^2)/det h = D_theta^2/det h,
// so sqrt(gamma^yy) = D_theta/J_c. Verified: max |lambda|_coord = alpha *
// D_theta/J_c, which is the same expression alpha_quad_y produces for the
// curved-Maxwell wave speed (since beta^y = 0).
GKYL_CU_DH double
lax_flux_y_2x_ser_p1_tetrad_roe(const double *dxv, const int theta_pole,
  const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *det_h,
  const double *flux_l, const double *flux_r,
  const double *field_con_l, const double *field_con_r,
  double* GKYL_RESTRICT conf_flux_surf)
{
  double dx11 = 2.0 / dxv[1];

  // Output goes to the dir=1 slice of conf_flux_surf (offset 16).
  double *out = &conf_flux_surf[16];
  double cflrate = 0.0;

  // Conserved-variable pointers.
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

  // h_ij and lapse/shift slicing per face GL node.
  const double *h_xx_n = &h_ij_nodal[0];   // h_rr
  const double *h_xy_n = &h_ij_nodal[2];   // h_rtheta = 0 in our metric
  const double *h_xz_n = &h_ij_nodal[4];   // h_rphi
  const double *h_yy_n = &h_ij_nodal[6];   // h_theta_theta
  const double *h_yz_n = &h_ij_nodal[8];   // h_thetaphi = 0 in our metric
  const double *h_zz_n = &h_ij_nodal[10];  // h_phiphi
  const double *shift_y_n = &shift_nodal[2];  // beta^theta (= 0 in Kerr-Schild)

  // Suppress unused-variable warnings for the metric components that vanish
  // in our metric; if the kernel is ever generalized to non-orthogonal
  // (theta) face metrics, h_xy_n and h_yz_n would be needed.
  (void)h_xy_n;
  (void)h_yz_n;

  for (int n = 0; n < 2; ++n) {
    // ---- Modal-to-nodal eval at face GL node n ----
    // Coefficients copied verbatim from dg_gr_maxwell_alpha_quad_y_2x_ser_p1.c
    // and lax_flux_y_2x_ser_p1.c. L cell is at the face's *upper* y-edge
    // (eta_y = +1); R cell is at the face's *lower* y-edge (eta_y = -1).
    double JDx_L, JDy_L, JDz_L, JBx_L, JBy_L, JBz_L;
    double JDx_R, JDy_R, JDz_R, JBx_R, JBy_R, JBz_R;
    if (n == 0) {
      JDx_L = -(0.8660254037844386)*JDx_l[3] + 0.8660254037844386*JDx_l[2] - 0.5*JDx_l[1] + 0.5*JDx_l[0];
      JDy_L = -(0.8660254037844386)*JDy_l[3] + 0.8660254037844386*JDy_l[2] - 0.5*JDy_l[1] + 0.5*JDy_l[0];
      JDz_L = -(0.8660254037844386)*JDz_l[3] + 0.8660254037844386*JDz_l[2] - 0.5*JDz_l[1] + 0.5*JDz_l[0];
      JBx_L = -(0.8660254037844386)*JBx_l[3] + 0.8660254037844386*JBx_l[2] - 0.5*JBx_l[1] + 0.5*JBx_l[0];
      JBy_L = -(0.8660254037844386)*JBy_l[3] + 0.8660254037844386*JBy_l[2] - 0.5*JBy_l[1] + 0.5*JBy_l[0];
      JBz_L = -(0.8660254037844386)*JBz_l[3] + 0.8660254037844386*JBz_l[2] - 0.5*JBz_l[1] + 0.5*JBz_l[0];

      JDx_R =   0.8660254037844386*JDx_r[3] - 0.8660254037844386*JDx_r[2] - 0.5*JDx_r[1] + 0.5*JDx_r[0];
      JDy_R =   0.8660254037844386*JDy_r[3] - 0.8660254037844386*JDy_r[2] - 0.5*JDy_r[1] + 0.5*JDy_r[0];
      JDz_R =   0.8660254037844386*JDz_r[3] - 0.8660254037844386*JDz_r[2] - 0.5*JDz_r[1] + 0.5*JDz_r[0];
      JBx_R =   0.8660254037844386*JBx_r[3] - 0.8660254037844386*JBx_r[2] - 0.5*JBx_r[1] + 0.5*JBx_r[0];
      JBy_R =   0.8660254037844386*JBy_r[3] - 0.8660254037844386*JBy_r[2] - 0.5*JBy_r[1] + 0.5*JBy_r[0];
      JBz_R =   0.8660254037844386*JBz_r[3] - 0.8660254037844386*JBz_r[2] - 0.5*JBz_r[1] + 0.5*JBz_r[0];
    } else {
      JDx_L =   0.8660254037844386*JDx_l[3] + 0.8660254037844386*JDx_l[2] + 0.5*JDx_l[1] + 0.5*JDx_l[0];
      JDy_L =   0.8660254037844386*JDy_l[3] + 0.8660254037844386*JDy_l[2] + 0.5*JDy_l[1] + 0.5*JDy_l[0];
      JDz_L =   0.8660254037844386*JDz_l[3] + 0.8660254037844386*JDz_l[2] + 0.5*JDz_l[1] + 0.5*JDz_l[0];
      JBx_L =   0.8660254037844386*JBx_l[3] + 0.8660254037844386*JBx_l[2] + 0.5*JBx_l[1] + 0.5*JBx_l[0];
      JBy_L =   0.8660254037844386*JBy_l[3] + 0.8660254037844386*JBy_l[2] + 0.5*JBy_l[1] + 0.5*JBy_l[0];
      JBz_L =   0.8660254037844386*JBz_l[3] + 0.8660254037844386*JBz_l[2] + 0.5*JBz_l[1] + 0.5*JBz_l[0];

      JDx_R = -(0.8660254037844386)*JDx_r[3] - 0.8660254037844386*JDx_r[2] + 0.5*JDx_r[1] + 0.5*JDx_r[0];
      JDy_R = -(0.8660254037844386)*JDy_r[3] - 0.8660254037844386*JDy_r[2] + 0.5*JDy_r[1] + 0.5*JDy_r[0];
      JDz_R = -(0.8660254037844386)*JDz_r[3] - 0.8660254037844386*JDz_r[2] + 0.5*JDz_r[1] + 0.5*JDz_r[0];
      JBx_R = -(0.8660254037844386)*JBx_r[3] - 0.8660254037844386*JBx_r[2] + 0.5*JBx_r[1] + 0.5*JBx_r[0];
      JBy_R = -(0.8660254037844386)*JBy_r[3] - 0.8660254037844386*JBy_r[2] + 0.5*JBy_r[1] + 0.5*JBy_r[0];
      JBz_R = -(0.8660254037844386)*JBz_r[3] - 0.8660254037844386*JBz_r[2] + 0.5*JBz_r[1] + 0.5*JBz_r[0];
    }

    double dDx = JDx_R - JDx_L;
    double dDy = JDy_R - JDy_L;
    double dDz = JDz_R - JDz_L;
    double dBx = JBx_R - JBx_L;
    double dBy = JBy_R - JBy_L;
    double dBz = JBz_R - JBz_L;

    // ---- Geometry at this node ----
    double h_rr  = h_xx_n[n];
    double h_rp  = h_xz_n[n];   // h_rphi
    double h_tt  = h_yy_n[n];   // h_theta_theta
    double h_pp  = h_zz_n[n];   // h_phiphi
    double Jc    = det_h[n];
    double alpha = lapse_nodal[n];
    double beta_t = shift_y_n[n];  // beta^theta (= 0 in Kerr-Schild for our preset)

    if (theta_pole) {
      // At the theta pole, the back-transform M^{-1} diverges (M^{-1}_{phi,phi}
      // = sqrt(h_rr)/D_theta, with D_theta proportional to |sin theta| -> 0).
      // Match the existing curved-LLF convention: use centered flux only,
      // zero dissipation. This produces a reflective-like behavior at the
      // coordinate singularity.
      out[n+0*2] = 0.5*(flux_l[0*2+n] + flux_r[0*2+n]);
      out[n+1*2] = 0.5*(flux_l[1*2+n] + flux_r[1*2+n]);
      out[n+2*2] = 0.5*(flux_l[2*2+n] + flux_r[2*2+n]);
      out[n+3*2] = 0.5*(flux_l[3*2+n] + flux_r[3*2+n]);
      out[n+4*2] = 0.5*(flux_l[4*2+n] + flux_r[4*2+n]);
      out[n+5*2] = 0.5*(flux_l[5*2+n] + flux_r[5*2+n]);
      continue;
    }

    // ---- Build M^a_i (y-face Gram-Schmidt tetrad) and tetrad-frame jumps ----
    // For our metric (h_rtheta = h_thetaphi = 0):
    //   M = [[ 0,           sqrt(h_tt),   0                ],   tetrad x = face normal e_theta
    //        [ sqrt(h_rr),  0,            h_rp/sqrt(h_rr)  ],   tetrad y = in-face e_r
    //        [ 0,           0,            D_theta/sqrt(h_rr)]]. tetrad z = in-face e_phi
    double D_theta = sqrt(h_rr * h_pp - h_rp * h_rp);
    double sqrt_h_rr = sqrt(h_rr);
    double sqrt_h_tt = sqrt(h_tt);
    double inv_sqrt_h_rr = 1.0 / sqrt_h_rr;
    double inv_D_theta = 1.0 / D_theta;

    // M^a_i:
    double M_xt = sqrt_h_tt;             // M[x_tet, theta]
    double M_yr = sqrt_h_rr;             // M[y_tet, r]
    double M_yp = h_rp * inv_sqrt_h_rr;  // M[y_tet, phi]
    double M_zp = D_theta * inv_sqrt_h_rr; // M[z_tet, phi]

    // Tetrad-frame jumps: dU_t^a = M^a_i dU^i.
    // For D and B independently:
    double dDx_t = M_xt * dDy;                       // D-component along tetrad x = M * dD where M only picks D^theta = dDy
    double dDy_t = M_yr * dDx + M_yp * dDz;          // tetrad y picks dD^r and dD^phi
    double dDz_t = M_zp * dDz;                        // tetrad z picks only dD^phi
    double dBx_t = M_xt * dBy;
    double dBy_t = M_yr * dBx + M_yp * dBz;
    double dBz_t = M_zp * dBz;

    // ---- Wave strengths from paper eq. 57 (tetrad-x face-normal flux) ----
    // Same eigenstructure as the x-face kernel; only the M (and thus the
    // tetrad jumps dD*_t, dB*_t) differ between directions.
    double alpha_m1 = (-dDy_t + dBz_t) * 0.5;  // r_minus_1
    double alpha_m2 = ( dDz_t + dBy_t) * 0.5;  // r_minus_2
    double alpha_01 =   dBx_t;                  // r_zero_1
    double alpha_02 =   dDx_t;                  // r_zero_2
    double alpha_p1 = ( dDy_t + dBz_t) * 0.5;  // r_plus_1
    double alpha_p2 = (-dDz_t + dBy_t) * 0.5;  // r_plus_2

    // ---- Modified eigenvalue magnitudes ----
    // tilde beta^x_face/alpha for the y-face: beta_tilde = M[x_tet, i] beta^i
    // = sqrt(h_tt) * beta^theta. For Kerr-Schild beta^theta = 0, so the
    // moving-interface correction vanishes.
    double beta_tilde = sqrt_h_tt * beta_t;
    double bxa = beta_tilde / alpha;
    double scale_t2c = alpha * D_theta / Jc;  // = alpha * sqrt(gamma^yy)
    double absL_m = fabs(-1.0 - bxa) * scale_t2c;
    double absL_0 = fabs(    -bxa)   * scale_t2c;
    double absL_p = fabs( 1.0 - bxa) * scale_t2c;

    if (absL_m > cflrate) cflrate = absL_m;
    if (absL_p > cflrate) cflrate = absL_p;
    if (absL_0 > cflrate) cflrate = absL_0;

    // ---- Roe dissipation in tetrad (paper eq. 57, same as x-kernel) ----
    double diss_Dx_t = 0.5 * absL_0 * alpha_02;
    double diss_Dy_t = 0.5 * (-absL_m * alpha_m1 + absL_p * alpha_p1);
    double diss_Dz_t = 0.5 * ( absL_m * alpha_m2 - absL_p * alpha_p2);
    double diss_Bx_t = 0.5 * absL_0 * alpha_01;
    double diss_By_t = 0.5 * ( absL_m * alpha_m2 + absL_p * alpha_p2);
    double diss_Bz_t = 0.5 * ( absL_m * alpha_m1 + absL_p * alpha_p1);

    // ---- Back-transform to coord basis: diss^i = (M^{-1})^i_a diss_t^a ----
    // M^{-1} = [[ 0,             1/sqrt(h_rr),    -h_rp/(sqrt(h_rr) D_theta) ],
    //           [ 1/sqrt(h_tt),  0,                0                          ],
    //           [ 0,             0,                sqrt(h_rr)/D_theta         ]].
    double Minv_yr_inv = inv_sqrt_h_rr;
    double Minv_xr_3   = -h_rp * inv_sqrt_h_rr * inv_D_theta;
    double Minv_xt_2   = 1.0 / sqrt_h_tt;
    double Minv_zp_3   = sqrt_h_rr * inv_D_theta;

    double diss_Dr  = Minv_yr_inv * diss_Dy_t + Minv_xr_3 * diss_Dz_t;
    double diss_Dt  = Minv_xt_2 * diss_Dx_t;
    double diss_Dp  = Minv_zp_3 * diss_Dz_t;
    double diss_Br  = Minv_yr_inv * diss_By_t + Minv_xr_3 * diss_Bz_t;
    double diss_Bt  = Minv_xt_2 * diss_Bx_t;
    double diss_Bp  = Minv_zp_3 * diss_Bz_t;

    // ---- Output: F* = (1/2)(F_L + F_R) - dissipation ----
    out[n + 0*2]  = 0.5*(flux_l[0*2+n] + flux_r[0*2+n]) - diss_Dr;  // D^r equation
    out[n + 1*2]  = 0.5*(flux_l[1*2+n] + flux_r[1*2+n]) - diss_Dt;  // D^theta equation
    out[n + 2*2]  = 0.5*(flux_l[2*2+n] + flux_r[2*2+n]) - diss_Dp;  // D^phi equation
    out[n + 3*2]  = 0.5*(flux_l[3*2+n] + flux_r[3*2+n]) - diss_Br;  // B^r equation
    out[n + 4*2]  = 0.5*(flux_l[4*2+n] + flux_r[4*2+n]) - diss_Bt;  // B^theta equation
    out[n + 5*2]  = 0.5*(flux_l[5*2+n] + flux_r[5*2+n]) - diss_Bp;  // B^phi equation
  }

  // Aux components (phi, psi divergence cleaning) decoupled; zero them in
  // the dir=1 slot.
  out[12] = 0.0;
  out[13] = 0.0;
  out[14] = 0.0;
  out[15] = 0.0;

  return 1.5 * dx11 * cflrate;
}
