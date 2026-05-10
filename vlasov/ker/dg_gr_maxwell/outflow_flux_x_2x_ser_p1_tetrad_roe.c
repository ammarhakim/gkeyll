#include <math.h>
#include <stdio.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Tetrad-frame characteristic-outflow flux at the +/- x boundary face for
// 2D serendipity p=1.
//
// Algorithm at each face GL node:
//   1. Read skin-cell modal coefficients; evaluate at face GL node (edge-
//      specific formula, same as the existing outflow_flux_x kernel).
//   2. Build M^a_i (tetrad transform) and M^{-1} (back-transform). Closed
//      form for our h_rtheta = h_thetaphi = 0 metric -- same as the LLF
//      tetrad-Roe kernel.
//   3. Transform U_skin to tetrad: Ũ^a = M^a_i U^i.
//   4. Wave decomposition using paper eq. 57 right eigenvectors of the flat-
//      Maxwell flux Jacobian along tetrad-x = face normal:
//        alpha_minus_1 = (-Ũ^y + Ũ^z_B)/2  [eigenvalue lambda_tetrad = -1]
//        alpha_minus_2 = ( Ũ^z + Ũ^y_B)/2  [-1]
//        alpha_zero_1  =   Ũ^x_B           [0]
//        alpha_zero_2  =   Ũ^x_D           [0]
//        alpha_plus_1  = ( Ũ^y + Ũ^z_B)/2  [+1]
//        alpha_plus_2  = (-Ũ^z + Ũ^y_B)/2  [+1]
//   5. Coord-frame eigenvalues (with moving-interface beta_tilde^x/alpha):
//        lambda_p^coord = (lambda_p_tetrad - beta_tilde^x/alpha) * alpha * sqrt(gamma^xx)
//      In our metric this equals the standard curved-Maxwell -beta^x +/- alpha*sqrt(gamma^xx)
//      eigenvalues (as verified analytically); we compute it via the tetrad
//      formulation here for consistency with the LLF tetrad-Roe kernel.
//   6. OUTGOING-only: keep lambda_p^coord whose sign matches the outgoing
//      direction. For edge=0 (lower-x boundary), outgoing means lambda < 0
//      (waves leaving the domain in -x); for edge=1 (upper-x), outgoing means
//      lambda > 0. Incoming waves get lambda zeroed -- their amplitude is
//      removed from the boundary flux.
//   7. Tetrad-frame outgoing flux: F_outflow^a = sum_p lambda_p^outgoing
//      alpha_p r_p^a. (NOTE: no 1/2 prefactor here -- this is the one-sided
//      outflow projection, not the centered LLF flux.)
//   8. Back-transform: F^i_coord = (M^{-1})^i_a F_outflow^a.
//
// Output: writes the F^i values to conf_flux_surf in the dir=0 slot, matching
// the layout of outflow_flux_x_2x_ser_p1.
//
// cflrate returns 1.5 * (2/dx) * |max lambda_outgoing|, same convention as the
// existing outflow_flux_x kernel.
GKYL_CU_DH double
outflow_flux_x_2x_ser_p1_tetrad_roe(const double *dxv, int edge,
  const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
  const double *det_h_nodal, const double *field_con_skin,
  double* GKYL_RESTRICT conf_flux_surf)
{
  double dx10 = 2.0 / dxv[0];

  const double *h_xz_nodal = &h_ij_nodal[4];
  const double *h_yy_nodal = &h_ij_nodal[6];
  const double *h_yz_nodal = &h_ij_nodal[8];
  const double *h_zz_nodal = &h_ij_nodal[10];
  const double *shift_x_n  = &shift_nodal[0];

  const double *JDx_s = &field_con_skin[0];
  const double *JDy_s = &field_con_skin[4];
  const double *JDz_s = &field_con_skin[8];
  const double *JBx_s = &field_con_skin[12];
  const double *JBy_s = &field_con_skin[16];
  const double *JBz_s = &field_con_skin[20];

  // Modal-to-nodal eval of U_skin at the 2 face quadrature nodes.
  // Edge=0 (lower-x face): skin sits to the RIGHT of the face, so evaluate
  // at the skin's xi=-1 (matches JD*_con_r_n in alpha_quad_x).
  // Edge=1 (upper-x face): skin to the LEFT, evaluate at skin's xi=+1
  // (matches JD*_con_l_n).
  double JDx_n[2], JDy_n[2], JDz_n[2], JBx_n[2], JBy_n[2], JBz_n[2];
  if (edge == 0) {
    JDx_n[0] = 0.8660254037844386*JDx_s[3] - 0.5*JDx_s[2] - 0.8660254037844386*JDx_s[1] + 0.5*JDx_s[0];
    JDy_n[0] = 0.8660254037844386*JDy_s[3] - 0.5*JDy_s[2] - 0.8660254037844386*JDy_s[1] + 0.5*JDy_s[0];
    JDz_n[0] = 0.8660254037844386*JDz_s[3] - 0.5*JDz_s[2] - 0.8660254037844386*JDz_s[1] + 0.5*JDz_s[0];
    JBx_n[0] = 0.8660254037844386*JBx_s[3] - 0.5*JBx_s[2] - 0.8660254037844386*JBx_s[1] + 0.5*JBx_s[0];
    JBy_n[0] = 0.8660254037844386*JBy_s[3] - 0.5*JBy_s[2] - 0.8660254037844386*JBy_s[1] + 0.5*JBy_s[0];
    JBz_n[0] = 0.8660254037844386*JBz_s[3] - 0.5*JBz_s[2] - 0.8660254037844386*JBz_s[1] + 0.5*JBz_s[0];
    JDx_n[1] = -(0.8660254037844386*JDx_s[3]) + 0.5*JDx_s[2] - 0.8660254037844386*JDx_s[1] + 0.5*JDx_s[0];
    JDy_n[1] = -(0.8660254037844386*JDy_s[3]) + 0.5*JDy_s[2] - 0.8660254037844386*JDy_s[1] + 0.5*JDy_s[0];
    JDz_n[1] = -(0.8660254037844386*JDz_s[3]) + 0.5*JDz_s[2] - 0.8660254037844386*JDz_s[1] + 0.5*JDz_s[0];
    JBx_n[1] = -(0.8660254037844386*JBx_s[3]) + 0.5*JBx_s[2] - 0.8660254037844386*JBx_s[1] + 0.5*JBx_s[0];
    JBy_n[1] = -(0.8660254037844386*JBy_s[3]) + 0.5*JBy_s[2] - 0.8660254037844386*JBy_s[1] + 0.5*JBy_s[0];
    JBz_n[1] = -(0.8660254037844386*JBz_s[3]) + 0.5*JBz_s[2] - 0.8660254037844386*JBz_s[1] + 0.5*JBz_s[0];
  } else {
    JDx_n[0] = -(0.8660254037844386*JDx_s[3]) - 0.5*JDx_s[2] + 0.8660254037844386*JDx_s[1] + 0.5*JDx_s[0];
    JDy_n[0] = -(0.8660254037844386*JDy_s[3]) - 0.5*JDy_s[2] + 0.8660254037844386*JDy_s[1] + 0.5*JDy_s[0];
    JDz_n[0] = -(0.8660254037844386*JDz_s[3]) - 0.5*JDz_s[2] + 0.8660254037844386*JDz_s[1] + 0.5*JDz_s[0];
    JBx_n[0] = -(0.8660254037844386*JBx_s[3]) - 0.5*JBx_s[2] + 0.8660254037844386*JBx_s[1] + 0.5*JBx_s[0];
    JBy_n[0] = -(0.8660254037844386*JBy_s[3]) - 0.5*JBy_s[2] + 0.8660254037844386*JBy_s[1] + 0.5*JBy_s[0];
    JBz_n[0] = -(0.8660254037844386*JBz_s[3]) - 0.5*JBz_s[2] + 0.8660254037844386*JBz_s[1] + 0.5*JBz_s[0];
    JDx_n[1] = 0.8660254037844386*JDx_s[3] + 0.5*JDx_s[2] + 0.8660254037844386*JDx_s[1] + 0.5*JDx_s[0];
    JDy_n[1] = 0.8660254037844386*JDy_s[3] + 0.5*JDy_s[2] + 0.8660254037844386*JDy_s[1] + 0.5*JDy_s[0];
    JDz_n[1] = 0.8660254037844386*JDz_s[3] + 0.5*JDz_s[2] + 0.8660254037844386*JDz_s[1] + 0.5*JDz_s[0];
    JBx_n[1] = 0.8660254037844386*JBx_s[3] + 0.5*JBx_s[2] + 0.8660254037844386*JBx_s[1] + 0.5*JBx_s[0];
    JBy_n[1] = 0.8660254037844386*JBy_s[3] + 0.5*JBy_s[2] + 0.8660254037844386*JBy_s[1] + 0.5*JBy_s[0];
    JBz_n[1] = 0.8660254037844386*JBz_s[3] + 0.5*JBz_s[2] + 0.8660254037844386*JBz_s[1] + 0.5*JBz_s[0];
  }

  double cflrate = 0.0;
  for (int n = 0; n < 2; ++n) {
    // ---- Geometry at this node ----
    double h_yy = h_yy_nodal[n];
    double h_yz = h_yz_nodal[n];
    double h_zz = h_zz_nodal[n];
    double h_xz = h_xz_nodal[n];
    double Jc   = det_h_nodal[n];
    double alpha   = lapse_nodal[n];
    double beta_x  = shift_x_n[n];

    // ---- M and M^{-1} (same closed form as the LLF tetrad-Roe kernel) ----
    double sqrt_gamma_xx     = sqrt(h_yy*h_zz - h_yz*h_yz) / Jc;
    double inv_sqrt_gamma_xx = 1.0 / sqrt_gamma_xx;
    double sqrt_h_yy = sqrt(h_yy);
    double sqrt_h_zz = sqrt(h_zz);

    double M_xx = inv_sqrt_gamma_xx;
    double M_yy = sqrt_h_yy;
    double M_zz = sqrt_h_zz;
    double M_zx = h_xz / sqrt_h_zz;

    double Minv_xx = sqrt_gamma_xx;
    double Minv_yy = 1.0 / sqrt_h_yy;
    double Minv_zx = -h_xz * sqrt_gamma_xx / h_zz;
    double Minv_zz = 1.0 / sqrt_h_zz;

    // ---- Tetrad-frame skin field ----
    double Ux_t = M_xx * JDx_n[n];
    double Uy_t = M_yy * JDy_n[n];
    double Uz_t = M_zx * JDx_n[n] + M_zz * JDz_n[n];
    double Vx_t = M_xx * JBx_n[n];
    double Vy_t = M_yy * JBy_n[n];
    double Vz_t = M_zx * JBx_n[n] + M_zz * JBz_n[n];

    // ---- Wave strengths (paper eq. 57; same formulas as LLF kernel) ----
    double alpha_m1 = (-Uy_t + Vz_t) * 0.5;
    double alpha_m2 = ( Uz_t + Vy_t) * 0.5;
    double alpha_01 =   Vx_t;
    double alpha_02 =   Ux_t;
    double alpha_p1 = ( Uy_t + Vz_t) * 0.5;
    double alpha_p2 = (-Uz_t + Vy_t) * 0.5;

    // ---- Signed coord-frame eigenvalues with moving-interface shift ----
    double bxa = beta_x * inv_sqrt_gamma_xx / alpha;
    double scale_t2c = alpha * sqrt_gamma_xx;
    double lam_m_coord = (-1.0 - bxa) * scale_t2c;  // r_minus_1, r_minus_2
    double lam_0_coord = (    -bxa)   * scale_t2c;  // r_zero_1, r_zero_2
    double lam_p_coord = ( 1.0 - bxa) * scale_t2c;  // r_plus_1, r_plus_2

    // ---- Outgoing-only eigenvalues ----
    double lam_m_out, lam_0_out, lam_p_out;
    if (edge == 0) {
      // Lower-x boundary: outgoing = negative coord-frame eigenvalue.
      lam_m_out = (lam_m_coord < 0.0) ? lam_m_coord : 0.0;
      lam_0_out = (lam_0_coord < 0.0) ? lam_0_coord : 0.0;
      lam_p_out = (lam_p_coord < 0.0) ? lam_p_coord : 0.0;
    } else {
      // Upper-x boundary: outgoing = positive coord-frame eigenvalue.
      lam_m_out = (lam_m_coord > 0.0) ? lam_m_coord : 0.0;
      lam_0_out = (lam_0_coord > 0.0) ? lam_0_coord : 0.0;
      lam_p_out = (lam_p_coord > 0.0) ? lam_p_coord : 0.0;
    }

    // CFL contribution: max magnitude of outgoing eigenvalue (note: not the
    // full eigenvalue spectrum, since the kernel only applies the outgoing
    // part).
    if (fabs(lam_m_out) > cflrate) cflrate = fabs(lam_m_out);
    if (fabs(lam_p_out) > cflrate) cflrate = fabs(lam_p_out);
    if (fabs(lam_0_out) > cflrate) cflrate = fabs(lam_0_out);

    // ---- One-sided outflow flux in tetrad: F = sum_p lambda_p^out alpha_p r_p ----
    // Same wave-vector layout as the LLF tetrad-Roe kernel's dissipation
    // formulas, BUT signed and without the 1/2 prefactor (this is the
    // one-sided projection, not the centered LLF flux).
    double Ft_Dx = lam_0_out * alpha_02;
    double Ft_Dy = (-lam_m_out * alpha_m1 + lam_p_out * alpha_p1);
    double Ft_Dz = ( lam_m_out * alpha_m2 - lam_p_out * alpha_p2);
    double Ft_Bx = lam_0_out * alpha_01;
    double Ft_By = ( lam_m_out * alpha_m2 + lam_p_out * alpha_p2);
    double Ft_Bz = ( lam_m_out * alpha_m1 + lam_p_out * alpha_p1);

    // ---- Back-transform to coord basis ----
    double F_Dx = Minv_xx * Ft_Dx;
    double F_Dy = Minv_yy * Ft_Dy;
    double F_Dz = Minv_zx * Ft_Dx + Minv_zz * Ft_Dz;
    double F_Bx = Minv_xx * Ft_Bx;
    double F_By = Minv_yy * Ft_By;
    double F_Bz = Minv_zx * Ft_Bx + Minv_zz * Ft_Bz;

    // Component-major layout: (comp 0..5) x (2 nodes), node index `n`.
    conf_flux_surf[0  + n] = F_Dx;
    conf_flux_surf[2  + n] = F_Dy;
    conf_flux_surf[4  + n] = F_Dz;
    conf_flux_surf[6  + n] = F_Bx;
    conf_flux_surf[8  + n] = F_By;
    conf_flux_surf[10 + n] = F_Bz;
  }

  // Cleaning slots: no outflow contribution.
  conf_flux_surf[12] = 0.0;
  conf_flux_surf[13] = 0.0;
  conf_flux_surf[14] = 0.0;
  conf_flux_surf[15] = 0.0;

  return 1.5 * dx10 * cflrate;
}
