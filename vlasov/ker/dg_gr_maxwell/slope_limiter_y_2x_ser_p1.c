#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Tetrad-frame characteristic slope limiter on the theta (y-direction) face,
// 2D serendipity p=1.
//
// Mirror of slope_limiter_x_2x_ser_p1.c but operating on the cell's q[2]
// (theta-slope mode), using the y-face Gram-Schmidt tetrad whose tetrad-x is
// the face normal e_theta. See lax_flux_y_2x_ser_p1_tetrad_roe.c for the
// full derivation of the M matrix.
//
// Theta-pole handling: if either the lower or upper theta neighbor lies
// across the pole, the inter-cell jump is dominated by the reflective ghost
// rather than a physical gradient. Rather than try to incorporate the
// reflection's sign pattern, we zero the theta-slope and cross-slope for
// pole-adjacent cells -- the cleanest conservative choice. theta_pole_lo
// = 1 means cell sits against the lower theta boundary; theta_pole_up = 1
// means the upper.
GKYL_CU_DH void
slope_limiter_y_2x_ser_p1(double limiter_fac, const int theta_pole_lo, const int theta_pole_up,
  const double *lapse_nodal_face, const double *shift_nodal_face,
  const double *h_ij_nodal_face, const double *det_h_face,
  const double *field_con_l, const double *field_con_c, const double *field_con_r,
  double *field_con_out)
{
  const double *JDx_l = &field_con_l[0];
  const double *JDy_l = &field_con_l[4];
  const double *JDz_l = &field_con_l[8];
  const double *JBx_l = &field_con_l[12];
  const double *JBy_l = &field_con_l[16];
  const double *JBz_l = &field_con_l[20];

  const double *JDx_c = &field_con_c[0];
  const double *JDy_c = &field_con_c[4];
  const double *JDz_c = &field_con_c[8];
  const double *JBx_c = &field_con_c[12];
  const double *JBy_c = &field_con_c[16];
  const double *JBz_c = &field_con_c[20];

  const double *JDx_r = &field_con_r[0];
  const double *JDy_r = &field_con_r[4];
  const double *JDz_r = &field_con_r[8];
  const double *JBx_r = &field_con_r[12];
  const double *JBy_r = &field_con_r[16];
  const double *JBz_r = &field_con_r[20];

  double *JDx_out = &field_con_out[0];
  double *JDy_out = &field_con_out[4];
  double *JDz_out = &field_con_out[8];
  double *JBx_out = &field_con_out[12];
  double *JBy_out = &field_con_out[16];
  double *JBz_out = &field_con_out[20];

  // Default identity copy (so q[0], q[1], q[3] pass through; q[2] overwritten).
  for (int m = 0; m < 4; ++m) {
    JDx_out[m] = JDx_c[m]; JDy_out[m] = JDy_c[m]; JDz_out[m] = JDz_c[m];
    JBx_out[m] = JBx_c[m]; JBy_out[m] = JBy_c[m]; JBz_out[m] = JBz_c[m];
  }

  // Theta-pole adjacency: kill the theta slope and cross slope outright.
  // This is the conservative choice -- no physical gradient information is
  // available across the pole, and the M^{-1} back-transform diverges at the
  // pole anyway.
  if (theta_pole_lo || theta_pole_up) {
    JDx_out[2] = 0.0; JDy_out[2] = 0.0; JDz_out[2] = 0.0;
    JBx_out[2] = 0.0; JBy_out[2] = 0.0; JBz_out[2] = 0.0;
    JDx_out[3] = 0.0; JDy_out[3] = 0.0; JDz_out[3] = 0.0;
    JBx_out[3] = 0.0; JBy_out[3] = 0.0; JBz_out[3] = 0.0;
    return;
  }

  // ---- Face-averaged geometry (2 GL nodes -> 1 representative) ----
  double h_rr = 0.5 * (h_ij_nodal_face[0] + h_ij_nodal_face[1]);   // h_xx
  double h_rp = 0.5 * (h_ij_nodal_face[4] + h_ij_nodal_face[5]);   // h_xz = h_rphi
  double h_tt = 0.5 * (h_ij_nodal_face[6] + h_ij_nodal_face[7]);   // h_yy = h_theta_theta
  double h_pp = 0.5 * (h_ij_nodal_face[10] + h_ij_nodal_face[11]); // h_zz = h_phiphi
  double Jc   = 0.5 * (det_h_face[0] + det_h_face[1]);
  (void)lapse_nodal_face; (void)shift_nodal_face;
  (void)Jc;  // not used for eigenvector projection

  double under = h_rr * h_pp - h_rp * h_rp;
  if (!(under > 0.0) || !(h_rr > 0.0) || !(h_tt > 0.0)) {
    return;  // degenerate geometry: identity no-op
  }
  double D_theta = sqrt(under);
  double sqrt_h_rr = sqrt(h_rr);
  double sqrt_h_tt = sqrt(h_tt);
  double inv_sqrt_h_rr = 1.0 / sqrt_h_rr;
  double inv_sqrt_h_tt = 1.0 / sqrt_h_tt;
  double inv_D_theta   = 1.0 / D_theta;

  // M^a_i (y-face tetrad):
  //   row x_tet:  (0,         sqrt(h_tt),  0                )
  //   row y_tet:  (sqrt(h_rr),0,           h_rp/sqrt(h_rr)  )
  //   row z_tet:  (0,         0,           D_theta/sqrt(h_rr))
  double M_xt = sqrt_h_tt;
  double M_yr = sqrt_h_rr;
  double M_yp = h_rp * inv_sqrt_h_rr;
  double M_zp = D_theta * inv_sqrt_h_rr;

  // (M^{-1})^i_a:
  //   row r:  (0,             1/sqrt(h_rr),    -h_rp/(sqrt(h_rr) D_theta))
  //   row t:  (1/sqrt(h_tt),  0,                0                          )
  //   row p:  (0,             0,                sqrt(h_rr)/D_theta         )
  double Minv_rt_inv = inv_sqrt_h_rr;
  double Minv_rp_3   = -h_rp * inv_sqrt_h_rr * inv_D_theta;
  double Minv_tx_2   = inv_sqrt_h_tt;
  double Minv_pp_3   = sqrt_h_rr * inv_D_theta;

  // ---- Build the 3 coord-basis "delta" 6-vectors ----
  double dDx_lc = limiter_fac * (JDx_c[0] - JDx_l[0]);
  double dDy_lc = limiter_fac * (JDy_c[0] - JDy_l[0]);
  double dDz_lc = limiter_fac * (JDz_c[0] - JDz_l[0]);
  double dBx_lc = limiter_fac * (JBx_c[0] - JBx_l[0]);
  double dBy_lc = limiter_fac * (JBy_c[0] - JBy_l[0]);
  double dBz_lc = limiter_fac * (JBz_c[0] - JBz_l[0]);

  // Cell's own theta-slope mode is q[2] (NOT q[1] like the x-limiter).
  double dDx_cc = JDx_c[2];
  double dDy_cc = JDy_c[2];
  double dDz_cc = JDz_c[2];
  double dBx_cc = JBx_c[2];
  double dBy_cc = JBy_c[2];
  double dBz_cc = JBz_c[2];

  double dDx_rc = limiter_fac * (JDx_r[0] - JDx_c[0]);
  double dDy_rc = limiter_fac * (JDy_r[0] - JDy_c[0]);
  double dDz_rc = limiter_fac * (JDz_r[0] - JDz_c[0]);
  double dBx_rc = limiter_fac * (JBx_r[0] - JBx_c[0]);
  double dBy_rc = limiter_fac * (JBy_r[0] - JBy_c[0]);
  double dBz_rc = limiter_fac * (JBz_r[0] - JBz_c[0]);

  // ---- Decompose each delta in tetrad basis (y-face M) ----
#define DECOMP6_Y(dDx, dDy, dDz, dBx, dBy, dBz, am1, am2, a01, a02, ap1, ap2) do { \
    double _dDx_t = M_xt * (dDy); \
    double _dDy_t = M_yr * (dDx) + M_yp * (dDz); \
    double _dDz_t = M_zp * (dDz); \
    double _dBx_t = M_xt * (dBy); \
    double _dBy_t = M_yr * (dBx) + M_yp * (dBz); \
    double _dBz_t = M_zp * (dBz); \
    (am1) = ( -_dDy_t + _dBz_t) * 0.5; \
    (am2) = (  _dDz_t + _dBy_t) * 0.5; \
    (a01) =    _dBx_t; \
    (a02) =    _dDx_t; \
    (ap1) = (  _dDy_t + _dBz_t) * 0.5; \
    (ap2) = ( -_dDz_t + _dBy_t) * 0.5; \
  } while (0)

  double aL_m1, aL_m2, aL_01, aL_02, aL_p1, aL_p2;
  double aC_m1, aC_m2, aC_01, aC_02, aC_p1, aC_p2;
  double aR_m1, aR_m2, aR_01, aR_02, aR_p1, aR_p2;

  DECOMP6_Y(dDx_lc, dDy_lc, dDz_lc, dBx_lc, dBy_lc, dBz_lc,
            aL_m1, aL_m2, aL_01, aL_02, aL_p1, aL_p2);
  DECOMP6_Y(dDx_cc, dDy_cc, dDz_cc, dBx_cc, dBy_cc, dBz_cc,
            aC_m1, aC_m2, aC_01, aC_02, aC_p1, aC_p2);
  DECOMP6_Y(dDx_rc, dDy_rc, dDz_rc, dBx_rc, dBy_rc, dBz_rc,
            aR_m1, aR_m2, aR_01, aR_02, aR_p1, aR_p2);

#undef DECOMP6_Y

#define MINMOD3(a, b, c) ( \
  ((a) > 0.0 && (b) > 0.0 && (c) > 0.0) ?  fmin(fmin((a),(b)),(c)) : \
  ((a) < 0.0 && (b) < 0.0 && (c) < 0.0) ? -fmin(fmin(-(a),-(b)),-(c)) : 0.0 )

  double lim_m1 = MINMOD3(aL_m1, aC_m1, aR_m1);
  double lim_m2 = MINMOD3(aL_m2, aC_m2, aR_m2);
  double lim_01 = MINMOD3(aL_01, aC_01, aR_01);
  double lim_02 = MINMOD3(aL_02, aC_02, aR_02);
  double lim_p1 = MINMOD3(aL_p1, aC_p1, aR_p1);
  double lim_p2 = MINMOD3(aL_p2, aC_p2, aR_p2);

#undef MINMOD3

  int clipped =
    (lim_m1 != aC_m1) || (lim_m2 != aC_m2) ||
    (lim_01 != aC_01) || (lim_02 != aC_02) ||
    (lim_p1 != aC_p1) || (lim_p2 != aC_p2);

  // Reconstruct in tetrad (paper eq. 57 component layout, same as x-kernel):
  double sDx_t = lim_02;
  double sDy_t = -lim_m1 + lim_p1;
  double sDz_t =  lim_m2 - lim_p2;
  double sBx_t = lim_01;
  double sBy_t =  lim_m2 + lim_p2;
  double sBz_t =  lim_m1 + lim_p1;

  // Back-transform to coord basis: sigma^i = (M^{-1})^i_a tilde sigma^a.
  //   row r: (0,           Minv_rt_inv,  Minv_rp_3)
  //   row t: (Minv_tx_2,   0,            0        )
  //   row p: (0,           0,            Minv_pp_3)
  double sDr = Minv_rt_inv * sDy_t + Minv_rp_3 * sDz_t;
  double sDt = Minv_tx_2 * sDx_t;
  double sDp = Minv_pp_3 * sDz_t;
  double sBr = Minv_rt_inv * sBy_t + Minv_rp_3 * sBz_t;
  double sBt = Minv_tx_2 * sBx_t;
  double sBp = Minv_pp_3 * sBz_t;

  // Write into q[2] of each conserved component (coord: r=x, theta=y, phi=z).
  JDx_out[2] = sDr; JDy_out[2] = sDt; JDz_out[2] = sDp;
  JBx_out[2] = sBr; JBy_out[2] = sBt; JBz_out[2] = sBp;

  if (clipped) {
    JDx_out[3] = 0.0; JDy_out[3] = 0.0; JDz_out[3] = 0.0;
    JBx_out[3] = 0.0; JBy_out[3] = 0.0; JBz_out[3] = 0.0;
  }
}
