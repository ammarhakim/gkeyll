#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Tetrad-frame characteristic slope limiter on the radial (x-direction) face,
// 2D serendipity p=1.
//
// Algorithm (per cell, in-place edit of field_con_c slope modes):
//
//   1. Read cell-average representatives (q[0]) for left, center, right cells.
//      Build "delta_l" = limiter_fac * (q_c - q_l) (cell-avg jump from left to
//      center, scaled to slope-mode units via limiter_fac = 1/sqrt(3) by
//      default), "delta_r" similarly, and "delta_c" = q[1]_c (the cell's own
//      r-slope mode, in slope-mode units).
//
//   2. At the cell's right face (already provided in the nodal arrays passed
//      in), average the metric across the 2 face GL nodes to get a single
//      representative h_ij, lapse, shift, det_h at that face. The slope
//      limiter doesn't need per-quadrature-node accuracy -- it's a comparator
//      between cell-averaged jumps.
//
//   3. Transform each of the 3 deltas (in coord basis) to tetrad basis via
//      M^a_i (paper eq. 57 closed form for h_xy = h_yz = 0):
//          M = [[ 1/sqrt(gamma^xx), 0,            0           ],
//               [ 0,                sqrt(h_yy),   0           ],
//               [ h_xz/sqrt(h_zz),  0,            sqrt(h_zz)  ]].
//
//   4. Decompose each tetrad delta into the 6 paper-eq.-57 wave strengths
//      alpha_p, just like the tetrad-Roe surface flux. Three (delta_l,
//      delta_c, delta_r) 6-vectors of characteristic-variable values.
//
//   5. Apply 3-argument minmod component-by-component:
//          alpha_p^lim = minmod(alpha_p^l, alpha_p^c, alpha_p^r)
//      with the convention minmod(a,b,c) = 0 if any signs disagree, else
//      sign(a) * min(|a|, |b|, |c|). Track whether the c-argument changed
//      ("clipped").
//
//   6. Reconstruct limited slope in tetrad: tilde sigma^a = sum_p alpha_p^lim
//      r_p^a. Back-transform to coord basis with M^{-1}.
//
//   7. Write back into field_con_c[1] of each conserved component (the
//      r-slope mode). If ANY of the 6 alphas was clipped, zero the cross
//      slope field_con_c[3] for every conserved component (Option B from
//      design doc Section 5).
//
// Cell-average q[0] is never touched -> the limiter is conservative.
//
// Aux components (phi, psi divergence cleaning) are decoupled from the EM
// system and not affected here; the caller should run a separate aux-limiter
// pass if those modes need bounding.
GKYL_CU_DH void
slope_limiter_x_2x_ser_p1(double limiter_fac, const int theta_pole_lo, const int theta_pole_up,
  const double *lapse_nodal_face, const double *shift_nodal_face,
  const double *h_ij_nodal_face, const double *det_h_face,
  const double *field_con_l, const double *field_con_c, const double *field_con_r,
  double *field_con_out)
{
  // Conserved-variable pointers (J*D and J*B; layout: 6 components x 4 modes).
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

  // Default output = identity on all 4 modes (so q[0], q[2], q[3] pass through
  // unchanged; q[1] gets overwritten below). The caller is allowed to alias
  // field_con_out == field_con_c, in which case these copies are no-ops.
  for (int m = 0; m < 4; ++m) {
    JDx_out[m] = JDx_c[m]; JDy_out[m] = JDy_c[m]; JDz_out[m] = JDz_c[m];
    JBx_out[m] = JBx_c[m]; JBy_out[m] = JBy_c[m]; JBz_out[m] = JBz_c[m];
  }

  // ---- Face-averaged geometry (2 GL nodes -> 1 representative) ----
  // Layout matches the curved-LLF / tetrad-Roe surface kernels: h_ij_nodal
  // is component-major with 2 nodes per face (6 components x 2 = 12 floats).
  double h_yy = 0.5 * (h_ij_nodal_face[6] + h_ij_nodal_face[7]);
  double h_yz = 0.5 * (h_ij_nodal_face[8] + h_ij_nodal_face[9]);
  double h_zz = 0.5 * (h_ij_nodal_face[10] + h_ij_nodal_face[11]);
  double h_xz = 0.5 * (h_ij_nodal_face[4] + h_ij_nodal_face[5]);
  double Jc   = 0.5 * (det_h_face[0] + det_h_face[1]);
  // lapse, shift not used by the limiter itself -- the eigenvectors of the
  // (instantaneous) flux Jacobian in the tetrad don't depend on alpha or
  // beta. Only eigenvalues do, and minmod doesn't use eigenvalues.
  (void)lapse_nodal_face; (void)shift_nodal_face;

  // Suppress unused warnings: h_xx, h_xy not needed for our (x-face,
  // h_xy = h_yz = 0) tetrad construction. h_yz is also 0 in spherical KS but
  // we read it for consistency with the surface flux kernel.
  (void)h_yz;

  // sqrt(gamma^xx) = sqrt(h_yy h_zz - h_yz^2)/Jc -- guard against negative
  // under sqrt (shouldn't happen for valid spatial 3-metrics, but if the
  // limiter is invoked before geometry is settled we want a sane no-op).
  double under = h_yy * h_zz - h_yz * h_yz;
  if (!(under > 0.0) || !(Jc > 0.0)) {
    return;  // identity: no-op slope limiter on degenerate geometry
  }
  double sqrt_gamma_xx     = sqrt(under) / Jc;
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

  // ---- Build the 3 coord-basis "delta" 6-vectors ----
  // delta_l = limiter_fac * (q_c[0] - q_l[0])
  // delta_c = q_c[1] (cell's own r-slope mode, already in slope-mode units)
  // delta_r = limiter_fac * (q_r[0] - q_c[0])
  // BC handling: ghost cells with theta-pole reflective at theta boundary do
  // not occur on the x-face neighbors, so we trust the (l, c, r) cells as
  // valid sources. The updater is responsible for not calling this kernel
  // when l or r is unphysical.
  (void)theta_pole_lo; (void)theta_pole_up;  // unused for x-direction

  double dDx_lc = limiter_fac * (JDx_c[0] - JDx_l[0]);
  double dDy_lc = limiter_fac * (JDy_c[0] - JDy_l[0]);
  double dDz_lc = limiter_fac * (JDz_c[0] - JDz_l[0]);
  double dBx_lc = limiter_fac * (JBx_c[0] - JBx_l[0]);
  double dBy_lc = limiter_fac * (JBy_c[0] - JBy_l[0]);
  double dBz_lc = limiter_fac * (JBz_c[0] - JBz_l[0]);

  double dDx_cc = JDx_c[1];
  double dDy_cc = JDy_c[1];
  double dDz_cc = JDz_c[1];
  double dBx_cc = JBx_c[1];
  double dBy_cc = JBy_c[1];
  double dBz_cc = JBz_c[1];

  double dDx_rc = limiter_fac * (JDx_r[0] - JDx_c[0]);
  double dDy_rc = limiter_fac * (JDy_r[0] - JDy_c[0]);
  double dDz_rc = limiter_fac * (JDz_r[0] - JDz_c[0]);
  double dBx_rc = limiter_fac * (JBx_r[0] - JBx_c[0]);
  double dBy_rc = limiter_fac * (JBy_r[0] - JBy_c[0]);
  double dBz_rc = limiter_fac * (JBz_r[0] - JBz_c[0]);

  // ---- Project each delta to tetrad and decompose into paper-eq.-57 waves ----
  // Helper macro: take a coord-basis (dDx, dDy, dDz, dBx, dBy, dBz) and emit
  // the 6 paper-eq.-57 wave amplitudes (alpha_m1, alpha_m2, alpha_01, alpha_02,
  // alpha_p1, alpha_p2). M is sparse (3-zero off-diagonal) so the transform
  // is closed form.
#define DECOMP6(dDx, dDy, dDz, dBx, dBy, dBz, am1, am2, a01, a02, ap1, ap2) do { \
    double _dDx_t = M_xx * (dDx); \
    double _dDy_t = M_yy * (dDy); \
    double _dDz_t = M_zx * (dDx) + M_zz * (dDz); \
    double _dBx_t = M_xx * (dBx); \
    double _dBy_t = M_yy * (dBy); \
    double _dBz_t = M_zx * (dBx) + M_zz * (dBz); \
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

  DECOMP6(dDx_lc, dDy_lc, dDz_lc, dBx_lc, dBy_lc, dBz_lc,
          aL_m1, aL_m2, aL_01, aL_02, aL_p1, aL_p2);
  DECOMP6(dDx_cc, dDy_cc, dDz_cc, dBx_cc, dBy_cc, dBz_cc,
          aC_m1, aC_m2, aC_01, aC_02, aC_p1, aC_p2);
  DECOMP6(dDx_rc, dDy_rc, dDz_rc, dBx_rc, dBy_rc, dBz_rc,
          aR_m1, aR_m2, aR_01, aR_02, aR_p1, aR_p2);

#undef DECOMP6

  // ---- 3-argument minmod across the 6 wave components ----
  // Returns 0 if any pair has opposite sign, else the smallest |.| with the
  // common sign. Standard generalized minmod.
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

  // Did any wave amplitude get clipped relative to the cell's own slope?
  // Used to decide whether to zero the cross slope q[3].
  int clipped =
    (lim_m1 != aC_m1) || (lim_m2 != aC_m2) ||
    (lim_01 != aC_01) || (lim_02 != aC_02) ||
    (lim_p1 != aC_p1) || (lim_p2 != aC_p2);

  // ---- Reconstruct limited slope in tetrad: tilde sigma^a = sum_p alpha_p r_p^a ----
  // From paper eq. 57:
  //   D_x slot: r_zero_2 entry = +1.
  //   D_y slot: r_minus_1 = -1, r_plus_1 = +1.
  //   D_z slot: r_minus_2 = +1, r_plus_2 = -1.
  //   B_x slot: r_zero_1 = +1.
  //   B_y slot: r_minus_2 = +1, r_plus_2 = +1.
  //   B_z slot: r_minus_1 = +1, r_plus_1 = +1.
  double sDx_t = lim_02;
  double sDy_t = -lim_m1 + lim_p1;
  double sDz_t =  lim_m2 - lim_p2;
  double sBx_t = lim_01;
  double sBy_t =  lim_m2 + lim_p2;
  double sBz_t =  lim_m1 + lim_p1;

  // ---- Back-transform to coord basis: sigma^i = (M^{-1})^i_a tilde sigma^a ----
  // Rows of M^{-1} (only non-zeros listed):
  //   row x:  (Minv_xx, 0,       0      )
  //   row y:  (0,       Minv_yy, 0      )
  //   row z:  (Minv_zx, 0,       Minv_zz)
  double sDx = Minv_xx * sDx_t;
  double sDy = Minv_yy * sDy_t;
  double sDz = Minv_zx * sDx_t + Minv_zz * sDz_t;
  double sBx = Minv_xx * sBx_t;
  double sBy = Minv_yy * sBy_t;
  double sBz = Minv_zx * sBx_t + Minv_zz * sBz_t;

  // ---- Write back limited slopes; zero cross-slope on any clip ----
  JDx_out[1] = sDx; JDy_out[1] = sDy; JDz_out[1] = sDz;
  JBx_out[1] = sBx; JBy_out[1] = sBy; JBz_out[1] = sBz;

  if (clipped) {
    JDx_out[3] = 0.0; JDy_out[3] = 0.0; JDz_out[3] = 0.0;
    JBx_out[3] = 0.0; JBy_out[3] = 0.0; JBz_out[3] = 0.0;
  }
}
