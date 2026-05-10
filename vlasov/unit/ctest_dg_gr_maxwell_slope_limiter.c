// Tetrad-frame characteristic slope limiter unit tests for the DG GR Maxwell
// scheme. Tests cover the design-doc Section 8.4 list:
//
//   1. flat_minmod_sanity      — flat metric reduces to standard component-
//                                wise minmod in characteristic variables.
//   2. zero_jump_no_change     — when all 3 deltas (left, center, right) are
//                                zero, the limiter leaves the cell untouched.
//   3. sign_disagreement_zeros — opposite-sign deltas zero the limited slope.
//   4. smooth_field_preserved  — linear field across 3 cells: slope unchanged.
//   5. conservation_q0         — q[0] (cell average) is never modified.
//   6. clip_zeros_q3           — clipping any wave amplitude zeros q[3].
//   7. curved_metric_smooth    — Kerr-Schild interior: smooth linear field's
//                                slope is preserved through the M /M^{-1}
//                                round-trip.
//   8. y_direction_smooth      — y-kernel: smooth field in theta preserved.
//   9. y_direction_theta_pole  — y-kernel at pole adjacency zeros q[2], q[3].
//  10. direction_decoupling    — x-pass doesn't touch q[2]; y-pass doesn't
//                                touch q[1].

#include <acutest.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_dg_gr_maxwell_kernels.h>

// ---- Test setup helpers ----

// Build a flat (Cartesian) face nodal geometry: h_ij = identity at both GL
// nodes, lapse = 1, shift = 0, det_h = 1. Layout matches the surface flux
// kernels: h_ij_nodal is component-major (xx, xy, xz, yy, yz, zz) with 2 GL
// nodes per component.
static void
set_flat_face_nodal(double *lapse, double *shift, double *h_ij, double *det_h)
{
  for (int n = 0; n < 2; ++n) {
    lapse[n] = 1.0;
    shift[n + 0] = 0.0; shift[n + 2] = 0.0; shift[n + 4] = 0.0;
    h_ij[n + 0] = 1.0;  h_ij[n + 2] = 0.0;  h_ij[n + 4] = 0.0;
    h_ij[n + 6] = 1.0;  h_ij[n + 8] = 0.0;  h_ij[n + 10] = 1.0;
    det_h[n] = 1.0;
  }
}

// Build a Kerr-Schild face geometry with given (h_rr, h_rp, h_tt, h_pp).
// Both face GL nodes share the same metric (we don't test sub-face variation
// here -- the limiter averages anyway).
static void
set_curved_face_nodal(double *lapse, double *shift, double *h_ij, double *det_h,
  double h_rr, double h_rp, double h_tt, double h_pp,
  double alpha, double beta_r)
{
  double det = h_rr * h_tt * h_pp - h_tt * h_rp * h_rp;
  double Jc = sqrt(det);
  for (int n = 0; n < 2; ++n) {
    lapse[n] = alpha;
    shift[n + 0] = beta_r; shift[n + 2] = 0.0; shift[n + 4] = 0.0;
    h_ij[n + 0] = h_rr;  h_ij[n + 2] = 0.0;  h_ij[n + 4] = h_rp;
    h_ij[n + 6] = h_tt;  h_ij[n + 8] = 0.0;  h_ij[n + 10] = h_pp;
    det_h[n] = Jc;
  }
}

// Initialize the 3 cells (24 floats each: 6 components x 4 modes) with all
// zeros, then let the caller set specific entries.
static void
zero_cells(double *fl, double *fc, double *fr)
{
  for (int k = 0; k < 24; ++k) { fl[k] = 0.0; fc[k] = 0.0; fr[k] = 0.0; }
}

// Set q[mode] of conserved-component `comp` (0..5: Dx Dy Dz Bx By Bz) to value.
// For each cell, the layout is component-major: comp 0 fills indices 0..3,
// comp 1 fills 4..7, etc.
static void
set_mode(double *cell, int comp, int mode, double value)
{
  cell[4*comp + mode] = value;
}

static double
get_mode(const double *cell, int comp, int mode)
{
  return cell[4*comp + mode];
}

// ---- Test 1: flat-metric minmod sanity ----
// In flat Cartesian, M = I and the tetrad transform is trivial; minmod on
// the 6 paper-eq.-57 waves reduces to minmod on coord-basis differences. We
// set up a single conserved component (D^z) with non-trivial q[1] and
// compatible cell-avg jumps to verify the standard minmod selection.
static void
test_flat_minmod_sanity(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // D^z cell averages: 1.0, 2.0, 4.0 (concave-up — slope increasing).
  // Cell's own r-slope on D^z: 0.5 (set in q[1] of comp 2 = D^z).
  set_mode(fl, /*Dz*/ 2, 0, 1.0);
  set_mode(fc, /*Dz*/ 2, 0, 2.0);
  set_mode(fr, /*Dz*/ 2, 0, 4.0);
  set_mode(fc, /*Dz*/ 2, 1, 0.5);  // own r-slope

  double fout[24];
  for (int k = 0; k < 24; ++k) fout[k] = fc[k];

  double limiter_fac = 0.5773502691896258;  // 1/sqrt(3)
  slope_limiter_x_2x_ser_p1(limiter_fac, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // Hand calc:
  //   ΔU_L = limiter_fac * (2 - 1) = 1/sqrt(3) ~ 0.5774
  //   σ_c = 0.5
  //   ΔU_R = limiter_fac * (4 - 2) = 2/sqrt(3) ~ 1.1547
  // All positive, so minmod picks the smallest in magnitude = 0.5 (σ_c).
  // For flat metric M = I, the D^z slot decomposes via waves r_minus_2 and
  // r_plus_2:
  //   alpha_m2 = (dDz + dBy)/2 = (dDz)/2  (dBy = 0 here)
  //   alpha_p2 = (-dDz + dBy)/2 = -dDz/2
  // Limited reconstructs as: D^z slot = alpha_m2 - alpha_p2 = dDz.
  // So minmod on alpha_m2 and alpha_p2 each selects the smallest |.| with
  // the common sign — alpha_m2 has same sign across all 3 (positive), alpha_p2
  // has same sign (negative) — both keep their σ_c values, reconstruct = 0.5.
  double q1_lim = get_mode(fout, /*Dz*/ 2, 1);
  double diff = fabs(q1_lim - 0.5);
  if (diff > 1e-12) {
    printf("  flat_minmod: expected 0.5, got %e (diff %e)\n", q1_lim, diff);
  }
  TEST_CHECK(diff < 1e-12);

  // q[0] preserved.
  TEST_CHECK(get_mode(fout, /*Dz*/ 2, 0) == 2.0);
}

// ---- Test 2: zero jump -> no change to limited slope ----
// If both ΔU_L and ΔU_R are zero, minmod sees argument 0 -> sign disagreement
// triggers a zero out. So the cell's own slope σ_c is clipped to zero.
// This is a feature of the limiter (cell is locally an extremum so its slope
// should be zeroed); we verify the kernel behaves this way.
static void
test_zero_jump_no_change(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // All cell avgs equal, but cell has a non-zero slope (local extremum).
  set_mode(fl, /*Dz*/ 2, 0, 5.0);
  set_mode(fc, /*Dz*/ 2, 0, 5.0);
  set_mode(fr, /*Dz*/ 2, 0, 5.0);
  set_mode(fc, /*Dz*/ 2, 1, 0.3);

  double fout[24];
  double limiter_fac = 0.5773502691896258;
  slope_limiter_x_2x_ser_p1(limiter_fac, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // ΔU_L = ΔU_R = 0, so minmod(0, σ_c, 0) = 0 -- sign disagreement triggers.
  double q1_lim = get_mode(fout, /*Dz*/ 2, 1);
  TEST_CHECK(fabs(q1_lim) < 1e-14);
  TEST_CHECK(get_mode(fout, /*Dz*/ 2, 0) == 5.0);  // q[0] preserved
}

// ---- Test 3: sign disagreement zeros slope ----
// Cell-avg jumps with opposite signs: left says "decreasing", right says
// "increasing". Limiter should zero the slope.
static void
test_sign_disagreement_zeros(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // Cell-avg pattern: 3, 1, 3 (valley) -- left jump negative, right positive.
  set_mode(fl, /*Dy*/ 1, 0, 3.0);
  set_mode(fc, /*Dy*/ 1, 0, 1.0);
  set_mode(fr, /*Dy*/ 1, 0, 3.0);
  set_mode(fc, /*Dy*/ 1, 1, 1.5);  // arbitrary non-zero slope

  double fout[24];
  double limiter_fac = 0.5773502691896258;
  slope_limiter_x_2x_ser_p1(limiter_fac, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // ΔU_L = (1/sqrt(3)) * (1-3) = -2/sqrt(3), ΔU_R = +2/sqrt(3). Sign
  // disagreement -> minmod -> 0 for every wave touching D_y.
  double q1_lim = get_mode(fout, /*Dy*/ 1, 1);
  TEST_CHECK(fabs(q1_lim) < 1e-14);
}

// ---- Test 4: smooth field preserves slope ----
// If the linear extrapolation matches the cell-avg jumps exactly, all three
// minmod arguments are equal and the slope is preserved.
static void
test_smooth_field_preserved(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // Linear field: q0_l = A - b, q0_c = A, q0_r = A + b. Cell's r-slope
  // mode is then q[1] = b/sqrt(3) (matches limiter_fac * b).
  double A = 5.0, b = 2.0;
  set_mode(fl, /*Dx*/ 0, 0, A - b);
  set_mode(fc, /*Dx*/ 0, 0, A);
  set_mode(fr, /*Dx*/ 0, 0, A + b);
  double limiter_fac = 0.5773502691896258;  // 1/sqrt(3)
  set_mode(fc, /*Dx*/ 0, 1, limiter_fac * b);

  double fout[24];
  slope_limiter_x_2x_ser_p1(limiter_fac, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // ΔU_L = limiter_fac * b, σ_c = limiter_fac * b, ΔU_R = limiter_fac * b
  // All identical -> minmod returns common value -> slope unchanged.
  double q1_lim = get_mode(fout, /*Dx*/ 0, 1);
  double diff = fabs(q1_lim - limiter_fac * b);
  if (diff > 1e-12) {
    printf("  smooth: expected %e, got %e (diff %e)\n", limiter_fac*b, q1_lim, diff);
  }
  TEST_CHECK(diff < 1e-12);
}

// ---- Test 5: conservation of q[0] ----
// q[0] must NEVER be modified by the limiter, regardless of how aggressive
// it is. Set wildly mismatched data and verify q[0] passes through.
static void
test_conservation_q0(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);
  // Random-ish nonzero data in q[0] for all 6 components and the slope modes.
  for (int comp = 0; comp < 6; ++comp) {
    set_mode(fl, comp, 0, 1.0 + 0.3 * comp);
    set_mode(fc, comp, 0, 2.0 - 0.2 * comp);
    set_mode(fr, comp, 0, 0.5 + 0.7 * comp);
    set_mode(fc, comp, 1, 1.5 * (comp - 2));
    set_mode(fc, comp, 2, -0.4);
    set_mode(fc, comp, 3, 0.9);
  }

  double fout[24];
  slope_limiter_x_2x_ser_p1(0.5773502691896258, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);
  for (int comp = 0; comp < 6; ++comp) {
    TEST_CHECK(get_mode(fout, comp, 0) == get_mode(fc, comp, 0));
  }
  // q[2] also not touched by the x-direction limiter.
  for (int comp = 0; comp < 6; ++comp) {
    TEST_CHECK(get_mode(fout, comp, 2) == get_mode(fc, comp, 2));
  }
}

// ---- Test 6: clipping zeros q[3] ----
// When any wave gets clipped, the cross-slope q[3] must be zeroed across
// all 6 conserved components (Option B from design Section 5).
static void
test_clip_zeros_q3(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // Set up the same sign-disagreement scenario as test 3 (valley shape),
  // but in conjunction with non-zero q[3] in all components.
  for (int comp = 0; comp < 6; ++comp) {
    set_mode(fl, comp, 0, 3.0);
    set_mode(fc, comp, 0, 1.0);
    set_mode(fr, comp, 0, 3.0);
    set_mode(fc, comp, 1, 1.0);
    set_mode(fc, comp, 3, 0.7);  // pre-clip cross slope value
  }

  double fout[24];
  slope_limiter_x_2x_ser_p1(0.5773502691896258, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // Some waves will be clipped (sign disagreement in dDx, dDy, etc.).
  // Verify q[3] zeroed for all 6 components.
  for (int comp = 0; comp < 6; ++comp) {
    double q3 = get_mode(fout, comp, 3);
    if (fabs(q3) > 1e-14) {
      printf("  clip_zeros_q3: comp=%d q3=%e (expected 0)\n", comp, q3);
    }
    TEST_CHECK(fabs(q3) < 1e-14);
  }
}

// ---- Test 7: curved-metric smooth field preserved through round-trip ----
// In Kerr-Schild interior, set up a smooth linear field; the limiter applies
// M (forward), minmod (identity since all 3 args equal), M^{-1} (back). The
// final slope should match the original within the round-trip tolerance.
static void
test_curved_metric_smooth(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  // Kerr-Schild interior: representative metric components in r=θ-natural
  // ordering. Numbers chosen to look like a mid-interior point.
  set_curved_face_nodal(lapse, shift, h_ij, det_h,
    /*h_rr=*/1.4, /*h_rp=*/-0.3, /*h_tt=*/4.0, /*h_pp=*/3.6,
    /*alpha=*/0.8, /*beta_r=*/0.2);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  // Linear field on one component: D^y (mid-component). q0 = A + B*x.
  double A = 5.0, b = 2.0;
  set_mode(fl, /*Dy*/ 1, 0, A - b);
  set_mode(fc, /*Dy*/ 1, 0, A);
  set_mode(fr, /*Dy*/ 1, 0, A + b);
  double limiter_fac = 0.5773502691896258;
  set_mode(fc, /*Dy*/ 1, 1, limiter_fac * b);

  double fout[24];
  slope_limiter_x_2x_ser_p1(limiter_fac, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // ΔU_L = ΔU_R = σ_c (in coord basis); after M-forward, minmod, M^{-1},
  // the result should equal σ_c. This tests that the M / M^{-1} round-trip
  // is exact (which we already separately verified in ctest_dg_gr_maxwell_tetrad
  // tests M_satisfies_metric_and_inverse). Limited slope should equal
  // limiter_fac * b on D^y.
  double q1_lim = get_mode(fout, /*Dy*/ 1, 1);
  double expected = limiter_fac * b;
  double diff = fabs(q1_lim - expected);
  if (diff > 1e-11) {
    printf("  curved_smooth: expected %e, got %e (diff %e)\n",
      expected, q1_lim, diff);
  }
  TEST_CHECK(diff < 1e-11);

  // Also verify all OTHER conserved components have q[1] = 0 (no slope
  // information was put into them, M is orthogonal-ish so they shouldn't pick
  // up spurious slope from D^y's data — at least not for the diagonal-only
  // case where dDy maps cleanly to tetrad-y-component slot).
  // For our M with h_rp != 0, dDz tetrad pulls from BOTH dDx and dDz, so
  // putting only dDy data should leave dDz_t = 0 and similarly the dDz coord
  // output = 0. Let's verify D^x and D^z slopes stay zero (numerically).
  TEST_CHECK(fabs(get_mode(fout, /*Dx*/ 0, 1)) < 1e-11);
  TEST_CHECK(fabs(get_mode(fout, /*Dz*/ 2, 1)) < 1e-11);
}

// ---- Test 8: y-direction smooth field preserved ----
static void
test_y_direction_smooth(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_curved_face_nodal(lapse, shift, h_ij, det_h,
    /*h_rr=*/1.4, /*h_rp=*/-0.3, /*h_tt=*/4.0, /*h_pp=*/3.6,
    /*alpha=*/0.8, /*beta_r=*/0.2);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);

  double A = 5.0, b = 1.5;
  set_mode(fl, /*Dx*/ 0, 0, A - b);
  set_mode(fc, /*Dx*/ 0, 0, A);
  set_mode(fr, /*Dx*/ 0, 0, A + b);
  double limiter_fac = 0.5773502691896258;
  // For y-direction, the cell's own slope mode is q[2].
  set_mode(fc, /*Dx*/ 0, 2, limiter_fac * b);

  double fout[24];
  slope_limiter_y_2x_ser_p1(limiter_fac, /*pole_lo=*/0, /*pole_up=*/0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  double q2_lim = get_mode(fout, /*Dx*/ 0, 2);
  double expected = limiter_fac * b;
  double diff = fabs(q2_lim - expected);
  if (diff > 1e-11) {
    printf("  y_smooth: expected %e, got %e (diff %e)\n", expected, q2_lim, diff);
  }
  TEST_CHECK(diff < 1e-11);
}

// ---- Test 9: y-direction theta-pole adjacency zeros q[2], q[3] ----
static void
test_y_direction_theta_pole(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  // Near-pole metric (h_tt small but positive; the limiter doesn't actually
  // do the M math when theta_pole flag is set).
  set_curved_face_nodal(lapse, shift, h_ij, det_h,
    /*h_rr=*/1.5, /*h_rp=*/0.0, /*h_tt=*/0.01, /*h_pp=*/0.01,
    /*alpha=*/0.7, /*beta_r=*/0.1);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);
  for (int comp = 0; comp < 6; ++comp) {
    set_mode(fc, comp, 0, 1.0);
    set_mode(fc, comp, 1, 0.5);
    set_mode(fc, comp, 2, 0.4);
    set_mode(fc, comp, 3, 0.3);
  }

  double fout[24];
  slope_limiter_y_2x_ser_p1(0.5773502691896258, /*pole_lo=*/1, /*pole_up=*/0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // Pole-adjacent: q[2] and q[3] zeroed; q[0] and q[1] preserved.
  for (int comp = 0; comp < 6; ++comp) {
    TEST_CHECK(fabs(get_mode(fout, comp, 2)) < 1e-14);
    TEST_CHECK(fabs(get_mode(fout, comp, 3)) < 1e-14);
    TEST_CHECK(get_mode(fout, comp, 0) == 1.0);
    TEST_CHECK(get_mode(fout, comp, 1) == 0.5);
  }
}

// ---- Test 10: x-pass doesn't touch q[2]; y-pass doesn't touch q[1] ----
static void
test_direction_decoupling(void)
{
  double lapse[2], shift[6], h_ij[12], det_h[2];
  set_flat_face_nodal(lapse, shift, h_ij, det_h);

  double fl[24], fc[24], fr[24];
  zero_cells(fl, fc, fr);
  // Set up some sign-disagreement so the x-limiter actually fires.
  for (int comp = 0; comp < 6; ++comp) {
    set_mode(fl, comp, 0, 3.0);
    set_mode(fc, comp, 0, 1.0);
    set_mode(fr, comp, 0, 3.0);
    set_mode(fc, comp, 1, 1.0);
    set_mode(fc, comp, 2, 0.7);  // pre-existing theta-slope
  }

  double fout[24];
  slope_limiter_x_2x_ser_p1(0.5773502691896258, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);

  // q[2] passed through unchanged on every component.
  for (int comp = 0; comp < 6; ++comp) {
    TEST_CHECK(get_mode(fout, comp, 2) == 0.7);
  }

  // Now reset and test y-limiter doesn't touch q[1].
  for (int comp = 0; comp < 6; ++comp) {
    set_mode(fc, comp, 1, 0.7);
    set_mode(fc, comp, 2, 1.0);
  }
  slope_limiter_y_2x_ser_p1(0.5773502691896258, 0, 0,
    lapse, shift, h_ij, det_h, fl, fc, fr, fout);
  for (int comp = 0; comp < 6; ++comp) {
    TEST_CHECK(get_mode(fout, comp, 1) == 0.7);
  }
}

TEST_LIST = {
  { "flat_minmod_sanity",       test_flat_minmod_sanity      },
  { "zero_jump_no_change",      test_zero_jump_no_change     },
  { "sign_disagreement_zeros",  test_sign_disagreement_zeros },
  { "smooth_field_preserved",   test_smooth_field_preserved  },
  { "conservation_q0",          test_conservation_q0         },
  { "clip_zeros_q3",            test_clip_zeros_q3           },
  { "curved_metric_smooth",     test_curved_metric_smooth    },
  { "y_direction_smooth",       test_y_direction_smooth      },
  { "y_direction_theta_pole",   test_y_direction_theta_pole  },
  { "direction_decoupling",     test_direction_decoupling    },
  { NULL, NULL },
};
