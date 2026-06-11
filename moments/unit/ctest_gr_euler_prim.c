// Primitive-recovery accuracy tests for wv_gr_euler_tetrad (Convention A:
// q[1..3] = √γ·γ_ij·ρhW²·v^j). Split out of ctest_wv_gr_euler_tetrad.c;
// see CTEST_SPLIT_PLAN.md.
//
// How accurately can primitives be reconstructed, across EOS closures
// (IDEAL / TM / RCC) and near floors:
//   (a) round_trip_*       — q → tetrad → q identity + triad invariants
//   (b) prim_consistency_* — curved-frame prim_vars vs SR Newton on q_tet
//   (c) floor_precision_*  — flux-jump residual scaling as τ → floor
//   (d) eos_*              — TM/RCC enthalpy vs the true Synge gas

#include <acutest.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

#include "ctest_gr_euler_test_util.h"

// (a) Round-trip identity: q_curved → q_tet → q_curved should recover the
//     original (modulo √γ factor — q_tet is undensitized in scalars,
//     curved q is densitized). Catches sign/index errors in M, M_inv
//     that flux-jump wouldn't see if they cancel on Σ s·w.
static void
run_round_trip(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  struct gkyl_gr_euler_eos eos)
{
  double gas_gamma = eos.gas_gamma;
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(spacetime, x, y, z, prods);
  if (prods[GKYL_GR_SP_EXCISION] < 0.0) return;

  double rhos[] = { 1.0, 0.3 };
  double pres[] = { 0.5, 2.0 };
  double vels[][3] = {
    { 0.10, 0.20, 0.30 },
    { -0.20, 0.05, 0.0 },
    { 0.0, 0.5, 0.0 },         // pure transverse v
    { 0.0, 0.0, 0.0 },         // static
  };

  // Build triad.
  double g_ij[3][3], inv_g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_ij[i][j]  = prods[GKYL_GR_SP_GIJ + 3*i + j];
      inv_g[i][j] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }
  double M[3][3], M_inv[3][3];
  gkyl_wave_spacetime_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);

  // Sanity: M·M_inv = I, M^T·γ·M = I.
  double max_id_err = 0.0;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double id_ab = M_inv[a][0]*M[0][b] + M_inv[a][1]*M[1][b] + M_inv[a][2]*M[2][b];
      double target = (a == b) ? 1.0 : 0.0;
      double e = fabs(id_ab - target);
      if (e > max_id_err) max_id_err = e;
    }
  }
  TEST_CHECK_( max_id_err < 1.0e-12,
    "[%s @ (%g,%g,%g)] M_inv·M ≠ I: max err = %.3e",
    label, x, y, z, max_id_err );

  // M^T γ M = I (orthonormality).
  double max_orth_err = 0.0;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double inner = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          inner += g_ij[i][j] * M[i][a] * M[j][b];
      double target = (a == b) ? 1.0 : 0.0;
      double e = fabs(inner - target);
      if (e > max_orth_err) max_orth_err = e;
    }
  }
  TEST_CHECK_( max_orth_err < 1.0e-12,
    "[%s @ (%g,%g,%g)] M^T γ M ≠ I: max err = %.3e",
    label, x, y, z, max_orth_err );

  // Round trip on conserved state. The forward+back identity holds at
  // the WAVE level (w_tet round-tripped through back-transform gives
  // back the input w_tet, modulo √γ). For the STATE itself, the
  // forward map undensitizes (q_tet[i] = q_GR[i]/√γ for scalars; with
  // contravariant transform for momentum). Back-transform on a "wave"
  // with all 5 slots equal to q_tet then re-densitizes.
  double max_rt_err = 0.0;
  for (size_t ir = 0; ir < sizeof(rhos)/sizeof(*rhos); ir++) {
    for (size_t ip = 0; ip < sizeof(pres)/sizeof(*pres); ip++) {
      for (size_t iv = 0; iv < sizeof(vels)/sizeof(*vels); iv++) {
        double rho = rhos[ir], p = pres[ip];
        const double *v = vels[iv];

        // Skip superluminal.
        double v_lo[3] = {0,0,0};
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            v_lo[i] += g_ij[i][j] * v[j];
        double vsq = 0.0;
        for (int i = 0; i < 3; i++) vsq += v_lo[i] * v[i];
        if (!(vsq < 1.0 - 1.0e-6)) continue;

        double q[5], q_tet[5], q_back[5];
        build_state_convA(eos, rho, v, p, prods, q);
        gkyl_gr_euler_tetrad_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);
        // The "wave_to_curved" maps tetrad-frame jumps to curved-frame
        // jumps. We use it on q_tet directly (which represents "the
        // full wave from 0 to q_tet"): the result should be the
        // densitized curved state.
        gkyl_gr_euler_tetrad_wave_to_curved_contra(q_tet, sqrt_det, M_inv, q_back);
        for (int i = 0; i < 5; i++) {
          double e = fabs(q_back[i] - q[i]);
          if (e > max_rt_err) max_rt_err = e;
        }
      }
    }
  }
  TEST_CHECK_( max_rt_err < 1.0e-11,
    "[%s @ (%g,%g,%g)] state round-trip max err = %.3e",
    label, x, y, z, max_rt_err );
}

// Each round-trip wrapper loops over [IDEAL, TM, RCC] — primitive
// recovery roundtrip is the most direct test of the EOS dispatch.
void test_round_trip_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_round_trip(st, "Mink", 0.3, 0.0, 0.0, eos_modes[ei]);
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_round_trip(st, "Schw", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_round_trip(st, "Schw", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_round_trip(st, "Schw", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}
void test_round_trip_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_round_trip(st, "Kerr", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_round_trip(st, "Kerr", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_round_trip(st, "Kerr", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}

// (b) Primitive-recovery consistency between curved-frame prim_vars and
//     SR Newton on q_tet. If both Newton solves are correct, they
//     should agree on (ρ, p, W) to machine precision. v_curved^i and
//     v_tet^a are related by the tetrad transform; v_tet^0 should be
//     v^x / √γ^{xx} by our construction.
static void
run_prim_consistency(struct gkyl_gr_spacetime *spacetime,
  const char *label, double x, double y, double z,
  struct gkyl_gr_euler_eos eos)
{
  double gas_gamma = eos.gas_gamma;
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(spacetime, x, y, z, prods);
  if (prods[GKYL_GR_SP_EXCISION] < 0.0) return;

  double g_ij[3][3], inv_g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_ij[i][j]  = prods[GKYL_GR_SP_GIJ + 3*i + j];
      inv_g[i][j] = prods[GKYL_GR_SP_INV_GIJ + 3*i + j];
    }
  double M[3][3], M_inv[3][3];
  gkyl_wave_spacetime_build_triad_contravariant_x(g_ij, inv_g, M, M_inv);
  double sqrt_det = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);

  double rho_in = 1.0, p_in = 1.5;
  double v_in[3] = { 0.12, 0.07, 0.03 };

  // Build state.
  double q[5];
  build_state_convA(eos, rho_in, v_in, p_in, prods, q);

  // Curved-frame recovery via shared helper.
  double D    = q[0] / sqrt_det;
  double momx = q[1] / sqrt_det;
  double momy = q[2] / sqrt_det;
  double momz = q[3] / sqrt_det;
  double tau  = q[4] / sqrt_det;
  struct gkyl_gr_euler_prim prim_curved;
  gkyl_gr_euler_recover_primitives(eos,
    D, momx, momy, momz, tau, inv_g, NULL, &prim_curved);

  // Tetrad-frame recovery: transform to tetrad, then run flat-space
  // Newton (replicate the body of sr_hll_minkowski).
  double q_tet[5];
  gkyl_gr_euler_tetrad_q_to_tetrad_contra(q, sqrt_det, inv_g, M_inv, q_tet);

  double D_t = q_tet[0], Sx_t = q_tet[1], Sy_t = q_tet[2], Sz_t = q_tet[3], tau_t = q_tet[4];
  // Flat-space Banyuls Newton (inline).
  double s_sq_t = ((tau_t + D_t)*(tau_t + D_t))
                - (Sx_t*Sx_t + Sy_t*Sy_t + Sz_t*Sz_t);
  double C, C0;
  if (s_sq_t < pow(10.0, -10.0)) {
    C  = D_t / sqrt(pow(10.0, -10.0));
    C0 = (D_t + tau_t) / sqrt(pow(10.0, -10.0));
  } else {
    C  = D_t / sqrt(s_sq_t);
    C0 = (D_t + tau_t) / sqrt(s_sq_t);
  }
  double a0 = -1.0 / (gas_gamma*gas_gamma);
  double a1 = -2.0 * C * ((gas_gamma - 1.0)/(gas_gamma*gas_gamma));
  double a2 = ((gas_gamma - 2.0)/gas_gamma) * ((C0*C0) - 1.0) + 1.0
            - (C*C)*((gas_gamma-1.0)/gas_gamma)*((gas_gamma-1.0)/gas_gamma);
  double a4 = (C0*C0) - 1.0;
  double eta = 2.0 * C * ((gas_gamma - 1.0)/gas_gamma);
  double g = 1.0;
  for (int it = 0; it < 100; it++) {
    double poly = (a4*g*g*g)*(g-eta) + (a2*g*g) + (a1*g) + a0;
    double dpoly = a1 + 2.0*a2*g + 4.0*a4*g*g*g - 3.0*eta*a4*g*g;
    double gn = g - poly/dpoly;
    if (fabs(g - gn) < pow(10.0, -14.0)) { g = gn; break; }
    g = gn;
  }
  double W_t = 0.5 * C0 * g * (1.0 + sqrt(1.0
    + 4.0*((gas_gamma-1.0)/gas_gamma)*((1.0 - C*g)/(C0*C0*g*g))));
  double h_t = 1.0 / (C * g);
  double rho_t = D_t / W_t;
  double p_t  = (rho_t * h_t * W_t * W_t) - D_t - tau_t;

  // (ρ, p, W, h) should match between curved-frame and tetrad-frame
  // recovery to round-off.
  TEST_CHECK_( fabs(prim_curved.rho - rho_t) < 1.0e-12,
    "[%s] ρ mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.rho, rho_t, fabs(prim_curved.rho - rho_t) );
  TEST_CHECK_( fabs(prim_curved.p - p_t) < 1.0e-12,
    "[%s] p mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.p, p_t, fabs(prim_curved.p - p_t) );
  TEST_CHECK_( fabs(prim_curved.W - W_t) < 1.0e-12,
    "[%s] W mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.W, W_t, fabs(prim_curved.W - W_t) );
  TEST_CHECK_( fabs(prim_curved.h - h_t) < 1.0e-12,
    "[%s] h mismatch curved=%g tetrad=%g, diff=%.3e",
    label, prim_curved.h, h_t, fabs(prim_curved.h - h_t) );

  // The tetrad-frame normal velocity v_tet^x should equal v^x/√γ^{xx}.
  double vx_tet_expected = v_in[0] / sqrt(inv_g[0][0]);
  double vx_tet_computed = Sx_t / (rho_t * h_t * W_t * W_t);
  TEST_CHECK_( fabs(vx_tet_expected - vx_tet_computed) < 1.0e-12,
    "[%s] v_tet^x: expected v^x/√γ^{xx} = %g, got %g (diff %.3e)",
    label, vx_tet_expected, vx_tet_computed,
    fabs(vx_tet_expected - vx_tet_computed) );
}

// Each prim-consistency wrapper loops over [IDEAL, TM, RCC].
void test_prim_consistency_minkowski(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++)
    run_prim_consistency(st, "Mink", 0.3, 0.0, 0.0, eos_modes[ei]);
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_schwarzschild(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_prim_consistency(st, "Schw @ (0.3,0.2,0)", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Schw @ (0.5,0.0,0)", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Schw @ (0.4,0.4,0)", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}
void test_prim_consistency_kerr(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  for (int ei = 0; ei < NUM_EOS_MODES; ei++) {
    run_prim_consistency(st, "Kerr @ (0.3,0.2,0)", 0.3, 0.2, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Kerr @ (0.5,0.0,0)", 0.5, 0.0, 0.0, eos_modes[ei]);
    run_prim_consistency(st, "Kerr @ (0.4,0.4,0)", 0.4, 0.4, 0.0, eos_modes[ei]);
  }
  gkyl_gr_spacetime_release(st);
}

// (d) Precision sweep: how does flux-jump residual scale with τ on the
//     low side? Report-only — diagnostic of catastrophic cancellation in
//     the Banyuls Newton when τ is near the admissibility boundary. Tells
//     us what the "safe" repair floor is from a finite-precision
//     standpoint — below this value, fluxes degrade and feed back into
//     larger downstream errors.
static void
run_floor_precision_sweep(struct gkyl_gr_spacetime *spacetime,
  enum gkyl_wv_gr_euler_tetrad_rp rp, const char *rp_name)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_euler_eos eos = gkyl_gr_euler_eos_ideal(gas_gamma);
  struct gkyl_range conf_range;
  int lower[1] = { 0 }, upper[1] = { 0 };
  gkyl_range_init(&conf_range, 1, lower, upper);
  struct gkyl_wv_eqn *eqn = make_eqn(eos, conf_range, rp);
  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, conf_range.volume);
  gkyl_gr_euler_tetrad_set_auxfields(eqn,
    (struct gkyl_wv_gr_euler_tetrad_auxfields){ .prods = prods });
  double *prods_row = gkyl_array_fetch(prods, 0);
  fill_prods_at(spacetime, 0.3, 0.2, 0.0, prods_row);
  if (prods_row[GKYL_GR_SP_EXCISION] < 0.0) {
    gkyl_array_release(prods); gkyl_wv_eqn_release(eqn);
    return;
  }

  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);

  double norm[3] = { 1.0, 0.0, 0.0 };
  double tau1[3] = { 0.0, 1.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };
  int idx[1] = { 0 };

  // Sweep pL (which sets τ_L in the cold limit). Other side is a typical
  // state. Same per-side primitives otherwise.
  double pL_values[] = {
    1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10,
    1.0e-12, 1.0e-14, 1.0e-16
  };
  double rho_L = 0.01;
  double v_L[3] = { 0.05, 0.10, 0.0 };
  double rho_R = 1.0, p_R = 0.5;
  double v_R[3] = { 0.10, 0.20, 0.0 };

  fprintf(stderr,
    "[floor-precision %s] p_L      |Σs·w − ΔF|_max  |amdq+apdq − ΔF|_max\n",
    rp_name);
  for (size_t i = 0; i < sizeof(pL_values)/sizeof(*pL_values); i++) {
    double p_L = pL_values[i];
    double qL_glob[5], qR_glob[5];
    build_state_convA(eos, rho_L, v_L, p_L, prods_row, qL_glob);
    build_state_convA(eos, rho_R, v_R, p_R, prods_row, qR_glob);

    eqn->set_interface_idx_func(eqn, idx, idx);
    double qL[5], qR[5];
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qL_glob, qL);
    eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qR_glob, qR);

    double delta[5];
    for (int k = 0; k < 5; k++) delta[k] = qR[k] - qL[k];
    double waves[3 * 5], speeds[3];
    eqn->waves_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      delta, qL, qR, 1.0, 1.0, waves, speeds);

    double dF[5];
    banyuls_delta_flux(eos, grm, qL, qR, dF);

    double max_fj = 0.0, max_fb = 0.0;
    int nw = eqn->num_waves;
    for (int k = 0; k < 5; k++) {
      double sw = 0.0;
      for (int w = 0; w < nw; w++) sw += speeds[w] * waves[w * 5 + k];
      double r = fabs(sw - dF[k]);
      if (r > max_fj) max_fj = r;
    }
    double amdq[5], apdq[5];
    eqn->qfluct_func(eqn, GKYL_WV_HIGH_ORDER_FLUX,
      qL, qR, 1.0, 1.0, waves, speeds, amdq, apdq);
    for (int k = 0; k < 5; k++) {
      double r = fabs(amdq[k] + apdq[k] - dF[k]);
      if (r > max_fb) max_fb = r;
    }

    fprintf(stderr, "  %-8s   %.0e  %.3e         %.3e\n",
      rp_name, p_L, max_fj, max_fb);
  }

  gkyl_array_release(prods);
  gkyl_wv_eqn_release(eqn);
}

void test_floor_precision_hllc(void)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_HLLC, "HLLC");
  gkyl_gr_spacetime_release(st);
}
void test_floor_precision_lax(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_LAX, "Schw-Lax");
  gkyl_gr_spacetime_release(st);
}
void test_floor_precision_hll(void)
{
  struct gkyl_gr_spacetime *st =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_floor_precision_sweep(st, WV_GR_EULER_TETRAD_RP_HLL, "Schw-HLL");
  gkyl_gr_spacetime_release(st);
}

// ---------------------------------------------------------------------------
// EOS-accuracy validation vs the true Synge gas.
//
// The Synge gas (single-component perfect relativistic gas) has the
// closed-form specific enthalpy
//   h_Synge(θ) = K3(1/θ) / K2(1/θ),     θ = p/ρ.
// where K_n is the modified Bessel function of the second kind. Both
// MATHEWS_TAUB (TM) and RYU_CHATTOPADHYAY (RCC) closures are designed
// to approximate this curve in the trans-relativistic regime
// θ ∈ [O(0.1), O(10)] with stated maximum relative errors:
//   TM (Mignone+ 2005)         ~ 2%   on h
//   RCC (Ryu+ 2006)            ~ 0.8% on h
//
// This test sweeps θ ∈ [0.1, 10], runs a primitive recovery roundtrip
// at each point under each EOS, then compares the EOS's enthalpy at
// the recovered primitives against the Synge value. We also assert
// the roundtrip recovery is accurate (the recovery should exactly
// invert the EOS-specific build_state_convA seeding to machine
// precision for the chosen EOS).
//
// Coverage:
//   TM no-flow (W = 1): tests the closed-form TM cubic recovery in
//                       the static-fluid M = 0 short-circuit.
//   RCC relativistic   (|v| = 0.9, W ≈ 2.29): tests the TM cubic
//                       warm-start → RC Newton refinement → paper-
//                       grounded (a)/(b)/(c) physicality check path.
// ---------------------------------------------------------------------------

// Modified Bessel functions K_0, K_1 via the Abramowitz & Stegun
// polynomial fits (9.8.1-9.8.8). Accurate to <2e-7 over all x > 0
// — far better than our test tolerances of 2% (TM) and 0.8% (RCC).
// libm provides these on Linux/glibc as XSI extensions but Apple's
// libSystem does NOT, so we implement them locally for portability.
static double bessel_i0(double x)
{
  if (x < 0.0) x = -x;
  if (x < 3.75) {
    double t = (x / 3.75) * (x / 3.75);
    return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
      + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
  }
  // x ≥ 3.75 — asymptotic form. Not needed for our K_0 path (we only
  // call I_0 from K_0's "small x" branch where x ≤ 2 < 3.75), but kept
  // for completeness.
  double t = 3.75 / x;
  return (exp(x) / sqrt(x)) * (0.39894228 + t * (0.01328592
    + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281
    + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633
    + t * 0.00392377))))))));
}
static double bessel_i1(double x)
{
  double ax = (x < 0.0) ? -x : x;
  double res;
  if (ax < 3.75) {
    double t = (ax / 3.75) * (ax / 3.75);
    res = ax * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934
      + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
  } else {
    double t = 3.75 / ax;
    res = (exp(ax) / sqrt(ax)) * (0.39894228 + t * (-0.03988024
      + t * (-0.00362018 + t * (0.00163801 + t * (-0.01031555
      + t * (0.02282967 + t * (-0.02895312 + t * (0.01787654
      + t * (-0.00420059)))))))));
  }
  return (x < 0.0) ? -res : res;
}
static double bessel_k0(double x)
{
  if (x <= 2.0) {
    double t = (x / 2.0) * (x / 2.0);
    return (-log(x / 2.0) * bessel_i0(x)) + (-0.57721566
      + t * (0.42278420 + t * (0.23069756 + t * (0.03488590
      + t * (0.00262698 + t * (0.00010750 + t * 0.00000740))))));
  }
  double t = 2.0 / x;
  return (exp(-x) / sqrt(x)) * (1.25331414 + t * (-0.07832358
    + t * (0.02189568 + t * (-0.01062446 + t * (0.00587872
    + t * (-0.00251540 + t * 0.00053208))))));
}
static double bessel_k1(double x)
{
  if (x <= 2.0) {
    double t = (x / 2.0) * (x / 2.0);
    return (log(x / 2.0) * bessel_i1(x)) + (1.0 / x) * (1.0
      + t * (0.15443144 + t * (-0.67278579 + t * (-0.18156897
      + t * (-0.01919402 + t * (-0.00110404 + t * (-0.00004686))))))) ;
  }
  double t = 2.0 / x;
  return (exp(-x) / sqrt(x)) * (1.25331414 + t * (0.23498619
    + t * (-0.03655620 + t * (0.01504268 + t * (-0.00780353
    + t * (0.00325614 + t * (-0.00068245)))))));
}

// K_2, K_3 via forward recurrence K_{n+1}(x) = K_{n-1}(x) + (2n/x)·K_n(x).
// Stable in the forward direction for K_n (unlike I_n).
//
// Numerical range used here: θ ∈ [0.1, 10] ⇒ x = 1/θ ∈ [0.1, 10].
//   - At x = 10: K_n(10) is small (~1e-5) but well above underflow;
//     the ratio K3/K2 is well-conditioned.
//   - At x = 0.1: K_n(0.1) is large (K0 ~ 2.4, K3 ~ 8e3); still finite,
//     ratio K3/K2 → 4/x = 4θ (ultra-relativistic limit).
// Both extremes evaluate cleanly in double precision.
static double bessel_k2_arg(double x)
{
  return bessel_k0(x) + (2.0 / x) * bessel_k1(x);
}
static double bessel_k3_arg(double x)
{
  return bessel_k1(x) + (4.0 / x) * bessel_k2_arg(x);
}
static double synge_enthalpy(double theta)
{
  double x = 1.0 / theta;
  return bessel_k3_arg(x) / bessel_k2_arg(x);
}

// Sweep θ ∈ [0.1, 10] (16 logspaced points) and validate that:
//   (a) primitive recovery roundtrip preserves (ρ, v, p) to round_trip_tol
//   (b) the EOS-specific enthalpy at the recovered (ρ, p) matches the
//       true Synge gas enthalpy to within max_rel_err.
static void
run_eos_synge_validation(struct gkyl_gr_euler_eos eos,
  const double v_co[3], double round_trip_tol, double max_rel_err,
  const char *label)
{
  struct gkyl_gr_spacetime *st = gkyl_gr_minkowski_new(false);
  double prods[GKYL_GR_SP_NCOMP_BASE];
  fill_prods_at(st, 0.3, 0.0, 0.0, prods);

  // Sanity-check velocity is subluminal at this metric (Minkowski → just
  // need |v|² < 1).
  double v_lo[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_lo[i] += prods[GKYL_GR_SP_GIJ + 3*i + j] * v_co[j];
  double vsq = 0.0;
  for (int i = 0; i < 3; i++) vsq += v_lo[i] * v_co[i];
  TEST_CHECK_( vsq < 1.0,
    "[%s] velocity is superluminal (γ_ij v^i v^j = %g)", label, vsq );

  // Logspaced θ scan over [0.1, 10] (16 points).
  static const double thetas[] = {
    0.10, 0.135, 0.18, 0.24, 0.32, 0.42, 0.56, 0.75,
    1.00, 1.33,  1.78, 2.37, 3.16, 4.22, 5.62, 10.0
  };
  const int n_theta = sizeof(thetas) / sizeof(*thetas);

  double max_rt_dr = 0.0, max_rt_dp = 0.0, max_rt_dv = 0.0;
  double max_h_err = 0.0;
  double max_h_err_theta = 0.0;

  for (int it = 0; it < n_theta; it++) {
    double theta = thetas[it];
    double rho = 1.0;
    double p = theta;  // → p/ρ = θ

    // Seed the conservative state using the EOS-specific enthalpy.
    double q[5];
    build_state_convA(eos, rho, v_co, p, prods, q);

    // Run the production recovery dispatch.
    double sd  = sqrt(prods[GKYL_GR_SP_SPATIAL_DET]);
    double D   = q[0] / sd;
    double Sx  = q[1] / sd;
    double Sy  = q[2] / sd;
    double Sz  = q[3] / sd;
    double tau = q[4] / sd;
    const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
    double inv_g[3][3] = {
      { ig[0], ig[1], ig[2] },
      { ig[3], ig[4], ig[5] },
      { ig[6], ig[7], ig[8] },
    };
    struct gkyl_gr_euler_prim prim;
    gkyl_gr_euler_recover_primitives(eos, D, Sx, Sy, Sz, tau, inv_g,
      NULL, &prim);

    // (a) Roundtrip: recovered primitives should match the inputs.
    double dr = fabs(prim.rho - rho);
    double dp = fabs(prim.p - p);
    double dv = 0.0;
    for (int i = 0; i < 3; i++) {
      double d = fabs(prim.v[i] - v_co[i]);
      if (d > dv) dv = d;
    }
    if (dr > max_rt_dr) max_rt_dr = dr;
    if (dp > max_rt_dp) max_rt_dp = dp;
    if (dv > max_rt_dv) max_rt_dv = dv;

    // (b) EOS enthalpy at recovered primitives vs Synge enthalpy at
    // the recovered θ. Using recovered values keeps the comparison
    // self-consistent if (a) has any FP slack.
    double theta_r = prim.p / prim.rho;
    double h_eos   = gkyl_gr_euler_eos_enthalpy(eos, prim.rho, prim.p);
    double h_synge = synge_enthalpy(theta_r);
    double rel_err = fabs(h_eos - h_synge) / h_synge;
    if (rel_err > max_h_err) {
      max_h_err = rel_err;
      max_h_err_theta = theta;
    }
  }

  TEST_CHECK_( max_rt_dr < round_trip_tol,
    "[%s] roundtrip ρ residual %.3e exceeds %.3e", label, max_rt_dr, round_trip_tol );
  TEST_CHECK_( max_rt_dp < round_trip_tol,
    "[%s] roundtrip p residual %.3e exceeds %.3e", label, max_rt_dp, round_trip_tol );
  TEST_CHECK_( max_rt_dv < round_trip_tol,
    "[%s] roundtrip v residual %.3e exceeds %.3e", label, max_rt_dv, round_trip_tol );
  TEST_CHECK_( max_h_err < max_rel_err,
    "[%s] h_EOS vs h_Synge max rel err %.3e exceeds %.3e (worst at θ=%.3g)",
    label, max_h_err, max_rel_err, max_h_err_theta );

  // Informational dump of the worst-case errors (visible only on fail
  // or via TEST_MSG always-shown channel for postmortem).
  TEST_MSG( "[%s] roundtrip residuals max: ρ=%.3e p=%.3e v=%.3e | "
            "max |h_EOS − h_Synge|/h_Synge = %.3e at θ=%.3g",
    label, max_rt_dr, max_rt_dp, max_rt_dv, max_h_err, max_h_err_theta );

  gkyl_gr_spacetime_release(st);
}

// TM (use_rcc=false) with no flow (v = 0, W = 1). Exercises the
// closed-form TM cubic recovery's M = 0 static-fluid short-circuit.
// Tolerance 2.1% rather than the literature 2%: a direct numerical
// sweep |h_TM − h_Synge|/h_Synge over θ ∈ [0.1, 10] peaks at ~2.011%
// near θ ≈ 0.52. Mignone+ 2005 quotes "≲ 2%" as a rounded-down
// description of this same peak; the extra 0.1% margin in the test
// captures the actual worst-case the bound is approximating.
void test_eos_tm_synge_no_flow(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = false };
  double v[3] = { 0.0, 0.0, 0.0 };
  run_eos_synge_validation(eos, v, 1.0e-12, 2.1e-2, "TM no-flow");
}

// RCC (use_rcc=true) with relativistic flow (|v|² = 0.81, W ≈ 2.29).
// Exercises the TM cubic warm-start → RC Newton refinement → paper-
// grounded physicality check path. Tolerance 0.8% per Ryu+ 2006.
void test_eos_rcc_synge_relativistic(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  // Components chosen so |v|² = (0.9)² with v isotropic across x/y/z.
  double v_each = 0.9 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-10, 8.0e-3, "RCC relativistic");
}

// RCC at higher Lorentz factor (|v| = 0.99 → W ≈ 7.09). Stresses the
// pressure recovery at large W where the (E²−M²)/E² ~ 1/W² ≈ 0.02
// conditioning measure starts to bite — historically the regime where
// TM cubic precision-loss appeared. RC Newton on the squared
// polynomial is supposed to converge robustly here without the TM
// precision-loss feedback loop. Round-trip and enthalpy tolerances
// stay tight (RCC's Synge accuracy is an enthalpy-formula property,
// independent of W; precision-loss would show up first in the
// roundtrip ρ/p residual).
void test_eos_rcc_synge_high_W(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  double v_each = 0.99 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-10, 8.0e-3, "RCC high-W");
}

// RCC at extreme Lorentz factor (|v| = 0.9999 → W ≈ 70.71). The
// (E²−M²)/E² ~ 1/W² ≈ 2e-4 conditioning measure is now squarely in
// the regime where the TM cubic's coefficients lose ~3 digits to
// cancellation. SESSION_NOTES_EOS_IMPROVEMENTS §2.2 catalogs cases
// where TM gave W = 5×10⁵ and RC corrected it to W = 50; this test
// ensures RC Newton's refinement + paper-grounded physicality check
// still locks onto the right root at extreme W.
//
// Round-trip tolerance loosened to 1e-6 to reflect realistic FP
// precision at extreme W:
//   - ρ residual scales as ULP × (D/ρ ratio) ~ 1e-12 × W ~ 7e-11
//     compounded with the W-solve precision; measured ≈ 7e-8.
//   - p residual is worse because p = ρhW² − τ − ρW is a catastrophic
//     cancellation: |ρhW²| ≈ |τ| ≈ ρ·4θ·W² ~ 2×10⁴ at the worst θ,
//     while p ≈ θ ~ 1, so p inherits relative precision degraded by
//     ρhW²/p ≈ 2×10⁴; measured ≈ 5e-7.
// What this test is really checking: RC Newton at extreme W does NOT
// diverge or fall back to TM with W = 5×10⁵ (a documented failure
// mode at this regime). Enthalpy accuracy stays at 0.8% — that's an
// algebraic property of the recovered θ and unaffected by W.
void test_eos_rcc_synge_extreme_W(void)
{
  struct gkyl_gr_euler_eos eos = {
    .type = GR_EULER_EOS_APPROXIMATE_SYNGE, .use_rcc = true };
  double v_each = 0.9999 / sqrt(3.0);
  double v[3] = { v_each, v_each, v_each };
  run_eos_synge_validation(eos, v, 1.0e-6, 8.0e-3, "RCC extreme-W");
}

TEST_LIST = {
  { "round_trip_minkowski",     test_round_trip_minkowski },
  { "round_trip_schwarzschild", test_round_trip_schwarzschild },
  { "round_trip_kerr",          test_round_trip_kerr },

  { "prim_consistency_minkowski",     test_prim_consistency_minkowski },
  { "prim_consistency_schwarzschild", test_prim_consistency_schwarzschild },
  { "prim_consistency_kerr",          test_prim_consistency_kerr },

  { "floor_precision_lax",  test_floor_precision_lax },
  { "floor_precision_hll",  test_floor_precision_hll },
  { "floor_precision_hllc", test_floor_precision_hllc },

  // EOS accuracy vs the true Synge gas (Bessel-function reference).
  { "eos_tm_synge_no_flow",          test_eos_tm_synge_no_flow },
  { "eos_rcc_synge_relativistic",    test_eos_rcc_synge_relativistic },
  { "eos_rcc_synge_high_W",          test_eos_rcc_synge_high_W },
  { "eos_rcc_synge_extreme_W",       test_eos_rcc_synge_extreme_W },

  { NULL, NULL },
};
