#pragma once

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

// Shared primitive-variable recovery + admissibility check + iterative
// state-repair helpers for the modular tetrad GR Euler equation
// (wv_gr_euler_tetrad_mod.c) and its source-term integrator
// (moment_spacetime_coupling.c). All call sites route the Banyuls Newton
// solve through gkyl_gr_euler_recover_primitives so the Convention-A
// bookkeeping (q[1..3] is genuine covariant momentum S_i; v^i is raised
// via γ^{ij}) lives in one place.
//
// Inputs are *undensitized* — the caller divides the packed conservatives
// by √γ before calling.

// Pressure / density / τ-repair floors used by the Newton recovery and
// the state-repair step. The three floors are kept equal so the
// regularized cold-flow state (ρ=floor, p=floor, τ≈floor) is internally
// consistent (τ = ρhW² − p − ρW ≈ 0 at floor values).
//
// Set to 1e-10 (down from a historical 1e-8). Reasoning:
//   - For BHL with high-Lorentz cells (W ~ 100), the dimensionless
//     pressure-precision quantity is p/(ρW²). At W²=10⁴ and ρ~1, a
//     floor of 1e-10 gives p/(ρW²) = 1e-14, still 2 digits above
//     machine epsilon (~1e-16). The 1e-8 floor was over-clamping cells
//     with physically meaningful p in [1e-10, 1e-8].
//   - For TM/RC cubic recovery, the pressure formula p = (e²−ρ²)/(3e)
//     suffers cancellation at e ≈ ρ. A lower floor lets the recovered
//     primitives carry more physical information about marginally-warm
//     cells (post-shock cooling, downstream wake) instead of artificially
//     declaring them "vacuum-cold-flow" at the 1e-8 boundary.
//   - Trade-off: less margin against re-triggering admissibility
//     violations on the next flux step. SESSION_NOTES record that 1e-14
//     gave 54× more wave_prop τ<0 fires than 1e-8 (insufficient margin
//     vs per-step flux scale). 1e-10 sits between — 2 digits below
//     the previously-validated 1e-8, 4 digits above the previously-
//     broken 1e-14.
//
// The source-step τ-limiter target follows automatically via
// TAU_TARGET = GR_EULER_TAU_REPAIR_FLOOR in moment_spacetime_coupling.c.
#define GR_EULER_PRESSURE_FLOOR   (1.0e-10)
#define GR_EULER_DENSITY_FLOOR    (1.0e-10)
#define GR_EULER_TAU_REPAIR_FLOOR (1.0e-10)

// ---------------------------------------------------------------------------
// Equation-of-state abstraction.
//
// All EOS-dependent math (specific enthalpy, sound speed, pressure recovery
// in the Newton iteration) routes through the three eos_* helpers below.
// The struct gkyl_gr_euler_eos is a POD bundle stored on the equation
// object and passed by value through the primitive solve, the SR Riemann
// cores, the flux helpers, and the source integrator.
//
//   IDEAL              — p = (γ-1)·ρ·ε  ⇒  h = 1 + γ/(γ-1)·p/ρ
//   MATHEWS_TAUB       — (h - θ)(h - 4θ) = 1, θ = p/ρ. Convex Taub-Mathews
//                        closure that interpolates between the non-rel
//                        monatomic limit (Γ=5/3, c_s²=5p/(3ρ)) and the
//                        ultra-rel limit (Γ=4/3, c_s²=1/3). Approximates
//                        the single-component perfect gas (Synge / RP) to
//                        within ~2% on enthalpy.
//   RYU_CHATTOPADHYAY  — h = 2(6θ² + 4θ + 1)/(3θ + 2). Ryu+ 2006 RC EOS.
//                        Same asymptotic limits as MATHEWS_TAUB but a
//                        tighter fit to the Synge gas (~0.8% on enthalpy
//                        vs MATHEWS_TAUB's ~2%). Recovery requires Newton
//                        iteration on a degree-8 polynomial (Ryu+ eq 29);
//                        not closed-form like the TM cubic.
// ---------------------------------------------------------------------------

enum gkyl_gr_euler_eos_type {
  GR_EULER_EOS_IDEAL = 0,
  GR_EULER_EOS_MATHEWS_TAUB = 1,
  GR_EULER_EOS_RYU_CHATTOPADHYAY = 2,
};

struct gkyl_gr_euler_eos {
  enum gkyl_gr_euler_eos_type type;
  double gas_gamma;  // adiabatic index (consulted only when type == IDEAL)
};

// Specific enthalpy h(ρ, p).
//   IDEAL:             h = 1 + γ/(γ-1)·p/ρ
//   MATHEWS_TAUB:      h(θ) = (5θ + sqrt(9θ² + 4))/2
//   RYU_CHATTOPADHYAY: h(θ) = 2(6θ² + 4θ + 1)/(3θ + 2)             [Ryu+ eq 15]
// where θ = p/ρ.
static inline double
gkyl_gr_euler_eos_enthalpy(struct gkyl_gr_euler_eos eos,
  double rho, double p)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return 1.0 + (p / rho) * (eos.gas_gamma / (eos.gas_gamma - 1.0));
  }
  double theta = p / rho;
  if (eos.type == GR_EULER_EOS_RYU_CHATTOPADHYAY) {
    return 2.0 * (6.0*theta*theta + 4.0*theta + 1.0) / (3.0*theta + 2.0);
  }
  return 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));
}

// Pressure recovery p(ρ, h). The EOS-specific inverse of h(ρ, p), used
// in primitive-variable recovery and in cell-interface state reconstruction.
//   IDEAL:             h = 1 + γ/(γ-1)·θ  ⇒  θ = (h-1)·(γ-1)/γ
//   MATHEWS_TAUB:      4θ² − 5hθ + (h² − 1) = 0
//                      ⇒  θ = (5h − sqrt(9h² + 16))/8
//   RYU_CHATTOPADHYAY: (3θ+2)·h = 12θ² + 8θ + 2
//                      ⇒  12θ² + (8-3h)θ + (2-2h) = 0
//                      ⇒  θ = (3h − 8 + sqrt(9h² + 48h − 32))/24      [Ryu+ eq 60]
static inline double
gkyl_gr_euler_eos_pressure_from_rho_h(struct gkyl_gr_euler_eos eos,
  double rho, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return rho * (h - 1.0) * (eos.gas_gamma - 1.0) / eos.gas_gamma;
  }
  if (eos.type == GR_EULER_EOS_RYU_CHATTOPADHYAY) {
    double theta = (3.0*h - 8.0 + sqrt(9.0*h*h + 48.0*h - 32.0)) / 24.0;
    return rho * theta;
  }
  double theta = 0.125 * (5.0 * h - sqrt(9.0 * h * h + 16.0));
  return rho * theta;
}

// Sound speed squared c_s²(ρ, p, h). Cold-flow limit (θ→0) is 5p/(3ρh)
// (γ=5/3 ideal); ultra-rel limit (θ→∞) is 1/3 (radiation fluid).
//   IDEAL:             c_s² = γp/(ρh)
//   MATHEWS_TAUB:      c_s² = θ(5h − 8θ) / (3·h·(h − θ))
//   RYU_CHATTOPADHYAY: c_s² = θ(3θ+2)(18θ²+24θ+5)
//                            / (3·(6θ²+4θ+1)·(9θ²+12θ+2))            [Ryu+ eq 16]
static inline double
gkyl_gr_euler_eos_cs2(struct gkyl_gr_euler_eos eos,
  double rho, double p, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return (eos.gas_gamma * p) / (rho * h);
  }
  double theta = p / rho;
  if (eos.type == GR_EULER_EOS_RYU_CHATTOPADHYAY) {
    double num = theta * (3.0*theta + 2.0) * (18.0*theta*theta + 24.0*theta + 5.0);
    double den = 3.0 * (6.0*theta*theta + 4.0*theta + 1.0)
                     * (9.0*theta*theta + 12.0*theta + 2.0);
    return num / den;
  }
  double num = theta * (5.0 * h - 8.0 * theta);
  double den = 3.0 * h * (h - theta);
  return num / den;
}

// Convenience for the legacy IDEAL-only call sites (packed tetrad path,
// older tests) that still hand in a bare gas_gamma. Constructs an IDEAL
// eos bundle inline.
static inline struct gkyl_gr_euler_eos
gkyl_gr_euler_eos_ideal(double gas_gamma)
{
  return (struct gkyl_gr_euler_eos){
    .type = GR_EULER_EOS_IDEAL,
    .gas_gamma = gas_gamma,
  };
}

// TM γ_eff matching: given recovered (ρ, p) from an EM Newton pass at
// some trial γ, returns the γ_eff such that h_IDEAL(γ_eff)(ρ, p) equals
// h_TM(ρ, p). This is the Picard fixed-point update γ_TM(EM(γ)) used by
// the TM branch of gkyl_gr_euler_recover_primitives.
//
//   θ = p/ρ
//   h_TM(θ) = (5θ + √(9θ² + 4)) / 2
//   f       = (h_TM − 1) / θ                  (f ∈ [5/2, 4) for TM)
//   γ_match = f / (f − 1)                     (γ_match ∈ (4/3, 5/3))
//
// θ → 0 (cold-flow) is the 0/0 limit where γ_match → 5/3 exactly; we
// short-circuit at θ < 1e-12 to avoid numerical noise. Density and
// pressure floors guard against negative inputs from a barely-recovered
// state.
static inline double
gkyl_gr_euler_tm_match_gamma_eff(double rho, double p)
{
  double rho_safe = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
  double p_safe   = (p   > GR_EULER_PRESSURE_FLOOR) ? p   : GR_EULER_PRESSURE_FLOOR;
  double theta = p_safe / rho_safe;
  if (theta < 1.0e-12) return 5.0 / 3.0;
  double h_tm = 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));
  double f = (h_tm - 1.0) / theta;
  return f / (f - 1.0);
}

// First-failing-constraint enum reported by check_admissibility. The
// repair_state callback uses this to decide which single projection to
// apply this iteration.
enum gkyl_gr_euler_admissibility_status {
  GR_EULER_ADM_OK = 0,
  GR_EULER_ADM_BAD_D,    // D ≤ 0
  GR_EULER_ADM_BAD_TAU,  // τ < 0
  GR_EULER_ADM_BAD_S2,   // s² = (D+τ)² − γ^{ij}·S_i·S_j ≤ 0
};

// Bitmask returned by the cascade-repair helper indicating which
// constraints were fixed in a single call. A single repair invocation
// may fix multiple constraints at once when the order is unambiguous
// (D and τ are independent; S² depends on the post-floor D+τ).
enum gkyl_gr_euler_repair_flags {
  GR_EULER_REPAIR_NONE = 0,
  GR_EULER_REPAIR_D    = 1u << 0,
  GR_EULER_REPAIR_TAU  = 1u << 1,
  GR_EULER_REPAIR_S2   = 1u << 2,
};

// Compute the Lorentz scalar |S|² = γ^{ij}·S_i·S_j with S covariant.
static inline double
gkyl_gr_euler_mom_sq(const double inv_g[3][3],
  double Sx, double Sy, double Sz)
{
  return inv_g[0][0]*Sx*Sx
       + inv_g[1][1]*Sy*Sy
       + inv_g[2][2]*Sz*Sz
       + 2.0*(inv_g[0][1]*Sx*Sy
            + inv_g[0][2]*Sx*Sz
            + inv_g[1][2]*Sy*Sz);
}

// Strict admissibility check for the Banyuls inversion.
//   D > 0
//   τ ≥ 0    (energy ≥ rest-mass: D + τ = ρhW² ≥ ρW = D, with h ≥ 1, W ≥ 1)
//   s² > 0
// Returns the FIRST failing constraint (in this order) so repair_state can
// project minimally on a single axis per iteration.
static inline enum gkyl_gr_euler_admissibility_status
gkyl_gr_euler_check_admissibility(
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3])
{
  if (!(D > 0.0))   return GR_EULER_ADM_BAD_D;
  if (!(tau >= 0.0)) return GR_EULER_ADM_BAD_TAU;
  double s_sq = ((D + tau) * (D + tau)) - gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  if (!(s_sq > 0.0)) return GR_EULER_ADM_BAD_S2;
  return GR_EULER_ADM_OK;
}

// Cascade state repair: in a single call, check each of the three
// admissibility constraints (D > 0, τ ≥ 0, s² > 0) and apply a minimal
// projection for each that is violated. Order is fixed and safe:
//   1. D-floor (independent of τ, S)
//   2. τ-floor (independent of D, S)
//   3. S-rescale (uses the post-floor D + τ for the bound)
//
// Because the order respects the dependencies, the post-call state is
// always admissible — no caller iteration needed.
//
// Repairs:
//   D ≤ 0   → D ← GR_EULER_DENSITY_FLOOR
//   τ < 0   → τ ← GR_EULER_PRESSURE_FLOOR. Cold-flow (W=1, h=1) inversion
//             of τ = ρhW² − p − D gives τ → p_floor, the safe lower
//             bound that keeps Banyuls recoverable. Newton then converges
//             to a physical (W, ρ, p) with p ≈ p_floor.
//   s² ≤ 0  → rescale (Sx, Sy, Sz) so |S|² = (1 − ε)·(D+τ)² with
//             ε = 1e-8, putting s² strictly in the interior of the
//             convex set so Newton converges away from the singularity.
//
// Returns a bitmask of GR_EULER_REPAIR_* flags indicating which
// constraints were fixed (zero means input was already admissible).
static inline unsigned int
gkyl_gr_euler_repair_admissibility_cascade(
  const double inv_g[3][3],
  double *D, double *Sx, double *Sy, double *Sz, double *tau)
{
  unsigned int fixed = GR_EULER_REPAIR_NONE;

  if (!(*D > 0.0)) {
    *D = GR_EULER_DENSITY_FLOOR;
    fixed |= GR_EULER_REPAIR_D;
  }

  if (!(*tau >= 0.0)) {
    // Land at GR_EULER_TAU_REPAIR_FLOOR (above p_floor) so the cell has
    // safe margin against the next flux step pushing it across zero.
    *tau = GR_EULER_TAU_REPAIR_FLOOR;
    fixed |= GR_EULER_REPAIR_TAU;
  }

  // S² check uses the (possibly updated) D and τ.
  double Dt = *D + *tau;
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, *Sx, *Sy, *Sz);
  double s_sq = (Dt * Dt) - mom_sq;
  if (!(s_sq > 0.0)) {
    // Repair into the INTERIOR of A_γ with a finite margin, so the next
    // flux step has breathing room before the cell touches the s²
    // boundary again. The original margin (1e-8) left cells essentially
    // ON the boundary — any subsequent Lax/HLL dissipation immediately
    // re-triggered the repair cascade. BHL production sweep results
    // (t_end=3.0, M=0.3 BH):
    //   1e-8 (orig) + 2x amax: 81 wave_prop s²<0 fires, 39 s wall
    //   1e-6           + 1x amax: 12 fires, 31.5 s wall  (best ratio)
    //   1e-4           + 1x amax: 12 fires, identical
    //   1e-2           + 1x amax: 14 fires (chaotic divergence)
    // 1e-6 sits at the plateau. Larger margin doesn't help; smaller
    // gets immediately overwhelmed by the next step's Lax dissipation.
    const double margin = 1.0e-6;
    if (mom_sq > 0.0) {
      double target = (1.0 - margin) * Dt * Dt;
      double scale = sqrt(target / mom_sq);
      *Sx *= scale;
      *Sy *= scale;
      *Sz *= scale;
    } else {
      // mom_sq is exactly zero (e.g., all S components zero) and s² ≤ 0
      // would require D + τ ≤ 0 — which the τ pass above already
      // handled. Defensive no-op.
      *Sx = *Sy = *Sz = 0.0;
    }
    fixed |= GR_EULER_REPAIR_S2;
  }

  return fixed;
}

struct gkyl_gr_euler_prim {
  double rho;     // rest-mass density (post-floor)
  double v[3];    // contravariant 3-velocity v^i
  double p;       // pressure (post-floor)
  double W;       // Lorentz factor
  double h;       // specific enthalpy (EOS-dependent closure)
  bool admissible; // true iff input lay in the strict admissibility set
                   // (D > 0, τ ≥ 0, s² > 0). check_inv consumes this; the
                   // floors below do NOT downgrade it.
};

// ---------------------------------------------------------------------------
// Recovery-iteration instrumentation.
//
// gkyl_gr_euler_recover_primitives accumulates per-call iteration counts for
// both the inner EM Newton (within a single Picard pass) and the outer
// Picard loop (TM only — IDEAL is one Picard pass by definition). The
// histograms tell us where the wall-clock cost lives: see
// SESSION_NOTES_PHASE_1 for the bucket definitions and analysis.
//
// All recording is process-global; the counters survive equation-object
// teardown and get printed by gkyl_gr_euler_tetrad_mod_free alongside the
// existing repair-cascade stats. The instrumentation function bodies live
// in moment_spacetime_coupling.c so the inline recover_primitives can
// remain header-only.
// ---------------------------------------------------------------------------

// Bucket boundaries for the inner-Newton histogram (upper edge, inclusive).
// Bin k counts calls with iter count ≤ GR_EULER_NEWTON_BIN_EDGES[k] and
// > GR_EULER_NEWTON_BIN_EDGES[k-1]. The final bin is "≥ last edge" overflow.
#define GR_EULER_NEWTON_HIST_NBINS 7
// 0–4, 5–8, 9–16, 17–32, 33–64, 65–99, 100+
extern const int gkyl_gr_euler_newton_bin_edges[GR_EULER_NEWTON_HIST_NBINS];

// Bucket boundaries for the outer-Picard histogram. Brent cap is 30, so
// the histogram splits the previous "20+" overflow into 20-29 and 30+
// to expose whether non-converging cells are dominantly hitting the cap
// or just slowly winding down through the bracket.
#define GR_EULER_PICARD_HIST_NBINS 9
// 1, 2, 3, 4, 5–9, 10–14, 15–19, 20–29, 30+
extern const int gkyl_gr_euler_picard_bin_edges[GR_EULER_PICARD_HIST_NBINS];

// Record one recovery call's inner-Newton iteration count. Called from
// gkyl_gr_euler_em_newton_at_gamma at the end of its iteration loop.
void gkyl_gr_euler_record_newton_iters(int n);

// Record one recovery call's outer-Picard iteration count. Called from
// gkyl_gr_euler_recover_primitives at the end of the Picard loop. For
// IDEAL (no Picard), this is invoked with n=1.
void gkyl_gr_euler_record_picard_iters(int n);

// Print the histograms (and totals) to stderr. Invoked by
// gkyl_gr_euler_tetrad_mod_free alongside repair-cascade stats.
void gkyl_gr_euler_print_recovery_stats(FILE *fp);

// Bounded snapshot of one EM-Newton-cap-hit event. Captured by the
// recovery helper when iter_done reaches 100 (no convergence). Logged in
// a process-global ring buffer so the post-run dump in
// gkyl_gr_euler_print_recovery_stats can show what the pathological
// states look like (which regime: cold/hot/near-floor/near-supraluminal).
struct gkyl_gr_euler_newton_capfail {
  double gas_gamma;          // γ being tried (after Picard / cold-restart logic)
  double D, Sx, Sy, Sz, tau; // undensitized conservatives
  double mom_sq;             // γ^{ij}·S_i·S_j
  double inv_g_diag[3];      // γ^{xx}, γ^{yy}, γ^{zz} (curvature indicator)
  double inv_g_off[3];       // γ^{xy}, γ^{xz}, γ^{yz}
  double rho, p, W, h;       // final (possibly non-converged) primitives
  double theta;              // p / ρ
  double s_sq_over_Dtau_sq;  // s²/(D+τ)² — how close to supraluminal
};

// Record one Newton-cap event. Stored in a bounded ring buffer (we cap
// the per-run number of snapshots so a buggy state doesn't flood log).
void gkyl_gr_euler_record_newton_capfail(
  const struct gkyl_gr_euler_newton_capfail *event);

// RC EOS comparison diagnostic. Called from the RC dispatch in
// gkyl_gr_euler_recover_primitives once per non-cold-flow cell. Records:
//   - rc_attempts:        total non-cold-flow cells routed through RC
//                         Newton
//   - rc_accepted:        cells where RC's recovered Γ passed the
//                         sanity check and replaced the TM result.
//   - rc_accepted_escape: cells accepted via the CONDITIONING ESCAPE
//                         (when TM's cubic is precision-suspect due to
//                         (E²−M²)/E² ≪ 1, i.e., near-luminal cells).
//                         Tracked separately so we can see how often
//                         the escape fires and whether it correlates
//                         with downstream physics improvements.
//   - dW_rel:             |W_rc − W_tm| / W_tm at this cell. Accumulated
//                         on accepted cells; max also tracked on
//                         rejected cells (those where neither the
//                         standard 10% threshold nor the conditioning
//                         escape rescued RC).
//
// Print via gkyl_gr_euler_print_recovery_stats alongside the other
// recovery counters.
void gkyl_gr_euler_record_rc_compare(
  bool accepted, bool conditioning_escape, double dW_rel);

// Bounded snapshot of a single RC-vs-TM disagreement event (rejected by
// the sanity check). Used to diagnose WHY RC was rejected:
//   - Is W_tm large? Suggests TM is losing precision from (E²−M²)²
//     coefficient cancellation — RC might actually be more accurate and
//     we'd be throwing away good refinement.
//   - Is W_tm small (~1-3) with reasonable conservatives? Suggests
//     Newton converged to a spurious root of the squared polynomial —
//     rejection is correct.
//   - Is θ_tm near the cold-flow boundary (just above 1e-4)? Suggests
//     RC's iterative scheme is failing near the regime where its
//     pressure formula starts to lose precision.
struct gkyl_gr_euler_rc_reject {
  double dW_rel;         // |W_rc − W_tm| / W_tm (>0.10 by construction)
  double W_tm, W_rc;     // Lorentz factors from both schemes
  double D, tau, mom_sq; // undensitized conservatives
  double s_sq_over_E2;   // (E_lab² − M²)/E_lab² — coefficient conditioning measure
  double theta_tm;       // p_tm / ρ_tm — EOS regime
  double p_tm, p_rc;     // pressures from each scheme
};

// Record one RC-rejection event. Stored in a bounded ring buffer (32
// entries) so a pathological run doesn't flood log; total count is
// preserved as s_rc_attempts − s_rc_accepted.
void gkyl_gr_euler_record_rc_reject(
  const struct gkyl_gr_euler_rc_reject *event);

// Recovery-call context — set by callers via gkyl_gr_euler_set_recovery_context
// before each gkyl_gr_euler_recover_primitives invocation so the floor-hit
// diagnostic can tag which call site is producing inadmissible cells.
// Crucial for distinguishing "real" floor hits (source step, prim_vars
// hook) from "transient" ones (HLL Riemann attempts that get rejected by
// the positivity sweep and replaced by Lax — these contribute false
// positives to the total floor count because the cells are corrected
// before being committed).
enum gkyl_gr_euler_recovery_context {
  GR_EULER_CTX_UNKNOWN = 0,
  GR_EULER_CTX_SOURCE  = 1,  // source-step recovery (committed)
  GR_EULER_CTX_PRIMS   = 2,  // prim_vars hook (committed)
  GR_EULER_CTX_HLL     = 3,  // HLL Riemann — may be rejected by positivity sweep
  GR_EULER_CTX_LAX     = 4,  // Lax Riemann — positivity-preserving fallback
  GR_EULER_CTX_HLLC    = 5,  // HLLC Riemann
  GR_EULER_CTX_COUNT   = 6,
};
void gkyl_gr_euler_set_recovery_context(int ctx);

// Floor-hit diagnostic. Records how often the outer recover_primitives
// floor clamps either ρ or p (separately tracked). Useful for floor-
// value sweeps (e.g., 1e-8 vs 1e-10 vs 1e-12) to see whether lowering
// the floor is "doing work" (catching cells previously clamped) or
// "irrelevant" (no cells reach the lower floor under the new value).
//
// Conservatives and raw recovered values are also passed in so the
// diagnostic can build a histogram of where the floored cells live in
// (D, mom_sq, E²−M²/E²) space — answers "are these all near-vacuum
// cells or do normal-density cells also get floored?"
//
// Note: this counts cells where the raw recovered value at the outer
// recover_primitives entry was BELOW the floor. It does NOT count the
// internal e_rest ≤ floor early-return in the TM/RC pressure formulas,
// since those return p = FLOOR exactly (which evaluates as "not below
// floor" at the outer check). For the floor-sweep experiment, the outer
// count is the primary signal.
void gkyl_gr_euler_record_floor_hit(
  bool rho_floored, bool p_floored,
  double D, double tau, double mom_sq,
  double rho_raw, double p_raw);

// One Eulderink-Mellema quartic Newton pass at fixed γ. Recovers
// (ρ, v^i, p, W, h) from undensitized (D, S_i, τ) and inv_g.
//
// The "guess" iterate is 1/(h·C), C = D/√s², which gives a 4th-order
// polynomial residual that converges quadratically from guess=1.0 for
// any admissible input. Outputs are pre-floor so callers can apply
// floors / dispatch after refinement.
//
// Returns the iteration count: < 100 = converged; 100 = hit cap. Picard
// callers use this to detect a bad warm-start γ_eff (rare, but seen in
// the BHL TM run when a cell's state changed enough between timesteps to
// invalidate the cached γ_eff without triggering cascade-repair).
static inline int
gkyl_gr_euler_em_newton_at_gamma(
  double gas_gamma, double mom_sq,
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  double *rho_out, double *vx_out, double *vy_out, double *vz_out,
  double *p_out, double *W_out, double *h_out)
{
  double s_sq = ((tau + D) * (tau + D)) - mom_sq;

  double C, C0;
  if (s_sq < pow(10.0, -10.0)) {
    C  = D / sqrt(pow(10.0, -10.0));
    C0 = (D + tau) / sqrt(pow(10.0, -10.0));
  } else {
    C  = D / sqrt(s_sq);
    C0 = (D + tau) / sqrt(s_sq);
  }

  double alpha0 = -1.0 / (gas_gamma * gas_gamma);
  double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
  double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0
                - (C*C) * ((gas_gamma - 1.0) / gas_gamma)
                        * ((gas_gamma - 1.0) / gas_gamma);
  double alpha4 = (C0*C0) - 1.0;
  double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);

  double guess = 1.0;
  int iter_done = 100;  // count actual iterations (1-100); 100 = no convergence
  for (int iter = 0; iter < 100; iter++) {
    double g2 = guess * guess;
    double g3 = g2 * guess;
    double poly = (alpha4 * g3 * (guess - eta))
                + (alpha2 * g2) + (alpha1 * guess) + alpha0;
    double poly_der = alpha1 + (2.0 * alpha2 * guess)
                    + (4.0 * alpha4 * g3) - (3.0 * eta * alpha4 * g2);
    double guess_new = guess - (poly / poly_der);
    if (fabs(guess - guess_new) < pow(10.0, -14.0)) {
      guess = guess_new;
      iter_done = iter + 1;
      break;
    }
    guess = guess_new;
  }
  gkyl_gr_euler_record_newton_iters(iter_done);

  double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0
    + (4.0 * ((gas_gamma - 1.0) / gas_gamma)
          * ((1.0 - (C * guess)) / ((C0*C0) * (guess*guess))))));
  double h = 1.0 / (C * guess);

  double rho = D / W;
  double rhohW2 = rho * h * W * W;
  double vx = (inv_g[0][0]*Sx + inv_g[0][1]*Sy + inv_g[0][2]*Sz) / rhohW2;
  double vy = (inv_g[1][0]*Sx + inv_g[1][1]*Sy + inv_g[1][2]*Sz) / rhohW2;
  double vz = (inv_g[2][0]*Sx + inv_g[2][1]*Sy + inv_g[2][2]*Sz) / rhohW2;
  double p  = rhohW2 - D - tau;

  *rho_out = rho; *vx_out = vx; *vy_out = vy; *vz_out = vz;
  *p_out = p; *W_out = W; *h_out = h;

  // Record cap-hit cells for postmortem. Cheap when the bin doesn't fire
  // (almost always); when it does, we snapshot the state for the post-run
  // diagnostic dump. The ring buffer caps total events so a bad state
  // doesn't fill memory or stderr.
  if (iter_done >= 100) {
    double Dtau = D + tau;
    struct gkyl_gr_euler_newton_capfail ev = {
      .gas_gamma = gas_gamma,
      .D = D, .Sx = Sx, .Sy = Sy, .Sz = Sz, .tau = tau,
      .mom_sq = mom_sq,
      .inv_g_diag = { inv_g[0][0], inv_g[1][1], inv_g[2][2] },
      .inv_g_off  = { inv_g[0][1], inv_g[0][2], inv_g[1][2] },
      .rho = rho, .p = p, .W = W, .h = h,
      .theta = (rho > 0.0) ? p / rho : 0.0,
      .s_sq_over_Dtau_sq = (Dtau > 0.0)
        ? (Dtau * Dtau - mom_sq) / (Dtau * Dtau)
        : 0.0,
    };
    gkyl_gr_euler_record_newton_capfail(&ev);
  }

  return iter_done;
}

// Direct (non-iterative) recovery of primitives under the Mathews-Taub
// closure, following Ryu, Chattopadhyay & Choi (2006, ApJS 166, 410) §3.2.
//
// The TM EOS (h − θ)(h − 4θ) = 1 combined with the Banyuls conservation
// form yields a CUBIC in W = Γ² − 1 whose physical root has a closed-form
// trigonometric solution (eqs 23–26 of the paper). This entirely replaces
// the Picard iteration over γ_eff used by earlier implementations
// (Aitken, Anderson, Brent variants — see SESSION_NOTES post-mortem).
// We were iterating an IDEAL-EOS quartic at varying γ_eff because the
// IDEAL recovery (Eulderink-Mellema) was what we had on hand; the TM
// EOS has its OWN analytic recovery via this cubic, so no iteration is
// required at all.
//
// Mapping to our Banyuls 3+1 convention:
//   D         → paper's D
//   τ + D     → paper's E (lab-frame total including rest-mass energy)
//   γⁱʲ Sᵢ Sⱼ → paper's M² (lab-frame |momentum|² as a curved-space scalar)
//
// Inputs are *undensitized* (caller has divided by √γ). Output uses the
// same Convention-A (contravariant vⁱ, post-floor ρ, p) as the EM Newton
// path.
//
// For any admissible (D, S, τ, γⁱʲ) — i.e., D > 0, τ ≥ 0, s² > 0 — the
// cubic has three real roots, two of which give Γ < 1 (unphysical) and
// one of which gives the physical Γ ≥ 1. The trig formula picks out the
// physical root directly.
static inline void
gkyl_gr_euler_recover_primitives_tm_cubic(
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  double *rho_out, double *vx_out, double *vy_out, double *vz_out,
  double *p_out, double *W_out, double *h_out)
{
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  if (mom_sq < 0.0) mom_sq = 0.0;

  double E_lab = tau + D;
  double M2 = mom_sq;
  double D2 = D * D;
  double E2 = E_lab * E_lab;

  // Cubic coefficients (Ryu+ 2006, eq 24a-c). Denominator (E² − M²)²
  // vanishes at the luminal boundary; admissibility (s² > 0) ⇒
  // (D+τ)² > M² ⇒ E² > M² strictly, so the only risk is floating
  // noise on a cell sitting essentially on the boundary. Defensive
  // floor on the squared denominator.
  double denom = E2 - M2;
  double denom_sq = denom * denom;
  if (denom_sq < 1.0e-30) denom_sq = 1.0e-30;

  double sum_E2_M2 = E2 + M2;
  double sum_M2_D2 = M2 + D2;
  double four_sum = 4.0 * sum_E2_M2 - sum_M2_D2;
  double M2E2 = M2 * E2;

  double c1 = (sum_E2_M2 * four_sum - 14.0 * M2E2) / (2.0 * denom_sq);
  double c2 = (four_sum * four_sum - 57.0 * M2E2) / (16.0 * denom_sq);
  double c3 = -9.0 * M2E2 / (16.0 * denom_sq);

  // Solve the depressed cubic w³ + (3J)w + (−2H) = 0  where W = w − c₁/3,
  //   J = (3·c₂ − c₁²) / 9              (depressed-cubic linear coef / 3)
  //   H = (9·c₁·c₂ − 27·c₃ − 2·c₁³) / 54 (depressed-cubic constant / −2)
  // The discriminant H² + J³ controls which root family applies:
  //   ≤ 0  → three real roots (Ryu+ 2006's trig form, eq 25)
  //   > 0  → one real root + complex pair (Cardano's form)
  //
  // The paper asserts the cubic always has three real roots for admissible
  // (D, M², E_lab), but in CURVED spacetime that's empirically false:
  // Schwarzschild and Kerr inputs can land on the one-real-root branch
  // (verified in unit tests). The unified solver below handles both.
  double J = (3.0 * c2 - c1 * c1) / 9.0;
  double H = (9.0 * c1 * c2 - 27.0 * c3 - 2.0 * c1 * c1 * c1) / 54.0;
  double disc = H * H + J * J * J;  // depressed-cubic discriminant

  double W_cubic;
  if (disc <= 0.0) {
    // Three real roots — trigonometric form (Ryu+ 2006 eq 25-26). The
    // physical root corresponds to k=0 (cos(ι/3)); the other two
    // (cos((ι ± 2π)/3)) give Γ < 1.
    double neg_J = -J;
    if (neg_J < 0.0) neg_J = 0.0;  // floating-noise guard
    double sqrt_neg_J = sqrt(neg_J);
    double neg_J_pow_3half = neg_J * sqrt_neg_J;  // (−J)^{3/2}
    double cos_iota;
    if (neg_J_pow_3half > 1.0e-60) {
      cos_iota = H / neg_J_pow_3half;
      if (cos_iota >  1.0) cos_iota =  1.0;
      if (cos_iota < -1.0) cos_iota = -1.0;
    } else {
      // J ≈ 0 degenerate: cubic has triple root at W = −c₁/3.
      cos_iota = 1.0;
    }
    double iota = acos(cos_iota);
    W_cubic = 2.0 * sqrt_neg_J * cos(iota / 3.0) - c1 / 3.0;
  } else {
    // One real root — Cardano's formula:
    //   w_real = ∛(H + √disc) + ∛(H − √disc)
    //   W = w_real − c₁/3
    // cbrt() correctly handles negative arguments (real cube root).
    double sqrt_disc = sqrt(disc);
    W_cubic = cbrt(H + sqrt_disc) + cbrt(H - sqrt_disc) - c1 / 3.0;
  }
  if (W_cubic < 0.0) W_cubic = 0.0;

  double Gamma = sqrt(W_cubic + 1.0);
  if (Gamma < 1.0) Gamma = 1.0;
  double rho = D / Gamma;

  // Velocity: v² = W / (W+1) = 1 − 1/Γ², direction parallel to the
  // contravariant momentum Sⁱ = γⁱʲ Sⱼ. The magnitude of S in the
  // spatial-metric inner product is √M² = √(γⁱʲ Sᵢ Sⱼ).
  double v2 = W_cubic / (W_cubic + 1.0);
  if (v2 < 0.0) v2 = 0.0;
  if (v2 >= 1.0) v2 = 1.0 - 1.0e-12;
  double v_mag = sqrt(v2);

  double Sx_up = inv_g[0][0]*Sx + inv_g[0][1]*Sy + inv_g[0][2]*Sz;
  double Sy_up = inv_g[1][0]*Sx + inv_g[1][1]*Sy + inv_g[1][2]*Sz;
  double Sz_up = inv_g[2][0]*Sx + inv_g[2][1]*Sy + inv_g[2][2]*Sz;
  double mom_mag = sqrt(M2);
  double vx, vy, vz;
  if (mom_mag > 1.0e-30) {
    double scale = v_mag / mom_mag;
    vx = Sx_up * scale;
    vy = Sy_up * scale;
    vz = Sz_up * scale;
  } else {
    vx = vy = vz = 0.0;
  }

  // Pressure (Ryu+ 2006, eq 28c). The combination E − Mᵢvⁱ collapses to
  // the rest-frame total energy density e = ρ + ε:
  //   Mᵢvⁱ = Γ²ρhv²,   E = Γ²ρh − p
  //   ⇒ E − Mᵢvⁱ = Γ²ρh(1 − v²) − p = ρh − p = e
  // Then the TM closure p/ρ = (1/3)(e/ρ − ρ/e) gives p = (e² − ρ²)/(3e).
  double e_rest = E_lab - mom_mag * v_mag;
  double p;
  if (e_rest > GR_EULER_DENSITY_FLOOR) {
    p = (e_rest * e_rest - rho * rho) / (3.0 * e_rest);
  } else {
    p = GR_EULER_PRESSURE_FLOOR;
  }

  // TM specific enthalpy (eq 13). Use floor-clamped (ρ, p) for the
  // theta computation so a barely-recovered floor state doesn't divide
  // by zero or take sqrt of a negative.
  double rho_safe = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
  double p_safe   = (p   > GR_EULER_PRESSURE_FLOOR) ? p   : GR_EULER_PRESSURE_FLOOR;
  double theta = p_safe / rho_safe;
  double h = 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));

  *rho_out = rho;
  *vx_out  = vx;
  *vy_out  = vy;
  *vz_out  = vz;
  *p_out   = p;
  *W_out   = Gamma;
  *h_out   = h;
}

// Newton-iterative recovery of primitives under the Ryu+ 2006 RC closure
// (better Synge fit than TM: ~0.8% vs ~2% h error). Unlike the TM cubic,
// RC's primitive-recovery equation (Ryu+ eq 29) does not factor through
// a low-degree analytic form — it's a degree-8 polynomial in Γ that
// requires numerical root-finding.
//
// The equation (Ryu+ eq 29):
//   M·√(Γ²−1)·G(Γ) = H(Γ)
// where
//   G(Γ) = 24EΓ³ − 8DΓ² − 3EΓ + 2D
//   H(Γ) = 12(M²+E²)Γ⁴ − 8DEΓ³
//          + (2D² − 3M² − 12E²)Γ² + 8DEΓ − 2D²
// and (D, M², E_lab) are the same scalars used by the TM cubic.
//
// To avoid the √(Γ²−1) singularity at Γ=1 (cold flow) we Newton-iterate
// on the SQUARED form:
//   P(Γ) = M²·(Γ²−1)·G(Γ)² − H(Γ)² = 0
// This is the degree-8 polynomial. Squaring introduces a spurious sign
// branch (LHS = −RHS solutions), but starting from a good initial guess
// in the basin of the physical root avoids these.
//
// Initial guess: the TM cubic's recovered Γ. RC and TM agree on Γ to ~1%
// at any admissible state (both approximate the Synge gas; TM's 2%
// enthalpy error translates to ~1% Γ error), so Newton from TM converges
// to RC's machine-precision root in 2-4 iterations.
//
// Cap at 50 iterations as a defensive bound — paper claims "a few
// iterations" robustly; we've never hit it in tests.
static inline void
gkyl_gr_euler_recover_primitives_rc_newton(
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  double *rho_out, double *vx_out, double *vy_out, double *vz_out,
  double *p_out, double *W_out, double *h_out)
{
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  if (mom_sq < 0.0) mom_sq = 0.0;

  double E_lab = tau + D;
  double M2 = mom_sq;
  double D2 = D * D;
  double E2 = E_lab * E_lab;
  double DE = D * E_lab;

  // Initial guess: TM cubic's Γ (excellent warm start; TM and RC differ
  // by ~1% on Γ for any admissible state).
  double rho_init, vx_init, vy_init, vz_init, p_init, Gamma_init, h_init;
  gkyl_gr_euler_recover_primitives_tm_cubic(
    D, Sx, Sy, Sz, tau, inv_g,
    &rho_init, &vx_init, &vy_init, &vz_init, &p_init, &Gamma_init, &h_init);

  double Gamma = Gamma_init;
  if (Gamma < 1.0) Gamma = 1.0;

  // Newton iteration on the squared polynomial.
  const int max_iter = 50;
  const double tol_rel = 1.0e-14;
  for (int iter = 0; iter < max_iter; iter++) {
    double G2_poly = Gamma * Gamma;
    double G3_poly = G2_poly * Gamma;
    double G4_poly = G3_poly * Gamma;

    double G  = 24.0*E_lab*G3_poly - 8.0*D*G2_poly - 3.0*E_lab*Gamma + 2.0*D;
    double Gp = 72.0*E_lab*G2_poly - 16.0*D*Gamma - 3.0*E_lab;
    double H  = 12.0*(M2 + E2)*G4_poly
              - 8.0*DE*G3_poly
              + (2.0*D2 - 3.0*M2 - 12.0*E2)*G2_poly
              + 8.0*DE*Gamma
              - 2.0*D2;
    double Hp = 48.0*(M2 + E2)*G3_poly
              - 24.0*DE*G2_poly
              + 2.0*(2.0*D2 - 3.0*M2 - 12.0*E2)*Gamma
              + 8.0*DE;

    double G2m1 = G2_poly - 1.0;
    if (G2m1 < 0.0) G2m1 = 0.0;

    double P  = M2*G2m1*G*G - H*H;
    double Pp = 2.0*M2*G*(Gamma*G + G2m1*Gp) - 2.0*H*Hp;

    if (fabs(Pp) < 1.0e-30) break;

    double dGamma = -P / Pp;
    Gamma += dGamma;
    if (Gamma < 1.0) Gamma = 1.0;
    if (fabs(dGamma) < tol_rel * fabs(Gamma)) break;
  }

  // Recover primitives from converged Γ.
  if (Gamma < 1.0) Gamma = 1.0;
  double rho = D / Gamma;

  // v² = (Γ²-1)/Γ² = 1 - 1/Γ². Same identity as TM.
  double v2 = 1.0 - 1.0 / (Gamma * Gamma);
  if (v2 < 0.0) v2 = 0.0;
  if (v2 >= 1.0) v2 = 1.0 - 1.0e-12;
  double v_mag = sqrt(v2);

  // Velocity direction parallel to contravariant momentum.
  double Sx_up = inv_g[0][0]*Sx + inv_g[0][1]*Sy + inv_g[0][2]*Sz;
  double Sy_up = inv_g[1][0]*Sx + inv_g[1][1]*Sy + inv_g[1][2]*Sz;
  double Sz_up = inv_g[2][0]*Sx + inv_g[2][1]*Sy + inv_g[2][2]*Sz;
  double mom_mag = sqrt(M2);
  double vx, vy, vz;
  if (mom_mag > 1.0e-30) {
    double scale = v_mag / mom_mag;
    vx = Sx_up * scale;
    vy = Sy_up * scale;
    vz = Sz_up * scale;
  } else {
    vx = vy = vz = 0.0;
  }

  // Pressure from Ryu+ eq 33c:
  //   p = [(E − Mᵢvⁱ) − 2ρ + √((E−Mᵢvⁱ)² + 4ρ(E−Mᵢvⁱ) − 4ρ²)] / 6
  // (E − Mᵢvⁱ) = ρh − p = e (rest-frame total energy density), same
  // identity as TM. For RC, the closing relation between p, ρ, e is
  // different than TM's (e² − ρ²)/(3e) — it's the quadratic root above.
  double e_minus_p = E_lab - mom_mag * v_mag;  // = e
  double disc = e_minus_p*e_minus_p + 4.0*rho*e_minus_p - 4.0*rho*rho;
  if (disc < 0.0) disc = 0.0;
  double p = (e_minus_p - 2.0*rho + sqrt(disc)) / 6.0;
  if (p < GR_EULER_PRESSURE_FLOOR) p = GR_EULER_PRESSURE_FLOOR;

  // RC specific enthalpy at recovered (ρ, p).
  double rho_safe = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
  double p_safe   = (p   > GR_EULER_PRESSURE_FLOOR) ? p   : GR_EULER_PRESSURE_FLOOR;
  double theta = p_safe / rho_safe;
  double h = 2.0 * (6.0*theta*theta + 4.0*theta + 1.0) / (3.0*theta + 2.0);

  *rho_out = rho;
  *vx_out  = vx;
  *vy_out  = vy;
  *vz_out  = vz;
  *p_out   = p;
  *W_out   = Gamma;
  *h_out   = h;
}

// Banyuls primitive-variable recovery under Convention A. Inputs are
// undensitized conservatives; inv_g is the 3×3 inverse spatial metric
// γ^{ij} for this cell.
//
// Dispatch by EOS:
//   IDEAL         → Single Eulderink-Mellema (EM) quartic Newton pass at
//                   eos.gas_gamma. Bit-identical to the pre-refactor
//                   recover_primitives (used by the production BHL run).
//
//   MATHEWS_TAUB  → Picard iteration of EM Newton, where each pass uses a
//                   γ_eff matched to the local h_TM. For any (ρ, p) there
//                   is a unique γ_eff with h_IDEAL(γ_eff)(ρ, p) = h_TM(ρ, p):
//                     f = (h_TM - 1)/θ,  θ = p/ρ
//                     γ_eff = f/(f - 1)
//                   γ_eff is monotone-decreasing in θ over its bounded
//                   range (5/3 at θ→0, 4/3 at θ→∞).
//
// gamma_eff_cell (optional, may be NULL): pointer to a per-cell cache slot
// for γ_eff. When non-NULL and the EOS is non-IDEAL, the Picard iteration
// is initialised from *gamma_eff_cell (the previous step's converged value)
// and the final converged γ_eff is written back. Warm-starting drops
// Picard iteration counts from ~10 to ~1–2 for cells in steady state. For
// IDEAL the pointer is ignored (no Picard loop). NULL preserves the cold-
// flow default initial guess of γ=5/3 (used by unit tests and any call
// site without per-cell context).
//
// Both paths recover the same convention (contravariant v^i, post-floor
// ρ, p) and populate the same out struct.
static inline void
gkyl_gr_euler_recover_primitives(
  struct gkyl_gr_euler_eos eos,
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  double *gamma_eff_cell,
  struct gkyl_gr_euler_prim *out)
{
  out->admissible =
    gkyl_gr_euler_check_admissibility(D, Sx, Sy, Sz, tau, inv_g) == GR_EULER_ADM_OK;

  // Lorentz scalar with covariant momentum: |S|² = γ^{ij} S_i S_j.
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  if (mom_sq < 0.0) mom_sq = 0.0;  // defensive against floating-point noise

  double rho, vx, vy, vz, p, W, h;

  if (eos.type == GR_EULER_EOS_IDEAL) {
    // Single EM Newton pass — bit-identical to pre-refactor.
    (void)gkyl_gr_euler_em_newton_at_gamma(
      eos.gas_gamma, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);
    gkyl_gr_euler_record_picard_iters(1);  // IDEAL is by definition 1 pass

  }
  else if (eos.type == GR_EULER_EOS_RYU_CHATTOPADHYAY) {
    // RC dispatch:
    //   θ_tm < 1e-6 or p_tm floored → EM Newton at γ=5/3 (cold-flow
    //                                  / TM-failure fallback)
    //   M = 0                       → static-fluid short-circuit, accept
    //                                  TM directly
    //   else                        → RC Newton, paper-grounded physicality
    //                                  check, accept if physical else fall
    //                                  back to TM
    //
    // The physicality check is built directly from Ryu+ 2006 §3.3:
    //   (a) W_rc > 1               (strict sub-luminal lower bound)
    //   (b) W_rc < Γ_u             (strict sub-luminal upper bound, eq 30)
    //   (c) sign(G(W_rc)) == sign(H(W_rc))   (eq 29 sign-branch check)
    //
    // All three inequalities are strict by design: this branch only runs
    // when M > 0, so Γ_u > 1 and Γ_l = 1 < Γ_phys < Γ_u with strict
    // separation. The M = 0 boundary (where Γ_l = Γ_u = 1 collapses the
    // range to a point and RC's squared polynomial has a double root)
    // is handled by the short-circuit above — RC has nothing to add
    // when v = 0.
    //
    // Together (a)-(c) uniquely identify the paper's physical root.
    // Squaring eq 29 to get the Newton polynomial doubles the algebraic
    // root count by introducing the LHS = −RHS sign-flipped roots;
    // check (c) filters those out. Check (b) rules out the paper's
    // "positive root larger than Γ_u" (unphysical supraluminal root);
    // see Ryu+ eq 30 + §3.3 second-to-last paragraph.
    //
    // Cancellation-free form of (b): W_rc < Γ_u  ⟺  W_rc² · s² < E²
    // where s² = (D+τ)² − M². Avoids dividing by (1 − M²/E²) which
    // can be very small for high-W cells.
    //
    // This replaces an earlier "trust RC if it agrees with TM within
    // 10%" heuristic plus a "conditioning escape" branch — both were
    // proxies for paper-grounded physicality before we wrote it out
    // directly. The new check is strictly more rigorous and removes
    // the need to compare against TM at all (TM is now purely the
    // Newton warm-start and the fallback when RC fails physicality).
    gkyl_gr_euler_recover_primitives_tm_cubic(
      D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);

    double rho_for_theta = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
    double theta_tm = p / rho_for_theta;
    // TM's pressure floor activated ⇒ TM failed to find a positive-p
    // solution. This is a more honest signal than θ_tm because TM's
    // precision-loss feedback loop near the admissibility boundary
    // inflates W ⇒ deflates ρ = D/W ⇒ artificially inflates
    // θ = p_floor/ρ above the cold-flow threshold even though the cell
    // is effectively pressureless. Route these directly to EM γ=5/3
    // alongside genuine cold-flow cells.
    bool tm_p_floored = (p <= GR_EULER_PRESSURE_FLOOR);

    if (theta_tm < 1.0e-6 || tm_p_floored) {
      (void)gkyl_gr_euler_em_newton_at_gamma(
        5.0/3.0, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
        &rho, &vx, &vy, &vz, &p, &W, &h);
    } else if (mom_sq == 0.0) {
      // Static-fluid short-circuit. v = 0 ⇒ W = 1 exactly and (ρ, p) are
      // determined from (D, τ) by TM's cubic (closed-form at M = 0).
      // RC's squared polynomial degenerates here: H(W) factors as
      // (W−1)(W+1)·2·(6E²W² − 4DE·W + D²) with the quadratic having
      // discriminant −8D²E² < 0, so the only real positive root of H²=0
      // is the double root W = 1. Newton on a double root converges
      // linearly and stops at W_rc = 1 ± √ε_machine, which makes the
      // strict-inequality physicality checks below trip on the boundary.
      // TM has already produced the exact answer (vx = vy = vz = 0,
      // W = 1) at this point — just take it.
    } else {
      double rho_rc, vx_rc, vy_rc, vz_rc, p_rc, W_rc, h_rc;
      gkyl_gr_euler_recover_primitives_rc_newton(
        D, Sx, Sy, Sz, tau, inv_g,
        &rho_rc, &vx_rc, &vy_rc, &vz_rc, &p_rc, &W_rc, &h_rc);

      // (a) Strict sub-luminal lower bound + NaN guard.
      bool rc_lower_bound = (W_rc > 1.0);

      // (b) Strict sub-luminal upper bound — Ryu+ eq 30 in
      //     cancellation-free form: W_rc² · (E² − M²) < E². The M = 0
      //     boundary case is handled by the static-fluid short-circuit
      //     above, so reaching this branch implies M > 0 and the strict
      //     inequality is the correct semantic.
      double E_lab_sq = (D + tau) * (D + tau);
      double s_sq_unfloored = E_lab_sq - mom_sq;
      bool rc_sub_luminal = (s_sq_unfloored > 0.0)
        && (W_rc * W_rc * s_sq_unfloored < E_lab_sq);

      // (c) Eq 29 strict sign-branch check. The un-squared eq 29 is
      //     M·√(Γ²−1)·G(Γ) = H(Γ). M > 0 here (M = 0 handled above),
      //     and Γ > 1 (lower bound), so sign(LHS) = sign(G). Physical
      //     root has sign(G) == sign(H) strictly.
      double E_lab = D + tau;
      double E2 = E_lab_sq;
      double M2 = mom_sq;
      double D2 = D * D;
      double DE = D * E_lab;
      double Wr2 = W_rc * W_rc;
      double Wr3 = Wr2 * W_rc;
      double Wr4 = Wr3 * W_rc;
      double G_at_rc = 24.0*E_lab*Wr3 - 8.0*D*Wr2 - 3.0*E_lab*W_rc + 2.0*D;
      double H_at_rc = 12.0*(M2 + E2)*Wr4 - 8.0*DE*Wr3
                     + (2.0*D2 - 3.0*M2 - 12.0*E2)*Wr2 + 8.0*DE*W_rc - 2.0*D2;
      bool rc_correct_branch = (G_at_rc * H_at_rc > 0.0);

      bool rc_ok = rc_lower_bound && rc_sub_luminal && rc_correct_branch;

      // Diagnostic: track acceptance / rejection magnitude. dW_rel
      // computed only for the counter; it does NOT participate in the
      // accept/reject decision anymore. `conditioning_escape` is kept
      // in the API signature as a vestigial bool (always false now);
      // simpler to leave the diagnostic API stable.
      double dW_rel = fabs(W_rc - W) / ((W > 0.0) ? W : 1.0);
      gkyl_gr_euler_record_rc_compare(rc_ok, false, dW_rel);

      if (rc_ok) {
        rho = rho_rc; vx = vx_rc; vy = vy_rc; vz = vz_rc;
        p   = p_rc;   W  = W_rc;  h  = h_rc;
      } else {
        // Capture diagnostic snapshot for the post-run dump.
        double s_sq_over_E2 = (E_lab_sq > 1.0e-30)
          ? s_sq_unfloored / E_lab_sq : 0.0;
        struct gkyl_gr_euler_rc_reject ev = {
          .dW_rel       = dW_rel,
          .W_tm         = W,
          .W_rc         = W_rc,
          .D            = D,
          .tau          = tau,
          .mom_sq       = mom_sq,
          .s_sq_over_E2 = s_sq_over_E2,
          .theta_tm     = theta_tm,
          .p_tm         = p,
          .p_rc         = p_rc,
        };
        gkyl_gr_euler_record_rc_reject(&ev);
      }
    }

    gkyl_gr_euler_record_picard_iters(1);

    if (gamma_eff_cell != NULL) {
      *gamma_eff_cell = gkyl_gr_euler_tm_match_gamma_eff(rho, p);
    }

  }
  else {
    // Direct (non-iterative) Mathews-Taub recovery via the Ryu+ 2006
    // cubic (see gkyl_gr_euler_recover_primitives_tm_cubic above) +
    // cold-flow fallback to EM Newton at γ=5/3 when θ = p/ρ is below
    // the EOS-equivalence threshold.
    //
    // Replaces the Picard machinery (Aitken / Anderson / Brent variants
    // from earlier sessions) — TM admits a closed-form analytic
    // recovery, derived independently from the IDEAL Eulderink-Mellema
    // quartic; we route to it for the bulk of cells.
    gkyl_gr_euler_recover_primitives_tm_cubic(
      D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);

    // COLD-FLOW FALLBACK. At low θ = p/ρ:
    //
    //   h_TM(θ)        = 1 + (5/2)θ + (9/8)θ² + O(θ⁴)
    //   h_IDEAL(γ=5/3) = 1 + (5/2)θ                 (exactly)
    //
    // so |h_TM − h_IDEAL|/h ~ (9/8)θ² — quadratic falloff. The
    // threshold trades off cubic-vs-EM-fallback rate against EOS
    // accuracy in the cold tail:
    //
    //   θ < 1e-8 → h agreement ~1e-16 (machine precision, bit-identical)
    //   θ < 1e-6 → h agreement ~1e-13
    //   θ < 1e-4 → h agreement ~1e-9   (still well below sim tolerance)
    //   θ < 1e-3 → h agreement ~1e-7
    //
    // Using 1e-4: catches the marginal-cubic-precision regime around
    // the compact object where (e − ρ)/ρ is small enough that the
    // cubic's p formula suffers cancellation, but generous enough to
    // include the "barely above the floor" cells whose cubic primitives
    // were inconsistent at the 1e-6 threshold (showed up in BHL as
    // residual ~10× excess repair counts vs Brent baseline).
    //
    // (Synge gas has same scaling with coefficient 15/8 — slightly
    // slower convergence but same regime.)
    //
    // At the same time, the cubic LOSES SELF-CONSISTENCY in this
    // regime because its pressure formula p = (e² − ρ²)/(3e) suffers
    // cancellation when e ≈ ρ (cold flow ⇔ h ≈ 1 ⇔ e ≈ ρ exactly).
    // The raw p comes out near or below the pressure floor; the outer
    // floor clamp then leaves the remaining primitives (W, ρ, h, vⁱ)
    // inconsistent with the floored p. Downstream, that inconsistency
    // manifests as repeated cascade-repair fires near the compact
    // object in BHL.
    //
    // EM Newton at γ=5/3 is the gold-standard cold-flow recovery:
    // robust convergence for p << 1, and gives a self-consistent
    // regularized primitive set that downstream tolerates cleanly.
    // This explicit fallback mirrors the Picard scheme's implicit
    // cold-flow shortcut (θ < 1e-12 → γ_eff = 5/3).
    double rho_for_theta = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
    double theta_check = p / rho_for_theta;
    if (theta_check < 1.0e-4) {
      (void)gkyl_gr_euler_em_newton_at_gamma(
        5.0/3.0, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
        &rho, &vx, &vy, &vz, &p, &W, &h);
    }

    // Histogram: 1 "pass" regardless of which branch fired — both the
    // cubic and the EM fallback are a single recovery call.
    gkyl_gr_euler_record_picard_iters(1);

    // Diagnostic γ_eff writeback (the IDEAL adiabatic index that would
    // give the same h as TM at the recovered (ρ, p)). Useful for
    // external monitors that want to see the per-cell effective regime;
    // not consumed as a warm-start anywhere.
    if (gamma_eff_cell != NULL) {
      *gamma_eff_cell = gkyl_gr_euler_tm_match_gamma_eff(rho, p);
    }
  }

  bool rho_floored = (rho < GR_EULER_DENSITY_FLOOR);
  bool p_floored   = (p   < GR_EULER_PRESSURE_FLOOR);
  // Pass raw values + conservatives to the diagnostic before clamping.
  gkyl_gr_euler_record_floor_hit(
    rho_floored, p_floored, D, tau, mom_sq, rho, p);
  if (rho_floored) rho = GR_EULER_DENSITY_FLOOR;
  if (p_floored)   p   = GR_EULER_PRESSURE_FLOOR;

  out->rho = rho;
  out->v[0] = vx;
  out->v[1] = vy;
  out->v[2] = vz;
  out->p = p;
  out->W = W;
  out->h = h;
}
