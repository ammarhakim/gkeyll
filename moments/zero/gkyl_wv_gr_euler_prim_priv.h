#pragma once

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

// Shared primitive-variable recovery + admissibility check + iterative
// state-repair helpers for the modular tetrad GR Euler equation
// (wv_gr_euler_tetrad.c) and its source-term integrator
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
// Two EOS types are supported:
//
//   IDEAL             — p = (γ-1)·ρ·ε  ⇒  h = 1 + γ/(γ-1)·p/ρ. Caller
//                       specifies gas_gamma; use_rcc is ignored.
//
//   APPROXIMATE_SYNGE — Single-component perfect relativistic gas
//                       approximation. Both branches share the same
//                       non-relativistic (Γ→5/3) and ultra-relativistic
//                       (Γ→4/3, c_s²→1/3) asymptotic limits and use the
//                       Ryu+ 2006 §3.2 closed-form TM cubic for the
//                       primary recovery. The `use_rcc` bool selects
//                       between the two enthalpy closures:
//                         use_rcc = false → Taub-Mathews (Mignone+ 2005)
//                                           h(θ) = (5θ + √(9θ²+4))/2
//                                           ~2% h error vs Synge.
//                         use_rcc = true  → Ryu-Chattopadhyay (Ryu+ 2006)
//                                           h(θ) = 2(6θ²+4θ+1)/(3θ+2)
//                                           ~0.8% h error vs Synge.
//                       For both, the recovery dispatch is: TM cubic →
//                       cold-flow EM γ=5/3 fallback when θ<1e-6 or TM
//                       floors. If use_rcc, additionally run RC Newton
//                       on the degree-8 polynomial (Ryu+ eq 29) and
//                       accept its result when it passes the strict
//                       paper-grounded physicality check; otherwise fall
//                       back to TM.
// ---------------------------------------------------------------------------

enum gkyl_gr_euler_eos_type {
  GR_EULER_EOS_IDEAL = 0,
  GR_EULER_EOS_APPROXIMATE_SYNGE = 1,
};

struct gkyl_gr_euler_eos {
  enum gkyl_gr_euler_eos_type type;
  double gas_gamma;  // adiabatic index (consulted only when type == IDEAL)
  bool use_rcc;      // RCC vs TM enthalpy closure within APPROXIMATE_SYNGE;
                     // ignored for IDEAL
};

// Specific enthalpy h(ρ, p). θ = p/ρ.
//   IDEAL:                       h = 1 + γ/(γ-1)·θ
//   APPROXIMATE_SYNGE (use_rcc=false, TM):
//                                h = (5θ + sqrt(9θ²+4))/2
//   APPROXIMATE_SYNGE (use_rcc=true, RCC):
//                                h = 2(6θ²+4θ+1)/(3θ+2)             [Ryu+ eq 15]
static inline double
gkyl_gr_euler_eos_enthalpy(struct gkyl_gr_euler_eos eos,
  double rho, double p)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return 1.0 + (p / rho) * (eos.gas_gamma / (eos.gas_gamma - 1.0));
  }
  double theta = p / rho;
  if (eos.use_rcc) {
    return 2.0 * (6.0*theta*theta + 4.0*theta + 1.0) / (3.0*theta + 2.0);
  }
  return 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));
}

// Pressure recovery p(ρ, h). EOS-specific inverse of h(ρ, p), used
// in cell-interface state reconstruction.
//   IDEAL:                       θ = (h-1)·(γ-1)/γ
//   APPROXIMATE_SYNGE (TM):      4θ² − 5hθ + (h² − 1) = 0
//                                ⇒ θ = (5h − sqrt(9h² + 16))/8
//   APPROXIMATE_SYNGE (RCC):     12θ² + (8-3h)θ + (2-2h) = 0
//                                ⇒ θ = (3h − 8 + sqrt(9h² + 48h − 32))/24
//                                                                    [Ryu+ eq 60]
static inline double
gkyl_gr_euler_eos_pressure_from_rho_h(struct gkyl_gr_euler_eos eos,
  double rho, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return rho * (h - 1.0) * (eos.gas_gamma - 1.0) / eos.gas_gamma;
  }
  if (eos.use_rcc) {
    double theta = (3.0*h - 8.0 + sqrt(9.0*h*h + 48.0*h - 32.0)) / 24.0;
    return rho * theta;
  }
  double theta = 0.125 * (5.0 * h - sqrt(9.0 * h * h + 16.0));
  return rho * theta;
}

// Sound speed squared c_s²(ρ, p, h). Cold-flow limit (θ→0) is 5p/(3ρh)
// (γ=5/3 ideal); ultra-rel limit (θ→∞) is 1/3 (radiation fluid).
//   IDEAL:                       c_s² = γp/(ρh)
//   APPROXIMATE_SYNGE (TM):      c_s² = θ(5h − 8θ) / (3·h·(h − θ))
//   APPROXIMATE_SYNGE (RCC):     c_s² = θ(3θ+2)(18θ²+24θ+5)
//                                       / (3·(6θ²+4θ+1)·(9θ²+12θ+2))  [Ryu+ eq 16]
static inline double
gkyl_gr_euler_eos_cs2(struct gkyl_gr_euler_eos eos,
  double rho, double p, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return (eos.gas_gamma * p) / (rho * h);
  }
  double theta = p / rho;
  if (eos.use_rcc) {
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
// Recovery-iteration instrumentation: per-callsite status structs.
//
// Modeled on gkyl_vlasov_lte_correct_status, the prim_status struct holds
// the running iteration counts, path-dispatch tallies, and floor-hit
// histograms for one logical callsite (e.g. the wave-prop hooks of a
// species, or the source-step recovery of that species). The app
// allocates one prim_status per (species, callsite-class) pair and
// passes a pointer to gkyl_gr_euler_recover_primitives so the recovery
// accumulates into it.
//
// Pass NULL for the status pointer in unit-test / probe callsites that
// don't care about instrumentation.
// ---------------------------------------------------------------------------

// Bucket boundaries for the Newton iteration histograms (upper edge,
// inclusive). Bin k counts calls with iter count ≤ NEWTON_BIN_EDGES[k]
// and > NEWTON_BIN_EDGES[k-1]. The final bin is "≥ last edge" overflow.
// Shared across the EM-Newton (cold-flow + IDEAL) and RC-Newton (RCC
// refinement) iteration trackers since both have similar convergence
// characteristics.
#define GR_EULER_NEWTON_HIST_NBINS 7
// 0–4, 5–8, 9–16, 17–32, 33–64, 65–99, 100+
extern const int gkyl_gr_euler_newton_bin_edges[GR_EULER_NEWTON_HIST_NBINS];

// Decade bins for the D-magnitude histogram of floored cells (answers
// "are floored cells near-vacuum or normal-density?")
#define GR_EULER_FLOOR_D_BINS 11
// <1e-10, 1e-10..1e-8, 1e-8..1e-6, 1e-6..1e-4, 1e-4..1e-2, 1e-2..1e-1,
// 1e-1..1, 1..10, 10..100, 100..1000, ≥1000

// Decade bins for the (E²−M²)/E² histogram of floored cells (Newton
// coefficient-conditioning measure ~1/W² for floor-typical cells —
// answers "are floored cells near-luminal or moderate W?")
#define GR_EULER_FLOOR_S2_BINS 8
// <1e-12, 1e-12..1e-10, 1e-10..1e-8, 1e-8..1e-6, 1e-6..1e-4, 1e-4..1e-2,
// 1e-2..1, ≥1

// Per-callsite recovery status. Accumulated by gkyl_gr_euler_recover_primitives.
//
// Newton iteration trackers: separate (calls, total iters, max, histogram)
// blocks for the EM-Newton inner solver (used by the IDEAL path and the
// cold-flow EM γ=5/3 fallback) and the RC-Newton outer solver (used only
// by APPROXIMATE_SYNGE with use_rcc = true). Lets the postmortem
// distinguish where wall-clock cost lives.
//
// Path counters: how often each dispatch branch fired. For IDEAL only
// path_ideal is ever non-zero. For APPROXIMATE_SYNGE the cells are
// distributed across (path_tm_cubic, path_cold_flow, path_rcc_accepted,
// path_rcc_rejected). path_rcc_* only non-zero when use_rcc = true.
//
// Floor stats: counts of recovery calls whose raw (pre-floor) ρ or p
// fell below the floor and was clamped. The decade-binned histograms
// characterize WHERE floored cells live in (D, conditioning) space.
struct gkyl_gr_euler_prim_status {
  // EM-Newton (IDEAL + cold-flow EM γ=5/3 fallback).
  uint64_t em_newton_calls;
  uint64_t em_newton_total_iters;
  int      em_newton_max_iters;
  uint64_t em_newton_iter_hist[GR_EULER_NEWTON_HIST_NBINS];
  uint64_t em_newton_cap_hits;   // calls that hit the 100-iter cap

  // RC-Newton (APPROXIMATE_SYNGE + use_rcc only).
  uint64_t rc_newton_calls;
  uint64_t rc_newton_total_iters;
  int      rc_newton_max_iters;
  uint64_t rc_newton_iter_hist[GR_EULER_NEWTON_HIST_NBINS];

  // Path counters — which dispatch branch ran per recovery call.
  uint64_t path_ideal;            // IDEAL: single EM Newton pass
  uint64_t path_tm_cubic;         // APPROXIMATE_SYNGE: TM cubic accepted
                                  //   (no cold-flow fallback, no RC refine)
  uint64_t path_cold_flow;        // APPROXIMATE_SYNGE: θ_tm<1e-6 or
                                  //   tm_p_floored → EM γ=5/3 fallback
  uint64_t path_rcc_accepted;     // APPROXIMATE_SYNGE + use_rcc: RC
                                  //   Newton accepted by (a)/(b)/(c) check
  uint64_t path_rcc_rejected;     // APPROXIMATE_SYNGE + use_rcc: RC
                                  //   rejected, fell back to TM cubic

  // Floor stats.
  uint64_t floor_calls;           // total recovery calls (one per invoke)
  uint64_t rho_floor_hits;        // ρ floored (raw < FLOOR before clamp)
  uint64_t p_floor_hits;          // p floored
  uint64_t both_floor_hits;       // BOTH ρ and p floored (subset)
  uint64_t floor_D_hist[GR_EULER_FLOOR_D_BINS];
  uint64_t floor_s2_hist[GR_EULER_FLOOR_S2_BINS];

  // RCC sanity diagnostics (only meaningful with use_rcc = true).
  // sum_dW_acc / max_dW_acc accumulate on path_rcc_accepted; max_dW_rej
  // on path_rcc_rejected. Lets the postmortem confirm that accepted
  // |ΔΓ|/Γ_tm is small (refinement is fine-grained) and rejected
  // |ΔΓ|/Γ_tm is large (rejections are real disagreements).
  double sum_dW_acc;
  double max_dW_acc;
  double max_dW_rej;

  // HLLC star-state diagnostics. Only populated by sr_hllc_minkowski (no-op
  // for the IDEAL+HLL/Roe/Lax paths — fields stay zero). Aggregated counter
  // + per-reason histogram for production postmortem; last_* fields capture
  // the most recent HLLC call for fine-grained inspection by single-state
  // tests (overwritten each call). Reason key:
  //   0 = no fallback (HLLC star-state used)
  //   1 = lam_diff ≈ 0 (degenerate Davis bracket — λ_L ≈ λ_R)
  //   2 = λ* not finite
  //   3 = |λ_L − λ*| < tol (would blow up 1/(λ_L − λ*) in U_L*)
  //   4 = |λ_R − λ*| < tol (would blow up 1/(λ_R − λ*) in U_R*)
  struct {
    uint64_t fallback_calls;             // HLLC calls that fell back to HLL
    uint64_t fallback_reason_hist[5];    // bin by reason index
    int      last_did_fallback;          // 0 or 1
    int      last_fallback_reason;       // 0..4
    double   last_lambda_L;
    double   last_lambda_R;
    double   last_lambda_star;
  } hllc;
};

// Per-callsite repair-cascade status. Tracks how often the repair
// cascade fixed each admissibility constraint, plus the source-step
// tau/s² limiter firing counts (the latter are only meaningful for the
// source-step callsite — wave-prop's status will see zeros there).
struct gkyl_gr_euler_repair_status {
  uint64_t bad_D_fixes;       // D≤0 floored to GR_EULER_DENSITY_FLOOR
  uint64_t bad_tau_fixes;     // τ<0 floored to GR_EULER_TAU_REPAIR_FLOOR
  uint64_t bad_s2_fixes;      // s²≤0 rescaled momentum to interior of A_γ
  uint64_t tau_limiter_fires; // source-step τ-positivity limiter (α<1)
  uint64_t s2_limiter_fires;  // source-step s²-positivity limiter
};

// Decade-bin selector for the floor D-magnitude histogram.
static inline int
gkyl_gr_euler_status_floor_D_bin(double D)
{
  if (D < 1.0e-10) return 0;
  if (D < 1.0e-8)  return 1;
  if (D < 1.0e-6)  return 2;
  if (D < 1.0e-4)  return 3;
  if (D < 1.0e-2)  return 4;
  if (D < 1.0e-1)  return 5;
  if (D < 1.0)     return 6;
  if (D < 10.0)    return 7;
  if (D < 100.0)   return 8;
  if (D < 1000.0)  return 9;
  return 10;
}

// Decade-bin selector for the floor (E²−M²)/E² histogram.
static inline int
gkyl_gr_euler_status_floor_s2_bin(double s_sq_over_E2)
{
  if (s_sq_over_E2 < 1.0e-12) return 0;
  if (s_sq_over_E2 < 1.0e-10) return 1;
  if (s_sq_over_E2 < 1.0e-8)  return 2;
  if (s_sq_over_E2 < 1.0e-6)  return 3;
  if (s_sq_over_E2 < 1.0e-4)  return 4;
  if (s_sq_over_E2 < 1.0e-2)  return 5;
  if (s_sq_over_E2 < 1.0)     return 6;
  return 7;
}

// Histogram bucket lookup (≤ semantic on the edge).
static inline int
gkyl_gr_euler_status_hist_bucket(int n, const int *edges, int nbins)
{
  for (int b = 0; b < nbins; b++) {
    if (n <= edges[b]) return b;
  }
  return nbins - 1;
}

// One Eulderink-Mellema quartic Newton pass at fixed γ. Recovers
// (ρ, v^i, p, W, h) from undensitized (D, S_i, τ) and inv_g.
//
// The "guess" iterate is 1/(h·C), C = D/√s², which gives a 4th-order
// polynomial residual that converges quadratically from guess=1.0 for
// any admissible input. Outputs are pre-floor so callers can apply
// floors / dispatch after refinement.
//
// Returns the iteration count: < 100 = converged; 100 = hit cap. The
// caller (gkyl_gr_euler_recover_primitives) accumulates this into the
// EM-Newton histogram block of its prim_status if non-NULL.
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
// Initial guess (Gamma_init): caller-supplied TM cubic Γ. RC and TM
// agree on Γ to ~1% at any admissible state (both approximate the Synge
// gas; TM's 2% enthalpy error translates to ~1% Γ error), so Newton
// from TM converges to RC's machine-precision root in 2-4 iterations.
// The caller in gkyl_gr_euler_recover_primitives already ran the TM
// cubic for dispatch decision; passing its Γ here avoids a redundant
// second cubic solve.
//
// Returns the iteration count (1..50). Caller accumulates this into
// the rc_newton histogram block of its prim_status if non-NULL. Cap at
// 50 iterations as a defensive bound — paper claims "a few iterations"
// robustly; we've never hit it in tests.
static inline int
gkyl_gr_euler_recover_primitives_rc_newton(
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  double Gamma_init,
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

  double Gamma = Gamma_init;
  if (Gamma < 1.0) Gamma = 1.0;

  // Newton iteration on the squared polynomial.
  const int max_iter = 50;
  const double tol_rel = 1.0e-14;
  int iter_done = max_iter;
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

    if (fabs(Pp) < 1.0e-30) {
      iter_done = iter + 1;
      break;
    }

    double dGamma = -P / Pp;
    Gamma += dGamma;
    if (Gamma < 1.0) Gamma = 1.0;
    if (fabs(dGamma) < tol_rel * fabs(Gamma)) {
      iter_done = iter + 1;
      break;
    }
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

  return iter_done;
}

// Banyuls primitive-variable recovery under Convention A. Inputs are
// undensitized conservatives; inv_g is the 3×3 inverse spatial metric
// γ^{ij} for this cell.
//
// Dispatch by EOS:
//   IDEAL              → Single Eulderink-Mellema (EM) quartic Newton
//                        pass at eos.gas_gamma.
//   APPROXIMATE_SYNGE  → (a) Closed-form TM cubic (Ryu+ 2006 §3.2).
//                        (b) Cold-flow / TM-failure fallback — EM Newton
//                            at γ=5/3 when θ_tm < 1e-6 or TM's pressure
//                            floor triggers. Cells in this regime have
//                            TM precision-corrupted primitives that
//                            downstream cannot consume cleanly; EM γ=5/3
//                            gives a self-consistent regularized state
//                            with h_TM ≈ h_IDEAL agreement ≲ 1e-13.
//                        (c) When use_rcc and mom_sq > 0, additionally
//                            run RC Newton on Ryu+ 2006 eq 29 squared
//                            polynomial (warm-started from TM cubic) and
//                            accept its result iff it passes the strict
//                            Ryu+ 2006 §3.3 physicality check (see
//                            below). Otherwise fall back to TM cubic.
//                        (d) When !use_rcc, the TM cubic from step (a)
//                            is the final answer.
//
// stat (optional, may be NULL): per-callsite prim_status to accumulate
// iteration counts, path counters, and floor stats into. NULL skips all
// instrumentation. Production callsites in the app pass per-species
// status objects; unit tests typically pass NULL.
//
// All paths recover the same convention (contravariant v^i, post-floor
// ρ, p) and populate the same out struct.
static inline void
gkyl_gr_euler_recover_primitives(
  struct gkyl_gr_euler_eos eos,
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  struct gkyl_gr_euler_prim_status *stat,
  struct gkyl_gr_euler_prim *out)
{
  out->admissible =
    gkyl_gr_euler_check_admissibility(D, Sx, Sy, Sz, tau, inv_g) == GR_EULER_ADM_OK;

  // Lorentz scalar with covariant momentum: |S|² = γ^{ij} S_i S_j.
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
  if (mom_sq < 0.0) mom_sq = 0.0;  // defensive against floating-point noise

  double rho, vx, vy, vz, p, W, h;

  if (eos.type == GR_EULER_EOS_IDEAL) {
    // Single EM Newton pass.
    int em_iters = gkyl_gr_euler_em_newton_at_gamma(
      eos.gas_gamma, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);
    if (stat) {
      stat->em_newton_calls       += 1;
      stat->em_newton_total_iters += (uint64_t)em_iters;
      if (em_iters > stat->em_newton_max_iters)
        stat->em_newton_max_iters = em_iters;
      stat->em_newton_iter_hist[gkyl_gr_euler_status_hist_bucket(
        em_iters, gkyl_gr_euler_newton_bin_edges,
        GR_EULER_NEWTON_HIST_NBINS)] += 1;
      if (em_iters >= 100) stat->em_newton_cap_hits += 1;
      stat->path_ideal += 1;
    }
  }
  else {
    // APPROXIMATE_SYNGE. Step (a): closed-form TM cubic.
    gkyl_gr_euler_recover_primitives_tm_cubic(
      D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);

    // Step (b): cold-flow / TM-failure fallback. tm_p_floored is more
    // honest than θ_tm alone — TM's precision-loss feedback loop near
    // the admissibility boundary inflates W ⇒ deflates ρ = D/W ⇒
    // artificially inflates θ = p_floor/ρ above the cold-flow
    // threshold even though the cell is effectively pressureless.
    // Route these directly to EM γ=5/3 alongside genuine cold-flow
    // cells.
    double rho_for_theta = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
    double theta_tm = p / rho_for_theta;
    bool tm_p_floored = (p <= GR_EULER_PRESSURE_FLOOR);
    if (theta_tm < 1.0e-6 || tm_p_floored) {
      int em_iters = gkyl_gr_euler_em_newton_at_gamma(
        5.0/3.0, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
        &rho, &vx, &vy, &vz, &p, &W, &h);
      if (stat) {
        stat->em_newton_calls       += 1;
        stat->em_newton_total_iters += (uint64_t)em_iters;
        if (em_iters > stat->em_newton_max_iters)
          stat->em_newton_max_iters = em_iters;
        stat->em_newton_iter_hist[gkyl_gr_euler_status_hist_bucket(
          em_iters, gkyl_gr_euler_newton_bin_edges,
          GR_EULER_NEWTON_HIST_NBINS)] += 1;
        if (em_iters >= 100) stat->em_newton_cap_hits += 1;
        stat->path_cold_flow += 1;
      }
    }
    // Step (c): RC Newton refinement. Only runs when use_rcc and the
    // momentum is non-zero (the M = 0 boundary case is handled by TM
    // exactly; RC's squared polynomial degenerates to a double root at
    // W = 1 there).
    //
    // The physicality check is built directly from Ryu+ 2006 §3.3:
    //   (a) W_rc > 1               (strict sub-luminal lower bound)
    //   (b) W_rc < Γ_u             (strict sub-luminal upper bound,
    //                               Ryu+ eq 30, cancellation-free form
    //                               W_rc² · (E² − M²) < E²)
    //   (c) sign(G(W_rc)) == sign(H(W_rc))   (eq 29 sign-branch check)
    // Together (a)-(c) uniquely identify the paper's physical root.
    // Squaring eq 29 to get the Newton polynomial doubles the algebraic
    // root count by introducing the LHS = −RHS sign-flipped roots;
    // check (c) filters those out. Check (b) rules out the paper's
    // positive root larger than Γ_u (supraluminal root).
    //
    // Step (d): when !use_rcc, the TM cubic is the final answer; this
    // branch is skipped.
    else if (eos.use_rcc && mom_sq > 0.0) {
      // Warm-start the RC Newton from the TM cubic's Γ we already have.
      double rho_rc, vx_rc, vy_rc, vz_rc, p_rc, W_rc, h_rc;
      int rc_iters = gkyl_gr_euler_recover_primitives_rc_newton(
        D, Sx, Sy, Sz, tau, inv_g, W,
        &rho_rc, &vx_rc, &vy_rc, &vz_rc, &p_rc, &W_rc, &h_rc);
      if (stat) {
        stat->rc_newton_calls       += 1;
        stat->rc_newton_total_iters += (uint64_t)rc_iters;
        if (rc_iters > stat->rc_newton_max_iters)
          stat->rc_newton_max_iters = rc_iters;
        stat->rc_newton_iter_hist[gkyl_gr_euler_status_hist_bucket(
          rc_iters, gkyl_gr_euler_newton_bin_edges,
          GR_EULER_NEWTON_HIST_NBINS)] += 1;
      }

      // (a) Strict sub-luminal lower bound + NaN guard.
      bool rc_lower_bound = (W_rc > 1.0);

      // (b) Strict sub-luminal upper bound. M > 0 here (M = 0 short-
      //     circuited above), so the strict inequality is the correct
      //     semantic.
      double E_lab_sq = (D + tau) * (D + tau);
      double s_sq_unfloored = E_lab_sq - mom_sq;
      bool rc_sub_luminal = (s_sq_unfloored > 0.0)
        && (W_rc * W_rc * s_sq_unfloored < E_lab_sq);

      // (c) Eq 29 strict sign-branch check.
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

      // |ΔΓ|/Γ_tm — diagnostic only, no role in accept/reject decision.
      double dW_rel = fabs(W_rc - W) / ((W > 0.0) ? W : 1.0);

      if (rc_ok) {
        rho = rho_rc; vx = vx_rc; vy = vy_rc; vz = vz_rc;
        p   = p_rc;   W  = W_rc;  h  = h_rc;
        if (stat) {
          stat->path_rcc_accepted += 1;
          stat->sum_dW_acc += dW_rel;
          if (dW_rel > stat->max_dW_acc) stat->max_dW_acc = dW_rel;
        }
      } else {
        if (stat) {
          stat->path_rcc_rejected += 1;
          if (dW_rel > stat->max_dW_rej) stat->max_dW_rej = dW_rel;
        }
      }
    }
    else if (stat) {
      // TM cubic accepted without RC refinement (use_rcc=false or M=0).
      stat->path_tm_cubic += 1;
    }
  }

  bool rho_floored = (rho < GR_EULER_DENSITY_FLOOR);
  bool p_floored   = (p   < GR_EULER_PRESSURE_FLOOR);

  // Floor-hit instrumentation (on raw pre-clamp values).
  if (stat) {
    stat->floor_calls += 1;
    if (rho_floored) stat->rho_floor_hits += 1;
    if (p_floored)   stat->p_floor_hits   += 1;
    if (rho_floored && p_floored) stat->both_floor_hits += 1;

    if (rho_floored || p_floored) {
      stat->floor_D_hist[gkyl_gr_euler_status_floor_D_bin(D)] += 1;
      double E_lab = tau + D;
      double s_sq = (D + tau) * (D + tau) - mom_sq;
      double s_sq_over_E2 = (E_lab * E_lab > 1.0e-300)
        ? s_sq / (E_lab * E_lab) : 0.0;
      if (s_sq_over_E2 < 0.0) s_sq_over_E2 = 0.0;
      stat->floor_s2_hist[gkyl_gr_euler_status_floor_s2_bin(s_sq_over_E2)] += 1;
    }
  }

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
