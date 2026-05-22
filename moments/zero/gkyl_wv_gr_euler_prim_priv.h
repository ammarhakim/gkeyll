#pragma once

#include <math.h>
#include <stdbool.h>

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

// Pressure floor used by the Newton recovery and the state-repair step.
// Matches the floor that the original inlined Newton solves apply.
#define GR_EULER_PRESSURE_FLOOR (1.0e-8)
#define GR_EULER_DENSITY_FLOOR  (1.0e-8)

// Repair-only τ floor. Used by:
//   (a) the cascade-repair callback in this header — restores τ to this
//       value when τ<0 is detected post-flux-step;
//   (b) the source-step τ limiter in moment_spacetime_coupling.c — the
//       target value the source step is scaled down to land at when
//       its δτ would drive τ below this floor.
// Set equal to GR_EULER_PRESSURE_FLOOR (= 1e-8). Earlier sessions ran
// at 1e-6 to give the cell ~100× margin against a re-trigger before
// the source-step limiter existed; once the limiter is in place
// (SESSION_NOTES_4 §9), it preempts the source step's τ-violation
// directly, so the cascade-repair almost never fires on the source
// path. Empirically the 1e-8 value gives slightly tighter physical
// accuracy (less artificial energy injection per repair) at no cost
// in residual fires under the BHL benchmark.
// Earlier wisdom (still relevant absent the limiter): lowering this to
// 1e-14 INCREASES wave_prop τ<0 fires 54× because the next half-CFL
// flux step can knock τ negative if the margin is below the per-step
// flux scale. 1e-8 is a tested floor that the post-limiter chain still
// handles cleanly.
#define GR_EULER_TAU_REPAIR_FLOOR (1.0e-8)

// ---------------------------------------------------------------------------
// Equation-of-state abstraction.
//
// All EOS-dependent math (specific enthalpy, sound speed, pressure recovery
// in the Newton iteration) routes through the three eos_* helpers below.
// The struct gkyl_gr_euler_eos is a POD bundle stored on the equation
// object and passed by value through the primitive solve, the SR Riemann
// cores, the flux helpers, and the source integrator.
//
//   IDEAL          — p = (γ-1)·ρ·ε  ⇒  h = 1 + γ/(γ-1)·p/ρ
//   MATHEWS_TAUB   — (h - θ)(h - 4θ) = 1, θ = p/ρ. Convex Taub-Mathews
//                    closure that interpolates between the non-relativistic
//                    monatomic limit (Γ=5/3, c_s²=5p/(3ρ)) and the ultra-
//                    relativistic limit (Γ=4/3, c_s²=1/3). Useful for
//                    BHL/shock-heated flows where the gas spans both
//                    regimes within one simulation.
// ---------------------------------------------------------------------------

enum gkyl_gr_euler_eos_type {
  GR_EULER_EOS_IDEAL = 0,
  GR_EULER_EOS_MATHEWS_TAUB = 1,
};

struct gkyl_gr_euler_eos {
  enum gkyl_gr_euler_eos_type type;
  double gas_gamma;  // adiabatic index (consulted only when type == IDEAL)
};

// Specific enthalpy h(ρ, p).
//   IDEAL: h = 1 + γ/(γ-1)·p/ρ
//   MATHEWS_TAUB: h(θ) = (5θ + sqrt(9θ² + 4))/2, θ = p/ρ
static inline double
gkyl_gr_euler_eos_enthalpy(struct gkyl_gr_euler_eos eos,
  double rho, double p)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return 1.0 + (p / rho) * (eos.gas_gamma / (eos.gas_gamma - 1.0));
  }
  double theta = p / rho;
  return 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));
}

// Pressure recovery p(ρ, h). This is the EOS-specific closing equation
// for the Newton-in-pressure iteration in gkyl_gr_euler_recover_primitives.
//   IDEAL: h = 1 + γ/(γ-1)·θ  ⇒  θ = (h-1)·(γ-1)/γ
//   MATHEWS_TAUB: 4θ² − 5hθ + (h² − 1) = 0  ⇒  θ = (5h − sqrt(9h² + 16))/8
static inline double
gkyl_gr_euler_eos_pressure_from_rho_h(struct gkyl_gr_euler_eos eos,
  double rho, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return rho * (h - 1.0) * (eos.gas_gamma - 1.0) / eos.gas_gamma;
  }
  double theta = 0.125 * (5.0 * h - sqrt(9.0 * h * h + 16.0));
  return rho * theta;
}

// Sound speed squared c_s²(ρ, p, h). For ideal this is γp/(ρh). For
// Mathews-Taub the closed form is c_s² = θ(5h − 8θ) / (3·h·(h − θ)) with
// θ = p/ρ. Verified limits: θ→0 ⇒ c_s² → 5p/(3ρ)/h (matches Γ=5/3 ideal),
// θ→∞ ⇒ c_s² → 1/3 (radiation-fluid asymptote).
static inline double
gkyl_gr_euler_eos_cs2(struct gkyl_gr_euler_eos eos,
  double rho, double p, double h)
{
  if (eos.type == GR_EULER_EOS_IDEAL) {
    return (eos.gas_gamma * p) / (rho * h);
  }
  double theta = p / rho;
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

// One Eulderink-Mellema quartic Newton pass at fixed γ. Recovers
// (ρ, v^i, p, W, h) from undensitized (D, S_i, τ) and inv_g.
//
// The "guess" iterate is 1/(h·C), C = D/√s², which gives a 4th-order
// polynomial residual that converges quadratically from guess=1.0 for
// any admissible input. Outputs are pre-floor so callers can apply
// floors / dispatch after refinement.
static inline void
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
//                   range (5/3 at θ→0, 4/3 at θ→∞), and is a smooth
//                   contraction map of (ρ, p), so Picard converges in 2-3
//                   outer iterations for any TM state. Each EM Newton pass
//                   stays well-conditioned (the same quartic that production
//                   IDEAL relies on), so TM inherits IDEAL's robust
//                   floor-region behavior — unlike a Newton-in-Z formulation
//                   whose iteration variable has poor scaling near low-p
//                   floors and near v² → 1.
//
// Both paths recover the same convention (contravariant v^i, post-floor
// ρ, p) and populate the same out struct.
static inline void
gkyl_gr_euler_recover_primitives(
  struct gkyl_gr_euler_eos eos,
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
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
    gkyl_gr_euler_em_newton_at_gamma(
      eos.gas_gamma, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
      &rho, &vx, &vy, &vz, &p, &W, &h);
  }
  else {
    // Picard iteration: outer loop refines γ_eff to match h_TM at the
    // current (ρ, p); inner loop is a full EM Newton at that γ_eff.
    //
    // Initial γ_eff = 5/3: this is the cold-flow TM limit and gives a
    // converged answer in 1 outer iteration for any cell with θ ≪ 1.
    // Hot cells converge in 2–3 outer iterations.
    double gamma_eff = 5.0 / 3.0;
    for (int picard = 0; picard < 20; picard++) {
      gkyl_gr_euler_em_newton_at_gamma(
        gamma_eff, mom_sq, D, Sx, Sy, Sz, tau, inv_g,
        &rho, &vx, &vy, &vz, &p, &W, &h);

      // Compute TM-matching γ_eff from the current (ρ, p). Guard against
      // θ → 0 (where the matching identity has a 0/0 limit) — in that
      // limit γ_eff = 5/3 exactly (TM cold-flow). Also guard against
      // non-physical p ≤ 0 (recovered just below the floor) — treat as
      // cold-flow.
      double rho_safe = (rho > GR_EULER_DENSITY_FLOOR) ? rho : GR_EULER_DENSITY_FLOOR;
      double p_safe   = (p   > GR_EULER_PRESSURE_FLOOR) ? p   : GR_EULER_PRESSURE_FLOOR;
      double theta = p_safe / rho_safe;
      double gamma_new;
      if (theta < 1.0e-12) {
        gamma_new = 5.0 / 3.0;
      } else {
        double h_tm = 0.5 * (5.0 * theta + sqrt(9.0 * theta * theta + 4.0));
        double f = (h_tm - 1.0) / theta;
        // f ∈ [5/2, 4) for TM ⇒ f − 1 ∈ [3/2, 3), bounded away from zero.
        gamma_new = f / (f - 1.0);
      }

      if (fabs(gamma_new - gamma_eff) < 1.0e-14) {
        gamma_eff = gamma_new;
        break;
      }
      gamma_eff = gamma_new;
    }
  }

  if (rho < GR_EULER_DENSITY_FLOOR)  rho = GR_EULER_DENSITY_FLOOR;
  if (p   < GR_EULER_PRESSURE_FLOOR) p   = GR_EULER_PRESSURE_FLOOR;

  out->rho = rho;
  out->v[0] = vx;
  out->v[1] = vy;
  out->v[2] = vz;
  out->p = p;
  out->W = W;
  out->h = h;
}
