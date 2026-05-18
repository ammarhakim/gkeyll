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
  double h;       // specific enthalpy h = 1 + γ·p / ((γ−1)·ρ)
  bool admissible; // true iff input lay in the strict admissibility set
                   // (D > 0, τ ≥ 0, s² > 0). check_inv consumes this; the
                   // floors below do NOT downgrade it.
};

// Banyuls primitive-variable recovery under Convention A. Inputs are
// undensitized conservatives; inv_g is the 3×3 inverse spatial metric
// γ^{ij} for this cell.
static inline void
gkyl_gr_euler_recover_primitives(
  double gas_gamma,
  double D, double Sx, double Sy, double Sz, double tau,
  const double inv_g[3][3],
  struct gkyl_gr_euler_prim *out)
{
  out->admissible =
    gkyl_gr_euler_check_admissibility(D, Sx, Sy, Sz, tau, inv_g) == GR_EULER_ADM_OK;

  // Lorentz scalar with covariant momentum: |S|² = γ^{ij} S_i S_j.
  double mom_sq = gkyl_gr_euler_mom_sq(inv_g, Sx, Sy, Sz);
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
