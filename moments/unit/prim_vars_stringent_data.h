#pragma once

// Shared parameter tables for the four GR Euler prim_vars stringent stress
// tests (ctest_wv_gr_euler.c, ctest_wv_gr_euler_mod.c,
// ctest_wv_gr_euler_tetrad.c, ctest_wv_gr_euler_tetrad.c). Keeping the
// data in one place ensures the four test runners exercise the same regime
// — pressure sweep, density sweep, velocity sweep, and combined extremes —
// at the same sample positions (near-horizon points where γ_ij departs
// strongly from δ_ij and β^i is large for spinning Kerr, mid-range, and
// far-field).

struct prim_state_strict {
  double rho, p;
  double vx, vy, vz;
};

static const struct prim_state_strict stringent_states[] = {
  // Pressure sweep over 6 orders of magnitude (mild velocity, ρ=1)
  { 1.0,    1.0e-3, 0.10, 0.05, 0.00 },  // cold gas
  { 1.0,    1.0e-1, 0.10, 0.05, 0.00 },
  { 1.0,    1.0,    0.10, 0.05, 0.00 },  // baseline (existing test regime)
  { 1.0,    1.0e+1, 0.10, 0.05, 0.00 },
  { 1.0,    1.0e+3, 0.10, 0.05, 0.00 },  // pressure-dominated (relativistic h)
  // Density sweep (4 orders of magnitude)
  { 1.0e-2, 1.0,    0.10, 0.05, 0.00 },  // low density
  { 1.0e+2, 1.0,    0.10, 0.05, 0.00 },  // high density
  // Velocity sweep (in flat magnitude; γ-norm may exceed when near horizon)
  { 1.0,    1.0,    0.30, 0.20, 0.10 },  // |v|_flat ≈ 0.37
  { 1.0,    1.0,    0.50, 0.30, 0.20 },  // |v|_flat ≈ 0.62
  { 1.0,    1.0,    0.60, 0.40, 0.30 },  // |v|_flat ≈ 0.78
  // Combined extremes
  { 0.1,    100.0,  0.20, 0.10, 0.05 },  // hot + low ρ + moderate v
  { 10.0,   0.01,   0.05, 0.02, 0.01 },  // cold + high ρ + slow
  { 0.01,   1000.0, 0.10, 0.10, 0.10 },  // pressure-dominated, low ρ, moderate v
};
#define N_STRINGENT_STATES \
  ((int)(sizeof(stringent_states) / sizeof(stringent_states[0])))

// Sample positions chosen to span:
//   - Near-horizon (Schwarzschild r_+ = 2M = 0.2 for M=0.1) where γ_xx ≫ 1
//     and α ≪ 1
//   - Off-axis near-horizon — γ acquires nonzero off-diagonal entries
//   - Mid-range — moderate curvature
//   - Far-field — nearly flat
static const double stringent_positions[][3] = {
  { 0.25, 0.00, 0.0 },   // just outside Schwarzschild horizon along +x
  { 0.30, 0.00, 0.0 },
  { 0.00, 0.25, 0.0 },   // along +y (rotational symmetry check)
  { 0.22, 0.22, 0.0 },   // diagonal near horizon (off-diagonal γ_xy)
  { 0.21, 0.10, 0.0 },   // very close on one axis
  { 0.50, 0.10, 0.0 },   // intermediate
  { 1.00, 1.00, 0.0 },   // mid-range
  { 5.00, 0.00, 0.0 },   // far-field
  { 10.0, 10.0, 0.0 },   // very far (nearly flat)
};
#define N_STRINGENT_POSITIONS \
  ((int)(sizeof(stringent_positions) / sizeof(stringent_positions[0])))

// Relative-error tolerance for the stringent stress-tests across all listed
// (state, position, spacetime) combinations.
#define STRINGENT_REL_TOL 1.0e-10

// Floor for the relative-error denominator: avoids amplifying noise on
// near-zero primitives. Matches the floor on ρ and p inside prim_vars.
#define STRINGENT_REL_FLOOR 1.0e-8

// ---------------------------------------------------------------------------
// Combined 3D ultra-stress sweep: full cross-product of W × ρ × (p/ρ),
// capped at W = 1e4 (the algorithmic-precision ceiling of the recovery
// polynomial — see the W²·ε analysis in the run_ultra_combined docstring).
// Systematically explores the joint parameter space and catches joint
// failure modes that single-axis slices would miss. Exercised in Minkowski,
// near-horizon Schwarzschild, and near-horizon Kerr.
// ---------------------------------------------------------------------------
static const double ultra_combined_W[] = {
  10.0, 100.0, 1.0e+3, 1.0e+4
};
#define N_ULTRA_COMBINED_W \
  ((int)(sizeof(ultra_combined_W) / sizeof(ultra_combined_W[0])))

static const double ultra_combined_rho[] = {
  1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e+2, 1.0e+4
};
#define N_ULTRA_COMBINED_RHO \
  ((int)(sizeof(ultra_combined_rho) / sizeof(ultra_combined_rho[0])))

static const double ultra_combined_p_over_rho[] = {
  1.0e-6, 1.0e-3, 1.0, 1.0e+3, 1.0e+6
};
#define N_ULTRA_COMBINED_PR \
  ((int)(sizeof(ultra_combined_p_over_rho) / sizeof(ultra_combined_p_over_rho[0])))
