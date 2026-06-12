#pragma once

// Shared extreme-state table for the GR-Euler positivity / convexity
// sweeps in ctest_gr_euler_positivity.c (registry, curved-Lax sweep,
// coarse-dx axis) and ctest_gr_euler_riemann.c (cache-equivalence
// sweep). Entries are PRIMITIVES
// (ρ, p, v^i), so the table is convention-agnostic — each consumer
// builds its own conservatives with its own convention and cell metric.
// (Historical: also consumed by ctest_wv_gr_euler_mod.c until the
// iface-flux experiment and the mod variant were retired.)
//
// Provenance:
//   - The mild / shock / vacuum block is the original tetrad-ctest
//     positivity table (g_positivity_cases).
//   - The hard block is the IFACE_FLUX_PLAN.md Phase-B sweep that exposed
//     the iface-flux Lax admissibility failures: near-floor, ultra-
//     relativistic, vacuum-vs-dense, and extreme-contrast regimes — the
//     regimes where the LOW_ORDER fallback fires in real BHL runs.
//
// Velocities are CONTRAVARIANT v^i. Tests must check γ_ij v^i v^j < 1
// against their own cell metrics and skip superluminal combinations:
// several entries are near-luminal in flat space and become superluminal
// in near-horizon metrics. That is intentional coverage, not an error —
// the skip belongs to the test harness, not the table.

struct gr_euler_pos_rp_case {
  const char *name;
  double rho_L, p_L; double v_L[3];
  double rho_R, p_R; double v_R[3];
};

static const struct gr_euler_pos_rp_case gr_euler_pos_cases[] = {
  // --- Mild cases — every scheme should preserve admissibility.
  { "uniform-mild",    1.0,    1.0,    { 0.0,   0.0,  0.0 },  1.0,    1.0,    { 0.0,   0.0, 0.0 } },
  { "equal-mid",       1.0,    1.5,    { 0.10,  0.0,  0.0 },  1.0,    1.5,    { 0.10,  0.0, 0.0 } },
  { "equal-cold",      1.0,    0.01,   { 0.10,  0.0,  0.0 },  1.0,    0.01,   { 0.10,  0.0, 0.0 } },
  { "smooth-shock",    2.0,    1.5,    { 0.10,  0.05, 0.0 },  1.0,    0.7,    { 0.0,   0.0, 0.0 } },
  { "smooth-rare",     1.0,    1.0,    { -0.10, 0.0,  0.0 },  1.0,    1.0,    { 0.10,  0.0, 0.0 } },

  // --- Shocks and contacts.
  { "strong-shock",    3.0,    5.0,    { 0.30,  0.0,  0.0 },  0.5,    0.5,    { -0.30, 0.0, 0.0 } },
  { "shock-mid",       10.0,   1.0,    { 0.30,  0.0,  0.0 },  0.1,    0.05,   { 0.10,  0.0, 0.0 } },
  { "hot-cold",        1.0,    100.0,  { 0.0,   0.0,  0.0 },  1.0,    0.001,  { 0.0,   0.0, 0.0 } },
  { "rel-shock",       5.0,    50.0,   { 0.50,  0.0,  0.0 },  0.1,    0.01,   { -0.50, 0.0, 0.0 } },
  { "strong-contrast", 5.0,    0.1,    { 0.50,  0.0,  0.0 },  0.05,   0.001,  { 0.10,  0.0, 0.0 } },

  // --- Rarefactions and collisions.
  { "strong-rare",     1.0,    1.0,    { -0.50, 0.0,  0.0 },  1.0,    1.0,    { 0.50,  0.0, 0.0 } },
  { "head-on",         1.0,    1.0,    { 0.85,  0.0,  0.0 },  1.0,    1.0,    { -0.85, 0.0, 0.0 } },
  { "ultrarel-head-on",1.0,    1.0,    { 0.95,  0.0,  0.0 },  1.0,    1.0,    { -0.95, 0.0, 0.0 } },
  { "equal-high-v",    1.0,    1.0,    { 0.60,  0.0,  0.0 },  1.0,    1.0,    { 0.60,  0.0, 0.0 } },

  // --- Near-vacuum / near-floor states.
  { "near-vacuum-R",   1.0,    1.0,    { 0.0,   0.0,  0.0 },  1.0e-3, 1.0e-3, { 0.0,   0.0, 0.0 } },
  { "near-vacuum-L",   1.0e-3, 1.0e-3, { 0.0,   0.0,  0.0 },  1.0,    1.0,    { 0.0,   0.0, 0.0 } },
  { "vac-L-vs-dense-R",1.0e-6, 1.0e-7, { 0.10,  0.0,  0.0 },  1.0,    1.0,    { 0.10,  0.0, 0.0 } },
  { "dense-L-vs-vac-R",1.0,    1.0,    { 0.10,  0.0,  0.0 },  1.0e-6, 1.0e-7, { 0.10,  0.0, 0.0 } },
  { "near-floor-equal",1.0e-5, 1.0e-7, { 0.50,  0.0,  0.0 },  1.0e-5, 1.0e-7, { 0.50,  0.0, 0.0 } },
  { "extra-cold",      0.01,   0.0001, { 0.20,  0.0,  0.0 },  0.01,   0.0001, { 0.20,  0.0, 0.0 } },
  { "extreme-contrast",100.0,  100.0,  { 0.30,  0.0,  0.0 },  0.001,  0.0001, { 0.10,  0.0, 0.0 } },

  // --- Cold / supersonic τ/D ≪ 1 (the canonical BHL pathology where the
  //     (τ + p)·v^x flux drives τ negative in the curved frame).
  { "cold-velocity-jump", 1.0,   0.001,  { 0.80, 0.0, 0.0 },  1.0, 0.001,  { 0.10,  0.0, 0.0 } },
  { "cold-supersonic",    1.0,   0.0001, { 0.95, 0.0, 0.0 },  1.0, 0.0001, { 0.10,  0.0, 0.0 } },
  { "cold-shock",         1.0,   1.0e-6, { 0.85, 0.0, 0.0 },  1.0, 1.0e-6, { 0.10,  0.0, 0.0 } },
  { "low-tau/D",          1.0e-3,1.0e-8, { 0.50, 0.0, 0.0 },  1.0, 1.0e-4, { 0.50,  0.0, 0.0 } },
  { "near-cold-v-jump",   1.0,   1.0e-5, { 0.90, 0.0, 0.0 },  1.0, 1.0e-5, { 0.10,  0.0, 0.0 } },
  { "near-cold-reverse",  1.0,   1.0e-5, { 0.90, 0.0, 0.0 },  1.0, 1.0e-5, { -0.10, 0.0, 0.0 } },
};

#define GR_EULER_POS_NUM_CASES \
  ((int)(sizeof(gr_euler_pos_cases) / sizeof(gr_euler_pos_cases[0])))

// Cone-marginal block: the production attractor regime. Cold flow rides
// the p≥0-cone boundary of the invariant domain (s² = D² ⟺ p = 0 —
// SESSION_NOTES_P_GE_0_CONE.md), and sign(s² − D²) = sign(p) means θ = p/ρ
// directly sets how far inside the cone a primitive-built state sits.
// The θ ladder spans fresh-repair landings (θ → 0⁺) up to the
// p-floor scale; the pairs model the BHL wake: uniform cold stream,
// cold velocity jumps, shear, and wake edges against warm/dense and
// near-vacuum neighbors. These rows are what "marginal state" means in
// production under the p≥0-cone domain.
static const struct gr_euler_pos_rp_case gr_euler_pos_cone_marginal_cases[] = {
  // Uniform cold streams at decreasing cone depth.
  { "cone-uniform-th1e-6",  1.0, 1.0e-6,  { 0.30, 0.0, 0.0 },  1.0, 1.0e-6,  { 0.30, 0.0, 0.0 } },
  { "cone-uniform-th1e-9",  1.0, 1.0e-9,  { 0.30, 0.0, 0.0 },  1.0, 1.0e-9,  { 0.30, 0.0, 0.0 } },
  { "cone-uniform-th1e-12", 1.0, 1.0e-12, { 0.30, 0.0, 0.0 },  1.0, 1.0e-12, { 0.30, 0.0, 0.0 } },
  // Cold velocity structure at fixed depth.
  { "cone-vjump-th1e-9",    1.0, 1.0e-9,  { 0.80, 0.0, 0.0 },  1.0, 1.0e-9,  { 0.10, 0.0, 0.0 } },
  { "cone-vjump-rev-th1e-9",1.0, 1.0e-9,  { 0.80, 0.0, 0.0 },  1.0, 1.0e-9,  { -0.10, 0.0, 0.0 } },
  { "cone-shear-th1e-9",    1.0, 1.0e-9,  { 0.50, 0.30, 0.0 }, 1.0, 1.0e-9,  { 0.50, -0.30, 0.0 } },
  { "cone-fast-th1e-9",     1.0, 1.0e-9,  { 0.90, 0.0, 0.0 },  1.0, 1.0e-9,  { 0.90, 0.0, 0.0 } },
  // Wake edges: marginal cell against warm-dense / near-vacuum drivers.
  { "cone-vs-warm",         1.0, 1.0e-9,  { 0.50, 0.0, 0.0 },  1.0, 1.0,     { 0.0, 0.0, 0.0 } },
  { "cone-vs-vacuum",       1.0, 1.0e-9,  { 0.50, 0.0, 0.0 },  1.0e-6, 1.0e-7, { 0.10, 0.0, 0.0 } },
  // Low-density marginal cell against a typical neighbor (post-floor
  // cell at the wake boundary).
  { "cone-low-rho-th1e-9",  1.0e-3, 1.0e-12, { 0.50, 0.0, 0.0 }, 1.0, 1.0e-4, { 0.50, 0.0, 0.0 } },
};

#define GR_EULER_POS_NUM_CONE_MARGINAL_CASES \
  ((int)(sizeof(gr_euler_pos_cone_marginal_cases) / \
         sizeof(gr_euler_pos_cone_marginal_cases[0])))
