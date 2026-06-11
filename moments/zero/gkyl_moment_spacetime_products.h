#pragma once

// Per-cell layout of the "spacetime products" gkyl_array. This array is filled
// by the moment_spacetime_coupling derive-products step from either an
// analytic gkyl_gr_spacetime callback set or the evolved Bona-Masso Einstein
// state. The modular GR fluid equations (wv_gr_euler_tetrad)
// acquire pointers to this array via auxfields and read the components below.
//
// Indices 0..66 reproduce the layout of the packed GR-Euler state q[5..66]
// re-indexed to start at 0, so source-term math ports across by mechanical
// substitution. Indices 62..71 (inv-metric, det) are precomputed here rather
// than recomputed inline at every source-update call.

enum {
  GKYL_GR_SP_LAPSE       =  0,  // 1   — α
  GKYL_GR_SP_SHIFT       =  1,  // 3   — β^i
  GKYL_GR_SP_GIJ         =  4,  // 9   — γ_ij row-major
  GKYL_GR_SP_KIJ         = 13,  // 9   — K_ij row-major
  GKYL_GR_SP_EXCISION    = 22,  // 1   — sign-coded flag (matches packed q[27])
  GKYL_GR_SP_DALPHA      = 23,  // 3   — ∂_i α
  GKYL_GR_SP_DBETA       = 26,  // 9   — ∂_j β^i row-major over (j,i)
  GKYL_GR_SP_DGIJ        = 35,  // 27  — ∂_k γ_ij row-major over (k,i,j)
  GKYL_GR_SP_INV_GIJ     = 62,  // 9   — γ^ij row-major (derived)
  GKYL_GR_SP_SPATIAL_DET = 71,  // 1   — det γ (derived)
  GKYL_GR_SP_CELLCTR     = 72,  // 3   — (x, y, z) cell center
  GKYL_GR_SP_NCOMP_BASE  = 75,

  // Tetrad add-on. Only allocated when a tetrad-mod species is registered
  // (gkyl_moment_spacetime_coupling input flag has_tetrad).
  GKYL_GR_SP_TETRAD       =  75, // 16 — e^μ_a row-major over (μ,a)
  GKYL_GR_SP_TETRAD_INV   =  91, // 16 — E^a_μ row-major over (a,μ)
  GKYL_GR_SP_NCOMP_TETRAD = 107
};
