#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad_mod.h>
#include <gkyl_wv_gr_euler_tetrad_mod_priv.h>

// ---------------------------------------------------------------------------
// Setters for cell-index communication from wave_prop. Identical to the
// regular mod variant; the per-instance idx slots populated below are read
// by downstream callbacks to fetch the appropriate row from the spacetime
// products array.
// ---------------------------------------------------------------------------

static void
gr_euler_tetrad_mod_set_interface_idx(const struct gkyl_wv_eqn *eqn,
  const int *idxl, const int *idxr)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = idxl[d];
    grm->cur_idxr[d] = idxr[d];
  }
  grm->rot_call_parity = 0;
}

static void
gr_euler_tetrad_mod_set_cell_idx(const struct gkyl_wv_eqn *eqn, const int *idx)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_cell_idx[d] = idx[d];
  }
}

void
gkyl_gr_euler_tetrad_mod_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_tetrad_mod_auxfields auxin)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  grm->auxfields.prods = auxin.prods;
}

void
gkyl_gr_euler_tetrad_mod_set_conf_range(const struct gkyl_wv_eqn *eqn,
  const struct gkyl_range *conf_range)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  grm->conf_range = *conf_range;
}

// ---------------------------------------------------------------------------
// Rotation of the spacetime products array from global to local
// tangent-normal frame. Shared verbatim with the regular mod variant — the
// tetrad equation reads the same spacetime block; rotation conventions are
// identical.
// ---------------------------------------------------------------------------

static inline void
rotate_rank1(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  out[0] = in[0]*norm[0] + in[1]*norm[1] + in[2]*norm[2];
  out[1] = in[0]*tau1[0] + in[1]*tau1[1] + in[2]*tau1[2];
  out[2] = in[0]*tau2[0] + in[1]*tau2[1] + in[2]*tau2[2];
}

static inline void
rotate_rank2(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  double tmp[3][3];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 3; i++)
      tmp[a][i] = in[3*i+0]*R[a][0] + in[3*i+1]*R[a][1] + in[3*i+2]*R[a][2];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      out[3*a + b] = tmp[a][0]*R[b][0] + tmp[a][1]*R[b][1] + tmp[a][2]*R[b][2];
}

static inline void
rotate_rank3(const double *tau1, const double *tau2, const double *norm,
  const double *in, double *out)
{
  const double *R[3] = { norm, tau1, tau2 };
  double tmp1[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        tmp1[a][j][k] =
          in[0*9 + 3*j + k]*R[a][0] +
          in[1*9 + 3*j + k]*R[a][1] +
          in[2*9 + 3*j + k]*R[a][2];
  double tmp2[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      for (int k = 0; k < 3; k++)
        tmp2[a][b][k] =
          tmp1[a][0][k]*R[b][0] +
          tmp1[a][1][k]*R[b][1] +
          tmp1[a][2][k]*R[b][2];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      for (int c = 0; c < 3; c++)
        out[9*a + 3*b + c] =
          tmp2[a][b][0]*R[c][0] +
          tmp2[a][b][1]*R[c][1] +
          tmp2[a][b][2]*R[c][2];
}

static inline void
rot_spacetime_to_local(const double *tau1, const double *tau2,
  const double *norm, const double *in_p, double *out_p)
{
  out_p[GKYL_GR_SP_LAPSE]        = in_p[GKYL_GR_SP_LAPSE];
  out_p[GKYL_GR_SP_EXCISION]     = in_p[GKYL_GR_SP_EXCISION];
  out_p[GKYL_GR_SP_SPATIAL_DET]  = in_p[GKYL_GR_SP_SPATIAL_DET];

  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_SHIFT],   &out_p[GKYL_GR_SP_SHIFT]);
  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_DALPHA],  &out_p[GKYL_GR_SP_DALPHA]);
  rotate_rank1(tau1, tau2, norm, &in_p[GKYL_GR_SP_CELLCTR], &out_p[GKYL_GR_SP_CELLCTR]);

  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_GIJ],     &out_p[GKYL_GR_SP_GIJ]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_KIJ],     &out_p[GKYL_GR_SP_KIJ]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_DBETA],   &out_p[GKYL_GR_SP_DBETA]);
  rotate_rank2(tau1, tau2, norm, &in_p[GKYL_GR_SP_INV_GIJ], &out_p[GKYL_GR_SP_INV_GIJ]);

  rotate_rank3(tau1, tau2, norm, &in_p[GKYL_GR_SP_DGIJ],    &out_p[GKYL_GR_SP_DGIJ]);
}

// ---------------------------------------------------------------------------
// Hydro helpers (prim_vars / flux / flux_correction / max_abs_speed).
// prim_vars and max_abs_speed match the regular mod variant — they don't
// see the tetrad split. The flux is what differs:
//   gkyl_gr_euler_tetrad_mod_flux:            flat-space SR flux (Cartesian W,
//                                             no lapse/shift/√γ)
//   gkyl_gr_euler_tetrad_mod_flux_correction: applies α·√γ and replaces
//                                             (vx, W_flat) with (vx - βˣ/α,
//                                             W_curved) to recover the
//                                             curved-space flux.
// Mirrors the packed factorization in wv_gr_euler_tetrad.c at lines 10-124.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_mod_prim_vars(double gas_gamma, const double q[5],
  const double *prods, double v[5])
{
  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];

    double D    = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);
    double Etot = q[4] / sqrt(spatial_det);

    // Convention A: q[1..3] is genuine covariant momentum S_i. Recovery
    // contracts |S|² with γ^{ij} and raises the velocity with γ^{ij};
    // both lookups go through the shared helper which is the single
    // source of truth for the Banyuls Newton solve.
    const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
    double inv_g[3][3] = {
      { ig[0], ig[1], ig[2] },
      { ig[3], ig[4], ig[5] },
      { ig[6], ig[7], ig[8] },
    };
    struct gkyl_gr_euler_prim prim;
    gkyl_gr_euler_recover_primitives(gas_gamma,
      D, momx, momy, momz, Etot, inv_g, &prim);

    v[0] = prim.rho;
    v[1] = prim.v[0];
    v[2] = prim.v[1];
    v[3] = prim.v[2];
    v[4] = prim.p;
  } else {
    for (int i = 0; i < 5; i++) v[i] = 0.0;
  }
}

void
gkyl_gr_euler_tetrad_mod_flux(double gas_gamma, const double q[5],
  const double *prods, double flux_sr[5])
{
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double v_dot = (vx*vx) + (vy*vy) + (vz*vz);
    double W = 1.0 / sqrt(1.0 - v_dot);
    if (v_dot > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(pow(10.0, -8.0));

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux_sr[0] = rho * W * vx;
    flux_sr[1] = (rho * h * (W*W) * (vx * vx)) + p;
    flux_sr[2] = rho * h * (W*W) * (vy * vx);
    flux_sr[3] = rho * h * (W*W) * (vz * vx);
    flux_sr[4] = ((rho * h * (W*W)) - (rho * W)) * vx;
  } else {
    for (int i = 0; i < 5; i++) flux_sr[i] = 0.0;
  }
}

void
gkyl_gr_euler_tetrad_mod_flux_correction(double gas_gamma, const double q[5],
  const double *prods, const double flux_sr[5], double flux_gr[5])
{
  // Mirrors gkyl_gr_euler_tetrad_flux_correction in the packed implementation:
  // scales the flat SR flux by the W_curved / W_flat ratios that arise when
  // mapping the tetrad-frame flux back into the coord frame, plus the α·√γ
  // densitization and shift correction. We deliberately keep this byte-for-
  // byte equivalent to packed so the *only* algorithmic difference between
  // packed and mod-tetrad lives in the Roe solve. See wv_gr_euler_tetrad.c
  // for the parallel implementation.
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  double lapse        = prods[GKYL_GR_SP_LAPSE];
  double shift_x      = prods[GKYL_GR_SP_SHIFT + 0];
  double spatial_det  = prods[GKYL_GR_SP_SPATIAL_DET];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * v[i+1] * v[j+1];

    double v_dot = (vx*vx) + (vy*vy) + (vz*vz);
    double W_flat = 1.0 / sqrt(1.0 - v_dot);
    if (v_dot > 1.0 - pow(10.0, -8.0)) W_flat = 1.0 / sqrt(pow(10.0, -8.0));

    double W_curved = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) W_curved = 1.0 / sqrt(pow(10.0, -8.0));

    // Numerical-safety clamp on vx for the divisions below. Mirrors packed.
    if (fabs(vx) < pow(10.0, -8.0)) {
      vx = (vx > 0.0) ? pow(10.0, -8.0) : -pow(10.0, -8.0);
    }

    // Convention A momentum flux: F^x[S_i] = α√γ·(S_i·v̂^x + p·δ_i^x)
    // with S_i·v̂^x = ρhW²·v_l[i]·v̂^x and v_l[i] = γ_ij·v^j (lowered
    // contravariant velocity). Reconstruct ρhW²_curved from flux_sr[1] =
    // (ρhW²_flat·v^x·v^x + p) so the D and τ slots can keep the existing
    // factorization while the momentum slots use the explicit Banyuls
    // form.
    double v_l[3];
    for (int i = 0; i < 3; i++) {
      v_l[i] = prods[GKYL_GR_SP_GIJ + 3*i + 0]*vx
             + prods[GKYL_GR_SP_GIJ + 3*i + 1]*vy
             + prods[GKYL_GR_SP_GIJ + 3*i + 2]*vz;
    }
    double rhohW2_c = ((flux_sr[1] - p) * (W_curved * W_curved))
                    / ((vx * vx) * (W_flat * W_flat));

    double prefac = lapse * sqrt(spatial_det);
    double vmsh   = vx - (shift_x / lapse);

    flux_gr[0] = prefac * ((flux_sr[0] * vmsh * W_curved) / (vx * W_flat));
    flux_gr[1] = prefac * (rhohW2_c * v_l[0] * vmsh + p);
    flux_gr[2] = prefac * (rhohW2_c * v_l[1] * vmsh);
    flux_gr[3] = prefac * (rhohW2_c * v_l[2] * vmsh);
    flux_gr[4] = prefac * (((((flux_sr[4] + (rho * vx * W_flat)) * (W_curved*W_curved))
                          / (vx * (W_flat*W_flat))) - p - (rho * W_curved)) * vmsh + (p * vx));
  } else {
    for (int i = 0; i < 5; i++) flux_gr[i] = 0.0;
  }
}

double
gkyl_gr_euler_tetrad_mod_max_abs_speed(double gas_gamma, const double q[5],
  const double *prods)
{
  double v[5];
  gkyl_gr_euler_tetrad_mod_prim_vars(gas_gamma, q, prods, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];

  double num = (gas_gamma * p) / rho;
  double den = 1.0 + ((p / rho) * gas_gamma / (gas_gamma - 1.0));
  double c_s = sqrt(num / den);

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double g_ij = prods[GKYL_GR_SP_GIJ + 3*i + j];
      if (i == j) {
        if (fabs(g_ij - 1.0) > pow(10.0, -8.0)) curved_spacetime = true;
      } else {
        if (fabs(g_ij)       > pow(10.0, -8.0)) curved_spacetime = true;
      }
    }
  }
  if (fabs(lapse - 1.0) > pow(10.0, -8.0) ||
      fabs(shift_x)     > pow(10.0, -8.0) ||
      fabs(shift_y)     > pow(10.0, -8.0) ||
      fabs(shift_z)     > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (in_excision_region) return pow(10.0, -8.0);

  if (curved_spacetime) {
    double vel[3]   = { vx, vy, vz };
    double shift[3] = { shift_x, shift_y, shift_z };
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];

    double max_eig = 0.0;
    for (int i = 0; i < 3; i++) {
      double inv_ii = prods[GKYL_GR_SP_INV_GIJ + 3*i + i];
      double material = (lapse * vel[i]) - shift[i];
      double common = lapse / (1.0 - (v_sq * (c_s*c_s)));
      double rad = (1.0 - v_sq) * (inv_ii * (1.0 - (v_sq * (c_s*c_s))) -
        (vel[i] * vel[i]) * (1.0 - (c_s*c_s)));
      double fast = common * ((vel[i] * (1.0 - (c_s*c_s))) + (c_s * sqrt(rad))) - shift[i];
      double slow = common * ((vel[i] * (1.0 - (c_s*c_s))) - (c_s * sqrt(rad))) - shift[i];

      if (fabs(material) > max_eig) max_eig = fabs(material);
      if (fabs(fast)     > max_eig) max_eig = fabs(fast);
      if (fabs(slow)     > max_eig) max_eig = fabs(slow);
    }

    return fabs(v_sq) + max_eig;
  }

  double v_mag = sqrt((vx*vx) + (vy*vy) + (vz*vz));
  return fabs(v_mag) + c_s;
}

// ---------------------------------------------------------------------------
// Modular tetrad-Roe pipeline. See header doc in
// gkyl_wv_gr_euler_tetrad_mod_priv.h for the contracts of each piece.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_mod_build_triad(const double g_ij[3][3],
  double L[3][3], double L_inv[3][3])
{
  // Cholesky factor of γ: γ = L L^T, L lower triangular. L_inv is the
  // inverse Cholesky factor. With these, the orthonormal triad ε^i_a =
  // (L^{-1})^T_{ia} satisfies γ_ij ε^i_a ε^j_b = δ_ab.
  //
  // We tried the alternative J.10 construction (Cholesky of γ^{-1}) — it
  // gives a different but equivalent tetrad with similar properties and
  // similar residual in off-diagonal-γ tests. Neither construction
  // eliminates the off-diagonal-γ flux-jump residual: that residual is a
  // structural limitation of 1D-sweep tetrad-Roe with non-diagonal coord-
  // frame γ. The SR Roe solves a 1D problem along the tetrad e_0 direction,
  // and the back-transformed flux F^0_coord depends on contributions from
  // tangent directions through γ_0j coupling, which a single 1D solve
  // cannot capture. See priv.h doc for the full analysis.
  L[0][0] = sqrt(g_ij[0][0]);
  L[1][0] = g_ij[1][0] / L[0][0];
  L[1][1] = sqrt(g_ij[1][1] - L[1][0]*L[1][0]);
  L[2][0] = g_ij[2][0] / L[0][0];
  L[2][1] = (g_ij[2][1] - L[2][0]*L[1][0]) / L[1][1];
  L[2][2] = sqrt(g_ij[2][2] - L[2][0]*L[2][0] - L[2][1]*L[2][1]);
  L[0][1] = L[0][2] = L[1][2] = 0.0;

  // Invert L by forward substitution (lower triangular).
  L_inv[0][0] = 1.0 / L[0][0];
  L_inv[1][1] = 1.0 / L[1][1];
  L_inv[2][2] = 1.0 / L[2][2];
  L_inv[1][0] = -L[1][0] * L_inv[0][0] / L[1][1];
  L_inv[2][1] = -L[2][1] * L_inv[1][1] / L[2][2];
  L_inv[2][0] = -(L[2][0] * L_inv[0][0] + L[2][1] * L_inv[1][0]) / L[2][2];
  L_inv[0][1] = L_inv[0][2] = L_inv[1][2] = 0.0;
}

// Build an orthonormal triad whose FIRST basis vector is aligned with the
// CONTRAVARIANT x-direction (∂^x = γ^{xj}·∂_j) rather than the coordinate
// x-direction. The other two basis vectors come from Gram-Schmidt-in-γ
// starting from (0,1,0) and (0,0,1). This eliminates the v_tet^x ↔ v^y,v^z
// mixing seen with Cholesky-on-γ for non-diagonal γ.
//
// Outputs:
//   M[i][a]    = e_a^i  (matrix of basis-vector components; columns = e_a)
//   M_inv[a][i] = (e^a)_i = γ_ij · M[j][a]  (dual co-vector components)
//
// Properties:
//   γ_ij·M[i][a]·M[j][b] = δ_ab          (orthonormality)
//   M_inv·M = M·M_inv = I (matrix product)
//   M[i][0] = γ^{xi}/√γ^{xx}              (contravariant x alignment)
//   M_inv[0][i] = δ^x_i/√γ^{xx} = (1/√γ^{xx}, 0, 0)
//
// For diagonal γ this reduces to the Cholesky construction.
void
gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(
  const double g_ij[3][3], const double inv_g[3][3],
  double M[3][3], double M_inv[3][3])
{
  // First basis vector: contravariant x-direction, normalized in γ.
  //   e_0^i = γ^{xi} / √γ^{xx}
  double sqrt_inv_gxx = sqrt(inv_g[0][0]);
  for (int i = 0; i < 3; i++) M[i][0] = inv_g[0][i] / sqrt_inv_gxx;

  // Helper: γ-inner-product of two contravariant vectors.
  // <u, v>_γ = γ_ij · u^i · v^j
  #define GAMMA_DOT(u, v)                                              \
    (g_ij[0][0]*(u)[0]*(v)[0] + g_ij[1][1]*(u)[1]*(v)[1]               \
   + g_ij[2][2]*(u)[2]*(v)[2]                                          \
   + g_ij[0][1]*((u)[0]*(v)[1] + (u)[1]*(v)[0])                        \
   + g_ij[0][2]*((u)[0]*(v)[2] + (u)[2]*(v)[0])                        \
   + g_ij[1][2]*((u)[1]*(v)[2] + (u)[2]*(v)[1]))

  // Second basis vector: Gram-Schmidt-in-γ of (0,1,0) against e_0.
  double e0[3] = { M[0][0], M[1][0], M[2][0] };
  double f1[3] = { 0.0, 1.0, 0.0 };
  double c01 = GAMMA_DOT(f1, e0);  // projection coefficient
  double e1[3] = { f1[0] - c01*e0[0], f1[1] - c01*e0[1], f1[2] - c01*e0[2] };
  double norm1 = sqrt(GAMMA_DOT(e1, e1));
  for (int i = 0; i < 3; i++) M[i][1] = e1[i] / norm1;

  // Third basis vector: Gram-Schmidt of (0,0,1) against e_0 and e_1.
  double e1_norm[3] = { M[0][1], M[1][1], M[2][1] };
  double f2[3] = { 0.0, 0.0, 1.0 };
  double c02 = GAMMA_DOT(f2, e0);
  double c12 = GAMMA_DOT(f2, e1_norm);
  double e2[3] = {
    f2[0] - c02*e0[0] - c12*e1_norm[0],
    f2[1] - c02*e0[1] - c12*e1_norm[1],
    f2[2] - c02*e0[2] - c12*e1_norm[2],
  };
  double norm2 = sqrt(GAMMA_DOT(e2, e2));
  for (int i = 0; i < 3; i++) M[i][2] = e2[i] / norm2;

  #undef GAMMA_DOT

  // M_inv = M^T · γ (since γ-orthonormality says M^T·γ·M = I, so
  // M^{-1} = M^T·γ).
  for (int a = 0; a < 3; a++) {
    for (int i = 0; i < 3; i++) {
      M_inv[a][i] = M[0][a]*g_ij[0][i]
                  + M[1][a]*g_ij[1][i]
                  + M[2][a]*g_ij[2][i];
    }
  }
}

// Forward transform of covariant momentum (Convention A) onto the
// contravariant-x-aligned tetrad. S_tet^a = M_inv[a][i]·γ^{ij}·S_j_curved.
// For a=0 this simplifies to (1/√γ^{xx})·S^x_curved, so v_tet^0 = v^x/√γ^{xx}
// (no metric mixing).
void
gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(const double q_GR[5],
  double sqrt_det, const double inv_g[3][3], const double M_inv[3][3],
  double q_tet[5])
{
  q_tet[0] = q_GR[0] / sqrt_det;
  q_tet[4] = q_GR[4] / sqrt_det;

  // Raise the index: S^i = γ^{ij}·S_j. Then project: S_tet^a = M_inv[a][i]·S^i.
  // We fold both steps to write S_tet^a = (M_inv · γ^{-1})[a][j] · S_j.
  // For efficiency, compute S^i first.
  double Sup[3];
  for (int i = 0; i < 3; i++) {
    Sup[i] = (inv_g[i][0]*q_GR[1] + inv_g[i][1]*q_GR[2]
            + inv_g[i][2]*q_GR[3]) / sqrt_det;
  }
  for (int a = 0; a < 3; a++) {
    q_tet[a+1] = M_inv[a][0]*Sup[0] + M_inv[a][1]*Sup[1] + M_inv[a][2]*Sup[2];
  }
}

// Back-transform of waves from the contravariant-x-aligned tetrad to
// curved-frame Convention A (covariant momentum). The wave amplitude on
// the i-th coord-momentum slot is:
//   w_GR[i+1] = √γ · γ_ij · M[j][a] · w_tet[a+1]
//             = √γ · M_inv^T[a][i] · w_tet[a+1]
// because M_inv[a][i] = γ_ij·M[j][a], so M_inv^T[i][a] = γ_ij·M[j][a].
//
// For the normal (i=0) momentum slot, M_inv^T[0][a] = M_inv[a][0]. And
// M_inv[0][0] = 1/√γ^{xx}, M_inv[a][0] for a>0 = 0 (by construction —
// the only "covariant-x" content lives in e_0, since e_1, e_2 are
// γ-orthogonal to e_0). So the back-transform for coord-x momentum is
//   w_GR[1] = √γ · (1/√γ^{xx}) · w_tet[1]
// clean: only the a=0 tetrad wave contributes to S_x. (This is the
// mirror of the forward-clean property in Cholesky for Convention A.)
void
gkyl_gr_euler_tetrad_mod_wave_to_curved_contra(const double w_tet[5],
  double sqrt_det, const double M_inv[3][3], double w_GR[5])
{
  w_GR[0] = sqrt_det * w_tet[0];
  w_GR[4] = sqrt_det * w_tet[4];

  // w_GR[i+1] = √γ · M_inv^T[a][i] · w_tet[a+1] = √γ · M_inv[a][i] · w_tet[a+1]
  // (with sum over a).
  for (int i = 0; i < 3; i++) {
    w_GR[i+1] = sqrt_det * (M_inv[0][i]*w_tet[1]
                          + M_inv[1][i]*w_tet[2]
                          + M_inv[2][i]*w_tet[3]);
  }
}

// Speed back-transform for the contravariant-x triad. The wave at tetrad
// speed s_tet propagates in the e_0 direction; the coord-frame x-speed
// is α·(e_0^x)·s_tet − β^x where e_0^x = M[0][0] = γ^{xx}/√γ^{xx} =
// √γ^{xx}.
//   s_coord = α · √γ^{xx} · s_tet − β^x
// For diagonal γ, √γ^{xx} = 1/√γ_xx = L_inv[0][0]_Cholesky, so the formula
// reduces to the Cholesky speed_to_curved.
double
gkyl_gr_euler_tetrad_mod_speed_to_curved_contra(double s_tet,
  double lapse, double shift_x, double inv_gxx)
{
  return (lapse * sqrt(inv_gxx) * s_tet) - shift_x;
}

void
gkyl_gr_euler_tetrad_mod_q_to_tetrad(const double q_GR[5],
  double sqrt_det, const double L_inv[3][3], double q_tet[5])
{
  // Convention A forward transform: project covariant momentum onto the
  // orthonormal triad basis ε^i_a = L_inv[a][i] (so γ_ij·ε^i_a·ε^j_b =
  // δ_ab). For a covariant vector S_i,
  //   S_a_tet = ε^i_a · S_i = L_inv[a][i] · (q_GR[i+1]/√γ)
  // and |S_tet|² = δ^{ab} S_a S_b = γ^{ij}·S_i·S_j matches the Lorentz
  // scalar in the curved frame.
  q_tet[0] = q_GR[0] / sqrt_det;
  q_tet[4] = q_GR[4] / sqrt_det;

  for (int a = 0; a < 3; a++) {
    q_tet[a+1] = (L_inv[a][0]*q_GR[1]
                + L_inv[a][1]*q_GR[2]
                + L_inv[a][2]*q_GR[3]) / sqrt_det;
  }
}

void
gkyl_gr_euler_tetrad_mod_wave_to_curved(const double w_tet[5],
  double sqrt_det, const double L[3][3], double w_GR[5])
{
  // Inverse of the forward transform S_a_tet = L_inv[a][i]·S_i:
  //   S_i = L[i][a]·S_a_tet   (sum over a)
  // so w_GR[i+1] = √γ·L[i][a]·w_tet[a+1]. With L lower triangular the
  // i=0 slot only sees tetrad-x; i=1 picks up tetrad-x and -y; i=2 all
  // three tetrad axes. This cross-coupling carries the off-diagonal-γ
  // contributions back into the covariant momentum slots.
  w_GR[0] = sqrt_det * w_tet[0];
  for (int i = 0; i < 3; i++) {
    w_GR[i+1] = sqrt_det * (L[i][0]*w_tet[1]
                          + L[i][1]*w_tet[2]
                          + L[i][2]*w_tet[3]);
  }
  w_GR[4] = sqrt_det * w_tet[4];
}

double
gkyl_gr_euler_tetrad_mod_speed_to_curved(double s_tet,
  double lapse, double shift_x, const double L_inv[3][3])
{
  // The tetrad x-axis basis vector in coord components has ε^0_0 =
  // L_inv[0][0] (= 1/√γ_xx for diagonal γ). A wave at tetrad x-velocity
  // s_tet has coord motion dx^0 = α · L_inv[0][0] · s_tet · dt − β^0 · dt.
  return (lapse * L_inv[0][0] * s_tet) - shift_x;
}

double
gkyl_gr_euler_tetrad_mod_sr_roe_minkowski(double gas_gamma,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[3 * 5], double speeds[3])
{
  // Pure Minkowski SR Roe. Inputs are tetrad-frame conserved variables
  // (Cartesian, no √γ). The metric is η^μν = diag(−1, +1, +1, +1).
  //
  // Recover primitives from each side. For Minkowski, prim_vars reduces
  // to the SR primitive-variable recovery (no √γ, no metric in v² norm).
  // We inline the recovery here rather than going through the auxfields
  // pipeline because q_tet doesn't correspond to any prods row.
  double D_l = ql_tet[0], D_r = qr_tet[0];
  double Sx_l = ql_tet[1], Sx_r = qr_tet[1];
  double Sy_l = ql_tet[2], Sy_r = qr_tet[2];
  double Sz_l = ql_tet[3], Sz_r = qr_tet[3];
  double tau_l = ql_tet[4], tau_r = qr_tet[4];

  // Newton recovery for primitives — same algorithm as
  // gkyl_gr_euler_tetrad_mod_prim_vars, but with √γ=1 and γ_ij=δ_ij.
  double rho_l, vx_l, vy_l, vz_l, p_l, W_l, h_l;
  double rho_r, vx_r, vy_r, vz_r, p_r, W_r, h_r;
  {
    double s_sq_l = ((tau_l + D_l)*(tau_l + D_l)) - (Sx_l*Sx_l + Sy_l*Sy_l + Sz_l*Sz_l);
    double C, C0;
    if (s_sq_l < pow(10.0, -8.0)) {
      C  = D_l / sqrt(pow(10.0, -8.0));
      C0 = (D_l + tau_l) / sqrt(pow(10.0, -8.0));
    } else {
      C  = D_l / sqrt(s_sq_l);
      C0 = (D_l + tau_l) / sqrt(s_sq_l);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_l = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_l = 1.0 / (C * guess);
    rho_l = D_l / W_l;
    vx_l = Sx_l / (rho_l * h_l * W_l*W_l);
    vy_l = Sy_l / (rho_l * h_l * W_l*W_l);
    vz_l = Sz_l / (rho_l * h_l * W_l*W_l);
    p_l  = (rho_l * h_l * W_l*W_l) - D_l - tau_l;
    if (rho_l < pow(10.0, -8.0)) rho_l = pow(10.0, -8.0);
    if (p_l   < pow(10.0, -8.0)) p_l   = pow(10.0, -8.0);
  }
  {
    double s_sq_r = ((tau_r + D_r)*(tau_r + D_r)) - (Sx_r*Sx_r + Sy_r*Sy_r + Sz_r*Sz_r);
    double C, C0;
    if (s_sq_r < pow(10.0, -8.0)) {
      C  = D_r / sqrt(pow(10.0, -8.0));
      C0 = (D_r + tau_r) / sqrt(pow(10.0, -8.0));
    } else {
      C  = D_r / sqrt(s_sq_r);
      C0 = (D_r + tau_r) / sqrt(s_sq_r);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_r = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_r = 1.0 / (C * guess);
    rho_r = D_r / W_r;
    vx_r = Sx_r / (rho_r * h_r * W_r*W_r);
    vy_r = Sy_r / (rho_r * h_r * W_r*W_r);
    vz_r = Sz_r / (rho_r * h_r * W_r*W_r);
    p_r  = (rho_r * h_r * W_r*W_r) - D_r - tau_r;
    if (rho_r < pow(10.0, -8.0)) rho_r = pow(10.0, -8.0);
    if (p_r   < pow(10.0, -8.0)) p_r   = pow(10.0, -8.0);
  }

  // Eulderink-Mellema Roe averages, Minkowski version (γ_ij = δ_ij so
  // D+τ+p = ρhW² and K = √(ρh)).
  double eps_l = p_l / (rho_l * h_l);
  double eps_r = p_r / (rho_r * h_r);
  double K_l = sqrt(D_l + tau_l + p_l) / W_l;
  double K_r = sqrt(D_r + tau_r + p_r) / W_r;
  double K_avg = 1.0 / (K_l + K_r);

  double v0 = ((K_l * W_l)        + (K_r * W_r))        * K_avg;
  double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * K_avg;
  double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * K_avg;
  double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * K_avg;
  double v4 = ((K_l * eps_l)      + (K_r * eps_r))      * K_avg;

  double c_minus = 1.0 - ((gas_gamma / (gas_gamma - 1.0)) * v4);
  double c_plus  = 1.0 + ((gas_gamma / (gas_gamma - 1.0)) * v4);

  double v_alpha_sq = -(v0*v0) + (v1*v1) + (v2*v2) + (v3*v3);
  double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_alpha_sq))
              - (0.5 * (gas_gamma - 1.0) * (1.0 + v_alpha_sq));
  double energy = (v0*v0) - (v1*v1);
  double y = sqrt(((1.0 - (gas_gamma * v4)) * energy) + s_sq);

  // Wave amplitudes (corrected basis with τ-slot rest-mass terms).
  double delta[5];
  for (int i = 0; i < 5; i++) delta[i] = qr_tet[i] - ql_tet[i];
  double sum04 = delta[0] + delta[4];
  double A_sum  = ((v0 * sum04) - (v1 * delta[1])) / energy;
  double B_diff = y * ((v0 * delta[1]) - (v1 * sum04)) / (sqrt(s_sq) * energy);

  double a4 = delta[2] - (v2 * A_sum);
  double a5 = delta[3] - (v3 * A_sum);
  double a3 = ((gas_gamma - 1.0) / s_sq) *
              (delta[0] - (A_sum * c_minus) + (c_plus * ((v2 * a4) + (v3 * a5))));
  double a1 = 0.5 * (A_sum - a3 - B_diff);
  double a2 = 0.5 * (A_sum - a3 + B_diff);

  for (int i = 0; i < 5 * 3; i++) waves_tet[i] = 0.0;

  double *wv;
  wv = &waves_tet[0 * 5];
  wv[0] = a1 * c_minus;
  wv[1] = a1 * (v1 - ((sqrt(s_sq) * v0) / y));
  wv[2] = a1 * v2;
  wv[3] = a1 * v3;
  wv[4] = a1 * (v0 - ((sqrt(s_sq) * v1) / y) - c_minus);
  speeds[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y))
            / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  wv = &waves_tet[1 * 5];
  wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0))))
        - (a4 * c_plus * v2) - (a5 * c_plus * v3);
  wv[1] = a3 * v1;
  wv[2] = (a3 * v2) + a4;
  wv[3] = (a3 * v3) + a5;
  wv[4] = (a3 * (v0 - c_minus - (s_sq / (gas_gamma - 1.0))))
        + (a4 * c_plus * v2) + (a5 * c_plus * v3);
  speeds[1] = v1 / v0;

  wv = &waves_tet[2 * 5];
  wv[0] = a2 * c_minus;
  wv[1] = a2 * (v1 + ((sqrt(s_sq) * v0) / y));
  wv[2] = a2 * v2;
  wv[3] = a2 * v3;
  wv[4] = a2 * (v0 + ((sqrt(s_sq) * v1) / y) - c_minus);
  speeds[2] = (((1.0 - (gas_gamma * v4)) * v0 * v1) + (sqrt(s_sq) * y))
            / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

  double max_s = fmax(fabs(speeds[0]), fmax(fabs(speeds[1]), fabs(speeds[2])));
  return max_s;
}

// Pure Minkowski SR HLL with Davis/Einfeldt wave-speed estimate. Inputs
// are tetrad-frame conserved variables (Cartesian, no √γ). The metric is
// η^μν = diag(−1, +1, +1, +1).
//
// Properties (provable in flat-space SR Euler):
//   - sl ≤ every characteristic ≤ sr (Davis precondition)
//   - q_HLL = (sr·qr − sl·ql + fl − fr)/(sr − sl) is admissible
//     (Mignone-Bodo 2005)
//   - cell update with CFL bound preserves admissibility
//
// Tetrad-first composition: this is the SR core; the curved-frame
// pipeline transforms states/fluxes into the tetrad frame, calls this,
// and back-transforms waves and speeds.
double
gkyl_gr_euler_tetrad_mod_sr_hll_minkowski(double gas_gamma,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[2 * 5], double speeds[2])
{
  double D_l = ql_tet[0], D_r = qr_tet[0];
  double Sx_l = ql_tet[1], Sx_r = qr_tet[1];
  double Sy_l = ql_tet[2], Sy_r = qr_tet[2];
  double Sz_l = ql_tet[3], Sz_r = qr_tet[3];
  double tau_l = ql_tet[4], tau_r = qr_tet[4];

  // Newton recovery of primitives (Banyuls flat). Mirrors the body of
  // gkyl_gr_euler_tetrad_mod_sr_roe_minkowski.
  double rho_l, vx_l, vy_l, vz_l, p_l, W_l, h_l;
  double rho_r, vx_r, vy_r, vz_r, p_r, W_r, h_r;
  {
    double s_sq_l = ((tau_l + D_l)*(tau_l + D_l)) - (Sx_l*Sx_l + Sy_l*Sy_l + Sz_l*Sz_l);
    double C, C0;
    if (s_sq_l < pow(10.0, -10.0)) {
      C  = D_l / sqrt(pow(10.0, -10.0));
      C0 = (D_l + tau_l) / sqrt(pow(10.0, -10.0));
    } else {
      C  = D_l / sqrt(s_sq_l);
      C0 = (D_l + tau_l) / sqrt(s_sq_l);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_l = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_l = 1.0 / (C * guess);
    rho_l = D_l / W_l;
    vx_l = Sx_l / (rho_l * h_l * W_l*W_l);
    vy_l = Sy_l / (rho_l * h_l * W_l*W_l);
    vz_l = Sz_l / (rho_l * h_l * W_l*W_l);
    p_l  = (rho_l * h_l * W_l*W_l) - D_l - tau_l;
    if (rho_l < pow(10.0, -8.0)) rho_l = pow(10.0, -8.0);
    if (p_l   < pow(10.0, -8.0)) p_l   = pow(10.0, -8.0);
  }
  {
    double s_sq_r = ((tau_r + D_r)*(tau_r + D_r)) - (Sx_r*Sx_r + Sy_r*Sy_r + Sz_r*Sz_r);
    double C, C0;
    if (s_sq_r < pow(10.0, -10.0)) {
      C  = D_r / sqrt(pow(10.0, -10.0));
      C0 = (D_r + tau_r) / sqrt(pow(10.0, -10.0));
    } else {
      C  = D_r / sqrt(s_sq_r);
      C0 = (D_r + tau_r) / sqrt(s_sq_r);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_r = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_r = 1.0 / (C * guess);
    rho_r = D_r / W_r;
    vx_r = Sx_r / (rho_r * h_r * W_r*W_r);
    vy_r = Sy_r / (rho_r * h_r * W_r*W_r);
    vz_r = Sz_r / (rho_r * h_r * W_r*W_r);
    p_r  = (rho_r * h_r * W_r*W_r) - D_r - tau_r;
    if (rho_r < pow(10.0, -8.0)) rho_r = pow(10.0, -8.0);
    if (p_r   < pow(10.0, -8.0)) p_r   = pow(10.0, -8.0);
  }

  // Sound speeds (relativistic ideal-gas).
  double c_sl = sqrt((gas_gamma * p_l) / (rho_l * h_l));
  double c_sr = sqrt((gas_gamma * p_r) / (rho_r * h_r));

  // Davis/Einfeldt wave-speed estimate in tetrad frame: per-side
  // relativistic acoustic eigenvalues with velocity addition, then
  // min/max across both sides. Guarantees sl ≤ smin, sr ≥ smax.
  double lambda_minus_l = (vx_l - c_sl) / (1.0 - vx_l * c_sl);
  double lambda_plus_l  = (vx_l + c_sl) / (1.0 + vx_l * c_sl);
  double lambda_minus_r = (vx_r - c_sr) / (1.0 - vx_r * c_sr);
  double lambda_plus_r  = (vx_r + c_sr) / (1.0 + vx_r * c_sr);
  double sl = fmin(lambda_minus_l, lambda_minus_r);
  double sr = fmax(lambda_plus_l, lambda_plus_r);

  // Banyuls fluxes in Minkowski (α=1, √γ=1, β=0):
  //   F[D]    = D · vx
  //   F[S_i]  = S_i · vx + p · δ_i^x  (with covariant S_i = Cartesian S^i in flat)
  //   F[τ]    = (τ + p) · vx
  double fl[5], fr[5];
  fl[0] = D_l  * vx_l;
  fl[1] = Sx_l * vx_l + p_l;
  fl[2] = Sy_l * vx_l;
  fl[3] = Sz_l * vx_l;
  fl[4] = (tau_l + p_l) * vx_l;
  fr[0] = D_r  * vx_r;
  fr[1] = Sx_r * vx_r + p_r;
  fr[2] = Sy_r * vx_r;
  fr[3] = Sz_r * vx_r;
  fr[4] = (tau_r + p_r) * vx_r;

  // HLL intermediate state: q_HLL = (sr·qR − sl·qL + fl − fr)/(sr − sl).
  // For Mignone-Bodo admissible inputs and Davis bracket, q_HLL is itself
  // admissible. Waves are conservative-state jumps from qL → q_HLL → qR.
  double qm[5];
  double denom = sr - sl;
  if (fabs(denom) < 1.0e-14) {
    // Degenerate (qL = qR or both supersonic in the same direction).
    // Fall back to zero waves; qfluct's central form recovers ΔF.
    for (int i = 0; i < 2 * 5; i++) waves_tet[i] = 0.0;
    speeds[0] = sl;
    speeds[1] = sr;
    return fmax(fabs(sl), fabs(sr));
  }
  for (int i = 0; i < 5; i++) {
    qm[i] = (sr * qr_tet[i] - sl * ql_tet[i] + fl[i] - fr[i]) / denom;
  }

  double *w0 = &waves_tet[0 * 5];
  double *w1 = &waves_tet[1 * 5];
  for (int i = 0; i < 5; i++) {
    w0[i] = qm[i] - ql_tet[i];
    w1[i] = qr_tet[i] - qm[i];
  }
  speeds[0] = sl;
  speeds[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

// Pure Minkowski SR Lax-Friedrichs. Symmetric envelope ±amax with
// amax = max over both sides of |λ|_max. Like HLL but with a symmetric
// (broader) speed bracket — more diffusive but still admissibility-
// preserving on admissible inputs.
//
// Tetrad-first composition: same role as sr_hll_minkowski but symmetric
// in speed. The curved-frame wave_lax pipeline transforms states into
// the tetrad frame, calls this, and back-transforms waves/speeds.
double
gkyl_gr_euler_tetrad_mod_sr_lax_minkowski(double gas_gamma,
  const double ql_tet[5], const double qr_tet[5],
  double waves_tet[2 * 5], double speeds[2])
{
  double D_l = ql_tet[0], D_r = qr_tet[0];
  double Sx_l = ql_tet[1], Sx_r = qr_tet[1];
  double Sy_l = ql_tet[2], Sy_r = qr_tet[2];
  double Sz_l = ql_tet[3], Sz_r = qr_tet[3];
  double tau_l = ql_tet[4], tau_r = qr_tet[4];

  // Banyuls Newton recovery (Minkowski flat) on each side.
  double rho_l, vx_l, vy_l, vz_l, p_l, W_l, h_l;
  double rho_r, vx_r, vy_r, vz_r, p_r, W_r, h_r;
  {
    double s_sq_l = ((tau_l + D_l)*(tau_l + D_l)) - (Sx_l*Sx_l + Sy_l*Sy_l + Sz_l*Sz_l);
    double C, C0;
    if (s_sq_l < pow(10.0, -10.0)) {
      C  = D_l / sqrt(pow(10.0, -10.0));
      C0 = (D_l + tau_l) / sqrt(pow(10.0, -10.0));
    } else {
      C  = D_l / sqrt(s_sq_l);
      C0 = (D_l + tau_l) / sqrt(s_sq_l);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_l = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_l = 1.0 / (C * guess);
    rho_l = D_l / W_l;
    vx_l = Sx_l / (rho_l * h_l * W_l*W_l);
    vy_l = Sy_l / (rho_l * h_l * W_l*W_l);
    vz_l = Sz_l / (rho_l * h_l * W_l*W_l);
    p_l  = (rho_l * h_l * W_l*W_l) - D_l - tau_l;
    if (rho_l < pow(10.0, -8.0)) rho_l = pow(10.0, -8.0);
    if (p_l   < pow(10.0, -8.0)) p_l   = pow(10.0, -8.0);
  }
  {
    double s_sq_r = ((tau_r + D_r)*(tau_r + D_r)) - (Sx_r*Sx_r + Sy_r*Sy_r + Sz_r*Sz_r);
    double C, C0;
    if (s_sq_r < pow(10.0, -10.0)) {
      C  = D_r / sqrt(pow(10.0, -10.0));
      C0 = (D_r + tau_r) / sqrt(pow(10.0, -10.0));
    } else {
      C  = D_r / sqrt(s_sq_r);
      C0 = (D_r + tau_r) / sqrt(s_sq_r);
    }
    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0*C0) - 1.0) + 1.0 -
      (C*C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0*C0) - 1.0;
    double eta = 2.0 * C * ((gas_gamma - 1.0) / gas_gamma);
    double guess = 1.0;
    for (int it = 0; it < 100; it++) {
      double poly = (alpha4 * guess*guess*guess) * (guess - eta)
                  + (alpha2 * guess*guess) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess)
                      + (4.0 * alpha4 * guess*guess*guess)
                      - (3.0 * eta * alpha4 * guess*guess);
      double guess_new = guess - poly/poly_der;
      if (fabs(guess - guess_new) < pow(10.0, -14.0)) { guess = guess_new; break; }
      guess = guess_new;
    }
    W_r = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0)/gas_gamma) *
      ((1.0 - C*guess) / (C0*C0 * guess*guess)))));
    h_r = 1.0 / (C * guess);
    rho_r = D_r / W_r;
    vx_r = Sx_r / (rho_r * h_r * W_r*W_r);
    vy_r = Sy_r / (rho_r * h_r * W_r*W_r);
    vz_r = Sz_r / (rho_r * h_r * W_r*W_r);
    p_r  = (rho_r * h_r * W_r*W_r) - D_r - tau_r;
    if (rho_r < pow(10.0, -8.0)) rho_r = pow(10.0, -8.0);
    if (p_r   < pow(10.0, -8.0)) p_r   = pow(10.0, -8.0);
  }

  // Sound speeds and per-side max-abs characteristic speed in tetrad.
  double c_sl = sqrt((gas_gamma * p_l) / (rho_l * h_l));
  double c_sr = sqrt((gas_gamma * p_r) / (rho_r * h_r));

  double lam_minus_l = (vx_l - c_sl) / (1.0 - vx_l * c_sl);
  double lam_plus_l  = (vx_l + c_sl) / (1.0 + vx_l * c_sl);
  double lam_minus_r = (vx_r - c_sr) / (1.0 - vx_r * c_sr);
  double lam_plus_r  = (vx_r + c_sr) / (1.0 + vx_r * c_sr);

  double max_l = fmax(fabs(lam_minus_l), fabs(lam_plus_l));
  double max_r = fmax(fabs(lam_minus_r), fabs(lam_plus_r));
  double amax = fmax(max_l, max_r);

  // Banyuls fluxes in flat tetrad: F[D] = D·vx, F[S_i] = S_i·vx + p·δ_i^x,
  // F[τ] = (τ+p)·vx.
  double fl[5], fr[5];
  fl[0] = D_l  * vx_l;
  fl[1] = Sx_l * vx_l + p_l;
  fl[2] = Sy_l * vx_l;
  fl[3] = Sz_l * vx_l;
  fl[4] = (tau_l + p_l) * vx_l;
  fr[0] = D_r  * vx_r;
  fr[1] = Sx_r * vx_r + p_r;
  fr[2] = Sy_r * vx_r;
  fr[3] = Sz_r * vx_r;
  fr[4] = (tau_r + p_r) * vx_r;

  // Lax wave decomposition: symmetric ±amax envelope.
  //   w_0 = 0.5·(Δq − ΔF/amax)
  //   w_1 = 0.5·(Δq + ΔF/amax)
  // Then Σ s·w = amax·(w_1 − w_0) = ΔF (flux jump exact in tetrad).
  if (!(amax > 0.0)) {
    for (int k = 0; k < 2 * 5; k++) waves_tet[k] = 0.0;
    speeds[0] = -pow(10.0, -8.0);
    speeds[1] =  pow(10.0, -8.0);
    return pow(10.0, -8.0);
  }
  double *w0 = &waves_tet[0 * 5];
  double *w1 = &waves_tet[1 * 5];
  for (int i = 0; i < 5; i++) {
    double dq = qr_tet[i] - ql_tet[i];
    double df = fr[i] - fl[i];
    w0[i] = 0.5 * (dq - df / amax);
    w1[i] = 0.5 * (dq + df / amax);
  }
  speeds[0] = -amax;
  speeds[1] = +amax;

  return amax;
}

// ---------------------------------------------------------------------------
// Riemann-variable conversions and Cartesian-frame rotations.
// ---------------------------------------------------------------------------

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *qin, double *wout)
{
  for (int i = 0; i < 5; i++) wout[i] = qin[i];
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn, const double *qstate,
  const double *win, double *qout)
{
  for (int i = 0; i < 5; i++) qout[i] = win[i];
}

static void
gr_euler_tetrad_mod_wall(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int i = 0; i < 5; i++) ghost[i] = skin[i];
  ghost[1] = -ghost[1];
}

static void
gr_euler_tetrad_mod_no_slip(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qglobal,
  double *GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (grm->auxfields.prods) {
    const int *idx = (grm->rot_call_parity == 0) ? grm->cur_idxl : grm->cur_idxr;
    double *out = (grm->rot_call_parity == 0) ? grm->prodl_local : grm->prodr_local;
    long cidx = gkyl_range_idx(&grm->conf_range, idx);
    const double *in_p = gkyl_array_cfetch(grm->auxfields.prods, cidx);
    rot_spacetime_to_local(tau1, tau2, norm, in_p, out);
  }
  grm->rot_call_parity ^= 1;
}

static inline void
rot_to_global(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qlocal,
  double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];
}

// ---------------------------------------------------------------------------
// Riemann-solver callbacks. Wave-speed estimates and Roe eigenstructure are
// shared with the regular mod variant; the tetrad-specific piece is the
// flat-flux + correction step used inside wave_lax / wave_hll / flux_jump.
// ---------------------------------------------------------------------------

// Tetrad-first Lax-Friedrichs. Mirrors the wave_hll pipeline:
//   1. Build Gram-Schmidt-on-γ⁻¹ triad from the interface metric.
//   2. Forward-transform qL, qR into the (flat) tetrad frame.
//   3. Run pure SR Lax in the tetrad: symmetric ±amax envelope.
//   4. Back-transform waves and speeds to the curved coord frame.
//
// With the Gram-Schmidt triad, flux jump (Σ s·w = ΔF_Banyuls) holds
// exactly in the curved frame and positivity is preserved. The pure
// wave-based qfluct (used by qfluct_lax) then ties conservation and
// flux-jump together.
static double
wave_lax(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  (void)delta;

  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (excise_l || excise_r) {
    for (int k = 0; k < 2 * 5; k++) waves[k] = 0.0;
    s[0] = -pow(10.0, -8.0);
    s[1] =  pow(10.0, -8.0);
    return pow(10.0, -8.0);
  }

  // Step 1: interface metric (arithmetic mean of L and R).
  double g_iface[3][3], inv_g_iface[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_GIJ + 3*i + j]
                           + grm->prodr_local[GKYL_GR_SP_GIJ + 3*i + j]);
      inv_g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_INV_GIJ + 3*i + j]
                               + grm->prodr_local[GKYL_GR_SP_INV_GIJ + 3*i + j]);
    }
  double alpha_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_LAPSE]
                            + grm->prodr_local[GKYL_GR_SP_LAPSE]);
  double shift_x_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_SHIFT + 0]
                              + grm->prodr_local[GKYL_GR_SP_SHIFT + 0]);
  double sqrt_det_l = sqrt(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_r = sqrt(grm->prodr_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_iface = 0.5 * (sqrt_det_l + sqrt_det_r);

  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(
    g_iface, inv_g_iface, M, M_inv);

  // Step 2: forward transform onto the contravariant-x triad.
  double ql_tet[5], qr_tet[5];
  gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(ql, sqrt_det_l, inv_g_iface, M_inv, ql_tet);
  gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(qr, sqrt_det_r, inv_g_iface, M_inv, qr_tet);

  // Step 3: SR Lax in tetrad with symmetric ±amax envelope.
  double waves_tet[2 * 5], speeds_tet[2];
  gkyl_gr_euler_tetrad_mod_sr_lax_minkowski(
    gas_gamma, ql_tet, qr_tet, waves_tet, speeds_tet);

  // Step 4: back-transform waves and speeds.
  double maxs_curved = 0.0;
  for (int k = 0; k < 2; k++) {
    gkyl_gr_euler_tetrad_mod_wave_to_curved_contra(&waves_tet[k * 5],
      sqrt_det_iface, M_inv, &waves[k * 5]);
    s[k] = gkyl_gr_euler_tetrad_mod_speed_to_curved_contra(
      speeds_tet[k], alpha_iface, shift_x_iface, inv_g_iface[0][0]);
    if (fabs(s[k]) > maxs_curved) maxs_curved = fabs(s[k]);
  }
  return maxs_curved;
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_roe(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  // True tetrad-Roe pipeline: compose the four modular helpers declared in
  // gkyl_wv_gr_euler_tetrad_mod_priv.h. Each piece is independently unit-
  // tested (see ctest_wv_gr_euler_tetrad_mod's modular helper tests).
  //
  //   1. Build orthonormal triad from the INTERFACE γ_ij (arithmetic mean
  //      of left and right side spatial metrics — symmetric under L↔R swap).
  //   2. Transform q_GR_L and q_GR_R to tetrad frame (using each side's own
  //      √γ to strip; momentum gets raised by interface γ^{-1} and rotated
  //      by interface L^{-1}). The result is two SR-conserved states in a
  //      common Minkowski frame.
  //   3. Run pure Minkowski SR Roe on the tetrad-frame states. The
  //      Eulderink-Mellema construction satisfies A_SR · ∆q_tet = ∆f_SR
  //      exactly, so wave-sum and flux-jump identities both hold to
  //      floating-point precision in the tetrad frame.
  //   4. Back-transform waves (×√γ + L) and speeds (×α/L[0][0] − β^1) to
  //      the curved coord frame for wave_prop's consumption. The `delta`
  //      parameter is unused — Roe reconstructs its own delta in the tetrad
  //      frame from qr_tet − ql_tet.

  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  // Step 1: interface metric = arithmetic mean of left and right γ_ij;
  // similarly mean of lapse, shift. √γ_interface from average of √γ.
  // (Cholesky-on-γ doesn't need γ^{-1}, so we don't compute it here.)
  double g_iface[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_GIJ + 3*i + j]
                           + grm->prodr_local[GKYL_GR_SP_GIJ + 3*i + j]);
  double alpha_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_LAPSE]
                            + grm->prodr_local[GKYL_GR_SP_LAPSE]);
  double shift_x_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_SHIFT + 0]
                              + grm->prodr_local[GKYL_GR_SP_SHIFT + 0]);
  // Average of the two sides' own √γ rather than √(det of mean γ) — both
  // are O(∆γ)-equivalent and the former matches how each side stripped its
  // √γ when forming q.
  double sqrt_det_l = sqrt(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_r = sqrt(grm->prodr_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_iface = 0.5 * (sqrt_det_l + sqrt_det_r);

  // Cholesky-of-γ construction: L L^T = γ, L lower triangular. Under
  // Convention A the slot stores genuine covariant momentum S_i, so the
  // forward transform projects with L_inv (S_a_tet = L_inv[a][i]·S_i) and
  // the back-transform of waves uses L (w_GR[i+1] = √γ·L[i][a]·w_tet[a+1]).
  double L[3][3], L_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad(g_iface, L, L_inv);

  // Step 2: forward transform each side to the (common) tetrad frame.
  double ql_tet[5], qr_tet[5];
  gkyl_gr_euler_tetrad_mod_q_to_tetrad(ql, sqrt_det_l, L_inv, ql_tet);
  gkyl_gr_euler_tetrad_mod_q_to_tetrad(qr, sqrt_det_r, L_inv, qr_tet);

  // Step 3: pure Minkowski SR Roe. The tetrad-frame max-speed return value
  // is discarded — CFL must use the curved-frame max computed in step 4,
  // because back-transform doesn't commute with abs-max when β^x ≠ 0.
  double waves_tet[3 * 5], speeds_tet[3];
  gkyl_gr_euler_tetrad_mod_sr_roe_minkowski(
    gas_gamma, ql_tet, qr_tet, waves_tet, speeds_tet);

  // Step 4: back-transform waves (using L) and speeds (using L_inv[0][0] =
  // ε^0_0 = 1/√γ_xx for diagonal γ, which is convention-independent for
  // the geometric basis vector).
  double maxs_curved = 0.0;
  for (int k = 0; k < 3; k++) {
    gkyl_gr_euler_tetrad_mod_wave_to_curved(&waves_tet[k * 5],
      sqrt_det_iface, L, &waves[k * 5]);
    s[k] = gkyl_gr_euler_tetrad_mod_speed_to_curved(
      speeds_tet[k], alpha_iface, shift_x_iface, L_inv);
    if (fabs(s[k]) > maxs_curved) maxs_curved = fabs(s[k]);
  }
  return maxs_curved;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  // Central-flux + Roe-dissipation form, parallel to qfluct_lax and qfluct_hll
  // (see wv_gr_euler_tetrad.c::qfluct_roe for the full docstring). amdq +
  // apdq = ΔF by construction even when the tetrad-Roe wave decomposition
  // does not strictly satisfy Σ s·w = ΔF in curved γ.
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql, grm->prodl_local, fl_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr, grm->prodr_local, fr_sr);
  double fl[5], fr[5];
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, ql, grm->prodl_local, fl_sr, fl);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qr, grm->prodr_local, fr_sr, fr);

  const double *w0 = &waves[0 * 5], *w1 = &waves[1 * 5], *w2 = &waves[2 * 5];
  double abs_s0 = fabs(s[0]), abs_s1 = fabs(s[1]), abs_s2 = fabs(s[2]);

  for (int i = 0; i < 5; i++) {
    double df   = fr[i] - fl[i];
    double diss = abs_s0 * w0[i] + abs_s1 * w1[i] + abs_s2 * w2[i];
    amdq[i] = 0.5 * (df - diss);
    apdq[i] = 0.5 * (df + diss);
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return wave_roe(eqn, delta, ql, qr, waves, s);
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// Tetrad-first HLL pipeline (mirrors wave_roe):
//   1. Build orthonormal triad ε^i_a from interface γ_ij via Cholesky.
//   2. Forward-transform qL, qR to the tetrad (Minkowski) frame.
//   3. Run pure-SR HLL with Davis bracket — provably admissible q_HLL.
//   4. Back-transform waves (via L) and speeds (via α·L_inv[0][0] − β^x).
//
// The flat-space SR HLL proof (Mignone-Bodo 2005) guarantees the
// tetrad-frame intermediate state q_HLL_tet is admissible. Because the
// tetrad transform is a valid change of orthonormal basis of the
// covariant momentum, admissibility (D > 0, τ ≥ 0, |S|_γ < D+τ)
// transfers to the curved frame state-by-state. This is the analog of
// wave_roe but exempt from the off-diagonal-γ structural residual,
// because HLL has only two waves and uses the curved-frame ΔF via
// qfluct's central+dissipation form for conservation.
static double
wave_hll(const struct gkyl_wv_eqn *eqn, const double *delta, const double *ql,
  const double *qr, double *waves, double *s)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  // Excision short-circuit (mirrors wave_roe).
  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (excise_l || excise_r) {
    for (int k = 0; k < 2 * 5; k++) waves[k] = 0.0;
    s[0] = 0.0;
    s[1] = 0.0;
    return pow(10.0, -8.0);
  }

  // Interface metric (arithmetic mean of L and R γ_ij, γ^{ij}).
  double g_iface[3][3], inv_g_iface[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_GIJ + 3*i + j]
                           + grm->prodr_local[GKYL_GR_SP_GIJ + 3*i + j]);
      inv_g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_INV_GIJ + 3*i + j]
                               + grm->prodr_local[GKYL_GR_SP_INV_GIJ + 3*i + j]);
    }
  double alpha_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_LAPSE]
                            + grm->prodr_local[GKYL_GR_SP_LAPSE]);
  double shift_x_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_SHIFT + 0]
                              + grm->prodr_local[GKYL_GR_SP_SHIFT + 0]);
  double sqrt_det_l = sqrt(grm->prodl_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_r = sqrt(grm->prodr_local[GKYL_GR_SP_SPATIAL_DET]);
  double sqrt_det_iface = 0.5 * (sqrt_det_l + sqrt_det_r);

  // Step 1: Gram-Schmidt-on-γ⁻¹ triad. First basis vector is the
  // contravariant x-direction (γ^{xi}/√γ^{xx}); the other two come from
  // Gram-Schmidt-in-γ. This eliminates the v_tet^x ↔ v^y, v^z mixing
  // seen with Cholesky for non-diagonal γ.
  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_mod_build_triad_contravariant_x(
    g_iface, inv_g_iface, M, M_inv);

  // Step 2: forward transform Convention-A covariant momentum onto the
  // contravariant-x triad. v_tet^x = v^x/√γ^{xx} (clean — no mixing).
  double ql_tet[5], qr_tet[5];
  gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(ql, sqrt_det_l, inv_g_iface, M_inv, ql_tet);
  gkyl_gr_euler_tetrad_mod_q_to_tetrad_contra(qr, sqrt_det_r, inv_g_iface, M_inv, qr_tet);

  // Step 3: pure Minkowski SR HLL with Davis bracket.
  double waves_tet[2 * 5], speeds_tet[2];
  gkyl_gr_euler_tetrad_mod_sr_hll_minkowski(
    gas_gamma, ql_tet, qr_tet, waves_tet, speeds_tet);

  // Step 4: back-transform waves and speeds with the new triad. The
  // coord-x momentum slot receives only the a=0 tetrad wave (mirror of
  // forward-clean property).
  double maxs_curved = 0.0;
  for (int k = 0; k < 2; k++) {
    gkyl_gr_euler_tetrad_mod_wave_to_curved_contra(&waves_tet[k * 5],
      sqrt_det_iface, M_inv, &waves[k * 5]);
    s[k] = gkyl_gr_euler_tetrad_mod_speed_to_curved_contra(
      speeds_tet[k], alpha_iface, shift_x_iface, inv_g_iface[0][0]);
    if (fabs(s[k]) > maxs_curved) maxs_curved = fabs(s[k]);
  }

  return maxs_curved;
}

// Pure wave-based qfluct: amdq = Σ min(0, s_k)·w_k, apdq = Σ max(0, s_k)·w_k.
// Conservation `amdq + apdq = Σ s·w` is identical to the flux-jump
// identity by construction. With the Gram-Schmidt-on-γ⁻¹ triad, flux
// jump `Σ s_curved · w_curved = ΔF_Banyuls` holds exactly in curved
// spacetime, so this form simultaneously gives strict conservation AND
// the convex-combination property needed for positivity preservation.
static void
qfluct_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  (void)eqn; (void)ql; (void)qr;
  const double *w0 = &waves[0 * 5], *w1 = &waves[1 * 5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = s0m * w0[i] + s1m * w1[i];
    apdq[i] = s0p * w0[i] + s1p * w1[i];
  }
}

// POSITIVITY_SWEEP fallback policy: HIGH_ORDER uses HLL (the production
// flux); LOW_ORDER (the second-pass sweep on cells flagged inadmissible)
// falls back to Lax. Even though tetrad-first HLL is admissibility-
// preserving in isolation, observation shows that HLL on the BHL
// bow-shock leaves substantial momentum overshoot (s² close to (D+τ)²)
// that propagates pathological local solution quality even when the
// admissibility constraints are nominally satisfied. Lax is more
// diffusive, smears these states, and leaves a much cleaner solution
// at the cost of additional viscosity. The mixed HIGH-then-LOW pattern
// is fine here: POSITIVITY_SWEEP resets qout to qin and recomputes
// fluctuations at flagged interfaces only; the unflagged interfaces
// retain their FIRST_SWEEP HLL values. The mixed cell-update is no
// worse than a CFL-constrained Lax pass for the flagged cells, and
// HLL elsewhere keeps the bulk solution sharp.
static double
wave_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return wave_hll(eqn, delta, ql, qr, waves, s);
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump_func(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, double *flux_jump)
{
  // Mirrors the regular mod variant's flux_jump_func but uses the corrected
  // (flat + GR-correction) flux to match packed wv_gr_euler_tetrad.c:1581.
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  double gas_gamma = grm->gas_gamma;

  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, ql, grm->prodl_local, fl_sr);
  gkyl_gr_euler_tetrad_mod_flux(gas_gamma, qr, grm->prodr_local, fr_sr);

  double fl[5], fr[5];
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, ql, grm->prodl_local, fl_sr, fl);
  gkyl_gr_euler_tetrad_mod_flux_correction(gas_gamma, qr, grm->prodr_local, fr_sr, fr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (!excise_l && !excise_r) {
    for (int m = 0; m < 5; m++) flux_jump[m] = fr[m] - fl[m];
  } else {
    for (int m = 0; m < 5; m++) flux_jump[m] = 0.0;
  }

  double amaxl = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, ql, grm->prodl_local);
  double amaxr = gkyl_gr_euler_tetrad_mod_max_abs_speed(gas_gamma, qr, grm->prodr_local);
  return fmax(amaxl, amaxr);
}

// Strict admissibility predicate (Convention A): D > 0, τ ≥ 0, and
// s² = (D+τ)² − γ^{ij}·S_i·S_j > 0. Returns true iff the cell lies in
// the convex set the Banyuls inversion can reach without floors firing.
//
// The previous test relied on `v[0] < 0 || v[4] < 0` AFTER floors had
// already pushed those positive — it never fired. Routing through
// gkyl_gr_euler_check_admissibility (which inspects the pre-floor
// conservatives directly) is what makes the test meaningful.
static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (!grm->auxfields.prods) return q[0] >= 0.0;

  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);

  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) return true;

  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];
  if (!(spatial_det > 0.0)) return false;

  const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
  double inv_g[3][3] = {
    { ig[0], ig[1], ig[2] },
    { ig[3], ig[4], ig[5] },
    { ig[6], ig[7], ig[8] },
  };
  double sd  = sqrt(spatial_det);
  double D   = q[0] / sd;
  double Sx  = q[1] / sd;
  double Sy  = q[2] / sd;
  double Sz  = q[3] / sd;
  double tau = q[4] / sd;

  return gkyl_gr_euler_check_admissibility(D, Sx, Sy, Sz, tau, inv_g)
    == GR_EULER_ADM_OK;
}

// One-iteration projection onto the convex admissibility set. Caller
// loops (check_inv → repair_state → check_inv) up to a small bounded
// number of passes; each call fixes the first failing constraint
// (D, then τ, then s²).
static void
repair_state(const struct gkyl_wv_eqn *eqn, double *q)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (!grm->auxfields.prods) return;

  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) return;

  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];
  if (!(spatial_det > 0.0)) return;

  const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
  double inv_g[3][3] = {
    { ig[0], ig[1], ig[2] },
    { ig[3], ig[4], ig[5] },
    { ig[6], ig[7], ig[8] },
  };
  double sd  = sqrt(spatial_det);
  double D   = q[0] / sd;
  double Sx  = q[1] / sd;
  double Sy  = q[2] / sd;
  double Sz  = q[3] / sd;
  double tau = q[4] / sd;

  unsigned int fixed =
    gkyl_gr_euler_repair_admissibility_cascade(inv_g, &D, &Sx, &Sy, &Sz, &tau);

  // Per-constraint independent tallies, split by call-site context. A
  // non-zero wave_prop count means the post-positivity-sweep low-order
  // (Lax/HLL) flux failed to produce an admissible state — which our
  // unit tests say should never happen. Source counts dominate when
  // the SSP-RK3 forward-Euler pushes already-admissible cells across
  // an admissibility boundary inside a stage.
  uint64_t *cnt = (eqn->cur_repair_ctx == 1)
    ? grm->repair_count_wave_prop
    : grm->repair_count_source;
  if (fixed & GR_EULER_REPAIR_D)   cnt[GR_EULER_ADM_BAD_D]   += 1;
  if (fixed & GR_EULER_REPAIR_TAU) cnt[GR_EULER_ADM_BAD_TAU] += 1;
  if (fixed & GR_EULER_REPAIR_S2)  cnt[GR_EULER_ADM_BAD_S2]  += 1;

  q[0] = D   * sd;
  q[1] = Sx  * sd;
  q[2] = Sy  * sd;
  q[3] = Sz  * sd;
  q[4] = tau * sd;
}

static double
max_speed_func(const struct gkyl_wv_eqn *eqn, const double *q)
{
  struct wv_gr_euler_tetrad_mod *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  if (!grm->auxfields.prods) return 1.0;
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);
  return gkyl_gr_euler_tetrad_mod_max_abs_speed(grm->gas_gamma, q, prods);
}

static inline void
gr_euler_tetrad_mod_cons_to_diag(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *diag)
{
  for (int i = 0; i < 5; i++) diag[i] = qin[i];
}

static inline void
gr_euler_tetrad_mod_source(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *sout)
{
  // Integrated by moment_spacetime_coupling.
  for (int i = 0; i < 5; i++) sout[i] = 0.0;
}

// ---------------------------------------------------------------------------
// Constructors and lifetime.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_mod_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_gr_euler_tetrad_mod *grm = container_of(base,
    struct wv_gr_euler_tetrad_mod, eqn);

  // Per-constraint independent counts, split by call site. A non-zero
  // wave_prop count would mean the low-order POSITIVITY_SWEEP flux did
  // not preserve admissibility — the unit tests assert this never
  // happens for Lax/HLL on admissible inputs, so any nonzero here is a
  // signal that the positivity-sweep mechanism itself is broken.
  uint64_t total_source =
      grm->repair_count_source[GR_EULER_ADM_BAD_D]
    + grm->repair_count_source[GR_EULER_ADM_BAD_TAU]
    + grm->repair_count_source[GR_EULER_ADM_BAD_S2];
  uint64_t total_wave_prop =
      grm->repair_count_wave_prop[GR_EULER_ADM_BAD_D]
    + grm->repair_count_wave_prop[GR_EULER_ADM_BAD_TAU]
    + grm->repair_count_wave_prop[GR_EULER_ADM_BAD_S2];
  if (total_source + total_wave_prop > 0) {
    fprintf(stderr,
      "[gr_euler_tetrad_mod] cascade-repair fix totals: %llu fixes "
      "(source: %llu, wave_prop: %llu)\n",
      (unsigned long long)(total_source + total_wave_prop),
      (unsigned long long)total_source,
      (unsigned long long)total_wave_prop);
    fprintf(stderr,
      "  source     — D<=0: %llu, tau<0: %llu, s^2<=0: %llu\n",
      (unsigned long long)grm->repair_count_source[GR_EULER_ADM_BAD_D],
      (unsigned long long)grm->repair_count_source[GR_EULER_ADM_BAD_TAU],
      (unsigned long long)grm->repair_count_source[GR_EULER_ADM_BAD_S2]);
    fprintf(stderr,
      "  wave_prop  — D<=0: %llu, tau<0: %llu, s^2<=0: %llu\n",
      (unsigned long long)grm->repair_count_wave_prop[GR_EULER_ADM_BAD_D],
      (unsigned long long)grm->repair_count_wave_prop[GR_EULER_ADM_BAD_TAU],
      (unsigned long long)grm->repair_count_wave_prop[GR_EULER_ADM_BAD_S2]);
  }

  gkyl_free(grm);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_new(double gas_gamma,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  return gkyl_wv_gr_euler_tetrad_mod_inew(
    &(struct gkyl_wv_gr_euler_tetrad_mod_inp) {
      .gas_gamma = gas_gamma,
      .conf_range = *conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_HLL,
      .use_gpu = use_gpu,
    });
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_inew(
  const struct gkyl_wv_gr_euler_tetrad_mod_inp *inp)
{
  struct wv_gr_euler_tetrad_mod *grm =
    gkyl_malloc(sizeof(struct wv_gr_euler_tetrad_mod));

  grm->eqn.type = GKYL_EQN_GR_EULER_TETRAD_MOD;
  grm->eqn.num_equations = 5;
  grm->eqn.num_diag = 5;

  grm->gas_gamma = inp->gas_gamma;
  grm->conf_range = inp->conf_range;
  grm->auxfields.prods = NULL;
  grm->rot_call_parity = 0;
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = 0;
    grm->cur_idxr[d] = 0;
    grm->cur_cell_idx[d] = 0;
  }
  for (int k = 0; k < 4; k++) {
    grm->repair_count_source[k]    = 0;
    grm->repair_count_wave_prop[k] = 0;
  }
  grm->eqn.cur_repair_ctx = 0;

  if (inp->rp_type == WV_GR_EULER_TETRAD_RP_LAX) {
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_lax_l;
    grm->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_GR_EULER_TETRAD_RP_ROE) {
    grm->eqn.num_waves = 3;
    grm->eqn.waves_func = wave_roe_l;
    grm->eqn.qfluct_func = qfluct_roe_l;
  }
  else {  // default: HLL
    grm->eqn.num_waves = 2;
    grm->eqn.waves_func = wave_hll_l;
    grm->eqn.qfluct_func = qfluct_hll_l;
  }

  grm->eqn.ffluct_func = NULL;
  grm->eqn.flux_jump = flux_jump_func;
  grm->eqn.check_inv_func = check_inv;
  grm->eqn.max_speed_func = max_speed_func;
  grm->eqn.rotate_to_local_func = rot_to_local;
  grm->eqn.rotate_to_global_func = rot_to_global;
  grm->eqn.wall_bc_func = gr_euler_tetrad_mod_wall;
  grm->eqn.no_slip_bc_func = gr_euler_tetrad_mod_no_slip;
  grm->eqn.cons_to_riem = cons_to_riem;
  grm->eqn.riem_to_cons = riem_to_cons;
  grm->eqn.cons_to_diag = gr_euler_tetrad_mod_cons_to_diag;
  grm->eqn.source_func = gr_euler_tetrad_mod_source;

  grm->eqn.set_interface_idx_func = gr_euler_tetrad_mod_set_interface_idx;
  grm->eqn.set_cell_idx_func = gr_euler_tetrad_mod_set_cell_idx;
  grm->eqn.repair_state_func = repair_state;

  grm->eqn.embed_geo = NULL;
  grm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(grm->eqn.flags);
  grm->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_tetrad_mod_free);
  grm->eqn.on_dev = &grm->eqn;

  return &grm->eqn;
}

double
gkyl_wv_gr_euler_tetrad_mod_gas_gamma(const struct gkyl_wv_eqn *eqn)
{
  const struct wv_gr_euler_tetrad_mod *grm = container_of(eqn,
    struct wv_gr_euler_tetrad_mod, eqn);
  return grm->gas_gamma;
}
