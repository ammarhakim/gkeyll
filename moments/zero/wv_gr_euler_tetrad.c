#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_wv_gr_euler_prim_priv.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

// ---------------------------------------------------------------------------
// Setters for cell-index communication from wave_prop. Identical to the
// regular mod variant; the per-instance idx slots populated below are read
// by downstream callbacks to fetch the appropriate row from the spacetime
// products array.
// ---------------------------------------------------------------------------

static void
gr_euler_tetrad_set_interface_idx(const struct gkyl_wv_eqn *eqn,
  const int *idxl, const int *idxr)
{
  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = idxl[d];
    grm->cur_idxr[d] = idxr[d];
  }
  grm->rot_call_parity = 0;
}

static void
gr_euler_tetrad_set_cell_idx(const struct gkyl_wv_eqn *eqn, const int *idx)
{
  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_cell_idx[d] = idx[d];
  }
}

void
gkyl_gr_euler_tetrad_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_tetrad_auxfields auxin)
{
  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
  grm->auxfields.prods                   = auxin.prods;
  grm->auxfields.prim_status_wave_prop   = auxin.prim_status_wave_prop;
  grm->auxfields.repair_status_wave_prop = auxin.repair_status_wave_prop;
  grm->auxfields.repair_status_source    = auxin.repair_status_source;
}

void
gkyl_gr_euler_tetrad_set_conf_range(const struct gkyl_wv_eqn *eqn,
  const struct gkyl_range *conf_range)
{
  struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
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
//   gkyl_gr_euler_tetrad_flux:            flat-space SR flux (Cartesian W,
//                                             no lapse/shift/√γ)
//   gkyl_gr_euler_tetrad_flux_correction: applies α·√γ and replaces
//                                             (vx, W_flat) with (vx - βˣ/α,
//                                             W_curved) to recover the
//                                             curved-space flux.
// Mirrors the packed factorization in wv_gr_euler_tetrad.c at lines 10-124.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_prim_vars(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat, double v[5])
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
    // source of truth for the Banyuls Newton solve. The EOS bundle
    // controls the closing equation (IDEAL or APPROXIMATE_SYNGE) in
    // the helper's dispatch.
    const double *ig = &prods[GKYL_GR_SP_INV_GIJ];
    double inv_g[3][3] = {
      { ig[0], ig[1], ig[2] },
      { ig[3], ig[4], ig[5] },
      { ig[6], ig[7], ig[8] },
    };
    struct gkyl_gr_euler_prim prim;
    gkyl_gr_euler_recover_primitives(eos,
      D, momx, momy, momz, Etot, inv_g, stat, &prim);

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
gkyl_gr_euler_tetrad_flux(struct gkyl_gr_euler_eos eos, const double q[5],
  const double *prods,
  struct gkyl_gr_euler_prim_status *stat, double flux_sr[5])
{
  double v[5];
  gkyl_gr_euler_tetrad_prim_vars(eos, q, prods, stat, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  bool in_excision_region = false;
  if (prods[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double v_dot = (vx*vx) + (vy*vy) + (vz*vz);
    double W = 1.0 / sqrt(1.0 - v_dot);
    if (v_dot > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(pow(10.0, -8.0));

    // EOS dispatch via the shared helper (IDEAL γ-law or
    // APPROXIMATE_SYNGE with the use_rcc closure selector).
    double h = gkyl_gr_euler_eos_enthalpy(eos, rho, p);

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
gkyl_gr_euler_tetrad_flux_correction(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat,
  const double flux_sr[5], double flux_gr[5])
{
  // Mirrors gkyl_gr_euler_tetrad_flux_correction in the packed implementation:
  // scales the flat SR flux by the W_curved / W_flat ratios that arise when
  // mapping the tetrad-frame flux back into the coord frame, plus the α·√γ
  // densitization and shift correction. We deliberately keep this byte-for-
  // byte equivalent to packed so the *only* algorithmic difference between
  // packed and mod-tetrad lives in the Roe solve. See wv_gr_euler_tetrad.c
  // for the parallel implementation.
  double v[5];
  gkyl_gr_euler_tetrad_prim_vars(eos, q, prods, stat, v);
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
gkyl_gr_euler_tetrad_max_abs_speed(struct gkyl_gr_euler_eos eos,
  const double q[5], const double *prods,
  struct gkyl_gr_euler_prim_status *stat)
{
  double v[5];
  gkyl_gr_euler_tetrad_prim_vars(eos, q, prods, stat, v);
  double rho = v[0], vx = v[1], vy = v[2], vz = v[3], p = v[4];

  double lapse   = prods[GKYL_GR_SP_LAPSE];
  double shift_x = prods[GKYL_GR_SP_SHIFT + 0];
  double shift_y = prods[GKYL_GR_SP_SHIFT + 1];
  double shift_z = prods[GKYL_GR_SP_SHIFT + 2];

  // Sound speed via EOS dispatch (see gkyl_gr_euler_eos_cs2 for the
  // per-EOS algebra).
  double h = gkyl_gr_euler_eos_enthalpy(eos, rho, p);
  double cs2 = gkyl_gr_euler_eos_cs2(eos, rho, p, h);
  if (cs2 < 0.0) cs2 = 0.0;
  double c_s = sqrt(cs2);

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
      // Clamp the radical to ≥ 0 to avoid NaN-poisoning when γ^{ii}·v²
      // approaches the relativistic limit. NaN > x is always false, so an
      // unclamped sqrt(rad<0) would silently drop the direction's
      // contribution to max_eig and under-bound the Lax envelope.
      if (rad < 0.0) rad = 0.0;
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
gkyl_gr_euler_tetrad_build_triad_contravariant_x(
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
gkyl_gr_euler_tetrad_q_to_tetrad_contra(const double q_GR[5],
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
gkyl_gr_euler_tetrad_wave_to_curved_contra(const double w_tet[5],
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
gkyl_gr_euler_tetrad_speed_to_curved_contra(double s_tet,
  double lapse, double shift_x, double inv_gxx)
{
  return (lapse * sqrt(inv_gxx) * s_tet) - shift_x;
}

// Detect the all-zero "vacuum" tetrad state used by the excision-boundary
// short-circuit in wave_lax / wave_hll / wave_hllc. The excised side is
// fed in as q = [0,0,0,0,0] (in coord frame) and back-transforms to
// [0,0,0,0,0] in the tetrad frame. The Banyuls Newton primitive recovery
// in sr_*_minkowski would divide by D = 0 and produce NaN / saturated
// floors. Detecting this state lets each SR core skip Newton on that
// side and substitute exact vacuum primitives (ρ=0, v=0, p=0, W=1, h=1).
// SR fluxes for a vacuum state are then identically zero, so the wave
// decomposition correctly reduces to the absorbing-BC case Σ s·w =
// F_R − F_L = ∓F(q_active). See SESSION_NOTES_3.md §11.10–§11.13.
//
// Tight tolerance: in production wave_prop the excised-side state is
// initialized to exactly zero and stays at zero (the source step zeros
// S_rate, cascade-repair never touches excised cells). So an exactly-
// equal-to-zero check would suffice; the small ε > 0 below is purely
// defensive against floating-point noise.
static inline bool
gkyl_gr_euler_tetrad_is_zero_state(const double q[5])
{
  const double eps = 1.0e-14;
  return fabs(q[0]) < eps && fabs(q[1]) < eps && fabs(q[2]) < eps
      && fabs(q[3]) < eps && fabs(q[4]) < eps;
}

double
gkyl_gr_euler_tetrad_sr_roe_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[3 * 5], double speeds[3])
{
  // Signature matches the unified gkyl_sr_riemann_kernel_t typedef. Roe
  // is IDEAL-only; the constructor rejects non-IDEAL eos.type, so we
  // can safely extract gas_gamma. The stat pointer is reserved for
  // future per-callsite instrumentation symmetry; Roe's inlined
  // primitive recovery does not yet write into it.
  double gas_gamma = eos.gas_gamma;

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
  // gkyl_gr_euler_tetrad_prim_vars, but with √γ=1 and γ_ij=δ_ij.
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
gkyl_gr_euler_tetrad_sr_hll_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[2 * 5], double speeds[2])
{
  double D_l = ql_tet[0], D_r = qr_tet[0];
  double Sx_l = ql_tet[1], Sx_r = qr_tet[1];
  double Sy_l = ql_tet[2], Sy_r = qr_tet[2];
  double Sz_l = ql_tet[3], Sz_r = qr_tet[3];
  double tau_l = ql_tet[4], tau_r = qr_tet[4];

  // Banyuls primitive recovery via the shared helper (IDEAL → Eulderink-
  // Mellema quartic Newton; APPROXIMATE_SYNGE → TM cubic + optional RC
  // Newton refinement). The Minkowski tetrad frame has γ_ij = δ_ij, so
  // inv_g = I and the shared helper's curved-frame contraction collapses
  // to the Cartesian dot product the original inline Newton used.
  //
  // Each side short-circuits to vacuum primitives when fed the all-zero
  // excision state (gkyl_gr_euler_tetrad_is_zero_state). SR fluxes
  // then identically zero on that side and the Davis bracket is bounded
  // by the active side's wave speeds.
  double inv_g_flat[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };
  double rho_l, vx_l, vy_l, vz_l, p_l, W_l, h_l;
  double rho_r, vx_r, vy_r, vz_r, p_r, W_r, h_r;
  if (gkyl_gr_euler_tetrad_is_zero_state(ql_tet)) {
    rho_l = 1.0e-30; vx_l = vy_l = vz_l = 0.0; p_l = 0.0; W_l = 1.0; h_l = 1.0;
  } else {
    struct gkyl_gr_euler_prim pl;
    gkyl_gr_euler_recover_primitives(eos, D_l, Sx_l, Sy_l, Sz_l, tau_l, inv_g_flat, stat, &pl);
    rho_l = pl.rho; vx_l = pl.v[0]; vy_l = pl.v[1]; vz_l = pl.v[2];
    p_l   = pl.p;   W_l  = pl.W;    h_l  = pl.h;
  }
  if (gkyl_gr_euler_tetrad_is_zero_state(qr_tet)) {
    rho_r = 1.0e-30; vx_r = vy_r = vz_r = 0.0; p_r = 0.0; W_r = 1.0; h_r = 1.0;
  } else {
    struct gkyl_gr_euler_prim pr;
    gkyl_gr_euler_recover_primitives(eos, D_r, Sx_r, Sy_r, Sz_r, tau_r, inv_g_flat, stat, &pr);
    rho_r = pr.rho; vx_r = pr.v[0]; vy_r = pr.v[1]; vz_r = pr.v[2];
    p_r   = pr.p;   W_r  = pr.W;    h_r  = pr.h;
  }

  // Sound speeds via EOS dispatch. Vacuum-state side gets c_s = 0 (p = 0
  // ⇒ c_s² = 0 for any reasonable EOS).
  double cs2_l = (p_l > 0.0) ? gkyl_gr_euler_eos_cs2(eos, rho_l, p_l, h_l) : 0.0;
  double cs2_r = (p_r > 0.0) ? gkyl_gr_euler_eos_cs2(eos, rho_r, p_r, h_r) : 0.0;
  if (cs2_l < 0.0) cs2_l = 0.0;
  if (cs2_r < 0.0) cs2_r = 0.0;
  double c_sl = sqrt(cs2_l);
  double c_sr = sqrt(cs2_r);

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
gkyl_gr_euler_tetrad_sr_lax_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[2 * 5], double speeds[2])
{
  double D_l = ql_tet[0], D_r = qr_tet[0];
  double Sx_l = ql_tet[1], Sx_r = qr_tet[1];
  double Sy_l = ql_tet[2], Sy_r = qr_tet[2];
  double Sz_l = ql_tet[3], Sz_r = qr_tet[3];
  double tau_l = ql_tet[4], tau_r = qr_tet[4];

  // Banyuls primitive recovery via the shared helper (eos dispatch). See
  // sr_hll_minkowski above for the rationale on identity inv_g and the
  // zero-state short-circuit.
  double inv_g_flat[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };
  double rho_l, vx_l, vy_l, vz_l, p_l, W_l, h_l;
  double rho_r, vx_r, vy_r, vz_r, p_r, W_r, h_r;
  if (gkyl_gr_euler_tetrad_is_zero_state(ql_tet)) {
    rho_l = 1.0e-30; vx_l = vy_l = vz_l = 0.0; p_l = 0.0; W_l = 1.0; h_l = 1.0;
  } else {
    struct gkyl_gr_euler_prim pl;
    gkyl_gr_euler_recover_primitives(eos, D_l, Sx_l, Sy_l, Sz_l, tau_l, inv_g_flat, stat, &pl);
    rho_l = pl.rho; vx_l = pl.v[0]; vy_l = pl.v[1]; vz_l = pl.v[2];
    p_l   = pl.p;   W_l  = pl.W;    h_l  = pl.h;
  }
  if (gkyl_gr_euler_tetrad_is_zero_state(qr_tet)) {
    rho_r = 1.0e-30; vx_r = vy_r = vz_r = 0.0; p_r = 0.0; W_r = 1.0; h_r = 1.0;
  } else {
    struct gkyl_gr_euler_prim pr;
    gkyl_gr_euler_recover_primitives(eos, D_r, Sx_r, Sy_r, Sz_r, tau_r, inv_g_flat, stat, &pr);
    rho_r = pr.rho; vx_r = pr.v[0]; vy_r = pr.v[1]; vz_r = pr.v[2];
    p_r   = pr.p;   W_r  = pr.W;    h_r  = pr.h;
  }

  // Sound speeds via EOS dispatch. Vacuum-state side gets c_s = 0.
  double cs2_l = (p_l > 0.0) ? gkyl_gr_euler_eos_cs2(eos, rho_l, p_l, h_l) : 0.0;
  double cs2_r = (p_r > 0.0) ? gkyl_gr_euler_eos_cs2(eos, rho_r, p_r, h_r) : 0.0;
  if (cs2_l < 0.0) cs2_l = 0.0;
  if (cs2_r < 0.0) cs2_r = 0.0;
  double c_sl = sqrt(cs2_l);
  double c_sr = sqrt(cs2_r);

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
// SR HLLC Riemann solver in Minkowski (tetrad-frame).
// Reference: Mignone & Bodo 2005, MNRAS 364, 126.
//
// Three-wave structure: λ_L (left acoustic), λ* (contact), λ_R (right
// acoustic). Star states give explicit τ-positivity from admissible
// inputs (see MB05 §3.1.2 and Appendix A) — this is the property HLL
// fails to deliver in the τ/D ≪ 1 regime (SESSION_NOTES_2.md §17).
//
// Conventions: this routine uses Banyuls (D, S_i, τ) conservatives.
// Mignone-Bodo writes E := τ + D; we convert internally.
//
// Waves are conservative-state jumps so that Σ s·w = ΔF in the tetrad
// frame:
//   w_L = U_L*  − U_L      (speed λ_L)
//   w_M = U_R* − U_L*      (speed λ*)
//   w_R = U_R   − U_R*     (speed λ_R)
//
// Production-grade fixes carried over from
// einstein/testing/gkeyll/moments/zero/wv_gr_euler_tetrad.c
// (wv_gr_euler_tetrad_sr_hllc_minkowski, lines 644-908):
//
//   1. Stable quadratic for λ* (citardauq form when −b > 0). The naive
//      (−b − √D)/(2a) suffers catastrophic cancellation when −b ≈ √D
//      (e.g. near-uniform states); the alternate 2c/(−b + √D) is well-
//      behaved. See line 784 in testing/gkeyll for the original note.
//
//   2. Sign-corrected p* formula:
//        p* = (A·λ* − B) / (1 − λ·λ*)
//      MB05 eq (17) as printed has (1 + λ·λ*) in the denominator; this
//      is a sign error. The corrected form satisfies p*_L = p*_R = p
//      in the trivial Riemann limit (ql = qr ⇒ U_L* = U_R* = U
//      exactly). See line 847 in testing/gkeyll for the derivation.
//
//   3. Cold-gas / degenerate-fan fallback to HLL: when c_s ≈ 0 the
//      acoustic waves λ_L, λ_R collapse and the λ* quadratic returns
//      either a root outside [λ_L, λ_R] or one arbitrarily close to
//      an outer wave (so 1/(λ_R − λ*) blows up the star state).
//      Detect either pathology and gracefully degrade to HLL: collapse
//      the contact wave to zero and let the two outer waves carry the
//      HLL-averaged jump. See line 810 in testing/gkeyll. CAVEAT: HLL
//      itself is NOT τ-positivity-preserving in this regime (see
//      SESSION_NOTES_2.md §17 and ctest test_small_tau_over_D_hll).
//      The fallback degrades gracefully but inherits HLL's known
//      weakness. Hopefully reaching the fallback is rare on actual
//      problems.
double
gkyl_gr_euler_tetrad_sr_hllc_minkowski(struct gkyl_gr_euler_eos eos,
  const double ql_tet[5], const double qr_tet[5],
  struct gkyl_gr_euler_prim_status *stat,
  double waves_tet[3 * 5], double speeds[3])
{
  // HLLC star-state diagnostics live on stat->hllc (if stat != NULL).
  // Per-call last_* fields overwrite on each call; aggregate counters
  // accumulate. Reset last_* at function entry so they reflect the
  // current Riemann problem only.
  if (stat) {
    stat->hllc.last_did_fallback = 0;
    stat->hllc.last_fallback_reason = 0;
    stat->hllc.last_lambda_L = 0.0;
    stat->hllc.last_lambda_R = 0.0;
    stat->hllc.last_lambda_star = 0.0 / 0.0;  // NaN until computed
  }
  // Pure Minkowski primitive recovery — γ_ij = δ_ij so the inv_g passed
  // to the shared Banyuls helper is identity. Convention A's covariant-
  // momentum bookkeeping is degenerate to contravariant in flat space.
  double D_l = ql_tet[0], Sx_l = ql_tet[1], Sy_l = ql_tet[2], Sz_l = ql_tet[3], tau_l = ql_tet[4];
  double D_r = qr_tet[0], Sx_r = qr_tet[1], Sy_r = qr_tet[2], Sz_r = qr_tet[3], tau_r = qr_tet[4];

  double inv_g_flat[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  struct gkyl_gr_euler_prim pl, pr;
  gkyl_gr_euler_recover_primitives(eos, D_l, Sx_l, Sy_l, Sz_l, tau_l, inv_g_flat, stat, &pl);
  gkyl_gr_euler_recover_primitives(eos, D_r, Sx_r, Sy_r, Sz_r, tau_r, inv_g_flat, stat, &pr);

  double rho_l = pl.rho, vx_l = pl.v[0], vy_l = pl.v[1], vz_l = pl.v[2], p_l = pl.p, W_l = pl.W, h_l = pl.h;
  double rho_r = pr.rho, vx_r = pr.v[0], vy_r = pr.v[1], vz_r = pr.v[2], p_r = pr.p, W_r = pr.W, h_r = pr.h;

  // Davis (1988) wave-speed estimate, eq (21–23) of MB05.
  // σ_s = c_s² / (γ²(1 − c_s²)), λ_± = (v_x ± √(σ_s(1 − v_x² + σ_s))) / (1 + σ_s).
  // c_s² via EOS dispatch (see gkyl_gr_euler_eos_cs2).
  double cs2_l = gkyl_gr_euler_eos_cs2(eos, pl.rho, p_l, h_l);
  double cs2_r = gkyl_gr_euler_eos_cs2(eos, pr.rho, p_r, h_r);
  if (cs2_l < 0.0) cs2_l = 0.0;
  if (cs2_r < 0.0) cs2_r = 0.0;
  // Cap c_s² strictly below 1 so σ_s is finite. For IDEAL with Γ ≤ 2 the
  // analytic bound is c_s² < Γ−1 ≤ 1; for APPROXIMATE_SYNGE the analytic
  // upper bound is the radiation-fluid limit c_s² → 1/3 as θ → ∞ (well
  // below 1). The clamp here catches numerical near-saturation in
  // degenerate configurations of the IDEAL branch (Γ → 2, ultra-rel).
  if (cs2_l > 1.0 - 1.0e-12) cs2_l = 1.0 - 1.0e-12;
  if (cs2_r > 1.0 - 1.0e-12) cs2_r = 1.0 - 1.0e-12;
  double sigma_l = cs2_l / (W_l*W_l * (1.0 - cs2_l));
  double sigma_r = cs2_r / (W_r*W_r * (1.0 - cs2_r));
  double rad_l = sqrt(fmax(sigma_l * (1.0 - vx_l*vx_l + sigma_l), 0.0));
  double rad_r = sqrt(fmax(sigma_r * (1.0 - vx_r*vx_r + sigma_r), 0.0));
  double lm_l = (vx_l - rad_l) / (1.0 + sigma_l);
  double lp_l = (vx_l + rad_l) / (1.0 + sigma_l);
  double lm_r = (vx_r - rad_r) / (1.0 + sigma_r);
  double lp_r = (vx_r + rad_r) / (1.0 + sigma_r);
  double lambda_L = fmin(lm_l, lm_r);
  double lambda_R = fmax(lp_l, lp_r);

  // Mignone-Bodo conserved E := τ + D. SR fluxes in the x-direction
  // (Banyuls form, no √γ since we're in flat tetrad):
  //   F^x[D]   = D · v_x
  //   F^x[S_i] = S_i · v_x + p · δ_i^x
  //   F^x[τ]   = (τ + p) · v_x = m_x − D · v_x   (since m_x = S_x = ρhW²v_x)
  //   F^x[E]   = F^x[τ] + F^x[D] = m_x
  double E_l = tau_l + D_l;
  double E_r = tau_r + D_r;

  double FD_l  = D_l  * vx_l;
  double Fmx_l = Sx_l * vx_l + p_l;
  double Fmy_l = Sy_l * vx_l;
  double Fmz_l = Sz_l * vx_l;
  double FE_l  = Sx_l;   // m_x

  double FD_r  = D_r  * vx_r;
  double Fmx_r = Sx_r * vx_r + p_r;
  double Fmy_r = Sy_r * vx_r;
  double Fmz_r = Sz_r * vx_r;
  double FE_r  = Sx_r;

  // HLL averages for the components needed by the λ* quadratic (MB05
  // eqs 9, 11). Pre-compute lam_diff and short-circuit the rare
  // degenerate λ_L ≈ λ_R case to avoid divide-by-zero.
  double lam_diff = lambda_R - lambda_L;
  if (stat) {
    stat->hllc.last_lambda_L = lambda_L;
    stat->hllc.last_lambda_R = lambda_R;
  }
  if (fabs(lam_diff) < 1.0e-14) {
    for (int i = 0; i < 3 * 5; i++) waves_tet[i] = 0.0;
    speeds[0] = lambda_L;
    speeds[1] = 0.5 * (lambda_L + lambda_R);
    speeds[2] = lambda_R;
    if (stat) {
      stat->hllc.last_did_fallback = 1;
      stat->hllc.last_fallback_reason = 1;
      stat->hllc.fallback_calls++;
      stat->hllc.fallback_reason_hist[1]++;
    }
    return fmax(fabs(lambda_L), fabs(lambda_R));
  }
  double mx_hll  = (lambda_R * Sx_r  - lambda_L * Sx_l  + Fmx_l - Fmx_r) / lam_diff;
  double E_hll   = (lambda_R * E_r   - lambda_L * E_l   + FE_l  - FE_r ) / lam_diff;
  double Fmx_hll = (lambda_R * Fmx_l - lambda_L * Fmx_r + lambda_R * lambda_L * (Sx_r - Sx_l)) / lam_diff;
  double FE_hll  = (lambda_R * FE_l  - lambda_L * FE_r  + lambda_R * lambda_L * (E_r  - E_l )) / lam_diff;

  // Quadratic for λ* (MB05 eq 18): F_E^hll·(λ*)² − (E^hll + F_mx^hll)·λ* + m_x^hll = 0.
  // The physical root is the minus-sign root (proven in MB05 Appendix A).
  // Stabilized via citardauq when −b > 0 to avoid catastrophic cancellation
  // — see comment block at the top of this function (fix #1).
  double a = FE_hll;
  double b = -(E_hll + Fmx_hll);
  double c = mx_hll;
  double lambda_star;
  if (fabs(a) < 1.0e-14) {
    // Linear case: −b·λ* + c = 0  ⇒  λ* = c / (−b).
    lambda_star = c / (-b);
  } else {
    double disc = b*b - 4.0*a*c;
    if (disc < 0.0) disc = 0.0;
    double sqrt_disc = sqrt(disc);
    if (-b >= 0.0) {
      lambda_star = (2.0 * c) / (-b + sqrt_disc);
    } else {
      lambda_star = (-b - sqrt_disc) / (2.0 * a);
    }
  }

  // Cold-gas / degenerate-fan fallback to HLL (fix #3 in header). Detect
  // λ* outside [λ_L, λ_R] (with a small clamp band) or non-finite, and
  // collapse the contact wave to zero. The two outer waves then carry the
  // entire HLL-averaged jump: w0 = U_HLL − U_L, w2 = U_R − U_HLL.
  // Fallback policy. Only catch *real* numerical pathologies that would
  // make the star-state computation unsafe:
  //   (1) λ* not finite (sqrt of negative discriminant after clamp, or
  //       analogous numerical blowup).
  //   (2) λ* arbitrarily close to λ_L or λ_R, which would blow up
  //       1/(λ_L − λ*) or 1/(λ_R − λ*) in the star-state RH formulas.
  //
  // We deliberately do NOT degrade when λ* is "outside the bracket"
  // [λ_L, λ_R] in the sense λ* < λ_L or λ* > λ_R. The MB05 Appendix A
  // proof of λ_L ≤ λ* ≤ λ_R explicitly assumes λ_L < 0 < λ_R (interface
  // sits inside the Riemann fan). For supersonic flow (entire fan on
  // one side of the interface) the algebraic minus root naturally lands
  // outside the bracket, but the wave decomposition still satisfies
  // Σ s·w = ΔF and wave_prop's bin-by-sign qfluct correctly handles the
  // supersonic upwinding via amdq=0 / apdq=ΔF (or vice versa).
  //
  // The earlier overly-aggressive "λ* < λ_L + clamp || λ* > λ_R − clamp"
  // check caused every BHL bow-shock interface to silently degrade to
  // HLL — see SESSION_NOTES_3.md (HLLC investigation) for the empirical
  // probe that uncovered this.
  double scale = fabs(lambda_R) + fabs(lambda_L) + 1.0;
  double dist_L = fabs(lambda_L - lambda_star);
  double dist_R = fabs(lambda_R - lambda_star);
  double tol = 1.0e-12 * scale;
  bool degrade_to_hll = !isfinite(lambda_star)
    || (dist_L < tol)
    || (dist_R < tol);
  if (stat) stat->hllc.last_lambda_star = lambda_star;
  if (degrade_to_hll) {
    int reason = !isfinite(lambda_star) ? 2
               : (dist_L < tol)         ? 3
               : 4;
    if (stat) {
      stat->hllc.last_did_fallback = 1;
      stat->hllc.last_fallback_reason = reason;
      stat->hllc.fallback_calls++;
      stat->hllc.fallback_reason_hist[reason]++;
    }
    double D_hll  = (lambda_R * D_r  - lambda_L * D_l  + FD_l  - FD_r ) / lam_diff;
    double Sy_hll = (lambda_R * Sy_r - lambda_L * Sy_l + Fmy_l - Fmy_r) / lam_diff;
    double Sz_hll = (lambda_R * Sz_r - lambda_L * Sz_l + Fmz_l - Fmz_r) / lam_diff;
    double tau_hll = E_hll - D_hll;
    double *w0 = &waves_tet[0 * 5];
    double *w1 = &waves_tet[1 * 5];
    double *w2 = &waves_tet[2 * 5];
    w0[0] = D_hll   - D_l;
    w0[1] = mx_hll  - Sx_l;
    w0[2] = Sy_hll  - Sy_l;
    w0[3] = Sz_hll  - Sz_l;
    w0[4] = tau_hll - tau_l;
    for (int k = 0; k < 5; k++) w1[k] = 0.0;
    w2[0] = D_r   - D_hll;
    w2[1] = Sx_r  - mx_hll;
    w2[2] = Sy_r  - Sy_hll;
    w2[3] = Sz_r  - Sz_hll;
    w2[4] = tau_r - tau_hll;
    speeds[0] = lambda_L;
    speeds[1] = 0.5 * (lambda_L + lambda_R);
    speeds[2] = lambda_R;
    return fmax(fabs(lambda_L), fabs(lambda_R));
  }

  // Success path: HLLC star-state proceeds without fallback. Bump the
  // no-fallback bin so the histogram totals match the number of calls.
  if (stat) stat->hllc.fallback_reason_hist[0]++;

  // p* from the left-state Rankine-Hugoniot relation. Sign-corrected vs
  // MB05 eq (17) as printed (fix #2 in header):
  //   p* = (A·λ* − B) / (1 − λ·λ*)
  // with A = λ·E − m_x, B = m_x·(λ − v_x) − p (using the LEFT state).
  double A_l = lambda_L * E_l - Sx_l;
  double B_l = Sx_l * (lambda_L - vx_l) - p_l;
  double p_star = (A_l * lambda_star - B_l) / (1.0 - lambda_L * lambda_star);

  // Star states from Rankine-Hugoniot (MB05 eq 16) in Banyuls (D, S_i, τ).
  // The fallback above ensures (lambda - lambda_star) is bounded away from
  // zero on both sides, so these inverses are safe.
  double inv_L_l = 1.0 / (lambda_L - lambda_star);
  double inv_L_r = 1.0 / (lambda_R - lambda_star);

  double Ds_l   = D_l  * (lambda_L - vx_l) * inv_L_l;
  double Sxs_l  = (Sx_l * (lambda_L - vx_l) + p_star - p_l) * inv_L_l;
  double Sys_l  = Sy_l * (lambda_L - vx_l) * inv_L_l;
  double Szs_l  = Sz_l * (lambda_L - vx_l) * inv_L_l;
  double Es_l   = (E_l  * (lambda_L - vx_l) + p_star * lambda_star - p_l * vx_l) * inv_L_l;
  double taus_l = Es_l - Ds_l;

  double Ds_r   = D_r  * (lambda_R - vx_r) * inv_L_r;
  double Sxs_r  = (Sx_r * (lambda_R - vx_r) + p_star - p_r) * inv_L_r;
  double Sys_r  = Sy_r * (lambda_R - vx_r) * inv_L_r;
  double Szs_r  = Sz_r * (lambda_R - vx_r) * inv_L_r;
  double Es_r   = (E_r  * (lambda_R - vx_r) + p_star * lambda_star - p_r * vx_r) * inv_L_r;
  double taus_r = Es_r - Ds_r;

  // Conservative-state jumps across the three Riemann fans.
  double *w0 = &waves_tet[0 * 5];  // λ_L: U_L  → U_L*
  double *w1 = &waves_tet[1 * 5];  // λ*:  U_L* → U_R*  (contact)
  double *w2 = &waves_tet[2 * 5];  // λ_R: U_R* → U_R

  w0[0] = Ds_l   - D_l;
  w0[1] = Sxs_l  - Sx_l;
  w0[2] = Sys_l  - Sy_l;
  w0[3] = Szs_l  - Sz_l;
  w0[4] = taus_l - tau_l;

  w1[0] = Ds_r   - Ds_l;
  w1[1] = Sxs_r  - Sxs_l;
  w1[2] = Sys_r  - Sys_l;
  w1[3] = Szs_r  - Szs_l;
  w1[4] = taus_r - taus_l;

  w2[0] = D_r   - Ds_r;
  w2[1] = Sx_r  - Sxs_r;
  w2[2] = Sy_r  - Sys_r;
  w2[3] = Sz_r  - Szs_r;
  w2[4] = tau_r - taus_r;

  speeds[0] = lambda_L;
  speeds[1] = lambda_star;
  speeds[2] = lambda_R;

  return fmax(fabs(lambda_L), fabs(lambda_R));
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
gr_euler_tetrad_wall(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int i = 0; i < 5; i++) ghost[i] = skin[i];
  ghost[1] = -ghost[1];
}

static void
gr_euler_tetrad_no_slip(const struct gkyl_wv_eqn *eqn, double t, int nc,
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

  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
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
// Riemann-solver callbacks. All four variants (Lax/HLL/HLLC/Roe) share the
// same Gorard+25 tetrad-first scaffolding; the per-variant SR kernel and
// configuration live on the equation object so a single dispatch entry
// point handles every variant — and the equation pointer is self-contained
// for GPU dispatch.
// ---------------------------------------------------------------------------

// Forward declarations of the curved-frame Lax LOW_ORDER fallback used by
// wave_tetrad_dispatch / qfluct_tetrad_dispatch. Definitions live further
// down in this file.
static double wave_lax_curved(const struct gkyl_wv_eqn *eqn, const double *delta,
  const double *ql, const double *qr, double *waves, double *s);
static void qfluct_lax_curved(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, const double *waves, const double *s, double *amdq,
  double *apdq);

// Tetrad-first HIGH_ORDER worker. Reads sr_kernel + num_waves +
// excision_policy from the equation object (set in the constructor based
// on rp_type).
static double
wave_tetrad_high_order(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, double *waves, double *s)
{
  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
  struct gkyl_gr_euler_eos eos = grm->eos;
  int num_waves = grm->num_waves;

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);

  // Both-excised, or HLLC any-side-excised: zero waves with sentinel speeds.
  if ((excise_l && excise_r) ||
      (grm->excision_policy == GKYL_TETRAD_EXCISION_SHORT_CIRCUIT && (excise_l || excise_r))) {
    for (int k = 0; k < num_waves * 5; k++) waves[k] = 0.0;
    if (num_waves == 2) {
      s[0] = -pow(10.0, -8.0);
      s[1] =  pow(10.0, -8.0);
    } else {
      s[0] = -pow(10.0, -8.0);
      s[1] = 0.0;
      s[2] =  pow(10.0, -8.0);
    }
    return pow(10.0, -8.0);
  }

  // Interface geometry. Active-cell-only for excision-adjacent faces;
  // otherwise the consistent g_iface / inv(g_iface) / sqrt(det) trio
  // (Phase 0 Fix 2).
  double g_iface[3][3], inv_g_iface[3][3];
  double alpha_iface, shift_x_iface, sqrt_det_iface;
  if (excise_l || excise_r) {
    const double *prods_active = excise_l ? grm->prodr_local : grm->prodl_local;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        g_iface[i][j]     = prods_active[GKYL_GR_SP_GIJ + 3*i + j];
        inv_g_iface[i][j] = prods_active[GKYL_GR_SP_INV_GIJ + 3*i + j];
      }
    alpha_iface    = prods_active[GKYL_GR_SP_LAPSE];
    shift_x_iface  = prods_active[GKYL_GR_SP_SHIFT + 0];
    sqrt_det_iface = sqrt(prods_active[GKYL_GR_SP_SPATIAL_DET]);
  } else {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        g_iface[i][j] = 0.5 * (grm->prodl_local[GKYL_GR_SP_GIJ + 3*i + j]
                             + grm->prodr_local[GKYL_GR_SP_GIJ + 3*i + j]);
    double det_iface = gkyl_gr_euler_tetrad_invert_metric_3x3(g_iface, inv_g_iface);
    alpha_iface   = 0.5 * (grm->prodl_local[GKYL_GR_SP_LAPSE]
                         + grm->prodr_local[GKYL_GR_SP_LAPSE]);
    shift_x_iface = 0.5 * (grm->prodl_local[GKYL_GR_SP_SHIFT + 0]
                         + grm->prodr_local[GKYL_GR_SP_SHIFT + 0]);
    sqrt_det_iface = sqrt(det_iface);
  }

  // Gram-Schmidt-on-γ⁻¹ triad (contravariant-x aligned).
  double M[3][3], M_inv[3][3];
  gkyl_gr_euler_tetrad_build_triad_contravariant_x(
    g_iface, inv_g_iface, M, M_inv);

  // Forward transform with sqrt_det_iface on BOTH sides (Phase 0 Fix 1).
  double ql_tet[5], qr_tet[5];
  if (excise_l) {
    for (int k = 0; k < 5; k++) ql_tet[k] = 0.0;
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(qr, sqrt_det_iface, inv_g_iface, M_inv, qr_tet);
  } else if (excise_r) {
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(ql, sqrt_det_iface, inv_g_iface, M_inv, ql_tet);
    for (int k = 0; k < 5; k++) qr_tet[k] = 0.0;
  } else {
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(ql, sqrt_det_iface, inv_g_iface, M_inv, ql_tet);
    gkyl_gr_euler_tetrad_q_to_tetrad_contra(qr, sqrt_det_iface, inv_g_iface, M_inv, qr_tet);
  }

  // SR Riemann kernel in the tetrad frame, dispatched via the eqn-stored
  // function pointer. Sized for the maximum num_waves (3).
  double waves_tet[3 * 5], speeds_tet[3];
  grm->sr_kernel(eos, ql_tet, qr_tet, grm->auxfields.prim_status_wave_prop,
    waves_tet, speeds_tet);

  // Back-transform waves and speeds; track curved-frame max-abs speed.
  double maxs_curved = 0.0;
  for (int k = 0; k < num_waves; k++) {
    gkyl_gr_euler_tetrad_wave_to_curved_contra(&waves_tet[k * 5],
      sqrt_det_iface, M_inv, &waves[k * 5]);
    s[k] = gkyl_gr_euler_tetrad_speed_to_curved_contra(
      speeds_tet[k], alpha_iface, shift_x_iface, inv_g_iface[0][0]);
    if (fabs(s[k]) > maxs_curved) maxs_curved = fabs(s[k]);
  }
  return maxs_curved;
}

// Unified HIGH_ORDER vs LOW_ORDER waves dispatcher. LOW_ORDER routes to
// direct curved-frame Lax for s²-positivity preservation on flagged cells
// (SESSION_NOTES_3 §13, §14).
static double
wave_tetrad_dispatch(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr,
  const double phil, const double phir, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_tetrad_high_order(eqn, ql, qr, waves, s);
  return wave_lax_curved(eqn, delta, ql, qr, waves, s);
}

// Unified qfluct dispatcher. HIGH_ORDER: bin-by-sign over the variant's
// num_waves (Σ s·w = ΔF_iface holds at machine ε for all four solvers).
// LOW_ORDER: direct curved-frame Lax fluctuation (positivity sweep).
static void
qfluct_tetrad_dispatch(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
      struct wv_gr_euler_tetrad, eqn);
    int num_waves = grm->num_waves;
    for (int i = 0; i < 5; i++) { amdq[i] = 0.0; apdq[i] = 0.0; }
    for (int k = 0; k < num_waves; k++) {
      double sm = fmin(0.0, s[k]);
      double sp = fmax(0.0, s[k]);
      const double *wk = &waves[k * 5];
      for (int i = 0; i < 5; i++) {
        amdq[i] += sm * wk[i];
        apdq[i] += sp * wk[i];
      }
    }
    return;
  }
  qfluct_lax_curved(eqn, ql, qr, waves, s, amdq, apdq);
}

// ---------------------------------------------------------------------------
// Direct curved-frame Lax-Friedrichs. Bypasses the tetrad-first wrap entirely:
// computes the GR Banyuls flux (α√γ · ...) directly on each side, takes Δq
// in curved-frame conservatives, and forms Lax waves with ±amax bracket.
//
// Motivation (SESSION_NOTES_3 §13): the tetrad-first Lax preserves the
// flat-SR admissibility set (D_tet, S_tet, τ_tet) by Einfeldt's convex-
// combination argument, but the back-transform with M_inv·γ does NOT
// preserve the curved-frame admissibility set
//   A_γ = {(D, S_i, τ) : D > 0, τ ≥ 0, (D+τ)² > γ^{ij} S_i S_j}
// which is what wave_prop's check_inv actually tests. Empirically the
// tetrad-first Lax can produce s²_curved < 0 from admissible inputs in
// curved spacetime (see test_direct_state_lax on bhl-repair-s2-#2).
//
// Direct curved-frame Lax works in the curved conservative-variable space
// throughout, so the convex-combination cell update stays in A_γ (which
// is convex; see §13.1 of the session notes for the Cauchy-Schwarz proof).
// Used as the POSITIVITY_SWEEP cleanup fallback for all three rp_types.
//
// Excision handling: flux+flux_correction already returns zero flux for
// excised cells and max_abs_speed_curved returns ~1e-8, so passing through
// the standard pipeline gives the correct absorbing-BC behavior without
// any special-case logic. Both-sides-excised is short-circuited cheaply.
static double
wave_lax_curved(const struct gkyl_wv_eqn *eqn, const double *delta,
  const double *ql, const double *qr, double *waves, double *s)
{

  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
  struct gkyl_gr_euler_eos eos = grm->eos;
  double gas_gamma = eos.gas_gamma;

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (excise_l && excise_r) {
    for (int k = 0; k < 2 * 5; k++) waves[k] = 0.0;
    s[0] = -pow(10.0, -8.0);
    s[1] =  pow(10.0, -8.0);
    return pow(10.0, -8.0);
  }

  // Per-side GR Banyuls fluxes. flux_correction zeros the flux on
  // excised cells, giving the absorbing-BC contribution automatically.
  struct gkyl_gr_euler_prim_status *stat = grm->auxfields.prim_status_wave_prop;
  double fl_sr[5], fr_sr[5], fl_gr[5], fr_gr[5];
  gkyl_gr_euler_tetrad_flux(eos, ql, grm->prodl_local, stat, fl_sr);
  gkyl_gr_euler_tetrad_flux(eos, qr, grm->prodr_local, stat, fr_sr);
  gkyl_gr_euler_tetrad_flux_correction(eos, ql, grm->prodl_local, stat, fl_sr, fl_gr);
  gkyl_gr_euler_tetrad_flux_correction(eos, qr, grm->prodr_local, stat, fr_sr, fr_gr);

  // amax — full-3D max-abs characteristic speed (max over x/y/z). The
  // x-only variant was too tight in metrics with off-diagonal γ_xy and
  // a general-orientation 3-velocity: with x-only, the Lax envelope
  // failed to dominate the full curved-frame flux Jacobian's spectral
  // radius, and the convex-combination argument broke. Reproducer:
  // test_direct_state_lax_curved_near_horizon at (sx,sy,sz)=(1,1,0)
  // on bhl-repair-s2-#2 showed s²_new=-1.45e-1 (x-only) → +1.16
  // (full-3D). See SESSION_NOTES_3 §17 for the analysis. Excision
  // short-circuits to 1e-8.
  double amaxl = gkyl_gr_euler_tetrad_max_abs_speed(eos, ql, grm->prodl_local, stat);
  double amaxr = gkyl_gr_euler_tetrad_max_abs_speed(eos, qr, grm->prodr_local, stat);
  double amax = fmax(amaxl, amaxr);
  if (!(amax > 0.0)) {
    for (int k = 0; k < 2 * 5; k++) waves[k] = 0.0;
    s[0] = -pow(10.0, -8.0);
    s[1] =  pow(10.0, -8.0);
    return pow(10.0, -8.0);
  }

  // Direct curved-frame Lax wave decomposition:
  //   Δq = q_R − q_L                (curved-frame densitized conservatives)
  //   ΔF = F_R_gr − F_L_gr          (curved-frame Banyuls flux)
  //   w_0 = 0.5·(Δq − ΔF/amax),    s_0 = −amax
  //   w_1 = 0.5·(Δq + ΔF/amax),    s_1 = +amax
  // Σ w_k = Δq, Σ s_k·w_k = ΔF (flux jump in curved frame, exact).
  // Cell update q^{n+1} = q − dt/dx · ΔF, which is a convex combination
  // of q_{i-1}, q_i, q_{i+1} for CFL ≤ 1 — preserves A_γ since A_γ is
  // convex.
  double *w0 = &waves[0 * 5];
  double *w1 = &waves[1 * 5];
  for (int i = 0; i < 5; i++) {
    double dq = qr[i] - ql[i];
    double df = fr_gr[i] - fl_gr[i];
    w0[i] = 0.5 * (dq - df / amax);
    w1[i] = 0.5 * (dq + df / amax);
  }
  s[0] = -amax;
  s[1] = +amax;

  return amax;
}

static void
qfluct_lax_curved(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  // Same bin-by-sign as qfluct_lax — Lax has the same wave-decomposition
  // shape regardless of which frame the waves were computed in.
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);
  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}


static double
flux_jump_func(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, double *flux_jump)
{
  // Mirrors the regular mod variant's flux_jump_func but uses the corrected
  // (flat + GR-correction) flux to match packed wv_gr_euler_tetrad.c:1581.
  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
  struct gkyl_gr_euler_eos eos = grm->eos;

  // flux_jump_func is the F-wave callback. Our Q-wave production setup
  // doesn't invoke it (see notes in priv.h).
  struct gkyl_gr_euler_prim_status *stat = grm->auxfields.prim_status_wave_prop;
  double fl_sr[5], fr_sr[5];
  gkyl_gr_euler_tetrad_flux(eos, ql, grm->prodl_local, stat, fl_sr);
  gkyl_gr_euler_tetrad_flux(eos, qr, grm->prodr_local, stat, fr_sr);

  double fl[5], fr[5];
  gkyl_gr_euler_tetrad_flux_correction(eos, ql, grm->prodl_local, stat, fl_sr, fl);
  gkyl_gr_euler_tetrad_flux_correction(eos, qr, grm->prodr_local, stat, fr_sr, fr);

  bool excise_l = grm->prodl_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  bool excise_r = grm->prodr_local[GKYL_GR_SP_EXCISION] < pow(10.0, -8.0);
  if (!excise_l && !excise_r) {
    for (int m = 0; m < 5; m++) flux_jump[m] = fr[m] - fl[m];
  } else {
    for (int m = 0; m < 5; m++) flux_jump[m] = 0.0;
  }

  double amaxl = gkyl_gr_euler_tetrad_max_abs_speed(eos, ql, grm->prodl_local, stat);
  double amaxr = gkyl_gr_euler_tetrad_max_abs_speed(eos, qr, grm->prodr_local, stat);
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
  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
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
  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
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

  // Per-constraint independent tallies. Wave-prop and source-step calls
  // route to separate repair_status buckets, selected by eqn->cur_repair_ctx
  // (set by the caller). A non-zero count on the wave-prop side means
  // the post-positivity-sweep low-order (Lax/HLL) flux failed to produce
  // an admissible state — should be rare. Source counts dominate when
  // the SSP-RK3 forward-Euler pushes already-admissible cells across an
  // admissibility boundary inside a stage.
  struct gkyl_gr_euler_repair_status *rstat = (eqn->cur_repair_ctx == 1)
    ? grm->auxfields.repair_status_wave_prop
    : grm->auxfields.repair_status_source;
  if (rstat) {
    if (fixed & GR_EULER_REPAIR_D)   rstat->bad_D_fixes   += 1;
    if (fixed & GR_EULER_REPAIR_TAU) rstat->bad_tau_fixes += 1;
    if (fixed & GR_EULER_REPAIR_S2)  rstat->bad_s2_fixes  += 1;
  }

  q[0] = D   * sd;
  q[1] = Sx  * sd;
  q[2] = Sy  * sd;
  q[3] = Sz  * sd;
  q[4] = tau * sd;
}

static double
max_speed_func(const struct gkyl_wv_eqn *eqn, const double *q)
{
  struct wv_gr_euler_tetrad *grm = container_of((struct gkyl_wv_eqn *)eqn,
    struct wv_gr_euler_tetrad, eqn);
  if (!grm->auxfields.prods) return 1.0;
  long cidx = gkyl_range_idx(&grm->conf_range, grm->cur_cell_idx);
  const double *prods = gkyl_array_cfetch(grm->auxfields.prods, cidx);
  return gkyl_gr_euler_tetrad_max_abs_speed(grm->eos, q, prods,
    grm->auxfields.prim_status_wave_prop);
}

static inline void
gr_euler_tetrad_cons_to_diag(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *diag)
{
  for (int i = 0; i < 5; i++) diag[i] = qin[i];
}

static inline void
gr_euler_tetrad_source(const struct gkyl_wv_eqn *eqn, const double *qin,
  double *sout)
{
  // Integrated by moment_spacetime_coupling.
  for (int i = 0; i < 5; i++) sout[i] = 0.0;
}

// ---------------------------------------------------------------------------
// Constructors and lifetime.
// ---------------------------------------------------------------------------

void
gkyl_gr_euler_tetrad_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_gr_euler_tetrad *grm = container_of(base,
    struct wv_gr_euler_tetrad, eqn);
  // Instrumentation now lives in app-owned per-species
  // {prim_status, repair_status} buckets — see
  // gkyl_moment_app_gr_euler_print_status for the opt-in dump.
  gkyl_free(grm);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_new(double gas_gamma,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  return gkyl_wv_gr_euler_tetrad_inew(
    &(struct gkyl_wv_gr_euler_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .conf_range = *conf_range,
      .rp_type = WV_GR_EULER_TETRAD_RP_HLL,
      .use_gpu = use_gpu,
    });
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_inew(
  const struct gkyl_wv_gr_euler_tetrad_inp *inp)
{
  struct wv_gr_euler_tetrad *grm =
    gkyl_malloc(sizeof(struct wv_gr_euler_tetrad));

  grm->eqn.type = GKYL_EQN_GR_EULER_TETRAD;
  grm->eqn.num_equations = 5;
  grm->eqn.num_diag = 5;

  // Resolve the EOS. Two input pathways:
  //   (a) New callers populate inp->eos directly (type IDEAL or
  //       APPROXIMATE_SYNGE, plus use_rcc for the latter).
  //   (b) Legacy callers populate only inp->gas_gamma; inp->eos is then
  //       zero-initialized (type = IDEAL, gas_gamma = 0.0) and we copy
  //       inp->gas_gamma into the eos bundle so downstream code sees IDEAL.
  grm->eos = inp->eos;
  if (grm->eos.type == GR_EULER_EOS_IDEAL && grm->eos.gas_gamma == 0.0) {
    grm->eos.gas_gamma = inp->gas_gamma;
  }

  // Roe Riemann solver is IDEAL-gas only — the SR-Roe eigenstructure uses
  // the Eulderink-Mellema ideal-gas Jacobian. Reject the combination at
  // construction time.
  assert(inp->rp_type != WV_GR_EULER_TETRAD_RP_ROE
         || grm->eos.type == GR_EULER_EOS_IDEAL);

  grm->conf_range = inp->conf_range;
  grm->auxfields.prods                   = NULL;
  grm->auxfields.prim_status_wave_prop   = NULL;
  grm->auxfields.repair_status_wave_prop = NULL;
  grm->auxfields.repair_status_source    = NULL;
  grm->eqn.cur_repair_ctx = 0;
  grm->rot_call_parity = 0;
  for (int d = 0; d < GKYL_MAX_DIM; d++) {
    grm->cur_idxl[d] = 0;
    grm->cur_idxr[d] = 0;
    grm->cur_cell_idx[d] = 0;
  }

  // Configure the tetrad-first scaffolding. sr_kernel, num_waves, and
  // excision_policy are read at every Riemann solve by
  // wave_tetrad_high_order / qfluct_tetrad_dispatch — all four variants
  // share the same dispatch entry point, so the equation pointer carries
  // everything needed for GPU dispatch.
  if (inp->rp_type == WV_GR_EULER_TETRAD_RP_LAX) {
    grm->sr_kernel = gkyl_gr_euler_tetrad_sr_lax_minkowski;
    grm->num_waves = 2;
    grm->excision_policy = GKYL_TETRAD_EXCISION_ZERO_VACUUM;
  }
  else if (inp->rp_type == WV_GR_EULER_TETRAD_RP_ROE) {
    grm->sr_kernel = gkyl_gr_euler_tetrad_sr_roe_minkowski;
    grm->num_waves = 3;
    grm->excision_policy = GKYL_TETRAD_EXCISION_ZERO_VACUUM;
  }
  else if (inp->rp_type == WV_GR_EULER_TETRAD_RP_HLLC) {
    grm->sr_kernel = gkyl_gr_euler_tetrad_sr_hllc_minkowski;
    grm->num_waves = 3;
    grm->excision_policy = GKYL_TETRAD_EXCISION_SHORT_CIRCUIT;
  }
  else {  // default: HLL
    grm->sr_kernel = gkyl_gr_euler_tetrad_sr_hll_minkowski;
    grm->num_waves = 2;
    grm->excision_policy = GKYL_TETRAD_EXCISION_ZERO_VACUUM;
  }
  grm->eqn.num_waves = grm->num_waves;
  grm->eqn.waves_func = wave_tetrad_dispatch;
  grm->eqn.qfluct_func = qfluct_tetrad_dispatch;

  grm->eqn.ffluct_func = NULL;
  grm->eqn.flux_jump = flux_jump_func;
  grm->eqn.check_inv_func = check_inv;
  grm->eqn.max_speed_func = max_speed_func;
  grm->eqn.rotate_to_local_func = rot_to_local;
  grm->eqn.rotate_to_global_func = rot_to_global;
  grm->eqn.wall_bc_func = gr_euler_tetrad_wall;
  grm->eqn.no_slip_bc_func = gr_euler_tetrad_no_slip;
  grm->eqn.cons_to_riem = cons_to_riem;
  grm->eqn.riem_to_cons = riem_to_cons;
  grm->eqn.cons_to_diag = gr_euler_tetrad_cons_to_diag;
  grm->eqn.source_func = gr_euler_tetrad_source;

  grm->eqn.set_interface_idx_func = gr_euler_tetrad_set_interface_idx;
  grm->eqn.set_cell_idx_func = gr_euler_tetrad_set_cell_idx;
  grm->eqn.repair_state_func = repair_state;

  grm->eqn.embed_geo = NULL;
  grm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(grm->eqn.flags);
  grm->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_tetrad_free);
  grm->eqn.on_dev = &grm->eqn;

  return &grm->eqn;
}

struct gkyl_gr_euler_eos
gkyl_wv_gr_euler_tetrad_eos(const struct gkyl_wv_eqn *eqn)
{
  const struct wv_gr_euler_tetrad *grm = container_of(eqn,
    struct wv_gr_euler_tetrad, eqn);
  return grm->eos;
}
