#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_mat.h>

gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp)
{
  gkyl_moment_em_coupling *mom_em = gkyl_malloc(sizeof(gkyl_moment_em_coupling));

  mom_em->grid = *(inp.grid);
  mom_em->ndim = mom_em->grid.ndim;
  mom_em->nfluids = inp.nfluids;
  for (int i = 0; i < mom_em->nfluids; i++) {
    mom_em->param[i] = inp.param[i];
  }

  mom_em->epsilon0 = inp.epsilon0;
  mom_em->mu0 = inp.mu0;

  mom_em->is_charged_species = false;
  for (int i = 0; i < mom_em->nfluids; i++) {
    if (inp.param[i].charge != 0.0) {
      mom_em->is_charged_species = true;
    }
  }

  mom_em->static_field = inp.static_field; 
  mom_em->t_ramp_E = inp.t_ramp_E;
  if (mom_em->t_ramp_E != 0.0) {
    mom_em->ramp_app_E = true;
  }
  else {
    mom_em->ramp_app_E = false;
  }
  mom_em->t_ramp_curr = inp.t_ramp_curr;
  if (mom_em->t_ramp_curr != 0.0) {
    mom_em->ramp_app_curr = true;
  }
  else {
    mom_em->ramp_app_curr = false;
  }

  mom_em->has_collision = inp.has_collision;
  mom_em->use_rel = inp.use_rel;

  if (mom_em->has_collision) {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = inp.nu_base[i][j];
      }
    }
  }
  else {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = 0.0;
      }
    }
  }

  mom_em->use_explicit_em_coupling = inp.use_explicit_em_coupling;

  mom_em->has_nT_sources = inp.has_nT_sources;

  // Source-family parameter bundles travel by struct assignment; each
  // family's parameters are read only when its 'enabled' flag is set.
  mom_em->friction = inp.friction;
  mom_em->volume_sources = inp.volume_sources;
  mom_em->reactivity = inp.reactivity;
  mom_em->einstein_medium = inp.einstein_medium;
  mom_em->gr_ultra_rel = inp.gr_ultra_rel;
  mom_em->gr_euler = inp.gr_euler;
  mom_em->gr_twofluid = inp.gr_twofluid;
  mom_em->vacuum_einstein = inp.vacuum_einstein;
  mom_em->vacuum_einstein_conformal = inp.vacuum_einstein_conformal;
  mom_em->gr_mhd = inp.gr_mhd;

  mom_em->has_gr_em_coupling = inp.has_gr_em_coupling;
  if (mom_em->has_gr_em_coupling) {
    for (int i = 0; i < mom_em->nfluids; i++) {
      mom_em->gr_em_eos[i] = inp.gr_em_eos[i];
    }
  }

  return mom_em;
}

void
gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
      p_rhs_s[i] = gkyl_array_cfetch(p_rhs[i], cell_idx);
      nT_sources_s[i] = gkyl_array_cfetch(nT_sources[i], cell_idx);
    }

    double *em_arr = em ? gkyl_array_fetch(em, cell_idx) : 0;
    const double *app_current_arr = app_current ? gkyl_array_cfetch(app_current, cell_idx) : 0;
    const double *ext_em_arr = ext_em ? gkyl_array_cfetch(ext_em, cell_idx) : 0;

    implicit_source_coupling_update(mom_em, t_curr, dt, fluid_s, app_accel_s, p_rhs_s, em_arr, app_current_arr, ext_em_arr, nT_sources_s);
  }
}

void
gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj *proj_app_curr, int nstrang)
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  double dt_local = 2.0 * dt;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
    }

    double *em_arr = em ? em_arr = gkyl_array_fetch(em, cell_idx) : 0;
    const double *app_current_arr = app_current ? gkyl_array_cfetch(app_current, cell_idx) : 0;
    const double *app_current1_arr = app_current1 ? gkyl_array_cfetch(app_current1, cell_idx) : 0;
    const double *app_current2_arr = app_current2 ? gkyl_array_cfetch(app_current2, cell_idx) : 0;
    const double *ext_em_arr = ext_em ? gkyl_array_cfetch(ext_em, cell_idx) : 0;

    if (mom_em->use_rel) {
      explicit_source_coupling_update(mom_em, t_curr, dt_local, fluid_s, app_accel_s, em_arr, app_current_arr, app_current1_arr, app_current2_arr,
        ext_em_arr, nstrang);
    }
  }
}

struct gkyl_moment_em_coupling_status
gkyl_moment_em_coupling_gr_em_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em, const struct gkyl_array* ext_em)
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];

  double omega_max = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
    }

    double *em_arr = gkyl_array_fetch(em, cell_idx);
    const double *ext_em_arr = ext_em ? gkyl_array_cfetch(ext_em, cell_idx) : 0;

    double cell_omega = explicit_gr_em_source_update(mom_em, t_curr, dt, fluid_s, em_arr, ext_em_arr);
    if (cell_omega > omega_max) {
      omega_max = cell_omega;
    }
  }

  // The explicit SSP-RK3 solve is conditionally stable; report whether dt
  // resolved the plasma/cyclotron frequencies and a stable dt for retry.
  // On failure the app discards this (mutated) sub-step from its backup copy
  // and retries the whole step with dt_suggested.
  struct gkyl_moment_em_coupling_status status;
  status.omega_max = omega_max;
  status.success = (omega_max * dt) <= GR_EM_EXPLICIT_OMEGA_DT_MAX;
  status.dt_suggested = (omega_max > 0.0) ? (GR_EM_EXPLICIT_OMEGA_DT_MAX / omega_max) : DBL_MAX;

  return status;
}

void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em)
{
  gkyl_free(mom_em);
}