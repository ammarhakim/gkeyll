#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_elasticity_rgfm.h>
#include <gkyl_wv_elasticity_rgfm_priv.h>

void
gkyl_elasticity_rgfm_prim_vars(int num_species, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, const double* q, double* v)
{
  double rho_total = q[0];
  double mom_x_total = q[1];
  double mom_y_total = q[2];
  double mom_z_total = q[3];
  
  double deformation_gradient_total_cons[3][3];
  deformation_gradient_total_cons[0][0] = q[4]; deformation_gradient_total_cons[0][1] = q[5]; deformation_gradient_total_cons[0][2] = q[6];
  deformation_gradient_total_cons[1][0] = q[7]; deformation_gradient_total_cons[1][1] = q[8]; deformation_gradient_total_cons[1][2] = q[9];
  deformation_gradient_total_cons[2][0] = q[10]; deformation_gradient_total_cons[2][1] = q[11]; deformation_gradient_total_cons[2][2] = q[12];

  double E_tot = q[13];
  double reinit_param = q[13 + (2 * num_species)];

  double *level_set_cons_s = gkyl_malloc(sizeof(double[num_species - 1]));
  for (int i = 0; i < num_species - 1; i++) {
    level_set_cons_s[i] = q[14 + i];
  }

  double *rho_cons_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_cons_s[i] = q[13 + num_species + i];
  }

  double vel_x_total = mom_x_total / rho_total;
  double vel_y_total = mom_y_total / rho_total;
  double vel_z_total = mom_z_total / rho_total;

  double deformation_gradient_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_total[i][j] = deformation_gradient_total_cons[i][j] / rho_total;
    }
  }

  double internal_energy_total = (E_tot / rho_total) - (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total)));

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = level_set_cons_s[i] / rho_total;
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = rho_cons_s[i] / level_set_s[i];
  }

  double *bulk_modulus_s = gkyl_malloc(sizeof(double[num_species]));
  double *shear_modulus_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    bulk_modulus_s[i] = (sound_speed_s[i] * sound_speed_s[i]) - ((4.0 / 3.0) * (shear_speed_s[i] * shear_speed_s[i]));
    shear_modulus_s[i] = shear_speed_s[i] * shear_speed_s[i];
  }

  double **inv_deformation_gradient_total = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_deformation_gradient_total[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_elasticity_rgfm_inv_deformation_gradient_total(q, &inv_deformation_gradient_total);

  double inv_deformation_gradient_total_transpose[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_total_transpose[i][j] = inv_deformation_gradient_total[j][i];
    }
  }

  double strain_tensor_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_total[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_total[i][j] += inv_deformation_gradient_total_transpose[i][k] * inv_deformation_gradient_total[k][j];
      }
    }
  }

  double strain_tensor_total_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_total_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_total_sq[i][j] += strain_tensor_total[i][k] * strain_tensor_total[k][j];
      }
    }
  }

  double strain_invariant1_total = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1_total += strain_tensor_total[i][i];
  }

  double sq_strain_trace_total = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace_total += strain_tensor_total_sq[i][i];
  }

  double strain_invariant2_total = 0.5 * ((strain_invariant1_total * strain_invariant1_total) - sq_strain_trace_total);

  double strain_invariant3_total = (strain_tensor_total[0][0] * ((strain_tensor_total[1][1] * strain_tensor_total[2][2]) - (strain_tensor_total[2][1] * strain_tensor_total[1][2]))) -
    (strain_tensor_total[0][1] * ((strain_tensor_total[1][0] * strain_tensor_total[2][2]) - (strain_tensor_total[1][2] * strain_tensor_total[2][0]))) +
    (strain_tensor_total[0][2] * ((strain_tensor_total[1][0] * strain_tensor_total[2][1]) - (strain_tensor_total[1][1] * strain_tensor_total[2][0])));

  double *specific_entropy_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    specific_entropy_s[i] = heat_capacity_s[i] * log((1.0 / (6.0 * (alpha_param_s[i] * alpha_param_s[i]) * heat_capacity_s[i] * T_ref_s[i])) * (pow(strain_invariant3_total, -0.5 * gamma_param_s[i]) *
    ((-3.0 * (pow(strain_invariant3_total, 0.5 * alpha_param_s[i]) - 1.0) * (pow(strain_invariant3_total, 0.5 * alpha_param_s[i]) - 1.0) * bulk_modulus_s[i]) +
    ((alpha_param_s[i] * alpha_param_s[i]) * ((6.0 * internal_energy_total) - (shear_modulus_s[i] * ((strain_invariant1_total * strain_invariant1_total) - (3.0 * strain_invariant2_total)) *
    pow(strain_invariant3_total, 0.5 * beta_param_s[i])) + (6.0 * heat_capacity_s[i] * pow(strain_invariant3_total, 0.5 * gamma_param_s[i]) * T_ref_s[i]))))));
  }

  double specific_entropy_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    specific_entropy_total += level_set_s[i] * specific_entropy_s[i];
  }

  v[0] = rho_total;
  v[1] = vel_x_total;
  v[2] = vel_y_total;
  v[3] = vel_z_total;

  v[4] = deformation_gradient_total[0][0]; v[5] = deformation_gradient_total[0][1]; v[6] = deformation_gradient_total[0][2];
  v[7] = deformation_gradient_total[1][0]; v[8] = deformation_gradient_total[1][1]; v[9] = deformation_gradient_total[1][2];
  v[10] = deformation_gradient_total[2][0]; v[11] = deformation_gradient_total[2][1]; v[12] = deformation_gradient_total[2][2];

  v[13] = specific_entropy_total;

  for (int i = 0; i < num_species - 1; i++) {
    v[14 + i] = level_set_s[i];
  }
  for (int i = 0; i < num_species; i++) {
    v[13 + num_species + i] = rho_s[i];
  }
  v[13 + (2 * num_species)] = reinit_param;

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_deformation_gradient_total[i]);
  }
  gkyl_free(inv_deformation_gradient_total);

  gkyl_free(level_set_cons_s);
  gkyl_free(rho_cons_s);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
  gkyl_free(bulk_modulus_s);
  gkyl_free(shear_modulus_s);
  gkyl_free(specific_entropy_s);
}

void 
gkyl_elasticity_rgfm_inv_deformation_gradient_total(const double* q, double ***inv_deformation_gradient_total)
{
  double rho_total = q[0];

  double deformation_gradient_total_cons[3][3];
  deformation_gradient_total_cons[0][0] = q[4]; deformation_gradient_total_cons[0][1] = q[5]; deformation_gradient_total_cons[0][2] = q[6];
  deformation_gradient_total_cons[1][0] = q[7]; deformation_gradient_total_cons[1][1] = q[8]; deformation_gradient_total_cons[1][2] = q[9];
  deformation_gradient_total_cons[2][0] = q[10]; deformation_gradient_total_cons[2][1] = q[11]; deformation_gradient_total_cons[2][2] = q[12];

  double deformation_gradient_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_total[i][j] = deformation_gradient_total_cons[i][j] / rho_total;
    }
  }

  double deformation_gradient_total_det = (deformation_gradient_total[0][0] * ((deformation_gradient_total[1][1] * deformation_gradient_total[2][2]) -
    (deformation_gradient_total[2][1] * deformation_gradient_total[1][2]))) -
    (deformation_gradient_total[0][1] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][2]) - (deformation_gradient_total[1][2] * deformation_gradient_total[2][0]))) +
    (deformation_gradient_total[0][2] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][1]) - (deformation_gradient_total[1][1] * deformation_gradient_total[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += deformation_gradient_total[i][i];
  }

  double deformation_gradient_total_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_total_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_total_sq[i][j] += deformation_gradient_total[i][k] * deformation_gradient_total[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += deformation_gradient_total_sq[i][i];
  }

  double identity_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        identity_tensor[i][j] = 1.0;
      }
      else {
        identity_tensor[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_deformation_gradient_total)[i][j] = (1.0 / deformation_gradient_total_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient_total[i][j]) + deformation_gradient_total_sq[i][j]);
    }
  }
}

void
gkyl_wv_elasticity_rgfm_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_elasticity_rgfm *elasticity_rgfm = container_of(base->on_dev, struct wv_elasticity_rgfm, eqn);
    gkyl_cu_free(elasticity_rgfm);
  }

  struct wv_elasticity_rgfm *elasticity_rgfm = container_of(base, struct wv_elasticity_rgfm, eqn);
  gkyl_free(elasticity_rgfm);
}

struct gkyl_wv_eqn*
gkyl_wv_elasticity_rgfm_new(int num_species, double* rho_ref_s, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, int reinit_freq, double surface_tension, bool use_gpu)
{
  return gkyl_wv_elasticity_rgfm_inew(&(struct gkyl_wv_elasticity_rgfm_inp) {
      .num_species = num_species,
      .rho_ref_s = rho_ref_s,
      .T_ref_s = T_ref_s,
      .sound_speed_s = sound_speed_s,
      .shear_speed_s = shear_speed_s,
      .heat_capacity_s = heat_capacity_s,
      .alpha_param_s = alpha_param_s,
      .beta_param_s = beta_param_s,
      .gamma_param_s = gamma_param_s,
      .reinit_freq = reinit_freq,
      .surface_tension = surface_tension,
      .rp_type = WV_ELASTICITY_RGFM_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_elasticity_rgfm_inew(const struct gkyl_wv_elasticity_rgfm_inp* inp)
{
  struct wv_elasticity_rgfm *elasticity_rgfm = gkyl_malloc(sizeof(struct wv_elasticity_rgfm));

  elasticity_rgfm->eqn.type = GKYL_EQN_ELASTICITY_RGFM;
  elasticity_rgfm->eqn.num_equations = 14 + (2 * inp->num_species);
  elasticity_rgfm->eqn.num_diag = 14 + (2 * inp->num_species);

  elasticity_rgfm->num_species = inp->num_species;

  elasticity_rgfm->rho_ref_s = inp->rho_ref_s;
  elasticity_rgfm->T_ref_s = inp->T_ref_s;
  elasticity_rgfm->sound_speed_s = inp->sound_speed_s;
  elasticity_rgfm->shear_speed_s = inp->shear_speed_s;
  elasticity_rgfm->heat_capacity_s = inp->heat_capacity_s;

  elasticity_rgfm->alpha_param_s = inp->alpha_param_s;
  elasticity_rgfm->beta_param_s = inp->beta_param_s;
  elasticity_rgfm->gamma_param_s = inp->gamma_param_s;

  elasticity_rgfm->reinit_freq = inp->reinit_freq;
  elasticity_rgfm->surface_tension = inp->surface_tension;

  elasticity_rgfm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(elasticity_rgfm->eqn.flags);
  elasticity_rgfm->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_elasticity_rgfm_free);
  elasticity_rgfm->eqn.on_dev = &elasticity_rgfm->eqn; // On the CPU, the equation object points to itself.

  elasticity_rgfm->eqn.embed_geo = NULL;

  return &elasticity_rgfm->eqn;
}

int
gkyl_wv_elasticity_rgfm_num_species(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  return num_species;
}

double*
gkyl_wv_elasticity_rgfm_rho_ref_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* rho_ref_s = elasticity_rgfm->rho_ref_s;

  return rho_ref_s;
}

double*
gkyl_wv_elasticity_rgfm_T_ref_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* T_ref_s = elasticity_rgfm->T_ref_s;

  return T_ref_s;
}

double*
gkyl_wv_elasticity_rgfm_sound_speed_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* sound_speed_s = elasticity_rgfm->sound_speed_s;

  return sound_speed_s;
}

double*
gkyl_wv_elasticity_rgfm_shear_speed_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* shear_speed_s = elasticity_rgfm->shear_speed_s;

  return shear_speed_s;
}

double*
gkyl_wv_elasticity_rgfm_heat_capacity_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* heat_capacity_s = elasticity_rgfm->heat_capacity_s;

  return heat_capacity_s;
}

double*
gkyl_wv_elasticity_rgfm_alpha_param_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* alpha_param_s = elasticity_rgfm->alpha_param_s;

  return alpha_param_s;
}

double*
gkyl_wv_elasticity_rgfm_beta_param_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* beta_param_s = elasticity_rgfm->beta_param_s;

  return beta_param_s;
}

double*
gkyl_wv_elasticity_rgfm_gamma_param_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double* gamma_param_s = elasticity_rgfm->gamma_param_s;

  return gamma_param_s;
}

int
gkyl_wv_elasticity_rgfm_reinit_freq(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int reinit_freq = elasticity_rgfm->reinit_freq;

  return reinit_freq;
}

double
gkyl_wv_elasticity_rgfm_surface_tension(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  double surface_tension = elasticity_rgfm->surface_tension;

  return surface_tension;
}