#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_elasticity.h>
#include <gkyl_wv_elasticity_priv.h>

void
gkyl_elasticity_prim_vars(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param,
  const double q[14], double v[14])
{
  double rho = q[0];
  double mom_x = q[1];
  double mom_y = q[2];
  double mom_z = q[3];
  
  double deformation_gradient_cons[3][3];
  deformation_gradient_cons[0][0] = q[4]; deformation_gradient_cons[0][1] = q[5]; deformation_gradient_cons[0][2] = q[6];
  deformation_gradient_cons[1][0] = q[7]; deformation_gradient_cons[1][1] = q[8]; deformation_gradient_cons[1][2] = q[9];
  deformation_gradient_cons[2][0] = q[10]; deformation_gradient_cons[2][1] = q[11]; deformation_gradient_cons[2][2] = q[12];

  double E_tot = q[13];

  double vel_x = mom_x / rho;
  double vel_y = mom_y / rho;
  double vel_z = mom_z / rho;

  double deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient[i][j] = deformation_gradient_cons[i][j] / rho;
    }
  }

  double internal_energy = (E_tot / rho) - (0.5 * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)));

  double bulk_modulus = (sound_speed * sound_speed) - ((4.0 / 3.0) * (shear_speed * shear_speed));
  double shear_modulus = shear_speed * shear_speed;

  double **inv_deformation_gradient = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_deformation_gradient[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_elasticity_inv_deformation_gradient(q, &inv_deformation_gradient);

  double inv_deformation_gradient_transpose[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose[i][j] = inv_deformation_gradient[j][i];
    }
  }

  double strain_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor[i][j] += inv_deformation_gradient_transpose[i][k] * inv_deformation_gradient[k][j];
      }
    }
  }

  double strain_tensor_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq[i][j] += strain_tensor[i][k] * strain_tensor[k][j];
      }
    }
  }

  double strain_invariant1 = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1 += strain_tensor[i][i];
  }

  double sq_strain_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace += strain_tensor_sq[i][i];
  }

  double strain_invariant2 = 0.5 * ((strain_invariant1 * strain_invariant1) - sq_strain_trace);

  double strain_invariant3 = (strain_tensor[0][0] * ((strain_tensor[1][1] * strain_tensor[2][2]) - (strain_tensor[2][1] * strain_tensor[1][2]))) -
    (strain_tensor[0][1] * ((strain_tensor[1][0] * strain_tensor[2][2]) - (strain_tensor[1][2] * strain_tensor[2][0]))) +
    (strain_tensor[0][2] * ((strain_tensor[1][0] * strain_tensor[2][1]) - (strain_tensor[1][1] * strain_tensor[2][0])));
  
  double specific_entropy = heat_capacity * log((1.0 / (6.0 * (alpha_param * alpha_param) * heat_capacity * T_ref)) * (pow(strain_invariant3, -0.5 * gamma_param) *
    ((-3.0 * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0) * bulk_modulus) +
    ((alpha_param * alpha_param) * ((6.0 * internal_energy) - (shear_modulus * ((strain_invariant1 * strain_invariant1) - (3.0 * strain_invariant2)) *
    pow(strain_invariant3, 0.5 * beta_param)) + (6.0 * heat_capacity * pow(strain_invariant3, 0.5 * gamma_param) * T_ref))))));
  
  v[0] = rho;
  v[1] = vel_x;
  v[2] = vel_y;
  v[3] = vel_z;

  v[4] = deformation_gradient[0][0]; v[5] = deformation_gradient[0][1]; v[6] = deformation_gradient[0][2];
  v[7] = deformation_gradient[1][0]; v[8] = deformation_gradient[1][1]; v[9] = deformation_gradient[1][2];
  v[10] = deformation_gradient[2][0]; v[11] = deformation_gradient[2][1]; v[12] = deformation_gradient[2][2];

  v[13] = specific_entropy;

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_deformation_gradient[i]);
  }
  gkyl_free(inv_deformation_gradient);
}

void 
gkyl_elasticity_inv_deformation_gradient(const double q[14], double ***inv_deformation_gradient)
{
  double rho = q[0];

  double deformation_gradient_cons[3][3];
  deformation_gradient_cons[0][0] = q[4]; deformation_gradient_cons[0][1] = q[5]; deformation_gradient_cons[0][2] = q[6];
  deformation_gradient_cons[1][0] = q[7]; deformation_gradient_cons[1][1] = q[8]; deformation_gradient_cons[1][2] = q[9];
  deformation_gradient_cons[2][0] = q[10]; deformation_gradient_cons[2][1] = q[11]; deformation_gradient_cons[2][2] = q[12];

  double deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient[i][j] = deformation_gradient_cons[i][j] / rho;
    }
  }

  double deformation_gradient_det = (deformation_gradient[0][0] * ((deformation_gradient[1][1] * deformation_gradient[2][2]) - (deformation_gradient[2][1] * deformation_gradient[1][2]))) -
    (deformation_gradient[0][1] * ((deformation_gradient[1][0] * deformation_gradient[2][2]) - (deformation_gradient[1][2] * deformation_gradient[2][0]))) +
    (deformation_gradient[0][2] * ((deformation_gradient[1][0] * deformation_gradient[2][1]) - (deformation_gradient[1][1] * deformation_gradient[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += deformation_gradient[i][i];
  }

  double deformation_gradient_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq[i][j] += deformation_gradient[i][k] * deformation_gradient[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += deformation_gradient_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_deformation_gradient)[i][j] = (1.0 / deformation_gradient_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * deformation_gradient[i][j]) + deformation_gradient_sq[i][j]);
    }
  }
}

void
gkyl_wv_elasticity_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_elasticity *elasticity = container_of(base->on_dev, struct wv_elasticity, eqn);
    gkyl_cu_free(elasticity);
  }

  struct wv_elasticity *elasticity = container_of(base, struct wv_elasticity, eqn);
  gkyl_free(elasticity);
}

struct gkyl_wv_eqn*
gkyl_wv_elasticity_new(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param, bool use_gpu)
{
  return gkyl_wv_elasticity_inew(&(struct gkyl_wv_elasticity_inp) {
      .T_ref = T_ref,
      .sound_speed = sound_speed,
      .shear_speed = shear_speed,
      .heat_capacity = heat_capacity,
      .alpha_param = alpha_param,
      .beta_param = beta_param,
      .gamma_param = gamma_param,
      .rp_type = WV_ELASTICITY_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_elasticity_inew(const struct gkyl_wv_elasticity_inp* inp)
{
  struct wv_elasticity *elasticity = gkyl_malloc(sizeof(struct wv_elasticity));

  elasticity->eqn.type = GKYL_EQN_ELASTICITY;
  elasticity->eqn.num_equations = 14;
  elasticity->eqn.num_diag = 14;

  elasticity->T_ref = inp->T_ref;
  elasticity->sound_speed = inp->sound_speed;
  elasticity->shear_speed = inp->shear_speed;
  elasticity->heat_capacity = inp->heat_capacity;

  elasticity->alpha_param = inp->alpha_param;
  elasticity->beta_param = inp->beta_param;
  elasticity->gamma_param = inp->gamma_param;

  elasticity->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(elasticity->eqn.flags);
  elasticity->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_elasticity_free);
  elasticity->eqn.on_dev = &elasticity->eqn; // On the CPU, the equation object points to itself.

  elasticity->eqn.embed_geo = NULL;

  return &elasticity->eqn;
}

double
gkyl_wv_elasticity_T_ref(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;

  return T_ref;
}

double
gkyl_wv_elasticity_sound_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double sound_speed = elasticity->sound_speed;

  return sound_speed;
}

double
gkyl_wv_elasticity_shear_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double shear_speed = elasticity->shear_speed;

  return shear_speed;
}

double
gkyl_wv_elasticity_heat_capacity(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double heat_capacity = elasticity->heat_capacity;

  return heat_capacity;
}

double
gkyl_wv_elasticity_alpha_param(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double alpha_param = elasticity->alpha_param;

  return alpha_param;
}

double
gkyl_wv_elasticity_beta_param(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double beta_param = elasticity->beta_param;

  return beta_param;
}

double
gkyl_wv_elasticity_gamma_param(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double gamma_param = elasticity->gamma_param;

  return gamma_param;
}