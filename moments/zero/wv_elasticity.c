#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_elasticity.h>
#include <gkyl_wv_elasticity_priv.h>

void
gkyl_elasticity_flux(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param,
  const double q[14], double flux[14])
{
  double v[14];
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q, v);
  
  double rho = v[0];
  double vel_x = v[1];
  double vel_y = v[2];
  double vel_z = v[3];

  double deformation_gradient[3][3];
  deformation_gradient[0][0] = v[4]; deformation_gradient[0][1] = v[5]; deformation_gradient[0][2] = v[6];
  deformation_gradient[1][0] = v[7]; deformation_gradient[1][1] = v[8]; deformation_gradient[1][2] = v[9];
  deformation_gradient[2][0] = v[10]; deformation_gradient[2][1] = v[11]; deformation_gradient[2][2] = v[12];

  double E_tot = q[13];
  double specific_entropy = v[13];

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
  
  double energy_deriv_invariant1 = (shear_modulus / 3.0) * strain_invariant1 * pow(strain_invariant3, 0.5 * beta_param);
  double energy_deriv_invariant2 = -(0.5 * shear_modulus) * pow(strain_invariant3, 0.5 * beta_param);

  double energy_deriv_invariant3 = (1.0 / (12.0 * alpha_param * strain_invariant3)) * (alpha_param * shear_modulus * beta_param *
    ((strain_invariant1 * strain_invariant1) - (3.0 * strain_invariant2)) * pow(strain_invariant3, 0.5 * beta_param));
  energy_deriv_invariant3 -= (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * pow(strain_invariant3, 0.5 * alpha_param) * bulk_modulus);
  energy_deriv_invariant3 += (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * pow(strain_invariant3, alpha_param) * bulk_modulus);
  energy_deriv_invariant3 += (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * alpha_param * heat_capacity * (exp(specific_entropy / heat_capacity) - 1.0) *
    gamma_param * pow(strain_invariant3, 0.5 * gamma_param) * T_ref);

  double invariant1_deriv_strain[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        invariant1_deriv_strain[i][j] = 1.0;
      }
      else {
        invariant1_deriv_strain[i][j] = 0.0;
      }
    }
  }

  double invariant2_deriv_strain[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      invariant2_deriv_strain[i][j] = (strain_invariant1 * invariant1_deriv_strain[i][j]) - strain_tensor[i][j];
    }
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

  double inv_strain_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_strain_tensor[i][j] = (1.0 / strain_invariant3) *
        ((0.5 * ((strain_invariant1 * strain_invariant1) - sq_strain_trace) * identity_tensor[i][j]) - (strain_invariant1 * strain_tensor[i][j]) + strain_tensor_sq[i][j]);
    }
  }

  double invariant3_deriv_strain[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      invariant3_deriv_strain[i][j] = strain_invariant3 * inv_strain_tensor[i][j];
    }
  }

  double stress_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      stress_tensor[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant1 * invariant1_deriv_strain[j][k]);
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant2 * invariant2_deriv_strain[j][k]);
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant3 * invariant3_deriv_strain[j][k]);
      }
    }
  }

  flux[0] = rho * vel_x;
  flux[1] = (rho * vel_x * vel_x) - stress_tensor[0][0];
  flux[2] = (rho * vel_y * vel_x) - stress_tensor[1][0];
  flux[3] = (rho * vel_z * vel_x) - stress_tensor[2][0];

  flux[4] = (rho * deformation_gradient[0][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_x);
  flux[5] = (rho * deformation_gradient[0][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_x);
  flux[6] = (rho * deformation_gradient[0][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_x);
  flux[7] = (rho * deformation_gradient[1][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_y);
  flux[8] = (rho * deformation_gradient[1][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_y);
  flux[9] = (rho * deformation_gradient[1][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_y);
  flux[10] = (rho * deformation_gradient[2][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_z);
  flux[11] = (rho * deformation_gradient[2][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_z);
  flux[12] = (rho * deformation_gradient[2][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_z);

  flux[13] = (E_tot * vel_x) - (vel_x * stress_tensor[0][0]) - (vel_y * stress_tensor[1][0]) - (vel_z * stress_tensor[2][0]);

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_deformation_gradient[i]);
  }
  gkyl_free(inv_deformation_gradient);
}

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
      (*inv_deformation_gradient)[i][j] = (1.0 / deformation_gradient_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient[i][j]) + deformation_gradient_sq[i][j]);
    }
  }
}

static inline double
gkyl_elasticity_max_abs_speed(double T_ref, double sound_speed, double shear_speed, double heat_capacity, double alpha_param, double beta_param, double gamma_param,
  const double q[14])
{
  double v[14] = { 0.0 };
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q, v);

  double vel_x = v[1];
  double vel_y = v[2];
  double vel_z = v[3];

  return sqrt((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)) + fmax(fabs(sound_speed), fabs(shear_speed));
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 14; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 14; i++) {
    qout[i] = win[i];
  }
}

static void
elasticity_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 14; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
elasticity_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];

  for (int i = 4; i < 14; i++) {
    ghost[i] = skin[i];
  }
}
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[4] * norm[0]) + (qglobal[5] * norm[1]) + (qglobal[6] * norm[2]);
  r1[1] = (qglobal[4] * tau1[0]) + (qglobal[5] * tau1[1]) + (qglobal[6] * tau1[2]);
  r1[2] = (qglobal[4] * tau2[0]) + (qglobal[5] * tau2[1]) + (qglobal[6] * tau2[2]);

  r2[0] = (qglobal[7] * norm[0]) + (qglobal[8] * norm[1]) + (qglobal[9] * norm[2]);
  r2[1] = (qglobal[7] * tau1[0]) + (qglobal[8] * tau1[1]) + (qglobal[9] * tau1[2]);
  r2[2] = (qglobal[7] * tau2[0]) + (qglobal[8] * tau2[1]) + (qglobal[9] * tau2[2]);

  r3[0] = (qglobal[10] * norm[0]) + (qglobal[11] * norm[1]) + (qglobal[12] * norm[2]);
  r3[1] = (qglobal[10] * tau1[0]) + (qglobal[11] * tau1[1]) + (qglobal[12] * tau1[2]);
  r3[2] = (qglobal[10] * tau2[0]) + (qglobal[11] * tau2[1]) + (qglobal[12] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[4] = v1[0]; qlocal[5] = v1[1]; qlocal[6] = v1[2];
  qlocal[7] = v2[0]; qlocal[8] = v2[1]; qlocal[9] = v2[2];
  qlocal[10] = v3[0]; qlocal[11] = v3[1]; qlocal[12] = v3[2];

  qlocal[13] = qglobal[13];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[4] * norm[0]) + (qlocal[5] * tau1[0]) + (qlocal[6] * tau2[0]);
  r1[1] = (qlocal[4] * norm[1]) + (qlocal[5] * tau1[1]) + (qlocal[6] * tau2[1]);
  r1[2] = (qlocal[4] * norm[2]) + (qlocal[5] * tau1[2]) + (qlocal[6] * tau2[2]);

  r2[0] = (qlocal[7] * norm[0]) + (qlocal[8] * tau1[0]) + (qlocal[9] * tau2[0]);
  r2[1] = (qlocal[7] * norm[1]) + (qlocal[8] * tau1[1]) + (qlocal[9] * tau2[1]);
  r2[2] = (qlocal[7] * norm[2]) + (qlocal[8] * tau1[2]) + (qlocal[9] * tau2[2]);

  r3[0] = (qlocal[10] * norm[0]) + (qlocal[11] * tau1[0]) + (qlocal[12] * tau2[0]);
  r3[1] = (qlocal[10] * norm[1]) + (qlocal[11] * tau1[1]) + (qlocal[12] * tau2[1]);
  r3[2] = (qlocal[10] * norm[2]) + (qlocal[11] * tau1[2]) + (qlocal[12] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[4] = v1[0]; qglobal[5] = v1[1]; qglobal[6] = v1[2];
  qglobal[7] = v2[0]; qglobal[8] = v2[1]; qglobal[9] = v2[2];
  qglobal[10] = v3[0]; qglobal[11] = v3[1]; qglobal[12] = v3[2];

  qglobal[13] = qlocal[13];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;
  double sound_speed = elasticity->sound_speed;
  double shear_speed = elasticity->shear_speed;
  double heat_capacity = elasticity->heat_capacity;

  double alpha_param = elasticity->alpha_param;
  double beta_param = elasticity->beta_param;
  double gamma_param = elasticity->gamma_param;

  double sl = gkyl_elasticity_max_abs_speed(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql);
  double sr = gkyl_elasticity_max_abs_speed(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr);
  double amax = fmax(sl, sr);

  double fl[14], fr[14];
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql, fl);
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[14];
  for (int i = 0; i < 14; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[14];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 14; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;
  double sound_speed = elasticity->sound_speed;
  double shear_speed = elasticity->shear_speed;
  double heat_capacity = elasticity->heat_capacity;

  double alpha_param = elasticity->alpha_param;
  double beta_param = elasticity->beta_param;
  double gamma_param = elasticity->gamma_param;

  double vl[14] = { 0.0 };
  double vr[14] = { 0.0 };
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql, vl);
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr, vr);

  double vx_l = vl[1];
  double vx_r = vr[1];

  double max_eig_l = fmax(fabs(sound_speed), fabs(shear_speed));
  double max_eig_r = fmax(fabs(sound_speed), fabs(shear_speed));

  double vx_avg = 0.5 * (vx_l + vx_r);
  double max_eig_avg = 0.5 * (max_eig_l + max_eig_r);

  double sl = fmin(vx_l - max_eig_l, vx_r - max_eig_r);
  double sr = fmax(vx_l + max_eig_l, vx_r + max_eig_r);

  double fl[14] = { 0.0 };
  double fr[14] = { 0.0 };
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql, fl);
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr, fr);

  double qm[14] = { 0.0 };
  double *w0 = &waves[0], *w1 = &waves[14];
  for (int i = 0; i < 14; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);

    w0[i] = qm[i] - ql[i];
    w1[i] = qr[i] - qm[i];
  }

  s[0] = sl;
  s[1] = sr;

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[14];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 14; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;
  double sound_speed = elasticity->sound_speed;
  double shear_speed = elasticity->shear_speed;
  double heat_capacity = elasticity->heat_capacity;

  double alpha_param = elasticity->alpha_param;
  double beta_param = elasticity->beta_param;
  double gamma_param = elasticity->gamma_param;

  double fr[14], fl[14];
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql, fl);
  gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr, fr);

  for (int m = 0; m < 14; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }
  
  double amaxl = gkyl_elasticity_max_abs_speed(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql);
  double amaxr = gkyl_elasticity_max_abs_speed(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;
  double sound_speed = elasticity->sound_speed;
  double shear_speed = elasticity->shear_speed;
  double heat_capacity = elasticity->heat_capacity;

  double alpha_param = elasticity->alpha_param;
  double beta_param = elasticity->beta_param;
  double gamma_param = elasticity->gamma_param;

  double v[14] = { 0.0 };
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q, v);

  if (v[0] < 0.0 || v[13] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_elasticity *elasticity = container_of(eqn, struct wv_elasticity, eqn);
  double T_ref = elasticity->T_ref;
  double sound_speed = elasticity->sound_speed;
  double shear_speed = elasticity->shear_speed;
  double heat_capacity = elasticity->heat_capacity;

  double alpha_param = elasticity->alpha_param;
  double beta_param = elasticity->beta_param;
  double gamma_param = elasticity->gamma_param;

  return gkyl_elasticity_max_abs_speed(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q);
}

static inline void
elasticity_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 14; i++) {
    diag[i] = qin[i];
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
      .rp_type = WV_ELASTICITY_RP_HLL,
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

  if (inp->rp_type == WV_ELASTICITY_RP_LAX) {
    elasticity->eqn.num_waves = 2;
    elasticity->eqn.waves_func = wave_lax_l;
    elasticity->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_ELASTICITY_RP_HLL) {
    elasticity->eqn.num_waves = 2;
    elasticity->eqn.waves_func = wave_hll_l;
    elasticity->eqn.qfluct_func = qfluct_hll_l;
  }

  elasticity->eqn.flux_jump = flux_jump;
  elasticity->eqn.check_inv_func = check_inv;
  elasticity->eqn.max_speed_func = max_speed;
  elasticity->eqn.rotate_to_local_func = rot_to_local;
  elasticity->eqn.rotate_to_global_func = rot_to_global;

  elasticity->eqn.wall_bc_func = elasticity_wall;
  elasticity->eqn.no_slip_bc_func = elasticity_no_slip;

  elasticity->eqn.cons_to_riem = cons_to_riem;
  elasticity->eqn.riem_to_cons = riem_to_cons;

  elasticity->eqn.cons_to_diag = elasticity_cons_to_diag;

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