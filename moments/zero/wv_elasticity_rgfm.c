#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_elasticity_rgfm.h>
#include <gkyl_wv_elasticity_rgfm_priv.h>

void
gkyl_elasticity_rgfm_flux(int num_species, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, const double* q, double* flux)
{
  double *v = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  gkyl_elasticity_rgfm_prim_vars(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q, v);
  
  double rho_total = v[0];
  double vel_x_total = v[1];
  double vel_y_total = v[2];
  double vel_z_total = v[3];

  double deformation_gradient_total[3][3];
  deformation_gradient_total[0][0] = v[4]; deformation_gradient_total[0][1] = v[5]; deformation_gradient_total[0][2] = v[6];
  deformation_gradient_total[1][0] = v[7]; deformation_gradient_total[1][1] = v[8]; deformation_gradient_total[1][2] = v[9];
  deformation_gradient_total[2][0] = v[10]; deformation_gradient_total[2][1] = v[11]; deformation_gradient_total[2][2] = v[12];

  double E_tot = q[13];

  double *level_set_s = gkyl_malloc(sizeof(double[num_species]));
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_s[i] = v[14 + i];
    level_set_total += level_set_s[i];
  }
  level_set_s[num_species - 1] = 1.0 - level_set_total;

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[13 + num_species + i];
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

  double sq_strain_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace += strain_tensor_total_sq[i][i];
  }

  double strain_invariant2_total = 0.5 * ((strain_invariant1_total * strain_invariant1_total) - sq_strain_trace);

  double strain_invariant3_total = (strain_tensor_total[0][0] * ((strain_tensor_total[1][1] * strain_tensor_total[2][2]) - (strain_tensor_total[2][1] * strain_tensor_total[1][2]))) -
    (strain_tensor_total[0][1] * ((strain_tensor_total[1][0] * strain_tensor_total[2][2]) - (strain_tensor_total[1][2] * strain_tensor_total[2][0]))) +
    (strain_tensor_total[0][2] * ((strain_tensor_total[1][0] * strain_tensor_total[2][1]) - (strain_tensor_total[1][1] * strain_tensor_total[2][0])));
  
  double *energy_deriv_invariant1_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    energy_deriv_invariant1_s[i] = (shear_modulus_s[i] / 3.0) * strain_invariant1_total * pow(strain_invariant3_total, 0.5 * beta_param_s[i]);
  }

  double *energy_deriv_invariant2_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    energy_deriv_invariant2_s[i] = -(0.5 * shear_modulus_s[i]) * pow(strain_invariant3_total, 0.5 * beta_param_s[i]);
  }

  double internal_energy_total = (E_tot / rho_total) - (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total)));
  double *specific_entropy_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    specific_entropy_s[i] = heat_capacity_s[i] * log((1.0 / (6.0 * (alpha_param_s[i] * alpha_param_s[i]) * heat_capacity_s[i] * T_ref_s[i])) * (pow(strain_invariant3_total, -0.5 * gamma_param_s[i]) *
      ((-3.0 * (pow(strain_invariant3_total, 0.5 * alpha_param_s[i]) - 1.0) * (pow(strain_invariant3_total, 0.5 * alpha_param_s[i]) - 1.0) * bulk_modulus_s[i]) +
      ((alpha_param_s[i] * alpha_param_s[i]) * ((6.0 * internal_energy_total) - (shear_modulus_s[i] * ((strain_invariant1_total * strain_invariant1_total) - (3.0 * strain_invariant2_total)) *
      pow(strain_invariant3_total, 0.5 * beta_param_s[i])) + (6.0 * heat_capacity_s[i] * pow(strain_invariant3_total, 0.5 * gamma_param_s[i]) * T_ref_s[i]))))));
  }

  double *energy_deriv_invariant3_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    energy_deriv_invariant3_s[i] = (1.0 / (12.0 * alpha_param_s[i] * strain_invariant3_total)) * (alpha_param_s[i] * shear_modulus_s[i] * beta_param_s[i] *
      ((strain_invariant1_total * strain_invariant1_total) - (3.0 * strain_invariant2_total)) * pow(strain_invariant3_total, 0.5 * beta_param_s[i]));
    
    energy_deriv_invariant3_s[i] -= (1.0 / (12.0 * alpha_param_s[i] * strain_invariant3_total)) * (6.0 * pow(strain_invariant3_total, 0.5 * alpha_param_s[i]) * bulk_modulus_s[i]);
    energy_deriv_invariant3_s[i] += (1.0 / (12.0 * alpha_param_s[i] * strain_invariant3_total)) * (6.0 * pow(strain_invariant3_total, alpha_param_s[i]) * bulk_modulus_s[i]);
    energy_deriv_invariant3_s[i] += (1.0 / (12.0 * alpha_param_s[i] * strain_invariant3_total)) * (6.0 * alpha_param_s[i] * heat_capacity_s[i] * (exp(specific_entropy_s[i] / heat_capacity_s[i]) - 1.0) *
      gamma_param_s[i] * pow(strain_invariant3_total, 0.5 * gamma_param_s[i]) * T_ref_s[i]);
  }

  double invariant1_deriv_strain_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        invariant1_deriv_strain_total[i][j] = 1.0;
      }
      else {
        invariant1_deriv_strain_total[i][j] = 0.0;
      }
    }
  }

  double invariant2_deriv_strain_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      invariant2_deriv_strain_total[i][j] = (strain_invariant1_total * invariant1_deriv_strain_total[i][j]) - strain_tensor_total[i][j];
    }
  }

  double inv_strain_tensor_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_strain_tensor_total[i][j] = (1.0 / strain_invariant3_total) *
        ((0.5 * ((strain_invariant1_total * strain_invariant1_total) - sq_strain_trace) * invariant1_deriv_strain_total[i][j]) - (strain_invariant1_total * strain_tensor_total[i][j]) + strain_tensor_total_sq[i][j]);
    }
  }

  double invariant3_deriv_strain_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      invariant3_deriv_strain_total[i][j] = strain_invariant3_total * inv_strain_tensor_total[i][j];
    }
  }

  double ***stress_tensor_s = gkyl_malloc(sizeof(double**[num_species]));
  for (int i = 0; i < num_species; i++) {
    stress_tensor_s[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      stress_tensor_s[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  for (int i = 0; i < num_species; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        stress_tensor_s[i][j][k] = 0.0;
      }
    }
  }

  for (int i = 0; i < num_species; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          stress_tensor_s[i][j][k] -= (2.0 * rho_total * strain_tensor_total[j][l]) * (energy_deriv_invariant1_s[i] * invariant1_deriv_strain_total[k][l]);
          stress_tensor_s[i][j][k] -= (2.0 * rho_total * strain_tensor_total[j][l]) * (energy_deriv_invariant2_s[i] * invariant2_deriv_strain_total[k][l]);
          stress_tensor_s[i][j][k] -= (2.0 * rho_total * strain_tensor_total[j][l]) * (energy_deriv_invariant3_s[i] * invariant3_deriv_strain_total[k][l]);
        }
      }
    }
  }

  double stress_tensor_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      stress_tensor_total[i][j] = 0.0;
    }
  }

  for (int i = 0; i < num_species; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        stress_tensor_total[j][k] += level_set_s[i] * stress_tensor_s[i][j][k];
      }
    }
  }

  flux[0] = rho_total * vel_x_total;
  flux[1] = (rho_total * vel_x_total * vel_x_total) - stress_tensor_total[0][0];
  flux[2] = (rho_total * vel_y_total * vel_x_total) - stress_tensor_total[1][0];
  flux[3] = (rho_total * vel_z_total * vel_x_total) - stress_tensor_total[2][0];

  flux[4] = (rho_total * deformation_gradient_total[0][0] * vel_x_total) - (rho_total * deformation_gradient_total[0][0] * vel_x_total);
  flux[5] = (rho_total * deformation_gradient_total[0][1] * vel_x_total) - (rho_total * deformation_gradient_total[0][1] * vel_x_total);
  flux[6] = (rho_total * deformation_gradient_total[0][2] * vel_x_total) - (rho_total * deformation_gradient_total[0][2] * vel_x_total);
  flux[7] = (rho_total * deformation_gradient_total[1][0] * vel_x_total) - (rho_total * deformation_gradient_total[0][0] * vel_y_total);
  flux[8] = (rho_total * deformation_gradient_total[1][1] * vel_x_total) - (rho_total * deformation_gradient_total[0][1] * vel_y_total);
  flux[9] = (rho_total * deformation_gradient_total[1][2] * vel_x_total) - (rho_total * deformation_gradient_total[0][2] * vel_y_total);
  flux[10] = (rho_total * deformation_gradient_total[2][0] * vel_x_total) - (rho_total * deformation_gradient_total[0][0] * vel_z_total);
  flux[11] = (rho_total * deformation_gradient_total[2][1] * vel_x_total) - (rho_total * deformation_gradient_total[0][1] * vel_z_total);
  flux[12] = (rho_total * deformation_gradient_total[2][2] * vel_x_total) - (rho_total * deformation_gradient_total[0][2] * vel_z_total);

  flux[13] = (E_tot * vel_x_total) - (vel_x_total * stress_tensor_total[0][0]) - (vel_y_total * stress_tensor_total[1][0]) - (vel_z_total * stress_tensor_total[2][0]);

  for (int i = 0; i < num_species - 1; i++) {
    flux[14 + i] = rho_total * (vel_x_total * level_set_s[i]);
  }
  for (int i = 0; i < num_species; i++) {
    flux[13 + num_species + i] = level_set_s[i] * (vel_x_total * rho_s[i]);
  }

  flux[13 + (2 * num_species)] = 0.0;

  gkyl_free(v);
  gkyl_free(level_set_s);
  gkyl_free(rho_s);
  gkyl_free(bulk_modulus_s);
  gkyl_free(shear_modulus_s);
  gkyl_free(energy_deriv_invariant1_s);
  gkyl_free(energy_deriv_invariant2_s);
  gkyl_free(energy_deriv_invariant3_s);
  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_deformation_gradient_total[i]);
  }
  gkyl_free(inv_deformation_gradient_total);

  for (int i = 0; i < num_species; i++) {
    for (int j = 0; j < 3; j++) {
      gkyl_free(stress_tensor_s[i][j]);
    }

    gkyl_free(stress_tensor_s[i]);
  }
  gkyl_free(stress_tensor_s);
}

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

static inline double
gkyl_elasticity_rgfm_max_abs_speed(int num_species, double* rho_ref_s, double* T_ref_s, double* sound_speed_s, double* shear_speed_s, double* heat_capacity_s,
  double* alpha_param_s, double* beta_param_s, double* gamma_param_s, const double* q)
{
  double *v = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  gkyl_elasticity_rgfm_prim_vars(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q, v);

  double vel_x_total = v[1];

  double max_eig = 0.0;
  for (int i = 0; i < num_species; i++) {
    if (fabs(sound_speed_s[i]) > max_eig) {
      max_eig = fabs(sound_speed_s[i]);
    }
    if (fabs(shear_speed_s[i]) > max_eig) {
      max_eig = fabs(shear_speed_s[i]);
    }
  }

  gkyl_free(v);

  return fabs(vel_x_total) + fabs(max_eig);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 14 + (2 * num_species); i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 14 + (2 * num_species); i++) {
    qout[i] = win[i];
  }
}

static void
elasticity_rgfm_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  for (int i = 0; i < 14 + (2 * num_species); i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
elasticity_rgfm_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];

  for (int i = 4; i < 14 + (2 * num_species); i++) {
    ghost[i] = skin[i];
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  for (int i = 0; i < 14 + (2 * num_species); i++) {
    qlocal[i] = qglobal[i];
  }

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

  // Rotate deformation gradient tensor to local coordinate frame.
  // NOTE: We are taking the transpose of the deformation gradient tensor here.
  qlocal[4] = v1[0]; qlocal[7] = v1[1]; qlocal[10] = v1[2];
  qlocal[5] = v2[0]; qlocal[8] = v2[1]; qlocal[11] = v2[2];
  qlocal[6] = v3[0]; qlocal[9] = v3[1]; qlocal[12] = v3[2];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  for (int i = 0; i < 14 + (2 * num_species); i++) {
    qglobal[i] = qlocal[i];
  }

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

  // Rotate deformation gradient tensor back to global coordinate frame.
  // NOTE: We are taking the transpose of the deformation gradient tensor here.
  qglobal[4] = v1[0]; qglobal[7] = v1[1]; qglobal[10] = v1[2];
  qglobal[5] = v2[0]; qglobal[8] = v2[1]; qglobal[11] = v2[2];
  qglobal[6] = v3[0]; qglobal[9] = v3[1]; qglobal[12] = v3[2];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  double *rho_ref_s = elasticity_rgfm->rho_ref_s;
  double *T_ref_s = elasticity_rgfm->T_ref_s;
  double *sound_speed_s = elasticity_rgfm->sound_speed_s;
  double *shear_speed_s = elasticity_rgfm->shear_speed_s;
  double *heat_capacity_s = elasticity_rgfm->heat_capacity_s;

  double *alpha_param_s = elasticity_rgfm->alpha_param_s;
  double *beta_param_s = elasticity_rgfm->beta_param_s;
  double *gamma_param_s = elasticity_rgfm->gamma_param_s;

  double sl = gkyl_elasticity_rgfm_max_abs_speed(num_species, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, ql);
  double sr = gkyl_elasticity_rgfm_max_abs_speed(num_species, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  double *fr = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  gkyl_elasticity_rgfm_flux(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, ql, fl);
  gkyl_elasticity_rgfm_flux(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[14 + ((2 * num_species))];
  for (int i = 0; i < 14 + (2 * num_species); i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  gkyl_free(fl);
  gkyl_free(fr);

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  const double *w0 = &waves[0], *w1 = &waves[14 + (2 * num_species)];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 14 + (2 * num_species); i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr,
  const double phil, const double phir, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir,
  const double* waves, const double* s, double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  double *rho_ref_s = elasticity_rgfm->rho_ref_s;
  double *T_ref_s = elasticity_rgfm->T_ref_s;
  double *sound_speed_s = elasticity_rgfm->sound_speed_s;
  double *shear_speed_s = elasticity_rgfm->shear_speed_s;
  double *heat_capacity_s = elasticity_rgfm->heat_capacity_s;

  double *alpha_param_s = elasticity_rgfm->alpha_param_s;
  double *beta_param_s = elasticity_rgfm->beta_param_s;
  double *gamma_param_s = elasticity_rgfm->gamma_param_s;

  double *fl = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  double *fr = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  gkyl_elasticity_rgfm_flux(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, ql, fl);
  gkyl_elasticity_rgfm_flux(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, qr, fr);

  for (int m = 0; m < 14 + (2 * num_species); m++) {
    flux_jump[m] = fr[m] - fl[m];
  }
  
  double amaxl = gkyl_elasticity_rgfm_max_abs_speed(num_species, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, ql);
  double amaxr = gkyl_elasticity_rgfm_max_abs_speed(num_species, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, qr);

  gkyl_free(fr);
  gkyl_free(fl);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  double *T_ref_s = elasticity_rgfm->T_ref_s;
  double *sound_speed_s = elasticity_rgfm->sound_speed_s;
  double *shear_speed_s = elasticity_rgfm->shear_speed_s;
  double *heat_capacity_s = elasticity_rgfm->heat_capacity_s;

  double *alpha_param_s = elasticity_rgfm->alpha_param_s;
  double *beta_param_s = elasticity_rgfm->beta_param_s;
  double *gamma_param_s = elasticity_rgfm->gamma_param_s;

  double *v = gkyl_malloc(sizeof(double[14 + (2 * num_species)]));
  gkyl_elasticity_rgfm_prim_vars(num_species, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q, v);

  for (int i = 0; i < num_species; i++) {
    if (v[14 + num_species + i] < 0.0) {
      gkyl_free(v);
      return false;
    }
  }
  
  double level_set_total = 0.0;
  for (int i = 0; i < num_species - 1; i++) {
    level_set_total += v[14 + i];
  }
  if (level_set_total > 1.0) {
    gkyl_free(v);
    return false;
  }

  if (v[0] < 0.0 || v[13] < 0.0) {
    gkyl_free(v);
    return false;
  }
  else {
    gkyl_free(v);
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  double *rho_ref_s = elasticity_rgfm->rho_ref_s;
  double *T_ref_s = elasticity_rgfm->T_ref_s;
  double *sound_speed_s = elasticity_rgfm->sound_speed_s;
  double *shear_speed_s = elasticity_rgfm->shear_speed_s;
  double *heat_capacity_s = elasticity_rgfm->heat_capacity_s;

  double *alpha_param_s = elasticity_rgfm->alpha_param_s;
  double *beta_param_s = elasticity_rgfm->beta_param_s;
  double *gamma_param_s = elasticity_rgfm->gamma_param_s;

  return gkyl_elasticity_rgfm_max_abs_speed(num_species, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q);
}

static inline void
elasticity_rgfm_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  const struct wv_elasticity_rgfm *elasticity_rgfm = container_of(eqn, struct wv_elasticity_rgfm, eqn);
  int num_species = elasticity_rgfm->num_species;

  for (int i = 0; i < 14 + (2 * num_species); i++) {
    diag[i] = qin[i];
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

  if (inp->rp_type == WV_ELASTICITY_RGFM_RP_LAX) {
    elasticity_rgfm->eqn.num_waves = 2;
    elasticity_rgfm->eqn.waves_func = wave_lax_l;
    elasticity_rgfm->eqn.qfluct_func = qfluct_lax_l;
  }

  elasticity_rgfm->eqn.flux_jump = flux_jump;
  elasticity_rgfm->eqn.check_inv_func = check_inv;
  elasticity_rgfm->eqn.max_speed_func = max_speed;
  elasticity_rgfm->eqn.rotate_to_local_func = rot_to_local;
  elasticity_rgfm->eqn.rotate_to_global_func = rot_to_global;

  elasticity_rgfm->eqn.wall_bc_func = elasticity_rgfm_wall;
  elasticity_rgfm->eqn.no_slip_bc_func = elasticity_rgfm_no_slip;

  elasticity_rgfm->eqn.cons_to_riem = cons_to_riem;
  elasticity_rgfm->eqn.riem_to_cons = riem_to_cons;

  elasticity_rgfm->eqn.cons_to_diag = elasticity_rgfm_cons_to_diag;

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