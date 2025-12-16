#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_elasticity_rgfm.h>
#include <gkyl_wv_elasticity_rgfm_priv.h>

void
test_elasticity_rgfm_twospecies_basic()
{
  double rho_ref1 = 8.93;
  double T_ref1 = 300.0;
  double sound_speed1 = 4.6;
  double shear_speed1 = 2.1;
  double heat_capacity1 = 3.9 * pow(10.0, -4.0);
  double alpha_param1 = 1.0;
  double beta_param1 = 3.0;
  double gamma_param1 = 2.0;

  double rho_ref2 = 2.71;
  double T_ref2 = 300.0;
  double sound_speed2 = 6.22;
  double shear_speed2 = 3.16;
  double heat_capacity2 = 9.0 * pow(10.0, -4.0);
  double alpha_param2 = 1.0;
  double beta_param2 = 3.577;
  double gamma_param2 = 2.088;

  double *rho_ref_s = gkyl_malloc(sizeof(double[2]));
  rho_ref_s[0] = rho_ref1;
  rho_ref_s[1] = rho_ref2;

  double *T_ref_s = gkyl_malloc(sizeof(double[2]));
  T_ref_s[0] = T_ref1;
  T_ref_s[1] = T_ref2;

  double *sound_speed_s = gkyl_malloc(sizeof(double[2]));
  sound_speed_s[0] = sound_speed1;
  sound_speed_s[1] = sound_speed2;
  
  double *shear_speed_s = gkyl_malloc(sizeof(double[2]));
  shear_speed_s[0] = shear_speed1;
  shear_speed_s[1] = shear_speed2;
  
  double *heat_capacity_s = gkyl_malloc(sizeof(double[2]));
  heat_capacity_s[0] = heat_capacity1;
  heat_capacity_s[1] = heat_capacity2;

  double *alpha_param_s = gkyl_malloc(sizeof(double[2]));
  alpha_param_s[0] = alpha_param1;
  alpha_param_s[1] = alpha_param2;

  double *beta_param_s = gkyl_malloc(sizeof(double[2]));
  beta_param_s[0] = beta_param1;
  beta_param_s[1] = beta_param2;

  double *gamma_param_s = gkyl_malloc(sizeof(double[2]));
  gamma_param_s[0] = gamma_param1;
  gamma_param_s[1] = gamma_param2;

  struct gkyl_wv_eqn *elasticity_rgfm = gkyl_wv_elasticity_rgfm_new(2, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, 0, 0.0, false);

  TEST_CHECK( elasticity_rgfm->num_equations == 18 );
  //TEST_CHECK( elasticity_rgfm->num_waves == 2 );

  double deformation_gradient_total[3][3];
  deformation_gradient_total[0][0] = 0.9; deformation_gradient_total[0][1] = 0.1; deformation_gradient_total[0][2] = 0.1;
  deformation_gradient_total[1][0] = 0.05; deformation_gradient_total[1][1] = 1.1; deformation_gradient_total[1][2] = 0.05;
  deformation_gradient_total[2][0] = 0.1; deformation_gradient_total[2][1] = 0.1; deformation_gradient_total[2][2] = 0.95;

  double phi1 = 0.75;
  double vel_x_total = 0.1, vel_y_total = 0.2, vel_z_total = 0.3;
  double specific_entropy_total = pow(10.0, -3.0);

  double deformation_gradient_total_det = (deformation_gradient_total[0][0] * ((deformation_gradient_total[1][1] * deformation_gradient_total[2][2]) -
    (deformation_gradient_total[2][1] * deformation_gradient_total[1][2]))) -
    (deformation_gradient_total[0][1] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][2]) - (deformation_gradient_total[1][2] * deformation_gradient_total[2][0]))) +
    (deformation_gradient_total[0][2] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][1]) - (deformation_gradient_total[1][1] * deformation_gradient_total[2][0])));
  double rho1 = rho_ref1 / deformation_gradient_total_det;
  double rho2 = rho_ref2 / deformation_gradient_total_det;

  double bulk_modulus1 = (sound_speed1 * sound_speed1) - ((4.0 / 3.0) * (shear_speed1 * shear_speed1));
  double bulk_modulus2 = (sound_speed2 * sound_speed2) - ((4.0 / 3.0) * (shear_speed2 * shear_speed2));
  double shear_modulus1 = shear_speed1 * shear_speed1;
  double shear_modulus2 = shear_speed2 * shear_speed2;

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
  
  double inv_deformation_gradient_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_total[i][j] = (1.0 / deformation_gradient_total_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient_total[i][j]) + deformation_gradient_total_sq[i][j]);
    }
  }

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

  double internal_energy1 = ((bulk_modulus1 / (2.0 * (alpha_param1 * alpha_param1))) * (pow(strain_invariant3_total, 0.5 * alpha_param1) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param1) - 1.0)) + (heat_capacity1 * T_ref1 * pow(strain_invariant3_total, 0.5 * gamma_param1) * (exp(specific_entropy_total / heat_capacity1) - 1.0)) +
    (0.5 * shear_modulus1 * pow(strain_invariant3_total, 0.5 * beta_param1) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));
  double internal_energy2 = ((bulk_modulus2 / (2.0 * (alpha_param2 * alpha_param2))) * (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0)) + (heat_capacity2 * T_ref2 * pow(strain_invariant3_total, 0.5 * gamma_param2) * (exp(specific_entropy_total / heat_capacity2) - 1.0)) +
    (0.5 * shear_modulus2 * pow(strain_invariant3_total, 0.5 * beta_param2) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));

  double rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2);
  double E1 = rho1 * (internal_energy1 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E2 = rho2 * (internal_energy2 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E_total = (phi1 * E1) + ((1.0 - phi1) * E2);

  double q[18];
  q[0] = rho_total;
  q[1] = rho_total * vel_x_total;
  q[2] = rho_total * vel_y_total;
  q[3] = rho_total * vel_z_total;

  q[4] = rho_total * deformation_gradient_total[0][0]; q[5] = rho_total * deformation_gradient_total[0][1]; q[6] = rho_total * deformation_gradient_total[0][2];
  q[7] = rho_total * deformation_gradient_total[1][0]; q[8] = rho_total * deformation_gradient_total[1][1]; q[9] = rho_total * deformation_gradient_total[1][2];
  q[10] = rho_total * deformation_gradient_total[2][0]; q[11] = rho_total * deformation_gradient_total[2][1]; q[12] = rho_total * deformation_gradient_total[2][2];

  q[13] = E_total;

  q[14] = rho_total * phi1;
  q[15] = phi1 * rho1;
  q[16] = (1.0 - phi1) * rho2;
  q[17] = 0.0;

  double prims[18];
  gkyl_elasticity_rgfm_prim_vars(2, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[1], vel_x_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[2], vel_y_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[3], vel_z_total, 1e-12) );

  TEST_CHECK( gkyl_compare(prims[4], deformation_gradient_total[0][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[5], deformation_gradient_total[0][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[6], deformation_gradient_total[0][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[7], deformation_gradient_total[1][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[8], deformation_gradient_total[1][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[9], deformation_gradient_total[1][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[10], deformation_gradient_total[2][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[11], deformation_gradient_total[2][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[12], deformation_gradient_total[2][2], 1e-12) );

  // For now, we check only that the reconstructed interface entropy is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical boundary conditions into the system.
  TEST_CHECK( gkyl_compare(prims[13], specific_entropy_total, 1e-2) );

  TEST_CHECK( gkyl_compare(prims[14], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[15], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[16], rho2, 1e-16) );
}

void
test_elasticity_rgfm_threespecies_basic()
{
  double rho_ref1 = 8.93;
  double T_ref1 = 300.0;
  double sound_speed1 = 4.6;
  double shear_speed1 = 2.1;
  double heat_capacity1 = 3.9 * pow(10.0, -4.0);
  double alpha_param1 = 1.0;
  double beta_param1 = 3.0;
  double gamma_param1 = 2.0;

  double rho_ref2 = 2.71;
  double T_ref2 = 300.0;
  double sound_speed2 = 6.22;
  double shear_speed2 = 3.16;
  double heat_capacity2 = 9.0 * pow(10.0, -4.0);
  double alpha_param2 = 1.0;
  double beta_param2 = 3.577;
  double gamma_param2 = 2.088;

  double rho_ref3 = 8.03;
  double T_ref3 = 300.0;
  double sound_speed3 = 5.68;
  double shear_speed3 = 3.1;
  double heat_capacity3 = 5.0 * pow(10.0, -4.0);
  double alpha_param3 = 0.596;
  double beta_param3 = 2.437;
  double gamma_param3 = 1.563;

  double *rho_ref_s = gkyl_malloc(sizeof(double[3]));
  rho_ref_s[0] = rho_ref1;
  rho_ref_s[1] = rho_ref2;
  rho_ref_s[2] = rho_ref3;

  double *T_ref_s = gkyl_malloc(sizeof(double[3]));
  T_ref_s[0] = T_ref1;
  T_ref_s[1] = T_ref2;
  T_ref_s[2] = T_ref3;

  double *sound_speed_s = gkyl_malloc(sizeof(double[3]));
  sound_speed_s[0] = sound_speed1;
  sound_speed_s[1] = sound_speed2;
  sound_speed_s[2] = sound_speed3;
  
  double *shear_speed_s = gkyl_malloc(sizeof(double[3]));
  shear_speed_s[0] = shear_speed1;
  shear_speed_s[1] = shear_speed2;
  shear_speed_s[2] = shear_speed3;
  
  double *heat_capacity_s = gkyl_malloc(sizeof(double[3]));
  heat_capacity_s[0] = heat_capacity1;
  heat_capacity_s[1] = heat_capacity2;
  heat_capacity_s[2] = heat_capacity3;

  double *alpha_param_s = gkyl_malloc(sizeof(double[3]));
  alpha_param_s[0] = alpha_param1;
  alpha_param_s[1] = alpha_param2;
  alpha_param_s[2] = alpha_param3;

  double *beta_param_s = gkyl_malloc(sizeof(double[3]));
  beta_param_s[0] = beta_param1;
  beta_param_s[1] = beta_param2;
  beta_param_s[2] = beta_param3;

  double *gamma_param_s = gkyl_malloc(sizeof(double[3]));
  gamma_param_s[0] = gamma_param1;
  gamma_param_s[1] = gamma_param2;
  gamma_param_s[2] = gamma_param3;

  struct gkyl_wv_eqn *elasticity_rgfm = gkyl_wv_elasticity_rgfm_new(3, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, 0, 0.0, false);

  TEST_CHECK( elasticity_rgfm->num_equations == 20 );
  //TEST_CHECK( elasticity_rgfm->num_waves == 2 );

  double deformation_gradient_total[3][3];
  deformation_gradient_total[0][0] = 0.9; deformation_gradient_total[0][1] = 0.1; deformation_gradient_total[0][2] = 0.1;
  deformation_gradient_total[1][0] = 0.05; deformation_gradient_total[1][1] = 1.1; deformation_gradient_total[1][2] = 0.05;
  deformation_gradient_total[2][0] = 0.1; deformation_gradient_total[2][1] = 0.1; deformation_gradient_total[2][2] = 0.95;

  double phi1 = 0.5, phi2 = 0.3;
  double vel_x_total = 0.1, vel_y_total = 0.2, vel_z_total = 0.3;
  double specific_entropy_total = pow(10.0, -3.0);

  double deformation_gradient_total_det = (deformation_gradient_total[0][0] * ((deformation_gradient_total[1][1] * deformation_gradient_total[2][2]) -
    (deformation_gradient_total[2][1] * deformation_gradient_total[1][2]))) -
    (deformation_gradient_total[0][1] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][2]) - (deformation_gradient_total[1][2] * deformation_gradient_total[2][0]))) +
    (deformation_gradient_total[0][2] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][1]) - (deformation_gradient_total[1][1] * deformation_gradient_total[2][0])));
  double rho1 = rho_ref1 / deformation_gradient_total_det;
  double rho2 = rho_ref2 / deformation_gradient_total_det;
  double rho3 = rho_ref3 / deformation_gradient_total_det;

  double bulk_modulus1 = (sound_speed1 * sound_speed1) - ((4.0 / 3.0) * (shear_speed1 * shear_speed1));
  double bulk_modulus2 = (sound_speed2 * sound_speed2) - ((4.0 / 3.0) * (shear_speed2 * shear_speed2));
  double bulk_modulus3 = (sound_speed3 * sound_speed3) - ((4.0 / 3.0) * (shear_speed3 * shear_speed3));
  double shear_modulus1 = shear_speed1 * shear_speed1;
  double shear_modulus2 = shear_speed2 * shear_speed2;
  double shear_modulus3 = shear_speed3 * shear_speed3;

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
  
  double inv_deformation_gradient_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_total[i][j] = (1.0 / deformation_gradient_total_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient_total[i][j]) + deformation_gradient_total_sq[i][j]);
    }
  }

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

  double internal_energy1 = ((bulk_modulus1 / (2.0 * (alpha_param1 * alpha_param1))) * (pow(strain_invariant3_total, 0.5 * alpha_param1) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param1) - 1.0)) + (heat_capacity1 * T_ref1 * pow(strain_invariant3_total, 0.5 * gamma_param1) * (exp(specific_entropy_total / heat_capacity1) - 1.0)) +
    (0.5 * shear_modulus1 * pow(strain_invariant3_total, 0.5 * beta_param1) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));
  double internal_energy2 = ((bulk_modulus2 / (2.0 * (alpha_param2 * alpha_param2))) * (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0)) + (heat_capacity2 * T_ref2 * pow(strain_invariant3_total, 0.5 * gamma_param2) * (exp(specific_entropy_total / heat_capacity2) - 1.0)) +
    (0.5 * shear_modulus2 * pow(strain_invariant3_total, 0.5 * beta_param2) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));
  double internal_energy3 = ((bulk_modulus3 / (2.0 * (alpha_param3 * alpha_param3))) * (pow(strain_invariant3_total, 0.5 * alpha_param3) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param3) - 1.0)) + (heat_capacity3 * T_ref3 * pow(strain_invariant3_total, 0.5 * gamma_param3) * (exp(specific_entropy_total / heat_capacity3) - 1.0)) +
    (0.5 * shear_modulus3 * pow(strain_invariant3_total, 0.5 * beta_param3) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));

  double rho_total = (phi1 * rho1) + (phi2 * rho2) + ((1.0 - (phi1 + phi2)) * rho3);
  double E1 = rho1 * (internal_energy1 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E2 = rho2 * (internal_energy2 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E3 = rho3 * (internal_energy3 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E_total = (phi1 * E1) + (phi2 * E2) + ((1.0 - (phi1 + phi2)) * E3);
  
  double q[20];
  q[0] = rho_total;
  q[1] = rho_total * vel_x_total;
  q[2] = rho_total * vel_y_total;
  q[3] = rho_total * vel_z_total;

  q[4] = rho_total * deformation_gradient_total[0][0]; q[5] = rho_total * deformation_gradient_total[0][1]; q[6] = rho_total * deformation_gradient_total[0][2];
  q[7] = rho_total * deformation_gradient_total[1][0]; q[8] = rho_total * deformation_gradient_total[1][1]; q[9] = rho_total * deformation_gradient_total[1][2];
  q[10] = rho_total * deformation_gradient_total[2][0]; q[11] = rho_total * deformation_gradient_total[2][1]; q[12] = rho_total * deformation_gradient_total[2][2];

  q[13] = E_total;

  q[14] = rho_total * phi1;
  q[15] = rho_total * phi2;
  q[16] = phi1 * rho1;
  q[17] = phi2 * rho2;
  q[18] = (1.0 - (phi1 + phi2)) * rho3;
  q[19] = 0.0;
  
  double prims[20];
  gkyl_elasticity_rgfm_prim_vars(3, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[1], vel_x_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[2], vel_y_total, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[3], vel_z_total, 1e-12) );

  TEST_CHECK( gkyl_compare(prims[4], deformation_gradient_total[0][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[5], deformation_gradient_total[0][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[6], deformation_gradient_total[0][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[7], deformation_gradient_total[1][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[8], deformation_gradient_total[1][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[9], deformation_gradient_total[1][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[10], deformation_gradient_total[2][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[11], deformation_gradient_total[2][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[12], deformation_gradient_total[2][2], 1e-12) );

  // For now, we check only that the reconstructed interface entropy is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical boundary conditions into the system.
  TEST_CHECK( gkyl_compare(prims[13], specific_entropy_total, 1e-2) );

  TEST_CHECK( gkyl_compare(prims[14], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[15], phi2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[16], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[17], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[18], rho3, 1e-16) );
}

TEST_LIST = {
  { "elasticity_rgfm_twospecies_basic", test_elasticity_rgfm_twospecies_basic },
  { "elasticity_rgfm_threespecies_basic", test_elasticity_rgfm_threespecies_basic },
  { NULL, NULL },
};