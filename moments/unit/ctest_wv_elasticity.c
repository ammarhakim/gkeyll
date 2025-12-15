#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_elasticity.h>
#include <gkyl_wv_elasticity_priv.h>

void
test_elasticity_basic()
{
  double rho_ref = 8.93;
  double T_ref = 300.0;
  double sound_speed = 4.6;
  double shear_speed = 2.1;
  double heat_capacity = 3.9 * pow(10.0, -4.0);
  double alpha_param = 1.0;
  double beta_param = 3.0;
  double gamma_param = 2.0;

  struct gkyl_wv_eqn *elasticity = gkyl_wv_elasticity_new(rho_ref, T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, false);

  TEST_CHECK( elasticity->num_equations == 14 );
  TEST_CHECK( elasticity->num_waves == 2 );

  double deformation_gradient[3][3];
  deformation_gradient[0][0] = 0.9; deformation_gradient[0][1] = 0.1; deformation_gradient[0][2] = 0.1;
  deformation_gradient[1][0] = 0.05; deformation_gradient[1][1] = 1.1; deformation_gradient[1][2] = 0.05;
  deformation_gradient[2][0] = 0.1; deformation_gradient[2][1] = 0.1; deformation_gradient[2][2] = 0.95;

  double vel_x = 0.1, vel_y = 0.2, vel_z = 0.3;
  double specific_entropy = pow(10.0, -3.0);

  double deformation_gradient_det = (deformation_gradient[0][0] * ((deformation_gradient[1][1] * deformation_gradient[2][2]) -
    (deformation_gradient[2][1] * deformation_gradient[1][2]))) -
    (deformation_gradient[0][1] * ((deformation_gradient[1][0] * deformation_gradient[2][2]) - (deformation_gradient[1][2] * deformation_gradient[2][0]))) +
    (deformation_gradient[0][2] * ((deformation_gradient[1][0] * deformation_gradient[2][1]) - (deformation_gradient[1][1] * deformation_gradient[2][0])));
  double rho = rho_ref / deformation_gradient_det;

  double bulk_modulus = (sound_speed * sound_speed) - ((4.0 / 3.0) * (shear_speed * shear_speed));
  double shear_modulus = shear_speed * shear_speed;

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
  
  double inv_deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient[i][j] = (1.0 / deformation_gradient_det) *
        ((0.5 * ((trace * trace) - sq_trace) * identity_tensor[i][j]) - (trace * deformation_gradient[i][j]) + deformation_gradient_sq[i][j]);
    }
  }

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
  
  double internal_energy = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0) *
    (pow(strain_invariant3, 0.5 * alpha_param) - 1.0)) + (heat_capacity * T_ref * pow(strain_invariant3, 0.5 * gamma_param) * (exp(specific_entropy / heat_capacity) - 1.0)) +
    (0.5 * shear_modulus * pow(strain_invariant3, 0.5 * beta_param) * (((strain_invariant1 * strain_invariant1) / 3.0) - strain_invariant2));
  
  double E_tot = rho * (internal_energy + (0.5 * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z))));
  
  double q[14];
  q[0] = rho;
  q[1] = rho * vel_x;
  q[2] = rho * vel_y;
  q[3] = rho * vel_z;

  q[4] = rho * deformation_gradient[0][0]; q[5] = rho * deformation_gradient[0][1]; q[6] = rho * deformation_gradient[0][2];
  q[7] = rho * deformation_gradient[1][0]; q[8] = rho * deformation_gradient[1][1]; q[9] = rho * deformation_gradient[1][2];
  q[10] = rho * deformation_gradient[2][0]; q[11] = rho * deformation_gradient[2][1]; q[12] = rho * deformation_gradient[2][2];

  q[13] = E_tot;

  double prims[14];
  gkyl_elasticity_prim_vars(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[1], vel_x, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[2], vel_y, 1e-12) );
  TEST_CHECK( gkyl_compare(prims[3], vel_z, 1e-12) );

  TEST_CHECK( gkyl_compare(prims[4], deformation_gradient[0][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[5], deformation_gradient[0][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[6], deformation_gradient[0][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[7], deformation_gradient[1][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[8], deformation_gradient[1][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[9], deformation_gradient[1][2], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[10], deformation_gradient[2][0], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[11], deformation_gradient[2][1], 1e-12) );
  TEST_CHECK( gkyl_compare(prims[12], deformation_gradient[2][2], 1e-12) );

  TEST_CHECK( gkyl_compare(prims[13], specific_entropy, 1e-12) );

  double energy_deriv_invariant1 = (shear_modulus / 3.0) * strain_invariant1 * pow(strain_invariant3, 0.5 * beta_param);
  double energy_deriv_invariant2 = -(0.5 * shear_modulus) * pow(strain_invariant3, 0.5 * beta_param);

  double energy_deriv_invariant3 = (1.0 / (12.0 * alpha_param * strain_invariant3)) * (alpha_param * shear_modulus * beta_param *
    ((strain_invariant1 * strain_invariant1) - (3.0 * strain_invariant2)) * pow(strain_invariant3, 0.5 * beta_param));
  energy_deriv_invariant3 -= (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * pow(strain_invariant3, 0.5 * alpha_param) * bulk_modulus);
  energy_deriv_invariant3 += (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * pow(strain_invariant3, alpha_param) * bulk_modulus);
  energy_deriv_invariant3 += (1.0 / (12.0 * alpha_param * strain_invariant3)) * (6.0 * alpha_param * heat_capacity * (exp(specific_entropy / heat_capacity) - 1.0) *
    gamma_param * pow(strain_invariant3, 0.5 * gamma_param) * T_ref);

  double invariant2_deriv_strain[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      invariant2_deriv_strain[i][j] = (strain_invariant1 * identity_tensor[i][j]) - strain_tensor[i][j];
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
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant1 * identity_tensor[j][k]);
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant2 * invariant2_deriv_strain[j][k]);
        stress_tensor[i][j] -= (2.0 * rho * strain_tensor[i][k]) * (energy_deriv_invariant3 * invariant3_deriv_strain[j][k]);
      }
    }
  }

  double fluxes[3][14] = {
    { rho * vel_x,
      (rho * vel_x * vel_x) - stress_tensor[0][0],
      (rho * vel_y * vel_x) - stress_tensor[1][0],
      (rho * vel_z * vel_x) - stress_tensor[2][0],
      (rho * deformation_gradient[0][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_x),
      (rho * deformation_gradient[0][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_x),
      (rho * deformation_gradient[0][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_x),
      (rho * deformation_gradient[1][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_y),
      (rho * deformation_gradient[1][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_y),
      (rho * deformation_gradient[1][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_y),
      (rho * deformation_gradient[2][0] * vel_x) - (rho * deformation_gradient[0][0] * vel_z),
      (rho * deformation_gradient[2][1] * vel_x) - (rho * deformation_gradient[0][1] * vel_z),
      (rho * deformation_gradient[2][2] * vel_x) - (rho * deformation_gradient[0][2] * vel_z),
      (E_tot * vel_x) - (vel_x * stress_tensor[0][0]) - (vel_y * stress_tensor[1][0]) - (vel_z * stress_tensor[2][0]) },
    { rho * vel_y,
      (rho * vel_x * vel_y) - stress_tensor[0][1],
      (rho * vel_y * vel_y) - stress_tensor[1][1],
      (rho * vel_z * vel_y) - stress_tensor[2][1],
      (rho * deformation_gradient[0][0] * vel_y) - (rho * deformation_gradient[1][0] * vel_x),
      (rho * deformation_gradient[0][1] * vel_y) - (rho * deformation_gradient[1][1] * vel_x),
      (rho * deformation_gradient[0][2] * vel_y) - (rho * deformation_gradient[1][2] * vel_x),
      (rho * deformation_gradient[1][0] * vel_y) - (rho * deformation_gradient[1][0] * vel_y),
      (rho * deformation_gradient[1][1] * vel_y) - (rho * deformation_gradient[1][1] * vel_y),
      (rho * deformation_gradient[1][2] * vel_y) - (rho * deformation_gradient[1][2] * vel_y),
      (rho * deformation_gradient[2][0] * vel_y) - (rho * deformation_gradient[1][0] * vel_z),
      (rho * deformation_gradient[2][1] * vel_y) - (rho * deformation_gradient[1][1] * vel_z),
      (rho * deformation_gradient[2][2] * vel_y) - (rho * deformation_gradient[1][2] * vel_z),
      (E_tot * vel_y) - (vel_x * stress_tensor[0][1]) - (vel_y * stress_tensor[1][1]) - (vel_z * stress_tensor[2][1]) },
    { rho * vel_z,
      (rho * vel_x * vel_z) - stress_tensor[0][2],
      (rho * vel_y * vel_z) - stress_tensor[1][2],
      (rho * vel_z * vel_z) - stress_tensor[2][2],
      (rho * deformation_gradient[0][0] * vel_z) - (rho * deformation_gradient[2][0] * vel_x),
      (rho * deformation_gradient[0][1] * vel_z) - (rho * deformation_gradient[2][1] * vel_x),
      (rho * deformation_gradient[0][2] * vel_z) - (rho * deformation_gradient[2][2] * vel_x),
      (rho * deformation_gradient[1][0] * vel_z) - (rho * deformation_gradient[2][0] * vel_y),
      (rho * deformation_gradient[1][1] * vel_z) - (rho * deformation_gradient[2][1] * vel_y),
      (rho * deformation_gradient[1][2] * vel_z) - (rho * deformation_gradient[2][2] * vel_y),
      (rho * deformation_gradient[2][0] * vel_z) - (rho * deformation_gradient[2][0] * vel_z),
      (rho * deformation_gradient[2][1] * vel_z) - (rho * deformation_gradient[2][1] * vel_z),
      (rho * deformation_gradient[2][2] * vel_z) - (rho * deformation_gradient[2][2] * vel_z),
      (E_tot * vel_z) - (vel_x * stress_tensor[0][2]) - (vel_y * stress_tensor[1][2]) - (vel_z * stress_tensor[2][2]) },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[14], flux_local[14], flux[14];

  for (int d = 0; d < 3; d++) {
    elasticity->rotate_to_local_func(elasticity, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, q_local, flux_local);
    elasticity->rotate_to_global_func(elasticity, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 14; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-12) );
    }
  }

  double q_l[14], q_g[14];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(elasticity, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 14; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[14], q1[14];
    elasticity->cons_to_riem(elasticity, q_local, q_local, w1);
    elasticity->riem_to_cons(elasticity, q_local, w1, q1);

    for (int i = 0; i < 14; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }
}

void
test_elasticity_waves()
{
  double rho_ref = 8.93;
  double T_ref = 300.0;
  double sound_speed = 4.6;
  double shear_speed = 2.1;
  double heat_capacity = 3.9 * pow(10.0, -4.0);
  double alpha_param = 1.0;
  double beta_param = 3.0;
  double gamma_param = 2.0;

  struct gkyl_wv_eqn *elasticity = gkyl_wv_elasticity_new(rho_ref, T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, false);

  double bulk_modulus = (sound_speed * sound_speed) - ((4.0 / 3.0) * (shear_speed * shear_speed));
  double shear_modulus = shear_speed * shear_speed;

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

  double deformation_gradient_l[3][3];
  deformation_gradient_l[0][0] = 0.9; deformation_gradient_l[0][1] = 0.1; deformation_gradient_l[0][2] = 0.1;
  deformation_gradient_l[1][0] = 0.05; deformation_gradient_l[1][1] = 1.1; deformation_gradient_l[1][2] = 0.05;
  deformation_gradient_l[2][0] = 0.1; deformation_gradient_l[2][1] = 0.1; deformation_gradient_l[2][2] = 0.95;

  double vel_x_l = 0.1, vel_y_l = 0.2, vel_z_l = 0.3;
  double specific_entropy_l = pow(10.0, -3.0);

  double deformation_gradient_det_l = (deformation_gradient_l[0][0] * ((deformation_gradient_l[1][1] * deformation_gradient_l[2][2]) -
    (deformation_gradient_l[2][1] * deformation_gradient_l[1][2]))) -
    (deformation_gradient_l[0][1] * ((deformation_gradient_l[1][0] * deformation_gradient_l[2][2]) - (deformation_gradient_l[1][2] * deformation_gradient_l[2][0]))) +
    (deformation_gradient_l[0][2] * ((deformation_gradient_l[1][0] * deformation_gradient_l[2][1]) - (deformation_gradient_l[1][1] * deformation_gradient_l[2][0])));
  double rho_l = rho_ref / deformation_gradient_det_l;

  double trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    trace_l += deformation_gradient_l[i][i];
  }

  double deformation_gradient_sq_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq_l[i][j] += deformation_gradient_l[i][k] * deformation_gradient_l[k][j];
      }
    }
  }

  double sq_trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace_l += deformation_gradient_sq_l[i][i];
  }
  
  double inv_deformation_gradient_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_l[i][j] = (1.0 / deformation_gradient_det_l) *
        ((0.5 * ((trace_l * trace_l) - sq_trace_l) * identity_tensor[i][j]) - (trace_l * deformation_gradient_l[i][j]) + deformation_gradient_sq_l[i][j]);
    }
  }

  double inv_deformation_gradient_transpose_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose_l[i][j] = inv_deformation_gradient_l[j][i];
    }
  }

  double strain_tensor_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_l[i][j] += inv_deformation_gradient_transpose_l[i][k] * inv_deformation_gradient_l[k][j];
      }
    }
  }

  double strain_tensor_sq_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq_l[i][j] += strain_tensor_l[i][k] * strain_tensor_l[k][j];
      }
    }
  }

  double strain_invariant1_l = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1_l += strain_tensor_l[i][i];
  }

  double sq_strain_trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace_l += strain_tensor_sq_l[i][i];
  }

  double strain_invariant2_l = 0.5 * ((strain_invariant1_l * strain_invariant1_l) - sq_strain_trace_l);

  double strain_invariant3_l = (strain_tensor_l[0][0] * ((strain_tensor_l[1][1] * strain_tensor_l[2][2]) - (strain_tensor_l[2][1] * strain_tensor_l[1][2]))) -
    (strain_tensor_l[0][1] * ((strain_tensor_l[1][0] * strain_tensor_l[2][2]) - (strain_tensor_l[1][2] * strain_tensor_l[2][0]))) +
    (strain_tensor_l[0][2] * ((strain_tensor_l[1][0] * strain_tensor_l[2][1]) - (strain_tensor_l[1][1] * strain_tensor_l[2][0])));
  
  double internal_energy_l = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3_l, 0.5 * alpha_param) - 1.0) *
    (pow(strain_invariant3_l, 0.5 * alpha_param) - 1.0)) + (heat_capacity * T_ref * pow(strain_invariant3_l, 0.5 * gamma_param) * (exp(specific_entropy_l / heat_capacity) - 1.0)) +
    (0.5 * shear_modulus * pow(strain_invariant3_l, 0.5 * beta_param) * (((strain_invariant1_l * strain_invariant1_l) / 3.0) - strain_invariant2_l));
  
  double E_tot_l = rho_l * (internal_energy_l + (0.5 * ((vel_x_l * vel_x_l) + (vel_y_l * vel_y_l) + (vel_z_l * vel_z_l))));
  
  double ql[14];
  ql[0] = rho_l;
  ql[1] = rho_l * vel_x_l;
  ql[2] = rho_l * vel_y_l;
  ql[3] = rho_l * vel_z_l;

  ql[4] = rho_l * deformation_gradient_l[0][0]; ql[5] = rho_l * deformation_gradient_l[0][1]; ql[6] = rho_l * deformation_gradient_l[0][2];
  ql[7] = rho_l * deformation_gradient_l[1][0]; ql[8] = rho_l * deformation_gradient_l[1][1]; ql[9] = rho_l * deformation_gradient_l[1][2];
  ql[10] = rho_l * deformation_gradient_l[2][0]; ql[11] = rho_l * deformation_gradient_l[2][1]; ql[12] = rho_l * deformation_gradient_l[2][2];

  ql[13] = E_tot_l;

  double deformation_gradient_r[3][3];
  deformation_gradient_r[0][0] = 0.8; deformation_gradient_r[0][1] = 0.05; deformation_gradient_r[0][2] = 0.2;
  deformation_gradient_r[1][0] = 0.02; deformation_gradient_r[1][1] = 1.2; deformation_gradient_r[1][2] = 0.02;
  deformation_gradient_r[2][0] = 0.2; deformation_gradient_r[2][1] = 0.05; deformation_gradient_r[2][2] = 0.9;

  double vel_x_r = 1.0, vel_y_r = 2.0, vel_z_r = 3.0;
  double specific_entropy_r = 2.0 * pow(10.0, -3.0);

  double deformation_gradient_det_r = (deformation_gradient_r[0][0] * ((deformation_gradient_r[1][1] * deformation_gradient_r[2][2]) -
    (deformation_gradient_r[2][1] * deformation_gradient_r[1][2]))) -
    (deformation_gradient_r[0][1] * ((deformation_gradient_r[1][0] * deformation_gradient_r[2][2]) - (deformation_gradient_r[1][2] * deformation_gradient_r[2][0]))) +
    (deformation_gradient_r[0][2] * ((deformation_gradient_r[1][0] * deformation_gradient_r[2][1]) - (deformation_gradient_r[1][1] * deformation_gradient_r[2][0])));
  double rho_r = rho_ref / deformation_gradient_det_r;

  double trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    trace_r += deformation_gradient_r[i][i];
  }

  double deformation_gradient_sq_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq_r[i][j] += deformation_gradient_r[i][k] * deformation_gradient_r[k][j];
      }
    }
  }

  double sq_trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace_r += deformation_gradient_sq_r[i][i];
  }
  
  double inv_deformation_gradient_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_r[i][j] = (1.0 / deformation_gradient_det_r) *
        ((0.5 * ((trace_r * trace_r) - sq_trace_r) * identity_tensor[i][j]) - (trace_r * deformation_gradient_r[i][j]) + deformation_gradient_sq_r[i][j]);
    }
  }

  double inv_deformation_gradient_transpose_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose_r[i][j] = inv_deformation_gradient_r[j][i];
    }
  }

  double strain_tensor_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_r[i][j] += inv_deformation_gradient_transpose_r[i][k] * inv_deformation_gradient_r[k][j];
      }
    }
  }

  double strain_tensor_sq_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq_r[i][j] += strain_tensor_r[i][k] * strain_tensor_r[k][j];
      }
    }
  }

  double strain_invariant1_r = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1_r += strain_tensor_r[i][i];
  }

  double sq_strain_trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace_r += strain_tensor_sq_r[i][i];
  }

  double strain_invariant2_r = 0.5 * ((strain_invariant1_r * strain_invariant1_r) - sq_strain_trace_r);

  double strain_invariant3_r = (strain_tensor_r[0][0] * ((strain_tensor_r[1][1] * strain_tensor_r[2][2]) - (strain_tensor_r[2][1] * strain_tensor_r[1][2]))) -
    (strain_tensor_r[0][1] * ((strain_tensor_r[1][0] * strain_tensor_r[2][2]) - (strain_tensor_r[1][2] * strain_tensor_r[2][0]))) +
    (strain_tensor_r[0][2] * ((strain_tensor_r[1][0] * strain_tensor_r[2][1]) - (strain_tensor_r[1][1] * strain_tensor_r[2][0])));
  
  double internal_energy_r = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3_r, 0.5 * alpha_param) - 1.0) *
    (pow(strain_invariant3_r, 0.5 * alpha_param) - 1.0)) + (heat_capacity * T_ref * pow(strain_invariant3_r, 0.5 * gamma_param) * (exp(specific_entropy_r / heat_capacity) - 1.0)) +
    (0.5 * shear_modulus * pow(strain_invariant3_r, 0.5 * beta_param) * (((strain_invariant1_r * strain_invariant1_r) / 3.0) - strain_invariant2_r));
  
  double E_tot_r = rho_r * (internal_energy_r + (0.5 * ((vel_x_r * vel_x_r) + (vel_y_r * vel_y_r) + (vel_z_r * vel_z_r))));
  
  double qr[14];
  qr[0] = rho_r;
  qr[1] = rho_r * vel_x_r;
  qr[2] = rho_r * vel_y_r;
  qr[3] = rho_r * vel_z_r;

  qr[4] = rho_r * deformation_gradient_r[0][0]; qr[5] = rho_r * deformation_gradient_r[0][1]; qr[6] = rho_r * deformation_gradient_r[0][2];
  qr[7] = rho_r * deformation_gradient_r[1][0]; qr[8] = rho_r * deformation_gradient_r[1][1]; qr[9] = rho_r * deformation_gradient_r[1][2];
  qr[10] = rho_r * deformation_gradient_r[2][0]; qr[11] = rho_r * deformation_gradient_r[2][1]; qr[12] = rho_r * deformation_gradient_r[2][2];

  qr[13] = E_tot_r;

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 14], waves_local[2 * 14];

    double ql_local[14], qr_local[14];
    gkyl_wv_eqn_rotate_to_local(elasticity, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(elasticity, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[14];
    for (int i = 0; i < 14; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(elasticity, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[14], amdq_local[14];
    gkyl_wv_eqn_qfluct(elasticity, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], &waves_local[i * 14], &waves[i * 14]);
    }

    double apdq[14], amdq[14];
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[14], fr_local[14];
    gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql_local, fl_local);
    gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr_local, fr_local);

    double fl[14], fr[14];
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 14; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }
}

void
test_elasticity_waves_2()
{
  double rho_ref = 8.93;
  double T_ref = 300.0;
  double sound_speed = 4.6;
  double shear_speed = 2.1;
  double heat_capacity = 3.9 * pow(10.0, -4.0);
  double alpha_param = 1.0;
  double beta_param = 3.0;
  double gamma_param = 2.0;

  struct gkyl_wv_eqn *elasticity = gkyl_wv_elasticity_new(rho_ref, T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, false);

  double bulk_modulus = (sound_speed * sound_speed) - ((4.0 / 3.0) * (shear_speed * shear_speed));
  double shear_modulus = shear_speed * shear_speed;

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

  double deformation_gradient_l[3][3];
  deformation_gradient_l[0][0] = 0.95; deformation_gradient_l[0][1] = 0.05; deformation_gradient_l[0][2] = 0.05;
  deformation_gradient_l[1][0] = 0.02; deformation_gradient_l[1][1] = 1.05; deformation_gradient_l[1][2] = 0.02;
  deformation_gradient_l[2][0] = 0.05; deformation_gradient_l[2][1] = 0.05; deformation_gradient_l[2][2] = 0.98;

  double vel_x_l = 1.0, vel_y_l = 2.0, vel_z_l = 3.0;
  double specific_entropy_l = 0.5 * pow(10.0, -3.0);

  double deformation_gradient_det_l = (deformation_gradient_l[0][0] * ((deformation_gradient_l[1][1] * deformation_gradient_l[2][2]) -
    (deformation_gradient_l[2][1] * deformation_gradient_l[1][2]))) -
    (deformation_gradient_l[0][1] * ((deformation_gradient_l[1][0] * deformation_gradient_l[2][2]) - (deformation_gradient_l[1][2] * deformation_gradient_l[2][0]))) +
    (deformation_gradient_l[0][2] * ((deformation_gradient_l[1][0] * deformation_gradient_l[2][1]) - (deformation_gradient_l[1][1] * deformation_gradient_l[2][0])));
  double rho_l = rho_ref / deformation_gradient_det_l;

  double trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    trace_l += deformation_gradient_l[i][i];
  }

  double deformation_gradient_sq_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq_l[i][j] += deformation_gradient_l[i][k] * deformation_gradient_l[k][j];
      }
    }
  }

  double sq_trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace_l += deformation_gradient_sq_l[i][i];
  }
  
  double inv_deformation_gradient_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_l[i][j] = (1.0 / deformation_gradient_det_l) *
        ((0.5 * ((trace_l * trace_l) - sq_trace_l) * identity_tensor[i][j]) - (trace_l * deformation_gradient_l[i][j]) + deformation_gradient_sq_l[i][j]);
    }
  }

  double inv_deformation_gradient_transpose_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose_l[i][j] = inv_deformation_gradient_l[j][i];
    }
  }

  double strain_tensor_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_l[i][j] += inv_deformation_gradient_transpose_l[i][k] * inv_deformation_gradient_l[k][j];
      }
    }
  }

  double strain_tensor_sq_l[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq_l[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq_l[i][j] += strain_tensor_l[i][k] * strain_tensor_l[k][j];
      }
    }
  }

  double strain_invariant1_l = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1_l += strain_tensor_l[i][i];
  }

  double sq_strain_trace_l = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace_l += strain_tensor_sq_l[i][i];
  }

  double strain_invariant2_l = 0.5 * ((strain_invariant1_l * strain_invariant1_l) - sq_strain_trace_l);

  double strain_invariant3_l = (strain_tensor_l[0][0] * ((strain_tensor_l[1][1] * strain_tensor_l[2][2]) - (strain_tensor_l[2][1] * strain_tensor_l[1][2]))) -
    (strain_tensor_l[0][1] * ((strain_tensor_l[1][0] * strain_tensor_l[2][2]) - (strain_tensor_l[1][2] * strain_tensor_l[2][0]))) +
    (strain_tensor_l[0][2] * ((strain_tensor_l[1][0] * strain_tensor_l[2][1]) - (strain_tensor_l[1][1] * strain_tensor_l[2][0])));
  
  double internal_energy_l = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3_l, 0.5 * alpha_param) - 1.0) *
    (pow(strain_invariant3_l, 0.5 * alpha_param) - 1.0)) + (heat_capacity * T_ref * pow(strain_invariant3_l, 0.5 * gamma_param) * (exp(specific_entropy_l / heat_capacity) - 1.0)) +
    (0.5 * shear_modulus * pow(strain_invariant3_l, 0.5 * beta_param) * (((strain_invariant1_l * strain_invariant1_l) / 3.0) - strain_invariant2_l));
  
  double E_tot_l = rho_l * (internal_energy_l + (0.5 * ((vel_x_l * vel_x_l) + (vel_y_l * vel_y_l) + (vel_z_l * vel_z_l))));
  
  double ql[14];
  ql[0] = rho_l;
  ql[1] = rho_l * vel_x_l;
  ql[2] = rho_l * vel_y_l;
  ql[3] = rho_l * vel_z_l;

  ql[4] = rho_l * deformation_gradient_l[0][0]; ql[5] = rho_l * deformation_gradient_l[0][1]; ql[6] = rho_l * deformation_gradient_l[0][2];
  ql[7] = rho_l * deformation_gradient_l[1][0]; ql[8] = rho_l * deformation_gradient_l[1][1]; ql[9] = rho_l * deformation_gradient_l[1][2];
  ql[10] = rho_l * deformation_gradient_l[2][0]; ql[11] = rho_l * deformation_gradient_l[2][1]; ql[12] = rho_l * deformation_gradient_l[2][2];

  ql[13] = E_tot_l;

  double deformation_gradient_r[3][3];
  deformation_gradient_r[0][0] = 0.7; deformation_gradient_r[0][1] = 0.0; deformation_gradient_r[0][2] = 0.4;
  deformation_gradient_r[1][0] = 0.0; deformation_gradient_r[1][1] = 1.4; deformation_gradient_r[1][2] = 0.0;
  deformation_gradient_r[2][0] = 0.3; deformation_gradient_r[2][1] = 0.0; deformation_gradient_r[2][2] = 0.8;

  double vel_x_r = 0.01, vel_y_r = 0.02, vel_z_r = 0.03;
  double specific_entropy_r = 5.0 * pow(10.0, -3.0);

  double deformation_gradient_det_r = (deformation_gradient_r[0][0] * ((deformation_gradient_r[1][1] * deformation_gradient_r[2][2]) -
    (deformation_gradient_r[2][1] * deformation_gradient_r[1][2]))) -
    (deformation_gradient_r[0][1] * ((deformation_gradient_r[1][0] * deformation_gradient_r[2][2]) - (deformation_gradient_r[1][2] * deformation_gradient_r[2][0]))) +
    (deformation_gradient_r[0][2] * ((deformation_gradient_r[1][0] * deformation_gradient_r[2][1]) - (deformation_gradient_r[1][1] * deformation_gradient_r[2][0])));
  double rho_r = rho_ref / deformation_gradient_det_r;

  double trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    trace_r += deformation_gradient_r[i][i];
  }

  double deformation_gradient_sq_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_sq_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        deformation_gradient_sq_r[i][j] += deformation_gradient_r[i][k] * deformation_gradient_r[k][j];
      }
    }
  }

  double sq_trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace_r += deformation_gradient_sq_r[i][i];
  }
  
  double inv_deformation_gradient_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_r[i][j] = (1.0 / deformation_gradient_det_r) *
        ((0.5 * ((trace_r * trace_r) - sq_trace_r) * identity_tensor[i][j]) - (trace_r * deformation_gradient_r[i][j]) + deformation_gradient_sq_r[i][j]);
    }
  }

  double inv_deformation_gradient_transpose_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient_transpose_r[i][j] = inv_deformation_gradient_r[j][i];
    }
  }

  double strain_tensor_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_r[i][j] += inv_deformation_gradient_transpose_r[i][k] * inv_deformation_gradient_r[k][j];
      }
    }
  }

  double strain_tensor_sq_r[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_tensor_sq_r[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        strain_tensor_sq_r[i][j] += strain_tensor_r[i][k] * strain_tensor_r[k][j];
      }
    }
  }

  double strain_invariant1_r = 0.0;
  for (int i = 0; i < 3; i++) {
    strain_invariant1_r += strain_tensor_r[i][i];
  }

  double sq_strain_trace_r = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_strain_trace_r += strain_tensor_sq_r[i][i];
  }

  double strain_invariant2_r = 0.5 * ((strain_invariant1_r * strain_invariant1_r) - sq_strain_trace_r);

  double strain_invariant3_r = (strain_tensor_r[0][0] * ((strain_tensor_r[1][1] * strain_tensor_r[2][2]) - (strain_tensor_r[2][1] * strain_tensor_r[1][2]))) -
    (strain_tensor_r[0][1] * ((strain_tensor_r[1][0] * strain_tensor_r[2][2]) - (strain_tensor_r[1][2] * strain_tensor_r[2][0]))) +
    (strain_tensor_r[0][2] * ((strain_tensor_r[1][0] * strain_tensor_r[2][1]) - (strain_tensor_r[1][1] * strain_tensor_r[2][0])));
  
  double internal_energy_r = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3_r, 0.5 * alpha_param) - 1.0) *
    (pow(strain_invariant3_r, 0.5 * alpha_param) - 1.0)) + (heat_capacity * T_ref * pow(strain_invariant3_r, 0.5 * gamma_param) * (exp(specific_entropy_r / heat_capacity) - 1.0)) +
    (0.5 * shear_modulus * pow(strain_invariant3_r, 0.5 * beta_param) * (((strain_invariant1_r * strain_invariant1_r) / 3.0) - strain_invariant2_r));
  
  double E_tot_r = rho_r * (internal_energy_r + (0.5 * ((vel_x_r * vel_x_r) + (vel_y_r * vel_y_r) + (vel_z_r * vel_z_r))));
  
  double qr[14];
  qr[0] = rho_r;
  qr[1] = rho_r * vel_x_r;
  qr[2] = rho_r * vel_y_r;
  qr[3] = rho_r * vel_z_r;

  qr[4] = rho_r * deformation_gradient_r[0][0]; qr[5] = rho_r * deformation_gradient_r[0][1]; qr[6] = rho_r * deformation_gradient_r[0][2];
  qr[7] = rho_r * deformation_gradient_r[1][0]; qr[8] = rho_r * deformation_gradient_r[1][1]; qr[9] = rho_r * deformation_gradient_r[1][2];
  qr[10] = rho_r * deformation_gradient_r[2][0]; qr[11] = rho_r * deformation_gradient_r[2][1]; qr[12] = rho_r * deformation_gradient_r[2][2];

  qr[13] = E_tot_r;

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 14], waves_local[2 * 14];

    double ql_local[14], qr_local[14];
    gkyl_wv_eqn_rotate_to_local(elasticity, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(elasticity, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[14];
    for (int i = 0; i < 14; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(elasticity, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[14], amdq_local[14];
    gkyl_wv_eqn_qfluct(elasticity, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], &waves_local[i * 14], &waves[i * 14]);
    }

    double apdq[14], amdq[14];
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[14], fr_local[14];
    gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, ql_local, fl_local);
    gkyl_elasticity_flux(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, qr_local, fr_local);

    double fl[14], fr[14];
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(elasticity, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 14; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-11) );
    }
  }
}

TEST_LIST = {
  { "elasticity_basic", test_elasticity_basic },
  { "elasticity_waves", test_elasticity_waves },
  { "elasticity_waves_2", test_elasticity_waves_2 },
  { NULL, NULL },
};
