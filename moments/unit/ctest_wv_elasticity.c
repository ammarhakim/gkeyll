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

  struct gkyl_wv_eqn *elasticity = gkyl_wv_elasticity_new(T_ref, sound_speed, shear_speed, heat_capacity, alpha_param, beta_param, gamma_param, false);

  TEST_CHECK( elasticity->num_equations == 14 );

  double deformation_gradient[3][3];
  deformation_gradient[0][0] = 0.9; deformation_gradient[0][1] = 0.1; deformation_gradient[0][2] = -0.1;
  deformation_gradient[1][0] = 0.05; deformation_gradient[1][1] = 1.1; deformation_gradient[1][2] = -0.05;
  deformation_gradient[2][0] = 0.1; deformation_gradient[2][1] = -0.1; deformation_gradient[2][2] = 0.95;

  double vel_x = 0.1, vel_y = 0.2, vel_z = 0.3;
  double specific_entropy = 0.05;

  double deformation_gradient_det = (deformation_gradient[0][0] * ((deformation_gradient[1][1] * deformation_gradient[2][2]) - (deformation_gradient[2][1] * deformation_gradient[1][2]))) -
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

  double inv_deformation_gradient[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_deformation_gradient[i][j] = (1.0 / deformation_gradient_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * deformation_gradient[i][j]) + deformation_gradient_sq[i][j]);
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

  double internal_energy = ((bulk_modulus / (2.0 * (alpha_param * alpha_param))) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0) * (pow(strain_invariant3, 0.5 * alpha_param) - 1.0)) +
    (heat_capacity * T_ref * pow(strain_invariant3, 0.5 * gamma_param) * (exp(specific_entropy / heat_capacity) - 1.0)) + (0.5 * shear_modulus * pow(strain_invariant3, 0.5 * beta_param) *
    (((strain_invariant1 * strain_invariant1) / 3.0) - strain_invariant2));
  
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
}

TEST_LIST = {
  { "elasticity_basic", test_elasticity_basic },
  { NULL, NULL },
};
