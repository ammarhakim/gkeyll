#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_elasticity_rgfm.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct elasticity_rgfm_shock_ctx
{
  // Physical constants (using normalized code units).
  double rho_ref1; // Reference density (unstressed configuration) of aluminum.
  double T_ref1; // Reference temperature (unstressed configuration) of aluminum.
  double sound_speed1; // Sound speed (speed of pressure waves) of aluminum.
  double shear_speed1; // Shear speed (speed of shear waves) of aluminum.
  double heat_capacity1; // Specific heat capacity at constant volume of aluminum.

  double alpha_param1; // Alpha parameter (nonlinear elasticity) of aluminum.
  double beta_param1; // Beta parameter (nonlinear elasticity) of aluminum.
  double gamma_param1; // Gamma parameter (nonlinear elasticity) of aluminum.

  double rho_ref2; // Reference density (unstressed configuration) of copper.
  double T_ref2; // Reference temperature (unstressed configuration) of copper.
  double sound_speed2; // Sound speed (speed of pressure waves) of copper.
  double shear_speed2; // Shear speed (speed of shear waves) of copper.
  double heat_capacity2; // Specific heat capacity at constant volume of copper.

  double alpha_param2; // Alpha parameter (nonlinear elasticity) of copper.
  double beta_param2; // Beta parameter (nonlinear elasticity) of copper.
  double gamma_param2; // Gamma parameter (nonlinear elasticity) of copper.

  double vel_x_l; // Left solid x-velocity.
  double vel_y_l; // Left solid y-velocity.
  double vel_z_l; // Left solid z-velocity.
  double Fxx_l; // Left solid deformation gradient tensor (xx-component).
  double Fxy_l; // Left solid deformation gradient tensor (xy-component).
  double Fxz_l; // Left solid deformation gradient tensor (xz-component).
  double Fyy_l; // Left solid deformation gradient tensor (yy-component).
  double Fzy_l; // Left solid deformation gradient tensor (zy-component).
  double Fzz_l; // Left solid deformation gradient tensor (zz-component).
  double entropy_l; // Left solid entropy.
  double phi1_l; // Left solid level set value (first species).

  double vel_x_r; // Right solid x-velocity.
  double vel_y_r; // Right solid y-velocity.
  double vel_z_r; // Right solid z-velocity.
  double Fxx_r; // Right solid deformation gradient tensor (xx-component).
  double Fxy_r; // Right solid deformation gradient tensor (xy-component).
  double Fxz_r; // Right solid deformation gradient tensor (xz-component).
  double Fyy_r; // Right solid deformation gradient tensor (yy-component).
  double Fzy_r; // Right solid deformation gradient tensor (zy-component).
  double Fzz_r; // Right solid deformation gradient tensor (zz-component).
  double entropy_r; // Right solid entropy.
  double phi1_r; // Right solid level set value (first species).

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.
  int reinit_freq; // Reinitialization frequency (for level set).
  double surface_tension; // Surface tension (for level set).

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct elasticity_rgfm_shock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double rho_ref1 = 2.71; // Reference density (unstressed configuration) of aluminum.
  double T_ref1 = 300.0; // Reference temperature (unstressed configuration) of aluminum.
  double sound_speed1 = 6.22; // Sound speed (speed of pressure waves) of aluminum.
  double shear_speed1 = 3.16; // Shear speed (speed of shear waves) of aluminum.
  double heat_capacity1 = 9.0 * pow(10.0, -4.0); // Specific heat capacity at constant volume of aluminum.

  double alpha_param1 = 1.0; // Alpha parameter (nonlinear elasticity) of aluminum.
  double beta_param1 = 3.577; // Beta parameter (nonlinear elasticity) of aluminum.
  double gamma_param1 = 2.088; // Gamma parameter (nonlinear elasticity) of aluminum.

  double rho_ref2 = 8.93; // Reference density (unstressed configuration) of copper.
  double T_ref2 = 300.0; // Reference temperature (unstressed configuration) of copper.
  double sound_speed2 = 4.6; // Sound speed (speed of pressure waves) of copper.
  double shear_speed2 = 2.1; // Shear speed (speed of shear waves) of copper.
  double heat_capacity2 = 3.9 * pow(10.0, -4.0); // Specific heat capacity at constant volume of copper.

  double alpha_param2 = 1.0; // Alpha parameter (nonlinear elasticity) of copper.
  double beta_param2 = 3.0; // Beta parameter (nonlinear elasticity) of copper.
  double gamma_param2 = 2.0; // Gamma parameter (nonlinear elasticity) of copper.

  double vel_x_l = 2.0; // Left solid x-velocity.
  double vel_y_l = 0.0; // Left solid y-velocity.
  double vel_z_l = 0.1; // Left solid z-velocity.
  double Fxx_l = 1.0; // Left solid deformation gradient tensor (xx-component).
  double Fxy_l = -0.01; // Left solid deformation gradient tensor (xy-component).
  double Fxz_l = -0.015; // Left solid deformation gradient tensor (xz-component).
  double Fyy_l = 0.95; // Left solid deformation gradient tensor (yy-component).
  double Fzy_l = 0.02; // Left solid deformation gradient tensor (zy-component).
  double Fzz_l = 0.9; // Left solid deformation gradient tensor (zz-component).
  double entropy_l = 0.0; // Left solid entropy.
  double phi1_l = 0.99999; // Left solid level set value (first species).

  double vel_x_r = 0.0; // Right solid x-velocity.
  double vel_y_r = -0.03; // Right solid y-velocity.
  double vel_z_r = -0.01; // Right solid z-velocity.
  double Fxx_r = 1.0; // Right solid deformation gradient tensor (xx-component).
  double Fxy_r = 0.015; // Right solid deformation gradient tensor (xy-component).
  double Fxz_r = -0.01; // Right solid deformation gradient tensor (xz-component).
  double Fyy_r = 0.95; // Right solid deformation gradient tensor (yy-component).
  double Fzy_r = 0.0; // Right solid deformation gradient tensor (zy-component).
  double Fzz_r = 0.9; // Right solid deformation gradient tensor (zz-component).
  double entropy_r = 0.0; // Right solid entropy.
  double phi1_r = 0.00001; // Right solid level set value (first species).

  // Simulation parameters.
  int Nx = 2048; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 0.8; // CFL coefficient.
  int reinit_freq = 3; // Reinitialization frequency (for level set).
  double surface_tension = 0.0; // Surface tension (for level set).

  double t_end = 0.05; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct elasticity_rgfm_shock_ctx ctx = {
    .rho_ref1 = rho_ref1,
    .T_ref1 = T_ref1,
    .sound_speed1 = sound_speed1,
    .shear_speed1 = shear_speed1,
    .heat_capacity1 = heat_capacity1,
    .alpha_param1 = alpha_param1,
    .beta_param1 = beta_param1,
    .gamma_param1 = gamma_param1,
    .rho_ref2 = rho_ref2,
    .T_ref2 = T_ref2,
    .sound_speed2 = sound_speed2,
    .shear_speed2 = shear_speed2,
    .heat_capacity2 = heat_capacity2,
    .alpha_param2 = alpha_param2,
    .beta_param2 = beta_param2,
    .gamma_param2 = gamma_param2,
    .vel_x_l = vel_x_l,
    .vel_y_l = vel_y_l,
    .vel_z_l = vel_z_l,
    .Fxx_l = Fxx_l,
    .Fxy_l = Fxy_l,
    .Fxz_l = Fxz_l,
    .Fyy_l = Fyy_l,
    .Fzy_l = Fzy_l,
    .Fzz_l = Fzz_l,
    .entropy_l = entropy_l,
    .phi1_l = phi1_l,
    .vel_x_r = vel_x_r,
    .vel_y_r = vel_y_r,
    .vel_z_r = vel_z_r,
    .Fxx_r = Fxx_r,
    .Fxy_r = Fxy_r,
    .Fxz_r = Fxz_r,
    .Fyy_r = Fyy_r,
    .Fzy_r = Fzy_r,
    .Fzz_r = Fzz_r,
    .entropy_r = entropy_r,
    .phi1_r = phi1_r,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .reinit_freq = reinit_freq,
    .surface_tension = surface_tension,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElasticityRGFMInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct elasticity_rgfm_shock_ctx *app = ctx;

  double rho_ref1 = app->rho_ref1;
  double T_ref1 = app->T_ref1;
  double sound_speed1 = app->sound_speed1;
  double shear_speed1 = app->shear_speed1;
  double heat_capacity1 = app->heat_capacity1;

  double alpha_param1 = app->alpha_param1;
  double beta_param1 = app->beta_param1;
  double gamma_param1 = app->gamma_param1;

  double rho_ref2 = app->rho_ref2;
  double T_ref2 = app->T_ref2;
  double sound_speed2 = app->sound_speed2;
  double shear_speed2 = app->shear_speed2;
  double heat_capacity2 = app->heat_capacity2;

  double alpha_param2 = app->alpha_param2;
  double beta_param2 = app->beta_param2;
  double gamma_param2 = app->gamma_param2;

  double vel_x_l = app->vel_x_l;
  double vel_y_l = app->vel_y_l;
  double vel_z_l = app->vel_z_l;
  double Fxx_l = app->Fxx_l;
  double Fxy_l = app->Fxy_l;
  double Fxz_l = app->Fxz_l;
  double Fyy_l = app->Fyy_l;
  double Fzy_l = app->Fzy_l;
  double Fzz_l = app->Fzz_l;
  double entropy_l = app->entropy_l;
  double phi1_l = app->phi1_l;

  double vel_x_r = app->vel_x_r;
  double vel_y_r = app->vel_y_r;
  double vel_z_r = app->vel_z_r;
  double Fxx_r = app->Fxx_r;
  double Fxy_r = app->Fxy_r;
  double Fxz_r = app->Fxz_r;
  double Fyy_r = app->Fyy_r;
  double Fzy_r = app->Fzy_r;
  double Fzz_r = app->Fzz_r;
  double entropy_r = app->entropy_r;
  double phi1_r = app->phi1_r;

  double phi1 = 0.0;
  double vel_x_total = 0.0;
  double vel_y_total = 0.0;
  double vel_z_total = 0.0;
  double Fxx_total = 0.0;
  double Fxy_total = 0.0;
  double Fxz_total = 0.0;
  double Fyy_total = 0.0;
  double Fzy_total = 0.0;
  double Fzz_total = 0.0;
  double entropy_total = 0.0;

  if (x < 0.5) {
    phi1 = phi1_l; // First species level set value (left).

    vel_x_total = vel_x_l; // Total solid x-velocity (left).
    vel_y_total = vel_y_l; // Total solid x-velocity (left).
    vel_z_total = vel_z_l; // Total solid x-velocity (left).

    Fxx_total = Fxx_l; // Total solid deformation gradient tensor (xx-component, left).
    Fxy_total = Fxy_l; // Total solid deformation gradient tensor (xy-component, left).
    Fxz_total = Fxz_l; // Total solid deformation gradient tensor (xz-component, left).
    Fyy_total = Fyy_l; // Total solid deformation gradient tensor (yy-component, left).
    Fzy_total = Fzy_l; // Total solid deformation gradient tensor (zy-component, left).
    Fzz_total = Fzz_l; // Total solid deformation gradient tensor (zz-component, left).

    entropy_total = entropy_l; // Total solid entropy (left).
  }
  else {
    phi1 = phi1_r; // First species level set value (right).

    vel_x_total = vel_x_r; // Total solid x-velocity (right).
    vel_y_total = vel_y_r; // Total solid x-velocity (right).
    vel_z_total = vel_z_r; // Total solid x-velocity (right).

    Fxx_total = Fxx_r; // Total solid deformation gradient tensor (xx-component, right).
    Fxy_total = Fxy_r; // Total solid deformation gradient tensor (xy-component, right).
    Fxz_total = Fxz_r; // Total solid deformation gradient tensor (xz-component, right).
    Fyy_total = Fyy_r; // Total solid deformation gradient tensor (yy-component, right).
    Fzy_total = Fzy_r; // Total solid deformation gradient tensor (zy-component, right).
    Fzz_total = Fzz_r; // Total solid deformation gradient tensor (zz-component, right).

    entropy_total = entropy_r; // Total solid entropy (right).
  }

  double deformation_gradient_total[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        deformation_gradient_total[i][j] = 1.0;
      }
      else {
        deformation_gradient_total[i][j] = 0.0;
      }
    }
  }

  deformation_gradient_total[0][0] = Fxx_total;
  deformation_gradient_total[0][1] = Fxy_total;
  deformation_gradient_total[0][2] = Fxz_total;
  deformation_gradient_total[1][1] = Fyy_total;
  deformation_gradient_total[2][1] = Fzy_total;
  deformation_gradient_total[2][2] = Fzz_total;

  double deformation_gradient_total_det = (deformation_gradient_total[0][0] * ((deformation_gradient_total[1][1] * deformation_gradient_total[2][2]) -
    (deformation_gradient_total[2][1] * deformation_gradient_total[1][2]))) -
    (deformation_gradient_total[0][1] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][2]) - (deformation_gradient_total[1][2] * deformation_gradient_total[2][0]))) +
    (deformation_gradient_total[0][2] * ((deformation_gradient_total[1][0] * deformation_gradient_total[2][1]) - (deformation_gradient_total[1][1] * deformation_gradient_total[2][0])));

  double rho1 = rho_ref1 / deformation_gradient_total_det;
  double rho2 = rho_ref2 / deformation_gradient_total_det;
  double rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2); // Total solid mass density.

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
    (pow(strain_invariant3_total, 0.5 * alpha_param1) - 1.0)) + (heat_capacity1 * T_ref1 * pow(strain_invariant3_total, 0.5 * gamma_param1) * (exp(entropy_total / heat_capacity1) - 1.0)) +
    (0.5 * shear_modulus1 * pow(strain_invariant3_total, 0.5 * beta_param1) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));
  double internal_energy2 = ((bulk_modulus2 / (2.0 * (alpha_param2 * alpha_param2))) * (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0) *
    (pow(strain_invariant3_total, 0.5 * alpha_param2) - 1.0)) + (heat_capacity2 * T_ref2 * pow(strain_invariant3_total, 0.5 * gamma_param2) * (exp(entropy_total / heat_capacity2) - 1.0)) +
    (0.5 * shear_modulus2 * pow(strain_invariant3_total, 0.5 * beta_param2) * (((strain_invariant1_total * strain_invariant1_total) / 3.0) - strain_invariant2_total));

  double mom_x_total = rho_total * vel_x_total; // Total solid momentum density (x-direction).
  double mom_y_total = rho_total * vel_y_total; // Total solid momentum density (y-direction).
  double mom_z_total = rho_total * vel_z_total; // Total solid momentum density (z-direction).

  double deformation_gradient_total_cons[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      deformation_gradient_total_cons[i][j] = rho_total * deformation_gradient_total[i][j]; // Solid deformation gradient tensor.
    }
  }

  double E1 = rho1 * (internal_energy1 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E2 = rho2 * (internal_energy2 + (0.5 * ((vel_x_total * vel_x_total) + (vel_y_total * vel_y_total) + (vel_z_total * vel_z_total))));
  double E_total = (phi1 * E1) + ((1.0 - phi1) * E2); // Total solid energy.

  double level_set1 = rho_total * phi1; // Conserved level set value (first species).
  double mass_frac1 = phi1 * rho1; // Conserved mass density (first species).
  double mass_frac2 = (1.0 - phi1) * rho2; // Conserved mass density (second species).

  // Set total solid mass density.
  fout[0] = rho_total;
  // Set total solid momentum density.
  fout[1] = mom_x_total; fout[2] = mom_y_total; fout[3] = mom_z_total;

  // Set total solid deformation gradient tensor
  fout[4] = deformation_gradient_total_cons[0][0]; fout[5] = deformation_gradient_total_cons[0][1]; fout[6] = deformation_gradient_total_cons[0][2];
  fout[7] = deformation_gradient_total_cons[1][0]; fout[8] = deformation_gradient_total_cons[1][1]; fout[9] = deformation_gradient_total_cons[1][2];
  fout[10] = deformation_gradient_total_cons[2][0]; fout[11] = deformation_gradient_total_cons[2][1]; fout[12] = deformation_gradient_total_cons[2][2];

  // Set total solid energy density.
  fout[13] = E_total;

  // Set conserved level set value (first species).
  fout[14] = level_set1;
  // Set conserved mass densities (first and second species).
  fout[15] = mass_frac1; fout[16] = mass_frac2;
  // Set reinitialization parameter (for level set).
  fout[17] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct elasticity_rgfm_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Fluid equations.
  double *rho_ref_s = gkyl_malloc(sizeof(double[2]));
  rho_ref_s[0] = ctx.rho_ref1;
  rho_ref_s[1] = ctx.rho_ref2;

  double *T_ref_s = gkyl_malloc(sizeof(double[2]));
  T_ref_s[0] = ctx.T_ref1;
  T_ref_s[1] = ctx.T_ref2;

  double *sound_speed_s = gkyl_malloc(sizeof(double[2]));
  sound_speed_s[0] = ctx.sound_speed1;
  sound_speed_s[1] = ctx.sound_speed2;
  
  double *shear_speed_s = gkyl_malloc(sizeof(double[2]));
  shear_speed_s[0] = ctx.shear_speed1;
  shear_speed_s[1] = ctx.shear_speed2;
  
  double *heat_capacity_s = gkyl_malloc(sizeof(double[2]));
  heat_capacity_s[0] = ctx.heat_capacity1;
  heat_capacity_s[1] = ctx.heat_capacity2;

  double *alpha_param_s = gkyl_malloc(sizeof(double[2]));
  alpha_param_s[0] = ctx.alpha_param1;
  alpha_param_s[1] = ctx.alpha_param2;

  double *beta_param_s = gkyl_malloc(sizeof(double[2]));
  beta_param_s[0] = ctx.beta_param1;
  beta_param_s[1] = ctx.beta_param2;

  double *gamma_param_s = gkyl_malloc(sizeof(double[2]));
  gamma_param_s[0] = ctx.gamma_param1;
  gamma_param_s[1] = ctx.gamma_param2;

  struct gkyl_wv_eqn *elasticity_rgfm = gkyl_wv_elasticity_rgfm_new(2, rho_ref_s, T_ref_s, sound_speed_s, shear_speed_s, heat_capacity_s, alpha_param_s, beta_param_s, gamma_param_s,
    ctx.reinit_freq, ctx.surface_tension, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "elasticity_rgfm",
    .equation = elasticity_rgfm,
    
    .init = evalElasticityRGFMInit,
    .ctx = &ctx,
    .force_low_order_flux = false,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "elasticity_rgfm_shock",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elasticity_rgfm);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
