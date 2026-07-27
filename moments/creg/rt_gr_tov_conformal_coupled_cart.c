// Coupled TOV star in Cartesian coordinates: GR Euler (gamma-law) fluid dynamically coupled to the conformal vacuum Einstein equations,
//  both evolved with the MP flux-form scheme.  Coupled through the MP scheme by:
//  - spacetime -> fluid: the live metric is synced into the fluid metric slots each RK stage (eqns 52/53 from the Multifluids paper)
//  - fluid -> spacetime: the fluid stress-energy adds the standard ADM matter sources  -8*pi*alpha*(T_ij - 1/2 gamma_ij T) to K_ij and +8*pi*alpha*T^0_k to V_k.
// The fluid still has the frozen-discrete well-balancing (toggle with gr_euler_disable_well_balanced). No WB for the spacetime for now.
// Both species have the same analytic TOV spacetime object for their initial data only.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_moment.h>
#include <gkyl_moment_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h> // gkyl_gr_euler_prim_vars, for the admissibility/floor diagnostic.
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_gr_tov.h>
#include "tov_solver.h"

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct gr_tov_coupled_ctx
{
  double gas_gamma; // Adiabatic index (evolution EOS).
  double K_poly; // Polytropic constant of the cold-IC EOS (p = K rho^gamma).
  double rho_c; // Central rest-mass density.
  double dr_tov; // TOV solver radial step.
  double rho_atm; // Atmosphere rest-mass density.
  double p_atm; // Atmosphere pressure floor.

  struct gkyl_tov *tov; // Frozen TOV table.
  struct gkyl_gr_spacetime *spacetime; // TOV spacetime object (shared for initial data only).
  double M_star;
  double R_star;

  enum gkyl_spacetime_gauge spacetime_gauge; // Fluid-side gauge tag (GR Euler eqn).
  int reinit_freq;

  // Conformal vacuum Einstein evolution parameters.
  double excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution;

  double perturb_frac; // x-momentum kick fraction in the inner star (0 => static test).

  int Nx, Ny;
  double Lx, Ly; // Domain [-Lx/2, Lx/2] x [-Ly/2, Ly/2], star centered at origin.
  double cfl_frac;

  double t_end;
  int num_frames;
  int field_energy_calcs;
  int integrated_mom_calcs;
  double dt_failure_tol;
  int num_failures_max;
};

struct gr_tov_coupled_ctx
create_ctx(void)
{
  double gas_gamma = 2.0;
  double K_poly = 100.0;
  double rho_c = 1.28e-3;
  double dr_tov = 0.01;

  double rho_atm = 1e-7 * rho_c;
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);
  double M_star = gkyl_tov_star_mass(tov);
  double R_star = gkyl_tov_star_radius(tov);
  printf("M_star = %e \n", M_star);
  printf("R_star = %e \n", R_star);
  printf("Compactness (2M_star / R_star) = %e \n", 2.0 * M_star / R_star);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, 0.0, 0.0, 0.0, false);

  struct gr_tov_coupled_ctx ctx = {
    .gas_gamma = gas_gamma,
    .K_poly = K_poly,
    .rho_c = rho_c,
    .dr_tov = dr_tov,
    .rho_atm = rho_atm,
    .p_atm = p_atm,
    .tov = tov,
    .spacetime = spacetime,
    .M_star = M_star,
    .R_star = R_star,
    .spacetime_gauge = GKYL_STATIC_GAUGE,
    .reinit_freq = 100,
    .excision_threshold = 0.3,
    .spacetime_slicing = GKYL_HARMONIC_SLICING,
    .spacetime_evolution = GKYL_EINSTEIN_EVOLUTION,
    .perturb_frac = 0.0, // set > 0 (e.g. 1.0e-3) for the momentum-perturbation test.
    .Nx = 200,
    .Ny = 200,
    .Lx = 24.0,
    .Ly = 24.0,
    .cfl_frac = 0.8,
    .t_end = 20.0,
    .num_frames = 40,
    .field_energy_calcs = INT_MAX,
    .integrated_mom_calcs = INT_MAX,
    .dt_failure_tol = 1.0e-4,
    .num_failures_max = 20,
  };

  return ctx;
}

void
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gr_tov_coupled_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double r = sqrt((x * x) + (y * y));
  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(app->tov, r, &bl);

  double p = fmax(bl.P, app->p_atm);
  double rho = pow(p / app->K_poly, 1.0 / gas_gamma);

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  double dl = pow(10.0, -8.0);
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &spatial_metric_der);

  // Once and for all with the velocity (I'll make a note for myself for later): In the Kerr-Schild slicing the static star is not at rest w.r.t. the Eulerian observer: 
  // a fluid element following the static Killing vector (zero coordinate velocity) has Eulerian 3-velocity v^i = beta^i / alpha (physical speed = 2m/r, subluminal for a regular star). 
  // v=0 would initialize a moving fluid. (For the areal gauge beta^i=0, so this reduces to v=0.)
  double vel[3] = { shift[0] / lapse, shift[1] / lapse, shift[2] / lapse };
  double v_sq = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];

  double W = 1.0 / sqrt(1.0 - v_sq);
  if (v_sq > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  double sqrt_det = sqrt(spatial_det);
  double cov_vel[3] = { 0.0 };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cov_vel[i] += spatial_metric[i][j] * vel[j];
    }
  }

  double rho_rel = sqrt_det * rho * W;
  double mom_x = sqrt_det * rho * h * (W * W) * cov_vel[0];
  double mom_y = sqrt_det * rho * h * (W * W) * cov_vel[1];
  double mom_z = sqrt_det * rho * h * (W * W) * cov_vel[2];
  double Etot = sqrt_det * ((rho * h * (W * W)) - p - (rho * W));

  fout[0] = rho_rel;
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  fout[4] = Etot;

  fout[5] = lapse;
  fout[6] = shift[0]; fout[7] = shift[1]; fout[8] = shift[2];

  fout[9]  = spatial_metric[0][0]; fout[10] = spatial_metric[0][1]; fout[11] = spatial_metric[0][2];
  fout[12] = spatial_metric[1][0]; fout[13] = spatial_metric[1][1]; fout[14] = spatial_metric[1][2];
  fout[15] = spatial_metric[2][0]; fout[16] = spatial_metric[2][1]; fout[17] = spatial_metric[2][2];

  fout[18] = extrinsic_curvature[0][0]; fout[19] = extrinsic_curvature[0][1]; fout[20] = extrinsic_curvature[0][2];
  fout[21] = extrinsic_curvature[1][0]; fout[22] = extrinsic_curvature[1][1]; fout[23] = extrinsic_curvature[1][2];
  fout[24] = extrinsic_curvature[2][0]; fout[25] = extrinsic_curvature[2][1]; fout[26] = extrinsic_curvature[2][2];

  fout[27] = in_excision_region ? -1.0 : 1.0;

  fout[28] = lapse_der[0]; fout[29] = lapse_der[1]; fout[30] = lapse_der[2];

  fout[31] = shift_der[0][0]; fout[32] = shift_der[0][1]; fout[33] = shift_der[0][2];
  fout[34] = shift_der[1][0]; fout[35] = shift_der[1][1]; fout[36] = shift_der[1][2];
  fout[37] = shift_der[2][0]; fout[38] = shift_der[2][1]; fout[39] = shift_der[2][2];

  fout[40] = spatial_metric_der[0][0][0]; fout[41] = spatial_metric_der[0][0][1]; fout[42] = spatial_metric_der[0][0][2];
  fout[43] = spatial_metric_der[0][1][0]; fout[44] = spatial_metric_der[0][1][1]; fout[45] = spatial_metric_der[0][1][2];
  fout[46] = spatial_metric_der[0][2][0]; fout[47] = spatial_metric_der[0][2][1]; fout[48] = spatial_metric_der[0][2][2];

  fout[49] = spatial_metric_der[1][0][0]; fout[50] = spatial_metric_der[1][0][1]; fout[51] = spatial_metric_der[1][0][2];
  fout[52] = spatial_metric_der[1][1][0]; fout[53] = spatial_metric_der[1][1][1]; fout[54] = spatial_metric_der[1][1][2];
  fout[55] = spatial_metric_der[1][2][0]; fout[56] = spatial_metric_der[1][2][1]; fout[57] = spatial_metric_der[1][2][2];

  fout[58] = spatial_metric_der[2][0][0]; fout[59] = spatial_metric_der[2][0][1]; fout[60] = spatial_metric_der[2][0][2];
  fout[61] = spatial_metric_der[2][1][0]; fout[62] = spatial_metric_der[2][1][1]; fout[63] = spatial_metric_der[2][1][2];
  fout[64] = spatial_metric_der[2][2][0]; fout[65] = spatial_metric_der[2][2][1]; fout[66] = spatial_metric_der[2][2][2];

  fout[67] = 0.0;
  fout[68] = x; fout[69] = y; fout[70] = 0.0;

  // Frozen-discrete well-balancing reference (num_equations = 73): t=0 conserved fluid in slots 71,72.
  fout[71] = rho_rel;
  fout[72] = Etot;

  // Optional momentum perturbation in the inner star
  if (app->perturb_frac > 0.0 && r < 0.5 * app->R_star) {
    fout[1] += app->perturb_frac * rho_rel;
  }

  if (in_excision_region) {
    for (int i = 0; i < 68; i++) fout[i] = 0.0;
    fout[27] = -1.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) { gkyl_free(spatial_metric[i]); gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]); }
  gkyl_free(spatial_metric); gkyl_free(extrinsic_curvature); gkyl_free(lapse_der); gkyl_free(shift_der);
  for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]); gkyl_free(spatial_metric_der[i]); }
  gkyl_free(spatial_metric_der);
}

void
evalVacuumEinsteinConformalInit(double t, const double* GKYL_RESTRICT xn,
  double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gr_tov_coupled_ctx *app = ctx;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double lapse, psi, chi;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **conformal_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) conformal_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **inv_conformal_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) inv_conformal_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***conformal_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    conformal_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) conformal_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  double *psi_der = gkyl_malloc(sizeof(double[3]));
  double *chi_der = gkyl_malloc(sizeof(double[3]));
  double **chi_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) chi_der2[i] = gkyl_malloc(sizeof(double[3]));

  double dl = 1.0e-8;
  double dl2 = 1.0e-6;
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &conformal_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_conformal_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &extrinsic_curvature);
  spacetime->conformal_factor_func(spacetime, 0.0, x, y, 0.0, &psi);
  spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, y, 0.0, &chi);
  spacetime->conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &psi_der);
  spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &chi_der);
  spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, y, 0.0,
    dl2, dl2, dl2, &chi_der2);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0,
    dl, dl, dl, &conformal_metric_der);

  double psi4 = psi * psi * psi * psi;
  for (int k = 0; k < 3; k++) psi_der[k] /= psi;

  // Convert gamma_ij and its inverse to tilde(gamma)_ij and
  // tilde(gamma)^ij. The derivative reduction variable is
  // tilde(D)_kij = 1/2 d_k tilde(gamma)_ij.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_metric[i][j] /= psi4;
      inv_conformal_metric[i][j] *= psi4;
      for (int k = 0; k < 3; k++) {
        conformal_metric_der[k][i][j] =
          (0.5 * conformal_metric_der[k][i][j] / psi4)
          - (2.0 * psi_der[k] * conformal_metric[i][j]);
      }
      shift_der[i][j] = 0.5 * shift_der[i][j];
    }
  }
  for (int i = 0; i < 3; i++) lapse_der[i] /= lapse;

  double conformal_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        conformal_metric_der_raised1[k][i][j] = 0.0;
        for (int l = 0; l < 3; l++)
          conformal_metric_der_raised1[k][i][j] +=
            inv_conformal_metric[k][l] * conformal_metric_der[l][i][j];
      }

  double conformal_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        conformal_metric_der_raised3[i][j][k] = 0.0;
        for (int l = 0; l < 3; l++)
          conformal_metric_der_raised3[i][j][k] +=
            inv_conformal_metric[l][k] * conformal_metric_der[i][j][l];
      }

  double conformal_aux_vect[3];
  for (int i = 0; i < 3; i++) {
    conformal_aux_vect[i] = -4.0 * psi_der[i];
    for (int s = 0; s < 3; s++) {
      conformal_aux_vect[i] += conformal_metric_der_raised3[i][s][s];
      conformal_aux_vect[i] -= conformal_metric_der_raised1[s][s][i];
    }
  }

  fout[0] = conformal_metric[0][0]; fout[1] = conformal_metric[0][1]; fout[2] = conformal_metric[0][2];
  fout[3] = conformal_metric[1][0]; fout[4] = conformal_metric[1][1]; fout[5] = conformal_metric[1][2];
  fout[6] = conformal_metric[2][0]; fout[7] = conformal_metric[2][1]; fout[8] = conformal_metric[2][2];

  fout[9] = lapse;

  fout[10] = extrinsic_curvature[0][0]; fout[11] = extrinsic_curvature[0][1]; fout[12] = extrinsic_curvature[0][2];
  fout[13] = extrinsic_curvature[1][0]; fout[14] = extrinsic_curvature[1][1]; fout[15] = extrinsic_curvature[1][2];
  fout[16] = extrinsic_curvature[2][0]; fout[17] = extrinsic_curvature[2][1]; fout[18] = extrinsic_curvature[2][2];

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fout[19 + (9 * k) + (3 * i) + j] = conformal_metric_der[k][i][j];
      }
    }
  }

  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];
  fout[49] = conformal_aux_vect[0]; fout[50] = conformal_aux_vect[1]; fout[51] = conformal_aux_vect[2];

  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];

  fout[55] = shift_der[0][0]; fout[56] = shift_der[0][1]; fout[57] = shift_der[0][2];
  fout[58] = shift_der[1][0]; fout[59] = shift_der[1][1]; fout[60] = shift_der[1][2];
  fout[61] = shift_der[2][0]; fout[62] = shift_der[2][1]; fout[63] = shift_der[2][2];

  fout[64] = chi;
  fout[65] = chi_der[0]; fout[66] = chi_der[1]; fout[67] = chi_der[2];
  fout[68] = chi_der2[0][0]; fout[69] = chi_der2[0][1]; fout[70] = chi_der2[0][2];
  fout[71] = chi_der2[1][0]; fout[72] = chi_der2[1][1]; fout[73] = chi_der2[1][2];
  fout[74] = chi_der2[2][0]; fout[75] = chi_der2[2][1]; fout[76] = chi_der2[2][2];

  if (in_excision_region) {
    for (int i = 0; i < 77; i++) fout[i] = 0.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) {
    gkyl_free(conformal_metric[i]); gkyl_free(inv_conformal_metric[i]);
    gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]);
    gkyl_free(chi_der2[i]);
    for (int j = 0; j < 3; j++) gkyl_free(conformal_metric_der[i][j]);
    gkyl_free(conformal_metric_der[i]);
  }
  gkyl_free(conformal_metric); gkyl_free(inv_conformal_metric); gkyl_free(extrinsic_curvature);
  gkyl_free(lapse_der); gkyl_free(shift_der); gkyl_free(conformal_metric_der);
  gkyl_free(psi_der); gkyl_free(chi_der); gkyl_free(chi_der2);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = (force_write) ? iot->curr : iot->curr - 1;
    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
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
#ifdef GKYL_HAVE_MPI
  if (argc > 1) {
    struct gkyl_app_args app_args_pre = parse_app_args(argc, argv);
    if (app_args_pre.use_mpi) MPI_Init(&argc, &argv);
  }
#endif

  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gr_tov_coupled_ctx ctx = create_ctx();

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // GR Euler fluid (well-balanced against the frozen TOV equilibrium).
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
    .gas_gamma = ctx.gas_gamma,
    .spacetime_gauge = ctx.spacetime_gauge,
    .reinit_freq = ctx.reinit_freq,
    .spacetime = ctx.spacetime,
    .rp_type = WV_GR_EULER_RP_LAX,
    .use_gpu = app_args.use_gpu,
    .tov_eq = ctx.tov,
    .p_atm = ctx.p_atm,
    .rho_atm = ctx.rho_atm,
  });

  // Conformal vacuum Einstein (Bona-Masso), evolved with MP and coupled to the fluid.
  struct gkyl_wv_eqn *vacuum_einstein_conformal = gkyl_wv_vacuum_einstein_conformal_new(ctx.excision_threshold,
    ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    .init = evalGREulerInit,
    .ctx = &ctx,

    .has_gr_euler = true,
    .gr_euler_gas_gamma = ctx.gas_gamma,
    //  Note to self:
    // Frozen-discrete WB (zero-momentum equilibrium) holds the fluid bit-statically in the areal gauge.
    // For the Kerr-Schild gauge the static star has v=beta/alpha != 0, so the frozen WB is invalid there
    // and the equilibrium-family WB (wb_family) is needed - which in turn needs the near-vacuum limiter.

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species einstein_conformal = {
    .name = "vacuum_einstein_conformal",
    .equation = vacuum_einstein_conformal,
    .init = evalVacuumEinsteinConformalInit,
    .ctx = &ctx,

    .has_vacuum_einstein_conformal = true,
    .vacuum_einstein_conformal_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_conformal_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_conformal_spacetime_evolution = ctx.spacetime_evolution,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int cells[] = { NX, NY };
  int dim = sizeof(cells) / sizeof(cells[0]);
  int cuts[dim];
  for (int d = 0; d < dim; d++) cuts[d] = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) for (int d = 0; d < dim; d++) cuts[d] = app_args.cuts[d];
#endif

  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) { .mpi_comm = MPI_COMM_WORLD } );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app_args.use_gpu } );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app_args.use_gpu } );
#endif

  struct gkyl_moment app_inp = {
    .name = "gr_tov_conformal_coupled_cart",

    .ndim = 2,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly },
    .upper = {  0.5 * ctx.Lx,  0.5 * ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_MP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { fluid, einstein_conformal },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { cuts[0], cuts[1] },
      .comm = comm,
    },
  };

  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end;
  int frame_curr = 0;
  gkyl_moment_app_apply_ic(app, t_curr);

  struct gkyl_tm_trigger io_trig = { .dt = t_end / ctx.num_frames, .tcurr = t_curr, .curr = frame_curr };
  write_data(&io_trig, app, t_curr, false);

  double dt = fmin(t_end - t_curr, gkyl_moment_app_max_dt(app));
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

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) dt_init = status.dt_actual;
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;
      gkyl_moment_app_cout(app, stdout, "WARNING: dt = %g below %g*dt_init; num_failures = %d\n",
        status.dt_actual, dt_failure_tol, num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: too many small steps. Aborting.\n");
        break;
      }
    }
    else num_failures = 0;

    step += 1;
  }

  write_data(&io_trig, app, t_curr, true);
  gkyl_moment_app_stat_write(app);

  gkyl_moment_app_release(app);
  gkyl_wv_eqn_release(gr_euler);
  gkyl_wv_eqn_release(vacuum_einstein_conformal);
  gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
  gkyl_tov_solution_release(ctx.tov);
  gkyl_comm_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Finalize();
#endif

  return 0;
}
