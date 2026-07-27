// Mixed-scheme coupled TOV star in Cartesian coordinates: GR Euler (gamma-law) fluid evolved with the MP flux-form scheme, vacuum Einstein (Bona-Masso) spacetime evolved with the
//  wave-propagation scheme, coupled each step.
//  Coupling once per step:
//  - spacetime -> fluid: the live metric is synced into the fluid metric slots at the start of the step (eqns 52/53 from the Multifluids paper);
//  - fluid -> spacetime: the fluid stress-energy adds the standard ADM matter sources
//    -8*pi*alpha*(T_ij - 1/2 gamma_ij T) to K_ij and +8*pi*alpha*T^0_k to V_k, folded into the spacetime's Strang-split half-step source integration.
// The fluid still has the frozen-discrete well-balancing (toggle with gr_euler_disable_well_balanced).
// Both species share the same analytic TOV spacetime object for their initial data only.

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_moment.h>
#include <gkyl_moment_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wv_vacuum_einstein.h>
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

  // Vacuum Einstein evolution parameters.
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

  // Gauge selection (env): GKYL_USE_KERR_SCHILD=1 -> Kerr-Schild geometry (horizon-penetrating, for collapse); default is the areal (static Schwarzschild-like) gauge. 
  // GKYL_SLICING = harmonic | 1pluslog | geodesic (default 1pluslog).
  bool use_kerr_schild = (getenv("GKYL_USE_KERR_SCHILD") != NULL);
  enum gkyl_spacetime_slicing slicing = GKYL_1PLUSLOG_SLICING;
  const char *slc = getenv("GKYL_SLICING");
  if (slc) {
    if (strcmp(slc, "harmonic") == 0) slicing = GKYL_HARMONIC_SLICING;
    else if (strcmp(slc, "geodesic") == 0) slicing = GKYL_GEODESIC_SLICING;
    else if (strcmp(slc, "1pluslog") == 0) slicing = GKYL_1PLUSLOG_SLICING;
  }
  printf("Gauge: use_kerr_schild=%d, slicing=%d\n", (int) use_kerr_schild, (int) slicing);
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, 0.0, 0.0, 0.0, use_kerr_schild);

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
    .spacetime_slicing = slicing, // env-selectable (default 1+log); see gauge sweep above
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

  // Frozen-discrete well-balancing reference: t=0 conserved equilibrium fluid.
  fout[71] = rho_rel;
  fout[72] = Etot;
  fout[73] = mom_x;
  fout[74] = mom_y;
  fout[75] = mom_z;

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
evalVacuumEinsteinInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct gr_tov_coupled_ctx *app = ctx;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

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
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, dl, dl, dl, &spatial_metric_der);

  // Bona-Masso conventions: 
  // D_kij = 1/2 d_k gamma_ij, (1/2) d_i beta^j, A_i = d_i ln(alpha)
  // (So what I calculate fits the vacuum einstein later)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k];
      }
      shift_der[i][j] = 0.5 * shift_der[i][j];
    }
  }
  for (int i = 0; i < 3; i++) {
    lapse_der[i] = lapse_der[i] / lapse;
  }

  double spatial_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised1[k][i][j] = 0.0;
        for (int l = 0; l < 3; l++)
          spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
      }

  double spatial_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised3[i][j][k] = 0.0;
        for (int l = 0; l < 3; l++)
          spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
      }

  double aux_vect[3];
  for (int i = 0; i < 3; i++) {
    aux_vect[i] = 0.0;
    for (int s = 0; s < 3; s++) {
      aux_vect[i] += spatial_metric_der_raised3[i][s][s];
      aux_vect[i] -= spatial_metric_der_raised1[s][s][i];
    }
  }

  fout[0] = spatial_metric[0][0]; fout[1] = spatial_metric[0][1]; fout[2] = spatial_metric[0][2];
  fout[3] = spatial_metric[1][0]; fout[4] = spatial_metric[1][1]; fout[5] = spatial_metric[1][2];
  fout[6] = spatial_metric[2][0]; fout[7] = spatial_metric[2][1]; fout[8] = spatial_metric[2][2];

  fout[9] = lapse;

  fout[10] = extrinsic_curvature[0][0]; fout[11] = extrinsic_curvature[0][1]; fout[12] = extrinsic_curvature[0][2];
  fout[13] = extrinsic_curvature[1][0]; fout[14] = extrinsic_curvature[1][1]; fout[15] = extrinsic_curvature[1][2];
  fout[16] = extrinsic_curvature[2][0]; fout[17] = extrinsic_curvature[2][1]; fout[18] = extrinsic_curvature[2][2];

  fout[19] = spatial_metric_der[0][0][0]; fout[20] = spatial_metric_der[0][0][1]; fout[21] = spatial_metric_der[0][0][2];
  fout[22] = spatial_metric_der[0][1][0]; fout[23] = spatial_metric_der[0][1][1]; fout[24] = spatial_metric_der[0][1][2];
  fout[25] = spatial_metric_der[0][2][0]; fout[26] = spatial_metric_der[0][2][1]; fout[27] = spatial_metric_der[0][2][2];

  fout[28] = spatial_metric_der[1][0][0]; fout[29] = spatial_metric_der[1][0][1]; fout[30] = spatial_metric_der[1][0][2];
  fout[31] = spatial_metric_der[1][1][0]; fout[32] = spatial_metric_der[1][1][1]; fout[33] = spatial_metric_der[1][1][2];
  fout[34] = spatial_metric_der[1][2][0]; fout[35] = spatial_metric_der[1][2][1]; fout[36] = spatial_metric_der[1][2][2];

  fout[37] = spatial_metric_der[2][0][0]; fout[38] = spatial_metric_der[2][0][1]; fout[39] = spatial_metric_der[2][0][2];
  fout[40] = spatial_metric_der[2][1][0]; fout[41] = spatial_metric_der[2][1][1]; fout[42] = spatial_metric_der[2][1][2];
  fout[43] = spatial_metric_der[2][2][0]; fout[44] = spatial_metric_der[2][2][1]; fout[45] = spatial_metric_der[2][2][2];

  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];

  fout[49] = aux_vect[0]; fout[50] = aux_vect[1]; fout[51] = aux_vect[2];

  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];

  fout[55] = shift_der[0][0]; fout[56] = shift_der[0][1]; fout[57] = shift_der[0][2];
  fout[58] = shift_der[1][0]; fout[59] = shift_der[1][1]; fout[60] = shift_der[1][2];
  fout[61] = shift_der[2][0]; fout[62] = shift_der[2][1]; fout[63] = shift_der[2][2];

  if (in_excision_region) {
    for (int i = 0; i < 64; i++) fout[i] = 0.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]); gkyl_free(inv_spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]);
    for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric); gkyl_free(inv_spatial_metric); gkyl_free(extrinsic_curvature);
  gkyl_free(lapse_der); gkyl_free(shift_der); gkyl_free(spatial_metric_der);
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

// [DIAG] Per-step mixed TOV tracker. Fluid quantities come from species 0; lapse comes from
// the live Einstein species 1, matching the collapse diagnostic.
static double diag_SumD0 = -1.0;
static void
diag_track(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  const struct gkyl_array *fluid = gkyl_moment_app_get_write_array_species(app, 0);
  const struct gkyl_array *einstein = gkyl_moment_app_get_write_array_species(app, 1);

  double SumD_local = 0.0;
  double Dmax_local = 0.0, r_Dmax_local = 0.0;
  double Dc_local = 0.0, lapse_c_local = 0.0, r_c_local = 1.0e300;
  double core_mom_local = 0.0, surf_mom_local = 0.0;
  double glob_mom_local = 0.0, r_glob_mom_local = 0.0;
  double lapse_min_local = 1.0e300, lapse_max_local = -1.0e300;
  int n_bad_local = 0;

  double R_core = 0.5 * ctx->R_star;
  double R_star = ctx->R_star;
  double surf_lo = R_star - 0.15 * R_star, surf_hi = R_star + 0.15 * R_star;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(fluid, loc);
    const double *es = gkyl_array_cfetch(einstein, loc);

    double xc[3];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1]);

    double D = fl[0];
    double lapse = es[9];
    double mom = sqrt(fl[1] * fl[1] + fl[2] * fl[2] + fl[3] * fl[3]);

    if (!isfinite(D) || !isfinite(lapse) || !isfinite(mom)) n_bad_local += 1;

    SumD_local += D;
    if (D > Dmax_local) { Dmax_local = D; r_Dmax_local = r; }

    if (r < r_c_local) { r_c_local = r; Dc_local = D; lapse_c_local = lapse; }

    if (r < R_core && mom > core_mom_local) core_mom_local = mom;
    if (r > surf_lo && r < surf_hi && mom > surf_mom_local) surf_mom_local = mom;
    if (mom > glob_mom_local) { glob_mom_local = mom; r_glob_mom_local = r; }

    if (r < 1.5 * R_star) {
      if (lapse < lapse_min_local) lapse_min_local = lapse;
      if (lapse > lapse_max_local) lapse_max_local = lapse;
    }
  }

  SumD_local *= app->grid.cellVolume;

  double SumD = 0.0;
  double Dmax = 0.0, core_mom = 0.0, surf_mom = 0.0, glob_mom = 0.0;
  double lapse_min = 0.0, lapse_max = 0.0, r_c = 0.0;
  int n_bad = 0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &SumD_local, &SumD);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Dmax_local, &Dmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &core_mom_local, &core_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &surf_mom_local, &surf_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &glob_mom_local, &glob_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &lapse_min_local, &lapse_min);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &lapse_max_local, &lapse_max);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_c_local, &r_c);
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_SUM, 1, &n_bad_local, &n_bad);

  double r_Dmax_candidate = 1.0e300;
  if (Dmax_local == Dmax) r_Dmax_candidate = r_Dmax_local;
  double r_glob_mom_candidate = 1.0e300;
  if (glob_mom_local == glob_mom) r_glob_mom_candidate = r_glob_mom_local;
  double Dc_candidate = 0.0, lapse_c_candidate = 0.0;
  if (r_c_local == r_c) {
    Dc_candidate = Dc_local;
    lapse_c_candidate = lapse_c_local;
  }

  double r_Dmax = 0.0, r_glob_mom = 0.0, Dc = 0.0, lapse_c = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_Dmax_candidate, &r_Dmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_glob_mom_candidate, &r_glob_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Dc_candidate, &Dc);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &lapse_c_candidate, &lapse_c);

  if (diag_SumD0 < 0.0) diag_SumD0 = SumD;
  double dSumD = SumD - diag_SumD0;

  int rank = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    printf("[DIAG] t=%.4f  Dc=%.6e lapse_c=%.6e  Dmax=%.6e @r=%.3f  core|mom|=%.3e surf|mom|=%.3e glob|mom|=%.3e @r=%.2f  lapse[%.4f,%.4f]  SumD=%.8e dSumD=%.3e bad=%d\n",
      t_curr, Dc, lapse_c, Dmax, r_Dmax, core_mom, surf_mom, glob_mom, r_glob_mom, lapse_min, lapse_max, SumD, dSumD, n_bad);
    fflush(stdout);
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

  // Stage switch (WB on for the quasi-static hold, off for the collapse): GKYL_WB_OFF disables the
  // frozen-discrete WB in ALL sites (reconstruction, flux jump, MoL source) consistently, via the eqn.
  // The frozen WB stores the FULL t=0 equilibrium (D_eq, Etot_eq, S_eq in slots 71-75), so it holds the
  // star bit-static in BOTH areal (v=0) AND Kerr-Schild (v=beta/alpha != 0) gauges - no family WB needed.
  bool wb_off = (getenv("GKYL_WB_OFF") != NULL);
  printf("WB: disable_well_balanced=%d\n", (int) wb_off);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
    .gas_gamma = ctx.gas_gamma,
    .spacetime_gauge = ctx.spacetime_gauge,
    .reinit_freq = ctx.reinit_freq,
    .spacetime = ctx.spacetime,
    .rp_type = WV_GR_EULER_RP_LAX,
    .use_gpu = app_args.use_gpu,
    .tov_eq = ctx.tov,
    .disable_well_balanced = wb_off,
    .p_atm = ctx.p_atm,
    .rho_atm = ctx.rho_atm,
  });

  // Vacuum Einstein (Bona-Masso), evolved with the MP scheme and coupled to the fluid.
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(ctx.excision_threshold,
    ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    .init = evalGREulerInit,
    .ctx = &ctx,

    .scheme_type = GKYL_MOMENT_MP, // fluid sector: MP + frozen-discrete WB

    .has_gr_euler = true,
    .gr_euler_gas_gamma = ctx.gas_gamma,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species einstein = {
    .name = "vacuum_einstein",
    .equation = vacuum_einstein,
    .init = evalVacuumEinsteinInit,
    .ctx = &ctx,

    .scheme_type = GKYL_MOMENT_WAVE_PROP, // spacetime sector: wave-propagation

    .has_vacuum_einstein = true,
    .vacuum_einstein_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_spacetime_evolution = ctx.spacetime_evolution,

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
    .name = "gr_tov_mixed_cart",

    .ndim = 2,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly },
    .upper = {  0.5 * ctx.Lx,  0.5 * ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_MIXED,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { fluid, einstein },

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
  diag_track(app, &ctx, t_curr);

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
    diag_track(app, &ctx, t_curr);

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
  gkyl_wv_eqn_release(vacuum_einstein);
  gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
  gkyl_tov_solution_release(ctx.tov);
  gkyl_comm_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Finalize();
#endif

  return 0;
}
