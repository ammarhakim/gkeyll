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
#include <gkyl_wave_prop.h>
#include <gkyl_wv_vacuum_einstein.h>
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

  // Vacuum Einstein evolution parameters.
  double excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution;

  double perturb_frac; // x-momentum kick fraction in the inner star (0 => static test).

  int Nx, Ny, Nz;
  double Lx, Ly, Lz; // Domain [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [-Lz/2,Lz/2], star centered at origin.
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

  double L = 24.0;
  const char *L_env = getenv("GKYL_LBOX");
  if (L_env != NULL) {
    L = atof(L_env);
  }

  double t_end = 20.0;
  const char *t_end_env = getenv("GKYL_T_END");
  if (t_end_env != NULL) {
    t_end = atof(t_end_env);
  }

  int num_frames = 40;
  const char *num_frames_env = getenv("GKYL_NUM_FRAMES");
  if (num_frames_env != NULL) {
    num_frames = atoi(num_frames_env);
  }

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
  enum gkyl_spacetime_evolution evolution = GKYL_RICCI_EVOLUTION;
  const char *evo = getenv("GKYL_EVOLUTION");
  if (evo && strcmp(evo, "einstein") == 0) evolution = GKYL_EINSTEIN_EVOLUTION;
  printf("Gauge: use_kerr_schild=%d, slicing=%d, evolution=%d\n",
    (int) use_kerr_schild, (int) slicing, (int) evolution);
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
    // RICCI (sigma=0): the K-flux is the standard Bona-Masso lambda^k_ij, whose
    // divergence reproduces -D_iD_j(alpha) + alpha R_ij exactly on the constraint
    // surface -- the correct partner of the Ricci-form matter source
    // -8 pi alpha (S_ij - 1/2 gamma_ij (S-E)). EINSTEIN (sigma=1) adds the
    // -(1/2) alpha V^k gamma_ij flux term, which shifts the trace sector by
    // -3/2 alpha c and, paired with the Ricci-form matter source, leaves a
    // continuum trK forcing at the density peak (the "lapse drift").
    // Env override: GKYL_EVOLUTION = ricci | einstein.
    .spacetime_evolution = evolution,
    .perturb_frac = 0.0, // set > 0 (e.g. 1.0e-3) for the momentum-perturbation test.
    .Nx = 48,
    .Ny = 48,
    .Nz = 48,
    .Lx = L,
    .Ly = L,
    .Lz = L,
    .cfl_frac = 0.8,
    .t_end = t_end,
    .num_frames = num_frames,
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
  double x = xn[0], y = xn[1], z = xn[2];
  struct gr_tov_coupled_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double r = sqrt((x * x) + (y * y) + (z * z));
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
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &spatial_metric_der);

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

  double rho_rel = sqrt_det * rho * W;
  double mom_x = sqrt_det * rho * h * (W * W) * vel[0];
  double mom_y = sqrt_det * rho * h * (W * W) * vel[1];
  double mom_z = sqrt_det * rho * h * (W * W) * vel[2];
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
  fout[68] = x; fout[69] = y; fout[70] = z;

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
  double x = xn[0], y = xn[1], z = xn[2];
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
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &spatial_metric_der);

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

// Conformal (77-comp, chi = psi^-2) vacuum-Einstein initial data, 3D. Ported from
// rt_gr_tov_conformal_mixed_cart.c with the third dimension carried through
// (z passed to every spacetime callback; 3D radius used by the fluid IC).
// Selected with GKYL_USE_CONFORMAL=1.
void
evalVacuumEinsteinConformalInit3d(double t, const double* GKYL_RESTRICT xn,
  double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
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
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &conformal_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_conformal_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->conformal_factor_func(spacetime, 0.0, x, y, z, &psi);
  spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, y, z, &chi);
  spacetime->conformal_factor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &psi_der);
  spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &chi_der);
  spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, y, z, dl2, dl2, dl2, &chi_der2);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &conformal_metric_der);

  double psi4 = psi * psi * psi * psi;
  for (int k = 0; k < 3; k++) psi_der[k] /= psi;

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

  double cmd_raised1[3][3][3], cmd_raised3[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        cmd_raised1[k][i][j] = 0.0;
        cmd_raised3[i][j][k] = 0.0;
        for (int l = 0; l < 3; l++) {
          cmd_raised1[k][i][j] += inv_conformal_metric[k][l] * conformal_metric_der[l][i][j];
          cmd_raised3[i][j][k] += inv_conformal_metric[l][k] * conformal_metric_der[i][j][l];
        }
      }

  double conformal_aux_vect[3];
  for (int i = 0; i < 3; i++) {
    conformal_aux_vect[i] = -4.0 * psi_der[i];
    for (int s = 0; s < 3; s++) {
      conformal_aux_vect[i] += cmd_raised3[i][s][s];
      conformal_aux_vect[i] -= cmd_raised1[s][s][i];
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fout[3 * i + j] = conformal_metric[i][j];
      fout[10 + 3 * i + j] = extrinsic_curvature[i][j];
    }
  }
  fout[9] = lapse;
  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        fout[19 + (9 * k) + (3 * i) + j] = conformal_metric_der[k][i][j];

  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];
  fout[49] = conformal_aux_vect[0]; fout[50] = conformal_aux_vect[1]; fout[51] = conformal_aux_vect[2];
  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      fout[55 + 3 * i + j] = shift_der[i][j];

  fout[64] = chi;
  fout[65] = chi_der[0]; fout[66] = chi_der[1]; fout[67] = chi_der[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      fout[68 + 3 * i + j] = chi_der2[i][j];

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
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);

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

static double
diag_sym3_inverse(const double g[3][3], double gi[3][3])
{
  double c00 = g[1][1] * g[2][2] - g[1][2] * g[2][1];
  double c01 = g[1][2] * g[2][0] - g[1][0] * g[2][2];
  double c02 = g[1][0] * g[2][1] - g[1][1] * g[2][0];
  double det = g[0][0] * c00 + g[0][1] * c01 + g[0][2] * c02;
  double idet = 0.0;
  if (det != 0.0) {
    idet = 1.0 / det;
  }

  gi[0][0] = c00 * idet;
  gi[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) * idet;
  gi[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) * idet;
  gi[0][1] = gi[1][0] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) * idet;
  gi[0][2] = gi[2][0] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) * idet;
  gi[1][2] = gi[2][1] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) * idet;

  return det;
}

// Set once at startup: the monitor reads the PHYSICAL metric in both formulations
// so the constraint diagnostics are directly comparable. Conformal state stores
// chi = psi^-2 at slot 64 (d_k chi at 65-67) and tilde(gamma)_ij at 0-8, with
// gamma_ij = chi^-2 tilde(gamma)_ij; K_ij is already physical at 10-18. The
// derivative reduction obeys d_k gamma_ij = chi^-2 (2 tilde(D)_kij - 2 (d_k chi/chi)
// tilde(gamma)_ij), i.e. D_kij = chi^-2 (tilde(D)_kij - (d_k chi/chi) tilde(gamma)_ij).
// (Matches sync_fluid_metric_from_einstein in mom_update_ssp_rk.c.)
static bool diag_is_conformal = false;

static void
diag_read_einstein_state(const double *q, double g[3][3], double gi[3][3],
  double K[3][3], double D[3][3][3], double A[3], double V[3])
{
  double psi4 = 1.0, dchi_over_chi[3] = { 0.0, 0.0, 0.0 };
  if (diag_is_conformal) {
    double chi = q[64];
    if (!(chi > 0.0)) chi = 1.0;
    double inv_chi = 1.0 / chi;
    psi4 = inv_chi * inv_chi;
    for (int k = 0; k < 3; k++) dchi_over_chi[k] = q[65 + k] * inv_chi;
  }

  double gtil[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gtil[i][j] = q[3 * i + j];
      g[i][j] = psi4 * gtil[i][j];
      K[i][j] = q[10 + 3 * i + j];
    }
  }
  diag_sym3_inverse(g, gi);

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double Dt = q[19 + 9 * k + 3 * i + j];
        D[k][i][j] = diag_is_conformal
          ? psi4 * (Dt - dchi_over_chi[k] * gtil[i][j])
          : Dt;
      }
    }
    A[k] = q[46 + k];
    V[k] = q[49 + k]; // Bona-Masso aux; C^V check is only meaningful non-conformally
  }
}

static void
diag_christoffel_from_state(const double *q, double Gamma[3][3][3])
{
  double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
  diag_read_einstein_state(q, g, gi, K, D, A, V);

  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      for (int c = 0; c < 3; c++) {
        Gamma[a][b][c] = 0.0;
        for (int m = 0; m < 3; m++) {
          Gamma[a][b][c] += gi[a][m] * (D[b][m][c] + D[c][m][b] - D[m][b][c]);
        }
      }
    }
  }
}

static void
diag_compute_Q(const double *q, double Q[3][3])
{
  double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
  diag_read_einstein_state(q, g, gi, K, D, A, V);

  double trK = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      trK += gi[i][j] * K[i][j];
    }
  }

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      Q[j][i] = 0.0;
      for (int a = 0; a < 3; a++) {
        Q[j][i] += gi[j][a] * K[a][i];
      }
      if (j == i) {
        Q[j][i] -= trK;
      }
    }
  }
}

static void
diag_matter_projections(double gas_gamma, const double *qf, const double g[3][3],
  double *E, double S_cov[3])
{
  double prim[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, qf, prim);
  double rho = prim[0];
  double p = prim[4];
  double vel[3] = { prim[1], prim[2], prim[3] };

  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += g[i][j] * vel[i] * vel[j];
    }
  }
  if (v_sq > 1.0 - 1.0e-12) {
    v_sq = 1.0 - 1.0e-12;
  }

  double W = 1.0 / sqrt(1.0 - v_sq);
  double h = gkyl_gr_euler_specific_enthalpy(gas_gamma, rho, p);
  *E = rho * h * W * W - p;

  for (int i = 0; i < 3; i++) {
    S_cov[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      S_cov[i] += g[i][j] * rho * h * W * W * vel[j];
    }
  }
}

// [CONS] ADM + first-order reduction monitor for the mixed 3D TOV run.
// This is diagnostic-only: it does not feed back into the update. It checks the
// quantities that would have caught the old Cartesian TOV metric bug immediately:
// C_kij = D_kij - 1/2 d_k gamma_ij, C_i = A_i - d_i log(alpha), and
// C^V_i = V_i - (D_is{}^s - D^s{}_si). It also reports ADM Hamiltonian and
// momentum residuals using the evolved first-order metric derivatives.
static void
constraint_monitor(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  static long call_count = 0;
  int every = 1;
  const char *every_env = getenv("GKYL_CONSTRAINT_DIAG_EVERY");
  if (every_env) {
    every = atoi(every_env);
    if (every < 1) {
      every = 1;
    }
  }
  call_count += 1;
  if ((call_count - 1) % every != 0) {
    return;
  }

  const struct gkyl_array *fluid = gkyl_moment_app_get_write_array_species(app, 0);
  const struct gkyl_array *einstein = gkyl_moment_app_get_write_array_species(app, 1);

  double Hmax_local = 0.0, Mmax_local = 0.0;
  double Cgmax_local = 0.0, Camax_local = 0.0, Cvmax_local = 0.0;
  double trKmax_local = 0.0, trKc_local = 0.0, r_c_local = 1.0e300;
  double r_H_local = 0.0, r_M_local = 0.0;
  int n_bad_local = 0;

  // Interior / exterior split: the global max|H| is typically dominated by the
  // star-surface kink and the outer-boundary fill error, hiding whether the CORE
  // is converging. Report the interior (r <= r_int, default just inside the star
  // surface; env GKYL_CONSTRAINT_RINT) maxima separately.
  double r_int = ctx->R_star - 1.0;
  const char *rint_env = getenv("GKYL_CONSTRAINT_RINT");
  if (rint_env) r_int = atof(rint_env);
  double Hint_local = 0.0, Mint_local = 0.0, Cgint_local = 0.0;
  // L2 (volume-weighted) interior norm of H. A max-norm is a poor convergence
  // diagnostic (it is grid-alignment sensitive and samples a single cell); the
  // L2 norm is the standard measure for a self-convergence test.
  double H2_local = 0.0, vol_local = 0.0;
  double Hnorm_int_local = 0.0, H2_noctr_local = 0.0, vol_noctr_local = 0.0;

  // [CONS-RAD] Radial profile of the constraint violation. This is the diagnostic that
  // separates the candidate causes: an ORIGIN/Cartesian-n_i representation problem loads
  // the innermost shell; a STAR-SURFACE problem loads the shell at R_star (~9.6); a
  // frozen-WB-vs-evolving-geometry (matter) problem loads the matter interior (E != 0)
  // and stops at the surface; a genuine formulation/ADM problem is spread everywhere
  // INCLUDING the vacuum exterior. Also bin the matter density so "where E lives" is
  // directly comparable to "where H lives".
  #define NRAD 8
  double rad_hi[NRAD] = { 2.0, 4.0, 6.0, 8.0, 9.6, 11.0, 14.0, 1.0e30 };
  double radH2[NRAD] = { 0.0 }, radVol[NRAD] = { 0.0 }, radE[NRAD] = { 0.0 };

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *q = gkyl_array_cfetch(einstein, loc);
    const double *qf = gkyl_array_cfetch(fluid, loc);

    double xc[3];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);

    double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
    diag_read_einstein_state(q, g, gi, K, D, A, V);

    bool bad_cell = false;
    if (!isfinite(q[9])) {
      bad_cell = true;
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (!isfinite(g[i][j]) || !isfinite(K[i][j])) {
          bad_cell = true;
        }
      }
    }
    if (bad_cell) {
      n_bad_local += 1;
      continue;
    }

    double Cg = 0.0, Ca = 0.0;
    for (int k = 0; k < 3; k++) {
      int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
      for (int d = 0; d < 3; d++) {
        idxp[d] = iter.idx[d];
        idxm[d] = iter.idx[d];
      }
      idxp[k] = iter.idx[k] + 1;
      idxm[k] = iter.idx[k] - 1;
      const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
      const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
      double dx = app->grid.dx[k];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          int gij = 3 * i + j;
          double half_dg = (qp[gij] - qm[gij]) / (4.0 * dx);
          double ckij = D[k][i][j] - half_dg;
          if (fabs(ckij) > Cg) {
            Cg = fabs(ckij);
          }
        }
      }

      double cp = qp[9], cm = qm[9];
      if (cp > 0.0 && cm > 0.0) {
        double dlog_alpha = (log(cp) - log(cm)) / (2.0 * dx);
        double ci = A[k] - dlog_alpha;
        if (fabs(ci) > Ca) {
          Ca = fabs(ci);
        }
      }
      else {
        n_bad_local += 1;
      }
    }

    double D_raised1[3][3][3], D_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          D_raised1[k][i][j] = 0.0;
          D_raised3[i][j][k] = 0.0;
          for (int l = 0; l < 3; l++) {
            D_raised1[k][i][j] += gi[k][l] * D[l][i][j];
            D_raised3[i][j][k] += gi[l][k] * D[i][j][l];
          }
        }
      }
    }

    double Cv = 0.0;
    for (int i = 0; i < 3; i++) {
      double aux_i = 0.0;
      for (int s = 0; s < 3; s++) {
        aux_i += D_raised3[i][s][s];
        aux_i -= D_raised1[s][s][i];
      }
      double cvi = V[i] - aux_i;
      if (fabs(cvi) > Cv) {
        Cv = fabs(cvi);
      }
    }

    double Gamma[3][3][3];
    diag_christoffel_from_state(q, Gamma);

    double Ricci[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Ricci[i][j] = 0.0;
        for (int k = 0; k < 3; k++) {
          int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
          for (int d = 0; d < 3; d++) {
            idxp[d] = iter.idx[d];
            idxm[d] = iter.idx[d];
          }
          idxp[k] = iter.idx[k] + 1;
          idxm[k] = iter.idx[k] - 1;
          const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
          const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
          double Gp[3][3][3], Gm[3][3][3];
          diag_christoffel_from_state(qp, Gp);
          diag_christoffel_from_state(qm, Gm);
          Ricci[i][j] += (Gp[k][i][j] - Gm[k][i][j]) / (2.0 * app->grid.dx[k]);
        }
        for (int k = 0; k < 3; k++) {
          int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
          for (int d = 0; d < 3; d++) {
            idxp[d] = iter.idx[d];
            idxm[d] = iter.idx[d];
          }
          idxp[j] = iter.idx[j] + 1;
          idxm[j] = iter.idx[j] - 1;
          const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
          const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
          double Gp[3][3][3], Gm[3][3][3];
          diag_christoffel_from_state(qp, Gp);
          diag_christoffel_from_state(qm, Gm);
          Ricci[i][j] -= (Gp[k][i][k] - Gm[k][i][k]) / (2.0 * app->grid.dx[j]);
        }
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            Ricci[i][j] += Gamma[k][i][j] * Gamma[l][k][l];
            Ricci[i][j] -= Gamma[l][i][k] * Gamma[k][j][l];
          }
        }
      }
    }

    double R = 0.0, trK = 0.0, KijKij = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        R += gi[i][j] * Ricci[i][j];
        trK += gi[i][j] * K[i][j];
      }
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double Kup = 0.0;
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            Kup += gi[i][a] * gi[j][b] * K[a][b];
          }
        }
        KijKij += K[i][j] * Kup;
      }
    }

    double E = 0.0, S_cov[3] = { 0.0, 0.0, 0.0 };
    diag_matter_projections(ctx->gas_gamma, qf, g, &E, S_cov);
    double H = R + trK * trK - KijKij - 16.0 * M_PI * E;
    // Normalized Hamiltonian: |H| divided by the sum of the magnitudes of its
    // own terms. This says whether H is a genuine O(1) imbalance or a small
    // relative residual (the standard "is the constraint violation actually big"
    // measure). And an interior L2 that EXCLUDES the innermost shell (r < 1.5 dx)
    // to separate a real distributed violation from a single center-cell FD spike.
    double Hscale = fabs(R) + trK*trK + KijKij + 16.0*M_PI*fabs(E) + 1.0e-30;
    double Hnorm = fabs(H)/Hscale;
    double dxmin = app->grid.dx[0];

    { // radial binning (see NRAD comment above)
      double cellvol = app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2];
      for (int b = 0; b < NRAD; ++b) {
        if (r <= rad_hi[b]) {
          radH2[b] += H*H*cellvol;
          radVol[b] += cellvol;
          radE[b] += fabs(E)*cellvol;
          break;
        }
      }
    }

    if (r <= r_int) {
      if (Hnorm > Hnorm_int_local) Hnorm_int_local = Hnorm;
      if (r > 1.5*dxmin) {
        double cellvol = app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2];
        H2_noctr_local += H*H*cellvol;
        vol_noctr_local += cellvol;
      }
    }

    double Mmax_cell = 0.0;
    for (int i = 0; i < 3; i++) {
      double divQ = 0.0;
      for (int j = 0; j < 3; j++) {
        int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
        for (int d = 0; d < 3; d++) {
          idxp[d] = iter.idx[d];
          idxm[d] = iter.idx[d];
        }
        idxp[j] = iter.idx[j] + 1;
        idxm[j] = iter.idx[j] - 1;
        const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
        const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
        double Qp[3][3], Qm[3][3];
        diag_compute_Q(qp, Qp);
        diag_compute_Q(qm, Qm);
        divQ += (Qp[j][i] - Qm[j][i]) / (2.0 * app->grid.dx[j]);
      }

      double Q[3][3];
      diag_compute_Q(q, Q);
      for (int j = 0; j < 3; j++) {
        for (int l = 0; l < 3; l++) {
          divQ += Gamma[j][j][l] * Q[l][i];
          divQ -= Gamma[l][j][i] * Q[j][l];
        }
      }

      double Mi = divQ - 8.0 * M_PI * S_cov[i];
      if (fabs(Mi) > Mmax_cell) {
        Mmax_cell = fabs(Mi);
      }
    }

    if (fabs(H) > Hmax_local) {
      Hmax_local = fabs(H);
      r_H_local = r;
    }
    if (Mmax_cell > Mmax_local) {
      Mmax_local = Mmax_cell;
      r_M_local = r;
    }
    if (r <= r_int) {
      if (fabs(H) > Hint_local) Hint_local = fabs(H);
      if (Mmax_cell > Mint_local) Mint_local = Mmax_cell;
      if (Cg > Cgint_local) Cgint_local = Cg;
      double cellvol = app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2];
      H2_local += H*H*cellvol;
      vol_local += cellvol;
    }
    if (Cg > Cgmax_local) Cgmax_local = Cg;
    if (Ca > Camax_local) Camax_local = Ca;
    if (Cv > Cvmax_local) Cvmax_local = Cv;
    if (fabs(trK) > fabs(trKmax_local)) trKmax_local = trK;
    if (r < r_c_local) {
      r_c_local = r;
      trKc_local = trK;
    }

    if (!isfinite(H) || !isfinite(Mmax_cell) || !isfinite(Cg) || !isfinite(Ca) || !isfinite(Cv)) {
      n_bad_local += 1;
    }
  }

  double Hmax = 0.0, Mmax = 0.0, Cgmax = 0.0, Camax = 0.0, Cvmax = 0.0, abs_trKmax = 0.0, r_c = 0.0;
  int n_bad = 0;
  double abs_trKmax_local = fabs(trKmax_local);
  double Hint = 0.0, Mint = 0.0, Cgint = 0.0, H2sum = 0.0, volsum = 0.0;
  double Hnorm_int = 0.0, H2noctr = 0.0, volnoctr = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hnorm_int_local, &Hnorm_int);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H2_noctr_local, &H2noctr);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &vol_noctr_local, &volnoctr);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H2_local, &H2sum);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &vol_local, &volsum);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hint_local, &Hint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Mint_local, &Mint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cgint_local, &Cgint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hmax_local, &Hmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Mmax_local, &Mmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cgmax_local, &Cgmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Camax_local, &Camax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cvmax_local, &Cvmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &abs_trKmax_local, &abs_trKmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_c_local, &r_c);
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_SUM, 1, &n_bad_local, &n_bad);

  double r_H_candidate = 1.0e300;
  if (Hmax_local == Hmax) r_H_candidate = r_H_local;
  double r_M_candidate = 1.0e300;
  if (Mmax_local == Mmax) r_M_candidate = r_M_local;
  double trKmax_candidate = -1.0e300;
  if (fabs(trKmax_local) == abs_trKmax) trKmax_candidate = trKmax_local;
  double trKc_candidate = -1.0e300;
  if (r_c_local == r_c) trKc_candidate = trKc_local;

  double r_H = 0.0, r_M = 0.0, trKmax = 0.0, trKc = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_H_candidate, &r_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_M_candidate, &r_M);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &trKmax_candidate, &trKmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &trKc_candidate, &trKc);

  int rank = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    printf("[CONS] t=%.4f  max|H|=%.6e @r=%.3f  max|M|=%.6e @r=%.3f  max|Cgamma|=%.6e max|CA|=%.6e max|CV|=%.6e  trK_c=%.6e max|trK|=%.6e  bad=%d\n",
      t_curr, Hmax, r_H, Mmax, r_M, Cgmax, Camax, Cvmax, trKc, trKmax, n_bad);
    double H_l2 = (volsum > 0.0) ? sqrt(H2sum/volsum) : 0.0;
    double H_l2_noctr = (volnoctr > 0.0) ? sqrt(H2noctr/volnoctr) : 0.0;
    printf("[CONS-INT] t=%.4f  r<=%.1f: max|H|=%.6e L2|H|=%.6e L2|H|(no-ctr)=%.6e maxHnorm=%.4f max|M|=%.6e max|Cgamma|=%.6e\n",
      t_curr, r_int, Hint, H_l2, H_l2_noctr, Hnorm_int, Mint, Cgint);
    printf("[CONS-RAD] t=%.4f  L2|H| per shell:", t_curr);
    double rlo = 0.0;
    for (int b = 0; b < NRAD; ++b) {
      double hb = (radVol[b] > 0.0) ? sqrt(radH2[b]/radVol[b]) : 0.0;
      double eb = (radVol[b] > 0.0) ? radE[b]/radVol[b] : 0.0;
      double rhi = (rad_hi[b] > 1.0e29) ? 99.0 : rad_hi[b];
      printf("  [%.0f-%.0f]=%.3e(E=%.1e)", rlo, rhi, hb, eb);
      rlo = rhi;
    }
    printf("\n");
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
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], ctx.Nz);

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
  // GKYL_USE_CONFORMAL=1 selects the conformal (77-comp, chi) formulation instead of
  // the non-conformal (64-comp) one, so the two can be compared under an identical
  // fluid sector, grid, gauge and constraint monitor.
  // Spacetime-sector limiter override. The default (GKYL_MONOTONIZED_CENTERED) is TVD
  // and therefore clips at extrema -- including the density peak at the star centre,
  // which is exactly where we measure the flux/source equilibrium residual. If that
  // residual is O(h) because of limiter clipping rather than the scheme, switching to
  // GKYL_NO_LIMITER should restore O(h^2). GKYL_ST_LIMITER = none|minmod|mc|vanleer.
  enum gkyl_wave_limiter st_limiter = GKYL_MONOTONIZED_CENTERED;
  const char *stlim = getenv("GKYL_ST_LIMITER");
  if (stlim) {
    if (strcmp(stlim, "none") == 0) st_limiter = GKYL_NO_LIMITER;
    else if (strcmp(stlim, "minmod") == 0) st_limiter = GKYL_MIN_MOD;
    else if (strcmp(stlim, "vanleer") == 0) st_limiter = GKYL_VAN_LEER;
  }
  printf("Spacetime limiter: %d (1=none,2=minmod,5=MC)\n", (int) st_limiter);

  bool use_conformal = (getenv("GKYL_USE_CONFORMAL") != NULL);
  diag_is_conformal = use_conformal;
  printf("Formulation: %s\n", use_conformal ? "CONFORMAL (77)" : "NON-CONFORMAL (64)");

  struct gkyl_wv_eqn *vacuum_einstein = use_conformal
    ? gkyl_wv_vacuum_einstein_conformal_new(ctx.excision_threshold,
        ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu)
    : gkyl_wv_vacuum_einstein_new(ctx.excision_threshold,
        ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);

  enum gkyl_species_bc_type einstein_bc = GKYL_SPECIES_RADIAL_FALLOFF;
  const char *sommerfeld_env = getenv("GKYL_USE_FALLOFF_BC");
  if (sommerfeld_env != NULL) {
    if (strcmp(sommerfeld_env, "0") == 0) {
      einstein_bc = GKYL_SPECIES_COPY;
    }
  }
  const char *einstein_bc_name = "Sommerfeld";
  if (einstein_bc == GKYL_SPECIES_COPY) {
    einstein_bc_name = "Copy";
  }
  printf("Einstein BC: %s\n", einstein_bc_name);

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
    .bcz = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species einstein = {
    .name = "vacuum_einstein",
    .equation = vacuum_einstein,
    .init = use_conformal ? evalVacuumEinsteinConformalInit3d : evalVacuumEinsteinInit,
    .ctx = &ctx,

    .scheme_type = GKYL_MOMENT_WAVE_PROP, // spacetime sector: wave-propagation
    .limiter = st_limiter,

    .has_vacuum_einstein = !use_conformal,
    .vacuum_einstein_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_spacetime_evolution = ctx.spacetime_evolution,

    .has_vacuum_einstein_conformal = use_conformal,
    .vacuum_einstein_conformal_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_conformal_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_conformal_spacetime_evolution = ctx.spacetime_evolution,

    .bcx = { einstein_bc, einstein_bc },
    .bcy = { einstein_bc, einstein_bc },
    .bcz = { einstein_bc, einstein_bc },
  };

  int cells[] = { NX, NY, NZ };
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
    .name = "gr_tov_mixed_cart3d",

    .ndim = 3,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly, -0.5 * ctx.Lz },
    .upper = {  0.5 * ctx.Lx,  0.5 * ctx.Ly,  0.5 * ctx.Lz },
    .cells = { NX, NY, NZ },

    .scheme_type = GKYL_MOMENT_MIXED,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { fluid, einstein },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { cuts[0], cuts[1], cuts[2] },
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
  constraint_monitor(app, &ctx, t_curr);

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
    constraint_monitor(app, &ctx, t_curr);

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
