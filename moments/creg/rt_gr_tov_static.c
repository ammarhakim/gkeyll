// 2D static TOV neutron star, for the general relativistic Euler equations.
// Uses gkyl_gr_tov_spacetime which provides the TOV interior metric for r < R_star
// and delegates to Schwarzschild (M = M_star) for r >= R_star.  No excision.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_moment_priv.h>
#include <gkyl_util.h>
//#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_polytropic_euler.h>
#include <gkyl_wv_gr_polytropic_euler_priv.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>
#include "tov_solver.h"
#include "gkyl_gr_tov_spacetime.h"

struct tov_static_ctx
{
  double pi;
  double gas_gamma;

  double K_poly;
  double rho_c;
  double dr_tov;

  double rho_atm;
  double p_atm;

  double pos_x;
  double pos_y;
  double pos_z;

  struct gkyl_tov *tov;
  struct gkyl_gr_spacetime *spacetime;

  //struct gkyl_gr_spacetime *bh_test;

  double M_star;
  double R_star;

  int Nx;
  int Ny;
  int Nz;
  double Lx;
  double Ly;
  double Lz;
  double cfl_frac;

  enum gkyl_spacetime_gauge spacetime_gauge;
  int reinit_freq;

  double t_end;
  int num_frames;
  int field_energy_calcs;
  int integrated_mom_calcs;
  double dt_failure_tol;
  int num_failures_max;
};

struct tov_static_ctx
create_ctx(void)
{
  double pi = M_PI;

  double gas_gamma = 2.0;
  double K_poly = 100.0; 

  /*Einstein Toolkit Linear TOV*/
  /*
  double rho_c = 1.28e-3;
  double dr_tov = 0.01;

  double rho_atm = 1e-2 * rho_c;
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);

  double M_star = gkyl_tov_star_mass(tov);
  double R_star = gkyl_tov_star_radius(tov);

  double Lx = 640.0 * M_star;
  double Ly = 640.0 * M_star;
  double Lz = 640.0 * M_star;
  double pos_x = 0.5 * Lx;
  double pos_y = 0.5 * Ly;
  double pos_z = 0.5 * Lz;

  // TOV spacetime: TOV interior for r < R_star, Schwarzschild for r >= R_star, no excision.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, pos_x, 0.0, 0.0);
  //struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  int Nx = 320;
  int Ny = 1;//320;
  int Nz = 1;//64;
  double cfl_frac = 0.5;
  */
  
  double rho_c = 5e-4;
  double dr_tov = 0.01;

  double rho_atm = 1e-2 * rho_c;
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);

  double M_star = gkyl_tov_star_mass(tov);
  printf("M_star = %e \n", M_star);
  double R_star = gkyl_tov_star_radius(tov);

  double Lx = 100.0;
  double Ly = 1.0;
  double Lz = 1.0;
  double pos_x = 0.5 * Lx;
  double pos_y = 0.0;
  double pos_z = 0.0;

  // TOV spacetime: TOV interior for r < R_star, Schwarzschild for r >= R_star, no excision.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, pos_x, 0.0, 0.0);
  //struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  int Nx = 512;
  const char *nx_env = getenv("GKYL_TOV_NX");
  if (nx_env) {
    int nx_in = atoi(nx_env);
    if (nx_in > 0) {
      Nx = nx_in;
    }
  }
  int Ny = 1;
  int Nz = 1;
  double cfl_frac = 0.5;
  

  enum gkyl_spacetime_gauge spacetime_gauge = GKYL_STATIC_GAUGE;
  int reinit_freq = INT_MAX;
  //printf("%i\n",INT_MAX);

  double t_end = 1000.0;
  int num_frames = 1000;
  int field_energy_calcs = INT_MAX;
  int integrated_mom_calcs = INT_MAX;
  double dt_failure_tol = 1.0e-4;
  int num_failures_max = 20;

  struct tov_static_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .K_poly    = K_poly,
    .rho_c     = rho_c,
    .dr_tov    = dr_tov,
    .rho_atm   = rho_atm,
    .p_atm     = p_atm,
    .pos_x     = pos_x,
    .pos_y     = pos_y,
    .pos_z     = pos_z,
    .tov       = tov,
    .spacetime = spacetime,
    //.bh_test = bh_test,
    .M_star    = M_star,
    .R_star    = R_star,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .cfl_frac = cfl_frac,
    .spacetime_gauge = spacetime_gauge,
    .reinit_freq     = reinit_freq,
    .t_end           = t_end,
    .num_frames      = num_frames,
    .field_energy_calcs   = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol   = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalGRPolytropicEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  
  double x = xn[0];
  struct tov_static_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double K_poly    = app->K_poly;
  double rho_atm   = app->rho_atm;
  double p_atm     = app->p_atm;
  double R_star    = app->R_star;

  struct gkyl_tov *tov = app->tov;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double pos_x = app->pos_x;
  double pos_y = app->pos_y;
  double pos_z = app->pos_z;

  // static int init_check = 0;
  // if (init_check < 3) {
  //   printf("[INIT] x=%.4f pos_x=%.4f dx=%.4f r=%.4f R_star=%.4f\n",
  //       x, pos_x, x-pos_x, fabs(x-pos_x), R_star);
  //   init_check++;
  // }

  double dx = x - pos_x;
  double r  = fabs(dx);

  struct tov_eval_fluid fluid = {0}; 

  double rho = 0.0;
  double p = 0.0;


  if (r <= R_star) {
    gkyl_tov_eval(tov, dx, 0.0, 0.0, &fluid);
    rho = fluid.rho;
    p = fluid.P;

    if (rho < rho_atm) {
      rho = rho_atm;
    }
    if (p < p_atm) {
      p = p_atm;
    }  
  }
  else {
    rho = rho_atm;
    p = p_atm;
  }

  bool is_atmosphere = (rho < rho_atm) || (p < p_atm) || (r > R_star);

  // Metric calc
  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++)
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++)
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++)
    shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++)
      spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->lapse_function_func(spacetime, 0.0, x, 0.0, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, 0.0, 0.0, &shift);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_metric);
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_det);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

  // double vel[3] = { shift[0]/lapse, shift[1]/lapse, shift[2]/lapse };
  double vel[3] = {0.0, 0.0, 0.0};  // w^x = 0 for static star
  double W = 1.0;  // since w=0 means v=β/α which gives the correct W

  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  // double W = 1.0 / (sqrt(1.0 - v_sq));
  // if (v_sq > 1.0 - pow(10.0, -8.0)) {
  //   W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  // }

  double eps = is_atmosphere ? (p / ((gas_gamma - 1.0) * rho)) : fluid.eps;
  double h = 1.0 + eps + p / rho;
  //double h = 1.0 + (gas_gamma / (gas_gamma - 1.0)) * (p / fmax(rho, 1.0e-15));
  // double cov_vel[3] = { 0.0, 0.0, 0.0 };
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     cov_vel[i] += spatial_metric[i][j] * vel[j];
  //   }
  // }

  double Etot = sqrt(spatial_det) * (rho * h * W * W - p); // Fluid total energy density.
  double mom_x = 0.0;
  double mom_y = 0.0;
  double mom_z = 0.0;
  // double mom_x = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[0]; // Fluid momentum density (x-direction).
  // double mom_y = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[1];
  // double mom_z = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[2];

  // Complement to inside Lapse velocity diagnostic
  // printf("r = %.6e vel[0]=%.6e, vel[1]=%.6e, vel[2]=%.6e\n", sqrt(dx*dx), vel[0], vel[1], vel[2]);
  // printf("v2=%.6e  W=%.6f (should be >1 inside star)\n", v_sq, 1.0/sqrt(1.0-v_sq));

  //Fluid diagnostic (from TOV)
  // double M_star = gkyl_tov_star_mass(tov);
  // //printf("mstar = %e", M_star);
  // R_star = gkyl_tov_star_radius(tov);

  //printf("FROM INIT r=%.6e  lapse=%.6f beta[0]=%.6f spatial_metric[0][0]=%.6f\n", r, lapse, shift[0], spatial_metric[0][0]);
  // r=10.772 inside
  // r = 11.10260 outside
  //static int tov_debug_printed = 0;
  //if (!tov_debug_printed) {

    // printf("[TOV] M_star=%.6e R_star=%.6e\n", M_star, R_star);
    // // compactness 0.1-0.2 for a realistic NS
    // printf("[TOV] compactness 2M/R = %.6f\n", 2.0 * M_star / R_star);

    // struct tov_eval_fluid out;
    // gkyl_tov_eval(tov, 0.0, 0.0, 0.0, &out);
    // printf("[CKS center] rho=%.6e P=%.6e eps=%.6e\n vel=%.6e", out.rho, out.P, out.eps, vel[0]);

    // gkyl_tov_eval(tov, 0.5 * R_star, 0.0, 0.0, &out);
    // printf("[CKS mid] rho=%.6e P=%.6e\n", out.rho, out.P);

    // gkyl_tov_eval(tov, 1.1 * R_star, 0.0, 0.0, &out);
    // printf("[CKS exterior] rho=%.6e P=%.6e (should be ~0)\n", out.rho, out.P);
    //tov_debug_printed = 1;
  //}

  

  // // BH exterior at r_test along x-axis
  // double lapse_ext, *beta_ext;
  // beta_ext = gkyl_malloc(sizeof(double[3]));

  // R_star = gkyl_tov_star_radius(tov);
  // double M_star = gkyl_tov_star_mass(tov);
  // double r_test = R_star + 0.1;
  // struct gkyl_gr_spacetime *bh = gkyl_gr_blackhole_new(false, M_star, 0.0, pos_x, 0.0, 0.0);
  

  // bh->lapse_function_func(bh, 0.0, r_test + pos_x, 0.0, 0.0, &lapse_ext);
  // bh->shift_vector_func(bh, 0.0, r_test + pos_x, 0.0, 0.0, &beta_ext);

  // // Expected values
  // double V = 2.0 * M_star / r_test;
  // double alpha_expected = 1.0 / sqrt(1.0 + V);
  // double beta_expected  = V / (1.0 + V);

  // printf("=== SURFACE CONTINUITY CHECK ===\n");
  // printf("R_star=%.6f  M_star=%.6f\n", R_star, M_star);
  // printf("EXT at r=R*+0.1: alpha=%.8f (expected %.8f)\n", lapse_ext, alpha_expected);
  // printf("EXT at r=R*+0.1: beta^x=%.8f (expected %.8f)\n", beta_ext[0], beta_expected);

  // // INT at r_test - 0.2 (just inside)
  // double lapse_int;
  // double *beta_int = gkyl_malloc(sizeof(double[3]));
  // double **gamma_int = gkyl_malloc(sizeof(double*[3]));
  // for (int i = 0; i < 3; i++)
  //   gamma_int[i] = gkyl_malloc(sizeof(double[3]));

  // //eval_cks_interior(tov, R_star - 0.1, 0.0, 0.0, &lapse_int, beta_int, gamma_int);
  // spacetime->lapse_function_func(spacetime, 0.0, r_test  + pos_x, 0.0, 0.0, &lapse_int);
  // spacetime->shift_vector_func(spacetime, 0.0, r_test  + pos_x, 0.0, 0.0, &beta_int);
  // spacetime->spatial_metric_tensor_func(spacetime, 0.0, r_test  + pos_x, 0.0, 0.0, &gamma_int);
  // printf("INT at r=R*-0.1: alpha=%.8f\n", lapse_int);
  // printf("INT at r=R*-0.1: beta^x=%.8f\n", beta_int[0]);
  // printf("================================\n");
  // gkyl_free(beta_ext);

  // // Check extrinsitc cuvrature
  // R_star = gkyl_tov_star_radius(tov);
  // double M_star = gkyl_tov_star_mass(tov);
  // double r_test = R_star + 0.1;
  
  // struct gkyl_gr_spacetime *bh = gkyl_gr_blackhole_new(false, M_star, 0.0, pos_x, 0.0, 0.0);

  // double **extrinsic_curvature_test = gkyl_malloc(sizeof(double*[3]));
  // for (int i = 0; i < 3; i++)
  //   extrinsic_curvature_test[i] = gkyl_malloc(sizeof(double[3]));

  // spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.5, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("TOV at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.5, **extrinsic_curvature_test);
  // spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.3, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("TOV at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.3, **extrinsic_curvature_test);
  // spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.2, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("TOV at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.2, **extrinsic_curvature_test);
  // spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.1, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("TOV at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.1, **extrinsic_curvature_test);
  
  // bh->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.5, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("BH at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.5, **extrinsic_curvature_test);
  // bh->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.3, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("BH at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.3, **extrinsic_curvature_test);
  // bh->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.2, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("BH at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.2, **extrinsic_curvature_test);
  // bh->extrinsic_curvature_tensor_func(spacetime, 0.0, r_test - 0.1 + pos_x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature_test);
  // printf("BH at r=%6e: extrinsic_curvature_test=%.8f\n", r_test - 0.1, **extrinsic_curvature_test);

  // printf("================================\n");


  // Pack output

  // Set fluid energy density
  fout[0] = Etot;
  // Set fluid momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;

  // Set lapse gauge variable.
  fout[4] = lapse;
  // Set shift gauge variables.
  fout[5] = shift[0]; fout[6] = shift[1]; fout[7] = shift[2];

  // Set spatial metric tensor.
  fout[8] = spatial_metric[0][0]; fout[9] = spatial_metric[0][1]; fout[10] = spatial_metric[0][2];
  fout[11] = spatial_metric[1][0]; fout[12] = spatial_metric[1][1]; fout[13] = spatial_metric[1][2];
  fout[14] = spatial_metric[2][0]; fout[15] = spatial_metric[2][1]; fout[16] = spatial_metric[2][2];

  // Set extrinsic curvature tensor.
  fout[17] = extrinsic_curvature[0][0]; fout[18] = extrinsic_curvature[0][1]; fout[19] = extrinsic_curvature[0][2];
  fout[20] = extrinsic_curvature[1][0]; fout[21] = extrinsic_curvature[1][1]; fout[22] = extrinsic_curvature[1][2];
  fout[23] = extrinsic_curvature[2][0]; fout[24] = extrinsic_curvature[2][1]; fout[25] = extrinsic_curvature[2][2];

  // Set excision boundary conditions.
  fout[26] = 1.0;

  // Set lapse function derivatives.
  fout[27] = lapse_der[0]; fout[28] = lapse_der[1]; fout[29] = lapse_der[2];
  // Set shift vector derivatives.
  fout[30] = shift_der[0][0]; fout[31] = shift_der[0][1]; fout[32] = shift_der[0][2];
  fout[33] = shift_der[1][0]; fout[34] = shift_der[1][1]; fout[35] = shift_der[1][2];
  fout[36] = shift_der[2][0]; fout[37] = shift_der[2][1]; fout[38] = shift_der[2][2];

  // Set spatial metric tensor derivatives.
  fout[39] = spatial_metric_der[0][0][0]; fout[40] = spatial_metric_der[0][0][1]; fout[41] = spatial_metric_der[0][0][2];
  fout[42] = spatial_metric_der[0][1][0]; fout[43] = spatial_metric_der[0][1][1]; fout[44] = spatial_metric_der[0][1][2];
  fout[45] = spatial_metric_der[0][2][0]; fout[46] = spatial_metric_der[0][2][1]; fout[47] = spatial_metric_der[0][2][2];

  fout[48] = spatial_metric_der[1][0][0]; fout[49] = spatial_metric_der[1][0][1]; fout[50] = spatial_metric_der[1][0][2];
  fout[51] = spatial_metric_der[1][1][0]; fout[52] = spatial_metric_der[1][1][1]; fout[53] = spatial_metric_der[1][1][2];
  fout[54] = spatial_metric_der[1][2][0]; fout[55] = spatial_metric_der[1][2][1]; fout[56] = spatial_metric_der[1][2][2];

  fout[57] = spatial_metric_der[2][0][0]; fout[58] = spatial_metric_der[2][0][1]; fout[59] = spatial_metric_der[2][0][2];
  fout[60] = spatial_metric_der[2][1][0]; fout[61] = spatial_metric_der[2][1][1]; fout[62] = spatial_metric_der[2][1][2];
  fout[63] = spatial_metric_der[2][2][0]; fout[64] = spatial_metric_der[2][2][1]; fout[65] = spatial_metric_der[2][2][2];

  // Set evolution parameter.
  fout[66] = 0.0;

  // Set spatial coordinates.
  fout[67] = x; fout[68] = 0.0; fout[69] = 0.0;
  // static int init_coord_check = 0;
  // if (init_coord_check < 5) {
  //   printf("[INIT COORD] x=%.6f fout[67]=%.6f fout[68]=%.6f\n", x, fout[67], fout[68]);
  //   init_coord_check++;
  // }
  

  // double v_check[70] = {0.0};
  // gkyl_gr_polytropic_euler_prim_vars(K_poly, gas_gamma, fout, v_check);

  // static int check_count = 0;
  // if (check_count < 5 && rho > 1e-6) {
  //   printf("[init] x=%.4f fout[67]=%.4f fout[68]=%.4f fout[69]=%.4f\n",
  //   x, fout[67], fout[68], fout[69]);
  //   // printf("INIT rho=%.8e, RECOVERED rho=%.8e, ratio=%.8e\n", rho, v_check[0], v_check[0]/rho);
  //   // printf("INIT p=%.8e, RECOVERED p=%.8e\n", p, K_poly * pow(v_check[0], gas_gamma));
  //   // printf("INIT Etot=%.8e, q[0]=%.8e, spatial_det=%.8e\n", Etot, fout[0], spatial_det);
  //   check_count++;
  // }

  // Free all tensorial quantities.
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
    gkyl_free(shift_der[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_metric_der[i][j]);
    }
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(extrinsic_curvature);
  gkyl_free(shift);
  gkyl_free(lapse_der);
  gkyl_free(shift_der);
  gkyl_free(spatial_metric_der);
}

static void
build_shifted_static_state(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  evalGRPolytropicEulerInit(t, xn, fout, ctx);

  struct tov_static_ctx *app = ctx;
  double x = xn[0];
  double dx = x - app->pos_x;
  double r = fabs(dx);

  if (r > app->R_star) {
    return;
  }

  double rho = 0.0, p = 0.0, eps = 0.0;
  struct tov_eval_fluid fluid = {0};
  gkyl_tov_eval(app->tov, dx, 0.0, 0.0, &fluid);
  rho = fmax(fluid.rho, app->rho_atm);
  p = fmax(fluid.P, app->p_atm);
  eps = fluid.eps;

  double lapse = fout[4];
  double shift[3] = { fout[5], fout[6], fout[7] };
  double spatial_metric[3][3];
  spatial_metric[0][0] = fout[8];  spatial_metric[0][1] = fout[9];  spatial_metric[0][2] = fout[10];
  spatial_metric[1][0] = fout[11]; spatial_metric[1][1] = fout[12]; spatial_metric[1][2] = fout[13];
  spatial_metric[2][0] = fout[14]; spatial_metric[2][1] = fout[15]; spatial_metric[2][2] = fout[16];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  // double vel[3] = {
  //   shift[0] / lapse,
  //   shift[1] / lapse,
  //   shift[2] / lapse
  // };
  double vel[3] = {0,0,0};


  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }
  //double W = 1.0 / sqrt(fmax(1.0 - v_sq, 1.0e-8));
  double W = 1.0;
  double h = 1.0 + eps + p / rho;
  double cov_vel[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cov_vel[i] += spatial_metric[i][j] * vel[j];
    }
  }

  fout[0] = sqrt(spatial_det) * (((rho * h) * (W * W)) - p);
  fout[1] = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[0];
  fout[2] = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[1];
  fout[3] = sqrt(spatial_det) * (rho * h) * (W * W) * cov_vel[2];
}

static inline void
calc_jump_diag(int n, const double *ql, const double *qr, double *jump)
{
  for (int d = 0; d < n; d++) {
    jump[d] = qr[d] - ql[d];
  }
}

static inline double
nearest_cell_center_x(const struct tov_static_ctx *ctx, double x)
{
  double dx = ctx->Lx / ctx->Nx;
  int ci = (int) floor(x / dx);

  if (ci < 0) {
    ci = 0;
  }
  if (ci > ctx->Nx - 1) {
    ci = ctx->Nx - 1;
  }

  return (ci + 0.5) * dx;
}

static void
compute_low_order_interface_fluct(const struct gkyl_wv_eqn *eqn, const double ql[70], const double qr[70],
  double amdq[70], double apdq[70])
{
  const double tau1[3] = { 0.0, 1.0, 0.0 };
  const double tau2[3] = { 0.0, 0.0, 1.0 };
  const double norm[3] = { 1.0, 0.0, 0.0 };
  const double phil = 1.0, phir = 1.0;

  double ql_local[70] = {0.0}, qr_local[70] = {0.0};
  double delta[70] = {0.0};
  double waves[140] = {0.0};
  double speeds[2] = {0.0};
  double amdq_local[70] = {0.0}, apdq_local[70] = {0.0};

  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, ql, ql_local);
  gkyl_wv_eqn_rotate_to_local(eqn, tau1, tau2, norm, qr, qr_local);

  calc_jump_diag(70, ql_local, qr_local, delta);
  gkyl_wv_eqn_waves(eqn, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, phil, phir, waves, speeds);
  gkyl_wv_eqn_qfluct(eqn, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, phil, phir, waves, speeds, amdq_local, apdq_local);

  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, amdq_local, amdq);
  gkyl_wv_eqn_rotate_to_global(eqn, tau1, tau2, norm, apdq_local, apdq);
}

static void
compute_low_order_fv_smx_divergence(const struct gkyl_wv_eqn *eqn, const double q_im1[70], const double q_i[70],
  const double q_ip1[70], double dx, double *div_smx_out)
{
  double amdq_l[70] = {0.0}, apdq_l[70] = {0.0};
  double amdq_r[70] = {0.0}, apdq_r[70] = {0.0};

  compute_low_order_interface_fluct(eqn, q_im1, q_i, amdq_l, apdq_l);
  compute_low_order_interface_fluct(eqn, q_i, q_ip1, amdq_r, apdq_r);

  *div_smx_out = (apdq_l[1] + amdq_r[1]) / dx;
}

static void
print_tov_mass_density_diag(gkyl_moment_app* app, const struct tov_static_ctx* ctx, double t_curr)
{
  struct gkyl_array *f = app->species[0].f[0];
  const struct gkyl_range *rng = &app->local;

  double x_lo = app->grid.lower[0];
  double dx = (app->grid.upper[0] - app->grid.lower[0]) / app->grid.cells[0];
  double dy = 1.0, dz = 1.0;
  if (app->ndim > 1) {
    dy = (app->grid.upper[1] - app->grid.lower[1]) / app->grid.cells[1];
  }
  if (app->ndim > 2) {
    dz = (app->grid.upper[2] - app->grid.lower[2]) / app->grid.cells[2];
  }
  double cell_volume = dx * dy * dz;

  double D_total_domain = 0.0;
  double D_total_star = 0.0;
  int idx[GKYL_MAX_DIM] = {0};

  for (int i = rng->lower[0]; i <= rng->upper[0]; i++) {
    idx[0] = i;
    const double *q = gkyl_array_cfetch(f, gkyl_range_idx(rng, idx));

    if (q[26] < 1.0e-8) {
      continue;
    }

    double v[70] = {0.0};
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q, v);

    double rho = v[0];
    double vx = v[1], vy = v[2], vz = v[3];

    double spatial_det =
        v[8]  * (v[12] * v[16] - v[13] * v[15])
      - v[9]  * (v[11] * v[16] - v[13] * v[14])
      + v[10] * (v[11] * v[15] - v[12] * v[14]);

    double v_sq = 0.0;
    v_sq += v[8]  * vx * vx + v[9]  * vx * vy + v[10] * vx * vz;
    v_sq += v[11] * vy * vx + v[12] * vy * vy + v[13] * vy * vz;
    v_sq += v[14] * vz * vx + v[15] * vz * vy + v[16] * vz * vz;

    double W = 1.0 / sqrt(fmax(1.0 - v_sq, 1.0e-14));
    double D_cell = sqrt(fmax(spatial_det, 0.0)) * rho * W * cell_volume;
    double x = q[67], y = q[68], z = q[69];
    double r = sqrt((x - ctx->pos_x) * (x - ctx->pos_x)
      + (y - ctx->pos_y) * (y - ctx->pos_y)
      + (z - ctx->pos_z) * (z - ctx->pos_z));

    D_total_domain += D_cell;
    if (r <= ctx->R_star + 0.5 * dx) {
      D_total_star += D_cell;
    }
  }

  int ic = (int) ((ctx->pos_x - x_lo) / dx) + rng->lower[0];
  if (ic < rng->lower[0]) {
    ic = rng->lower[0];
  }
  if (ic > rng->upper[0]) {
    ic = rng->upper[0];
  }

  idx[0] = ic;
  const double *q_center = gkyl_array_cfetch(f, gkyl_range_idx(rng, idx));
  double v_center[70] = {0.0};
  gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_center, v_center);
  double rho_c_num = v_center[0];

  struct tov_eval_fluid fl = {0};
  gkyl_tov_eval(ctx->tov, 0.0, 0.0, 0.0, &fl);
  double rho_c_exact = fl.rho;
  double rho_rel = fabs(rho_c_num - rho_c_exact) / fmax(rho_c_exact, 1.0e-30);

  printf("t=%.4f  rho_c=%.10e  rho_c_exact=%.10e  rho_rel=%.3e  D_domain=%.10e  D_star=%.10e\n",
    t_curr, rho_c_num, rho_c_exact, rho_rel, D_total_domain, D_total_star);
}

static void
print_tov_force_balance_diag(gkyl_moment_app* app, const struct tov_static_ctx* ctx)
{
  struct gkyl_array *f = app->species[0].f[0];
  const struct gkyl_range *rng = &app->local;
  const struct gkyl_wv_eqn *eqn = app->species[0].equation;

  double x_lo = app->grid.lower[0];
  double dx = (app->grid.upper[0] - app->grid.lower[0]) / app->grid.cells[0];

  // Probe the center and a few interior offsets so we can see whether the
  // force balance is good throughout the star, not just at one location.
  int ic = (int) ((ctx->pos_x - x_lo) / dx) + rng->lower[0];
  int sample_offsets[] = { 0, 8, 16, 24, 32 };
  int nsamp = sizeof(sample_offsets) / sizeof(sample_offsets[0]);

  printf("\n=== TOV force-balance diagnostic at t=0 ===\n");
  printf("# dP_dr      : centered pressure gradient proxy from neighboring cells\n");
  printf("# -divF_momx : actual low-order FV x-momentum flux divergence\n");
  printf("# src_raw_x  : momentum source before multiplying by lapse*sqrt(det(gamma))\n");
  printf("# src_rhs_x  : conservative RHS source term = sqrt(det(gamma))*lapse*src_raw_x\n");
  printf("# ratioP     : dP_dr / src_raw_x, should be near -1 for static balance\n");
  printf("# ratioFV    : (-divF_momx) / src_rhs_x, should be near -1 for discrete balance\n");

  for (int s = 0; s < nsamp; s++) {
    int i = ic + sample_offsets[s];
    if (i <= rng->lower[0] || i >= rng->upper[0]) {
      continue;
    }

    int idx_m[GKYL_MAX_DIM] = {0};
    int idx_0[GKYL_MAX_DIM] = {0};
    int idx_p[GKYL_MAX_DIM] = {0};
    idx_m[0] = i - 1;
    idx_0[0] = i;
    idx_p[0] = i + 1;

    const double *q_im1 = gkyl_array_cfetch(f, gkyl_range_idx(rng, idx_m));
    const double *q_i   = gkyl_array_cfetch(f, gkyl_range_idx(rng, idx_0));
    const double *q_ip1 = gkyl_array_cfetch(f, gkyl_range_idx(rng, idx_p));

    if (q_i[26] < 1.0e-8 || q_im1[26] < 1.0e-8 || q_ip1[26] < 1.0e-8) {
      continue;
    }

    double v_im1[70] = {0.0};
    double v_i[70] = {0.0};
    double v_ip1[70] = {0.0};
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_im1, v_im1);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_i,   v_i);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_ip1, v_ip1);

    double p_im1 = ctx->K_poly * pow(v_im1[0], ctx->gas_gamma);
    double p_i   = ctx->K_poly * pow(v_i[0],   ctx->gas_gamma);
    double p_ip1 = ctx->K_poly * pow(v_ip1[0], ctx->gas_gamma);

    // Pressure-gradient proxy: this is the intuitive hydrostatic term you
    // would expect in a static star in the Newtonian limit.
    double dP_dr = (p_ip1 - p_im1) / (2.0 * dx);

    // Actual discrete finite-volume momentum-flux divergence used by the
    // wave-propagation solver.
    double div_smx = 0.0;
    compute_low_order_fv_smx_divergence(eqn, q_im1, q_i, q_ip1, dx, &div_smx);
    double minus_divF_momx = -div_smx;

    double spatial_det =
        v_i[8]  * (v_i[12] * v_i[16] - v_i[13] * v_i[15])
      - v_i[9]  * (v_i[11] * v_i[16] - v_i[13] * v_i[14])
      + v_i[10] * (v_i[11] * v_i[15] - v_i[12] * v_i[14]);
    double sqrt_det = sqrt(fmax(spatial_det, 0.0));
    double lapse = v_i[4];

    double src[70] = {0.0};
    gkyl_wv_eqn_source(eqn, q_i, src);

    // gkyl_wv_eqn_source returns lapse*src_raw_x; the explicit source update
    // later multiplies by sqrt(det(gamma))*dt before adding to conserved q.
    double src_raw_x = (fabs(lapse) > 1.0e-14) ? (src[1] / lapse) : 0.0;
    double src_rhs_x = sqrt_det * src[1];

    double x = q_i[67];
    double y = q_i[68];
    double z = q_i[69];
    double r = sqrt((x - ctx->pos_x) * (x - ctx->pos_x)
      + (y - ctx->pos_y) * (y - ctx->pos_y)
      + (z - ctx->pos_z) * (z - ctx->pos_z));

    double ratio_p = (fabs(src_raw_x) > 1.0e-30) ? dP_dr / src_raw_x : 0.0;
    double ratio_fv = (fabs(src_rhs_x) > 1.0e-30) ? minus_divF_momx / src_rhs_x : 0.0;

    printf("r=%.4f  P=%.8e  dP_dr=%.8e  -divF_momx=%.8e  src_raw_x=%.8e  src_rhs_x=%.8e  ratioP=%.6f  ratioFV=%.6f\n",
      r, p_i, dP_dr, minus_divF_momx, src_raw_x, src_rhs_x, ratio_p, ratio_fv);
  }

  printf("=== end force-balance diagnostic ===\n\n");
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
    print_tov_mass_density_diag(app, app->species[0].ctx, t_curr);
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

static void
dump_lineout(gkyl_moment_app* app, int nx)
{
  char path[128];
  snprintf(path, sizeof(path), "/tmp/rt_gr_tov_static_lineout_%d.dat", nx);

  FILE *fp = fopen(path, "w");
  if (!fp) {
    return;
  }

  struct gkyl_array *fld = app->species[0].f[0];
  struct gkyl_range *rng = &app->local;
  double x_lo = app->grid.lower[0];
  double dx = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

  fprintf(fp, "# nx=%d time=%.16e\n", nx, app->tcurr);
  fprintf(fp, "# i x Etot mom_x lapse\n");

  int idx[GKYL_MAX_DIM] = {0};
  for (int ci = rng->lower[0]; ci <= rng->upper[0]; ci++) {
    idx[0] = ci;
    const double *q = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));
    double x = x_lo + (ci - rng->lower[0] + 0.5) * dx;
    fprintf(fp, "%d %.16e %.16e %.16e %.16e\n", ci, x, q[0], q[1], q[4]);
  }

  fclose(fp);
}

static void
probe_tov_cks_metric_identities(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);

  printf("\n=== TOV/CKS 4-METRIC IDENTITY PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "gtt_BL", "gtt_3p1", "gtr_BL", "gtr_3p1", "grr_BL", "grr_3p1");

  for (int n = 0; n < ns; n++) {
    double r = sample_r[n];

    struct tov_eval_bl bl = {0.0};
    if (!gkyl_tov_eval_bl(ctx->tov, r, &bl)) {
      continue;
    }

    double alpha_BL = exp(bl.Phi);
    double f = 1.0 - 2.0 * bl.m / r;
    if (f <= 1e-12) {
      continue;
    }

    double A_BL = 1.0 / f;
    double h_prime = 2.0 * bl.m * A_BL / r;

    /* BL -> KS transformed metric components */
    double gtt_BL = -alpha_BL * alpha_BL;
    double gtr_BL = alpha_BL * alpha_BL * h_prime;
    double grr_BL = A_BL - alpha_BL * alpha_BL * h_prime * h_prime;

    /* 3+1 quantities from spacetime object */
    double alpha = 0.0;
    double *shift = gkyl_malloc(sizeof(double[3]));
    double **gam = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      gam[i] = gkyl_malloc(sizeof(double[3]));
    }

    double x = ctx->pos_x + r;
    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &alpha);
    ctx->spacetime->shift_vector_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &shift);
    ctx->spacetime->spatial_metric_tensor_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &gam);

    double beta_r_contra = shift[0];
    double gamma_rr = gam[0][0];

    double gtt_3p1 = -alpha * alpha + gamma_rr * beta_r_contra * beta_r_contra;
    double gtr_3p1 = gamma_rr * beta_r_contra;
    double grr_3p1 = gamma_rr;

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, gtt_BL, gtt_3p1, gtr_BL, gtr_3p1, grr_BL, grr_3p1);

    printf("           dgtt=% .3e  dgtr=% .3e  dgrr=% .3e\n",
      gtt_3p1 - gtt_BL, gtr_3p1 - gtr_BL, grr_3p1 - grr_BL);

    gkyl_free(shift);
    for (int i = 0; i < 3; i++) {
      gkyl_free(gam[i]);
    }
    gkyl_free(gam);
  }
}


static void
probe_tov_spacetime_balance(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);

  const double h_fd = 1e-4;

  printf("\n=== TOV/SPACETIME BALANCE PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "alpha_tov", "alpha_st", "a_ratio", "dp_dr", "rhs_alpha", "p_ratio");

  for (int n = 0; n < ns; n++) {
    double r = sample_r[n];

    struct tov_eval_bl bl = {0.0};
    gkyl_tov_eval_bl(ctx->tov, r, &bl);

    double alpha_tov = exp(bl.Phi);

    double alpha_st = 0.0;
    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0,
      ctx->pos_x + r, 0.0, 0.0, &alpha_st);

    struct tov_eval_fluid fm = {0.0}, fp = {0.0}, fc = {0.0};
    gkyl_tov_eval(ctx->tov, r - h_fd, 0.0, 0.0, &fm);
    gkyl_tov_eval(ctx->tov, r + h_fd, 0.0, 0.0, &fp);
    gkyl_tov_eval(ctx->tov, r,        0.0, 0.0, &fc);

    if (fc.rho <= 0.0 || fc.P <= 0.0) {
      continue;
    }

    double dp_dr = (fp.P - fm.P) / (2.0 * h_fd);

    double alpha_m = 0.0, alpha_p = 0.0;
    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0,
      ctx->pos_x + (r - h_fd), 0.0, 0.0, &alpha_m);
    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0,
      ctx->pos_x + (r + h_fd), 0.0, 0.0, &alpha_p);

    double dlogalpha_dr = (log(alpha_p) - log(alpha_m)) / (2.0 * h_fd);

    double h = 1.0 + fc.eps + fc.P / fc.rho;
    double rhs_alpha = -(fc.rho * h) * dlogalpha_dr;

    double a_ratio = (fabs(alpha_tov) > 0.0) ? alpha_st / alpha_tov : 0.0;
    double p_ratio = (fabs(rhs_alpha) > 0.0) ? dp_dr / rhs_alpha : 0.0;

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, alpha_tov, alpha_st, a_ratio, dp_dr, rhs_alpha, p_ratio);
  }
}

static void
probe_extrinsic_curvature_consistency(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double h_fd = 1e-4;

  printf("\n=== EXTRINSIC CURVATURE CONSISTENCY PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "Kxx", "Lbeta_gxx", "pred_Kxx", "resid", "relerr");

  for (int n = 0; n < ns; n++) {
    double r = sample_r[n];
    double x = ctx->pos_x + r;

    double alpha = 0.0;
    double *shift = gkyl_malloc(sizeof(double[3]));
    double **gam = gkyl_malloc(sizeof(double*[3]));
    double **K = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      gam[i] = gkyl_malloc(sizeof(double[3]));
      K[i] = gkyl_malloc(sizeof(double[3]));
    }

    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &alpha);
    ctx->spacetime->shift_vector_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &shift);
    ctx->spacetime->spatial_metric_tensor_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &gam);
    ctx->spacetime->extrinsic_curvature_tensor_func(ctx->spacetime, 0.0, x, 0.0, 0.0,
      h_fd, h_fd, h_fd, &K);

    double *shift_m = gkyl_malloc(sizeof(double[3]));
    double *shift_p = gkyl_malloc(sizeof(double[3]));
    double **gam_m = gkyl_malloc(sizeof(double*[3]));
    double **gam_p = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      gam_m[i] = gkyl_malloc(sizeof(double[3]));
      gam_p[i] = gkyl_malloc(sizeof(double[3]));
    }

    ctx->spacetime->shift_vector_func(ctx->spacetime, 0.0, x - h_fd, 0.0, 0.0, &shift_m);
    ctx->spacetime->shift_vector_func(ctx->spacetime, 0.0, x + h_fd, 0.0, 0.0, &shift_p);
    ctx->spacetime->spatial_metric_tensor_func(ctx->spacetime, 0.0, x - h_fd, 0.0, 0.0, &gam_m);
    ctx->spacetime->spatial_metric_tensor_func(ctx->spacetime, 0.0, x + h_fd, 0.0, 0.0, &gam_p);

    double dbetax_dx = (shift_p[0] - shift_m[0]) / (2.0 * h_fd);
    double dgamxx_dx = (gam_p[0][0] - gam_m[0][0]) / (2.0 * h_fd);

    double lie_beta_gxx = shift[0] * dgamxx_dx + 2.0 * gam[0][0] * dbetax_dx;
    double pred_Kxx = -(0.5 / alpha) * lie_beta_gxx;
    double resid = lie_beta_gxx + 2.0 * alpha * K[0][0];
    double relerr = fabs(pred_Kxx - K[0][0]) / fmax(fabs(K[0][0]), 1e-14);

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, K[0][0], lie_beta_gxx, pred_Kxx, resid, relerr);

    gkyl_free(shift);
    gkyl_free(shift_m);
    gkyl_free(shift_p);
    for (int i = 0; i < 3; i++) {
      gkyl_free(gam[i]);
      gkyl_free(K[i]);
      gkyl_free(gam_m[i]);
      gkyl_free(gam_p[i]);
    }
    gkyl_free(gam);
    gkyl_free(K);
    gkyl_free(gam_m);
    gkyl_free(gam_p);
  }
}

static void
compute_polytrope_explicit_smx_densitized(const struct tov_static_ctx *ctx, const double q[70], double *s_mx_out)
{
  double v[70] = {0.0};
  gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q, v);

  if (v[26] < pow(10.0, -8.0)) {
    *s_mx_out = 0.0;
    return;
  }

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = ctx->K_poly * pow(rho, ctx->gas_gamma);

  double lapse = v[4];
  double shift[3] = { v[5], v[6], v[7] };

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[8];  spatial_metric[0][1] = v[9];  spatial_metric[0][2] = v[10];
  spatial_metric[1][0] = v[11]; spatial_metric[1][1] = v[12]; spatial_metric[1][2] = v[13];
  spatial_metric[2][0] = v[14]; spatial_metric[2][1] = v[15]; spatial_metric[2][2] = v[16];

  double lapse_der[3] = { v[27], v[28], v[29] };

  double shift_der[3][3];
  shift_der[0][0] = v[30]; shift_der[0][1] = v[31]; shift_der[0][2] = v[32];
  shift_der[1][0] = v[33]; shift_der[1][1] = v[34]; shift_der[1][2] = v[35];
  shift_der[2][0] = v[36]; shift_der[2][1] = v[37]; shift_der[2][2] = v[38];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = v[39]; spatial_metric_der[0][0][1] = v[40]; spatial_metric_der[0][0][2] = v[41];
  spatial_metric_der[0][1][0] = v[42]; spatial_metric_der[0][1][1] = v[43]; spatial_metric_der[0][1][2] = v[44];
  spatial_metric_der[0][2][0] = v[45]; spatial_metric_der[0][2][1] = v[46]; spatial_metric_der[0][2][2] = v[47];

  spatial_metric_der[1][0][0] = v[48]; spatial_metric_der[1][0][1] = v[49]; spatial_metric_der[1][0][2] = v[50];
  spatial_metric_der[1][1][0] = v[51]; spatial_metric_der[1][1][1] = v[52]; spatial_metric_der[1][1][2] = v[53];
  spatial_metric_der[1][2][0] = v[54]; spatial_metric_der[1][2][1] = v[55]; spatial_metric_der[1][2][2] = v[56];

  spatial_metric_der[2][0][0] = v[57]; spatial_metric_der[2][0][1] = v[58]; spatial_metric_der[2][0][2] = v[59];
  spatial_metric_der[2][1][0] = v[60]; spatial_metric_der[2][1][1] = v[61]; spatial_metric_der[2][1][2] = v[62];
  spatial_metric_der[2][2][0] = v[63]; spatial_metric_der[2][2][1] = v[64]; spatial_metric_der[2][2][2] = v[65];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  double vel[3] = { vx, vy, vz };
  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / sqrt(fmax(1.0 - v_sq, 1.0e-8));
  double h = 1.0 + (ctx->gas_gamma / (ctx->gas_gamma - 1.0)) * (p / fmax(rho, 1.0e-15));

  double **stress_energy = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    stress_energy[i] = gkyl_malloc(sizeof(double[4]));
  }
  gkyl_gr_polytropic_euler_stress_energy_tensor(ctx->K_poly, ctx->gas_gamma, q, &stress_energy);

  double mom[3];
  mom[0] = (rho * h) * (W * W) * vx;
  mom[1] = (rho * h) * (W * W) * vy;
  mom[2] = (rho * h) * (W * W) * vz;

  double s_mx = -stress_energy[0][0] * lapse * lapse_der[0];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      s_mx += 0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[0][k][l];
      s_mx += 0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[0][k][l];
    }
    s_mx += (mom[k] / lapse) * shift_der[0][k];
    for (int i = 0; i < 3; i++) {
      s_mx += stress_energy[0][i + 1] * shift[k] * spatial_metric_der[0][i][k];
    }
  }

  *s_mx_out = sqrt(spatial_det) * lapse * s_mx;

  for (int i = 0; i < 4; i++) {
    gkyl_free(stress_energy[i]);
  }
  gkyl_free(stress_energy);
}

static void
compute_polytrope_explicit_smx_terms_densitized(const struct tov_static_ctx *ctx, const double q[70],
  double *s_lapse, double *s_beta_beta_dg, double *s_spatial_dg, double *s_shift_der, double *s_mixed, double *s_total)
{
  double v[70] = {0.0};
  gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q, v);

  if (v[26] < pow(10.0, -8.0)) {
    *s_lapse = *s_beta_beta_dg = *s_spatial_dg = *s_shift_der = *s_mixed = *s_total = 0.0;
    return;
  }

  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = ctx->K_poly * pow(rho, ctx->gas_gamma);

  double lapse = v[4];
  double shift[3] = { v[5], v[6], v[7] };

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[8];  spatial_metric[0][1] = v[9];  spatial_metric[0][2] = v[10];
  spatial_metric[1][0] = v[11]; spatial_metric[1][1] = v[12]; spatial_metric[1][2] = v[13];
  spatial_metric[2][0] = v[14]; spatial_metric[2][1] = v[15]; spatial_metric[2][2] = v[16];

  double lapse_der[3] = { v[27], v[28], v[29] };

  double shift_der[3][3];
  shift_der[0][0] = v[30]; shift_der[0][1] = v[31]; shift_der[0][2] = v[32];
  shift_der[1][0] = v[33]; shift_der[1][1] = v[34]; shift_der[1][2] = v[35];
  shift_der[2][0] = v[36]; shift_der[2][1] = v[37]; shift_der[2][2] = v[38];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = v[39]; spatial_metric_der[0][0][1] = v[40]; spatial_metric_der[0][0][2] = v[41];
  spatial_metric_der[0][1][0] = v[42]; spatial_metric_der[0][1][1] = v[43]; spatial_metric_der[0][1][2] = v[44];
  spatial_metric_der[0][2][0] = v[45]; spatial_metric_der[0][2][1] = v[46]; spatial_metric_der[0][2][2] = v[47];

  spatial_metric_der[1][0][0] = v[48]; spatial_metric_der[1][0][1] = v[49]; spatial_metric_der[1][0][2] = v[50];
  spatial_metric_der[1][1][0] = v[51]; spatial_metric_der[1][1][1] = v[52]; spatial_metric_der[1][1][2] = v[53];
  spatial_metric_der[1][2][0] = v[54]; spatial_metric_der[1][2][1] = v[55]; spatial_metric_der[1][2][2] = v[56];

  spatial_metric_der[2][0][0] = v[57]; spatial_metric_der[2][0][1] = v[58]; spatial_metric_der[2][0][2] = v[59];
  spatial_metric_der[2][1][0] = v[60]; spatial_metric_der[2][1][1] = v[61]; spatial_metric_der[2][1][2] = v[62];
  spatial_metric_der[2][2][0] = v[63]; spatial_metric_der[2][2][1] = v[64]; spatial_metric_der[2][2][2] = v[65];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

  double vel[3] = { vx, vy, vz };
  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / sqrt(fmax(1.0 - v_sq, 1.0e-8));
  double h = 1.0 + (ctx->gas_gamma / (ctx->gas_gamma - 1.0)) * (p / fmax(rho, 1.0e-15));

  double mom[3];
  mom[0] = (rho * h) * (W * W) * vx;
  mom[1] = (rho * h) * (W * W) * vy;
  mom[2] = (rho * h) * (W * W) * vz;

  double **stress_energy = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    stress_energy[i] = gkyl_malloc(sizeof(double[4]));
  }
  gkyl_gr_polytropic_euler_stress_energy_tensor(ctx->K_poly, ctx->gas_gamma, q, &stress_energy);

  *s_lapse = -stress_energy[0][0] * lapse * lapse_der[0];
  *s_beta_beta_dg = 0.0;
  *s_spatial_dg = 0.0;
  *s_shift_der = 0.0;
  *s_mixed = 0.0;

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      *s_beta_beta_dg += 0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[0][k][l];
      *s_spatial_dg += 0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[0][k][l];
    }
    *s_shift_der += (mom[k] / lapse) * shift_der[0][k];
    for (int i = 0; i < 3; i++) {
      *s_mixed += stress_energy[0][i + 1] * shift[k] * spatial_metric_der[0][i][k];
    }
  }

  *s_lapse *= sqrt(spatial_det) * lapse;
  *s_beta_beta_dg *= sqrt(spatial_det) * lapse;
  *s_spatial_dg *= sqrt(spatial_det) * lapse;
  *s_shift_der *= sqrt(spatial_det) * lapse;
  *s_mixed *= sqrt(spatial_det) * lapse;
  *s_total = *s_lapse + *s_beta_beta_dg + *s_spatial_dg + *s_shift_der + *s_mixed;

  for (int i = 0; i < 4; i++) {
    gkyl_free(stress_energy[i]);
  }
  gkyl_free(stress_energy);
}

static void
probe_stress_energy_consistency(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);

  printf("\n=== STRESS-ENERGY CONSISTENCY PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "E_exp", "E_wv", "E_loc", "Sx_exp", "Sx_wv", "Sx_loc");

  for (int n = 0; n < ns; n++) {
    double r = sample_r[n];
    double xn[3] = { ctx->pos_x + r, 0.0, 0.0 };
    double q[70] = {0.0};
    double v[70] = {0.0};

    evalGRPolytropicEulerInit(0.0, xn, q, (void*) ctx);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q, v);

    if (v[26] < pow(10.0, -8.0)) {
      continue;
    }

    double rho = v[0];
    double vx = v[1];
    double vy = v[2];
    double vz = v[3];
    double p = ctx->K_poly * pow(rho, ctx->gas_gamma);

    double lapse = v[4];
    double shift[3] = { v[5], v[6], v[7] };

    double spatial_metric[3][3];
    spatial_metric[0][0] = v[8];  spatial_metric[0][1] = v[9];  spatial_metric[0][2] = v[10];
    spatial_metric[1][0] = v[11]; spatial_metric[1][1] = v[12]; spatial_metric[1][2] = v[13];
    spatial_metric[2][0] = v[14]; spatial_metric[2][1] = v[15]; spatial_metric[2][2] = v[16];

    double vel[3] = { vx, vy, vz };
    double v_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }
    double W = 1.0 / sqrt(fmax(1.0 - v_sq, 1.0e-8));
    double h = 1.0 + (ctx->gas_gamma / (ctx->gas_gamma - 1.0)) * (p / fmax(rho, 1.0e-15));

    double cov_vel[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }

    double E_exp = (rho * h * W * W) - p;
    double Sx_exp = (rho * h * W * W) * cov_vel[0];

    double **T_wv = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      T_wv[i] = gkyl_malloc(sizeof(double[4]));
    }
    gkyl_gr_polytropic_euler_stress_energy_tensor(ctx->K_poly, ctx->gas_gamma, q, &T_wv);

    double beta_cov_x = 0.0;
    for (int i = 0; i < 3; i++) {
      beta_cov_x += spatial_metric[0][i] * shift[i];
    }
    double E_wv = lapse * lapse * T_wv[0][0];
    double Sx_wv = lapse * (beta_cov_x * T_wv[0][0]
      + spatial_metric[0][0] * T_wv[0][1]
      + spatial_metric[0][1] * T_wv[0][2]
      + spatial_metric[0][2] * T_wv[0][3]);

    double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    }
    gkyl_gr_polytropic_euler_inv_spatial_metric(q, &inv_spatial_metric);

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift[0] * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift[1] * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift[2] * (W / lapse));

    double inv_spacetime_metric[4][4] = {0.0};
    inv_spacetime_metric[0][0] = -(1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i + 1] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i + 1][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i + 1][j + 1] = inv_spatial_metric[i][j]
          - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double T_loc[4][4] = {0.0};
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        T_loc[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
      }
    }
    double E_loc = lapse * lapse * T_loc[0][0];
    double Sx_loc = lapse * (beta_cov_x * T_loc[0][0]
      + spatial_metric[0][0] * T_loc[0][1]
      + spatial_metric[0][1] * T_loc[0][2]
      + spatial_metric[0][2] * T_loc[0][3]);

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, E_exp, E_wv, E_loc, Sx_exp, Sx_wv, Sx_loc);
    printf("           dE_wv=% .3e dE_loc=% .3e dS_wv=% .3e dS_loc=% .3e\n",
      E_wv - E_exp, E_loc - E_exp, Sx_wv - Sx_exp, Sx_loc - Sx_exp);

    for (int i = 0; i < 4; i++) {
      gkyl_free(T_wv[i]);
    }
    gkyl_free(T_wv);
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_spatial_metric[i]);
    }
    gkyl_free(inv_spatial_metric);
  }
}

static void
probe_flux_source_balance(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double dx = ctx->Lx / ctx->Nx;

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  printf("\n=== FV / EQN-SOURCE / EXPL-SOURCE PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "fvDiv_Sx", "Sx_eqn", "Sx_expl", "FV/eqn", "FV/expl");

  for (int n = 0; n < ns; n++) {
    double x_center = nearest_cell_center_x(ctx, ctx->pos_x + sample_r[n]);
    double r = fabs(x_center - ctx->pos_x);

    double xn_im1[3] = { x_center - dx, 0.0, 0.0 };
    double xn_i[3]   = { x_center,      0.0, 0.0 };
    double xn_ip1[3] = { x_center + dx, 0.0, 0.0 };

    double q_im1[70] = {0.0}, q_i[70] = {0.0}, q_ip1[70] = {0.0};
    double s_0[70] = {0.0};

    evalGRPolytropicEulerInit(0.0, xn_im1, q_im1, (void*) ctx);
    evalGRPolytropicEulerInit(0.0, xn_i, q_i, (void*) ctx);
    evalGRPolytropicEulerInit(0.0, xn_ip1, q_ip1, (void*) ctx);

    gkyl_wv_eqn_source(eqn, q_i, s_0);

    double fv_div_sx = 0.0;
    compute_low_order_fv_smx_divergence(eqn, q_im1, q_i, q_ip1, dx, &fv_div_sx);

    double spatial_metric[3][3];
    spatial_metric[0][0] = q_i[8];  spatial_metric[0][1] = q_i[9];  spatial_metric[0][2] = q_i[10];
    spatial_metric[1][0] = q_i[11]; spatial_metric[1][1] = q_i[12]; spatial_metric[1][2] = q_i[13];
    spatial_metric[2][0] = q_i[14]; spatial_metric[2][1] = q_i[15]; spatial_metric[2][2] = q_i[16];
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
    double Sx_eqn = sqrt(spatial_det) * s_0[1];
    double Sx_expl = 0.0;
    compute_polytrope_explicit_smx_densitized(ctx, q_i, &Sx_expl);
    double ratio_eqn = (fabs(Sx_eqn) > 0.0) ? fv_div_sx / Sx_eqn : 0.0;
    double ratio_expl = (fabs(Sx_expl) > 0.0) ? fv_div_sx / Sx_expl : 0.0;

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, fv_div_sx, Sx_eqn, Sx_expl, ratio_eqn, ratio_expl);
  }

  gkyl_wv_eqn_release(eqn);
}

static void
probe_shifted_static_balance(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.5,
    0.25 * ctx->R_star,
    0.5 * ctx->R_star,
    0.75 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double dx = ctx->Lx / ctx->Nx;

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  printf("\n=== SHIFTED STATIC BALANCE PROBE ===\n");
  printf("%-10s %-12s %-12s %-12s %-12s %-12s\n",
    "r", "mom_x", "fvDiv_Sx", "Sx_eqn", "Sx_expl", "FV/expl");

  for (int n = 0; n < ns; n++) {
    double x_center = nearest_cell_center_x(ctx, ctx->pos_x + sample_r[n]);
    double r = fabs(x_center - ctx->pos_x);
    double xn_im1[3] = { x_center - dx, 0.0, 0.0 };
    double xn_i[3]   = { x_center,      0.0, 0.0 };
    double xn_ip1[3] = { x_center + dx, 0.0, 0.0 };

    double q_im1[70] = {0.0}, q_i[70] = {0.0}, q_ip1[70] = {0.0}, s_0[70] = {0.0};

    build_shifted_static_state(0.0, xn_im1, q_im1, (void*) ctx);
    build_shifted_static_state(0.0, xn_i, q_i, (void*) ctx);
    build_shifted_static_state(0.0, xn_ip1, q_ip1, (void*) ctx);

    gkyl_wv_eqn_source(eqn, q_i, s_0);

    double spatial_metric[3][3];
    spatial_metric[0][0] = q_i[8];  spatial_metric[0][1] = q_i[9];  spatial_metric[0][2] = q_i[10];
    spatial_metric[1][0] = q_i[11]; spatial_metric[1][1] = q_i[12]; spatial_metric[1][2] = q_i[13];
    spatial_metric[2][0] = q_i[14]; spatial_metric[2][1] = q_i[15]; spatial_metric[2][2] = q_i[16];
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double fv_div_sx = 0.0;
    compute_low_order_fv_smx_divergence(eqn, q_im1, q_i, q_ip1, dx, &fv_div_sx);

    double Sx_eqn = sqrt(spatial_det) * s_0[1];
    double Sx_expl = 0.0;
    compute_polytrope_explicit_smx_densitized(ctx, q_i, &Sx_expl);
    double ratio_expl = (fabs(Sx_expl) > 0.0) ? fv_div_sx / Sx_expl : 0.0;

    printf("%-10.6f %-12.6e %-12.6e %-12.6e %-12.6e %-12.6e\n",
      r, q_i[1], fv_div_sx, Sx_eqn, Sx_expl, ratio_expl);
  }

  gkyl_wv_eqn_release(eqn);
}

static void
probe_low_order_face_flux_detail(const struct tov_static_ctx *ctx, bool shifted_static)
{
  double sample_r[] = {
    0.5 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double dx = ctx->Lx / ctx->Nx;

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  printf("\n=== %s FACE FLUX DETAIL ===\n",
    shifted_static ? "SHIFTED STATIC" : "INITIAL");

  for (int n = 0; n < ns; n++) {
    double x_left = nearest_cell_center_x(ctx, ctx->pos_x + sample_r[n]);
    double x_right = x_left + dx;
    double r_face = fabs((0.5 * (x_left + x_right)) - ctx->pos_x);

    double xn_l[3] = { x_left, 0.0, 0.0 };
    double xn_r[3] = { x_right, 0.0, 0.0 };
    double ql[70] = {0.0}, qr[70] = {0.0};
    double fl[70] = {0.0}, fr[70] = {0.0};
    double flux_jump[70] = {0.0};
    double amdq[70] = {0.0}, apdq[70] = {0.0};

    if (shifted_static) {
      build_shifted_static_state(0.0, xn_l, ql, (void*) ctx);
      build_shifted_static_state(0.0, xn_r, qr, (void*) ctx);
    }
    else {
      evalGRPolytropicEulerInit(0.0, xn_l, ql, (void*) ctx);
      evalGRPolytropicEulerInit(0.0, xn_r, qr, (void*) ctx);
    }

    gkyl_gr_polytropic_euler_flux(ctx->K_poly, ctx->gas_gamma, ql, fl);
    gkyl_gr_polytropic_euler_flux(ctx->K_poly, ctx->gas_gamma, qr, fr);
    compute_low_order_interface_fluct(eqn, ql, qr, amdq, apdq);
    gkyl_wv_eqn_flux_jump(eqn, ql, qr, flux_jump);

    double vl[70] = {0.0}, vr[70] = {0.0};
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, ql, vl);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, qr, vr);

    double pl = ctx->K_poly * pow(vl[0], ctx->gas_gamma);
    double pr = ctx->K_poly * pow(vr[0], ctx->gas_gamma);

    double gaml[3][3], gamr[3][3];
    gaml[0][0] = vl[8];  gaml[0][1] = vl[9];  gaml[0][2] = vl[10];
    gaml[1][0] = vl[11]; gaml[1][1] = vl[12]; gaml[1][2] = vl[13];
    gaml[2][0] = vl[14]; gaml[2][1] = vl[15]; gaml[2][2] = vl[16];
    gamr[0][0] = vr[8];  gamr[0][1] = vr[9];  gamr[0][2] = vr[10];
    gamr[1][0] = vr[11]; gamr[1][1] = vr[12]; gamr[1][2] = vr[13];
    gamr[2][0] = vr[14]; gamr[2][1] = vr[15]; gamr[2][2] = vr[16];

    double detl = (gaml[0][0] * ((gaml[1][1] * gaml[2][2]) - (gaml[2][1] * gaml[1][2]))) -
      (gaml[0][1] * ((gaml[1][0] * gaml[2][2]) - (gaml[1][2] * gaml[2][0]))) +
      (gaml[0][2] * ((gaml[1][0] * gaml[2][1]) - (gaml[1][1] * gaml[2][0])));
    double detr = (gamr[0][0] * ((gamr[1][1] * gamr[2][2]) - (gamr[2][1] * gamr[1][2]))) -
      (gamr[0][1] * ((gamr[1][0] * gamr[2][2]) - (gamr[1][2] * gamr[2][0]))) +
      (gamr[0][2] * ((gamr[1][0] * gamr[2][1]) - (gamr[1][1] * gamr[2][0])));

    double press_l = vl[4] * sqrt(detl) * pl;
    double press_r = vr[4] * sqrt(detr) * pr;
    double press_jump = press_r - press_l;
    double full_jump = flux_jump[1];
    double fluct_sum = amdq[1] + apdq[1];
    double diffusive_residual = full_jump - press_jump;

    printf("face_r=%.6f xL=%.6f xR=%.6f\n", r_face, x_left, x_right);
    printf("  flux_jump=%.6e  amdq+apdq=%.6e  diff=%.6e\n",
      full_jump, fluct_sum, fluct_sum - full_jump);
    printf("  pressure_jump=%.6e  nonpressure_jump=%.6e\n",
      press_jump, diffusive_residual);
    printf("  qL_mom=%.6e qR_mom=%.6e  pL=%.6e pR=%.6e\n",
      ql[1], qr[1], pl, pr);
  }

  gkyl_wv_eqn_release(eqn);
}

static void
probe_well_balanced_source_detail(const struct tov_static_ctx *ctx, bool shifted_static)
{
  double sample_r[] = {
    0.5 * ctx->R_star,
    0.9 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double dx = ctx->Lx / ctx->Nx;

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  printf("\n=== %s WELL-BALANCED SOURCE DETAIL ===\n",
    shifted_static ? "SHIFTED STATIC" : "INITIAL");

  for (int n = 0; n < ns; n++) {
    double x_center = nearest_cell_center_x(ctx, ctx->pos_x + sample_r[n]);
    double r = fabs(x_center - ctx->pos_x);
    double xn_im1[3] = { x_center - dx, 0.0, 0.0 };
    double xn_i[3]   = { x_center,      0.0, 0.0 };
    double xn_ip1[3] = { x_center + dx, 0.0, 0.0 };

    double q_im1[70] = {0.0}, q_i[70] = {0.0}, q_ip1[70] = {0.0};
    if (shifted_static) {
      build_shifted_static_state(0.0, xn_im1, q_im1, (void*) ctx);
      build_shifted_static_state(0.0, xn_i, q_i, (void*) ctx);
      build_shifted_static_state(0.0, xn_ip1, q_ip1, (void*) ctx);
    }
    else {
      evalGRPolytropicEulerInit(0.0, xn_im1, q_im1, (void*) ctx);
      evalGRPolytropicEulerInit(0.0, xn_i, q_i, (void*) ctx);
      evalGRPolytropicEulerInit(0.0, xn_ip1, q_ip1, (void*) ctx);
    }

    double fv_div_sx = 0.0;
    compute_low_order_fv_smx_divergence(eqn, q_im1, q_i, q_ip1, dx, &fv_div_sx);

    double s_lapse = 0.0, s_beta_beta_dg = 0.0, s_spatial_dg = 0.0;
    double s_shift_der = 0.0, s_mixed = 0.0, s_total = 0.0;
    compute_polytrope_explicit_smx_terms_densitized(ctx, q_i,
      &s_lapse, &s_beta_beta_dg, &s_spatial_dg, &s_shift_der, &s_mixed, &s_total);

    double s_rest = s_total - s_lapse;
    double s_wb_total = fv_div_sx + s_rest;
    double ratio_old = fabs(s_total) > 0.0 ? fv_div_sx / s_total : 0.0;
    double ratio_wb = fabs(s_wb_total) > 0.0 ? fv_div_sx / s_wb_total : 0.0;

    printf("r=%.6f x=%.6f\n", r, x_center);
    printf("  fv_div=%.6e old_total=%.6e lapse=%.6e rest=%.6e\n",
      fv_div_sx, s_total, s_lapse, s_rest);
    printf("  wb_total=%.6e wb_static=%.6e ratios old=%.6e wb=%.6e\n",
      s_wb_total, fv_div_sx, ratio_old, ratio_wb);
  }

  gkyl_wv_eqn_release(eqn);
}

static void
probe_shifted_static_surface_detail(const struct tov_static_ctx *ctx)
{
  double sample_r[] = {
    0.80 * ctx->R_star,
    0.85 * ctx->R_star,
    0.9 * ctx->R_star,
    0.94 * ctx->R_star,
    0.96 * ctx->R_star,
    0.98 * ctx->R_star
  };
  int ns = sizeof(sample_r) / sizeof(sample_r[0]);
  const double h_fd = 1e-4;

  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  printf("\n=== SHIFTED STATIC SURFACE DETAIL ===\n");

  for (int n = 0; n < ns; n++) {
    double r = sample_r[n];
    double xn_m[3] = { ctx->pos_x + (r - h_fd), 0.0, 0.0 };
    double xn_0[3] = { ctx->pos_x + r,          0.0, 0.0 };
    double xn_p[3] = { ctx->pos_x + (r + h_fd), 0.0, 0.0 };

    double q_m[70] = {0.0}, q_0[70] = {0.0}, q_p[70] = {0.0};
    double f_m[70] = {0.0}, f_0[70] = {0.0}, f_p[70] = {0.0};
    double s_0[70] = {0.0};
    double v_0[70] = {0.0};

    build_shifted_static_state(0.0, xn_m, q_m, (void*) ctx);
    build_shifted_static_state(0.0, xn_0, q_0, (void*) ctx);
    build_shifted_static_state(0.0, xn_p, q_p, (void*) ctx);

    gkyl_gr_polytropic_euler_flux(ctx->K_poly, ctx->gas_gamma, q_m, f_m);
    gkyl_gr_polytropic_euler_flux(ctx->K_poly, ctx->gas_gamma, q_0, f_0);
    gkyl_gr_polytropic_euler_flux(ctx->K_poly, ctx->gas_gamma, q_p, f_p);
    gkyl_wv_eqn_source(eqn, q_0, s_0);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_0, v_0);

    double rho = v_0[0];
    double p = ctx->K_poly * pow(rho, ctx->gas_gamma);
    double lapse = v_0[4];
    double shift_x = v_0[5];
    double spatial_metric[3][3];
    spatial_metric[0][0] = v_0[8]; spatial_metric[0][1] = v_0[9]; spatial_metric[0][2] = v_0[10];
    spatial_metric[1][0] = v_0[11]; spatial_metric[1][1] = v_0[12]; spatial_metric[1][2] = v_0[13];
    spatial_metric[2][0] = v_0[14]; spatial_metric[2][1] = v_0[15]; spatial_metric[2][2] = v_0[16];
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double pressure_flux_m = f_m[4];
    double pressure_flux_0 = f_0[4];
    double pressure_flux_p = f_p[4];
    (void) pressure_flux_0;
    double pressure_part_m = q_m[0]*0.0;
    double pressure_part_p = q_p[0]*0.0;
    double v_m[70] = {0.0}, v_p[70] = {0.0};
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_m, v_m);
    gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, q_p, v_p);
    double p_m = ctx->K_poly * pow(v_m[0], ctx->gas_gamma);
    double p_p = ctx->K_poly * pow(v_p[0], ctx->gas_gamma);
    double spatial_metric_m[3][3], spatial_metric_p[3][3];
    spatial_metric_m[0][0] = v_m[8]; spatial_metric_m[0][1] = v_m[9]; spatial_metric_m[0][2] = v_m[10];
    spatial_metric_m[1][0] = v_m[11]; spatial_metric_m[1][1] = v_m[12]; spatial_metric_m[1][2] = v_m[13];
    spatial_metric_m[2][0] = v_m[14]; spatial_metric_m[2][1] = v_m[15]; spatial_metric_m[2][2] = v_m[16];
    spatial_metric_p[0][0] = v_p[8]; spatial_metric_p[0][1] = v_p[9]; spatial_metric_p[0][2] = v_p[10];
    spatial_metric_p[1][0] = v_p[11]; spatial_metric_p[1][1] = v_p[12]; spatial_metric_p[1][2] = v_p[13];
    spatial_metric_p[2][0] = v_p[14]; spatial_metric_p[2][1] = v_p[15]; spatial_metric_p[2][2] = v_p[16];
    double spatial_det_m = (spatial_metric_m[0][0] * ((spatial_metric_m[1][1] * spatial_metric_m[2][2]) - (spatial_metric_m[2][1] * spatial_metric_m[1][2]))) -
      (spatial_metric_m[0][1] * ((spatial_metric_m[1][0] * spatial_metric_m[2][2]) - (spatial_metric_m[1][2] * spatial_metric_m[2][0]))) +
      (spatial_metric_m[0][2] * ((spatial_metric_m[1][0] * spatial_metric_m[2][1]) - (spatial_metric_m[1][1] * spatial_metric_m[2][0])));
    double spatial_det_p = (spatial_metric_p[0][0] * ((spatial_metric_p[1][1] * spatial_metric_p[2][2]) - (spatial_metric_p[2][1] * spatial_metric_p[1][2]))) -
      (spatial_metric_p[0][1] * ((spatial_metric_p[1][0] * spatial_metric_p[2][2]) - (spatial_metric_p[1][2] * spatial_metric_p[2][0]))) +
      (spatial_metric_p[0][2] * ((spatial_metric_p[1][0] * spatial_metric_p[2][1]) - (spatial_metric_p[1][1] * spatial_metric_p[2][0])));
    pressure_part_m = v_m[4] * sqrt(spatial_det_m) * p_m;
    pressure_part_p = v_p[4] * sqrt(spatial_det_p) * p_p;

    double dFdx_total = (f_p[1] - f_m[1]) / (2.0 * h_fd);
    double dFdx_press = (pressure_part_p - pressure_part_m) / (2.0 * h_fd);
    double dFdx_adv = dFdx_total - dFdx_press;

    double s_lapse, s_beta_beta_dg, s_spatial_dg, s_shift_der, s_mixed, s_total;
    compute_polytrope_explicit_smx_terms_densitized(ctx, q_0,
      &s_lapse, &s_beta_beta_dg, &s_spatial_dg, &s_shift_der, &s_mixed, &s_total);

    double ratio_total = fabs(s_total) > 0.0 ? dFdx_total / s_total : 0.0;
    double ratio_press = fabs(s_lapse) > 0.0 ? dFdx_press / s_lapse : 0.0;
    double ratio_adv = fabs(s_total - s_lapse) > 0.0 ? dFdx_adv / (s_total - s_lapse) : 0.0;

    printf("r=%.6f rho=%.6e p=%.6e lapse=%.6e shift_x=%.6e mom_x=%.6e\n",
      r, rho, p, lapse, shift_x, q_0[1]);
    printf("  flux: total=%.6e pressure=%.6e advective=%.6e\n",
      dFdx_total, dFdx_press, dFdx_adv);
    printf("  source: total=%.6e lapse=%.6e beta_beta_dg=%.6e spatial_dg=%.6e shift_der=%.6e mixed=%.6e\n",
      s_total, s_lapse, s_beta_beta_dg, s_spatial_dg, s_shift_der, s_mixed);
    printf("  ratios: total/source=%.6e pressure/lapse=%.6e adv/(source-lapse)=%.6e\n",
      ratio_total, ratio_press, ratio_adv);
  }

  gkyl_wv_eqn_release(eqn);
}


static void
report_first_bad_cell(gkyl_moment_app *app, const struct tov_static_ctx *ctx, double t_curr)
{
  static bool reported = false;
  if (reported) {
    return;
  }

  struct gkyl_array *fld = app->species[0].f[0];
  struct gkyl_range *rng = &app->local;
  double x_lo = app->grid.lower[0];
  double dx_cell = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

  int idx[GKYL_MAX_DIM] = {0};

  for (int ci = rng->lower[0]; ci <= rng->upper[0]; ci++) {
    idx[0] = ci;
    const double *q = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));

    double x = x_lo + (ci - rng->lower[0] + 0.5) * dx_cell;
    double r = fabs(x - ctx->pos_x);

    bool bad = false;
    for (int k = 0; k < 5; k++) {
      if (isnan(q[k]) || isinf(q[k])) {
        bad = true;
      }
    }
    if (fabs(q[0]) > 1e6 || fabs(q[1]) > 1e6) {
      bad = true;
    }

    if (bad) {
      printf("[BAD-CELL] t=%.16e i=%d x=%.16e r=%.16e Etot=%.16e mom_x=%.16e lapse=%.16e\n",
        t_curr, ci, x, r, q[0], q[1], q[4]);

      printf("[BAD-CELL-NEIGHBORHOOD]\n");
      printf("%-6s %-14s %-14s %-14s %-14s %-14s\n",
        "i", "x", "Etot", "mom_x", "rho_rec", "p_rec");
      for (int cj = ci - 2; cj <= ci + 2; cj++) {
        if (cj < rng->lower[0] || cj > rng->upper[0]) {
          continue;
        }

        idx[0] = cj;
        const double *qj = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));
        double xj = x_lo + (cj - rng->lower[0] + 0.5) * dx_cell;

        double v[70] = {0.0};
        gkyl_gr_polytropic_euler_prim_vars(ctx->K_poly, ctx->gas_gamma, qj, v);

        double rho_rec = v[0];
        double p_rec = ctx->K_poly * pow(fmax(rho_rec, 0.0), ctx->gas_gamma);

        printf("%-6d %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n",
          cj, xj, qj[0], qj[1], rho_rec, p_rec);
      }
      reported = true;
      return;
    }
  }
}


struct conv_norms {
  double l1;
  double l2;
  double linf;
};

static void
append_energy_convergence_record(const char *kind, int nx, int nx_fine, double time,
  struct conv_norms etot)
{
  FILE *fp = fopen("gr_tov_static_energy_conv.dat", "a");
  if (!fp) {
    return;
  }

  fprintf(fp, "%s %d %d %.16e %.16e %.16e %.16e\n",
    kind, nx, nx_fine, time, etot.l1, etot.l2, etot.linf);
  fclose(fp);
}

static void
compute_exact_error_norms(gkyl_moment_app* app, const struct tov_static_ctx* ctx)
{
  struct gkyl_array *fld = app->species[0].f[0];
  struct gkyl_range *rng = &app->local;
  double x_lo = app->grid.lower[0];
  double dx_cell = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

  double l1_etot = 0.0, l2_etot = 0.0, linf_etot = 0.0;
  double l1_momx = 0.0, l2_momx = 0.0, linf_momx = 0.0;
  double l1_lapse = 0.0, l2_lapse = 0.0, linf_lapse = 0.0;
  int ncell = 0;

  int idx[GKYL_MAX_DIM] = {0};

  for (int ci = rng->lower[0]; ci <= rng->upper[0]; ci++) {
    idx[0] = ci;

    const double *q = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));
    double x = x_lo + (ci - rng->lower[0] + 0.5) * dx_cell;
    double dx = x - ctx->pos_x;
    double r = fabs(dx);

    struct tov_eval_fluid fluid = {0.0};

    double rho_exact = 0.0;
    double p_exact = 0.0;
    double eps_exact = 0.0;

    if (r <= ctx->R_star) {
      gkyl_tov_eval(ctx->tov, dx, 0.0, 0.0, &fluid);

      rho_exact = fluid.rho;
      p_exact = fluid.P;
      eps_exact = fluid.eps;

      if (rho_exact < ctx->rho_atm) {
        rho_exact = ctx->rho_atm;
      }
      if (p_exact < ctx->p_atm) {
        p_exact = ctx->p_atm;
      }

      eps_exact = p_exact / ((ctx->gas_gamma - 1.0) * rho_exact);
    }
    else {
      rho_exact = ctx->rho_atm;
      p_exact = ctx->p_atm;
      eps_exact = p_exact / ((ctx->gas_gamma - 1.0) * rho_exact);
    }

    double lapse_exact = 0.0;
    double spatial_det = 0.0;

    ctx->spacetime->lapse_function_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &lapse_exact);
    ctx->spacetime->spatial_metric_det_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &spatial_det);

    // double *shift_exact = gkyl_malloc(sizeof(double[3]));
    // double **spatial_metric_exact = gkyl_malloc(sizeof(double*[3]));
    // for (int i = 0; i < 3; i++) {
    //   spatial_metric_exact[i] = gkyl_malloc(sizeof(double[3]));
    // }

    // ctx->spacetime->shift_vector_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &shift_exact);
    // ctx->spacetime->spatial_metric_tensor_func(ctx->spacetime, 0.0, x, 0.0, 0.0, &spatial_metric_exact);

    double h_exact = 1.0 + eps_exact + p_exact / rho_exact;
    double Etot_exact = sqrt(spatial_det) * ((rho_exact * h_exact) - p_exact);
    double mom_x_exact = 0.0;
    // double vel_exact[3] = {0.0, 0.0, 0.0};
    // if (r <= ctx->R_star) {
    //   vel_exact[0] = shift_exact[0] / lapse_exact;
    //   vel_exact[1] = shift_exact[1] / lapse_exact;
    //   vel_exact[2] = shift_exact[2] / lapse_exact;
    // }
    // double v_sq_exact = 0.0;
    // for (int i = 0; i < 3; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     v_sq_exact += spatial_metric_exact[i][j] * vel_exact[i] * vel_exact[j];
    //   }
    // }

    // double W_exact = 1.0 / sqrt(1.0 - v_sq_exact);
    // double h_exact = 1.0 + eps_exact + p_exact / rho_exact;

    // double Etot_exact = sqrt(spatial_det) * ((rho_exact * h_exact * W_exact * W_exact) - p_exact);
    // double mom_x_exact = sqrt(spatial_det) * (rho_exact * h_exact * W_exact * W_exact) * vel_exact[0];
    // double cov_vel_exact[3] = {0.0, 0.0, 0.0};
    // for (int i = 0; i < 3; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     cov_vel_exact[i] += spatial_metric_exact[i][j] * vel_exact[j];
    //   }
    // }

    // double mom_x_exact = sqrt(spatial_det) * (rho_exact * h_exact * W_exact * W_exact) * cov_vel_exact[0];


    double err_etot = fabs(q[0] - Etot_exact);
    double err_momx = fabs(q[1] - mom_x_exact);
    double err_lapse = fabs(q[4] - lapse_exact);

    l1_etot += err_etot;
    l2_etot += err_etot * err_etot;
    if (err_etot > linf_etot) linf_etot = err_etot;

    l1_momx += err_momx;
    l2_momx += err_momx * err_momx;
    if (err_momx > linf_momx) linf_momx = err_momx;

    l1_lapse += err_lapse;
    l2_lapse += err_lapse * err_lapse;
    if (err_lapse > linf_lapse) linf_lapse = err_lapse;

    ncell++;
  }

  if (ncell > 0) {
    l1_etot /= ncell;
    l2_etot = sqrt(l2_etot / ncell);

    l1_momx /= ncell;
    l2_momx = sqrt(l2_momx / ncell);

    l1_lapse /= ncell;
    l2_lapse = sqrt(l2_lapse / ncell);
  }

  printf("[EXACT-CONV] nx=%d time=%.16e\n", app->grid.cells[0], app->tcurr);
  printf("[EXACT-CONV] Etot  L1=%.16e L2=%.16e Linf=%.16e\n", l1_etot, l2_etot, linf_etot);
  printf("[EXACT-CONV] mom_x L1=%.16e L2=%.16e Linf=%.16e\n", l1_momx, l2_momx, linf_momx);
  printf("[EXACT-CONV] lapse L1=%.16e L2=%.16e Linf=%.16e\n", l1_lapse, l2_lapse, linf_lapse);

  append_energy_convergence_record("exact", app->grid.cells[0], 0, app->tcurr,
    (struct conv_norms) { .l1 = l1_etot, .l2 = l2_etot, .linf = linf_etot });
}


static int
compute_self_convergence(const char *coarse_path, const char *fine_path)
{
  FILE *fc = fopen(coarse_path, "r");
  FILE *ff = fopen(fine_path, "r");
  if (!fc || !ff) {
    if (fc) fclose(fc);
    if (ff) fclose(ff);
    return 0;
  }

  int coarse_n = 0, fine_n = 0;
  char line[512];

  while (fgets(line, sizeof(line), fc)) {
    if (line[0] != '#') coarse_n++;
  }
  while (fgets(line, sizeof(line), ff)) {
    if (line[0] != '#') fine_n++;
  }

  if (coarse_n <= 0 || fine_n != 2 * coarse_n) {
    fclose(fc);
    fclose(ff);
    return 0;
  }

  rewind(fc);
  rewind(ff);

  double *coarse_x = gkyl_malloc(sizeof(double[coarse_n]));
  double *coarse_etot = gkyl_malloc(sizeof(double[coarse_n]));
  double *coarse_momx = gkyl_malloc(sizeof(double[coarse_n]));
  double *coarse_lapse = gkyl_malloc(sizeof(double[coarse_n]));

  double *fine_x = gkyl_malloc(sizeof(double[fine_n]));
  double *fine_etot = gkyl_malloc(sizeof(double[fine_n]));
  double *fine_momx = gkyl_malloc(sizeof(double[fine_n]));
  double *fine_lapse = gkyl_malloc(sizeof(double[fine_n]));

  int i = 0;
  while (fgets(line, sizeof(line), fc)) {
    if (line[0] == '#') continue;
    int cell_idx;
    if (sscanf(line, "%d %lf %lf %lf %lf",
        &cell_idx, &coarse_x[i], &coarse_etot[i], &coarse_momx[i], &coarse_lapse[i]) == 5) {
      i++;
    }
  }

  i = 0;
  while (fgets(line, sizeof(line), ff)) {
    if (line[0] == '#') continue;
    int cell_idx;
    if (sscanf(line, "%d %lf %lf %lf %lf",
        &cell_idx, &fine_x[i], &fine_etot[i], &fine_momx[i], &fine_lapse[i]) == 5) {
      i++;
    }
  }

  fclose(fc);
  fclose(ff);

  double l1_etot = 0.0, l2_etot = 0.0, linf_etot = 0.0;
  double l1_momx = 0.0, l2_momx = 0.0, linf_momx = 0.0;
  double l1_lapse = 0.0, l2_lapse = 0.0, linf_lapse = 0.0;

  for (int c = 0; c < coarse_n; ++c) {
    int f0 = 2 * c;
    int f1 = 2 * c + 1;

    double fine_etot_restricted = 0.5 * (fine_etot[f0] + fine_etot[f1]);
    double fine_momx_restricted = 0.5 * (fine_momx[f0] + fine_momx[f1]);
    double fine_lapse_restricted = 0.5 * (fine_lapse[f0] + fine_lapse[f1]);

    double err_etot = fabs(coarse_etot[c] - fine_etot_restricted);
    double err_momx = fabs(coarse_momx[c] - fine_momx_restricted);
    double err_lapse = fabs(coarse_lapse[c] - fine_lapse_restricted);

    l1_etot += err_etot;
    l2_etot += err_etot * err_etot;
    if (err_etot > linf_etot) linf_etot = err_etot;

    l1_momx += err_momx;
    l2_momx += err_momx * err_momx;
    if (err_momx > linf_momx) linf_momx = err_momx;

    l1_lapse += err_lapse;
    l2_lapse += err_lapse * err_lapse;
    if (err_lapse > linf_lapse) linf_lapse = err_lapse;
  }

  l1_etot /= coarse_n;
  l2_etot = sqrt(l2_etot / coarse_n);

  l1_momx /= coarse_n;
  l2_momx = sqrt(l2_momx / coarse_n);

  l1_lapse /= coarse_n;
  l2_lapse = sqrt(l2_lapse / coarse_n);

  printf("[SELF-CONV] coarse=%s fine=%s\n", coarse_path, fine_path);
  printf("[SELF-CONV] Etot  L1=%.16e L2=%.16e Linf=%.16e\n", l1_etot, l2_etot, linf_etot);
  printf("[SELF-CONV] mom_x L1=%.16e L2=%.16e Linf=%.16e\n", l1_momx, l2_momx, linf_momx);
  printf("[SELF-CONV] lapse L1=%.16e L2=%.16e Linf=%.16e\n", l1_lapse, l2_lapse, linf_lapse);

  int coarse_nx = 0, fine_nx = 0;
  if (sscanf(coarse_path, "/tmp/rt_gr_tov_static_lineout_%d.dat", &coarse_nx) == 1 &&
      sscanf(fine_path, "/tmp/rt_gr_tov_static_lineout_%d.dat", &fine_nx) == 1) {
    append_energy_convergence_record("self", coarse_nx, fine_nx, 0.0,
      (struct conv_norms) { .l1 = l1_etot, .l2 = l2_etot, .linf = linf_etot });
  }

  gkyl_free(coarse_x);
  gkyl_free(coarse_etot);
  gkyl_free(coarse_momx);
  gkyl_free(coarse_lapse);
  gkyl_free(fine_x);
  gkyl_free(fine_etot);
  gkyl_free(fine_momx);
  gkyl_free(fine_lapse);

  return 1;
}

static void
probe_initial_source_residual(const struct tov_static_ctx *ctx)
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_gr_polytropic_euler_new(
    ctx->K_poly, ctx->gas_gamma, ctx->spacetime_gauge, ctx->reinit_freq,
    ctx->spacetime, false
  );

  double x_probe[5] = {
    ctx->pos_x,
    ctx->pos_x + 0.5 * ctx->R_star,
    ctx->pos_x + 0.9 * ctx->R_star,
    ctx->pos_x + 0.99 * ctx->R_star,
    ctx->pos_x + 1.01 * ctx->R_star
  };

  printf("\n=== INITIAL SOURCE RESIDUAL PROBE ===\n");
  printf("%-12s %-12s %-12s %-12s %-12s %-12s\n",
    "x", "Etot", "mom_x", "lapse", "S_E", "S_mx");

  for (int p = 0; p < 5; p++) {
    double x = x_probe[p];
    double xn[3] = { x, 0.0, 0.0 };
    double q[70] = { 0.0 };
    double s[70] = { 0.0 };

    evalGRPolytropicEulerInit(0.0, xn, q, (void*) ctx);
    gkyl_wv_eqn_source(eqn, q, s);

    printf("%-12.6f %-12.4e %-12.4e %-12.6f %-12.4e %-12.4e\n",
      x, q[0], q[1], q[4], s[0], s[1]);
  }

  gkyl_wv_eqn_release(eqn);
}

static void
probe_geometry_derivatives(const struct tov_static_ctx *ctx)
{
  double x_probe[4] = {
    ctx->pos_x + 0.5 * ctx->R_star,
    ctx->pos_x + 0.9 * ctx->R_star,
    ctx->pos_x + 0.99 * ctx->R_star,
    ctx->pos_x + 1.01 * ctx->R_star
  };

  double h_probe[3] = { 1.0e-4, 1.0e-6, 1.0e-8 };

  printf("\n=== GEOMETRY DERIVATIVE PROBE ===\n");

  for (int p = 0; p < 4; p++) {
    double x = x_probe[p];
    printf("x = %.8f\n", x);

    for (int m = 0; m < 3; m++) {
      double h = h_probe[m];

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      double ***metric_der = gkyl_malloc(sizeof(double**[3]));
      double **K = gkyl_malloc(sizeof(double*[3]));

      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
        metric_der[i] = gkyl_malloc(sizeof(double*[3]));
        K[i] = gkyl_malloc(sizeof(double[3]));
        for (int j = 0; j < 3; j++) {
          metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      ctx->spacetime->lapse_function_der_func(ctx->spacetime, 0.0, x, 0.0, 0.0, h, h, h, &lapse_der);
      ctx->spacetime->shift_vector_der_func(ctx->spacetime, 0.0, x, 0.0, 0.0, h, h, h, &shift_der);
      ctx->spacetime->spatial_metric_tensor_der_func(ctx->spacetime, 0.0, x, 0.0, 0.0, h, h, h, &metric_der);
      ctx->spacetime->extrinsic_curvature_tensor_func(ctx->spacetime, 0.0, x, 0.0, 0.0, h, h, h, &K);

      printf("  h=%.1e dalpha_dx=%.6e dbetax_dx=%.6e dgamxx_dx=%.6e Kxx=%.6e\n",
        h, lapse_der[0], shift_der[0][0], metric_der[0][0][0], K[0][0]);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          gkyl_free(metric_der[i][j]);
        }
        gkyl_free(shift_der[i]);
        gkyl_free(metric_der[i]);
        gkyl_free(K[i]);
      }
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(metric_der);
      gkyl_free(K);
    }
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

  struct tov_static_ctx ctx = create_ctx();

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  struct gkyl_wv_eqn *gr_polytropic_euler = gkyl_wv_gr_polytropic_euler_new(ctx.K_poly, ctx.gas_gamma, ctx.spacetime_gauge, ctx.reinit_freq, ctx.spacetime, app_args.use_gpu);
  // struct gkyl_wv_eqn *gr_polytropic_euler = gkyl_wv_gr_polytropic_euler_new(
  //   ctx.K_poly, ctx.gas_gamma, ctx.spacetime_gauge, ctx.reinit_freq,
  //   ctx.spacetime, ctx.tov, ctx.pos_x, ctx.pos_y, ctx.pos_z,
  //   app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "gr_polytropic_euler",
    .equation = gr_polytropic_euler,
    .init = evalGRPolytropicEulerInit,
    .force_low_order_flux = true,
    .ctx = &ctx,
    .has_gr_polytropic = true,
    .gr_polytropic_gas_gamma = ctx.gas_gamma,
    .gr_polytropic_K_poly = ctx.K_poly,
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int cells[] = { NX };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) cuts[d] = app_args.cuts[d];
    else cuts[d] = 1;
  }
#else
  for (int d = 0; d < dim; d++) cuts[d] = 1;
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

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) ncuts *= cuts[d];

  if (ncuts != comm_size) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    goto mpifinalize;
  }

  struct gkyl_moment app_inp = {
    .name = "gr_tov_static",
    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,
    .cfl_frac = ctx.cfl_frac,
    .num_species = 1,
    .species = { fluid },
    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end;

  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);
    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }
    frame_curr = status.frame;
    t_curr = status.stime;
    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d at time = %g\n", frame_curr, t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
    print_tov_force_balance_diag(app, &ctx);
    //probe_initial_source_residual(&ctx);
    // probe_geometry_derivatives(&ctx);
    // probe_tov_spacetime_balance(&ctx);
    // probe_tov_cks_metric_identities(&ctx);
    // probe_extrinsic_curvature_consistency(&ctx);
    // probe_stress_energy_consistency(&ctx);
    //probe_flux_source_balance(&ctx);
    //probe_shifted_static_balance(&ctx);
    //probe_low_order_face_flux_detail(&ctx, false);
    //probe_well_balanced_source_detail(&ctx, false);
    // probe_low_order_face_flux_detail(&ctx, true);
    // probe_well_balanced_source_detail(&ctx, true);
    // probe_shifted_static_surface_detail(&ctx);
    //apply_atmosphere_reset(app, &ctx);
  }

  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };
  calc_field_energy(&fe_trig, app, t_curr, false);

  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };
  calc_integrated_mom(&im_trig, app, t_curr, false);

  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };
  write_data(&io_trig, app, t_curr, false);

  double dt = fmin(gkyl_moment_app_max_dt(app), t_end - t_curr); //t_end - t_curr; //
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
    //apply_atmosphere_reset(app, &ctx);
    dt = status.dt_suggested;

    report_first_bad_cell(app, &ctx, t_curr);

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;
      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g is below %g*dt_init ... num_failures = %d\n",
        status.dt_actual, dt_failure_tol, num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init %d consecutive times. Aborting simulation ....\n",
          dt_failure_tol, num_failures_max);
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
  
  //dump_lineout(app, NX);
  //compute_exact_error_norms(app, &ctx);

  // if (NX == 256) {
  // compute_self_convergence("/tmp/rt_gr_tov_static_lineout_128.dat",
  //                           "/tmp/rt_gr_tov_static_lineout_256.dat");
  // }
  // if (NX == 512) {
  //   compute_self_convergence("/tmp/rt_gr_tov_static_lineout_256.dat",
  //                           "/tmp/rt_gr_tov_static_lineout_512.dat");
  // }

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
  gkyl_wv_eqn_release(gr_polytropic_euler);
  gkyl_gr_spacetime_release(ctx.spacetime);
  gkyl_tov_solution_release(ctx.tov);
  //gkyl_gr_spacetime_release(ctx.bh_test);
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
