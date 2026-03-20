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
  double R_star = gkyl_tov_star_radius(tov);

  double Lx = 200.0;
  double Ly = 1.0;
  double Lz = 1.0;
  double pos_x = 0.5 * Lx;
  double pos_y = 0.0;
  double pos_z = 0.0;

  // TOV spacetime: TOV interior for r < R_star, Schwarzschild for r >= R_star, no excision.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, pos_x, 0.0, 0.0);
  //struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  int Nx = 512;
  // const char *nx_env = getenv("GKYL_TOV_NX");
  // if (nx_env) {
  //   int nx_in = atoi(nx_env);
  //   if (nx_in > 0) {
  //     Nx = nx_in;
  //   }
  // }
  int Ny = 1;
  int Nz = 1;
  double cfl_frac = 0.5;
  

  enum gkyl_spacetime_gauge spacetime_gauge = GKYL_STATIC_GAUGE;
  int reinit_freq = INT_MAX;
  //printf("%i\n",INT_MAX);

  double t_end = 10.0;
  int num_frames = 10;
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

  double vel[3] = {0.0, 0.0, 0.0};

  // // TOV solution relative to star center in Spherical Kerr-Schild/Schwarzschild coordinates (they share areal radius)
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
     vel[0] = vel[1] = vel[2] = 0.0;
  }

  bool is_atmosphere = (rho < rho_atm) || (p < p_atm) || (r > R_star);

  if (is_atmosphere) {
    rho = rho_atm;
    p = p_atm;
  }

  // Smooth weight: w ~ 1 inside star, w ~ 0 outside.
  // double dcell = app->Lx / app->Nx;
  // double dr_smooth = 2.0 * dcell;
  // dr_smooth = 4.0 * fmax(app->Lx / app->Nx, app->Ly / app->Ny);
  // double w = 0.5 * (1.0 - tanh((r - R_star) / dr_smooth));

  // rho = w * rho + (1.0 - w) * rho_atm;
  // p = w * p + (1.0 - w) * p_atm;

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
  // static int lapse_init_check = 0;
  // if (lapse_init_check < 3 && r > 7.0 && r < 9.0) {
  //     printf("[init lapse] x=%.4f r=%.4f lapse=%.6f\n", x, r, lapse);
  //     lapse_init_check++;
  // }
  // static int lapse_init_check2 = 0;
  // if (lapse_init_check2 < 3 && x > pos_x && r > 7.0 && r < 9.0) {
  //     printf("[init lapse RIGHT] x=%.4f r=%.4f lapse=%.6f\n", x, r, lapse);
  //     lapse_init_check2++;
  // }
  spacetime->shift_vector_func(spacetime, 0.0, x, 0.0, 0.0, &shift);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_metric);
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_det);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, 0.0, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / (sqrt(1.0 - v_sq));
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double eps = is_atmosphere ? (p / ((gas_gamma - 1.0) * rho)) : fluid.eps;
  double h = 1.0 + eps + p / rho;

  double Etot = sqrt(spatial_det) * (((rho * h) * (W * W)) - p); // Fluid total energy density.
  double mom_x = sqrt(spatial_det) * (rho * h) * (W * W) * vel[0]; // Fluid momentum density (x-direction).
  double mom_y = sqrt(spatial_det) * (rho * h) * (W * W) * vel[1];
  double mom_z = sqrt(spatial_det) * (rho * h) * (W * W) * vel[2];
  

  //Fluid diagnostic (from TOV)
  double M_star = gkyl_tov_star_mass(tov);
  R_star = gkyl_tov_star_radius(tov);
  static int tov_debug_printed = 0;
  if (!tov_debug_printed) {
    printf("[TOV] M_star=%.6e R_star=%.6e\n", M_star, R_star);
    // compactness 0.1-0.2 for a realistic NS
    printf("[TOV] compactness 2M/R = %.6f\n", 2.0 * M_star / R_star);

    struct tov_eval_fluid out;
    gkyl_tov_eval(tov, 0.0, 0.0, 0.0, &out);
    printf("[CKS center] rho=%.6e P=%.6e eps=%.6e\n", out.rho, out.P, out.eps);

    gkyl_tov_eval(tov, 0.5 * R_star, 0.0, 0.0, &out);
    printf("[CKS mid] rho=%.6e P=%.6e\n", out.rho, out.P);

    gkyl_tov_eval(tov, 1.1 * R_star, 0.0, 0.0, &out);
    printf("[CKS exterior] rho=%.6e P=%.6e (should be ~0)\n", out.rho, out.P);
    tov_debug_printed = 1;
  }

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

// static void
// apply_atmosphere_reset(gkyl_moment_app* app, const struct tov_static_ctx* ctx)
// {
//   if (app->num_species < 1) {
//     return;
//   }

//   struct gkyl_array *fld = app->species[0].f[0];
//   struct gkyl_range *rng = &app->local;
//   double x_lo = app->grid.lower[0];
//   double dx = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

//   double gas_gamma = ctx->gas_gamma;
//   double K_poly = ctx->K_poly;
//   double rho_atm = ctx->rho_atm;
//   double p_atm = ctx->p_atm;
//   double R_star = ctx->R_star;
//   double pos_x = ctx->pos_x;
//   double eps_atm = p_atm / ((gas_gamma - 1.0) * rho_atm);
//   double h_atm = 1.0 + eps_atm + (p_atm / rho_atm);

//   int idx[GKYL_MAX_DIM] = {0};
//   for (int ci = rng->lower[0]; ci <= rng->upper[0]; ci++) {
//     idx[0] = ci;

//     double x = x_lo + (ci - rng->lower[0] + 0.5) * dx;
//     double r = fabs(x - pos_x);

//     double *q = gkyl_array_fetch(fld, gkyl_range_idx(rng, idx));
//     double v[70] = {0.0};
//     gkyl_gr_polytropic_euler_prim_vars(K_poly, gas_gamma, q, v);

//     double rho = v[0];
//     double p = K_poly * pow(fmax(rho, 0.0), gas_gamma);
//     bool bad_state = isnan(rho) || isnan(p) || rho < rho_atm || p < p_atm;
//     bool exterior = r > R_star;

//     if (exterior || bad_state) {
//       double spatial_metric[3][3];
//       spatial_metric[0][0] = q[8];  spatial_metric[0][1] = q[9];  spatial_metric[0][2] = q[10];
//       spatial_metric[1][0] = q[11]; spatial_metric[1][1] = q[12]; spatial_metric[1][2] = q[13];
//       spatial_metric[2][0] = q[14]; spatial_metric[2][1] = q[15]; spatial_metric[2][2] = q[16];

//       double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2])))
//         - (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0])))
//         + (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

//       q[0] = sqrt(spatial_det) * ((rho_atm * h_atm) - p_atm);
//       q[1] = 0.0;
//       q[2] = 0.0;
//       q[3] = 0.0;
//     }
//   }
// }

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

  double dt = t_end - t_curr;// fmin(gkyl_moment_app_max_dt(app), t_end - t_curr);
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
