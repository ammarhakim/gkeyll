#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_gr_ultra_rel_euler.h>
#include <gkyl_gr_minkowski.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct coupled_einstein_plane_shock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  //double tau0; // Reparameterized initial time.
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor in the Einstein field equations.

  double exp_2a; // Exponential appearing in dt and dx metric terms.
  
  // Derived physical quantities (using normalized code units).
  double rhol; // Left fluid mass density.
  double rhor; // Right fluid mass density.

  double Etot_l; // Left fluid total energy density.
  double Etot_r; // Right fluid total energy density.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Evolution parameters.
  double excision_threshold; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution; // Spacetime evolution system.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct coupled_einstein_plane_shock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  //double tau0 = 9.8753205829098; // Reparameterized initial time.
  double gas_gamma = 4.0 / 3.0; // Adiabatic index.

  double kappa = 8.0 * pi; // Stress-energy prefactor in the Einstein field equations.

  //double a_init = 0.0;
  //double exp_2a = exp(2.0 * a_init);
  double exp_2a = exp(-4.0); // Exponential appearing in dt and dx metric terms.

  // Derived physical quantities (using normalized code units).
  double rhol = 100.0 / kappa; // Left fluid mass density.
  double rhor = 1.0 / kappa; // Right fluid mass density.

  double Etot_l = rhol; // Left fluid total energy density.
  double Etot_r = rhor; // Right fluid total energy density.

  // Pointer to spacetime metric (placeholder; actual metric is set analytically in init functions).
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  // Evolution parameters.
  double excision_threshold = 0.3; // Excision threshold (lapse).
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_HARMONIC_SLICING; // Spacetime slicing condition.
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION; // Spacetime evolution system.

  enum gkyl_spacetime_gauge spacetime_gauge = GKYL_STATIC_GAUGE; // Spacetime gauge choice.
  int reinit_freq = 100; // Spacetime reinitialization frequency.

  // Simulation parameters.
  int Nx = 800; // Cell count (x-direction).
  double Lx = 2.0; // Domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.5; // Final simulation time.
  int num_frames = 100; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct coupled_einstein_plane_shock_ctx ctx = {
    .pi = pi,
    //.tau0 = tau0,
    .gas_gamma = gas_gamma,
    .kappa = kappa,
    .exp_2a = exp_2a,
    .rhol = rhol,
    .rhor = rhor,
    .Etot_l = Etot_l,
    .Etot_r = Etot_r,
    .spacetime = spacetime,
    .spacetime_gauge = spacetime_gauge,
    .reinit_freq = reinit_freq,
    .excision_threshold = excision_threshold,
    .spacetime_slicing = spacetime_slicing,
    .spacetime_evolution = spacetime_evolution,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
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
evalVacuumEinsteinInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct coupled_einstein_plane_shock_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double kappa = app->kappa;
  double exp_2a = app->exp_2a;
  double exp_a = sqrt(exp_2a);
  double rhol = app->rhol;
  double rhor = app->rhor;
  double Etot_l = app->Etot_l;
  double Etot_r = app->Etot_r;

  double Etot = 0.0;


  if (x < 0.0) {
    Etot = Etot_l; // Fluid total energy density (left).
  }
  else {
    Etot = Etot_r; // Fluid total energy density (right).
  }

  
  double f = kappa * exp_2a * Etot;
  double b = (2.0 / 3.0) * log(fabs(cos(0.5 * x * sqrt(3.0 * f))));
  double c = 0.0;
  double exp_2bpc = exp(2 * ( b + c));
  double exp_2bmc = exp(2 * ( b - c));

  double a_dt = 0.0; // Time derivative of metric term a.
  double a_dx = 0.0; // Space derivative of metric term a.
  double b_dt = 0.0; // Time derivative of metric term b.
  double b_dx = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * x * sqrt(3.0 * kappa * exp_2a * Etot))); // Space derivative of metric term b.
  double c_dt = 0.0; // Time derivative of metric term c.
  double c_dx = 0.0; // Space derivative of metric term c.

  double b_dx_plus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x + (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));
  double b_dx_minus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x - (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));

  double a_dt_dx = 0.0; // Mixed space-time derivative of metric term a.
  double a_dx_dx = 0.0; // Second space derivative of metric term a.
  double b_dt_dx = 0.0; // Mixed space-time derivative of metric term b.
  double b_dx_dx = (b_dx_plus - b_dx_minus) / pow(10.0, -8.0); // Second space derivative of metric term b.
  double c_dt_dx = 0.0; // Mixed space-time derivative of metric term c.
  double c_dx_dx = 0.0; // Second space derivative of metric term c.

  double mom_x = 0.0; // Fluid momentum (x-direction).

  // Compute metric quantities.
  double lapse;

  double shift[3];
  double spatial_metric[3][3], inv_spatial_metric[3][3];
  double extrinsic_curvature[3][3];
  double lapse_der[3];
  double shift_der[3][3];
  double spatial_metric_der[3][3][3];
  double aux_vect[3];

  lapse = exp_a;
  lapse_der[0] = exp_a * a_dt; lapse_der[1] = 0.0; lapse_der[2] = 0.0;

  shift[0] = 0.0; shift[1] = 0.0; shift[2] = 0.0;
  shift_der[0][0] = 0.0; shift_der[0][1] = 0.0; shift_der[0][2] = 0.0;
  shift_der[1][0] = 0.0; shift_der[1][1] = 0.0; shift_der[1][2] = 0.0;
  shift_der[2][0] = 0.0; shift_der[2][1] = 0.0; shift_der[2][2] = 0.0;

  spatial_metric[0][0] = exp_2a;   spatial_metric[0][1] = 0.0;      spatial_metric[0][2] = 0.0;
  spatial_metric[1][0] = 0.0;      spatial_metric[1][1] = exp_2bpc;  spatial_metric[1][2] = 0.0;
  spatial_metric[2][0] = 0.0;      spatial_metric[2][1] = 0.0;       spatial_metric[2][2] = exp_2bmc;

  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        spatial_metric_der[k][i][j] = 0.0;
      }
    }
  }
  spatial_metric_der[0][0][0] = 2 * spatial_metric[0][0] * a_dx;
  spatial_metric_der[0][1][1] = 2 * spatial_metric[1][1] * (b_dx + c_dx);
  spatial_metric_der[0][2][2] = 2 * spatial_metric[2][2] * (b_dx - c_dx);

  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
      extrinsic_curvature[k][i] = 0.0;
    }
  }
  extrinsic_curvature[0][0] = exp_a * a_dt;
  extrinsic_curvature[1][1] = exp_2bpc * (b_dt + c_dt);
  extrinsic_curvature[2][2] = exp_2bmc * (b_dt - c_dt);

  aux_vect[0] = 2 * b_dx;
  aux_vect[1] = 0.0;
  aux_vect[2] = 0.0;

  // Pack vacuum einstein q[64].
  // Spatial metric tensor.
  fout[0] = spatial_metric[0][0]; fout[1] = spatial_metric[0][1]; fout[2] = spatial_metric[0][2];
  fout[3] = spatial_metric[1][0]; fout[4] = spatial_metric[1][1]; fout[5] = spatial_metric[1][2];
  fout[6] = spatial_metric[2][0]; fout[7] = spatial_metric[2][1]; fout[8] = spatial_metric[2][2];

  // Lapse.
  fout[9] = lapse;

  // Extrinsic curvature tensor.
  fout[10] = extrinsic_curvature[0][0]; fout[11] = extrinsic_curvature[0][1]; fout[12] = extrinsic_curvature[0][2];
  fout[13] = extrinsic_curvature[1][0]; fout[14] = extrinsic_curvature[1][1]; fout[15] = extrinsic_curvature[1][2];
  fout[16] = extrinsic_curvature[2][0]; fout[17] = extrinsic_curvature[2][1]; fout[18] = extrinsic_curvature[2][2];

  // Spatial metric tensor derivatives.
  fout[19] = spatial_metric_der[0][0][0]; fout[20] = spatial_metric_der[0][0][1]; fout[21] = spatial_metric_der[0][0][2];
  fout[22] = spatial_metric_der[0][1][0]; fout[23] = spatial_metric_der[0][1][1]; fout[24] = spatial_metric_der[0][1][2];
  fout[25] = spatial_metric_der[0][2][0]; fout[26] = spatial_metric_der[0][2][1]; fout[27] = spatial_metric_der[0][2][2];

  fout[28] = spatial_metric_der[1][0][0]; fout[29] = spatial_metric_der[1][0][1]; fout[30] = spatial_metric_der[1][0][2];
  fout[31] = spatial_metric_der[1][1][0]; fout[32] = spatial_metric_der[1][1][1]; fout[33] = spatial_metric_der[1][1][2];
  fout[34] = spatial_metric_der[1][2][0]; fout[35] = spatial_metric_der[1][2][1]; fout[36] = spatial_metric_der[1][2][2];

  fout[37] = spatial_metric_der[2][0][0]; fout[38] = spatial_metric_der[2][0][1]; fout[39] = spatial_metric_der[2][0][2];
  fout[40] = spatial_metric_der[2][1][0]; fout[41] = spatial_metric_der[2][1][1]; fout[42] = spatial_metric_der[2][1][2];
  fout[43] = spatial_metric_der[2][2][0]; fout[44] = spatial_metric_der[2][2][1]; fout[45] = spatial_metric_der[2][2][2];

  // Lapse function derivatives.
  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];

  // Auxiliary vector.
  fout[49] = aux_vect[0]; fout[50] = aux_vect[1]; fout[51] = aux_vect[2];

  // Shift gauge variables.
  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];

  // Shift vector derivatives.
  fout[55] = shift_der[0][0]; fout[56] = shift_der[0][1]; fout[57] = shift_der[0][2];
  fout[58] = shift_der[1][0]; fout[59] = shift_der[1][1]; fout[60] = shift_der[1][2];
  fout[61] = shift_der[2][0]; fout[62] = shift_der[2][1]; fout[63] = shift_der[2][2];
}

void
evalGRUltraRelEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct coupled_einstein_plane_shock_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double kappa = app->kappa;
  double exp_2a = app->exp_2a;
  double exp_a = sqrt(exp_2a);
  double rhol = app->rhol;
  double rhor = app->rhor;
  double Etot_l = app->Etot_l;
  double Etot_r = app->Etot_r;

  double Etot = 0.0;


  if (x < 0.0) {
    Etot = Etot_l; // Fluid total energy density (left).
  }
  else {
    Etot = Etot_r; // Fluid total energy density (right).
  }
  
  double f = kappa * exp_2a * Etot;
  double b = (2.0 / 3.0) * log(fabs(cos(0.5 * x * sqrt(3.0 * f))));
  double c = 0.0;
  double exp_2bpc = exp(2 * ( b + c));
  double exp_2bmc = exp(2 * ( b - c));

  double a_dt = 0.0; // Time derivative of metric term a.
  double a_dx = 0.0; // Space derivative of metric term a.
  double b_dt = 0.0; // Time derivative of metric term b.
  double b_dx = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * x * sqrt(3.0 * kappa * exp_2a * Etot))); // Space derivative of metric term b.
  double c_dt = 0.0; // Time derivative of metric term c.
  double c_dx = 0.0; // Space derivative of metric term c.

  double b_dx_plus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x + (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));
  double b_dx_minus = -sqrt((kappa * exp_2a * Etot) / 3.0) * tan((0.5 * (x - (0.5 * pow(10.0, -8.0))) * sqrt(3.0 * kappa * exp_2a * Etot)));

  double a_dt_dx = 0.0; // Mixed space-time derivative of metric term a.
  double a_dx_dx = 0.0; // Second space derivative of metric term a.
  double b_dt_dx = 0.0; // Mixed space-time derivative of metric term b.
  double b_dx_dx = (b_dx_plus - b_dx_minus) / pow(10.0, -8.0); // Second space derivative of metric term b.
  double c_dt_dx = 0.0; // Mixed space-time derivative of metric term c.
  double c_dx_dx = 0.0; // Second space derivative of metric term c.

  double mom_x = 0.0; // Fluid momentum (x-direction).

  // Compute metric quantities.
  double lapse;

  double shift[3];
  double spatial_metric[3][3], inv_spatial_metric[3][3];
  double extrinsic_curvature[3][3];
  double lapse_der[3];
  double shift_der[3][3];
  double spatial_metric_der[3][3][3];
  double aux_vect[3];

  lapse = exp_a;
  lapse_der[0] = exp_a * a_dt; lapse_der[1] = 0.0; lapse_der[2] = 0.0;

  shift[0] = 0.0; shift[1] = 0.0; shift[2] = 0.0;
  shift_der[0][0] = 0.0; shift_der[0][1] = 0.0; shift_der[0][2] = 0.0;
  shift_der[1][0] = 0.0; shift_der[1][1] = 0.0; shift_der[1][2] = 0.0;
  shift_der[2][0] = 0.0; shift_der[2][1] = 0.0; shift_der[2][2] = 0.0;

  spatial_metric[0][0] = exp_2a;   spatial_metric[0][1] = 0.0;      spatial_metric[0][2] = 0.0;
spatial_metric[1][0] = 0.0;      spatial_metric[1][1] = exp_2bpc;  spatial_metric[1][2] = 0.0;
spatial_metric[2][0] = 0.0;      spatial_metric[2][1] = 0.0;       spatial_metric[2][2] = exp_2bmc;

  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        spatial_metric_der[k][i][j] = 0.0;
      }
    }
  }
  spatial_metric_der[0][0][0] = 2 * spatial_metric[0][0] * a_dx;
  spatial_metric_der[0][1][1] = 2 * spatial_metric[1][1] * (b_dx + c_dx);
  spatial_metric_der[0][2][2] = 2 * spatial_metric[2][2] * (b_dx - c_dx);

  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
      extrinsic_curvature[k][i] = 0.0;
    }
  }
  extrinsic_curvature[0][0] = exp_a * a_dt;
  extrinsic_curvature[1][1] = exp_2bpc * (b_dt + c_dt);
  extrinsic_curvature[2][2] = exp_2bmc * (b_dt - c_dt);

  aux_vect[0] = 2 * b_dx;
  aux_vect[1] = 0.0;
  aux_vect[2] = 0.0;



  // double spatial_det = spatial_metric[0][0] * spatial_metric[1][1] * spatial_metric[2][2];
  // double p = (gas_gamma - 1.0) * rho;
  // double W = 1.0;
  // double Etot = sqrt(spatial_det) * (((rho + p) * (W * W)) - p) // Fluid total energy density.
  fout[0] = Etot;
  fout[1] = mom_x;  
  fout[2] = 0; 
  fout[3] = 0;

  // Lapse gauge variable.
  fout[4] = lapse;
  // Shift gauge variables.
  fout[5] = shift[0]; fout[6] = shift[1]; fout[7] = shift[2];

  // Spatial metric tensor.
  fout[8]  = spatial_metric[0][0]; fout[9]  = spatial_metric[0][1]; fout[10] = spatial_metric[0][2];
  fout[11] = spatial_metric[1][0]; fout[12] = spatial_metric[1][1]; fout[13] = spatial_metric[1][2];
  fout[14] = spatial_metric[2][0]; fout[15] = spatial_metric[2][1]; fout[16] = spatial_metric[2][2];

  // Extrinsic curvature tensor.
  fout[17] = extrinsic_curvature[0][0]; fout[18] = extrinsic_curvature[0][1]; fout[19] = extrinsic_curvature[0][2];
  fout[20] = extrinsic_curvature[1][0]; fout[21] = extrinsic_curvature[1][1]; fout[22] = extrinsic_curvature[1][2];
  fout[23] = extrinsic_curvature[2][0]; fout[24] = extrinsic_curvature[2][1]; fout[25] = extrinsic_curvature[2][2];

  // Excision boundary conditions.
  fout[26] = 1.0; // Not in excision region.

  // Lapse function derivatives.
  fout[27] = lapse_der[0]; fout[28] = lapse_der[1]; fout[29] = lapse_der[2];
  // Shift vector derivatives.
  fout[30] = shift_der[0][0]; fout[31] = shift_der[0][1]; fout[32] = shift_der[0][2];
  fout[33] = shift_der[1][0]; fout[34] = shift_der[1][1]; fout[35] = shift_der[1][2];
  fout[36] = shift_der[2][0]; fout[37] = shift_der[2][1]; fout[38] = shift_der[2][2];

  // Spatial metric tensor derivatives.
  fout[39] = spatial_metric_der[0][0][0]; fout[40] = spatial_metric_der[0][0][1]; fout[41] = spatial_metric_der[0][0][2];
  fout[42] = spatial_metric_der[0][1][0]; fout[43] = spatial_metric_der[0][1][1]; fout[44] = spatial_metric_der[0][1][2];
  fout[45] = spatial_metric_der[0][2][0]; fout[46] = spatial_metric_der[0][2][1]; fout[47] = spatial_metric_der[0][2][2];

  fout[48] = spatial_metric_der[1][0][0]; fout[49] = spatial_metric_der[1][0][1]; fout[50] = spatial_metric_der[1][0][2];
  fout[51] = spatial_metric_der[1][1][0]; fout[52] = spatial_metric_der[1][1][1]; fout[53] = spatial_metric_der[1][1][2];
  fout[54] = spatial_metric_der[1][2][0]; fout[55] = spatial_metric_der[1][2][1]; fout[56] = spatial_metric_der[1][2][2];

  fout[57] = spatial_metric_der[2][0][0]; fout[58] = spatial_metric_der[2][0][1]; fout[59] = spatial_metric_der[2][0][2];
  fout[60] = spatial_metric_der[2][1][0]; fout[61] = spatial_metric_der[2][1][1]; fout[62] = spatial_metric_der[2][1][2];
  fout[63] = spatial_metric_der[2][2][0]; fout[64] = spatial_metric_der[2][2][1]; fout[65] = spatial_metric_der[2][2][2];

  // Evolution parameter.
  fout[66] = 0.0;

  // Spatial coordinates.
  fout[67] = x; fout[68] = 0.0; fout[69] = 0.0;
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

  struct coupled_einstein_plane_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Einstein equations.
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(ctx.excision_threshold, ctx.spacetime_slicing, ctx.spacetime_evolution, app_args.use_gpu);
  // Fluid equations.
  struct gkyl_wv_eqn *gr_ultra_rel_euler = gkyl_wv_gr_ultra_rel_euler_new(ctx.gas_gamma, ctx.spacetime_gauge, ctx.reinit_freq, ctx.spacetime, app_args.use_gpu);
  // Wire the vacuum einstein equation into the spacetime object.
  ctx.spacetime->vacuum_einstein = vacuum_einstein;

  struct gkyl_moment_species einstein = {
    .name = "vacuum_einstein",
    .equation = vacuum_einstein,
    
    .init = evalVacuumEinsteinInit,
    .force_low_order_flux = false, // Use HLL fluxes.
    .ctx = &ctx,

    .has_vacuum_einstein = true,
    .vacuum_einstein_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_spacetime_evolution = ctx.spacetime_evolution,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species fluid = {
    .name = "gr_ultra_rel_euler",
    .equation = gr_ultra_rel_euler,
    
    .init = evalGRUltraRelEulerInit,
    .force_low_order_flux = true, // Use Lax fluxes.
    .ctx = &ctx,

    .has_gr_ultra_rel = true,
    .gr_ultra_rel_gas_gamma = ctx.gas_gamma,

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
    .name = "coupled_einstein_plane_shock",

    .ndim = 1,
    .lower = { -0.5 * ctx.Lx },
    .upper = { 0.5 * ctx.Lx },
    .cells = { NX },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { fluid, einstein },
    .spacetime = ctx.spacetime,

    .num_periodic_dir = 0, 

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
  gkyl_wv_eqn_release(vacuum_einstein);
  gkyl_wv_eqn_release(gr_ultra_rel_euler);
  gkyl_gr_spacetime_release(ctx.spacetime);
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
