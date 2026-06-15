// GR Euler Ultra Rel in Spherical Coorindates in 1D for a static TOV
// Euler eqns from: https://arxiv.org/pdf/gr-qc/9904052

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
#include <gkyl_wv_gr_tov.h>
#include "tov_solver.h"

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPIh
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct gr_tov_static_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double kappa; // Stress-energy prefactor in the Einstein field equations.

  double K_poly;
  double rho_c;
  double dr_tov;

  double rho_atm;
  double p_atm;

  struct gkyl_tov *tov;
  double M_star;
  double R_star;

  // Simulation parameters.
  int Nx; // Cell count (r-direction).
  double Lx; // Domain size (r-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gr_tov_static_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 2.0;//4.0 / 3.0; // Adiabatic index.
  double kappa = 8.0 * pi; // Stress-energy prefactor in the Einstein field equations.

  double K_poly = 100.0; 
  double rho_c = 1.28e-3;//5e-8;
  double dr_tov = 0.01;

  double rho_atm = 1e-10 * rho_c;
  //double p_atm = (gas_gamma - 1.0) * rho_atm;
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);

  double M_star = gkyl_tov_star_mass(tov);
  printf("M_star = %e \n", M_star);
  double R_star = gkyl_tov_star_radius(tov);
  printf("R_star = %e \n", R_star);
  printf("Compactness (2M_star / R_star) = %e \n", 2 * M_star / R_star);

  // Simulation parameters.
  int Nx = 4096; // Cell count (r-direction).
  double Lx = 20.0; // Domain size (r-direction).
  double cfl_frac = 0.8; // CFL coefficient.

  double t_end = 6000; // Final simulation time.
  int num_frames = 1000; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gr_tov_static_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .kappa = kappa,
    .K_poly = K_poly,
    .rho_c = rho_c,
    .dr_tov = dr_tov,
    .rho_atm = rho_atm,
    .p_atm = p_atm,
    .tov = tov,
    .M_star = M_star,
    .R_star = R_star,
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
evalGRTovInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double r = xn[0]; // r_relative
  struct gr_tov_static_ctx *app = ctx;

  struct gkyl_tov *tov = app->tov;
  
  struct tov_eval_fluid fluid = {0}; 

  struct tov_eval_bl bl;
  gkyl_tov_eval_bl(tov, r, &bl);

  double lapse = exp(bl.Phi);
  double a = 1.0 / sqrt(1.0 - (2.0 * bl.m / r));

  double v_con[3] = {0.0, 0.0, 0.0}; //fluid coordinate velocity

  double v_sq = 0.0; // For a static TOV
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     v_sq += spatial_metric[i][j] * v_con[i] * v_con[j];
  //   }
  // }

  double W = 1.0 / (sqrt(1.0 - v_sq)); // Lorentz gamma function
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double u_con[4] = {0.0}; //fluid four-velocity
  double shift[3] = {0.0, 0.0, 0.0};
  u_con[0] = W / lapse;
  for (int i = 1; i < 4; i ++) {
    u_con[i] = 0.0; // for a static TOV
  }

  double p = fmax(bl.P, app->p_atm);
  //double rho = p / (app->gas_gamma - 1.0); // fluid-frame total energy density
  double rho = pow(p / app->K_poly, 1.0 / app->gas_gamma);
  double e = rho + p / (app->gas_gamma - 1.0);

  double D = rho * W; // baryonic mass density;
  double Etot = ((e + p) * W * W) - p; // Eulerian-frame (conserved) energy dneisty 
  double mom_r = ((e + p) * W * W) * v_con[0]; // radial momenutm density

  //Evolved conservative variables: D, tau, momentum_r, lapse

  double r2 = r * r;
  fout[0] = r2 * D;
  fout[1] = r2 * (Etot - D); // tau = Etot - D (Valencia energy, rest mass removed)
  fout[2] = r2 * mom_r;
  fout[3] = bl.Phi; // Set Phi appearing in dt metric term
  fout[4] = bl.m;
  fout[5] = r;
  // Frozen discrete well-balanced reference: we store the projected t=0 cell state (r^2*p and r^2*D) in the conserved slots q[6], q[7]. 
  // These ride with the cell and are r^2-averaged by the FV projection IDENTICALLY to q[0] and q[1] (by r^2 over the cell, not the center), so the recovered equilibrium (p_eq=q[6]/r^2, D_eq=q[7]/r^2) matches the current state exactly at t=0. 
  fout[6] = r2 * p;
  fout[7] = r2 * D;

  // for (int i = 0; i < 3; i++) {
  //   gkyl_free(spatial_metric[i]);
  // }
  // gkyl_free(spatial_metric);
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

static double
recovery_r2(double r)
{
  double r2 = r * r;
  if (r2 < 1.0e-300) {
    r2 = 1.0e-300;
  }
  return r2;
}

static void
recover_idealgas_prim(double gas_gamma, const double q[8], double *rho, double *p, double *vel)
{
  double r = q[5];
  double r2 = recovery_r2(r); // cell-center r^2, floored

  double D = q[0] / r2;
  double tau = q[1] / r2;   // tau = Etot - D (Valencia energy)
  double mom_r = q[2] / r2;

  double rho_floor = 1.0e-15;
  double p_floor = 1.0e-300;
  double D_safe = fmax(D, rho_floor);
  double p_new = fmax((gas_gamma - 1.0) * tau, p_floor);
  double p_min = fmax(p_floor, fabs(mom_r) - tau - D_safe + 1.0e-16);
  if (p_new < p_min) {
    p_new = p_min;
  }
  double v = 0.0;
  double rho_new = D_safe;

  for (int iter = 0; iter < 100; iter++) {
    double Q = tau + D_safe + p_new;
    v = mom_r / Q;
    if (v * v > 1.0 - 1.0e-12) {
      v = copysign(sqrt(1.0 - 1.0e-12), v);
    }

    double W2 = 1.0 / (1.0 - v * v);
    double W = sqrt(W2);
    rho_new = D_safe / W;
    double eps = (tau - D_safe * (W - 1.0) - p_new * (W2 - 1.0)) / (D_safe * W);
    if (eps < 0.0) {
      eps = 0.0;
    }

    double p_eos = (gas_gamma - 1.0) * rho_new * eps;
    double f_val = p_eos - p_new;
    double h = 1.0 + eps + p_new / rho_new;
    double cs2 = gas_gamma * p_new / (rho_new * h);
    if (cs2 < 0.0) {
      cs2 = 0.0;
    }
    if (cs2 > 1.0 - 1.0e-12) {
      cs2 = 1.0 - 1.0e-12;
    }
    double f_prime = v * v * cs2 - 1.0;

    double dp = -f_val / f_prime;
    p_new += dp;
    if (p_new < p_min) {
      p_new = p_min;
    }
    if (fabs(dp) < 1.0e-12 * fmax(fabs(p_new), p_floor)) {
      break;
    }
  }

  double denom = tau + D_safe + p_new;
  v = mom_r / denom;
  if (v * v > 1.0 - 1.0e-12) {
    v = copysign(sqrt(1.0 - 1.0e-12), v);
  }

  double W = 1.0 / sqrt(1.0 - v * v);
  *rho = D_safe / W;
  *p = p_new;
  *vel = v;
}

static double
compactness_diagnostic(gkyl_moment_app *app, double t_curr, bool print_crossing)
{
  struct gkyl_array *q = gkyl_moment_app_get_write_array_species(app, 0);

  double max_compact_local = -DBL_MAX;
  double r_at_max_local = 0.0;
  double m_at_max_local = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qcell = gkyl_array_cfetch(q, loc);

    double m = qcell[4];
    double r = qcell[5];
    if (r <= 0.0) {
      continue;
    }

    double compact = 2.0 * m / r;
    if (compact > max_compact_local) {
      max_compact_local = compact;
      r_at_max_local = r;
      m_at_max_local = m;
    }
  }

  double max_compact_global = max_compact_local;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &max_compact_local, &max_compact_global);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0 && print_crossing) {
    printf("Compactness diagnostic at t = %.16e: max(2m/r) = %.16e\n",
       t_curr, max_compact_global);
    //  printf("Compactness diagnostic at t = %.16e: max(2m/r) = %.16e, r = %.16e, m = %.16e\n",
    //    t_curr, max_compact_global, r_at_max_local, m_at_max_local);
    if (max_compact_global > 0.95) {
      printf("Collapsed! max(2m/r) = %.16e crossed 0.95 at t = %.16e\n",
        max_compact_global, t_curr);
    }
  }

  return max_compact_global;
}

static void
polytropic_state_health_diagnostic(const struct gr_tov_static_ctx *ctx, gkyl_moment_app *app, double t_curr, bool force_print)
{
  struct gkyl_array *q = gkyl_moment_app_get_write_array_species(app, 0);

  double min_D = DBL_MAX, r_min_D = 0.0;
  double min_enthalpy = DBL_MAX, r_min_enthalpy = 0.0;
  double min_invariant = DBL_MAX, r_min_invariant = 0.0;
  double max_abs_vel = 0.0, r_max_abs_vel = 0.0, vel_at_max = 0.0;
  double max_abs_vel_int = 0.0, r_max_abs_vel_int = 0.0; // interior only (r < 0.85 R_star)
  double min_p = DBL_MAX, r_min_p = 0.0;
  double R_star = gkyl_tov_star_radius(ctx->tov);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qcell = gkyl_array_cfetch(q, loc);

    double rho, p, vel;
    recover_idealgas_prim(ctx->gas_gamma, qcell, &rho, &p, &vel);

    double r = qcell[5];
    double r2 = r * r;
    if (r2 < 1.0e-300) {
      r2 = 1.0e-300;
    }

    double D = qcell[0] / r2;
    double tau = qcell[1] / r2;
    double Etot = tau + D; // physical energy, restored from tau
    double mom_r = qcell[2] / r2;
    double enthalpy = Etot + p;
    double invariant = (Etot * Etot) - (mom_r * mom_r);

    if (D < min_D) {
      min_D = D;
      r_min_D = r;
    }
    if (p < min_p) {
      min_p = p;
      r_min_p = r;
    }
    if (enthalpy < min_enthalpy) {
      min_enthalpy = enthalpy;
      r_min_enthalpy = r;
    }
    if (invariant < min_invariant) {
      min_invariant = invariant;
      r_min_invariant = r;
    }
    if (fabs(vel) > max_abs_vel) {
      max_abs_vel = fabs(vel);
      vel_at_max = vel;
      r_max_abs_vel = r;
    }
    if (r > 0.0 && r < 0.85 * R_star && fabs(vel) > max_abs_vel_int) {
      max_abs_vel_int = fabs(vel);
      r_max_abs_vel_int = r;
    }
  }

  bool should_print = force_print || max_abs_vel > 0.4 || min_D < 0.0 || min_enthalpy < 0.0 || min_invariant < 0.0;
  if (should_print) {
    printf("Polytropic state health at t = %.16e:\n", t_curr);
    printf("  max |v| = %.16e (v = %.16e) at r = %.16e\n", max_abs_vel, vel_at_max, r_max_abs_vel);
    printf("  max |v| interior (r < 0.85 R_star = %.4e) = %.16e at r = %.16e\n", 0.85 * R_star, max_abs_vel_int, r_max_abs_vel_int);
    printf("  min D = %.16e at r = %.16e\n", min_D, r_min_D);
    printf("  min p = %.16e at r = %.16e\n", min_p, r_min_p);
    printf("  min Etot+p = %.16e at r = %.16e\n", min_enthalpy, r_min_enthalpy);
    printf("  min Etot^2-S^2 = %.16e at r = %.16e\n", min_invariant, r_min_invariant);
  }
}

static void
write_energy_exact_error_norms(const struct gr_tov_static_ctx *ctx, gkyl_moment_app *app, double t_curr, const char *fname)
{
  struct gkyl_array *q = gkyl_moment_app_get_write_array_species(app, 0);
  double dx = app->grid.dx[0];

  double l1_local = 0.0;
  double l2sq_local = 0.0;
  double linf_local = 0.0;
  double l1_int_local = 0.0;
  double l2sq_int_local = 0.0;
  double linf_int_local = 0.0;
  double r_int_cut = 0.8 * ctx->R_star;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qnum = gkyl_array_cfetch(q, loc);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    double qexact[8] = { 0.0 };
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    evalGRTovInit(0.0, xc, qexact, (void*) ctx);

    double r = qnum[5];
    double r2 = r * r;
    if (r2 < 1.0e-300) {
      r2 = 1.0e-300;
    }
    double rexact = qexact[5];
    double rexact2 = rexact * rexact;
    if (rexact2 < 1.0e-300) {
      rexact2 = 1.0e-300;
    }

    double tau_num = qnum[1] / r2;
    double tau_exact = qexact[1] / rexact2;
    double err = tau_num - tau_exact;
    double abs_err = fabs(err);

    l1_local += dx * abs_err;
    l2sq_local += dx * err * err;
    linf_local = fmax(linf_local, abs_err);
    if (xc[0] <= r_int_cut) {
      l1_int_local += dx * abs_err;
      l2sq_int_local += dx * err * err;
      linf_int_local = fmax(linf_int_local, abs_err);
  }
  }

  double sum_local[2] = { l1_local, l2sq_local };
  double sum_global[2] = { 0.0, 0.0 };
  double linf_global = 0.0;
  double sum_int_local[2] = { l1_int_local, l2sq_int_local };
  double sum_int_global[2] = { 0.0, 0.0 };
  double linf_int_global = 0.0;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2, sum_local, sum_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &linf_local, &linf_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2, sum_int_local, sum_int_global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &linf_int_local, &linf_int_global);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    bool write_header = false;
    FILE *check = fopen(fname, "r");
    if (check == NULL) {
      write_header = true;
    }
    else {
      fclose(check);
    }

    FILE *fp = fopen(fname, "a");
    if (fp == NULL) {
      fprintf(stderr, "Could not open %s for convergence diagnostics.\n", fname);
      return;
    }

    if (write_header) {
      fprintf(fp, "# kind nx nx_fine time L1 L2 Linf\n");
    }
    fprintf(fp, "exact %d 0 %.16e %.16e %.16e %.16e\n",
      app->grid.cells[0], t_curr, sum_global[0], sqrt(sum_global[1]), linf_global);
    fprintf(fp, "interior_exact %d 0 %.16e %.16e %.16e %.16e\n",
      app->grid.cells[0], t_curr, sum_int_global[0], sqrt(sum_int_global[1]), linf_int_global);
    fclose(fp);

    printf("\nFinal Etot exact-error norms written to %s:\n", fname);
    printf("  NX = %d, t = %.16e\n", app->grid.cells[0], t_curr);
    printf("  Global:   L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_global[0], sqrt(sum_global[1]), linf_global);
    printf("  Interior r <= 0.8 R_star = %.16e:\n", r_int_cut);
    printf("            L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_int_global[0], sqrt(sum_int_global[1]), linf_int_global);
    printf("  L1 = %.16e, L2 = %.16e, Linf = %.16e\n",
      sum_global[0], sqrt(sum_global[1]), linf_global);
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

  struct gr_tov_static_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  double lower = 0.5 * (ctx.Lx / NX);
  double upper = ctx.Lx + 0.5 * (ctx.Lx / NX);

  // Fluid equations.
  // HLL is the default here. HLLC (WV_GR_TOV_RP_HLLC) is implemented and flux
  // conservative, but for the *static* star it destabilizes the center: removing
  // HLL's contact dissipation lets the unbalanced flux/source residual at r -> 0
  // drive v -> 1. The kick is a well-balancing problem, addressed by equilibrium-
  // subtracted reconstruction, not by the Riemann solver.
  struct gkyl_wv_eqn *gr_tov = gkyl_wv_gr_tov_inew(&(struct gkyl_wv_gr_tov_inp) {
    .gas_gamma = ctx.gas_gamma,
    .kappa = ctx.kappa,
    .rp_type = WV_GR_TOV_RP_HLL_WB,
    .use_gpu = app_args.use_gpu,
  });

  struct gkyl_moment_species fluid = {
    .name = "gr_tov",
    .equation = gr_tov,
    
    .init = evalGRTovInit,
    .ctx = &ctx,

    .has_gr_tov = true,
    .tov_gas_gamma = ctx.gas_gamma,
    .tov_kappa = ctx.kappa,
    .tov_p_atm = ctx.p_atm,
    .has_dynamic_lapse = true,

    .force_low_order_flux = false,
    .limiter = GKYL_MIN_MOD,

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
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
    .name = "gr_einstein_tov_static",

    .ndim = 1,
    .lower = { lower },
    .upper = { upper },
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
  bool collapse_happened = false;
  double max_compactness_seen = compactness_diagnostic(app, t_curr, true);
  polytropic_state_health_diagnostic(&ctx, app, t_curr, true);
  if (max_compactness_seen > 0.95) {
    collapse_happened = true;
  }

  // Compute initial guess of maximum stable time-step.
  double dt = fmin(t_end - t_curr, gkyl_moment_app_max_dt(app)); //double dt = t_end - t_curr;

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
      polytropic_state_health_diagnostic(&ctx, app, t_curr, true);
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);
    double max_compactness = compactness_diagnostic(app, t_curr, !collapse_happened);
    polytropic_state_health_diagnostic(&ctx, app, t_curr, false);
    max_compactness_seen = fmax(max_compactness_seen, max_compactness);
    if (max_compactness > 0.95) {
      collapse_happened = true;
    }

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
        double max_compactness = compactness_diagnostic(app, t_curr, !collapse_happened);
        polytropic_state_health_diagnostic(&ctx, app, t_curr, true);
        max_compactness_seen = fmax(max_compactness_seen, max_compactness);
        if (max_compactness > 0.95) {
          collapse_happened = true;
        }

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
  double max_compactness = compactness_diagnostic(app, t_curr, !collapse_happened);
  polytropic_state_health_diagnostic(&ctx, app, t_curr, true);
  max_compactness_seen = fmax(max_compactness_seen, max_compactness);
  if (max_compactness > 0.95) {
    collapse_happened = true;
  }
  write_energy_exact_error_norms(&ctx, app, t_curr, "gr_tov_static_energy_conv.dat");
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);
  if (collapse_happened) {
    gkyl_moment_app_cout(app, stdout, "Collapsed! max(2m/r) crossed 0.95 during the run. max seen = %.16e\n",
      max_compactness_seen);
  }
  else {
    gkyl_moment_app_cout(app, stdout, "No collapse detected. max(2m/r) seen = %.16e\n",
      max_compactness_seen);
  }

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(gr_tov);
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
