// Static TOV  star in Cartesian cooridnates: GR Euler (gamma-law) fluid on a fixed analytic TOV spacetime background 
// provided by the gr_tov spacetime object (interior TOV + exterior Schwarzschild, regular at the origin), and is not evolved (GKYL_STATIC_GAUGE). 

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
#include <gkyl_gr_tov.h>
#include "tov_solver.h"

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>


struct gr_tov_cart_ctx
{
  double gas_gamma; // Adiabatic index (evolution EOS).
  double K_poly; // Polytropic constant of the cold-IC EOS (p = K rho^gamma) - later gamma-law EOS
  double rho_c; // Central rest-mass density.
  double dr_tov; // TOV solver radial step.
  double rho_atm; // Atmosphere rest-mass density.
  double p_atm; // Atmosphere pressure floor.

  struct gkyl_tov *tov; // Frozen TOV table
  struct gkyl_gr_spacetime *spacetime; // TOV spacetime object (fixed background for now)
  double M_star;
  double R_star;

  enum gkyl_spacetime_gauge spacetime_gauge;
  int reinit_freq;

  int Nx, Ny;
  double Lx, Ly; // Domain is [-Lx/2, Lx/2] x [-Ly/2, Ly/2], star centered at origin!!! Important not to get thtat confused for laters
  double cfl_frac;

  double t_end;
  int num_frames;
  int field_energy_calcs;
  int integrated_mom_calcs;
  double dt_failure_tol;
  int num_failures_max;
};

struct gr_tov_cart_ctx
create_ctx(void)
{
  double gas_gamma = 2.0;
  double K_poly = 100.0;
  double rho_c = 1.28e-3;
  double dr_tov = 0.01;

  double rho_atm = 1e-7 * rho_c; // ETK rho_rel_min = 1e-7 (= 1.28e-10 for this canonical TOV, with this rho_c)
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);
  double M_star = gkyl_tov_star_mass(tov);
  double R_star = gkyl_tov_star_radius(tov);
  printf("M_star = %e \n", M_star);
  printf("R_star = %e \n", R_star);
  printf("Compactness (2M_star / R_star) = %e \n", 2.0 * M_star / R_star);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, 0.0, 0.0, 0.0);

  struct gr_tov_cart_ctx ctx = {
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
    .Nx = 240,
    .Ny = 240,
    .Lx = 24.0,
    .Ly = 24.0,
    .cfl_frac = 0.8,
    .t_end = 50.0,
    .num_frames = 10,
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
  struct gr_tov_cart_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;


  double r = sqrt((x * x) + (y * y)); // this r is measured from the ORIGIN. 
  // The spacetime object measures r from its center `pos` (gkyl_gr_tov_spacetime_new(..., pos_x,y,z)). 
  // They agree ONLY because the star is at pos = (0,0,0). If we ever offset the star (pos != 0), 
  // we must subtract the same pos here -> r = sqrt((x-pos_x)^2 + (y-pos_y)^2), or the fluid and the metric will reference different centers.
  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(app->tov, r, &bl);

  double p = fmax(bl.P, app->p_atm);
  double rho = pow(p / app->K_poly, 1.0 / gas_gamma); // cold polytrope inversion
  double u = 0.0; // static star -> zero coordinate velocity

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

  //(Valencia, v = 0 -> W = 1)
  double vel[3] = { u, 0.0, 0.0 };
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
  fout[68] = x; fout[69] = y; fout[70] = 0.0;

  // Frozen discrete well-balancing reference (num_equations = 73): store the t=0 conserved fluid in slots 71,72. 
  fout[71] = rho_rel; // = fout[0]
  fout[72] = Etot;    // = fout[4]

  // Perturbation away from the equilibrium: small x-momentum kick in the stellar interior. 
  // The slots [71],[72] above keep the unperturbed equilibrium, so w = q - q_eq = (0, dS, 0, 0, 0) is a genuine deviation the WB must evolve (not absorb). Pulsation should stay bounded.
  // if (r < 0.5 * app->R_star) {
  //   fout[1] += 1.0e-3 * rho_rel; // ~1e-3 coordinate velocity perturbation
  // }

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

// --- Milestone tracker: is the origin stable + are the conserved totals flat? ---
// At static equilibrium the momentum density should stay ~0 everywhere. We report the
// max |mom| (sqrt(mom_x^2+mom_y^2+mom_z^2)) and WHERE it is (radius), separating a deep-
// core zone (r < 1) from the rest, plus the relative drift of SUM(sqrt(g) D) and
// SUM(Etot). Cartesian success = max|mom| stays tiny and the sums stay ~machine flat.
static double g_sum_D0 = -1.0, g_sum_E0 = -1.0;
static void
milestone_tracker(gkyl_moment_app *app, double t_curr, long step)
{
  struct gkyl_array *q = gkyl_moment_app_get_write_array_species(app, 0);

  double core_mom = 0.0, r_core = 0.0;   // r < 1 (near origin)
  double glob_mom = 0.0, r_glob = 0.0;   // anywhere
  double sum_D = 0.0, sum_E = 0.0;
  double min_inv = DBL_MAX, r_inv = 0.0;  // conserved physicality invariant (D+tau)^2 - S^2 (<0 => unphysical)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qc = gkyl_array_cfetch(q, loc);
    if (qc[27] < 0.0) continue; // excised cell

    sum_D += qc[0];
    sum_E += qc[4];

    double mom = sqrt((qc[1] * qc[1]) + (qc[2] * qc[2]) + (qc[3] * qc[3]));
    double inv = ((qc[0] + qc[4]) * (qc[0] + qc[4])) - (mom * mom); // (D+tau)^2 - S^2 (conserved)
    double r = sqrt((qc[68] * qc[68]) + (qc[69] * qc[69]));
    if (mom > glob_mom) { glob_mom = mom; r_glob = r; }
    if (r < 1.0 && mom > core_mom) { core_mom = mom; r_core = r; }
    if (inv < min_inv) { min_inv = inv; r_inv = r; }
  }
  if (g_sum_D0 < 0.0) { g_sum_D0 = sum_D; g_sum_E0 = sum_E; }

  printf("[CART] step=%5ld t=%.4e | core|mom|=%.3e | glob|mom|=%.3e@%.2f | min(D+t)^2-S^2=%+.3e@%.2f | dSUM_D=%+.3e dSUM_E=%+.3e\n",
    step, t_curr, core_mom, glob_mom, r_glob, min_inv, r_inv,
    (sum_D - g_sum_D0) / g_sum_D0, (sum_E - g_sum_E0) / g_sum_E0);
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

  struct gr_tov_cart_ctx ctx = create_ctx();

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Well-balanced against the frozen TOV equilibrium: pass tov_eq to switch WB on, plus the
  // atmosphere rho/p so the C2P recovery floor matches the IC. Other GR-Euler problems leave
  // tov_eq unset (NULL) and get plain GR Euler.
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
    .gas_gamma = ctx.gas_gamma,
    .spacetime_gauge = ctx.spacetime_gauge,
    .reinit_freq = ctx.reinit_freq,
    .spacetime = ctx.spacetime,
    .rp_type = WV_GR_EULER_RP_LAX, // the Riemann solver the MP scheme uses at each reconstructed edge
    .use_gpu = app_args.use_gpu,
    .tov_eq = ctx.tov, // turns on the static-TOV well-balancing (MP flux-form WB under GKYL_MOMENT_MP)
    .p_atm = ctx.p_atm,
    .rho_atm = ctx.rho_atm, // C2P recovery floor matches the IC atmosphere (ET rho_rel_min=1e-7)
  });

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    .init = evalGREulerInit,
    //.force_low_order_flux is a wave-prop-only flag; we use the MP scheme
    .ctx = &ctx,

    .has_gr_euler = true,
    .gr_euler_gas_gamma = ctx.gas_gamma,

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
    .name = "gr_tov_static_cart",

    .ndim = 2,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly },
    .upper = {  0.5 * ctx.Lx,  0.5 * ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_MP, // MP (Suresh-Huynh) flux-form scheme + well-balanced deviation reconstruction
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

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
  milestone_tracker(app, t_curr, 0);

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
    if (step <= 20 || step % 25 == 0) milestone_tracker(app, t_curr, step); // whole run now

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
  gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
  gkyl_tov_solution_release(ctx.tov);
  gkyl_comm_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Finalize();
#endif

  return 0;
}
