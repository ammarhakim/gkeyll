#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Parameters passed to the run_simulation function
// to describe the current phase
struct phase_params {
  double alpha;
  bool static_field;
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type;
  bool positivity;
  bool update_species;
};

// Structure to hold restart information for POA
struct poa_restart_info {
  int cycle;              // Which cycle we're in (0-based)
  bool is_oap_phase;      // true if in OAP phase, false if in FDP phase
  bool is_final_fdp;      // true if in final FDP phase with extra frames
  double phase_start_time; // Start time of current phase
  double phase_end_time;   // End time of current phase
  int phase_start_frame;   // Start frame of current phase
  int phase_end_frame;     // End frame of current phase
};

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p;
  double beta;
  double tau;
  double Ti0;
  double kperpRhos;
  // Parameters controlling initial conditions.
  double alim;
  double nuFrac;
  // Electron-electron collision freq.
  double logLambdaElc;
  double nuElc;
  double elc_nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  // Thermal speeds.
  double vti;
  double vte;
  double c_s;
  // Gyrofrequencies and gyroradii.
  double omega_ci;
  double rho_s;
  double kperp; // Perpendicular wavenumber in SI units.
  double RatZeq0; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_min;
  double psi_eval;
  double psi_max;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double t_end;
  int num_frames;
  double write_phase_freq; // Frequency of writing phase-space data.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude;
  double ion_source_sigma;
  double ion_source_temp;

  // Initial conditions reading
  double *f_dist_ion;
  double *f_dist_elc;
  double *phi_vals;
  double *psi_grid;
  double *z_grid;
  double *v_grid;
  double *theta_grid;
  double *B_grid;
  int *dims;
  int rank;

  // POA parameters
  int oap_first; // Whether to start with OAP phase (1) or FDP phase (0)
  double alpha_oap; // Alpha parameter for OAP phase
  double alpha_fdp; // Alpha parameter for FDP phase
  double tau_oap; // Duration of OAP phase in microseconds
  double tau_fdp; // Duration of FDP phase in microseconds
  double tau_fdp_extra; // Extra FDP duration in final cycle in microseconds
  bool static_field_oap; // Static field setting for OAP phase
  bool static_field_fdp; // Static field setting for FDP phase
  bool update_species_oap; // Species update setting for OAP phase
  bool update_species_fdp; // Species update setting for FDP phase
  double mask_loss_fac_oap; // Mask loss factor for OAP phase
  double mask_loss_fac_fdp; // Mask loss factor for FDP phase
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_oap; // fdot multiplier type for OAP phase
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_fdp; // fdot multiplier type for FDP phase
  bool positivity_oap; // Positivity hack setting for OAP phase
  bool positivity_fdp; // Positivity hack setting for FDP phase
  int num_cycles; // Number of OAP/FDP cycles
  int frames_per_phase; // Number of frames to output per phase
};

// Evaluate collision frequencies
void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = 1e17;
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->Ti0;
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double src_amp = app->ion_source_amplitude;
  double z_src = 0.0;
  double src_sigma = app->ion_source_sigma;
  double src_amp_floor = src_amp*1e-2;
  if (fabs(z) <= 1.0)
  {
    // fout[0] = fmax(src_amp_floor, (src_amp / sqrt(2.0 * M_PI * pow(src_sigma, 2))) *
      // exp(-1 * pow((z - z_src), 2) / (2.0 * pow(src_sigma, 2))));
    
      // cubic polynomial drop of to the edge
    fout[0] = src_amp * (1 - pow(fabs(z), 6));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double TSrc0 = app->ion_source_temp;
  double Tfloor = TSrc0*1e-2;
  if (fabs(z) <= 1.0)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./6.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_ion*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_ion*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;
  double kperpRhos = 0.1;

  // Parameters controlling initial conditions.
  double alim = 0.125;
  double alphaIC0 = 2;
  double alphaIC1 = 10;

  double nuFrac = 1.0;
  double elc_nuFrac = 1/5.489216862238348;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Te0 / eV);
  double nuElc = elc_nuFrac * nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
                 (6. * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(me) * pow(Te0, 3. / 2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Perpendicular wavenumber in SI units:
  double kperp = kperpRhos / rho_s;

  // Geometry parameters.
  double z_min = -2.0;
  double z_max =  2.0;
  double psi_min = 1e-6; // Go smaller. 1e-4 might be too small
  double psi_eval= 1e-3;
  double psi_max = 3e-3; // aim for 2e-2

  // Grid parameters
  double vpar_max_elc = 30 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 30 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nx = 16;
  int Nz = 288;
  int Nvpar = 32; // 96 uniform
  int Nmu = 32;  // 192 uniform
  int poly_order = 1;
  double write_phase_freq = 1;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude = 1.e20;
  double ion_source_sigma = 0.5;
  double ion_source_temp = 5000. * eV;

  // POA parameters  
  int oap_first = 1; // Start with OAP phase
  double alpha_oap = 0.000005; // Alpha for OAP phase
  double alpha_fdp = 1.0; // Alpha for FDP phase
  double tau_oap = 3e-4; // OAP duration in microseconds
  double tau_fdp = 2e-9; // FDP duration in microseconds
  double tau_fdp_extra = 0.0; // Extra FDP duration in final cycle
  bool static_field_oap = true; // Dynamic field for OAP
  bool static_field_fdp = false; // Dynamic field for FDP
  bool update_species_oap = true; // Update species during OAP
  bool update_species_fdp = true; // Update species during FDP
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_oap = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE; // Default fdot multiplier for OAP
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_fdp = GKYL_GK_FDOT_MULTIPLIER_NONE; // Default fdot multiplier for FDP
  bool positivity_oap = false; // No positivity hack for OAP
  bool positivity_fdp = true; // Yes positivity hack for FDP
  int num_cycles = 100; // Number of OAP/FDP cycles
  int frames_per_phase = 5; // Frames per phase

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te0 = Te0,
    .n0 = n0,
    .B_p = B_p,
    .beta = beta,
    .tau = tau,
    .Ti0 = Ti0,
    .kperpRhos = kperpRhos,
    .alim = alim,
    .nuFrac = nuFrac,
    .logLambdaElc = logLambdaElc,
    .nuElc = nuElc,
    .elc_nuFrac = elc_nuFrac,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .kperp = kperp,
    .z_min = z_min,
    .z_max = z_max,
    .psi_min = psi_min,
    .psi_eval = psi_eval,
    .psi_max = psi_max,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .write_phase_freq = write_phase_freq,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .ion_source_amplitude = ion_source_amplitude,
    .ion_source_sigma = ion_source_sigma,
    .ion_source_temp = ion_source_temp,
    // POA parameters
    .oap_first = oap_first,
    .alpha_oap = alpha_oap,
    .alpha_fdp = alpha_fdp,
    .tau_oap = tau_oap,
    .tau_fdp = tau_fdp,
    .tau_fdp_extra = tau_fdp_extra,
    .static_field_oap = static_field_oap,
    .static_field_fdp = static_field_fdp,
    .update_species_oap = update_species_oap,
    .update_species_fdp = update_species_fdp,
    .fdot_mult_type_oap = fdot_mult_type_oap,
    .fdot_mult_type_fdp = fdot_mult_type_fdp,
    .positivity_oap = positivity_oap,
    .positivity_fdp = positivity_fdp,
    .num_cycles = num_cycles,
    .frames_per_phase = frames_per_phase,
  };
  
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

// Function to determine restart phase and cycle info
struct poa_restart_info
determine_restart_phase(struct gk_mirror_ctx* ctx, double t_curr, int frame_curr)
{
  struct poa_restart_info info = {0};
  
  // Use frame-based calculation for more reliable phase detection
  // Each cycle has 2 phases, each with frames_per_phase frames
  int frames_per_cycle = 2 * ctx->frames_per_phase;
  
  // Check if we're in the final extra FDP phase
  int regular_frames = ctx->num_cycles * frames_per_cycle;
  if (frame_curr >= regular_frames) {
    info.cycle = ctx->num_cycles - 1; // Last cycle
    info.is_oap_phase = false;
    info.is_final_fdp = true;
    
    double total_regular_time = ctx->num_cycles * (ctx->tau_oap + ctx->tau_fdp);
    info.phase_start_time = total_regular_time;
    info.phase_end_time = total_regular_time + ctx->tau_fdp_extra;
    info.phase_start_frame = regular_frames;
    int extra_frames = (int)((ctx->tau_fdp_extra / ctx->tau_fdp) * ctx->frames_per_phase);
    info.phase_end_frame = info.phase_start_frame + extra_frames;
    return info;
  }
  
  // We're in a regular cycle - use frame-based calculation
  info.cycle = frame_curr / frames_per_cycle;
  int frame_in_cycle = frame_curr % frames_per_cycle;
  info.is_final_fdp = false;
  
  // Calculate phase timing
  double total_phase_time = ctx->tau_oap + ctx->tau_fdp;
  
  // Determine which phase within the cycle based on frames
  if (ctx->oap_first) {
    if (frame_in_cycle < ctx->frames_per_phase) {
      // In OAP phase
      info.is_oap_phase = true;
      info.phase_start_time = info.cycle * total_phase_time;
      info.phase_end_time = info.phase_start_time + ctx->tau_oap;
      info.phase_start_frame = info.cycle * frames_per_cycle;
      info.phase_end_frame = info.phase_start_frame + ctx->frames_per_phase;
    } else {
      // In FDP phase
      info.is_oap_phase = false;
      info.phase_start_time = info.cycle * total_phase_time + ctx->tau_oap;
      info.phase_end_time = info.phase_start_time + ctx->tau_fdp;
      info.phase_start_frame = info.cycle * frames_per_cycle + ctx->frames_per_phase;
      info.phase_end_frame = info.phase_start_frame + ctx->frames_per_phase;
    }
  } else {
    if (frame_in_cycle < ctx->frames_per_phase) {
      // In FDP phase
      info.is_oap_phase = false;
      info.phase_start_time = info.cycle * total_phase_time;
      info.phase_end_time = info.phase_start_time + ctx->tau_fdp;
      info.phase_start_frame = info.cycle * frames_per_cycle;
      info.phase_end_frame = info.phase_start_frame + ctx->frames_per_phase;
    } else {
      // In OAP phase
      info.is_oap_phase = true;
      info.phase_start_time = info.cycle * total_phase_time + ctx->tau_fdp;
      info.phase_end_time = info.phase_start_time + ctx->tau_oap;
      info.phase_start_frame = info.cycle * frames_per_cycle + ctx->frames_per_phase;
      info.phase_end_frame = info.phase_start_frame + ctx->frames_per_phase;
    }
  }
  
  return info;
}

bool
run_simulation_phase(const char* phase_name, double phase_duration, 
                    struct phase_params params, int frames_per_phase,
                    double* tfinal, int* frame_end_phase, 
                    gkyl_gyrokinetic_app* app, struct gk_mirror_ctx* ctx,
                    struct gkyl_app_args* app_args, double* t_curr, int my_rank)
{
  // Calculate phase end time and frames
  *tfinal += phase_duration;
  int frame_start_phase = *frame_end_phase;
  *frame_end_phase += frames_per_phase;
  
  if (my_rank == 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "%s: tfinal=%.9e s | frame_start_phase=%d | frame_end_phase=%d\n", 
           phase_name, *tfinal, frame_start_phase, *frame_end_phase);
  }
  
  // Update simulation parameters for this phase
  gkyl_gyrokinetic_app_reset_fdot_mult(app, 0, params.alpha, params.fdot_mult_type);
  gkyl_gyrokinetic_app_reset_enforce_positivity(app, params.positivity);
  gkyl_gyrokinetic_app_reset_update_field(app, !params.static_field);
  gkyl_gyrokinetic_app_reset_update_species(app, 0, params.update_species);

  // For the triggers, use the current time as the starting point
  double phase_start_time = *t_curr;  // Start from current time
  
  // Create triggers for this phase
  struct gkyl_tm_trigger trig_write_conf = { 
    .dt = phase_duration/frames_per_phase, 
    .tcurr = phase_start_time, 
    .curr = frame_start_phase 
  };
  struct gkyl_tm_trigger trig_write_phase = { 
    .dt = phase_duration/(ctx->write_phase_freq * frames_per_phase), 
    .tcurr = phase_start_time, 
    .curr = frame_start_phase
  };
  struct gkyl_tm_trigger trig_calc_intdiag = { 
    .dt = phase_duration/GKYL_MAX2(frames_per_phase, (ctx->int_diag_calc_num/ctx->num_frames)*frames_per_phase),
    .tcurr = phase_start_time, 
    .curr = frame_start_phase 
  };
  
  // Run this phase
  double phase_t_end = *tfinal;
  double phase_t_curr = *t_curr;  // Start from current time
  
  // Initialize small time-step check for this phase
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;
  
  long step = 1, num_steps = app_args->num_steps;
  while ((phase_t_curr < phase_t_end) && (step <= num_steps)) {
    if (step == 1 || step % 1 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, phase_t_curr);

    double dt = phase_t_end - phase_t_curr;
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step == 1 || step % 1 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      return false;
    }
    phase_t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, phase_t_curr, phase_t_curr >= phase_t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, phase_t_curr, phase_t_curr >= phase_t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, phase_t_curr, true);
        write_data(&trig_write_conf, &trig_write_phase, app, phase_t_curr, true);
        return false;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  
  *t_curr = phase_t_curr;  // Update global time
  return true;
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mirror_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);


  int my_rank = 0;
  int comm_sz = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi){
    gkyl_comm_get_rank(comm, &my_rank);
    int comm_sz;
    gkyl_comm_get_size(comm, &comm_sz);
  }
#endif

  // Calculate total simulation time and frames from POA parameters
  double total_phase_time = ctx.tau_oap + ctx.tau_fdp;
  ctx.t_end = ctx.num_cycles * total_phase_time;
  ctx.num_frames = ctx.num_cycles * 2 * ctx.frames_per_phase;

  ctx.int_diag_calc_num = ctx.num_frames*100;

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .no_by = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = eval_density_ion,
      .ctx_density = &ctx,
      .upar = eval_upar_ion,
      .ctx_upar = &ctx,
      .temp = eval_temp_ion,
      .ctx_temp = &ctx,
    },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless_scale_factor = ctx.alpha_oap,
    .time_rate_multiplier = {
      .type = ctx.fdot_mult_type_oap,
      .cellwise_const = true,
    },

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Ti0,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .write_diagnostics = true,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      },
      .diagnostics = {
        .num_diag_moments = 6,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
      },
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP},
    },
  };
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .is_static = false,
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "./core/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_MIRROR_GRID_GEN_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_gk app_inp = {  // GK app
    .name = "gk_wham",
    .cdim = ctx.cdim ,  .vdim = ctx.vdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .enforce_positivity = true,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {ctx.psi_eval, 0.0},
      .mirror_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {ion},

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  struct poa_restart_info restart_info = {0};
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    // Determine which phase and cycle we're restarting from
    restart_info = determine_restart_phase(&ctx, t_curr, frame_curr);

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, "Restart info: cycle %d, %s phase", 
                             restart_info.cycle, restart_info.is_oap_phase ? "OAP" : "FDP");
    if (restart_info.is_final_fdp) {
      gkyl_gyrokinetic_app_cout(app, stdout, " (final FDP with extra frames)");
    }
    gkyl_gyrokinetic_app_cout(app, stdout, "\n");
    
    // Set up simulation parameters for the current phase
    if (restart_info.is_oap_phase) {
      gkyl_gyrokinetic_app_reset_fdot_mult(app, 0, ctx.alpha_oap, ctx.fdot_mult_type_oap);
      gkyl_gyrokinetic_app_reset_enforce_positivity(app, ctx.positivity_oap);
      gkyl_gyrokinetic_app_reset_update_field(app, !ctx.static_field_oap);
      gkyl_gyrokinetic_app_reset_update_species(app, 0, ctx.update_species_oap);
    } else {
      gkyl_gyrokinetic_app_reset_fdot_mult(app, 0, ctx.alpha_fdp, ctx.fdot_mult_type_fdp);
      gkyl_gyrokinetic_app_reset_enforce_positivity(app, ctx.positivity_fdp);
      gkyl_gyrokinetic_app_reset_update_field(app, !ctx.static_field_fdp);
      gkyl_gyrokinetic_app_reset_update_species(app, 0, ctx.update_species_fdp);
    }
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  
  // POA simulation structure similar to jobscript_POA_single
  // For restart, we need to set tfinal and frame_end_phase based on restart info
  double tfinal = app_args.is_restart ? t_curr : 0.0;
  int frame_end_phase = app_args.is_restart ? restart_info.phase_start_frame : 0;

  // Write out ICs
  struct gkyl_tm_trigger trig_write_conf_ic = { .dt = 1.0, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase_ic = { .dt = 1.0, .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag_ic = { .dt = 1.0, .tcurr = t_curr, .curr = frame_curr };
  calc_integrated_diagnostics(&trig_calc_intdiag_ic, app, t_curr, false);
  write_data(&trig_write_conf_ic, &trig_write_phase_ic, app, t_curr, false);

  // POA Loop Structure (similar to jobscript_POA_single)
  bool simulation_success = true;
  
  // Determine starting point for restart or normal run
  int start_cycle = app_args.is_restart ? restart_info.cycle : 0;
  bool start_in_oap = app_args.is_restart ? restart_info.is_oap_phase : ctx.oap_first;
  bool skip_to_final_fdp = app_args.is_restart && restart_info.is_final_fdp;
  
  if (skip_to_final_fdp) {
    // We're restarting in the final FDP phase
    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting in final FDP phase with extra frames\n");
    int extra_frames = (int)((ctx.tau_fdp_extra / ctx.tau_fdp) * ctx.frames_per_phase);
    frame_end_phase += extra_frames - ctx.frames_per_phase; // Adjust for extra frames
    struct phase_params fdp_extra_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                    ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
    
    // Calculate remaining time for this phase
    double remaining_time = restart_info.phase_end_time - t_curr;
    simulation_success = run_simulation_phase("FDP_EXTRA", remaining_time, fdp_extra_params, 
                                             extra_frames, &tfinal, &frame_end_phase, app, &ctx, 
                                             &app_args, &t_curr, my_rank);
  }
  else if (ctx.oap_first) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Running OAP first%s\n", app_args.is_restart ? " (restarting)" : "");
    
    // If restarting and we're in the middle of a phase, complete current phase first
    if (app_args.is_restart && t_curr > restart_info.phase_start_time) {
      const char* phase_name = restart_info.is_oap_phase ? "OAP" : "FDP";
      double remaining_time = restart_info.phase_end_time - t_curr;
      
      gkyl_gyrokinetic_app_cout(app, stdout, "Completing current %s phase (%.6e s remaining)\n", phase_name, remaining_time);
      
      if (restart_info.is_oap_phase) {
        struct phase_params oap_params = {ctx.alpha_oap, ctx.static_field_oap, 
                      ctx.fdot_mult_type_oap, ctx.positivity_oap, ctx.update_species_oap};
        simulation_success = run_simulation_phase("OAP", remaining_time, oap_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        if (!simulation_success) goto cleanup;
        
        // Now run the FDP phase for this cycle
        struct phase_params fdp_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                      ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
        simulation_success = run_simulation_phase("FDP", ctx.tau_fdp, fdp_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        start_cycle++; // Move to next cycle
      } else {
        struct phase_params fdp_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                      ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
        simulation_success = run_simulation_phase("FDP", remaining_time, fdp_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        start_cycle++; // Move to next cycle
      }
      if (!simulation_success) goto cleanup;
    }
    
    // Main loop for alternating OAP-FDP cycles
    for (int i = start_cycle; i < ctx.num_cycles && simulation_success; i++) {
      // OAP phase
      struct phase_params oap_params = {ctx.alpha_oap, ctx.static_field_oap, 
                    ctx.fdot_mult_type_oap, ctx.positivity_oap, ctx.update_species_oap};
      simulation_success = run_simulation_phase("OAP", ctx.tau_oap, oap_params, ctx.frames_per_phase,
                                               &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
      if (!simulation_success) break;
      
      // FDP phase
      struct phase_params fdp_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                    ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
      simulation_success = run_simulation_phase("FDP", ctx.tau_fdp, fdp_params, ctx.frames_per_phase,
                                               &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
    }
    
    // Final FDP phase with extra frames
    if (simulation_success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "Final FDP phase with extra frames\n");
      int extra_frames = (int)((ctx.tau_fdp_extra / ctx.tau_fdp) * ctx.frames_per_phase);
      frame_end_phase += extra_frames - ctx.frames_per_phase; // Adjust for extra frames
      struct phase_params fdp_extra_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                      ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
      simulation_success = run_simulation_phase("FDP_EXTRA", ctx.tau_fdp_extra, fdp_extra_params, 
                                               extra_frames, &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
    }
    
  } else {
    gkyl_gyrokinetic_app_cout(app, stdout, "Running FDP first%s\n", app_args.is_restart ? " (restarting)" : "");
    
    // If restarting and we're in the middle of a phase, complete current phase first
    if (app_args.is_restart && t_curr > restart_info.phase_start_time) {
      const char* phase_name = restart_info.is_oap_phase ? "OAP" : "FDP";
      double remaining_time = restart_info.phase_end_time - t_curr;
      
      gkyl_gyrokinetic_app_cout(app, stdout, "Completing current %s phase (%.6e s remaining)\n", phase_name, remaining_time);
      
      if (restart_info.is_oap_phase) {
        struct phase_params oap_params = {ctx.alpha_oap, ctx.static_field_oap, 
                      ctx.fdot_mult_type_oap, ctx.positivity_oap, ctx.update_species_oap};
        simulation_success = run_simulation_phase("OAP", remaining_time, oap_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        start_cycle++; // Move to next cycle
      } else {
        struct phase_params fdp_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                      ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
        simulation_success = run_simulation_phase("FDP", remaining_time, fdp_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        if (!simulation_success) goto cleanup;
        
        // Now run the OAP phase for this cycle
        struct phase_params oap_params = {ctx.alpha_oap, ctx.static_field_oap, 
                      ctx.fdot_mult_type_oap, ctx.positivity_oap, ctx.update_species_oap};
        simulation_success = run_simulation_phase("OAP", ctx.tau_oap, oap_params, ctx.frames_per_phase,
                                                 &tfinal, &frame_end_phase, app, &ctx, 
                                                 &app_args, &t_curr, my_rank);
        start_cycle++; // Move to next cycle
      }
      if (!simulation_success) goto cleanup;
    }
    
    // Main loop for alternating FDP-OAP cycles
    for (int i = start_cycle; i < ctx.num_cycles && simulation_success; i++) {
      // FDP phase
      struct phase_params fdp_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                    ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
      simulation_success = run_simulation_phase("FDP", ctx.tau_fdp, fdp_params, ctx.frames_per_phase,
                                               &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
      if (!simulation_success) break;
      
      // OAP phase
      struct phase_params oap_params = {ctx.alpha_oap, ctx.static_field_oap, 
                    ctx.fdot_mult_type_oap, ctx.positivity_oap, ctx.update_species_oap};
      simulation_success = run_simulation_phase("OAP", ctx.tau_oap, oap_params, ctx.frames_per_phase,
                                               &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
    }
    
    // Final FDP phase with extra frames
    if (simulation_success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "Final FDP phase with extra frames\n");
      int extra_frames = (int)((ctx.tau_fdp_extra / ctx.tau_fdp) * ctx.frames_per_phase);
      frame_end_phase += extra_frames - ctx.frames_per_phase; // Adjust for extra frames
      struct phase_params fdp_extra_params = {ctx.alpha_fdp, ctx.static_field_fdp, 
                      ctx.fdot_mult_type_fdp, ctx.positivity_fdp, ctx.update_species_fdp};
      simulation_success = run_simulation_phase("FDP_EXTRA", ctx.tau_fdp_extra, fdp_extra_params, 
                                               extra_frames, &tfinal, &frame_end_phase, app, &ctx, 
                                               &app_args, &t_curr, my_rank);
    }
  }

cleanup:

  gkyl_gyrokinetic_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Simulation completed at t = %g s\n", t_curr);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
