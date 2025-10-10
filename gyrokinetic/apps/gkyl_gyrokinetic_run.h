#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

struct gkyl_gyrokinetic_timings
{
  double t_end; // End time.
  int num_frames; // Number of output frames.
  int write_phase_freq; // Frequency (in multiples of num_frames) of writing phase-space data.
  int int_diag_calc_num; // Number of integrated diagnostic calculations to do.
  double dt_failure_tol; // Tolerance for small time-step failures.
  int num_failures_max; // Maximum number of consecutive small time-step failures before aborting simulation.
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Frame to restart from.
  int num_steps; // Maximum number of time-steps to take.
};

struct gkyl_gyrokinetic_run_inp
{
  struct gkyl_gk* app_inp; // Input parameters for the gyrokinetic application.
  struct gkyl_gyrokinetic_multib* multib_app_inp; // Input parameters for the multi-species gyrokinetic application.
  struct gkyl_gyrokinetic_timings timing; // Timing parameters for the simulation.
};

/**
 * Run the gyrokinetic simulation.
 * 
 * @param app_inp Input parameters for the gyrokinetic application.
 * @param timing Timing parameters for the simulation.
 * @param app_args Command-line arguments for the application.
 */
void gkyl_gyrokinetic_run_simulation(struct gkyl_gyrokinetic_run_inp inp);