#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_multib.h>

struct gkyl_gyrokinetic_timings {
  double t_end; // End time for the simulation
  int num_frames; // Number of output frames. Output every (t_end/num_frames) time units.
  int write_phase_freq; // Frequency (in multiples of num_frames) of writing phase-space data. e.g. 0.2 means write every 5 frames.
  int int_diag_calc_num; // Number of integrated diagnostic calculations to do. 100*num_frames means 100 calculations per frame.
  double dt_failure_tol; // Tolerance for small time-step failures. Typical value is 1e-4
  int num_failures_max; // Maximum number of consecutive small time-step failures before aborting simulation. Typical value is 20.
  // Provided by app_args
  bool is_restart; // Is this a restarted simulation?
  int restart_frame; // Frame to restart from.
  int num_steps; // Maximum number of time-steps to take.
};

struct gkyl_gyrokinetic_run_verbosity {
  bool enabled; // Is verbosity enabled? Prints information every time step. Defaults false
  double frequency; // Print information with given frequency. Defaults to 1.0
  bool est_completion_time; // Estimate completion time based on current progress. Defaults false
};

enum gkyl_gyrokinetic_run_app_type {
    GKYL_GK_SINGLEB, // Single-block simulation. Default
    GKYL_GK_MULTIB, // Multi-block simulation
};

struct gkyl_gyrokinetic_run_inp {
  enum gkyl_gyrokinetic_run_app_type app_type; // Type of gyrokinetic application to run.
  union{
    struct gkyl_gk app_inp; // Single-block application input struct.
    struct gkyl_gyrokinetic_multib multib_app_inp; // Multi-block application input.
  };
  struct gkyl_gyrokinetic_timings timing; // Timing parameters for the simulation.
  struct gkyl_gyrokinetic_run_verbosity print_verbosity; // Verbosity settings for the simulation.
};

/**
 * Run the gyrokinetic simulation.
 * 
 * @param inp Input parameters for the simulation.
 */
void
gkyl_gyrokinetic_run_simulation(struct gkyl_gyrokinetic_run_inp* inp);