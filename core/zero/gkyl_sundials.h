#pragma once

/**
 *
 * SUNDIALS NVECTOR wrapper of Gkeyll data structure and methods.
 *
 */

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_comm.h>

// Time stepping method options.
enum gkyl_sundials_rk_method {
  GKYL_SUNDIALS_METHOD_NONE = 0,
  GKYL_RK_METHOD_SSP_3_3, // Gkeyll's native 3rd order 3-stage SSP RK method (for testing only).
  GKYL_SUNDIALS_LSRK_METHOD_RKC_2, // 2nd order Runge-Kutta-Chebyshev (RKC).
  GKYL_SUNDIALS_LSRK_METHOD_RKL_2, // 2nd order Runge-Kutta-Legendre (RKL).
  GKYL_SUNDIALS_LSRK_METHOD_SSP_S_2, // Optimal 2nd order s-stage SSP RK method.
  GKYL_SUNDIALS_LSRK_METHOD_SSP_S_3, // Optimal 3rd order s-stage SSP RK method.
  GKYL_SUNDIALS_LSRK_METHOD_SSP_10_4, // Optimal 4th order 10-stage SSP RK method.
};

// Context for functions that are app-specific and/or
// need to handle the app pointer as a void *.
struct gkyl_sundials_app_ctx {
  void *app_ptr; // Gkeyll app.
  // Function that computes df/dt.
  double (*dfdt_func)(void *app, double t_curr,
    const struct gkyl_array *fin[], struct gkyl_array *fout[], struct gkyl_array **bflux_out[],
    const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[]);
  // Function that reduces a local dt across MPI processes.
  double (*reduce_dt_func)(void *app, double t_curr, double dt_local);
  // Function that computes weight for the error norm.
  int (*error_wgt_func)(void *app, const struct gkyl_array *xarr,
    struct gkyl_array *wgt, struct gkyl_range *local_range);
  int num_species; // Number of species.
  int num_neut_species; // Number of neutral species.
  // Objects below are private.
  double dt_local; // CFL constrained time step in local MPI process.
  double dt_global; // Reduction of dt_local over all MPI processes.
};

// Sundials inputs specified in input file.
struct gkyl_sundials_stepper_inp {
  double rel_tol; // Relative tolerance.
  double abs_tol; // Absolute tolerance.
  long max_steps; // Maximum number of steps (Default: 1e5).
  int num_SSP_stages; // Number of stages in SSP RK method.
  enum gkyl_sundials_rk_method rk_method; // Time stepping method (Default: GKYL_RK_METHOD_SSP_3_3).
  struct gkyl_sundials_nvec *gsnv; // Input NVECTOR.
  double t_curr; // Current simulation time.
  struct gkyl_sundials_app_ctx *app_ctx; // Context with app-specific data and functions.
};

//
// SUNDIALS module methods.
//

typedef struct gkyl_sundials gkyl_sundials; // Sundials object.

/**
 * Create a new SUNDIALS object (e.g. a SUNDIALS context).
 *
 * @param use_gpu Whether data wrapped by SUNDIALS NVECTORS lives on the GPU.
 * @return A new SUNDIALS object.
 */
struct gkyl_sundials* gkyl_sundials_new(bool use_gpu);

/**
 * Initialize the SSP RK time stepper.
 *
 * @param gksun SUNDIALS object.
 * @param inp SUNDIALS inputs for time stepper.
 */
void gkyl_sundials_stepper_init_ssp_rk(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp);

/**
 * Reset the ARKODE time stepper given a new time and state vector.
 *
 * @param gksun SUNDIALS object.
 * @param time New simulation time.
 * @param gsmanynv New state vectors (as ManyNVector).
 * @param gsmanynv_buff Buffer array, as big as state vector array.
 */
void gkyl_sundials_arkode_reset(struct gkyl_sundials *gksun, double time,
  struct gkyl_sundials_nvec *gsmanynv, struct gkyl_sundials_nvec *gsmanynv_buff);

/**
 * Evolve the solution contained in a given Nvector from
 * the current to the new time.
 *
 * @param gksun SUNDIALS object.
 * @param t_new New time to evolve to.
 * @param gsnv Gkeyll wrapper of SUNDIALS NVECTOR.
 * @param t_curr Current time.
 */
int gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double *t_curr);

/**
 * Free resources associates with a SUNDIALS object.
 *
 * @param gksun SUNDIALS object to be freed.
 */
void gkyl_sundials_release(struct gkyl_sundials *gksun);

//
// SUNDIALS NVECTOR methods.

typedef struct gkyl_sundials_nvec gkyl_sundials_nvec; // Sundials nvector.

/**
 * Create a new NVECTOR wrapping for a given gkyl_array.
 *
 * @param gksun SUNDIALS object.
 * @param arr Gkeyll array to wrap.
 * @param comm Gkeyll communicator object.
 * @param local_range Local range to loop over arr.
 * @return A new NVECTOR.
 */
struct gkyl_sundials_nvec* gkyl_sundials_nvec_new(struct gkyl_sundials *gksun,
  struct gkyl_array *arr, struct gkyl_comm *comm, struct gkyl_range *local_range);

/**
 * Fetch the Gkeyll array wrapped by an NVECTOR.
 *
 * @param nvin Input NVECTOR.
 * @return Gkeyll array wrapped by NVECTOR.
 */
struct gkyl_array* gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *nvin);

/**
 * Destroy SUNDIALS NVECTOR.
 *
 * @param nvin NVECTOR to be destroyed.
 */
void gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *nvin);

/**
 * Create a new ManyNVECTOR wrapping of all our Nvectors (one for each
 * gkyl_array that is evolved via time integration of a PDE).
 *
 * @param gksun SUNDIALS object.
 * @param num_nvector Number of NVectors to wrap.
 * @param gsnv_arr Array of Gkeyll nvector objects.
 * @return A new NVECTOR, which holds the ManyVector.
 */
struct gkyl_sundials_nvec* gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun,
  int num_nvector, struct gkyl_sundials_nvec* gsnv_arr[]);

/**
 * Destroy SUNDIALS ManyNVECTOR.
 *
 * @param gsmanynv ManyNVECTOR to be destroyed.
 */
void gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv);

