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
enum gkyl_sundials_lsrk_method {
  GKYL_SUNDIALS_LSRK_METHOD_NONE = 0,
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
  int (*dfdt_func)(void *app, double t_curr,
    const struct gkyl_array *fin[], struct gkyl_array *fout[], struct gkyl_array **bflux_out[],
    const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[]);
  // Function that computes weight for the error norm.
  int (*error_wgt_func)(void *app, const struct gkyl_array *xarr,
    struct gkyl_array *wgt, struct gkyl_range *local_range);
};

// Sundials inputs specified in input file.
struct gkyl_sundials_stepper_inp {
  double rel_tol; // Relative tolerance.
  double abs_tol; // Absolute tolerance.
  long max_steps; // Maximum number of steps.
  int num_SSP_stages; // Number of stages in SSP RK method.
  struct gkyl_sundials_nvec *gsnv; // Input NVECTOR.
  double t_curr; // Current simulation time.
  enum gkyl_sundials_lsrk_method method; // Time stepping method.
  struct gkyl_sundials_app_ctx *app_ctx; // Context with app-specific data and functions.
};

//
// SUNDIALS module methods.
//

typedef struct gkyl_sundials gkyl_sundials; // Sundials object.

/**
 * Create a new SUNDIALS object (e.g. a SUNDIALS context).
 *
 * @param ncomp Number of components per cell in the state vector to be wrapped
 *              by SUNDIALS NVECTORS. Used to allocate memory for reductions.
 * @param use_gpu Whether data wrapped by SUNDIALS NVECTORS lives on the GPU.
 * @return A new SUNDIALS object.
 */
struct gkyl_sundials* gkyl_sundials_new(int ncomp, bool use_gpu);

/**
 * Initialize the SSP RK time stepper.
 *
 * @param gksun SUNDIALS object.
 * @param inp SUNDIALS inputs for time stepper.
 */
void gkyl_sundials_stepper_init_ssp_rk(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp);

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
  struct gkyl_sundials_nvec *gsnv, double t_curr);

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

