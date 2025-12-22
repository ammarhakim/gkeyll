#pragma once

// Private header file for SUNDIALS nvector wrapper.

#include <gkyl_sundials.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_types.h> // Definition of type sunrealtype.
#include <arkode/arkode_lsrkstep.h> // Access to LSRKStep.
#include <arkode/arkode_erkstep.h> // Access to ERKStep.
#include <sundomeigest/sundomeigest_power.h> // Power iteration dominant eigen value estimator.

// Give access to the Gkeyll vector from within the NVECTOR.
#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

struct gkyl_sundials_reduce_mem
{
  double *red_local, *red_global; // Result of local and global reductions.
  double *red_global_ho; // Result of global reduction on the host.
};

// Gkeyll implementation of N_Vector.
struct _N_VectorContent_Gkeyll
{
  struct gkyl_array *arr; // Array wrapped by this vector.
  sunbooleantype own_vector; // Ownership of the gkyl_array.
  sunbooleantype use_gpu; // Whether data is on GPU.
  struct gkyl_comm *comm; // Communicator.
  struct gkyl_range *local_range; // Local range.
  struct gkyl_sundials_reduce_mem *red_mem; // Memory for reductions.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

struct gkyl_sundials_nvec
{
  N_Vector nvec; // SUNDIALS NVector.
  struct gkyl_sundials_reduce_mem red_mem; // Memory for reductions.
};

struct gkyl_sundials
{
  SUNContext sunctx; // Sundials context.
  enum gkyl_sundials_rk_method rk_method; // RK method.
  bool use_gpu; // Whether to run on GPU.
  void *arkode_mem; // Memory for ARKODE.
  struct gkyl_sundials_app_ctx *app_ctx; // App-specific context.
  SUNDomEigEstimator dom_eig_est; // Dominant eigenvalue estimator.
};

/**
 * Check function return value.
 *   opt == 0 means function allocates memory so check if returned NULL pointer.
 *   opt == 1 means function returns a flag so check if flag >= 0.
 */
static void
sundials_check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  if (opt == 0) {
    if (flagvalue == NULL) {
      // Check if function returned NULL pointer - no memory allocated.
      fprintf(stderr, "\nError: %s() failed - returned NULL pointer\n", funcname);
      assert(false);
    }
  }
  else if (opt == 1) {
    // Check if flag != 0.
    errflag = (int*)flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nError %s() failed with flag = %d\n", funcname, *errflag);
      assert(false);
    }
  }
}
