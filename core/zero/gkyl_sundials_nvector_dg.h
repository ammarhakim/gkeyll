#pragma once

#ifdef GKYL_HAVE_SUNDIALS
/**
 *
 * SUNDIALS NVECTOR wrapper of Gkeyll data structure and methods.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h> /* definition of type sunrealtype          */

// Gkeyll header files.
#include <gkyl_array.h>
#include <gkyl_range.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#include <mpi.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

// Give access to the Gkeyll vector from within the NVECTOR.
#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

// Gkeyll implementation of N_Vector.
struct _N_VectorContent_Gkeyll
{
  sunbooleantype own_vector; // Ownership Gkeyll vector.
  sunbooleantype use_gpu;    // Whether data is on GPU.
  // Gkeyll pointers.
  struct gkyl_comm* comm;         // Communicator.
  struct gkyl_range* local_range; // Local range.
  struct gkyl_array* dataptr;     // Data of this array/vector.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

/**
 * Construct a new SUNDIALS NVECTOR without the underlying Gkeyll array.
 *
 * @param sunctx SUNDIALS context object.
 * @return A new nvector.
 */
N_Vector gkyl_sundials_nvec_vnew_empty(SUNContext sunctx);

/**
 * Create a new NVECTOR wrapping the provided gkyl_array.
 *
 * @param garr Gkeyll array to wrap.
 * @param use_gpu Whether array data lives on the GPU.
 * @param comm Gkeyll communicatr object.
 * @param local_range Local range to loop over garr.
 * @param sunctx SUNDIALS context object.
 * @return A new NVECTOR.
 */
N_Vector gkyl_sundials_nvec_make(struct gkyl_array* garr, sunbooleantype use_gpu,
  struct gkyl_comm* comm, struct gkyl_range* local_range, SUNContext sunctx);

/**
 * Fetch the Gkeyll array wrapped by an NVECTOR.
 *
 * @param nvin Input NVECTOR.
 * @return Gkeyll array wrapped by NVECTOR.
 */
struct gkyl_array* gkyl_sundials_nvec_get_vector(N_Vector nvin);

/**
 * Make an empty copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
N_Vector gkyl_sundials_nvec_clone_empty(N_Vector nvin);

/**
 * Make an copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
N_Vector gkyl_sundials_nvec_clone(N_Vector nvin);

// Vector operations.
void N_VLinearSum_Gkeyll(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
  N_Vector z);
void N_VConst_Gkeyll(sunrealtype c, N_Vector z);
void N_VScale_Gkeyll(sunrealtype c, N_Vector x, N_Vector z);
sunrealtype N_VWrmsNorm_abs_comp_Gkeyll(N_Vector x, N_Vector w);
sunrealtype N_VWrmsNorm_cell_norm_Gkeyll(N_Vector x, N_Vector w);
sunrealtype N_VDotProd_Gkeyll(N_Vector x, N_Vector y);
void N_VSpace_Gkeyll(N_Vector v, sunindextype* x, sunindextype* y);
void N_VDiv_Gkeyll(N_Vector u, N_Vector v, N_Vector w);
void N_VAbs_Gkeyll(N_Vector u, N_Vector v);
void N_VInv_Gkeyll(N_Vector u, N_Vector v);
sunrealtype N_VMaxnorm_Gkeyll(N_Vector u);
void N_VAddconst_Gkeyll(N_Vector u, sunrealtype x, N_Vector v);

/**
 * Destroy SUNDIALS NVECTOR.
 *
 * @param nvec NVECTOR to be destroyed.
 */
void gkyl_sundials_nvec_destroy(N_Vector nvec);

#ifdef __cplusplus
}
#endif

#endif
