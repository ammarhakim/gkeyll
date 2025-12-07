#pragma once

/**
 *
 * SUNDIALS NVECTOR wrapper of Gkeyll data structure and methods.
 *
 */

#ifndef _NVECTOR_GKYLZERO_H
#define _NVECTOR_GKYLZERO_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h> /* definition of type sunrealtype          */

/* gkylzero header files -- ADD MORE AS NECESSARY */
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#include <mpi.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Gkeyll implementation of N_Vector
 * -----------------------------------------------------------------*/

struct _N_VectorContent_Gkeyll
{
  sunbooleantype own_vector; /* ownership Gkeyll vector */
  sunbooleantype use_gpu;    /* where to reside */
  // Gkeyll pointers.
  struct gkyl_comm* comm;         // Communicator.
  struct gkyl_range* local_range; // Local range.
  struct gkyl_array* dataptr;     // Data of this array/vector.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

/* -----------------------------------------------------------------
 * Functions exported by nvector_gkylzero
 * -----------------------------------------------------------------*/
struct gkyl_array* mkarr(bool on_gpu, long nc, long size);

N_Vector N_VNewEmpty_Gkeyll(SUNContext sunctx);
N_Vector N_VMake_Gkeyll(struct gkyl_array* x, sunbooleantype use_gpu,
                          struct gkyl_comm* comm,
                          struct gkyl_range* local_range, SUNContext sunctx);
struct gkyl_array* N_VGetVector_Gkeyll(N_Vector v);

N_Vector N_VCloneEmpty_Gkeyll(N_Vector w);
N_Vector N_VClone_Gkeyll(N_Vector w);
void N_VDestroy_Gkeyll(N_Vector v);

/* vector operations -- DELETE ALL THAT WILL BE UNUSED BY LSRKSTEP */
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

#ifdef __cplusplus
}
#endif

#endif
