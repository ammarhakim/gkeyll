#pragma once

// Private header file for SUNDIALS nvector wrapper.

#include <gkyl_sundials.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h> // Definition of type sunrealtype.

// Give access to the Gkeyll vector from within the NVECTOR.
#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

struct gkyl_sundials_op_mem
{
  int ncomp; // Number of components per cell of vectors wrapped.
  double *red_local, *red_global; // Result of local and global reductions.
  double *red_global_ho; // Result of global reduction on the host.
  struct gkyl_array *arr_buff; // Auxiliary array buffer.
};

// Gkeyll implementation of N_Vector.
struct _N_VectorContent_Gkeyll
{
  struct gkyl_array *arr; // Array wrapped by this vector.
  sunbooleantype own_vector; // Ownership of the gkyl_array.
  sunbooleantype use_gpu; // Whether data is on GPU.
  struct gkyl_comm *comm; // Communicator.
  struct gkyl_range *local_range; // Local range.
  struct gkyl_sundials_op_mem *op_mem; // Memory for reductions.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

struct gkyl_sundials_nvec
{
  N_Vector nvec; // SUNDIALS NVector.
};

struct gkyl_sundials
{
  SUNContext sunctx; // Sundials context.
  struct gkyl_sundials_op_mem op_mem; // Memory for reductions.
  bool use_gpu;
};
