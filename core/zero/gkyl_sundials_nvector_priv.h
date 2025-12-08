#pragma once

// Private header file for SUNDIALS nvector wrapper.

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h> // Definition of type sunrealtype.

// Give access to the Gkeyll vector from within the NVECTOR.
#define NV_CONTENT_GKZ(v) ((N_VectorContent_Gkeyll)(v->content))

// Gkeyll implementation of N_Vector.
struct _N_VectorContent_Gkeyll
{
  sunbooleantype own_vector; // Ownership Gkeyll vector.
  sunbooleantype use_gpu;    // Whether data is on GPU.
  // Gkeyll pointers.
  struct gkyl_comm *comm;         // Communicator.
  struct gkyl_range *local_range; // Local range.
  struct gkyl_array *arr;         // Array of this vector.
};

typedef struct _N_VectorContent_Gkeyll* N_VectorContent_Gkeyll;

struct gkyl_sundials_nvec
{
  SUNContext sunctx;
  N_Vector nvec;
};
