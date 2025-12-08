#pragma once

#ifdef GKYL_HAVE_SUNDIALS
/**
 *
 * SUNDIALS NVECTOR wrapper of Gkeyll data structure and methods.
 *
 */

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_comm.h>

typedef struct gkyl_sundials_nvec gkyl_sundials_nvec;

/**
 * Create a new NVECTOR wrapping the provided gkyl_array.
 *
 * @param arr Gkeyll array to wrap.
 * @param use_gpu Whether array data lives on the GPU.
 * @param comm Gkeyll communicatr object.
 * @param local_range Local range to loop over arr.
 * @return A new NVECTOR.
 */
struct gkyl_sundials_nvec* gkyl_sundials_nvec_new(struct gkyl_array* arr, bool use_gpu,
  struct gkyl_comm* comm, struct gkyl_range* local_range);

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

#endif
