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

typedef struct gkyl_sundials gkyl_sundials; // Sundials object.
typedef struct gkyl_sundials_nvec gkyl_sundials_nvec; // Sundials nvector.

/**
 * Create a new SUNDIALS object (e.g. a SUNDIALS context).
 *
 * @param ncomp Number of components per cell in the state vector to be wrapped
 *              by SUNDIALS NVECTORS. Used to allocate memory for reductions.
 * @param local_ext_vol Volume of local extended range Gkeyll arrays are
 *                      allocated in. Used to allocate gkyl_array for dot
 *                      product (should be rewritten so this isn't needed).
 * @param use_gpu Whether data wrapped by SUNDIALS NVECTORS lives on the GPU.
 * @return A new SUNDIALS object.
 */
struct gkyl_sundials* gkyl_sundials_new(int ncomp, long local_ext_vol, bool use_gpu);

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
 * Free resources associates with a SUNDIALS object.
 *
 * @param gksun SUNDIALS object to be freed.
 */
void gkyl_sundials_release(struct gkyl_sundials *gksun);

#endif
