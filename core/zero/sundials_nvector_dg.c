/** 
 *
 * NVECTOR wrapper of the Gkeyll data structure.
 *
 */

#include <gkyl_sundials_nvector_dg.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <string.h> // for memcpy.

/**
 * Create a new Gkeyll array.
 *
 * @param use_gpu Whether to allocate it on the GPU.
 * @param nc Number of components per cell.
 * @param size Number of cells.
 * @return A new gkyl_array pointer.
 */
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

N_Vector
gkyl_sundials_nvec_new_empty(SUNContext sunctx)
{
  N_Vector nvec;
  N_VectorContent_Gkeyll content;

  // Create an empty vector object.
  nvec = NULL;
  nvec = N_VNewEmpty(sunctx);
  if (nvec == NULL)
    return (NULL);

  // Attach operations.
  // Constructors, destructors, and utility operations.
  nvec->ops->nvclone      = gkyl_sundials_nvec_clone;
  nvec->ops->nvcloneempty = gkyl_sundials_nvec_clone_empty;
  nvec->ops->nvdestroy    = gkyl_sundials_nvec_destroy;

//  // Data operations.
//  nvec->ops->nvlinearsum = N_VLinearSum_Gkeyll;
//  nvec->ops->nvconst     = N_VConst_Gkeyll;
//  nvec->ops->nvscale     = N_VScale_Gkeyll;
//  nvec->ops->nvwrmsnorm  = N_VWrmsNorm_abs_comp_Gkeyll;
//  nvec->ops->nvdotprod   = N_VDotProd_Gkeyll;
//  nvec->ops->nvspace     = N_VSpace_Gkeyll;
//  nvec->ops->nvdiv       = N_VDiv_Gkeyll;
//  nvec->ops->nvabs       = N_VAbs_Gkeyll;
//  nvec->ops->nvinv       = N_VInv_Gkeyll;
//  nvec->ops->nvmaxnorm   = N_VMaxnorm_Gkeyll;
//  nvec->ops->nvaddconst  = N_VAddconst_Gkeyll;

  // Create content.
  content = NULL;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == NULL) {
    N_VDestroy(nvec);
    return (NULL);
  }

  // Attach content.
  nvec->content = content;

  // Attach gkyl_array.
  content->own_vector  = SUNFALSE;
  content->use_gpu     = SUNFALSE;
  content->dataptr     = NULL;
  content->comm        = NULL;
  content->local_range = NULL;

  return nvec;
}

N_Vector
gkyl_sundials_nvec_make(struct gkyl_array* garr, sunbooleantype use_gpu,
  struct gkyl_comm* comm, struct gkyl_range* local_range, SUNContext sunctx)
{
  N_Vector nvout;
  nvout = NULL;
  nvout = gkyl_sundials_nvec_new_empty(sunctx);
  if (nvout == NULL)
    return (NULL);

  NV_CONTENT_GKZ(nvout)->own_vector  = SUNFALSE;
  NV_CONTENT_GKZ(nvout)->use_gpu     = use_gpu;
  NV_CONTENT_GKZ(nvout)->comm        = comm;
  NV_CONTENT_GKZ(nvout)->local_range = local_range;
  NV_CONTENT_GKZ(nvout)->dataptr     = garr;
  return nvout;
}

struct gkyl_array*
gkyl_sundials_nvec_get_vector(N_Vector nvin)
{
  return NV_CONTENT_GKZ(nvin)->dataptr;
}

N_Vector
gkyl_sundials_nvec_clone_empty(N_Vector nvin)
{
  N_Vector nvout;
  N_VectorContent_Gkeyll content;
  if (nvin == NULL)
    return (NULL);

  // Create vector.
  nvout = NULL;
  nvout = N_VNewEmpty(nvin->sunctx);
  if (nvout == NULL)
    return (NULL);

  // Attach operations.
  if (N_VCopyOps(nvin, nvout)) {
    N_VDestroy(nvout);
    return (NULL);
  }

  // Create content.
  content = NULL;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == NULL) {
    N_VDestroy(nvout);
    return (NULL);
  }

  // Attach content.
  nvout->content = content;

  // Initialize content.
  content->own_vector = SUNFALSE;
  // gkyl_sundials_nvec_clone_empty must be called only in
  // gkyl_sundials_nvec_clone to ensure use_gpu has the correct flag.
  // Otherwise, use_gpu flag will be false even if it must be true.

  // TO DO: Check to verify if this function is called separately
  content->use_gpu     = SUNFALSE;
  content->dataptr     = NULL;
  content->comm        = NULL;
  content->local_range = NULL;

  return nvout;
}

N_Vector
gkyl_sundials_nvec_clone(N_Vector nvin)
{
  N_Vector nvout;
  struct gkyl_array* nvout_dptr;
  struct gkyl_array* nvin_dptr = NV_CONTENT_GKZ(nvin)->dataptr;

  nvout = NULL;
  nvout = gkyl_sundials_nvec_clone_empty(nvin);
  if (nvout == NULL)
    return (NULL);

  nvout_dptr = mkarr(NV_CONTENT_GKZ(nvin)->use_gpu, nvin_dptr->ncomp, nvin_dptr->size);

  NV_CONTENT_GKZ(nvout)->dataptr     = nvout_dptr;
  NV_CONTENT_GKZ(nvout)->use_gpu     = NV_CONTENT_GKZ(nvin)->use_gpu;
  NV_CONTENT_GKZ(nvout)->comm        = NV_CONTENT_GKZ(nvin)->comm;
  NV_CONTENT_GKZ(nvout)->local_range = NV_CONTENT_GKZ(nvin)->local_range;
  NV_CONTENT_GKZ(nvout)->own_vector  = SUNTRUE;
  return nvout;
}

void
gkyl_sundials_nvec_destroy(N_Vector nvin)
{
  if (nvin == NULL)
    return;

  // Free content.
  if (nvin->content != NULL) {
    // Free the gkeyll parvector if it's owned by the vector wrapper.
    if (NV_CONTENT_GKZ(nvin)->own_vector && NV_CONTENT_GKZ(nvin)->dataptr != NULL) {
      gkyl_array_release(NV_CONTENT_GKZ(nvin)->dataptr);
      NV_CONTENT_GKZ(nvin)->dataptr = NULL;
    }
    free(nvin->content);
    nvin->content = NULL;
  }

  // Free ops and vector.
  if (nvin->ops != NULL) {
    free(nvin->ops);
    nvin->ops = NULL;
  }
  free(nvin);
  nvin = NULL;
}
