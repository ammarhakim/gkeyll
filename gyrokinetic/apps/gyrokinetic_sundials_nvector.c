/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds and Mustafa Aggul @ UMBC
 * -----------------------------------------------------------------
 * This is the implementation file for an NVECTOR wrapper of the
 * Gkeyll data structure.
 * -----------------------------------------------------------------*/

#include <gkyl_gyrokinetic_sundials_nvector.h>
#include <gkyl_alloc.h>
#include <string.h> // for memcpy.

/* -----------------------------------------------------------------
 * Simplifying macro: NV_CONTENT_GKZ
 *
 * This gives access to the Gkeyll vector from within the NVECTOR.
 * -----------------------------------------------------------------*/

/* -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------*/

struct gkyl_array* mkarr(bool on_gpu, long nc, long size)
{
  // Allocate array (filled with zeros).
  struct gkyl_array* a;
  if (on_gpu) a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

/* Construct a new Gkeyll vector without the underlying Gkeyll vector */
N_Vector N_VNewEmpty_Gkeyll(SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_Gkeyll content;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations -- DELETE ANY THAT WERE DELETED FROM nvector_gkylzero.h */

//  /* constructors, destructors, and utility operations */
//  v->ops->nvclone      = N_VClone_Gkeyll;
//  v->ops->nvcloneempty = N_VCloneEmpty_Gkeyll;
//  v->ops->nvdestroy    = N_VDestroy_Gkeyll;
//
//  /* vector operations */
//  v->ops->nvlinearsum = N_VLinearSum_Gkeyll;
//  v->ops->nvconst     = N_VConst_Gkeyll;
//  v->ops->nvscale     = N_VScale_Gkeyll;
//  v->ops->nvwrmsnorm  = N_VWrmsNorm_abs_comp_Gkeyll;
//  v->ops->nvdotprod   = N_VDotProd_Gkeyll;
//  v->ops->nvspace     = N_VSpace_Gkeyll;
//  v->ops->nvdiv       = N_VDiv_Gkeyll;
//  v->ops->nvabs       = N_VAbs_Gkeyll;
//  v->ops->nvinv       = N_VInv_Gkeyll;
//  v->ops->nvmaxnorm   = N_VMaxnorm_Gkeyll;
//  v->ops->nvaddconst  = N_VAddconst_Gkeyll;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Gkeyll)malloc(sizeof *content);
  if (content == NULL)
  {
    N_VDestroy(v);
    return (NULL);
  }

  /* Attach content */
  v->content = content;

  /* Attach gkyl_array */
  content->own_vector  = SUNFALSE;
  content->use_gpu     = SUNFALSE;
  content->dataptr     = NULL;
  content->comm        = NULL;
  content->local_range = NULL;

  return (v);
}
