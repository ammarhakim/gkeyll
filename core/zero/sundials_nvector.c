/** 
 *
 * NVECTOR wrapper of the Gkeyll data structure.
 *
 */

#include <gkyl_sundials_nvector.h>
#include <gkyl_sundials_nvector_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <string.h> // for memcpy.

//// Vector operations.
//sunrealtype N_VWrmsNorm_cell_norm_Gkeyll(N_Vector x, N_Vector w);
//sunrealtype N_VDotProd_Gkeyll(N_Vector x, N_Vector y);
//void N_VSpace_Gkeyll(N_Vector v, sunindextype* x, sunindextype* y);
//void N_VDiv_Gkeyll(N_Vector u, N_Vector v, N_Vector w);
//void N_VAbs_Gkeyll(N_Vector u, N_Vector v);
//void N_VInv_Gkeyll(N_Vector u, N_Vector v);
//sunrealtype N_VMaxnorm_Gkeyll(N_Vector u);
//void N_VAddconst_Gkeyll(N_Vector u, sunrealtype x, N_Vector v);

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

/**
 * Make an empty copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
static N_Vector
snvec_clone_empty(N_Vector nvin)
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
  content->arr         = NULL;
  content->comm        = NULL;
  content->local_range = NULL;

  return nvout;
}

/**
 * Make an copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
static N_Vector
snvec_clone(N_Vector nvin)
{
  N_Vector nvout;

  struct gkyl_array* nvout_arr;
  struct gkyl_array* nvin_arr = NV_CONTENT_GKZ(nvin)->arr;

  nvout = NULL;
  nvout = snvec_clone_empty(nvin);
  if (nvout == NULL)
    return (NULL);

  nvout_arr = mkarr(NV_CONTENT_GKZ(nvin)->use_gpu, nvin_arr->ncomp, nvin_arr->size);

  NV_CONTENT_GKZ(nvout)->arr         = nvout_arr;
  NV_CONTENT_GKZ(nvout)->use_gpu     = NV_CONTENT_GKZ(nvin)->use_gpu;
  NV_CONTENT_GKZ(nvout)->comm        = NV_CONTENT_GKZ(nvin)->comm;
  NV_CONTENT_GKZ(nvout)->local_range = NV_CONTENT_GKZ(nvin)->local_range;
  NV_CONTENT_GKZ(nvout)->own_vector  = SUNTRUE;

  return nvout;
}

/**
 * Destroy SUNDIALS NVECTOR.
 *
 * @param nvin NVECTOR to be destroyed.
 */
static void
snvec_destroy(N_Vector nvin)
{
  if (nvin == 0)
    return;

  // Free content.
  if (nvin->content != 0) {
    // Free the gkeyll parvector if it's owned by the vector wrapper.
    if (NV_CONTENT_GKZ(nvin)->own_vector && NV_CONTENT_GKZ(nvin)->arr != 0) {
      gkyl_array_release(NV_CONTENT_GKZ(nvin)->arr);
      NV_CONTENT_GKZ(nvin)->arr = 0;
    }
    free(nvin->content);
    nvin->content = 0;
  }

  // Free ops and vector.
  if (nvin->ops != 0) {
    free(nvin->ops);
    nvin->ops = 0;
  }
  free(nvin);
  nvin = 0;
}

/**
 * Perform the linear sum z = a*x + b*y where z, x and y are SUNDIALS NVECTORS,
 * and a and b are scalars.
 *
 * @param a Scalar to multiply `x`.
 * @param x Nvector to multiply by `a` and add to `z`.
 * @param b Scalar to multiply `y`.
 * @param y Nvector to multiply by `b` and add to `z`.
 * @param z Output Nvector.
 */
void
snvec_linear_sum(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y, N_Vector z)
{
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *yarr = NV_CONTENT_GKZ(y)->arr;
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;

  gkyl_array_comp_op_range(zarr, GKYL_AXPBY, a, xarr, b, yarr, local_range);
}

/**
 * Set the whole SUNDIALS NVECTOR to a constant value.
 *
 * @param c Scalar value.
 * @param z Nvector to assign `c` value to.
 */
void
snvec_const(sunrealtype c, N_Vector z)
{
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(z)->local_range;

  gkyl_array_clear_range(zarr, c, local_range);
}

/**
 * Scale a given SUNDIALS NVECTOR x and assign it to another
 * NVECTOR z: z = c * x.
 *
 * @param c Scalar value.
 * @param x Nvector to multiply by c.
 * @param z Output Nvector.
 */
void
snvec_scale(sunrealtype c, N_Vector x, N_Vector z)
{
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;

  gkyl_array_set_range(zarr, c, xarr, local_range);
}

/**
 * Compute the WRMS norm of a Nvector x (typically the error in the solution
 * vector) given the weight w. The norm is
 *   out = sqrt( (1/N) * sum_i^N w_i |x_i|_C^2 ),
 * where
 *   |x_i|_C^2 = (1/N_b) * sum_{k=1}^{N_b} x_i^{(k)},
 *   N = number of cells,
 *   N_b = number of basis functions per cell,
 *   x_i^{(k)} = k-th coefficient in the i-th cell.
 * The weight w is computed externally, e.g.
 *   w_i = 1/[eps_rel*|y_i|_C+eps_abs]
 * where eps_rel and eps_abs are the relative and absolute tolerance,
 * respectively, and y is the solution vector at some stage.
 *
 * @param c Scalar value.
 * @param x Nvector to multiply by c.
 * @param z Output Nvector.
 */
sunrealtype
snvec_wrms_norm(N_Vector nvx, N_Vector nvwgt)
{
  struct gkyl_array* nvx_arr     = NV_CONTENT_GKZ(nvx)->arr;
  struct gkyl_array* nvwgt_arr   = NV_CONTENT_GKZ(nvwgt)->arr;
  struct gkyl_comm* comm         = NV_CONTENT_GKZ(nvwgt)->comm;
  struct gkyl_range* local_range = NV_CONTENT_GKZ(nvx)->local_range;
  bool use_gpu                   = NV_CONTENT_GKZ(nvx)->use_gpu;

  // TODO: change code so these allocations only happen once.
  int ncomp      = nvx_arr->ncomp;
  double* red_ho = gkyl_malloc(ncomp * sizeof(double));
  double *red_local, *red_global;
  if (use_gpu) {
    red_local  = gkyl_cu_malloc(ncomp * sizeof(double));
    red_global = gkyl_cu_malloc(ncomp * sizeof(double));
  }
  else {
    red_local  = gkyl_malloc(ncomp * sizeof(double));
    red_global = gkyl_malloc(ncomp * sizeof(double));
  }

  // Reduce over cells.
  gkyl_array_reduce_weighted_range(red_local, nvx_arr, nvwgt_arr, GKYL_RMS, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, ncomp, red_local, red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_ho, red_global, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_ho, red_global, ncomp * sizeof(double));

  // Reduce over components.
  //  sunrealtype red_out = 0.0;
  //  for (sunindextype i = 0; i < ncomp; ++i) red_out += red_ho[i];
  // Use the 0th component because each component should have the same result.
  sunrealtype red_out = red_ho[0];

  red_out = SUNRsqrt(red_out / local_range->volume);

  gkyl_free(red_ho);
  if (use_gpu) {
    gkyl_cu_free(red_local);
    gkyl_cu_free(red_global);
  }
  else {
    gkyl_free(red_local);
    gkyl_free(red_global);
  }

  return red_out;
}

/**
 * Construct a new SUNDIALS NVECTOR without the underlying Gkeyll array.
 *
 * @param sunctx SUNDIALS context object.
 * @return A new nvector.
 */
static N_Vector
snvec_new_empty(SUNContext sunctx)
{
  N_Vector nvec;
  N_VectorContent_Gkeyll content;

  // Create an empty vector object.
  nvec = 0;
  nvec = N_VNewEmpty(sunctx);
  if (nvec == 0)
    return 0;

  // Attach operations.
  // Constructors, destructors, and utility operations.
  nvec->ops->nvclone      = snvec_clone;
  nvec->ops->nvcloneempty = snvec_clone_empty;
  nvec->ops->nvdestroy    = snvec_destroy;

  // Data operations.
  nvec->ops->nvlinearsum = snvec_linear_sum;
  nvec->ops->nvconst     = snvec_const;
  nvec->ops->nvscale     = snvec_scale;
  nvec->ops->nvwrmsnorm  = snvec_wrms_norm;
//  nvec->ops->nvdotprod   = N_VDotProd_Gkeyll;
//  nvec->ops->nvspace     = N_VSpace_Gkeyll;
//  nvec->ops->nvdiv       = N_VDiv_Gkeyll;
//  nvec->ops->nvabs       = N_VAbs_Gkeyll;
//  nvec->ops->nvinv       = N_VInv_Gkeyll;
//  nvec->ops->nvmaxnorm   = N_VMaxnorm_Gkeyll;
//  nvec->ops->nvaddconst  = N_VAddconst_Gkeyll;

  // Create content.
  content = 0;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == 0) {
    N_VDestroy(nvec);
    return 0;
  }

  // Attach content.
  nvec->content = content;

  // Attach gkyl_array.
  content->own_vector  = SUNFALSE;
  content->use_gpu     = SUNFALSE;
  content->arr         = 0;
  content->comm        = 0;
  content->local_range = 0;

  return nvec;
}

struct gkyl_sundials *
gkyl_sundials_new()
{
  struct gkyl_sundials *gksun = gkyl_malloc(sizeof(*gksun));

  // Create the SUNDIALS context object.
  SUNContext_Create(SUN_COMM_NULL, &gksun->sunctx);

  return gksun;
}

struct gkyl_sundials_nvec*
gkyl_sundials_nvec_new(struct gkyl_sundials *gksun, struct gkyl_array *arr,
  struct gkyl_comm *comm, struct gkyl_range *local_range)
{
  struct gkyl_sundials_nvec *gsnv = gkyl_malloc(sizeof(*gsnv));

  N_Vector nvout = gsnv->nvec;
  nvout = 0;
  nvout = snvec_new_empty(gksun->sunctx);
  if (nvout == 0)
    return 0;

  NV_CONTENT_GKZ(nvout)->own_vector  = SUNFALSE;
  NV_CONTENT_GKZ(nvout)->use_gpu     = gkyl_array_is_cu_dev(arr);
  NV_CONTENT_GKZ(nvout)->comm        = comm;
  NV_CONTENT_GKZ(nvout)->local_range = local_range;
  NV_CONTENT_GKZ(nvout)->arr         = arr;

  return gsnv;
}

struct gkyl_array*
gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec* gsnv)
{
  N_Vector nvin = gsnv->nvec;
  return NV_CONTENT_GKZ(nvin)->arr;
}

void
gkyl_sundials_nvec_release(struct gkyl_sundials_nvec* gsnv)
{
  N_Vector nvin = gsnv->nvec;
  snvec_destroy(nvin);
  gkyl_free(gsnv);
}

void
gkyl_sundials_release(struct gkyl_sundials *gksun)
{
  gkyl_free(gksun);
}
