/** 
 *
 * NVECTOR wrapper of the Gkeyll data structure.
 *
 */

#ifdef GKYL_HAVE_SUNDIALS

#include <gkyl_sundials.h>
#include <gkyl_sundials_priv.h>
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
  content->is_passive  = SUNFALSE;
  content->red_mem     = NULL;

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
  NV_CONTENT_GKZ(nvout)->is_passive  = NV_CONTENT_GKZ(nvin)->is_passive;
  NV_CONTENT_GKZ(nvout)->red_mem     = NV_CONTENT_GKZ(nvin)->red_mem;
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
 * Set some entries 0. Needed by SUNDIALS but may be removed in
 * future SUNDIALS releases.
 *
 * @param v Input Nvector.
 * @param x Input index.
 * @param y Input index.
 */
static void
snvec_space(N_Vector v, sunindextype* x, sunindextype* y)
{
  *x = 0;
  *y = 0;
}

/**
 * Return the number of DOF wrapped by the nvector, excluding ghost cells.
 * future SUNDIALS releases.
 *
 * @param v Input Nvector.
 * @return Number of DOF wrapped by v (excluding ghosts).
 */
static sunindextype
snvec_get_length(N_Vector x)
{
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;
  return xarr->ncomp * local_range->volume;
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
static void
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
static void
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
static void
snvec_scale(sunrealtype c, N_Vector x, N_Vector z)
{
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;

  gkyl_array_set_range(zarr, c, xarr, local_range);
}

/**
 * Component-wise division of two SUNDIALS NVECTORS.
 *   w = u / v.
 *
 * @param u Numerator Nvector.
 * @param v Denominator Nvector.
 * @param w Output Nvector.
 */
static void
snvec_div(N_Vector u, N_Vector v, N_Vector w)
{
  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_array *w_arr = NV_CONTENT_GKZ(w)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(w_arr, GKYL_DIV, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(0.0), v_arr, local_range);
}

/**
 * Component-wise absolute value.
 *   v_k = |u_k|.
 *
 * @param u Input Nvector.
 * @param v Output Nvector.
 */
static void
snvec_abs(N_Vector u, N_Vector v)
{
  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(v_arr, GKYL_ABS, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(1.0), v_arr, local_range);
}

/**
 * Component-wise reciprocal:
 *   v_k = 1/u_k.
 *
 * @param u Input Nvector.
 * @param v Output Nvector.
 */
static void
snvec_inv(N_Vector u, N_Vector v)
{
  struct gkyl_array *u_arr       = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr       = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(v_arr, GKYL_INV, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(1.0), v_arr, local_range);
}

/**
 * Add a constant component-wise:
 *   v_k = u_k + x.
 *
 * @param u Input Nvector.
 * @param x Scala value to add.
 * @param v Output Nvector.
 */
static void
snvec_add_const(N_Vector u, sunrealtype x, N_Vector v)
{
  struct gkyl_array *u_arr       = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr       = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_copy_range(v_arr, u_arr, local_range);

  sunindextype ncomp = u_arr->ncomp;

  for (sunindextype i = 0; i < ncomp; ++i)
    gkyl_array_shiftc_range(v_arr, x, i, local_range);
}

/**
 * Compute the component-wise dot product between two SUNDIALS NVECTORS.
 *
 * @param x Input Nvector.
 * @param y Input Nvector.
 * @return Dot product of x and y.
 */
static sunrealtype
snvec_dot_product(N_Vector x, N_Vector y)
{
  bool is_passive = NV_CONTENT_GKZ(x)->is_passive;
  if (is_passive)
    return 0.0;

  struct gkyl_array *x_arr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *y_arr = NV_CONTENT_GKZ(y)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(y)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(x)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(y)->use_gpu;

  // Sum reduce (component-wise) x_i^{(k)} * y_i^{(k)}.
  int ncomp = x_arr->ncomp;
  gkyl_array_reduce_weighted_range(red_mem->red_local, x_arr, y_arr, GKYL_SUM, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp*sizeof(double));

  sunrealtype dot_prod = 0.0;
  for (sunindextype i = 0; i < ncomp; ++i)
    dot_prod += red_mem->red_global_ho[i];

  return dot_prod;
}

/**
 * Compute the absolute max norm of a SUNDIALS NVECTOR (component-wise).
 *
 * @param u Input Nvector.
 * @return Absolute max norm.
 */
static sunrealtype
snvec_abs_max_norm(N_Vector u)
{
  bool is_passive = NV_CONTENT_GKZ(u)->is_passive;
  if (is_passive)
    return 0.0;

  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(u)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(u)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(u)->use_gpu;

  int ncomp = u_arr->ncomp;

  gkyl_array_reduce_range(red_mem->red_local, u_arr, GKYL_ABS_MAX, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_MAX, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double));

  sunrealtype u_abs_max = -1.0;
  for (sunindextype i = 0; i < ncomp; ++i)
    u_abs_max = fmax(u_abs_max, red_mem->red_global_ho[i]);

  return u_abs_max;
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
static sunrealtype
snvec_wrms_norm(N_Vector nvx, N_Vector nvwgt)
{
  bool is_passive = NV_CONTENT_GKZ(nvx)->is_passive;
  if (is_passive)
    return 0.0;

  struct gkyl_array *nvx_arr = NV_CONTENT_GKZ(nvx)->arr;
  struct gkyl_array *nvwgt_arr = NV_CONTENT_GKZ(nvwgt)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(nvwgt)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(nvx)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(nvx)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(nvx)->use_gpu;

  int ncomp = nvx_arr->ncomp;

  // Reduce over cells.
  gkyl_array_reduce_weighted_range(red_mem->red_local, nvx_arr, nvwgt_arr, GKYL_RMS, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double));

  // Reduce over components.
  //  sunrealtype red_out = 0.0;
  //  for (sunindextype i = 0; i < ncomp; ++i) red_out += red_ho[i];
  // Use the 0th component because each component should have the same result.
  sunrealtype red_out = red_mem->red_global_ho[0];

  red_out = SUNRsqrt(red_out / local_range->volume);

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
  if (nvec == 0) {
    fprintf(stderr, "Error: Creating empty nvector.\n");
    assert(false);
  }

  // Attach operations.
  // Constructors, destructors, and utility operations.
  nvec->ops->nvclone      = snvec_clone;
  nvec->ops->nvcloneempty = snvec_clone_empty;
  nvec->ops->nvdestroy    = snvec_destroy;

  // Data operations.
  nvec->ops->nvspace     = snvec_space;
  nvec->ops->nvgetlength = snvec_get_length;
  nvec->ops->nvlinearsum = snvec_linear_sum;
  nvec->ops->nvconst     = snvec_const;
  nvec->ops->nvscale     = snvec_scale;
  nvec->ops->nvdiv       = snvec_div;
  nvec->ops->nvabs       = snvec_abs;
  nvec->ops->nvinv       = snvec_inv;
  nvec->ops->nvaddconst  = snvec_add_const;
  nvec->ops->nvdotprod   = snvec_dot_product;
  nvec->ops->nvmaxnorm   = snvec_abs_max_norm;
  nvec->ops->nvwrmsnorm  = snvec_wrms_norm;

  // Create content.
  content = 0;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == 0) {
    fprintf(stderr, "Error: Allocating content for nvector.\n");
    assert(false);
  }

  // Attach content.
  nvec->content = content;

  // Attach gkyl_array.
  content->own_vector  = SUNFALSE;
  content->use_gpu     = SUNFALSE;
  content->arr         = 0;
  content->comm        = 0;
  content->local_range = 0;
  content->is_passive  = SUNFALSE;
  content->red_mem     = 0;

  return nvec;
}

struct gkyl_array*
snvec_get_array(N_Vector nvin)
{
  return NV_CONTENT_GKZ(nvin)->arr;
}

struct gkyl_array*
smanynvec_get_array(N_Vector manynvin, int nvidx)
{
  N_Vector nvin = N_VGetSubvector_ManyVector(manynvin, nvidx);
  return snvec_get_array(nvin);
}

struct gkyl_sundials *
gkyl_sundials_new(bool use_gpu)
{
  struct gkyl_sundials *gksun = gkyl_malloc(sizeof(*gksun));

  gksun->use_gpu = use_gpu;

  // Create the SUNDIALS context object.
  SUNContext_Create(SUN_COMM_NULL, &gksun->sunctx);

  return gksun;
}

/**
 * Translate the enum used to indicate LSRK method in Gkeyll,
 * to the LSRK method types in SUNDIALS.
 *
 * @param gk_lsrk_method LSRK method Gkeyll enum.
 * @return SUNDIALS LSRK method flag.
 */
static ARKODE_LSRKMethodType
translate_gk_to_sundials_rk_method(enum gkyl_sundials_rk_method gk_rk_method)
{
  switch (gk_rk_method) {
    case GKYL_RK_METHOD_SSP_3_3:
      assert(false); // Shouldn't be necessary.
      break;
    case GKYL_SUNDIALS_LSRK_METHOD_SSP_S_2: // Optimal 2nd order s-stage SSP RK method.
      return ARKODE_LSRK_SSP_S_2;
      break;
    case GKYL_SUNDIALS_LSRK_METHOD_SSP_S_3: // Optimal 3rd order s-stage SSP RK method.
      return ARKODE_LSRK_SSP_S_3;
      break;
    case GKYL_SUNDIALS_LSRK_METHOD_SSP_10_4: // Optimal 4th order 10-stage SSP RK method.
      return ARKODE_LSRK_SSP_10_4;
      break;
    case GKYL_SUNDIALS_LSRK_METHOD_RKC_2: // 2nd order Runge-Kutta-Chebyshev (RKC).
      return ARKODE_LSRK_RKC_2;
      break;
    case GKYL_SUNDIALS_LSRK_METHOD_RKL_2: // 2nd order Runge-Kutta-Legendre (RKL).
      return ARKODE_LSRK_RKL_2;
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

static void
gkyl_sundials_stepper_init_ssp_rk33(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize Gkeyll's native 3rd order 3-stage SSP RK but through SUNDIALS.
  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\nError: gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call ERKStepCreate to initialize the ARK timestepper module and
  // specify the right-hand side function in y'=f(t,y), the initial time
  // T0, and the initial dependent variable vector y.
  gksun->arkode_mem = ERKStepCreate(gksun->dfdt_func, inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem, "ERKStepCreate", 0);

  // Set routines.
  inp->app_ctx->arkode_mem = gksun->arkode_mem;
  flag = ARKodeSetUserData(gksun->arkode_mem, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Set the initial step size.
  sunrealtype dt_stab = 0.1; // Recalculated and passed with _reset method.
  flag = ARKodeSetInitStep(gksun->arkode_mem, dt_stab);
  sundials_check_flag(&flag, "ARKodeSetInitStep", 1);

  // Set CFL stable time step size function.
  flag = ARKodeSetStabilityFn(gksun->arkode_mem, gksun->cfl_stable_dt_func, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetStabilityFn", 1);

  // Set CFL safety factor one to ensure exact use of the stable time step size.
  flag = ARKodeSetCFLFraction(gksun->arkode_mem, SUN_RCONST(1.0));
  sundials_check_flag(&flag, "ARKodeSetCFLFraction", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  // Insert the SSP33 Butcher tableau
  ARKodeButcherTable B_ssp33;
  B_ssp33 = ARKodeButcherTable_Alloc(3, SUNTRUE);
  sundials_check_flag((void*)(B_ssp33), "ARKodeButcherTable_Alloc", 0);
  B_ssp33->A[1][0] = SUN_RCONST(1.0);
  B_ssp33->A[2][0] = SUN_RCONST(0.25);
  B_ssp33->A[2][1] = SUN_RCONST(0.25);
  B_ssp33->b[0]    = SUN_RCONST(1.0/6.0);
  B_ssp33->b[1]    = SUN_RCONST(1.0/6.0);
  B_ssp33->b[2]    = SUN_RCONST(2.0/3.0);
  B_ssp33->c[1]    = SUN_RCONST(1.0);
  B_ssp33->c[2]    = SUN_RCONST(0.5);
  B_ssp33->q       = 3;
  B_ssp33->p       = 3;    // dummy embedding order
  B_ssp33->d[0]    = SUN_RCONST(1.0/6.0); // dummy embedding
  B_ssp33->d[1]    = SUN_RCONST(1.0/6.0);
  B_ssp33->d[2]    = SUN_RCONST(2.0/3.0);

  flag = ERKStepSetTable(gksun->arkode_mem, B_ssp33);
  sundials_check_flag(&flag, "ERKStepSetTable", 1);

  // Set pre/post processing methods in arkode mem.
  flag = ARKodeSetPreprocessRHSFn(gksun->arkode_mem, gksun->pre_process_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPreprocessRHSFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem, gksun->post_process_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPostprocessStepFailFn(gksun->arkode_mem, gksun->post_process_failed_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStepFailFn", 1);

  // Free the Butcher tableau.
  ARKodeButcherTable_Free(B_ssp33);
}

static void
gkyl_sundials_stepper_init_ssp_rk(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize a SSP RK time stepper (not 3-stage SSP RK3).
  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\nError: gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call LSRKStepCreateSSP to initialize the ARK timestepper module and
  // specify the right-hand side function in dfdt, the initial time
  // t_curr, and the initial dependent variable vector nvin.
  gksun->arkode_mem = LSRKStepCreateSSP(gksun->dfdt_func, inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem, "LSRKStepCreateSSP", 0);

  // Set user data (app pointer).
  inp->app_ctx->arkode_mem = gksun->arkode_mem;
  flag = ARKodeSetUserData(gksun->arkode_mem, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Specify tolerances.
  flag = ARKodeSStolerances(gksun->arkode_mem, inp->rel_tol, inp->abs_tol);
  sundials_check_flag(&flag, "ARKStepSStolerances", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  // Specify the SSP method.
  flag = LSRKStepSetSSPMethod(gksun->arkode_mem, translate_gk_to_sundials_rk_method(inp->rk_method));
  sundials_check_flag(&flag, "LSRKStepSetSSPMethod", 1);

  // Specify the number of SSP stages.
  flag = LSRKStepSetNumSSPStages(gksun->arkode_mem, inp->num_stages);
  sundials_check_flag(&flag, "LSRKStepSetNumSSPStages", 1);

  // Attach the error function.
  flag = ARKodeWFtolerances(gksun->arkode_mem, gksun->snvec_efun_cell_norm_func);
  sundials_check_flag(&flag, "ARKodeWFtolerances", 1);

  // Set pre/post processing methods in arkode mem.
  flag = ARKodeSetPreprocessRHSFn(gksun->arkode_mem, gksun->pre_process_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPreprocessRHSFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem, gksun->post_process_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPostprocessStepFailFn(gksun->arkode_mem, gksun->post_process_failed_rk_stage_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStepFailFn", 1);
}

static void
gkyl_sundials_stepper_init_sts(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize the Super Time Stepping integrator.

  // Set default values if user didn't provide a value.
  if (inp->dee_max_iter == 0)
    inp->dee_max_iter = 1e3;
  
  if (fabs(inp->dee_rel_tol) < 1e-16)
    inp->dee_rel_tol = 0.01;
  
  if (inp->dee_num_init_warmups == 0)
    inp->dee_num_init_warmups = 1e2;
  
  if (inp->dee_num_succ_warmups != 0)
    inp->dee_num_succ_warmups = inp->dee_num_succ_warmups;
  
  if (inp->dee_frequency == 0)
    inp->dee_frequency = 10;
  // Finished setting default values.

  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\nError: gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call LSRKStepCreateSTS to initialize the ARK timestepper module and
  // specify the right-hand side function in dfdt, the initial time
  // t_curr, and the initial dependent variable vector nvin.
  gksun->arkode_mem = LSRKStepCreateSTS(gksun->dfdt_sts_func, inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem, "LSRKStepCreateSTS", 0);

  // Set user data (app pointer).
  inp->app_ctx->arkode_mem = gksun->arkode_mem;
  flag = ARKodeSetUserData(gksun->arkode_mem, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Specify tolerances.
  flag = ARKodeSStolerances(gksun->arkode_mem, inp->rel_tol, inp->abs_tol);
  sundials_check_flag(&flag, "ARKStepSStolerances", 1);

  gksun->dom_eig_est = 0;
  if (inp->dee_by_gkeyll) {
    // Gkeyll estimates the dominant eigenvalue.
    flag = LSRKStepSetDomEigFn(gksun->arkode_mem, gksun->sts_dom_eig_func);
    sundials_check_flag(&flag, "LSRKStepSetDomEigFn", 1);
  }
  else {
    // Set the initial eigenvector for DEE. Should reset after ICs are set.
    gksun->dom_eig_est = SUNDomEigEstimator_Power(nvin, inp->dee_max_iter,
      inp->dee_rel_tol, nvin->sunctx);
    sundials_check_flag(&flag, "SUNDomEigEstimator_Power", 0);
  
    flag = LSRKStepSetDomEigEstimator(gksun->arkode_mem, gksun->dom_eig_est);
    sundials_check_flag(&flag, "LSRKStepSetDomEigEstimator", 1);
  
    flag = LSRKStepSetNumDomEigEstInitPreprocessIters(gksun->arkode_mem, inp->dee_num_init_warmups);
    sundials_check_flag(&flag, "LSRKStepSetNumDomEigEstInitPreprocessIters", 1);
  
    flag = LSRKStepSetNumDomEigEstPreprocessIters(gksun->arkode_mem, inp->dee_num_succ_warmups);
    sundials_check_flag(&flag, "LSRKStepSetNumDomEigEstPreprocessIters", 1);
  }

  // Specify after how many successful steps dom_eig is recomputed.
  // Note that nsteps = 0 refers to constant dominant eigenvalue.
  flag = LSRKStepSetDomEigFrequency(gksun->arkode_mem, inp->dee_frequency);
  sundials_check_flag(&flag, "LSRKStepSetDomEigFrequency", 1);

  // Specify max number of stages allowed.
  flag = LSRKStepSetMaxNumStages(gksun->arkode_mem, inp->max_num_stages);
  sundials_check_flag(&flag, "LSRKStepSetMaxNumStages", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem, inp->max_steps);
  sundials_check_flag(&flag, "LSRKStepSetMaxNumStages", 1);

  // Specify safety factor for user provided dom_eig.
  flag = LSRKStepSetDomEigSafetyFactor(gksun->arkode_mem, inp->dee_safety_fac);
  sundials_check_flag(&flag, "LSRKStepSetDomEigSafetyFactor", 1);

  // Specify the STS method.
  flag = LSRKStepSetSTSMethod(gksun->arkode_mem, translate_gk_to_sundials_rk_method(inp->rk_method));
  sundials_check_flag(&flag, "LSRKStepSetSTSMethod", 1);
}


void
gkyl_sundials_stepper_init(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Set default values if user didn't provide a value.
  if ( (inp->rk_method == GKYL_SUNDIALS_METHOD_NONE) ||
       ((inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_SSP_S_3) && (inp->num_stages == 3)) )
    inp->rk_method = GKYL_RK_METHOD_SSP_3_3;
  
  if (inp->max_steps == 0)
    inp->max_steps = 100000;
  // Finished setting default values.

  gksun->stepper_inp = inp; // Store stepper inputs.
  gksun->app_ctx = inp->app_ctx; // Copy pointer to app pointer.

  if (inp->rk_method == GKYL_RK_METHOD_SSP_3_3) {
    // Gkeyll's native 3rd order 3-stage SSP RK.
    gkyl_sundials_stepper_init_ssp_rk33(gksun, inp);
  }
  else if ( (inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_SSP_S_2) ||
            (inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_SSP_S_3) ||
            (inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_SSP_10_4) ) {
    if (inp->num_stages == 3) {
      // Gkeyll's native 3rd order 3-stage SSP RK, without embedding,
      // adapting dt using Gkeyll's CFL constraint.
      inp->rk_method = GKYL_RK_METHOD_SSP_3_3;
      gkyl_sundials_stepper_init_ssp_rk33(gksun, inp);
    }
    else {
      // Other SSP RK methods (with embedding to adapt dt).
      gkyl_sundials_stepper_init_ssp_rk(gksun, inp);
    }
  }
  else {
    // Super time stepping.
    gkyl_sundials_stepper_init_sts(gksun, inp);
  }

}

void
gkyl_sundials_arkode_reset(struct gkyl_sundials *gksun, double time,
  struct gkyl_sundials_nvec *gsmanynv, struct gkyl_sundials_nvec *gsmanynv_buff)
{
  N_Vector manynvin = gsmanynv->nvec;

  int flag;

  flag = ARKodeReset(gksun->arkode_mem, time, manynvin);
  sundials_check_flag(&flag, "ARKodeReset", 1);

  if (gksun->stepper_inp->rk_method == GKYL_RK_METHOD_SSP_3_3) {
    // Estimate the CFL stable dt.
    // Create a temporary NVector.
    N_Vector manynvbuff = gsmanynv_buff->nvec;
  
    // Compute dt.
    double dt_init;
    flag = gksun->dfdt_func(time, manynvin, manynvbuff, gksun->app_ctx);
    flag = gksun->cfl_stable_dt_func(manynvin, time, &dt_init, gksun->app_ctx);

    // Set the initial step size.
    flag = ARKodeSetInitStep(gksun->arkode_mem, dt_init);
    sundials_check_flag(&flag, "ARKodeSetInitStep", 1);
  }
  else if ( (gksun->stepper_inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_RKC_2) ||
            (gksun->stepper_inp->rk_method == GKYL_SUNDIALS_LSRK_METHOD_RKL_2) ) {
    if ( !(gksun->stepper_inp->dee_by_gkeyll) ) {
      // Pass ICs to the eigenvalue estimate.
      flag = SUNDomEigEstimator_SetInitialGuess(gksun->dom_eig_est, manynvin);
      sundials_check_flag(&flag, "SUNDomEigEstimator_SetInitialGuess", 1);
    }
  }

}

int
gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double *t_curr)
{
  N_Vector nvin = gsnv->nvec;

  // Call integrator to evolve the solution to time t_new.
  int flag = ARKodeEvolve(gksun->arkode_mem, t_new, nvin, t_curr, ARK_ONE_STEP);
  sundials_check_flag(&flag, "ARKodeEvolve", 1);

  return flag;
}

long
gkyl_sundials_get_num_error_test_failures(struct gkyl_sundials *gksun)
{
  long num_fail;
  int flag = ARKodeGetNumErrTestFails(gksun->arkode_mem, &num_fail); 
  sundials_check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  return num_fail;
}

long
gkyl_sundials_get_num_rhs_evals(struct gkyl_sundials *gksun)
{
  long num_evals;
  int flag = ARKodeGetNumRhsEvals(gksun->arkode_mem, 0, &num_evals); 
  sundials_check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  return num_evals;
}

void
gkyl_sundials_release(struct gkyl_sundials *gksun)
{
  gkyl_free(gksun);
}

struct gkyl_sundials_nvec*
gkyl_sundials_nvec_new(struct gkyl_sundials *gksun, struct gkyl_array *arr,
  struct gkyl_comm *comm, struct gkyl_range *local_range, bool is_passive)
{
  struct gkyl_sundials_nvec *gsnv = gkyl_malloc(sizeof(*gsnv));

  // Allocate memory needed for reductions
  bool use_gpu = gkyl_array_is_cu_dev(arr);
  int ncomp = arr->ncomp;
  if (use_gpu) {
    gsnv->red_mem.red_local  = gkyl_cu_malloc(ncomp * sizeof(double));
    gsnv->red_mem.red_global = gkyl_cu_malloc(ncomp * sizeof(double));
  }
  else {
    gsnv->red_mem.red_local  = gkyl_malloc(ncomp * sizeof(double));
    gsnv->red_mem.red_global = gkyl_malloc(ncomp * sizeof(double));
  }
  gsnv->red_mem.red_global_ho = gkyl_malloc(ncomp * sizeof(double));

  // Allocate the Sundials NVector.
  N_Vector nvout;
  nvout = 0;
  nvout = snvec_new_empty(gksun->sunctx);
  if (nvout == 0) {
    fprintf(stderr, "Error: Creating new empty nvector.\n");
    assert(false);
  }

  NV_CONTENT_GKZ(nvout)->own_vector = SUNFALSE;
  NV_CONTENT_GKZ(nvout)->use_gpu = gkyl_array_is_cu_dev(arr);
  NV_CONTENT_GKZ(nvout)->comm = comm;
  NV_CONTENT_GKZ(nvout)->local_range = local_range;
  NV_CONTENT_GKZ(nvout)->is_passive = is_passive;
  NV_CONTENT_GKZ(nvout)->arr = arr;
  NV_CONTENT_GKZ(nvout)->red_mem = &gsnv->red_mem;

  gsnv->nvec = nvout;
  return gsnv;
}

struct gkyl_array*
gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *gsnv)
{
  N_Vector nvin = gsnv->nvec;
  return snvec_get_array(nvin);
}

void
gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *gsnv)
{
  N_Vector nvin = gsnv->nvec;
  struct gkyl_array *arr = snvec_get_array(nvin);
  bool use_gpu = gkyl_array_is_cu_dev(arr);

  if (use_gpu) {
    gkyl_cu_free(gsnv->red_mem.red_local);
    gkyl_cu_free(gsnv->red_mem.red_global);
  }
  else {
    gkyl_free(gsnv->red_mem.red_local);
    gkyl_free(gsnv->red_mem.red_global);
  }
  gkyl_free(gsnv->red_mem.red_global_ho);

  snvec_destroy(nvin);
  gkyl_free(gsnv);
}

struct gkyl_sundials_nvec*
gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun, int num_nvector,
  struct gkyl_sundials_nvec *gsnv_arr[])
{
  struct gkyl_sundials_nvec *gsmanynv = gkyl_malloc(sizeof(*gsmanynv));
  
  // Create array of NVectors.
  N_Vector nvarr[num_nvector];
  for (int i=0; i<num_nvector; i++)
    nvarr[i] = gsnv_arr[i]->nvec;

  // Allocate the Sundials ManyNVector.
  N_Vector nvout = 0;
  nvout = N_VNew_ManyVector(num_nvector, nvarr, gksun->sunctx);
  if (nvout == 0) {
    fprintf(stderr, "Error: Creating new many nvector.\n");
    assert(false);
  }

  gsmanynv->nvec = nvout;

  return gsmanynv;
}

int 
gkyl_sundials_many_nvec_get_num_subvec(struct gkyl_sundials *gksun, struct gkyl_sundials_nvec* gsmanynv)
{
  return N_VGetNumSubvectors_ManyVector(gsmanynv->nvec);
}

void
gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv)
{
  N_Vector nvin = gsmanynv->nvec;
  N_VDestroy(nvin);
  gkyl_free(gsmanynv);
}

#else

struct gkyl_sundials*
gkyl_sundials_new(int ncomp, bool use_gpu)
{
  fprintf(stderr, "\nGkeyll was not built with SUNDIALS. Use Gkeyll-native steppers or build Gkeyll with SUNDIALS.\n");
}

void
gkyl_sundials_stepper_init_ssp_rk(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Do nothing.
}

int
gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double t_curr)
{
  // Do nothing.
}

void
gkyl_sundials_release(struct gkyl_sundials *gksun)
{
  // Do nothing.
}

struct gkyl_sundials_nvec*
gkyl_sundials_nvec_new(struct gkyl_sundials *gksun, struct gkyl_array *arr,
  struct gkyl_comm *comm, struct gkyl_range *local_range, bool is_passive)
{
  // Do nothing.
}

struct gkyl_array*
gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *gsnv)
{
  // Do nothing.
}

void
gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *gsnv)
{
  // Do nothing.
}

struct gkyl_sundials_nvec*
gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun, int num_nvector,
  struct gkyl_sundials_nvec *gsnv_arr[])
{
  // Do nothing.
}


int 
gkyl_sundials_many_nvec_get_num_subvec(struct gkyl_sundials *gksun, struct gkyl_sundials_nvec* gsmanynv)
{
  return 0;
}

void
gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv)
{
  // Do nothing.
}
#endif
