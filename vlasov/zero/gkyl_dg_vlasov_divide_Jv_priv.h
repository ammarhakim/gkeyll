// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*divide_Jv_t)(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J);
typedef void (*rescale_Jv_t)(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static const struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { divide_Jv_t kernels[4]; } gkyl_dg_divide_Jv_kern_list;
typedef struct { rescale_Jv_t kernels[4]; } gkyl_dg_rescale_Jv_kern_list;

// Velocity-space Jacobian division kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_divide_Jv_kern_list ser_divide_Jv_kernels[] = {
  // 1x kernels
  { NULL, vlasov_divide_Jv_1x1v_ser_p1, vlasov_divide_Jv_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_divide_Jv_1x2v_ser_p1, vlasov_divide_Jv_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_divide_Jv_1x3v_ser_p1, vlasov_divide_Jv_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_divide_Jv_2x1v_ser_p1, vlasov_divide_Jv_2x1v_ser_p2, NULL }, // 3
  { NULL, vlasov_divide_Jv_2x2v_ser_p1, vlasov_divide_Jv_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_divide_Jv_2x3v_ser_p1, vlasov_divide_Jv_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_divide_Jv_3x3v_ser_p1, NULL, NULL }, // 6
};

// Velocity-space Jacobian division kernels (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_divide_Jv_kern_list tensor_divide_Jv_kernels[] = {
  // 1x kernels
  { NULL, NULL, vlasov_divide_Jv_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, vlasov_divide_Jv_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, vlasov_divide_Jv_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, vlasov_divide_Jv_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, vlasov_divide_Jv_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, vlasov_divide_Jv_2x3v_tensor_p2, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Velocity-space Jacobian rescale kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_rescale_Jv_kern_list ser_rescale_Jv_kernels[] = {
  // 1x kernels
  { NULL, vlasov_rescale_Jv_1x1v_ser_p1, vlasov_rescale_Jv_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_rescale_Jv_1x2v_ser_p1, vlasov_rescale_Jv_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_rescale_Jv_1x3v_ser_p1, vlasov_rescale_Jv_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_rescale_Jv_2x1v_ser_p1, vlasov_rescale_Jv_2x1v_ser_p2, NULL }, // 3
  { NULL, vlasov_rescale_Jv_2x2v_ser_p1, vlasov_rescale_Jv_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_rescale_Jv_2x3v_ser_p1, vlasov_rescale_Jv_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_rescale_Jv_3x3v_ser_p1, NULL, NULL }, // 6
};

// Velocity-space Jacobian rescale kernels (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_rescale_Jv_kern_list tensor_rescale_Jv_kernels[] = {
  // 1x kernels
  { NULL, NULL, vlasov_rescale_Jv_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, vlasov_rescale_Jv_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, vlasov_rescale_Jv_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, vlasov_rescale_Jv_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, vlasov_rescale_Jv_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, vlasov_rescale_Jv_2x3v_tensor_p2, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static divide_Jv_t
choose_ser_divide_Jv_kern(int cdim, int vdim, int poly_order)
{
  return ser_divide_Jv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static divide_Jv_t
choose_tensor_divide_Jv_kern(int cdim, int vdim, int poly_order)
{
  return tensor_divide_Jv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static rescale_Jv_t
choose_ser_rescale_Jv_kern(int cdim, int vdim, int poly_order)
{
  return ser_rescale_Jv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static rescale_Jv_t
choose_tensor_rescale_Jv_kern(int cdim, int vdim, int poly_order)
{
  return tensor_rescale_Jv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}