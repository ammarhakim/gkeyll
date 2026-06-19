// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_differentiate_ser.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Function pointer type for differentiation kernels.
// Computes the DG derivative of f in a single direction and at a given
// derivative order, writing the result into df.  Both arrays have
// num_basis entries.
typedef void (*differentiate_op_t)(const double *f, double *df);

// for use in kernel tables
typedef struct { differentiate_op_t kernels[4]; } differentiate_op_kern_list;

// Serendipity differentiation kernels.
// [ndim-1][dir][order-1].kernels[poly_order]
//   ndim  : 1..3
//   dir   : 0..ndim-1
//   order : 1..poly_order  (order > poly_order gives identically zero output)
GKYL_CU_D
static const differentiate_op_kern_list ser_differentiate_list[3][3][3] = {
  // ndim = 1
  {
    // dir0
    {
      { NULL, dg_differentiate_1d_ser_p1_dir0_ord1, dg_differentiate_1d_ser_p2_dir0_ord1, dg_differentiate_1d_ser_p3_dir0_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_1d_ser_p2_dir0_ord2, dg_differentiate_1d_ser_p3_dir0_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_1d_ser_p3_dir0_ord3 }, // ord3
    },
    // dir1: not applicable for 1D
    {
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
    },
    // dir2: not applicable for 1D
    {
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
    },
  },
  // ndim = 2
  {
    // dir0
    {
      { NULL, dg_differentiate_2d_ser_p1_dir0_ord1, dg_differentiate_2d_ser_p2_dir0_ord1, dg_differentiate_2d_ser_p3_dir0_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_2d_ser_p2_dir0_ord2, dg_differentiate_2d_ser_p3_dir0_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_2d_ser_p3_dir0_ord3 }, // ord3
    },
    // dir1
    {
      { NULL, dg_differentiate_2d_ser_p1_dir1_ord1, dg_differentiate_2d_ser_p2_dir1_ord1, dg_differentiate_2d_ser_p3_dir1_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_2d_ser_p2_dir1_ord2, dg_differentiate_2d_ser_p3_dir1_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_2d_ser_p3_dir1_ord3 }, // ord3
    },
    // dir2: not applicable for 2D
    {
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL },
    },
  },
  // ndim = 3
  {
    // dir0
    {
      { NULL, dg_differentiate_3d_ser_p1_dir0_ord1, dg_differentiate_3d_ser_p2_dir0_ord1, dg_differentiate_3d_ser_p3_dir0_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_3d_ser_p2_dir0_ord2, dg_differentiate_3d_ser_p3_dir0_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_3d_ser_p3_dir0_ord3 }, // ord3
    },
    // dir1
    {
      { NULL, dg_differentiate_3d_ser_p1_dir1_ord1, dg_differentiate_3d_ser_p2_dir1_ord1, dg_differentiate_3d_ser_p3_dir1_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_3d_ser_p2_dir1_ord2, dg_differentiate_3d_ser_p3_dir1_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_3d_ser_p3_dir1_ord3 }, // ord3
    },
    // dir2
    {
      { NULL, dg_differentiate_3d_ser_p1_dir2_ord1, dg_differentiate_3d_ser_p2_dir2_ord1, dg_differentiate_3d_ser_p3_dir2_ord1 }, // ord1
      { NULL, NULL,                                  dg_differentiate_3d_ser_p2_dir2_ord2, dg_differentiate_3d_ser_p3_dir2_ord2 }, // ord2
      { NULL, NULL,                                  NULL,                                 dg_differentiate_3d_ser_p3_dir2_ord3 }, // ord3
    },
  },
};

GKYL_CU_D
static differentiate_op_t
choose_ser_differentiate_kern(int ndim, int dir, int poly_order, int order)
{
  assert(ndim >= 1 && ndim <= 3);
  assert(dir >= 0 && dir < ndim);
  assert(order >= 1 && order <= 3);
  return ser_differentiate_list[ndim-1][dir][order-1].kernels[poly_order];
}
