// Private header for dg_eval_at_coord_proj updater. Not for direct use.
#pragma once

#include <assert.h>
#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dg_eval_at_coord_proj_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Function pointer type for projection kernels.
// coords[i] is the reference coordinate for the i-th evaluated direction.
typedef void (*eval_at_coord_t)(const double *coords, const double *fdo, double *ftar);

// For use in kernel tables.
typedef struct { eval_at_coord_t kernels[3]; } eval_at_coord_kern_list;  // [poly_order-1]

// Struct with double array for passing coordinates to CUDA kernels.
typedef struct { double c[GKYL_MAX_DIM]; } dg_evproj_struct_double_t;

// Struct with int array for passing cell index to CUDA kernels.
typedef struct { int c[GKYL_MAX_DIM]; } dg_evproj_struct_int_t;

// Serendipity kernels.
// Indexed as [ndim_do-1][dir_bitmask-1].kernels[poly_order-1]
//   ndim_do = 1 -> index 0   (1 valid mask:  {x}=1)
//   ndim_do = 2 -> index 1   (3 valid masks: {x}=1, {y}=2, {xy}=3)
//   ndim_do = 3 -> index 2   (7 valid masks: 1..7)
GKYL_CU_D
static const eval_at_coord_kern_list ser_eval_at_coord_list[3][7] = {
  // ndim_do = 1
  {
    // mask=1: eval_x (1x -> scalar)
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0, gkyl_dg_eval_at_coord_proj_1x_ser_p2_eval_dirs_0, NULL },
    // masks 2..7: not applicable for 1D
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 2
  {
    // mask=1: eval_x (2x -> 1D in y)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0 , NULL },
    // mask=2: eval_y (2x -> 1D in x)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1 , NULL },
    // mask=3: eval_xy (2x -> scalar)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01, gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01, NULL },
    // masks 4..7: not applicable for 2D
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 3
  {
    // mask=1: eval_x (3x -> 2D in y,z)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0  , NULL, NULL },
    // mask=2: eval_y (3x -> 2D in x,z)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1  , NULL, NULL },
    // mask=3: eval_xy (3x -> 1D in z)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01 , NULL, NULL },
    // mask=4: eval_z (3x -> 2D in x,y)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2  , NULL, NULL },
    // mask=5: eval_xz (3x -> 1D in y)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02 , NULL, NULL },
    // mask=6: eval_yz (3x -> 1D in x)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12 , NULL, NULL },
    // mask=7: eval_xyz (3x -> scalar)
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012, NULL, NULL },
  },
};

// Tensor kernels.
GKYL_CU_D
static const eval_at_coord_kern_list ten_eval_at_coord_list[3][7] = {
  // ndim_do = 1 (tensor p1 = ser p1, no tensor p2 for 1D)
  {
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 2
  {
    // mask=1: eval_x  (p1 reuses ser; p2 uses tensor kernel)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0 , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0 , NULL },
    // mask=2: eval_y
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1 , gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1 , NULL },
    // mask=3: eval_xy
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01, gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
    { NULL, NULL, NULL },
  },
  // ndim_do = 3
  {
    // mask=1: eval_x
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0  , NULL, NULL },
    // mask=2: eval_y
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1  , NULL, NULL },
    // mask=3: eval_xy
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01 , NULL, NULL },
    // mask=4: eval_z
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2  , NULL, NULL },
    // mask=5: eval_xz
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02 , NULL, NULL },
    // mask=6: eval_yz
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12 , NULL, NULL },
    // mask=7: eval_xyz
    { gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012, NULL, NULL },
  },
};

GKYL_CU_DH static void
eval_at_coord_get_idx_do(const bool *is_eval, int ndim_do,
  const int *idx_tar, const int *cell_idx, int *idx_do)
{
  // Translate a target index to the corresponding donor index.
  int c = 0;
  for (int d=0; d<ndim_do; d++) {
    if (!is_eval[d])
      idx_do[d] = idx_tar[c++];
    else
      idx_do[d] = cell_idx[d];
  }
}

GKYL_CU_D static eval_at_coord_t
choose_ser_eval_at_coord_kern(int ndim_do, int dir_mask, int poly_order)
{
  assert(ndim_do >= 1 && ndim_do <= 3);
  assert(dir_mask >= 1 && dir_mask < (1 << ndim_do));
  assert(poly_order >= 1 && poly_order <= 3);
  return ser_eval_at_coord_list[ndim_do-1][dir_mask-1].kernels[poly_order-1];
}

GKYL_CU_D static eval_at_coord_t
choose_ten_eval_at_coord_kern(int ndim_do, int dir_mask, int poly_order)
{
  assert(ndim_do >= 1 && ndim_do <= 3);
  assert(dir_mask >= 1 && dir_mask < (1 << ndim_do));
  assert(poly_order >= 1 && poly_order <= 3);
  return ten_eval_at_coord_list[ndim_do-1][dir_mask-1].kernels[poly_order-1];
}

// Primary struct for this updater.
struct gkyl_dg_eval_at_coord_proj {
  int ndim_do;                 // Number of dimensions in donor field.
  int ndim_tar;                // Number of dimensions in target field.
  int num_basis_do;            // Number of monomials in donor basis.
  int num_basis_tar;           // Number of monomials in target basis.
  int num_eval_dirs;           // Number of directions being evaluated.
  int eval_dirs[GKYL_MAX_DIM]; // Which directions are evaluated.
  bool is_eval[GKYL_MAX_DIM];  // Is direction evaluated.
  int dir_mask;                // Bitmask encoding evaluated directions.
  bool use_gpu;                // Whether to run on GPU.
  eval_at_coord_t kernel;      // Projection kernel.
};

#ifdef GKYL_HAVE_CUDA
// Host-side wrapper for the GPU implementation of gkyl_dg_eval_at_coord_proj_advance.
void gkyl_dg_eval_at_coord_proj_advance_cu(struct gkyl_dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar);
#endif
