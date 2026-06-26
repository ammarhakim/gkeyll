// Private header for dg_eval_at_coord_proj updater. Not for direct use.
#pragma once

#include <assert.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dg_eval_at_coord_proj_kernels.h>
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

// Struct with bool array for passing is_dir_eval to CUDA kernels.
typedef struct { bool c[GKYL_MAX_DIM]; } dg_evproj_struct_bool_t;

// Serendipity kernels.
// Indexed as [ndim_do-1][dir_bitmask-1].kernels[poly_order-1]
//   ndim_do = 1 -> index 0   (1 valid mask:  {x}=1)
//   ndim_do = 2 -> index 1   (3 valid masks: 1..3)
//   ndim_do = 3 -> index 2   (7 valid masks: 1..7)
//   ndim_do = 4 -> index 3   (15 valid masks: 1..15; p1 and p2)
//   ndim_do = 5 -> index 4   (31 valid masks: 1..31; p1 only)
//   ndim_do = 6 -> index 5   (63 valid masks: 1..63; p1 only)
// Unlisted entries are zero-initialized (NULL kernels).
GKYL_CU_D
static const eval_at_coord_kern_list ser_eval_at_coord_list[6][63] = {
  // ndim_do = 1: mask=1 valid; masks 2..63 zero-initialized
  {
    // mask=1: eval_x (1x -> scalar)
    { gkyl_dg_eval_at_coord_proj_1x_ser_p1_eval_dirs_0, gkyl_dg_eval_at_coord_proj_1x_ser_p2_eval_dirs_0, NULL },
  },
  // ndim_do = 2: masks 1..3 valid; masks 4..63 zero-initialized
  {
    // mask=1: eval_x (2x -> 1D in y)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0 , NULL },
    // mask=2: eval_y (2x -> 1D in x)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1 , gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1 , NULL },
    // mask=3: eval_xy (2x -> scalar)
    { gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01, gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01, NULL },
  },
  // ndim_do = 3: masks 1..7 valid, p1 only; masks 8..63 zero-initialized
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
  // ndim_do = 4: masks 1..15, p1 and p2; masks 16..63 zero-initialized
  {
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3   , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3   , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23  , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23  , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123 , gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123 , NULL },
    { gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123, gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123, NULL },
  },
  // ndim_do = 5: masks 1..31, p1 only; masks 32..63 zero-initialized
  {
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234, NULL, NULL },
  },
  // ndim_do = 6: masks 1..63, p1 only
  {
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5     , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45    , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345   , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345  , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345 , NULL, NULL },
    { gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345, NULL, NULL },
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

struct dg_ev_proj_kernels {
  eval_at_coord_t ev_ker;      // Projection kernel.
};

// Primary struct for this updater.
struct gkyl_dg_eval_at_coord_proj {
  int ndim_do;                 // Number of dimensions in donor field.
  int ndim_tar;                // Number of dimensions in target field.
  int num_basis_do;            // Number of monomials in donor basis.
  int num_basis_tar;           // Number of monomials in target basis.
  int num_eval_dirs;           // Number of directions being evaluated.
  int eval_dirs[GKYL_MAX_DIM]; // Which directions are evaluated.
  bool is_eval[GKYL_MAX_DIM];  // Is direction evaluated.
  bool use_gpu;                // Whether to run on GPU.
  struct dg_ev_proj_kernels *kers; // Projection kernel.
};

#ifdef GKYL_HAVE_CUDA
// Host-side wrapper for the GPU implementation of gkyl_dg_eval_at_coord_proj_advance.
void gkyl_dg_eval_at_coord_proj_advance_cu(struct gkyl_dg_eval_at_coord_proj *up, const double *eval_coords,
  const struct gkyl_rect_grid *grid, const bool *pick_lower, const int *known_index,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *fdo, struct gkyl_array *ftar);

struct dg_ev_proj_kernels* dg_eval_at_coord_choose_ker_cu(int cdim, int ndim,
  const struct gkyl_basis *basis, int num_eval_dirs, const int *eval_dirs);
#endif

GKYL_CU_DH static int
eval_dirs_to_mask(int num_eval_dirs, const int *eval_dirs)
{
  // Encode eval_dirs as a bitmask: bit d is set iff direction d is evaluated.
  // Used to index the dispatch table as [ndim_do-1][mask-1].
  int mask = 0;
  for (int i=0; i<num_eval_dirs; i++)
    mask |= (1 << eval_dirs[i]);
  return mask;
}

GKYL_CU_D static struct dg_ev_proj_kernels*
dg_eval_at_coord_choose_ker(bool use_gpu, int cdim, int ndim, const struct gkyl_basis *basis, int num_eval_dirs,
  const int *eval_dirs)
{
  // Choose the projection kernel.

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return dg_eval_at_coord_choose_ker_cu(cdim, ndim, basis, num_eval_dirs, eval_dirs);
#endif

  struct dg_ev_proj_kernels *kers = (struct dg_ev_proj_kernels *) gkyl_calloc(1, sizeof(struct dg_ev_proj_kernels));

  int dir_mask = eval_dirs_to_mask(num_eval_dirs, eval_dirs);

  int poly_order = basis->poly_order;

  assert(ndim >= 1 && ndim <= 6);
  assert(dir_mask >= 1 && dir_mask < (1 << ndim));
  assert(poly_order >= 1 && poly_order <= 3);

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->ev_ker = ser_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      kers->ev_ker = ten_eval_at_coord_list[ndim-1][dir_mask-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
  assert(kers->ev_ker);

  return kers;
}
