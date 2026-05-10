// Private header: not for direct use.
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>

// Kernel signature: takes the limiter_fac, theta-pole flags for this cell on
// each side, face nodal geometry (single-face: lapse, shift, h_ij, det_h),
// and the three-cell (l, c, r) field_con. Writes the limited slope back into
// field_con_out (which may alias field_con_c for in-place).
typedef void (*slope_limiter_kern_t)(double limiter_fac,
  const int theta_pole_lo, const int theta_pole_up,
  const double *lapse_nodal_face, const double *shift_nodal_face,
  const double *h_ij_nodal_face, const double *det_h_face,
  const double *field_con_l, const double *field_con_c, const double *field_con_r,
  double *field_con_out);

typedef struct { slope_limiter_kern_t kernels[4]; } gkyl_slope_limiter_kern_list;

struct gkyl_dg_gr_maxwell_slope_limiter {
  struct gkyl_rect_grid conf_grid;
  int cdim;
  double limiter_fac;
  slope_limiter_kern_t kernels[3]; // one per direction
  int theta_pole_lo[GKYL_MAX_CDIM];
  int theta_pole_up[GKYL_MAX_CDIM];
  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_slope_limiter *on_dev;
};

// Kernel table: serendipity, indexed by [cdim-1].kernels[poly_order].
// Currently only 2D ser p=1 implemented; other entries NULL.
static const gkyl_slope_limiter_kern_list ser_slope_limiter_x_kernels[] = {
  { NULL, NULL, NULL, NULL },                              // cdim=1 (unused)
  { NULL, slope_limiter_x_2x_ser_p1, NULL, NULL },         // cdim=2: p0,p1,p2,p3
  { NULL, NULL, NULL, NULL },                              // cdim=3 (unused)
};

static const gkyl_slope_limiter_kern_list ser_slope_limiter_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, slope_limiter_y_2x_ser_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

static const gkyl_slope_limiter_kern_list ser_slope_limiter_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_DH
static slope_limiter_kern_t
choose_slope_limiter_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0) return ser_slope_limiter_x_kernels[cdim-1].kernels[poly_order];
      if (dir == 1) return ser_slope_limiter_y_kernels[cdim-1].kernels[poly_order];
      if (dir == 2) return ser_slope_limiter_z_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
  }
  return NULL;
}
