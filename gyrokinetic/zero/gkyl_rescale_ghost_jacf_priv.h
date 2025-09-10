#pragma once

// Private header for rescale_ghost_jacf updater, not for direct use in user code.

#include <gkyl_rescale_ghost_jacf.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_deflate_surf_kernels.h>
#include <gkyl_inflate_surf_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*deflate_surf_op_t)(const double *fld, double* deflated_fld);
typedef void (*inflate_surf_op_t)(const double *deflated_fld, double* fld);

typedef struct { deflate_surf_op_t kernels[2]; } deflate_surf_kern_list;  // For use in kernel tables.
typedef struct { deflate_surf_kern_list dirlist[3]; } dir_deflate_surf_kern_list;
typedef struct { dir_deflate_surf_kern_list edgedlist[4]; } edged_deflate_surf_kern_list;

typedef struct { inflate_surf_op_t kernels[2]; } inflate_surf_kern_list;  // For use in kernel tables.
typedef struct { inflate_surf_kern_list dirlist[4]; } edged_inflate_surf_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_deflate_surf_kern_list ser_deflate_surf_conf_list[] = {
  {
    .edgedlist =
    {
      {
        .dirlist =
        {
          { deflate_surfx_lower_1x_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_lower_2x_ser_p1, NULL },
          { deflate_surfy_lower_2x_ser_p1, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_lower_3x_ser_p1, NULL },
          { deflate_surfy_lower_3x_ser_p1, NULL },
          { deflate_surfz_lower_3x_ser_p1, NULL },
        }
      },
      {
        .dirlist =
        {
          { NULL, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
    }
  },
  {
    .edgedlist =
    {
      {
        .dirlist =
        {
          { deflate_surfx_upper_1x_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_upper_2x_ser_p1, NULL },
          { deflate_surfy_upper_2x_ser_p1, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_upper_3x_ser_p1, NULL },
          { deflate_surfy_upper_3x_ser_p1, NULL },
          { deflate_surfz_upper_3x_ser_p1, NULL },
        }
      },
      {
        .dirlist =
        {
          { NULL, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
    }
  }
};

GKYL_CU_D
static const edged_deflate_surf_kern_list ser_deflate_surf_phase_list[] = {
  {
    .edgedlist =
    {
      {
        .dirlist =
        {
          { deflate_surfx_lower_1x1v_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_lower_1x2v_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_lower_2x2v_ser_p1, NULL },
          { deflate_surfy_lower_2x2v_ser_p1, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_lower_3x2v_ser_p1, NULL },
          { deflate_surfy_lower_3x2v_ser_p1, NULL },
          { deflate_surfz_lower_3x2v_ser_p1, NULL },
        }
      }
    }
  },
  {
    .edgedlist =
    {
      {
        .dirlist =
        {
          { deflate_surfx_upper_1x1v_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_upper_1x2v_ser_p1, NULL },
          { NULL, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_upper_2x2v_ser_p1, NULL },
          { deflate_surfy_upper_2x2v_ser_p1, NULL },
          { NULL, NULL },
        }
      },
      {
        .dirlist =
        {
          { deflate_surfx_upper_3x2v_ser_p1, NULL },
          { deflate_surfy_upper_3x2v_ser_p1, NULL },
          { deflate_surfz_upper_3x2v_ser_p1, NULL },
        }
      }
    }
  }
};

GKYL_CU_D
static const edged_inflate_surf_kern_list ser_inflate_surf_phase_list[] = {
  {
    .dirlist =
    {
      { inflate_surfx_1x1v_ser_p1, NULL },
      { NULL, NULL },
      { NULL, NULL },
    }
  },
  {
    .dirlist =
    {
      { inflate_surfx_1x2v_ser_p1, NULL },
      { NULL, NULL },
      { NULL, NULL },
    }
  },
  {
    .dirlist =
    {
      { inflate_surfx_2x2v_ser_p1, NULL },
      { inflate_surfy_2x2v_ser_p1, NULL },
      { NULL, NULL },
    }
  },
  {
    .dirlist =
    {
      { inflate_surfx_3x2v_ser_p1, NULL },
      { inflate_surfy_3x2v_ser_p1, NULL },
      { inflate_surfz_3x2v_ser_p1, NULL },
    }
  }
};

// Need a special kernel for cdim=1:
GKYL_CU_D
static void
conf_inv_op_0x(const double *A, double *A_inv)
{
  A_inv[0] = 1.0/A[0];
}

GKYL_CU_D
static void
conf_mul_op_0x(const double *f, const double *g, double *fg)
{
  fg[0] = f[0]*g[0];
};

GKYL_CU_D
static void
conf_mul_op_0x_1x1v_gkhybrid(const double *f, const double *g, double *fg)
{
  fg[0] = f[0]*g[0];
  fg[1] = f[0]*g[1];
  fg[2] = f[0]*g[2];
  fg[3] = f[0]*g[3];
  fg[4] = f[0]*g[4];
  fg[5] = f[0]*g[5];
};

GKYL_CU_D
static void
conf_mul_op_0x_1x2v_gkhybrid(const double *f, const double *g, double *fg)
{
  fg[0] = f[0]*g[0];
  fg[1] = f[0]*g[1];
  fg[2] = f[0]*g[2];
  fg[3] = f[0]*g[3];
  fg[4] = f[0]*g[4];
  fg[5] = f[0]*g[5];
  fg[6] = f[0]*g[6];
  fg[7] = f[0]*g[7];
  fg[8] = f[0]*g[8];
  fg[9] = f[0]*g[9];
  fg[10] = f[0]*g[10];
  fg[11] = f[0]*g[11];
};

struct gkyl_rescale_ghost_jacf_kernels {
  deflate_surf_op_t deflate_conf_ghost_op; // Project conf-field onto plane at lower/upper surface.
  deflate_surf_op_t deflate_conf_skin_op; // Project conf-field onto plane at lower/upper surface.
  deflate_surf_op_t deflate_phase_ghost_op; // Project phase-field onto plane at lower/upper surface.
  inflate_surf_op_t inflate_phase_ghost_op; // Inflate phase-field surface field to volume.
  inv_op_t conf_inv_op; // Conf-space weak inversion (1/A) kernel (p=1 only).
  mul_op_t conf_mul_op; // Conf-space weak multiplication kernel.
  mul_op_t conf_phase_mul_op; // Conf-phase weak multiplication kernel.
};

// Primary struct in this updater.
struct gkyl_rescale_ghost_jacf {
  int dir; // Direction perpendicular to the boundary.
  enum gkyl_edge_loc edge; // Boundary edge (lower/upper).
  bool use_gpu; // Whether to run on the GPU or not.
  struct gkyl_rescale_ghost_jacf_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void
rescale_ghost_jacf_choose_kernel_cu(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

void
gkyl_rescale_ghost_jacf_advance_cu(const struct gkyl_rescale_ghost_jacf *up,
  const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf);
#endif

GKYL_CU_D
static void rescale_ghost_jacf_choose_kernel(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    rescale_ghost_jacf_choose_kernel_cu(kernels, dir, edge, cbasis, pbasis);
    return;
  }
#endif

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  enum gkyl_basis_type cbasis_type = cbasis->b_type, pbasis_type = pbasis->b_type;
  int poly_order = pbasis->poly_order;

  enum gkyl_edge_loc ghost_edge = edge == GKYL_LOWER_EDGE? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE;

  switch (pbasis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->deflate_phase_ghost_op = ser_deflate_surf_phase_list[ghost_edge].edgedlist[pdim-2].dirlist[dir].kernels[poly_order-1];
      kernels->inflate_phase_ghost_op = ser_inflate_surf_phase_list[pdim-2].dirlist[dir].kernels[poly_order-1];
      if (cdim == 1) {
        if (vdim == 1)
          kernels->conf_phase_mul_op = conf_mul_op_0x_1x1v_gkhybrid;
        else if (vdim == 2)
          kernels->conf_phase_mul_op = conf_mul_op_0x_1x2v_gkhybrid;
      }
      else
        kernels->conf_phase_mul_op = choose_mul_conf_phase_kern(pbasis_type, cdim-1, vdim, poly_order);
      break;
    default:
      assert(false);
      break;
  }

  switch (cbasis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->deflate_conf_skin_op = ser_deflate_surf_conf_list[edge].edgedlist[cdim-1].dirlist[dir].kernels[poly_order-1];
      kernels->deflate_conf_ghost_op = ser_deflate_surf_conf_list[ghost_edge].edgedlist[cdim-1].dirlist[dir].kernels[poly_order-1];
      kernels->conf_inv_op = cdim==1? conf_inv_op_0x : choose_ser_inv_kern(cdim-1, poly_order);
      kernels->conf_mul_op = cdim==1? conf_mul_op_0x : choose_ser_mul_kern(cdim-1, poly_order);
      break;
    default:
      assert(false);
      break;
  }

}
