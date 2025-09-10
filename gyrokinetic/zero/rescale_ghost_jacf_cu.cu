/* -*- c++ -*- */

extern "C" {
#include <gkyl_rescale_ghost_jacf.h>
#include <gkyl_rescale_ghost_jacf_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_rescale_ghost_jacf_set_cu_ker_ptrs(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  int dir, enum gkyl_edge_loc edge, struct gkyl_basis cbasis, struct gkyl_basis pbasis)
{
  int cdim = cbasis.ndim, pdim = pbasis.ndim;
  int vdim = pdim-cdim;
  enum gkyl_basis_type cbasis_type = cbasis.b_type, pbasis_type = pbasis.b_type;
  int poly_order = cbasis.poly_order;

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
};

void
rescale_ghost_jacf_choose_kernel_cu(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  gkyl_rescale_ghost_jacf_set_cu_ker_ptrs<<<1,1>>>(kernels, dir, edge, *cbasis, *pbasis);
}

// CUDA kernel to copy ghost cell values to the adjacent skin (boundary) cells on the GPU.
__global__ static void
gkyl_rescale_ghost_jacf_advance_cu_ker(struct gkyl_rescale_ghost_jacf_kernels *kers,
  int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range conf_skin_r, const struct gkyl_range conf_ghost_r,
  const struct gkyl_range phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf)
{
  int sidx[GKYL_MAX_DIM]; // skin idx
  int gidx[GKYL_MAX_DIM]; // ghost idx

  int cdim = conf_skin_r.ndim;

  const int num_basis_conf_surf = 8; // MF 2024/11/12: Hardcoded to 2x p2 for now.
  const int num_basis_phase_surf = 24; // MF 2024/11/12: Hardcoded to 2x2v GkHybrid for now.

  // Loop over all points in the skin range using CUDA threads.
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x;
       linc < phase_ghost_r.volume; linc += blockDim.x * gridDim.x) {

    // Convert the linear index to a multi-dimensional index in the skin range.
    gkyl_sub_range_inv_idx(&phase_ghost_r, linc, gidx);

    // Get skin cell corresponding to this ghost cell.
    gkyl_copy_int_arr(cdim, gidx, sidx);
    sidx[dir] = edge == GKYL_LOWER_EDGE? gidx[dir]+1 : gidx[dir]-1;

    // Compute the linear indices for both skin and ghost locations.
    long clinidx_skin = gkyl_range_idx(&conf_skin_r, sidx);
    long clinidx_ghost = gkyl_range_idx(&conf_ghost_r, gidx);

    const double *jacskin_c = (const double *) gkyl_array_cfetch(jac, clinidx_skin);
    const double *jacghost_c = (const double *) gkyl_array_cfetch(jac, clinidx_ghost);

    // Deflate the jacobian in the ghost cell and compute its reciprocal.
    double jacghost_surf_c[num_basis_conf_surf], jacghost_surf_inv_c[num_basis_conf_surf];
    kers->deflate_conf_ghost_op(jacghost_c, jacghost_surf_c);
    kers->conf_inv_op(jacghost_surf_c, jacghost_surf_inv_c);

    // Deflate the jacobian in the skin cell.
    double jacskin_surf_c[num_basis_conf_surf];
    kers->deflate_conf_skin_op(jacskin_c, jacskin_surf_c);

    long plinidx_ghost = gkyl_range_idx(&phase_ghost_r, gidx);
    double *jf_c = (double *) gkyl_array_fetch(jf, plinidx_ghost);

    // Deflate Jf in the ghost cell.
    double jf_surf_c[num_basis_phase_surf];
    kers->deflate_phase_ghost_op(jf_c, jf_surf_c);

    // Multiply the surface Jf by the surface 1/J_ghost and by the surface J_skin.
    kers->conf_phase_mul_op(jacghost_surf_inv_c, jf_surf_c, jf_surf_c);
    kers->conf_phase_mul_op(jacskin_surf_c, jf_surf_c, jf_surf_c);

    // Inflate Jf.
    for (int k=0; k<jf->ncomp; k++)
      jf_c[k] = 0.0;
    kers->inflate_phase_ghost_op(jf_surf_c, jf_c);
  }
}

// Function to launch the CUDA kernel that performs the ghost-to-skin value transfer on the GPU.
void
gkyl_rescale_ghost_jacf_advance_cu(const struct gkyl_rescale_ghost_jacf *up,
  const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf)
{
  // Only proceed if the skin range has a non-zero volume (i.e., there are skin cells to update).
  int nblocks = phase_ghost_r->nblocks, nthreads = phase_ghost_r->nthreads; // CUDA grid configuration.

  // Launch the CUDA kernel to advance the ghost-to-skin update.
  gkyl_rescale_ghost_jacf_advance_cu_ker<<<nblocks, nthreads>>>(up->kernels, up->dir, up->edge,
    *conf_skin_r, *conf_ghost_r, *phase_ghost_r, jac->on_dev, jf->on_dev);
}
