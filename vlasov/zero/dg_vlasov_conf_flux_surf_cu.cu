/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_conf_flux_surf_priv.h>
#include <gkyl_util.h> 
}

#include <cassert>

__global__ void
gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range, struct gkyl_range phase_range_ext, 
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM]; 
  int idx_hamil[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);  
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    } 
    long vidx = gkyl_range_idx(&up->vel_range, idx_vel); 

    for (int i=0; i<up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset+i];
    } 
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil); 

    // Grab the cell center location for NC bracket calculation 
    double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    const double *f_c = (const double*) gkyl_array_cfetch(fin, pidx); 
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, pidx);
    const double *poisson_tensor_conf_d = (const double*) gkyl_array_cfetch(poisson_tensor_conf, cidx);
    const double *hamil_d = (const double*) gkyl_array_cfetch(hamil, hidx);
    double *flux = (double*) gkyl_array_fetch(conf_flux_surf, pidx); 
    
    // Each cell owns *lower* fluxes in each configuration-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower configuration-space
    // edge, we call ghost cells
    for (int dir = 0; dir<cdim; ++dir) {

      // Create an index for the left cell (which may be a ghost cell)
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long pidx_l = gkyl_range_idx(&phase_range, idx_l); 
      const double* f_l = (const double*) gkyl_array_cfetch(fin, pidx_l);

      // Which face to evalute the hamiltonian / pt on
      int hamil_pt_edge = -1;

      cflrate_d[0] += up->conf_flux_surf(up, dir, xcC, up->phase_grid.dx, 
        hamil_pt_edge, poisson_tensor_conf_d, hamil_d, f_l, f_c, flux);    

      // If at the right boundary compute flux owned by the point in the ghost cell
      if (idx[dir] == phase_range.upper[dir]) {
        
        // Index the right cell (ghost cell)
        gkyl_copy_int_arr(pdim, idx, idx_r);
        idx_r[dir] = idx_r[dir]+1;
        long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r); 
        
        const double *f_r = (const double*) gkyl_array_cfetch(fin, pidx_r);
        double *flux_r = (double*) gkyl_array_fetch(conf_flux_surf, pidx_r); 
        double *cflrate_d_r = (double*) gkyl_array_fetch(cflrate, pidx_r);

        /* As a concequence of not having ghost cells for PT/Hamil, they are shifted here
          and evaluated in the kernels at the upper boundary +1. This is allowed by continuity of hamil/pt */
        hamil_pt_edge = 1;

        gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xcR);
        cflrate_d_r[0] += up->conf_flux_surf(up, dir, xcR, up->phase_grid.dx, hamil_pt_edge,
          poisson_tensor_conf_d, hamil_d, f_c, f_r, flux_r); 
      }
    }    
  }
}

void 
gkyl_dg_vlasov_conf_flux_surf_advance_cu(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, 
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_vlasov_conf_flux_surf_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, *phase_range_ext, poisson_tensor_conf->on_dev,
    hamil->on_dev, fin->on_dev, cflrate->on_dev, conf_flux_surf->on_dev);  
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs(struct gkyl_dg_vlasov_conf_flux_surf *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_hamil_id hamil_id, bool use_lo)
{
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (hamil_id == GKYL_HAMIL_VEL_SPARSE);

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( use_lo ) {
        up->lax_flux_nodal_to_modal[0] = ser_lax_flux_nodal_to_modal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[1] = ser_lax_flux_nodal_to_modal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[2] = ser_lax_flux_nodal_to_modal_z_kernels[kernel_index].kernels[poly_order];
      }
      else {
        up->lax_flux_nodal_to_modal[0] = ser_ho_lax_flux_nodal_to_modal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[1] = ser_ho_lax_flux_nodal_to_modal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[2] = ser_ho_lax_flux_nodal_to_modal_z_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry. 
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Only triads are supported currently for configuration space nodal fluxes.
        assert(false);
      }
      else if (model_id == GKYL_MODEL_TRIAD) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_hamil_vel_sparse_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_hamil_vel_sparse_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order] :
            ser_hamil_vel_dense_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (model_id == GKYL_MODEL_TRIAD_GR) {
        if ( use_lo ) {
          up->hamil_alpha_quad[0] = ser_hamil_phase_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = ser_hamil_phase_ho_alpha_quad_x_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_ho_alpha_quad_y_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_ho_alpha_quad_z_kernels[kernel_index].kernels[poly_order];
        }
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      if (use_lo) {
        up->lax_flux_nodal_to_modal[0] = tensor_lax_flux_nodal_to_modal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[1] = tensor_lax_flux_nodal_to_modal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[2] = tensor_lax_flux_nodal_to_modal_z_kernels[kernel_index].kernels[poly_order];
      } 
      else {
        up->lax_flux_nodal_to_modal[0] = tensor_ho_lax_flux_nodal_to_modal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[1] = tensor_ho_lax_flux_nodal_to_modal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal_to_modal[2] = tensor_ho_lax_flux_nodal_to_modal_z_kernels[kernel_index].kernels[poly_order];
      }
      
      // Only have Hamiltonian forces in general geometry. 
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        assert(false);
      }  

      break;      

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->conf_flux_surf = conf_flux_surf_kernels[kernel_index].kernels[poly_order];
}

gkyl_dg_vlasov_conf_flux_surf*
gkyl_dg_vlasov_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_conf_flux_surf *up = (struct gkyl_dg_vlasov_conf_flux_surf*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = true; 

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    up->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  }
  else {
    up->skip_cell_thresh = -1.0;
  }  

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  up->hamil_range = *inp->hamil_range; 
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR || inp->model_id == GKYL_MODEL_TRIAD) {
    up->hamil_dim = vdim; 
    up->hamil_offset = cdim; 
  }
  else {
    up->hamil_dim = pdim; 
    up->hamil_offset = 0; 
  }
  up->vel_range = *inp->vel_range; 
  
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_vlasov_conf_flux_surf *up_cu = (struct gkyl_dg_vlasov_conf_flux_surf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_conf_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_conf_flux_surf_set_cu_dev_ptrs<<<1,1>>>(up_cu, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->hamil_id, inp->use_lo);  

  // set parent on_dev pointer
  up->on_dev = up_cu;
 
  return up;
}
