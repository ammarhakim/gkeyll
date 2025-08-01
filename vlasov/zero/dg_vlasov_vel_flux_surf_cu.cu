/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range, 
  const struct gkyl_array *jacob_vel, const struct gkyl_array *hamil, 
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM], idx_hamil[GKYL_MAX_DIM];
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

    const double *hamil_d = (const double*) gkyl_array_cfetch(hamil, hidx);
    const double *qmem_d = (const double*) gkyl_array_cfetch(qmem, cidx); 
    const double *pot_tot_d = (const double*) gkyl_array_cfetch(pot_tot, cidx); 
    const double *f_c = (const double*) gkyl_array_cfetch(fin, pidx); 
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, pidx);
    double *flux = (double*) gkyl_array_fetch(vel_flux_surf, pidx); 
    
    // Each cell owns *lower* fluxes in each velocity-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower velocity-space
    // edge, we call special kernels because we have *no* ghost cells in velocity
    // space, so we cannot index the distribution function in the lower direction
    // along that velocity-space edge. 
    for (int dir = 0; dir<vdim; ++dir) {
      if (idx[cdim+dir] == phase_range.lower[cdim+dir]) {
        cflrate_d[0] += up->vel_flux_surf_edge(up, dir, up->phase_grid.dx, 
          jacob_vel ? (const double*) gkyl_array_cfetch(jacob_vel, vidx) : 0,
          hamil_d, qmem_d, pot_tot_d, f_c, flux); 
      }
      else {
        gkyl_copy_int_arr(pdim, idx, idx_l);
        idx_l[cdim+dir] = idx_l[cdim+dir]-1;
        long pidx_l = gkyl_range_idx(&phase_range, idx_l); 
        const double* f_l = (const double*) gkyl_array_cfetch(fin, pidx_l);  
        cflrate_d[0] += up->vel_flux_surf(up, dir, up->phase_grid.dx, 
          jacob_vel ? (const double*) gkyl_array_cfetch(jacob_vel, vidx) : 0,
          hamil_d, qmem_d, pot_tot_d, f_l, f_c, flux);      
      }
    }    
  }
}

void 
gkyl_dg_vlasov_vel_flux_surf_advance_cu(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel, const struct gkyl_array *hamil, 
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_vlasov_vel_flux_surf_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, jacob_vel ? jacob_vel->on_dev : 0,
    hamil->on_dev, qmem->on_dev, pot_tot->on_dev, fin->on_dev, cflrate->on_dev, vel_flux_surf->on_dev);  
}

// CUDA kernel to set device pointers to canonical pb vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs(struct gkyl_dg_vlasov_vel_flux_surf *up,
  enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, bool has_qmem, bool has_phi)
{
  // By default, we have no forces from Hamiltonian, E, B, or phi. 
  for (int d=0; d<vdim; ++d) {
    up->hamil_alpha_quad[d] = no_hamil_alpha_quad; 
    up->EB_alpha_quad[d] = no_EB_alpha_quad;
    up->phi_alpha_quad[d] = no_phi_alpha_quad; 
  } 

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->lax_flux_nodal_to_modal[0] = ser_lax_flux_nodal_to_modal_vx_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[1] = ser_lax_flux_nodal_to_modal_vy_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[2] = ser_lax_flux_nodal_to_modal_vz_kernels[kernel_index].kernels[poly_order];

      // Only have Hamiltonian forces in general geometry. 
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        up->hamil_alpha_quad[0] = ser_hamil_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[1] = ser_hamil_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[2] = ser_hamil_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (has_qmem) {
        up->EB_alpha_quad[0] = ser_EB_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->EB_alpha_quad[1] = ser_EB_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->EB_alpha_quad[2] = ser_EB_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (has_phi) {
        up->phi_alpha_quad[0] = ser_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[1] = ser_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[2] = ser_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      up->lax_flux_nodal_to_modal[0] = tensor_lax_flux_nodal_to_modal_vx_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[1] = tensor_lax_flux_nodal_to_modal_vy_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[2] = tensor_lax_flux_nodal_to_modal_vz_kernels[kernel_index].kernels[poly_order];
      
      // Only have Hamiltonian forces in general geometry. 
      if (model_id == GKYL_MODEL_CANONICAL_PB || model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        up->hamil_alpha_quad[0] = tensor_hamil_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[1] = tensor_hamil_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[2] = tensor_hamil_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (has_qmem) {
        up->EB_alpha_quad[0] = tensor_EB_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->EB_alpha_quad[1] = tensor_EB_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->EB_alpha_quad[2] = tensor_EB_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (has_phi) {
        up->phi_alpha_quad[0] = tensor_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[1] = tensor_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[2] = tensor_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }
      
      break;      

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->vel_flux_surf = vel_flux_surf_kernels[kernel_index].kernels[poly_order];
  // Currently only support zero-flux boundary, so edge velocity flux an empty function (flux = 0.0). 
  up->vel_flux_surf_edge = no_vel_flux_surf_edge;   
}

gkyl_dg_vlasov_vel_flux_surf*
gkyl_dg_vlasov_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp)
{
  struct gkyl_dg_vlasov_vel_flux_surf *up = (struct gkyl_dg_vlasov_vel_flux_surf*) gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim; 
  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index. 
  up->hamil_range = *inp->hamil_range; 
  if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
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

  struct gkyl_dg_vlasov_vel_flux_surf *up_cu = (struct gkyl_dg_vlasov_vel_flux_surf*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_vlasov_vel_flux_surf), GKYL_CU_MEMCPY_H2D);

  gkyl_dg_vlasov_vel_flux_surf_set_cu_dev_ptrs<<<1,1>>>(up_cu, inp->conf_basis->b_type, 
    cdim, vdim, poly_order, inp->model_id, inp->has_qmem, inp->has_phi);  

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}