#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf_priv.h>
#include <gkyl_util.h>

gkyl_dg_vlasov_vel_flux_surf*
gkyl_dg_vlasov_vel_flux_surf_inew(const struct gkyl_dg_vlasov_vel_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if(inp->use_gpu) {
    return gkyl_dg_vlasov_vel_flux_surf_cu_dev_inew(inp);
  } 
#endif  
  struct gkyl_dg_vlasov_vel_flux_surf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu; 

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

  // By default, we have no forces from Hamiltonian, E, B, phi, or radiation. 
  for (int d=0; d<vdim; ++d) {
    up->hamil_alpha_quad[d] = no_hamil_alpha_quad; 
    up->E_alpha_quad[d] = no_E_alpha_quad;
    up->phi_alpha_quad[d] = no_phi_alpha_quad; 
    up->B_alpha_quad[d] = no_B_alpha_quad;
    up->rad_alpha_quad[d] = no_rad_alpha_quad; 
  } 

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      up->lax_flux_nodal_to_modal[0] = ser_lax_flux_nodal_to_modal_vx_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[1] = ser_lax_flux_nodal_to_modal_vy_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[2] = ser_lax_flux_nodal_to_modal_vz_kernels[kernel_index].kernels[poly_order];

      // Only have Hamiltonian forces in general geometry. 
      if (inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        up->hamil_alpha_quad[0] = ser_hamil_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[1] = ser_hamil_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[2] = ser_hamil_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }
      else if (inp->model_id == GKYL_MODEL_TRIAD) {
        up->hamil_alpha_quad[0] = ser_nc_hamil_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[1] = ser_nc_hamil_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->hamil_alpha_quad[2] = ser_nc_hamil_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_E) {
        up->E_alpha_quad[0] = ser_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[1] = ser_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[2] = ser_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_phi) {
        up->phi_alpha_quad[0] = ser_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[1] = ser_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[2] = ser_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_B) {
        up->B_alpha_quad[0] = ser_B_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[1] = ser_B_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[2] = ser_B_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_rad) {
        up->rad_alpha_quad[0] = ser_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[1] = ser_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[2] = ser_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }      

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      up->lax_flux_nodal_to_modal[0] = tensor_lax_flux_nodal_to_modal_vx_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[1] = tensor_lax_flux_nodal_to_modal_vy_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal_to_modal[2] = tensor_lax_flux_nodal_to_modal_vz_kernels[kernel_index].kernels[poly_order];
      
      // Only have Hamiltonian forces in general geometry. 
      if (inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        gkyl_exit("dg_vlasov_vel_flux_surf: Tensor basis and general Hamiltonian, GKYL_MODEL_CAN_PB or GKYL_MODEL_CANONICAL_PB_GR not yet supported!"); 
      }

      if (inp->has_E) {
        up->E_alpha_quad[0] = tensor_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[1] = tensor_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[2] = tensor_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_phi) {
        up->phi_alpha_quad[0] = tensor_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[1] = tensor_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[2] = tensor_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_B) {
        up->B_alpha_quad[0] = tensor_B_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[1] = tensor_B_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[2] = tensor_B_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_rad) {
        up->rad_alpha_quad[0] = tensor_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[1] = tensor_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[2] = tensor_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
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

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) {
    assert(up->lax_flux_nodal_to_modal[i]);
    assert(up->hamil_alpha_quad[i]);
    assert(up->E_alpha_quad[i]);
    assert(up->phi_alpha_quad[i]);
    assert(up->B_alpha_quad[i]);
    assert(up->rad_alpha_quad[i]);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;  
}

void gkyl_dg_vlasov_vel_flux_surf_advance(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_surf, const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, const struct gkyl_array *rad, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vel_flux_surf)) {
    return gkyl_dg_vlasov_vel_flux_surf_advance_cu(up, conf_range, phase_range, 
      jacob_vel_surf, poisson_tensor_conf, hamil, qmem, pot_tot, rad, fin, cflrate, vel_flux_surf);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM], idx_hamil[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long cidx = gkyl_range_idx(conf_range, idx);
    long pidx = gkyl_range_idx(phase_range, idx);

    for (int i=0; i<vdim; ++i) {
      idx_vel[i] = idx[cdim+i];
    } 
    long vidx = gkyl_range_idx(&up->vel_range, idx_vel); 

    for (int i=0; i<up->hamil_dim; ++i) {
      idx_hamil[i] = idx[up->hamil_offset+i];
    } 
    long hidx = gkyl_range_idx(&up->hamil_range, idx_hamil); 

    // Grab the cell center location for NC bracket calculation 
    double xcC[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    const double *poisson_tensor_conf_d = gkyl_array_cfetch(poisson_tensor_conf, cidx);
    const double *hamil_d = gkyl_array_cfetch(hamil, hidx);
    const double *qmem_d = gkyl_array_cfetch(qmem, cidx); 
    const double *pot_tot_d = gkyl_array_cfetch(pot_tot, cidx); 
    const double *rad_d = gkyl_array_cfetch(rad, vidx); 
    const double *f_c = gkyl_array_cfetch(fin, pidx); 
    double *cflrate_d = gkyl_array_fetch(cflrate, pidx);
    double *flux = gkyl_array_fetch(vel_flux_surf, pidx); 

    // Each cell owns *lower* fluxes in each velocity-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower velocity-space
    // edge, we call special kernels because we have *no* ghost cells in velocity
    // space, so we cannot index the distribution function in the lower direction
    // along that velocity-space edge. 
    for (int dir = 0; dir<vdim; ++dir) {
      if (idx[cdim+dir] == phase_range->lower[cdim+dir]) {
        cflrate_d[0] += up->vel_flux_surf_edge(up, dir, xcC, up->phase_grid.dx, 
          jacob_vel_surf ? gkyl_array_cfetch(jacob_vel_surf, vidx) : 0, poisson_tensor_conf_d,
          hamil_d, qmem_d, pot_tot_d, rad_d, f_c, flux); 
      }
      else {
        gkyl_copy_int_arr(pdim, iter.idx, idx_l);
        idx_l[cdim+dir] = idx_l[cdim+dir]-1;
        long pidx_l = gkyl_range_idx(phase_range, idx_l); 
        const double *f_l = gkyl_array_cfetch(fin, pidx_l);  
        cflrate_d[0] += up->vel_flux_surf(up, dir, xcC, up->phase_grid.dx, 
          jacob_vel_surf ? gkyl_array_cfetch(jacob_vel_surf, vidx) : 0, poisson_tensor_conf_d,
          hamil_d, qmem_d, pot_tot_d, rad_d, f_l, f_c, flux);      
      }
    }
  }
}

void
gkyl_dg_vlasov_vel_flux_surf_release(struct gkyl_dg_vlasov_vel_flux_surf* up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  gkyl_free(up);
}
