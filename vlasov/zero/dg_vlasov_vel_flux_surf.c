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
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    up->hamil_dim = pdim;
    up->hamil_offset = 0;
  }
  else {
    up->hamil_dim = vdim;
    up->hamil_offset = cdim;
  }
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE);
  // The velocity map is required: it provides the velocity-space Jacobian at
  // surface quadrature points and the velocity-space range used to index it.
  assert(inp->vel_map);
  up->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  up->vel_range = inp->vel_map->local_vel;
  // Borrowed pointer; kept alive by the acquired vel_map.
  up->jacob_vel_surf = inp->vel_map->jacob_vel_surf;
  // Position map: provides the (per-conf-cell constant) Jacobian used to
  // transform the -grad(phi) force to the mapped grid. Borrowed pointer.
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  // By default, we have no forces from Hamiltonian, E, B, phi, or radiation.
  for (int d=0; d<vdim; ++d) {
    up->hamil_alpha_quad[d] = no_hamil_alpha_quad; 
    up->E_alpha_quad[d] = no_E_alpha_quad;
    up->phi_alpha_quad[d] = no_phi_alpha_quad; 
    up->B_alpha_quad[d] = no_B_alpha_quad;
    up->rad_alpha_quad[d] = no_rad_alpha_quad; 
    up->hamil_alpha_quad_arr[d] = no_hamil_alpha_quad_arr; 
    up->E_alpha_quad_arr[d] = no_E_alpha_quad_arr;
    up->phi_alpha_quad_arr[d] = no_phi_alpha_quad_arr; 
    up->B_alpha_quad_arr[d] = no_B_alpha_quad_arr;
    up->rad_alpha_quad_arr[d] = no_rad_alpha_quad_arr; 
  } 

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( inp->use_lo ) {
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[0] = ser_lax_flux_nodal_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[1] = ser_lax_flux_nodal_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[2] = ser_lax_flux_nodal_vz_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      } 
      else {
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[0] = ser_ho_lax_flux_nodal_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[1] = ser_ho_lax_flux_nodal_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_arr[2] = ser_ho_lax_flux_nodal_vz_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry. 
      if (inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        if ( inp->use_lo ) {
          up->hamil_alpha_quad[0] = ser_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = ser_hamil_phase_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = ser_hamil_phase_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = ser_hamil_phase_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }
        else {
          up->hamil_alpha_quad[0] = ser_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = ser_hamil_phase_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = ser_hamil_phase_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = ser_hamil_phase_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (inp->model_id == GKYL_MODEL_TRIAD) {
        if ( inp->use_lo ) {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = hamil_sparse ?
            ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }
      }
      else if (inp->model_id == GKYL_MODEL_TRIAD_GR) {
        if ( inp->use_lo ) {
          up->hamil_alpha_quad[0] = ser_nc_hamil_phase_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = ser_nc_hamil_phase_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_nc_hamil_phase_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = ser_nc_hamil_phase_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_nc_hamil_phase_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = ser_nc_hamil_phase_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        } 
        else {
          up->hamil_alpha_quad[0] = ser_nc_hamil_phase_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[0] = ser_nc_hamil_phase_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[1] = ser_nc_hamil_phase_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[1] = ser_nc_hamil_phase_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad[2] = ser_nc_hamil_phase_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->hamil_alpha_quad_arr[2] = ser_nc_hamil_phase_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }
      }

      if ( inp->use_lo ) {
        if (inp->has_E) {
          up->E_alpha_quad[0] = ser_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[0] = ser_E_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = ser_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[1] = ser_E_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = ser_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[2] = ser_E_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_phi) {
          up->phi_alpha_quad[0] = ser_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[0] = ser_phi_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = ser_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[1] = ser_phi_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = ser_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[2] = ser_phi_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_B) {
          up->B_alpha_quad[0] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[0] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[1] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[2] = hamil_sparse ?
            ser_B_sparse_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_rad) {
          up->rad_alpha_quad[0] = ser_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[0] = ser_rad_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = ser_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[1] = ser_rad_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = ser_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[2] = ser_rad_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }    
      }
      else {
        if (inp->has_E) {
          up->E_alpha_quad[0] = ser_E_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[0] = ser_E_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[1] = ser_E_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[1] = ser_E_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad[2] = ser_E_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->E_alpha_quad_arr[2] = ser_E_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_phi) {
          up->phi_alpha_quad[0] = ser_phi_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[0] = ser_phi_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[1] = ser_phi_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[1] = ser_phi_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad[2] = ser_phi_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->phi_alpha_quad_arr[2] = ser_phi_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_B) {
          up->B_alpha_quad[0] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[0] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[1] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[1] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad[2] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->B_alpha_quad_arr[2] = hamil_sparse ?
            ser_B_sparse_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order] :
            ser_B_dense_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }

        if (inp->has_rad) {
          up->rad_alpha_quad[0] = ser_rad_ho_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[0] = ser_rad_ho_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[1] = ser_rad_ho_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[1] = ser_rad_ho_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad[2] = ser_rad_ho_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
          up->rad_alpha_quad_arr[2] = ser_rad_ho_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
        }    
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      
      up->lax_flux_nodal[0] = tensor_lax_flux_nodal_vx_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_arr[0] = tensor_lax_flux_nodal_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_vx_cfl_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal[1] = tensor_lax_flux_nodal_vy_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_arr[1] = tensor_lax_flux_nodal_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_vy_cfl_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_nodal[2] = tensor_lax_flux_nodal_vz_kernels[kernel_index].kernels[poly_order];
      up->lax_flux_arr[2] = tensor_lax_flux_nodal_vz_arr_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_vz_cfl_kernels[kernel_index].kernels[poly_order];


      // Only have Hamiltonian forces in general geometry. 
      if (inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        gkyl_exit("dg_vlasov_vel_flux_surf: Tensor basis and general Hamiltonian, GKYL_MODEL_CAN_PB or GKYL_MODEL_CANONICAL_PB_GR not yet supported!"); 
      }

      if (inp->has_E) {
        up->E_alpha_quad[0] = tensor_E_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad_arr[0] = tensor_E_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[1] = tensor_E_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad_arr[1] = tensor_E_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad[2] = tensor_E_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        up->E_alpha_quad_arr[2] = tensor_E_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_phi) {
        up->phi_alpha_quad[0] = tensor_phi_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad_arr[0] = tensor_phi_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[1] = tensor_phi_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad_arr[1] = tensor_phi_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad[2] = tensor_phi_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        up->phi_alpha_quad_arr[2] = tensor_phi_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_B) {
        up->B_alpha_quad[0] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vx_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad_arr[0] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[1] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vy_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad_arr[1] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad[2] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vz_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        up->B_alpha_quad_arr[2] = hamil_sparse ?
          tensor_B_sparse_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order] :
          tensor_B_dense_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
      }

      if (inp->has_rad) {
        up->rad_alpha_quad[0] = tensor_rad_alpha_quad_vx_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad_arr[0] = tensor_rad_alpha_quad_vx_arr_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[1] = tensor_rad_alpha_quad_vy_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad_arr[1] = tensor_rad_alpha_quad_vy_arr_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad[2] = tensor_rad_alpha_quad_vz_kernels[kernel_index].kernels[poly_order];
        up->rad_alpha_quad_arr[2] = tensor_rad_alpha_quad_vz_arr_kernels[kernel_index].kernels[poly_order];
      }      

      break;      

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->vel_flux_surf = vel_flux_surf_arrays;
  // Surface node counts for the per-node dispatch: (p+1) points per direction,
  // p+2 for the higher-order (anti-aliasing) and tensor (cubic-map) kernels.
  int nq = poly_order + 1;
  if ((poly_order > 1) && !inp->use_lo) nq = poly_order + 2;
  if (inp->conf_basis->b_type == GKYL_BASIS_MODAL_TENSOR) nq = poly_order + 2;
  up->num_nodes_conf = 1;
  for (int d=0; d<cdim; ++d) up->num_nodes_conf *= nq;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim-1; ++d) up->num_nodes_vel *= nq;
  // Currently only support zero-flux boundary, so edge velocity flux an empty function (flux = 0.0). 
  up->vel_flux_surf_edge = no_vel_flux_surf_edge;   

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) {
    assert(up->lax_flux_nodal[i]);
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
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil,
  const struct gkyl_array *qmem, const struct gkyl_array *pot_tot, const struct gkyl_array *rad,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vel_flux_surf)) {
    return gkyl_dg_vlasov_vel_flux_surf_advance_cu(up, conf_range, phase_range,
      poisson_tensor_conf, hamil, qmem, pot_tot, rad, fin, cflrate, vel_flux_surf);
  }
#endif
  const struct gkyl_array *jacob_vel_surf = up->jacob_vel_surf;
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
        // Velocity-space Jacobian of the lower neighbor in this direction, for
        // the minimum-Jacobian time-step estimate of the C^0 linear map.
        int idx_vel_l[GKYL_MAX_DIM];
        for (int i=0; i<vdim; ++i) {
          idx_vel_l[i] = idx_l[cdim+i];
        }
        long vidx_l = gkyl_range_idx(&up->vel_range, idx_vel_l);
        cflrate_d[0] += up->vel_flux_surf(up, dir, xcC, up->phase_grid.dx,
          up->jacob_pos ? gkyl_array_cfetch(up->jacob_pos, cidx) : 0,
          gkyl_array_cfetch(jacob_vel_surf, vidx_l),
          gkyl_array_cfetch(jacob_vel_surf, vidx), poisson_tensor_conf_d,
          hamil_d, qmem_d, pot_tot_d, rad_d, f_l, f_c, flux);
      }
    }
  }
}

void
gkyl_dg_vlasov_vel_flux_surf_release(struct gkyl_dg_vlasov_vel_flux_surf* up)
{
  // Release memory associated with this updater.
  gkyl_vlasov_velocity_map_release(up->vel_map);
  gkyl_vlasov_position_map_release(up->pos_map);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  gkyl_free(up);
}
