#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_vlasov *vlasov = container_of(base->on_dev, struct dg_vlasov, eqn);
    gkyl_cu_free(vlasov);
  }
  
  struct dg_vlasov *vlasov = container_of(base, struct dg_vlasov, eqn);
  gkyl_array_release(vlasov->poisson_tensor_conf);
  gkyl_array_release(vlasov->hamil);
  gkyl_array_release(vlasov->qmem);
  gkyl_array_release(vlasov->pot_tot);
  gkyl_array_release(vlasov->rad);
  if (vlasov->vel_map) {
    gkyl_vlasov_velocity_map_release(vlasov->vel_map);
  }
  if (vlasov->pos_map) {
    gkyl_vlasov_position_map_release(vlasov->pos_map);
  }
  if (vlasov->use_conf_flux_surf) {
    gkyl_array_release(vlasov->conf_flux_surf);
  }
  gkyl_array_release(vlasov->vel_flux_surf); 
  gkyl_array_release(vlasov->f_no_J); 
  gkyl_free(vlasov);
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_inew(const struct gkyl_dg_vlasov_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if(inp->use_gpu) {
    return gkyl_dg_vlasov_cu_dev_inew(inp);
  } 
#endif
  struct dg_vlasov *vlasov = gkyl_malloc(sizeof(*vlasov));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = inp->conf_basis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;  

  // Are we skipping cells with small phase space density?
  if (inp->skip_cell_thresh > 0.0) {
    vlasov->skip_cell_thresh = inp->skip_cell_thresh * pow(sqrt(2.0), pdim);
  }
  else {
    vlasov->skip_cell_thresh = -1.0;
  }

  // Determine Hamiltonian dimensionality and index offset for indexing Hamiltonian
  // from an input phase space index.
  if (inp->hamil_id == GKYL_HAMIL_PHASE) {
    vlasov->hamil_dim = pdim;
    vlasov->hamil_offset = 0;
  }
  else {
    vlasov->hamil_dim = vdim;
    vlasov->hamil_offset = cdim;
  }
  vlasov->hamil_range = *inp->hamil_range;
  vlasov->conf_range = *inp->conf_range;
  vlasov->phase_range = *inp->phase_range;

  // The velocity map is required: it provides the velocity-space range used
  // to index per-velocity-cell quantities (Jacobian, radiation drag).
  assert(inp->vel_map);
  vlasov->vel_map = gkyl_vlasov_velocity_map_acquire(inp->vel_map);
  vlasov->vel_range = inp->vel_map->local_vel;
  vlasov->use_vmap = inp->vel_map->is_mapped;
  vlasov->jacob_vel = 0;
  if (vlasov->use_vmap) {
    // Borrowed pointer; kept alive by the acquired vel_map.
    vlasov->jacob_vel = inp->vel_map->jacob_vel;
  }
  // The position map is required: it provides the (per-conf-cell constant)
  // configuration-space Jacobian used to transform the streaming term.
  assert(inp->pos_map);
  vlasov->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  vlasov->jacob_pos = inp->pos_map->jacob_pos; // Borrowed; kept alive by the acquired pos_map.
  vlasov->poisson_tensor_conf = gkyl_array_acquire(inp->poisson_tensor_conf); 
  vlasov->hamil = gkyl_array_acquire(inp->hamil); 
  vlasov->qmem = gkyl_array_acquire(inp->qmem); 
  vlasov->pot_tot = gkyl_array_acquire(inp->pot_tot); 
  vlasov->rad = gkyl_array_acquire(inp->rad);
  vlasov->use_conf_flux_surf = false;
  if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR){
    vlasov->use_conf_flux_surf = true;
    vlasov->conf_flux_surf = gkyl_array_acquire(inp->conf_flux_surf);
  }
  vlasov->vel_flux_surf = gkyl_array_acquire(inp->vel_flux_surf); 
  vlasov->f_no_J = gkyl_array_acquire(inp->f_no_J); 

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = vlasov_vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;
  vlasov->eqn.boundary_diag_term = boundary_diag;

  // By default, we have no forces from E, B, phi, or radiation. 
  vlasov->E_vol = no_E_vol; 
  vlasov->Bx_vol = no_B_vol; 
  vlasov->By_vol = no_B_vol; 
  vlasov->Bz_vol = no_B_vol; 
  vlasov->phi_vol = no_phi_vol; 
  vlasov->rad_vol = no_rad_vol; 

  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, 
    *stream_surf_y_kernels, 
    *stream_surf_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels, 
    *stream_boundary_surf_y_kernels,
    *stream_boundary_surf_z_kernels;

  const gkyl_dg_vlasov_stream_surf_from_flux_kern_list *stream_surf_from_flux_x_kernels, 
    *stream_surf_from_flux_y_kernels, 
    *stream_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list *stream_boundary_surf_from_flux_x_kernels, 
    *stream_boundary_surf_from_flux_y_kernels,
    *stream_boundary_surf_from_flux_z_kernels;

  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, 
    *accel_surf_vy_kernels, 
    *accel_surf_vz_kernels;

  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, 
    *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  int kernel_index = cv_index[cdim].vdim[vdim]; 
  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Set function pointers for individual pieces of the volume update.
      if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
        if (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = ser_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol = ser_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = ser_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = ser_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = ser_stream_hamil_vel_sparse_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_vel_sparse_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_vel_sparse_surf_z_kernels;
          stream_boundary_surf_x_kernels = ser_stream_hamil_vel_sparse_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_vel_sparse_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_vel_sparse_boundary_surf_z_kernels;
        }
        else {
          vlasov->hamil_vol = ser_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol = ser_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = ser_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = ser_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = ser_stream_hamil_vel_dense_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_vel_dense_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_vel_dense_surf_z_kernels;
          stream_boundary_surf_x_kernels = ser_stream_hamil_vel_dense_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_vel_dense_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_vel_dense_boundary_surf_z_kernels;
        }
        if (inp->has_rad) vlasov->rad_vol = ser_rad_vol_kernels[kernel_index].kernels[poly_order];

      }
      else if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {

        if (inp->model_id == GKYL_MODEL_TRIAD) {
          vlasov->hamil_vol = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) ?
            ser_nc_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order] :
            ser_nc_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
        }
        else if (inp->model_id == GKYL_MODEL_TRIAD_GR) {
          vlasov->hamil_vol = ser_nc_hamil_phase_vol_kernels[kernel_index].kernels[poly_order];
        }

        if ( inp->use_lo ) {
          stream_surf_from_flux_x_kernels = ser_stream_nc_hamil_gen_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_nc_hamil_gen_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_nc_hamil_gen_surf_z_kernels;
          
          stream_boundary_surf_from_flux_x_kernels = ser_stream_nc_hamil_gen_boundary_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_nc_hamil_gen_boundary_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_nc_hamil_gen_boundary_surf_z_kernels; 
        } 
        else {
          stream_surf_from_flux_x_kernels = ser_stream_nc_hamil_gen_ho_surf_x_kernels;
          stream_surf_from_flux_y_kernels = ser_stream_nc_hamil_gen_ho_surf_y_kernels;
          stream_surf_from_flux_z_kernels = ser_stream_nc_hamil_gen_ho_surf_z_kernels;
          
          stream_boundary_surf_from_flux_x_kernels = ser_stream_nc_hamil_gen_boundary_ho_surf_x_kernels;
          stream_boundary_surf_from_flux_y_kernels = ser_stream_nc_hamil_gen_boundary_ho_surf_y_kernels;
          stream_boundary_surf_from_flux_z_kernels = ser_stream_nc_hamil_gen_boundary_ho_surf_z_kernels; 
        }
      }
      else {
        vlasov->hamil_vol = ser_hamil_gen_vol_kernels[kernel_index].kernels[poly_order];

        if ( inp->use_lo ) {
          stream_surf_x_kernels = ser_stream_hamil_gen_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_gen_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_gen_surf_z_kernels;
          
          stream_boundary_surf_x_kernels = ser_stream_hamil_gen_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_gen_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_gen_boundary_surf_z_kernels;    
        } 
        else {
          stream_surf_x_kernels = ser_stream_hamil_gen_ho_surf_x_kernels;
          stream_surf_y_kernels = ser_stream_hamil_gen_ho_surf_y_kernels;
          stream_surf_z_kernels = ser_stream_hamil_gen_ho_surf_z_kernels;
          
          stream_boundary_surf_x_kernels = ser_stream_hamil_gen_boundary_ho_surf_x_kernels;
          stream_boundary_surf_y_kernels = ser_stream_hamil_gen_boundary_ho_surf_y_kernels;
          stream_boundary_surf_z_kernels = ser_stream_hamil_gen_boundary_ho_surf_z_kernels;    
        }     
      }
      if (inp->has_E) vlasov->E_vol = ser_E_vol_kernels[kernel_index].kernels[poly_order];
      if (inp->has_phi) vlasov->phi_vol = ser_phi_vol_kernels[kernel_index].kernels[poly_order];


      if ( inp->use_lo ) {
        accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
        accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
        accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
        accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      }
      else {
        accel_surf_vx_kernels = ser_accel_ho_surf_vx_kernels;
        accel_surf_vy_kernels = ser_accel_ho_surf_vy_kernels;
        accel_surf_vz_kernels = ser_accel_ho_surf_vz_kernels;
        accel_boundary_surf_vx_kernels = ser_accel_boundary_ho_surf_vx_kernels;
        accel_boundary_surf_vy_kernels = ser_accel_boundary_ho_surf_vy_kernels;
        accel_boundary_surf_vz_kernels = ser_accel_boundary_ho_surf_vz_kernels;
      }
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      // Set function pointers for individual pieces of the volume update.
      if (inp->model_id == GKYL_MODEL_DEFAULT || inp->model_id == GKYL_MODEL_SR) {
        if (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE) {
          vlasov->hamil_vol = tensor_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol = tensor_Bx_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = tensor_By_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = tensor_Bz_hamil_vel_sparse_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_sparse_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_sparse_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_sparse_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_sparse_boundary_surf_z_kernels;
        }
        else {
          vlasov->hamil_vol = tensor_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];

          if (inp->has_B) {
            vlasov->Bx_vol = tensor_Bx_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->By_vol = tensor_By_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
            vlasov->Bz_vol = tensor_Bz_hamil_vel_dense_vol_kernels[kernel_index].kernels[poly_order];
          }

          stream_surf_x_kernels = tensor_stream_hamil_vel_dense_surf_x_kernels;
          stream_surf_y_kernels = tensor_stream_hamil_vel_dense_surf_y_kernels;
          stream_surf_z_kernels = tensor_stream_hamil_vel_dense_surf_z_kernels;

          stream_boundary_surf_x_kernels = tensor_stream_hamil_vel_dense_boundary_surf_x_kernels;
          stream_boundary_surf_y_kernels = tensor_stream_hamil_vel_dense_boundary_surf_y_kernels;
          stream_boundary_surf_z_kernels = tensor_stream_hamil_vel_dense_boundary_surf_z_kernels;
        }
        if (inp->has_rad) vlasov->rad_vol = tensor_rad_vol_kernels[kernel_index].kernels[poly_order];
      }
      else if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {
        gkyl_exit("dg_vlasov: Tensor basis and general Hamiltonian, GKYL_MODEL_TRIAD not yet supported!"); 
      }
      else {
        gkyl_exit("dg_vlasov: Tensor basis and general Hamiltonian, GKYL_MODEL_CAN_PB or GKYL_MODEL_CANONICAL_PB_GR not yet supported!"); 
      }
      if (inp->has_E) vlasov->E_vol = tensor_E_vol_kernels[kernel_index].kernels[poly_order];
      if (inp->has_phi) vlasov->phi_vol = tensor_phi_vol_kernels[kernel_index].kernels[poly_order];

      accel_surf_vx_kernels = tensor_accel_surf_vx_kernels;
      accel_surf_vy_kernels = tensor_accel_surf_vy_kernels;
      accel_surf_vz_kernels = tensor_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = tensor_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = tensor_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = tensor_accel_boundary_surf_vz_kernels;
      break;      

    default:
      assert(false);
      break;    
  }


  if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {
    vlasov->stream_surf_from_flux[0] = stream_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[1] = stream_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf_from_flux[2] = stream_surf_from_flux_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf_from_flux[0] = stream_boundary_surf_from_flux_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[1] = stream_boundary_surf_from_flux_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf_from_flux[2] = stream_boundary_surf_from_flux_z_kernels[kernel_index].kernels[poly_order];    
  } 
  else {
    vlasov->stream_surf[0] = stream_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[1] = stream_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_surf[2] = stream_surf_z_kernels[kernel_index].kernels[poly_order];

    vlasov->stream_boundary_surf[0] = stream_boundary_surf_x_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[1] = stream_boundary_surf_y_kernels[kernel_index].kernels[poly_order];
    vlasov->stream_boundary_surf[2] = stream_boundary_surf_z_kernels[kernel_index].kernels[poly_order]; 
  }   

  vlasov->accel_surf[0] = accel_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[1] = accel_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_surf[2] = accel_surf_vz_kernels[kernel_index].kernels[poly_order];

  vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[kernel_index].kernels[poly_order];
  vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[kernel_index].kernels[poly_order]; 

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) {
    if (inp->model_id == GKYL_MODEL_TRIAD || inp->model_id == GKYL_MODEL_TRIAD_GR) {
      assert(vlasov->stream_surf_from_flux[i]);
      assert(vlasov->stream_boundary_surf_from_flux[i]);
    } else {
      assert(vlasov->stream_surf[i]);
      assert(vlasov->stream_boundary_surf[i]);
    }
  }
  for (int i=0; i<vdim; ++i) {
    assert(vlasov->accel_surf[i]);
    assert(vlasov->accel_boundary_surf[i]);
  }
  
  vlasov->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vlasov->eqn.flags);

  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);
  vlasov->eqn.on_dev = &vlasov->eqn; // CPU eqn obj points to itself
  
  return &vlasov->eqn;
}
