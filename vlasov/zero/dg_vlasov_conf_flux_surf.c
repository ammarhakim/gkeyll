#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_conf_flux_surf_priv.h> 
#include <gkyl_util.h> 

gkyl_dg_vlasov_conf_flux_surf*
gkyl_dg_vlasov_conf_flux_surf_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if(inp->use_gpu) {
    return gkyl_dg_vlasov_conf_flux_surf_cu_dev_inew(inp);
  } 
#endif  
  struct gkyl_dg_vlasov_conf_flux_surf *up = gkyl_malloc(sizeof(*up));

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
  if (inp->model_id == GKYL_MODEL_TRIAD) {
    up->hamil_dim = vdim; 
    up->hamil_offset = cdim; 
  }
  else if (inp->model_id == GKYL_MODEL_TRIAD_GR
    || inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
    // Full phase-space Hamiltonian.
    up->hamil_dim = pdim;
    up->hamil_offset = 0;
  }
  else {
    assert(false); // Should not be here for other models
  }
  up->vel_range = *inp->vel_range; 
  // Sparse (separable) vs. dense velocity-space Hamiltonian kernel selection.
  bool hamil_sparse = (inp->hamil_id == GKYL_HAMIL_VEL_SPARSE);

  int kernel_index = cv_index[cdim].vdim[vdim];   
  switch (inp->conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( inp->use_lo ) {
        up->lax_flux_nodal[0] = ser_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } 
      else {
        up->lax_flux_nodal[0] = ser_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = ser_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = ser_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = ser_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = ser_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = ser_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      }

      // Only have Hamiltonian forces in general geometry.
      if (inp->model_id == GKYL_MODEL_TRIAD) {
        if ( inp->use_lo ) {
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
      else if (inp->model_id == GKYL_MODEL_TRIAD_GR
        || inp->model_id == GKYL_MODEL_CANONICAL_PB || inp->model_id == GKYL_MODEL_CANONICAL_PB_GR) {
        // Full phase-space Hamiltonian: alpha_dir = P . grad_v H evaluated at
        // the surface nodes. Canonical-PB models supply the identity Poisson
        // tensor, reducing this to the canonical streaming speed dH/dv_dir.
        if ( inp->use_lo ) {
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
      if ( inp->use_lo ) {
        up->lax_flux_nodal[0] = tensor_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } 
      else {
        up->lax_flux_nodal[0] = tensor_ho_lax_flux_nodal_x_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[0] = tensor_ho_lax_flux_nodal_x_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[1] = tensor_ho_lax_flux_nodal_y_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[1] = tensor_ho_lax_flux_nodal_y_cfl_kernels[kernel_index].kernels[poly_order];
        up->lax_flux_nodal[2] = tensor_ho_lax_flux_nodal_z_kernels[kernel_index].kernels[poly_order];
        up->lax_cfl[2] = tensor_ho_lax_flux_nodal_z_cfl_kernels[kernel_index].kernels[poly_order];
      } 

      break;      

    default:
      assert(false);
      break;    
  } 
  // Set assembly functions for computing fluxes. 
  up->conf_flux_surf = conf_flux_surf_nodes;
  // Surface node counts for the per-node dispatch: (p+1) points per direction,
  // p+2 for the higher-order (anti-aliasing) kernels.
  int nq = poly_order + 1;
  if ((poly_order > 1) && !inp->use_lo) nq = poly_order + 2;
  up->num_nodes_conf = 1;
  for (int d=0; d<cdim-1; ++d) up->num_nodes_conf *= nq;
  up->num_nodes_vel = 1;
  for (int d=0; d<vdim; ++d) up->num_nodes_vel *= nq;

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) {
    assert(up->lax_flux_nodal[i]);
    assert(up->hamil_alpha_quad[i]);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;  
}

void gkyl_dg_vlasov_conf_flux_surf_advance(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(conf_flux_surf)) {
    return gkyl_dg_vlasov_conf_flux_surf_advance_cu(up, conf_range, phase_range, phase_range_ext, 
      poisson_tensor_conf, hamil, fin, cflrate, conf_flux_surf);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int vdim = pdim - cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM], idx_vel[GKYL_MAX_DIM]; 
  int idx_hamil[GKYL_MAX_DIM];
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
    double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xcC);

    const double *f_c = gkyl_array_cfetch(fin, pidx); 
    double *cflrate_d = gkyl_array_fetch(cflrate, pidx);
    const double *hamil_d = gkyl_array_cfetch(hamil, hidx); 
    const double *poisson_tensor_conf_d = gkyl_array_cfetch(poisson_tensor_conf, cidx); 
    double *flux = gkyl_array_fetch(conf_flux_surf, pidx); 

    // Each cell owns *lower* fluxes in each configuration-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower configuration-space
    // edge, we call ghost cells
    for (int dir = 0; dir<cdim; ++dir) {

      // Create an index for the left cell (which may be a ghost cell) 
      gkyl_copy_int_arr(pdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long pidx_l = gkyl_range_idx(phase_range, idx_l); 
      const double *f_l = gkyl_array_cfetch(fin, pidx_l);

      // Which face to evalute the hamiltonian / pt on
      int hamil_pt_edge = -1;

      cflrate_d[0] += up->conf_flux_surf(up, dir, xcC, up->phase_grid.dx, hamil_pt_edge,
        poisson_tensor_conf_d, hamil_d, f_l, f_c, flux);     

      // If at the right boundary compute flux owned by the point in the ghost cell
      if (idx[dir] == phase_range->upper[dir]) {

        // Index the right cell (ghost cell)
        gkyl_copy_int_arr(pdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir]+1;
        long pidx_r = gkyl_range_idx(phase_range_ext, idx_r); 

        const double *f_r = gkyl_array_cfetch(fin, pidx_r);
        double *flux_r = gkyl_array_fetch(conf_flux_surf, pidx_r); 
        double *cflrate_d_r = gkyl_array_fetch(cflrate, pidx_r);

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
gkyl_dg_vlasov_conf_flux_surf_release(struct gkyl_dg_vlasov_conf_flux_surf* up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif
  gkyl_free(up);
}

