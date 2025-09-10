/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_calc_gyrokinetic_vars_priv.h>
#include <gkyl_util.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads and can parallelize over components too
  dimBlock->y = ncomp;
  dimGrid->y = 1;
  dimBlock->x = gkyl_int_div_up(252, ncomp);
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ void
gkyl_dg_calc_gyrokinetic_vars_alpha_surf_cu_kernel(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range phase_ext_range, const struct gkyl_array *phi, 
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{ 
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM], idx_vel[2];
  double xc[GKYL_MAX_DIM];

  // 2D thread grid
  // linc2 = c where c is the component index (from 0 to cdim + 1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);
    const double *jacobtot_inv_d = (const double*) gkyl_array_cfetch(up->gk_geom->jacobtot_inv, loc_conf);
    const double *cmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->cmag, loc_conf);
    const double *b_i_d = (const double*) gkyl_array_cfetch(up->gk_geom->b_i, loc_conf);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double* alpha_surf_d = (double*) gkyl_array_fetch(alpha_surf, loc_phase);
    double* sgn_alpha_surf_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc_phase);
    int* const_sgn_alpha_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc_phase);

    for (int dir = 0; dir<cdim+1; ++dir) {
      // Each thread in linc2 thread grid handles a different component
      if (linc2 == dir) {
        const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->phase_grid.dx, 
          vmap_d, vmapSq_d, up->charge, up->mass, bmag_d, jacobtot_inv_d,
          cmag_d, b_i_d, phi_d, alpha_surf_d, sgn_alpha_surf_d);

        // If the phase space index is at the local configuration space upper value, we
        // we are at the configuration space upper edge and we also need to evaluate 
        // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
        // where it is not defined when computing the final surface alpha we need
        // (since the surface alpha array stores only the *lower* surface expansion)
        if (dir < cdim && idx[dir] == conf_range.upper[dir]) {
          gkyl_copy_int_arr(pdim, idx, idx_edge);
          idx_edge[dir] = idx_edge[dir]+1;
          long loc_phase_ext = gkyl_range_idx(&phase_ext_range, idx_edge);

          double* alpha_surf_ext_d = (double*) gkyl_array_fetch(alpha_surf, loc_phase_ext);
          double* sgn_alpha_surf_ext_d = (double*) gkyl_array_fetch(sgn_alpha_surf, loc_phase_ext);
          int* const_sgn_alpha_ext_d = (int*) gkyl_array_fetch(const_sgn_alpha, loc_phase_ext);
          const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->phase_grid.dx, 
            vmap_d, vmapSq_d, up->charge, up->mass, bmag_d, jacobtot_inv_d,
            cmag_d, b_i_d, phi_d, alpha_surf_ext_d, sgn_alpha_surf_ext_d);
        }  
      }
    }
  }
}

// Host-side wrapper for gyrokinetic surface alpha calculation
void gkyl_dg_calc_gyrokinetic_vars_alpha_surf_cu(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
  dim3 dimGrid, dimBlock;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, up->cdim+1);
  gkyl_dg_calc_gyrokinetic_vars_alpha_surf_cu_kernel<<<dimGrid, dimBlock>>>(up->on_dev, 
    *conf_range, *phase_range, *phase_ext_range, phi->on_dev,
    alpha_surf->on_dev, sgn_alpha_surf->on_dev, const_sgn_alpha->on_dev);
}

// CUDA kernel to set device pointers to gyrokinetic vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_gyrokinetic_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  int cdim, int vdim, int poly_order, enum gkyl_gkmodel_id gkmodel_id)
{
  if (gkmodel_id == GKYL_GK_MODEL_NO_BY) {
    for (int d=0; d<cdim; ++d) {
      up->alpha_surf[d] = choose_gyrokinetic_alpha_no_by_surf_conf_kern(d, cdim, vdim, poly_order);
      up->alpha_edge_surf[d] = choose_gyrokinetic_alpha_no_by_edge_surf_conf_kern(d, cdim, vdim, poly_order);
    }
    up->alpha_surf[cdim] = choose_gyrokinetic_alpha_no_by_surf_vpar_kern(cdim, vdim, poly_order);
  }
  else {
    for (int d=0; d<cdim; ++d) {
      up->alpha_surf[d] = choose_gyrokinetic_alpha_surf_conf_kern(d, cdim, vdim, poly_order);
      up->alpha_edge_surf[d] = choose_gyrokinetic_alpha_edge_surf_conf_kern(d, cdim, vdim, poly_order);
    }
    up->alpha_surf[cdim] = choose_gyrokinetic_alpha_surf_vpar_kern(cdim, vdim, poly_order);
  }
}

gkyl_dg_calc_gyrokinetic_vars*
gkyl_dg_calc_gyrokinetic_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass, enum gkyl_gkmodel_id gkmodel_id, 
  const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map)
{
  struct gkyl_dg_calc_gyrokinetic_vars *up = (struct gkyl_dg_calc_gyrokinetic_vars*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  up->gk_geom = geom_ho->on_dev;
  up->vel_map = vel_map_ho->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_gyrokinetic_vars *up_cu = (struct gkyl_dg_calc_gyrokinetic_vars*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_gyrokinetic_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_gyrokinetic_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order, gkmodel_id);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // Updater should store host pointers.
  up->gk_geom = geom_ho; 
  up->vel_map = vel_map_ho; 
  
  return up;
}
