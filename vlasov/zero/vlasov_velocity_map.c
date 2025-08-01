#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_velocity_map_priv.h>

void
gkyl_vlasov_velocity_map_new(const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange, 
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM], 
  struct gkyl_array *vmap, struct gkyl_array *jacob_vel, 
  struct gkyl_array *vmap_pgkyl, struct gkyl_array *jacob_vel_pgkyl, 
  struct gkyl_array *vmap_avg_pgkyl, struct gkyl_array *jacob_vel_avg_pgkyl, 
  struct gkyl_array *jacob_vel_gauss)
{
  int vdim = vgrid->ndim;
  struct gkyl_array *v_nodal[3];
  struct gkyl_array *v_cubic[3];
  struct gkyl_dg_basis_op_mem *mem[3];

  // Make 1D cubic basis for constructing C^1 expansion
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  // 1D ranges for indexing 1D cubic bases
  struct gkyl_range local[3], local_ext[3];

  // Loop over number of dimensions and construct 1D mappings
  for (int i=0; i<vdim; ++i) {
    double lower[] = { vgrid->lower[i] }, upper[] = { vgrid->upper[i] };
    int cells[] = { vgrid->cells[i] };

    struct gkyl_rect_grid grid_1d;
    gkyl_rect_grid_init(&grid_1d, 1, lower, upper, cells); 

    // nodal grid used for constructing physical coordinates
    double nc_lower[] = { lower[0] - 0.5*grid_1d.dx[0] };
    double nc_upper[] = { upper[0] + 0.5*grid_1d.dx[0] };
    int nc_cells[] = { cells[0] + 1 };
    struct gkyl_rect_grid nc_grid;
    gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

    int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
    gkyl_create_grid_ranges(&grid_1d, nghost, &local_ext[i], &local[i]);

    struct gkyl_range nc_local, nc_local_ext;
    gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

    v_nodal[i] = gkyl_array_new(GKYL_DOUBLE, 1, cells[0]+1);
    v_cubic[i] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext[i].volume);
    mem[i] = gkyl_dg_alloc_cubic_1d(cells[0]);
    double xn[1];

    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid_1d, iter.idx, xn);

      double *pn = gkyl_array_fetch(v_nodal[i], nidx);
      inp_vmap[i].eval_vmap(0.0, xn, pn, inp_vmap[i].ctx);
    }
    
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem[i], cells[0], grid_1d.dx[0],
      v_nodal[i], v_cubic[i]);
  }

  // initialize the mapping
  vmap_cubic_t vmap_op = choose_vmap_kern(vdim); 
  const double *v_cubic_dir[3]; // 1D cubic in each direction
  int vidx_1D[1]; // 1D index for indexing correct cubic mapping

  struct gkyl_range_iter iter_vmap;
  gkyl_range_iter_init(&iter_vmap, vrange);
  while (gkyl_range_iter_next(&iter_vmap)) {
    long loc_vel = gkyl_range_idx(vrange, iter_vmap.idx);

    for (int i=0; i<vdim; ++i) {
      vidx_1D[0] = iter_vmap.idx[i];
      long loc_vel_1D = gkyl_range_idx(&local[i], vidx_1D);
      v_cubic_dir[i] = gkyl_array_cfetch(v_cubic[i], loc_vel_1D);
    }
    double *vmap_d = gkyl_array_fetch(vmap, loc_vel);
    double *jacob_vel_d = gkyl_array_fetch(jacob_vel, loc_vel);
    double *vmap_pgkyl_d = gkyl_array_fetch(vmap_pgkyl, loc_vel);
    double *jacob_vel_pgkyl_d = gkyl_array_fetch(jacob_vel_pgkyl, loc_vel);
    double *vmap_avg_pgkyl_d = gkyl_array_fetch(vmap_avg_pgkyl, loc_vel);
    double *jacob_vel_avg_pgkyl_d = gkyl_array_fetch(jacob_vel_avg_pgkyl, loc_vel);
    double *jacob_vel_gauss_d = gkyl_array_fetch(jacob_vel_gauss, loc_vel);
    
    vmap_op(vgrid->dx, v_cubic_dir, vmap_d, jacob_vel_d, 
      vmap_pgkyl_d, jacob_vel_pgkyl_d, 
      vmap_avg_pgkyl_d, jacob_vel_avg_pgkyl_d, 
      jacob_vel_gauss_d);
  }

  // free temporary memory
  for (int i=0; i<vdim; ++i) {
    gkyl_array_release(v_nodal[i]);
    gkyl_array_release(v_cubic[i]);
    gkyl_dg_basis_op_mem_release(mem[i]);
  }
}