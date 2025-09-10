#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_deflate_geo_priv.h>



struct gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis,const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, bool use_gpu){

  gkyl_deflate_geo *up = gkyl_malloc(sizeof(gkyl_deflate_geo));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;
  up->deflated_basis = deflated_cbasis;
  up->rem_dirs = gkyl_malloc(3*sizeof(int));
  for(int i=0; i<3; i++)
    up->rem_dirs[i] = rem_dirs[i]; // 1 to remove
  up->kernel = deflate_geo_choose_kernel(up->rem_dirs, up->deflated_grid->ndim, cbasis->b_type, cbasis->poly_order);

  return up;
}


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp){
  // Inflated grid will always be 1 long in other directions
  // So deflate reange to the cell at index 1 in the ignored directions (z or x and y). 1 and not 0 because it is a local range 

  struct gkyl_range update_range;
  int loc_dir[3] = {0.};

  for(int i = 0; i<3; i++){
    loc_dir[i] = 1;
    if(up->rem_dirs[i] == 1){
      loc_dir[i] = 1;
    }
  }

  int do_idx[3];
  for(int i = 0; i<3; i++)
    do_idx[i] = loc_dir[i];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_range);

  int idx_deflated[3];
  while(gkyl_range_iter_next(&iter)){
      int count = 0;
      for(int i = 0; i<3; i++){
        if(up->rem_dirs[i]==0){
          do_idx[i] = iter.idx[count];
          count += 1;
        }
      }

      long loc = gkyl_range_idx(range, do_idx);
      const double *fld = gkyl_array_cfetch(field, loc);

      long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
      double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
      for(int c = 0; c<ncomp; c++){
        up->kernel(&fld[c*up->basis->num_basis], &fld_deflated[c*up->deflated_basis->num_basis]);
      }
 }

}

void gkyl_deflate_geo_release(gkyl_deflate_geo* up){
  gkyl_free(up->rem_dirs);
  gkyl_free(up);
}


