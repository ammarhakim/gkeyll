#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_alloc_flags_priv.h>
#include <assert.h>
#include <float.h>

struct gk_geometry*
gkyl_gk_geometry_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp, bool use_gpu)
{

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_gk_geometry_cu_dev_new(geo_host, geometry_inp);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;
  up->geqdsk_sign_convention = geo_host->geqdsk_sign_convention;
  up->has_LCFS = geometry_inp->has_LCFS;
  if (up->has_LCFS) {
    up->x_LCFS = geometry_inp->x_LCFS;
    // Check that the split happens within the domain.
    assert((up->grid.lower[0] <= up->x_LCFS) && (up->x_LCFS <= up->grid.upper[0]));
    // Check that the split happens at a cell boundary;
    double needint = (up->x_LCFS - up->grid.lower[0])/up->grid.dx[0];
    double rem_floor = fabs(needint-floor(needint));
    double rem_ceil = fabs(needint-ceil(needint));
    if (rem_floor < 1.0e-12) {
      up->idx_LCFS_lo = (int) floor(needint);
    }
    else if (rem_ceil < 1.0e-12) {
      up->idx_LCFS_lo = (int) ceil(needint);
    }
    else {
      fprintf(stderr, "x_LCFS = %.9e must be at a cell boundary.\n", up->x_LCFS);
      assert(false);
    }
  }

  // bmag, metrics and derived geo quantities
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxxj = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxyj = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gyyj = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxzj = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->eps2 = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

struct gkyl_rect_grid gkyl_gk_geometry_augment_grid(struct gkyl_rect_grid grid, struct gkyl_gk_geometry_inp geometry)
{
  struct gkyl_rect_grid augmented_grid;
  int cells[3];
  double lower[3];
  double upper[3];

  if (grid.ndim==1) {
    cells[0] = 1;
    cells[1] = 1;
    cells[2] = grid.cells[0];

    lower[0] = geometry.world[0] - (geometry.world[0]>1e-14? fmin(1e-5, geometry.world[0]*0.1) : 1e-5);
    lower[1] = geometry.world[1] - 1e-1;
    lower[2] = grid.lower[0];

    upper[0] = geometry.world[0] + (geometry.world[0]>1e-14? fmin(1e-5, geometry.world[0]*0.1) : 1e-5);
    upper[1] = geometry.world[1] + 1e-1;
    upper[2] = grid.upper[0];
  }
  else if (grid.ndim==2) {
    cells[0] = grid.cells[0];
    cells[1] = 1;
    cells[2] = grid.cells[1];

    lower[0] = grid.lower[0];
    lower[1] = geometry.world[0] - 1e-5;
    lower[2] = grid.lower[1];

    upper[0] = grid.upper[0];
    upper[1] = geometry.world[0] + 1e-5;
    upper[2] = grid.upper[1];
  }

  gkyl_rect_grid_init(&augmented_grid, 3, lower, upper, cells);
  return augmented_grid;
}


void gkyl_gk_geometry_augment_local(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  if (inrange->ndim == 2) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = inrange->lower[0]-nghost[0];
    upper_ext[0] = inrange->upper[0]+nghost[0];
    lower[0] = inrange->lower[0];
    upper[0] = inrange->upper[0];

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[1]-nghost[1];
    upper_ext[2] = inrange->upper[1]+nghost[1];
    lower[2] = inrange->lower[1];
    upper[2] = inrange->upper[1];


    gkyl_range_init(ext_range, inrange->ndim+1, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
  else if (inrange->ndim == 1) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = 1 - 1;
    upper_ext[0] = 1 + 1;
    lower[0] = 1;
    upper[0] = 1;

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[0]-nghost[0];
    upper_ext[2] = inrange->upper[0]+nghost[0];
    lower[2] = inrange->lower[0];
    upper[2] = inrange->upper[0];



    gkyl_range_init(ext_range, inrange->ndim+2, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
}

double
gkyl_gk_geometry_reduce_bmag(struct gk_geometry* up, enum gkyl_array_op op)
{
  int cdim = up->grid.ndim;
  double b_m;
  if (op == GKYL_MIN)
    b_m = DBL_MAX;
  else if (op == GKYL_MAX)
    b_m = -DBL_MAX;
  else
    assert(false);

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cdim, up->basis.num_basis);
  up->basis.node_list(gkyl_array_fetch(nodes, 0));

  struct gkyl_array *bmag_ho = gkyl_array_new(GKYL_DOUBLE, up->bmag->ncomp, up->bmag->size);
  gkyl_array_copy(bmag_ho, up->bmag);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->local, iter.idx);
    double *b_d = gkyl_array_fetch(bmag_ho, linidx);
    double blog[cdim];
    for (int n = 0; n < up->basis.num_basis; n++) {
      const double *blog = gkyl_array_cfetch(nodes,n);
      double b = up->basis.eval_expand(blog, b_d);
      if (op == GKYL_MIN)
        b_m = GKYL_MIN2(b_m, b);
      else if (op == GKYL_MAX)
        b_m = GKYL_MAX2(b_m, b);
    }
  }

  gkyl_array_release(nodes);
  gkyl_array_release(bmag_ho);

  return b_m;
}

void
gkyl_gk_geometry_init_nodal_range( struct gkyl_range *nrange, struct gkyl_range *range, int poly_order)
{
    int nodes[GKYL_MAX_DIM];
    if (poly_order == 1) {
      for (int d=0; d<range->ndim; ++d)
        nodes[d] = gkyl_range_shape(range, d) + 1;
    }
    if (poly_order == 2) {
      for (int d=0; d<range->ndim; ++d)
        nodes[d] = 2*gkyl_range_shape(range, d) + 1;
    }
    gkyl_range_init_from_shape(nrange, range->ndim, nodes);

}

void
gkyl_gk_geometry_init_nodal_grid(struct gkyl_rect_grid *ngrid, struct gkyl_rect_grid *grid, struct gkyl_range *nrange)
{
    double lower[GKYL_MAX_DIM];
    double upper[GKYL_MAX_DIM];
    int cells[GKYL_MAX_DIM];
    for (int i=0; i<nrange->ndim; ++i) {
      lower[i] = grid->lower[i];
      upper[i] = grid->upper[i];
      cells[i] = gkyl_range_shape(nrange, i);
    }
    gkyl_rect_grid_init(ngrid, nrange->ndim, lower, upper, cells);
}

struct gk_geometry*
gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->grid = geometry_inp->grid;
  up->geqdsk_sign_convention = up_3d->geqdsk_sign_convention;
  up->has_LCFS = up_3d->has_LCFS;
  up->x_LCFS = up_3d->x_LCFS;
  up->idx_LCFS_lo = up_3d->idx_LCFS_lo;

  // bmag, metrics and derived geo quantities
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);

  // Now fill the arrays by deflation
  int rem_dirs[3] = {0};
  if (up->grid.ndim==1) {
    rem_dirs[0] = 1;
    rem_dirs[1] = 1;
  }
  else if (up->grid.ndim==2) {
    rem_dirs[1] = 1;
  }
  struct gkyl_deflate_geo* deflator = gkyl_deflate_geo_new(&up_3d->basis, &up->basis, &up_3d->grid, &up->grid, rem_dirs, false);

  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->mc2p, up->mc2p, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->mc2nu_pos, up->mc2nu_pos, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag, up->bmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->g_ij, up->g_ij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->g_ij_neut, up->g_ij_neut, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->dxdz, up->dxdz, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->dzdx, up->dzdx, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->dualmag, up->dualmag, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->normals, up->normals, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobgeo, up->jacobgeo, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobgeo_inv, up->jacobgeo_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gij, up->gij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gij_neut, up->gij_neut, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->b_i, up->b_i, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bcart, up->bcart, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->cmag, up->cmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobtot, up->jacobtot, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobtot_inv, up->jacobtot_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag_inv, up->bmag_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag_inv_sq, up->bmag_inv_sq, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxxj, up->gxxj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxyj, up->gxyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gyyj, up->gyyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxzj, up->gxzj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->eps2, up->eps2, 1);
  // Done deflating
  gkyl_deflate_geo_release(deflator);

  if (up->grid.ndim==1) {
    // In 1D geometry, make mapc2p a function of only Z and mc2nu_pos only a function of length along field line
    gkyl_array_set_offset(up->mc2p_deflated, 1.0, up->mc2p, 1 * up->basis.num_basis);
    gkyl_array_set_offset(up->mc2nu_pos_deflated, 1.0, up->mc2nu_pos, 2 * up->basis.num_basis);
  }
  else if (up->grid.ndim==2) {
    // In 2D geometry, make mapc2p a function of only R and Z 
    // and mc2nu_pos only a function of psi and length along field line
    struct gkyl_array *temp = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    gkyl_array_set_offset(up->mc2p_deflated, 1.0, up->mc2p, 0 * up->basis.num_basis); // 
    gkyl_array_set_offset(temp,              1.0, up->mc2p, 1 * up->basis.num_basis);
    gkyl_array_set_offset(up->mc2p_deflated, 1.0, temp,     1 * up->basis.num_basis);
    gkyl_array_set_offset(up->mc2nu_pos_deflated, 1.0, up->mc2nu_pos, 0 * up->basis.num_basis);
    gkyl_array_set_offset(temp,                   1.0, up->mc2nu_pos, 2 * up->basis.num_basis);
    gkyl_array_set_offset(up->mc2nu_pos_deflated, 1.0, temp,          1 * up->basis.num_basis);
    gkyl_array_release(temp);
  }
  else if (up->grid.ndim==3) {
    // In 3D geoemtry, these are identical
    gkyl_array_copy(up->mc2p_deflated, up->mc2p);
    gkyl_array_copy(up->mc2nu_pos_deflated, up->mc2nu_pos);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gk_geometry *up = container_of(ref, struct gk_geometry, ref_count);
  gkyl_array_release(up->mc2p);
  gkyl_array_release(up->mc2p_deflated);
  gkyl_array_release(up->mc2nu_pos);
  gkyl_array_release(up->mc2nu_pos_deflated);
  gkyl_array_release(up->bmag);
  gkyl_array_release(up->g_ij);
  gkyl_array_release(up->g_ij_neut);
  gkyl_array_release(up->jacobgeo);
  gkyl_array_release(up->jacobgeo_inv);
  gkyl_array_release(up->dxdz);
  gkyl_array_release(up->dzdx);
  gkyl_array_release(up->dualmag);
  gkyl_array_release(up->normals);
  gkyl_array_release(up->gij);
  gkyl_array_release(up->gij_neut);
  gkyl_array_release(up->b_i);
  gkyl_array_release(up->bcart);
  gkyl_array_release(up->cmag);
  gkyl_array_release(up->jacobtot);
  gkyl_array_release(up->jacobtot_inv);
  gkyl_array_release(up->bmag_inv);
  gkyl_array_release(up->bmag_inv_sq);
  gkyl_array_release(up->gxxj);
  gkyl_array_release(up->gxyj);
  gkyl_array_release(up->gyyj);
  gkyl_array_release(up->gxzj);
  gkyl_array_release(up->eps2);
  if (gkyl_gk_geometry_is_cu_dev(up)) 
    gkyl_cu_free(up->on_dev); 

  gkyl_free(up);
}

struct gk_geometry*
gkyl_gk_geometry_acquire(const struct gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



