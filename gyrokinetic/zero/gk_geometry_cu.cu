/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>
#include <assert.h>
}

__global__ static void
gk_geometry_set_corn_cu_kernel(struct gk_geometry *gk_geom,
   struct gkyl_array *mc2p, struct gkyl_array *mc2nu_pos, struct gkyl_array *bmag)
{
  gk_geom->geo_corn.mc2p = mc2p;
  gk_geom->geo_corn.mc2nu_pos = mc2nu_pos;
  gk_geom->geo_corn.bmag = bmag;
}

__global__ static void
gk_geometry_set_int_cu_kernel(struct gk_geometry *gk_geom,
   struct gkyl_array *mc2p, struct gkyl_array *bmag, struct gkyl_array *g_ij, struct gkyl_array *g_ij_neut,
   struct gkyl_array *dxdz, struct gkyl_array *dzdx, struct gkyl_array *dualmag, struct gkyl_array *normals,
   struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_ghost, struct gkyl_array *jacobgeo_inv, struct gkyl_array *gij,
   struct gkyl_array *gij_neut, struct gkyl_array *b_i, struct gkyl_array *bcart, struct gkyl_array *cmag,
   struct gkyl_array *jacobtot, struct gkyl_array *jacobtot_inv, struct gkyl_array *bmag_inv, struct gkyl_array *bmag_inv_sq,
   struct gkyl_array *gxxj, struct gkyl_array *gxyj, struct gkyl_array *gyyj, struct gkyl_array *gxzj,
   struct gkyl_array *eps2, struct gkyl_array *rtg33inv, struct gkyl_array *dualcurlbhatoverB, struct gkyl_array *bioverJB, struct gkyl_array *B3, struct gkyl_array *dualcurlbhat)
{
  gk_geom->geo_int.mc2p = mc2p;
  gk_geom->geo_int.bmag = bmag;
  gk_geom->geo_int.g_ij = g_ij;
  gk_geom->geo_int.g_ij_neut = g_ij_neut;
  gk_geom->geo_int.dxdz = dxdz;
  gk_geom->geo_int.dzdx = dzdx;
  gk_geom->geo_int.dualmag = dualmag;
  gk_geom->geo_int.normals = normals;
  gk_geom->geo_int.jacobgeo = jacobgeo;
  gk_geom->geo_int.jacobgeo_ghost = jacobgeo_ghost;
  gk_geom->geo_int.jacobgeo_inv = jacobgeo_inv;
  gk_geom->geo_int.gij = gij;
  gk_geom->geo_int.gij_neut = gij_neut;
  gk_geom->geo_int.b_i = b_i;
  gk_geom->geo_int.bcart = bcart;
  gk_geom->geo_int.cmag = cmag;
  gk_geom->geo_int.jacobtot = jacobtot;
  gk_geom->geo_int.jacobtot_inv = jacobtot_inv;
  gk_geom->geo_int.bmag_inv = bmag_inv;
  gk_geom->geo_int.bmag_inv_sq = bmag_inv_sq;
  gk_geom->geo_int.gxxj = gxxj;
  gk_geom->geo_int.gxyj = gxyj;
  gk_geom->geo_int.gyyj = gyyj;
  gk_geom->geo_int.gxzj = gxzj;
  gk_geom->geo_int.eps2 = eps2;
  gk_geom->geo_int.rtg33inv = rtg33inv;
  gk_geom->geo_int.dualcurlbhatoverB = dualcurlbhatoverB;
  gk_geom->geo_int.dualcurlbhat = dualcurlbhat;
  gk_geom->geo_int.bioverJB = bioverJB;
  gk_geom->geo_int.B3 = B3;
}

__global__ static void
gk_geometry_set_surf_cu_kernel(struct gk_geometry *gk_geom, int dir,
   struct gkyl_array *bmag,  struct gkyl_array *jacobgeo,  struct gkyl_array *jacobgeo_sync, 
   struct gkyl_array *b_i,  struct gkyl_array *cmag,  struct gkyl_array *jacobtot_inv,
   struct gkyl_array *B3, struct gkyl_array *normcurlbhat, struct gkyl_array *normals, struct gkyl_array *lenr)
{
  gk_geom->geo_surf[dir].bmag = bmag;
  gk_geom->geo_surf[dir].jacobgeo = jacobgeo;
  gk_geom->geo_surf[dir].jacobgeo_sync = jacobgeo_sync;
  gk_geom->geo_surf[dir].b_i = b_i;
  gk_geom->geo_surf[dir].cmag = cmag;
  gk_geom->geo_surf[dir].jacobtot_inv = jacobtot_inv;
  gk_geom->geo_surf[dir].B3 = B3;
  gk_geom->geo_surf[dir].normcurlbhat = normcurlbhat;
  gk_geom->geo_surf[dir].normals = normals;
  gk_geom->geo_surf[dir].lenr = lenr;
}

// Host-side wrapper for set_corn_cu_kernel
void
gkyl_geometry_set_corn_cu(struct gk_geometry *gk_geom, struct gk_geom_corn *geo_corn)
{
  gk_geometry_set_corn_cu_kernel<<<1,1>>>(gk_geom, 
        geo_corn->mc2p->on_dev, geo_corn->mc2nu_pos->on_dev, geo_corn->bmag->on_dev);
}

// Host-side wrapper for set_int_cu_kernel
void
gkyl_geometry_set_int_cu(struct gk_geometry *gk_geom, struct gk_geom_int *geo_int)
{
  gk_geometry_set_int_cu_kernel<<<1,1>>>(gk_geom,
   geo_int->mc2p->on_dev, geo_int->bmag->on_dev, geo_int->g_ij->on_dev, geo_int->g_ij_neut->on_dev,
   geo_int->dxdz->on_dev, geo_int->dzdx->on_dev, geo_int->dualmag->on_dev, geo_int->normals->on_dev,
   geo_int->jacobgeo->on_dev, geo_int->jacobgeo_ghost->on_dev, geo_int->jacobgeo_inv->on_dev, geo_int->gij->on_dev,
   geo_int->gij_neut->on_dev, geo_int->b_i->on_dev, geo_int->bcart->on_dev, geo_int->cmag->on_dev,
   geo_int->jacobtot->on_dev, geo_int->jacobtot_inv->on_dev, geo_int->bmag_inv->on_dev, geo_int->bmag_inv_sq->on_dev,
   geo_int->gxxj->on_dev, geo_int->gxyj->on_dev, geo_int->gyyj->on_dev, geo_int->gxzj->on_dev,
   geo_int->eps2->on_dev, geo_int->rtg33inv->on_dev, geo_int->dualcurlbhatoverB->on_dev, geo_int->bioverJB->on_dev, geo_int->B3->on_dev, geo_int->dualcurlbhat->on_dev);
}

// Host-side wrapper for set_surf_cu_kernel
void
gkyl_geometry_set_surf_cu(struct gk_geometry *gk_geom, struct gk_geom_surf *geo_surf, int dir)
{
  gk_geometry_set_surf_cu_kernel<<<1,1>>>(gk_geom, dir,
    geo_surf->bmag->on_dev, geo_surf->jacobgeo->on_dev, geo_surf->jacobgeo_sync->on_dev, 
    geo_surf->b_i->on_dev, geo_surf->cmag->on_dev, geo_surf->jacobtot_inv->on_dev, geo_surf->B3->on_dev, geo_surf->normcurlbhat->on_dev, geo_surf->normals->on_dev, geo_surf->lenr->on_dev);
}

struct gk_geom_corn*
gk_geometry_corn_cu_dev_alloc(struct gk_geom_corn up_corn_host)
{
  struct gk_geom_corn *up_corn_dev = (struct gk_geom_corn*) gkyl_malloc(sizeof(struct gk_geom_corn));
  up_corn_dev->mc2p = gkyl_array_cu_dev_new(up_corn_host.mc2p->type, up_corn_host.mc2p->ncomp, up_corn_host.mc2p->size);
  up_corn_dev->mc2nu_pos = gkyl_array_cu_dev_new(up_corn_host.mc2nu_pos->type, up_corn_host.mc2nu_pos->ncomp, up_corn_host.mc2nu_pos->size);
  up_corn_dev->bmag = gkyl_array_cu_dev_new(up_corn_host.bmag->type, up_corn_host.bmag->ncomp, up_corn_host.bmag->size);
  up_corn_dev->mc2p_deflated = gkyl_array_cu_dev_new(up_corn_host.mc2p_deflated->type, up_corn_host.mc2p_deflated->ncomp, up_corn_host.mc2p_deflated->size);
  up_corn_dev->mc2nu_pos_deflated = gkyl_array_cu_dev_new(up_corn_host.mc2nu_pos_deflated->type, up_corn_host.mc2nu_pos_deflated->ncomp, up_corn_host.mc2nu_pos_deflated->size);
  return up_corn_dev;
}

struct gk_geom_int*
gk_geometry_int_cu_dev_alloc(struct gk_geom_int up_int_host)
{
  struct gk_geom_int *up_int_dev = (struct gk_geom_int*) gkyl_malloc(sizeof(struct gk_geom_int));
  up_int_dev->mc2p = gkyl_array_cu_dev_new(up_int_host.mc2p->type, up_int_host.mc2p->ncomp, up_int_host.mc2p->size);
  up_int_dev->bmag = gkyl_array_cu_dev_new(up_int_host.bmag->type, up_int_host.bmag->ncomp, up_int_host.bmag->size);
  up_int_dev->g_ij = gkyl_array_cu_dev_new(up_int_host.g_ij->type, up_int_host.g_ij->ncomp, up_int_host.g_ij->size);
  up_int_dev->g_ij_neut = gkyl_array_cu_dev_new(up_int_host.g_ij_neut->type, up_int_host.g_ij_neut->ncomp, up_int_host.g_ij_neut->size);
  up_int_dev->dxdz = gkyl_array_cu_dev_new(up_int_host.dxdz->type, up_int_host.dxdz->ncomp, up_int_host.dxdz->size);
  up_int_dev->dzdx = gkyl_array_cu_dev_new(up_int_host.dzdx->type, up_int_host.dzdx->ncomp, up_int_host.dzdx->size);
  up_int_dev->dualmag = gkyl_array_cu_dev_new(up_int_host.dualmag->type, up_int_host.dualmag->ncomp, up_int_host.dualmag->size);
  up_int_dev->normals = gkyl_array_cu_dev_new(up_int_host.normals->type, up_int_host.normals->ncomp, up_int_host.normals->size);
  up_int_dev->jacobgeo = gkyl_array_cu_dev_new(up_int_host.jacobgeo->type, up_int_host.jacobgeo->ncomp, up_int_host.jacobgeo->size);
  up_int_dev->jacobgeo_ghost = gkyl_array_cu_dev_new(up_int_host.jacobgeo_ghost->type, up_int_host.jacobgeo_ghost->ncomp, up_int_host.jacobgeo_ghost->size);
  up_int_dev->jacobgeo_inv = gkyl_array_cu_dev_new(up_int_host.jacobgeo_inv->type, up_int_host.jacobgeo_inv->ncomp, up_int_host.jacobgeo_inv->size);
  up_int_dev->gij = gkyl_array_cu_dev_new(up_int_host.gij->type, up_int_host.gij->ncomp, up_int_host.gij->size);
  up_int_dev->gij_neut = gkyl_array_cu_dev_new(up_int_host.gij_neut->type, up_int_host.gij_neut->ncomp, up_int_host.gij_neut->size);
  up_int_dev->b_i = gkyl_array_cu_dev_new(up_int_host.b_i->type, up_int_host.b_i->ncomp, up_int_host.b_i->size);
  up_int_dev->bcart = gkyl_array_cu_dev_new(up_int_host.bcart->type, up_int_host.bcart->ncomp, up_int_host.bcart->size);
  up_int_dev->cmag = gkyl_array_cu_dev_new(up_int_host.cmag->type, up_int_host.cmag->ncomp, up_int_host.cmag->size);
  up_int_dev->jacobtot = gkyl_array_cu_dev_new(up_int_host.jacobtot->type, up_int_host.jacobtot->ncomp, up_int_host.jacobtot->size);
  up_int_dev->jacobtot_inv = gkyl_array_cu_dev_new(up_int_host.jacobtot_inv->type, up_int_host.jacobtot_inv->ncomp, up_int_host.jacobtot_inv->size);
  up_int_dev->bmag_inv = gkyl_array_cu_dev_new(up_int_host.bmag_inv->type, up_int_host.bmag_inv->ncomp, up_int_host.bmag_inv->size);
  up_int_dev->bmag_inv_sq = gkyl_array_cu_dev_new(up_int_host.bmag_inv_sq->type, up_int_host.bmag_inv_sq->ncomp, up_int_host.bmag_inv_sq->size);
  up_int_dev->gxxj = gkyl_array_cu_dev_new(up_int_host.gxxj->type, up_int_host.gxxj->ncomp, up_int_host.gxxj->size);
  up_int_dev->gxyj = gkyl_array_cu_dev_new(up_int_host.gxyj->type, up_int_host.gxyj->ncomp, up_int_host.gxyj->size);
  up_int_dev->gyyj = gkyl_array_cu_dev_new(up_int_host.gyyj->type, up_int_host.gyyj->ncomp, up_int_host.gyyj->size);
  up_int_dev->gxzj = gkyl_array_cu_dev_new(up_int_host.gxzj->type, up_int_host.gxzj->ncomp, up_int_host.gxzj->size);
  up_int_dev->eps2 = gkyl_array_cu_dev_new(up_int_host.eps2->type, up_int_host.eps2->ncomp, up_int_host.eps2->size);
  up_int_dev->rtg33inv = gkyl_array_cu_dev_new(up_int_host.rtg33inv->type, up_int_host.rtg33inv->ncomp, up_int_host.rtg33inv->size);
  up_int_dev->dualcurlbhatoverB = gkyl_array_cu_dev_new(up_int_host.dualcurlbhatoverB->type, up_int_host.dualcurlbhatoverB->ncomp, up_int_host.dualcurlbhatoverB->size);
  up_int_dev->bioverJB = gkyl_array_cu_dev_new(up_int_host.bioverJB->type, up_int_host.bioverJB->ncomp, up_int_host.bioverJB->size);
  up_int_dev->B3 = gkyl_array_cu_dev_new(up_int_host.B3->type, up_int_host.B3->ncomp, up_int_host.B3->size);
  up_int_dev->dualcurlbhat = gkyl_array_cu_dev_new(up_int_host.dualcurlbhat->type, up_int_host.dualcurlbhat->ncomp, up_int_host.dualcurlbhat->size);
  return up_int_dev;
}

struct gk_geom_surf*
gk_geometry_surf_cu_dev_alloc(struct gk_geom_surf up_surf_host)
{
  struct gk_geom_surf *up_surf_dev = (struct gk_geom_surf*) gkyl_malloc(sizeof(struct gk_geom_surf));
  up_surf_dev->bmag = gkyl_array_cu_dev_new(up_surf_host.bmag->type, up_surf_host.bmag->ncomp, up_surf_host.bmag->size);
  up_surf_dev->jacobgeo = gkyl_array_cu_dev_new(up_surf_host.jacobgeo->type, up_surf_host.jacobgeo->ncomp, up_surf_host.jacobgeo->size);
  up_surf_dev->jacobgeo_sync = gkyl_array_cu_dev_new(up_surf_host.jacobgeo_sync->type, up_surf_host.jacobgeo_sync->ncomp, up_surf_host.jacobgeo_sync->size);
  up_surf_dev->b_i = gkyl_array_cu_dev_new(up_surf_host.b_i->type, up_surf_host.b_i->ncomp, up_surf_host.b_i->size);
  up_surf_dev->cmag = gkyl_array_cu_dev_new(up_surf_host.cmag->type, up_surf_host.cmag->ncomp, up_surf_host.cmag->size);
  up_surf_dev->jacobtot_inv = gkyl_array_cu_dev_new(up_surf_host.jacobtot_inv->type, up_surf_host.jacobtot_inv->ncomp, up_surf_host.jacobtot_inv->size);
  up_surf_dev->B3 = gkyl_array_cu_dev_new(up_surf_host.B3->type, up_surf_host.B3->ncomp, up_surf_host.B3->size);
  up_surf_dev->normcurlbhat = gkyl_array_cu_dev_new(up_surf_host.normcurlbhat->type, up_surf_host.normcurlbhat->ncomp, up_surf_host.normcurlbhat->size);
  up_surf_dev->normals = gkyl_array_cu_dev_new(up_surf_host.normals->type, up_surf_host.normals->ncomp, up_surf_host.normals->size);
  up_surf_dev->lenr = gkyl_array_cu_dev_new(up_surf_host.lenr->type, up_surf_host.lenr->ncomp, up_surf_host.lenr->size);
  return up_surf_dev;
}

// CPU interface to create and track a GPU object
struct gk_geometry* 
gkyl_gk_geometry_cu_dev_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up =(struct gk_geometry*) gkyl_malloc(sizeof(struct gk_geometry));

  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;
  if (up->grid.ndim > 1) {
    gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
    up->num_surf_basis = up->surf_basis.num_basis;
  }
  else {
    up->num_surf_basis = 1;
  }
  up->geqdsk_sign_convention = geo_host->geqdsk_sign_convention;
  up->has_LCFS = geo_host->has_LCFS;
  if (up->has_LCFS) {
    up->x_LCFS = geo_host->x_LCFS;
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

  struct gk_geom_corn *geo_corn_dev = gk_geometry_corn_cu_dev_alloc(geo_host->geo_corn);
  struct gk_geom_int *geo_int_dev = gk_geometry_int_cu_dev_alloc(geo_host->geo_int);
  struct gk_geom_surf *geo_surf_dev[up->grid.ndim];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    geo_surf_dev[dir] = gk_geometry_surf_cu_dev_alloc(geo_host->geo_surf[dir]);


  gkyl_array_copy(geo_corn_dev->mc2p, geo_host->geo_corn.mc2p);
  gkyl_array_copy(geo_corn_dev->mc2nu_pos, geo_host->geo_corn.mc2nu_pos);
  gkyl_array_copy(geo_corn_dev->bmag, geo_host->geo_corn.bmag);

  gkyl_array_copy(geo_int_dev->mc2p, geo_host->geo_int.mc2p);
  gkyl_array_copy(geo_int_dev->bmag, geo_host->geo_int.bmag);
  gkyl_array_copy(geo_int_dev->g_ij, geo_host->geo_int.g_ij);
  gkyl_array_copy(geo_int_dev->g_ij_neut, geo_host->geo_int.g_ij_neut);
  gkyl_array_copy(geo_int_dev->dxdz, geo_host->geo_int.dxdz);
  gkyl_array_copy(geo_int_dev->dzdx, geo_host->geo_int.dzdx);
  gkyl_array_copy(geo_int_dev->dualmag, geo_host->geo_int.dualmag);
  gkyl_array_copy(geo_int_dev->normals, geo_host->geo_int.normals);
  gkyl_array_copy(geo_int_dev->jacobgeo , geo_host->geo_int.jacobgeo);
  gkyl_array_copy(geo_int_dev->jacobgeo_ghost , geo_host->geo_int.jacobgeo_ghost);
  gkyl_array_copy(geo_int_dev->jacobgeo_inv, geo_host->geo_int.jacobgeo_inv);
  gkyl_array_copy(geo_int_dev->gij, geo_host->geo_int.gij);
  gkyl_array_copy(geo_int_dev->gij_neut, geo_host->geo_int.gij_neut);
  gkyl_array_copy(geo_int_dev->b_i, geo_host->geo_int.b_i);
  gkyl_array_copy(geo_int_dev->bcart, geo_host->geo_int.bcart);
  gkyl_array_copy(geo_int_dev->cmag, geo_host->geo_int.cmag);
  gkyl_array_copy(geo_int_dev->jacobtot, geo_host->geo_int.jacobtot);
  gkyl_array_copy(geo_int_dev->jacobtot_inv, geo_host->geo_int.jacobtot_inv);
  gkyl_array_copy(geo_int_dev->bmag_inv, geo_host->geo_int.bmag_inv);
  gkyl_array_copy(geo_int_dev->bmag_inv_sq, geo_host->geo_int.bmag_inv_sq);
  gkyl_array_copy(geo_int_dev->gxxj, geo_host->geo_int.gxxj);
  gkyl_array_copy(geo_int_dev->gxyj, geo_host->geo_int.gxyj);
  gkyl_array_copy(geo_int_dev->gyyj, geo_host->geo_int.gyyj);
  gkyl_array_copy(geo_int_dev->gxzj, geo_host->geo_int.gxzj);
  gkyl_array_copy(geo_int_dev->eps2, geo_host->geo_int.eps2);
  gkyl_array_copy(geo_int_dev->rtg33inv, geo_host->geo_int.rtg33inv);
  gkyl_array_copy(geo_int_dev->dualcurlbhatoverB, geo_host->geo_int.dualcurlbhatoverB);
  gkyl_array_copy(geo_int_dev->bioverJB, geo_host->geo_int.bioverJB);
  gkyl_array_copy(geo_int_dev->B3, geo_host->geo_int.B3);
  gkyl_array_copy(geo_int_dev->dualcurlbhat, geo_host->geo_int.dualcurlbhat);

  for (int dir=0; dir<up->grid.ndim; ++dir) {
   gkyl_array_copy(geo_surf_dev[dir]->bmag, geo_host->geo_surf[dir].bmag);
   gkyl_array_copy(geo_surf_dev[dir]->jacobgeo, geo_host->geo_surf[dir].jacobgeo);
   gkyl_array_copy(geo_surf_dev[dir]->jacobgeo_sync, geo_host->geo_surf[dir].jacobgeo_sync);
   gkyl_array_copy(geo_surf_dev[dir]->b_i, geo_host->geo_surf[dir].b_i);
   gkyl_array_copy(geo_surf_dev[dir]->cmag, geo_host->geo_surf[dir].cmag);
   gkyl_array_copy(geo_surf_dev[dir]->jacobtot_inv, geo_host->geo_surf[dir].jacobtot_inv);
   gkyl_array_copy(geo_surf_dev[dir]->B3, geo_host->geo_surf[dir].B3);
   gkyl_array_copy(geo_surf_dev[dir]->normcurlbhat, geo_host->geo_surf[dir].normcurlbhat);
   gkyl_array_copy(geo_surf_dev[dir]->normals, geo_host->geo_surf[dir].normals);
   gkyl_array_copy(geo_surf_dev[dir]->lenr, geo_host->geo_surf[dir].lenr);
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);

  // Initialize the device geometry object
  struct gk_geometry *up_cu = (struct gk_geometry*) gkyl_cu_malloc(sizeof(struct gk_geometry));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gk_geometry), GKYL_CU_MEMCPY_H2D);
  gkyl_geometry_set_corn_cu(up_cu, geo_corn_dev);
  gkyl_geometry_set_int_cu(up_cu, geo_int_dev);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_geometry_set_surf_cu(up_cu, geo_surf_dev[dir], dir);
  up->on_dev = up_cu;

  // geometry object should store host pointer
  up->geo_corn.mc2p  = geo_corn_dev->mc2p;
  up->geo_corn.mc2nu_pos  = geo_corn_dev->mc2nu_pos;
  up->geo_corn.bmag  = geo_corn_dev->bmag;
  up->geo_corn.mc2p_deflated  = geo_corn_dev->mc2p_deflated;
  up->geo_corn.mc2nu_pos_deflated  = geo_corn_dev->mc2nu_pos_deflated;

  up->geo_int.mc2p  = geo_int_dev->mc2p;
  up->geo_int.bmag  = geo_int_dev->bmag;
  up->geo_int.g_ij  = geo_int_dev->g_ij;
  up->geo_int.g_ij_neut  = geo_int_dev->g_ij_neut;
  up->geo_int.dxdz  = geo_int_dev->dxdz;
  up->geo_int.dzdx  = geo_int_dev->dzdx;
  up->geo_int.dualmag  = geo_int_dev->dualmag;
  up->geo_int.normals  = geo_int_dev->normals;
  up->geo_int.jacobgeo  = geo_int_dev->jacobgeo;
  up->geo_int.jacobgeo_ghost  = geo_int_dev->jacobgeo_ghost;
  up->geo_int.jacobgeo_inv = geo_int_dev->jacobgeo_inv;
  up->geo_int.gij  = geo_int_dev->gij;
  up->geo_int.gij_neut  = geo_int_dev->gij_neut;
  up->geo_int.b_i  = geo_int_dev->b_i;
  up->geo_int.bcart  = geo_int_dev->bcart;
  up->geo_int.cmag  =  geo_int_dev->cmag;
  up->geo_int.jacobtot  = geo_int_dev->jacobtot;
  up->geo_int.jacobtot_inv = geo_int_dev->jacobtot_inv;
  up->geo_int.bmag_inv  = geo_int_dev->bmag_inv;
  up->geo_int.bmag_inv_sq = geo_int_dev->bmag_inv_sq;
  up->geo_int.gxxj  = geo_int_dev->gxxj;
  up->geo_int.gxyj  = geo_int_dev->gxyj;
  up->geo_int.gyyj  = geo_int_dev->gyyj;
  up->geo_int.gxzj  = geo_int_dev->gxzj;
  up->geo_int.eps2  = geo_int_dev->eps2;
  up->geo_int.rtg33inv  = geo_int_dev->rtg33inv;
  up->geo_int.dualcurlbhatoverB  = geo_int_dev->dualcurlbhatoverB;
  up->geo_int.bioverJB  = geo_int_dev->bioverJB;
  up->geo_int.B3  = geo_int_dev->B3;
  up->geo_int.dualcurlbhat  = geo_int_dev->dualcurlbhat;

  for (int dir=0; dir<up->grid.ndim; ++dir) {
   up->geo_surf[dir].bmag = geo_surf_dev[dir]->bmag;
   up->geo_surf[dir].jacobgeo = geo_surf_dev[dir]->jacobgeo;
   up->geo_surf[dir].jacobgeo_sync = geo_surf_dev[dir]->jacobgeo_sync;
   up->geo_surf[dir].b_i = geo_surf_dev[dir]->b_i;
   up->geo_surf[dir].cmag = geo_surf_dev[dir]->cmag;
   up->geo_surf[dir].jacobtot_inv = geo_surf_dev[dir]->jacobtot_inv;
   up->geo_surf[dir].B3 = geo_surf_dev[dir]->B3;
   up->geo_surf[dir].normcurlbhat = geo_surf_dev[dir]->normcurlbhat;
   up->geo_surf[dir].normals = geo_surf_dev[dir]->normals;
   up->geo_surf[dir].lenr = geo_surf_dev[dir]->lenr;
   gkyl_free(geo_surf_dev[dir]);
  }
  
  return up;
}

