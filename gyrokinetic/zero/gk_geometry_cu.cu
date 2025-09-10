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

  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *mc2p_dev = gkyl_array_cu_dev_new(geo_host->mc2p->type, geo_host->mc2p->ncomp, geo_host->mc2p->size);
  struct gkyl_array *mc2p_delfated_dev = gkyl_array_cu_dev_new(geo_host->mc2p_deflated->type, geo_host->mc2p_deflated->ncomp, geo_host->mc2p_deflated->size);
  struct gkyl_array *mc2nu_pos_dev = gkyl_array_cu_dev_new(geo_host->mc2nu_pos->type, geo_host->mc2nu_pos->ncomp, geo_host->mc2nu_pos->size);
  struct gkyl_array *mc2nu_pos_deflated_dev = gkyl_array_cu_dev_new(geo_host->mc2nu_pos_deflated->type, geo_host->mc2nu_pos_deflated->ncomp, geo_host->mc2nu_pos_deflated->size);
  struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(geo_host->bmag->type, geo_host->bmag->ncomp, geo_host->bmag->size);
  struct gkyl_array *g_ij_dev = gkyl_array_cu_dev_new(geo_host->g_ij->type, geo_host->g_ij->ncomp, geo_host->g_ij->size);
  struct gkyl_array *g_ij_neut_dev = gkyl_array_cu_dev_new(geo_host->g_ij_neut->type, geo_host->g_ij_neut->ncomp, geo_host->g_ij_neut->size);
  struct gkyl_array *dxdz_dev = gkyl_array_cu_dev_new(geo_host->dxdz->type, geo_host->dxdz->ncomp, geo_host->dxdz->size);
  struct gkyl_array *dzdx_dev = gkyl_array_cu_dev_new(geo_host->dzdx->type, geo_host->dzdx->ncomp, geo_host->dzdx->size);
  struct gkyl_array *dualmag_dev = gkyl_array_cu_dev_new(geo_host->dualmag->type, geo_host->dualmag->ncomp, geo_host->dualmag->size);
  struct gkyl_array *normals_dev = gkyl_array_cu_dev_new(geo_host->normals->type, geo_host->normals->ncomp, geo_host->normals->size);
  struct gkyl_array *jacobgeo_dev = gkyl_array_cu_dev_new(geo_host->jacobgeo->type, geo_host->jacobgeo->ncomp, geo_host->jacobgeo->size);
  struct gkyl_array *jacobgeo_inv_dev = gkyl_array_cu_dev_new(geo_host->jacobgeo_inv->type, geo_host->jacobgeo_inv->ncomp, geo_host->jacobgeo_inv->size);
  struct gkyl_array *gij_dev = gkyl_array_cu_dev_new(geo_host->gij->type, geo_host->gij->ncomp, geo_host->gij->size);
  struct gkyl_array *gij_neut_dev = gkyl_array_cu_dev_new(geo_host->gij_neut->type, geo_host->gij_neut->ncomp, geo_host->gij_neut->size);
  struct gkyl_array *b_i_dev = gkyl_array_cu_dev_new(geo_host->b_i->type, geo_host->b_i->ncomp, geo_host->b_i->size);
  struct gkyl_array *bcart_dev = gkyl_array_cu_dev_new(geo_host->bcart->type, geo_host->bcart->ncomp, geo_host->bcart->size);
  struct gkyl_array *cmag_dev = gkyl_array_cu_dev_new(geo_host->cmag->type, geo_host->cmag->ncomp, geo_host->cmag->size);
  struct gkyl_array *jacobtot_dev = gkyl_array_cu_dev_new(geo_host->jacobtot->type, geo_host->jacobtot->ncomp, geo_host->jacobtot->size);
  struct gkyl_array *jacobtot_inv_dev = gkyl_array_cu_dev_new(geo_host->jacobtot_inv->type, geo_host->jacobtot_inv->ncomp, geo_host->jacobtot_inv->size);
  struct gkyl_array *bmag_inv_dev = gkyl_array_cu_dev_new(geo_host->bmag_inv->type, geo_host->bmag_inv->ncomp, geo_host->bmag_inv->size);
  struct gkyl_array *bmag_inv_sq_dev = gkyl_array_cu_dev_new(geo_host->bmag_inv_sq->type, geo_host->bmag_inv_sq->ncomp, geo_host->bmag_inv_sq->size);
  struct gkyl_array *gxxj_dev = gkyl_array_cu_dev_new(geo_host->gxxj->type, geo_host->gxxj->ncomp, geo_host->gxxj->size);
  struct gkyl_array *gxyj_dev = gkyl_array_cu_dev_new(geo_host->gxyj->type, geo_host->gxyj->ncomp, geo_host->gxyj->size);
  struct gkyl_array *gyyj_dev = gkyl_array_cu_dev_new(geo_host->gyyj->type, geo_host->gyyj->ncomp, geo_host->gyyj->size);
  struct gkyl_array *gxzj_dev = gkyl_array_cu_dev_new(geo_host->gxzj->type, geo_host->gxzj->ncomp, geo_host->gxzj->size);
  struct gkyl_array *eps2_dev = gkyl_array_cu_dev_new(geo_host->eps2->type, geo_host->eps2->ncomp, geo_host->eps2->size);

  gkyl_array_copy(mc2p_dev, geo_host->mc2p);
  gkyl_array_copy(mc2p_delfated_dev, geo_host->mc2p_deflated);
  gkyl_array_copy(mc2nu_pos_dev, geo_host->mc2nu_pos);
  gkyl_array_copy(mc2nu_pos_deflated_dev, geo_host->mc2nu_pos_deflated);
  gkyl_array_copy(bmag_dev, geo_host->bmag);
  gkyl_array_copy(g_ij_dev, geo_host->g_ij);
  gkyl_array_copy(g_ij_neut_dev, geo_host->g_ij_neut);
  gkyl_array_copy(dxdz_dev, geo_host->dxdz);
  gkyl_array_copy(dzdx_dev, geo_host->dzdx);
  gkyl_array_copy(dualmag_dev, geo_host->dualmag);
  gkyl_array_copy(normals_dev, geo_host->normals);
  gkyl_array_copy(jacobgeo_dev , geo_host->jacobgeo);
  gkyl_array_copy(jacobgeo_inv_dev, geo_host->jacobgeo_inv);
  gkyl_array_copy(gij_dev, geo_host->gij);
  gkyl_array_copy(gij_neut_dev, geo_host->gij_neut);
  gkyl_array_copy(b_i_dev, geo_host->b_i);
  gkyl_array_copy(bcart_dev, geo_host->bcart);
  gkyl_array_copy(cmag_dev, geo_host->cmag);
  gkyl_array_copy(jacobtot_dev, geo_host->jacobtot);
  gkyl_array_copy(jacobtot_inv_dev, geo_host->jacobtot_inv);
  gkyl_array_copy(bmag_inv_dev, geo_host->bmag_inv);
  gkyl_array_copy(bmag_inv_sq_dev, geo_host->bmag_inv_sq);
  gkyl_array_copy(gxxj_dev, geo_host->gxxj);
  gkyl_array_copy(gxyj_dev, geo_host->gxyj);
  gkyl_array_copy(gyyj_dev, geo_host->gyyj);
  gkyl_array_copy(gxzj_dev, geo_host->gxzj);
  gkyl_array_copy(eps2_dev, geo_host->eps2);

  // this is for the memcpy below
  up->mc2p = mc2p_dev->on_dev;
  up->mc2p_deflated = mc2p_delfated_dev->on_dev;
  up->mc2nu_pos = mc2nu_pos_dev->on_dev;
  up->mc2nu_pos_deflated = mc2nu_pos_deflated_dev->on_dev;
  up->bmag = bmag_dev->on_dev;
  up->g_ij = g_ij_dev->on_dev;
  up->g_ij_neut = g_ij_neut_dev->on_dev;
  up->dxdz = dxdz_dev->on_dev;
  up->dzdx = dzdx_dev->on_dev;
  up->dualmag = dualmag_dev->on_dev;
  up->normals = normals_dev->on_dev;
  up->jacobgeo = jacobgeo_dev->on_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev->on_dev;
  up->gij = gij_dev->on_dev;
  up->gij_neut = gij_neut_dev->on_dev;
  up->b_i = b_i_dev->on_dev;
  up->bcart = bcart_dev->on_dev;
  up->cmag = cmag_dev->on_dev;
  up->jacobtot = jacobtot_dev->on_dev;
  up->jacobtot_inv = jacobtot_inv_dev->on_dev;
  up->bmag_inv = bmag_inv_dev->on_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev->on_dev;
  up->gxxj = gxxj_dev->on_dev;
  up->gxyj = gxyj_dev->on_dev;
  up->gyyj = gyyj_dev->on_dev;
  up->gxzj = gxzj_dev->on_dev;
  up->eps2 = eps2_dev->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);

  // Initialize the device geometry object
  struct gk_geometry *up_cu = (struct gk_geometry*) gkyl_cu_malloc(sizeof(struct gk_geometry));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gk_geometry), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;

  // geometry object should store host pointer
  up->mc2p = mc2p_dev;
  up->mc2p_deflated = mc2p_delfated_dev;
  up->mc2nu_pos  = mc2nu_pos_dev;
  up->mc2nu_pos_deflated = mc2nu_pos_deflated_dev;
  up->bmag = bmag_dev;
  up->g_ij = g_ij_dev;
  up->g_ij_neut = g_ij_neut_dev;
  up->dxdz = dxdz_dev;
  up->dzdx = dzdx_dev;
  up->dualmag = dualmag_dev;
  up->normals = normals_dev;
  up->jacobgeo = jacobgeo_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev;
  up->gij = gij_dev;
  up->gij_neut = gij_neut_dev;
  up->b_i = b_i_dev;
  up->bcart = bcart_dev;
  up->cmag = cmag_dev;
  up->jacobtot = jacobtot_dev;
  up->jacobtot_inv = jacobtot_inv_dev;
  up->bmag_inv = bmag_inv_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev;
  up->gxxj = gxxj_dev;
  up->gxyj = gxyj_dev;
  up->gyyj = gyyj_dev;
  up->gxzj = gxzj_dev;
  up->eps2 = eps2_dev;
  
  return up;
}

