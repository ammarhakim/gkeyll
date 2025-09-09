#include <gkyl_array.h>
#include <gkyl_mirror_geo_dg.h>
#include <gkyl_gk_geometry.h>

void
gk_geometry_mirror_array_acquire(struct gk_geometry *up, struct gkyl_mirror_geo_dg *geo_dg)
{
    up->mc2p = gkyl_array_acquire(geo_dg->mapc2p);
    up->mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->mc2nu_pos = gkyl_array_acquire(geo_dg->mc2nu_pos);
    up->mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->dxdz = gkyl_array_acquire(geo_dg->tang);
    up->dzdx = gkyl_array_acquire(geo_dg->dual);
    up->dualmag = gkyl_array_acquire(geo_dg->dualmag);
    up->normals = gkyl_array_acquire(geo_dg->normals);
    up->g_ij = gkyl_array_acquire(geo_dg->metric_covar);
    up->g_ij_neut = gkyl_array_acquire(geo_dg->metric_covar_neut);
    up->gij = gkyl_array_acquire(geo_dg->metric_contr);
    up->gij_neut = gkyl_array_acquire(geo_dg->metric_contr_neut);
    up->gxxj = gkyl_array_acquire(geo_dg->gxxj);
    up->gxyj = gkyl_array_acquire(geo_dg->gxyj);
    up->gyyj = gkyl_array_acquire(geo_dg->gyyj);
    up->gxzj = gkyl_array_acquire(geo_dg->gxzj);
    up->jacobgeo = gkyl_array_acquire(geo_dg->Jc);
    up->jacobgeo_inv = gkyl_array_acquire(geo_dg->Jc_inv);
    up->jacobtot = gkyl_array_acquire(geo_dg->JB);
    up->jacobtot_inv = gkyl_array_acquire(geo_dg->JB_inv);
    up->b_i = gkyl_array_acquire(geo_dg->b_covar);
    up->bcart = gkyl_array_acquire(geo_dg->b_cart);
    up->bmag = gkyl_array_acquire(geo_dg->Bmag);
    up->bmag_inv = gkyl_array_acquire(geo_dg->Bmag_inv);
    up->bmag_inv_sq = gkyl_array_acquire(geo_dg->Bmag_inv_sq);
    up->cmag = gkyl_array_acquire(geo_dg->C);
    up->eps2 = gkyl_array_acquire(geo_dg->eps2);
}

void 
gk_geometry_mirror_array_new(struct gk_geometry *up)
{
    up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
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
}

void
gk_geometry_mirror_array_copy(struct gk_geometry *up, struct gkyl_mirror_geo_dg *mirror_geo_dg)
{
    struct gkyl_range global_sub_range;
    int intersect = gkyl_sub_range_intersect(&global_sub_range, &up->global, &up->local);
    gkyl_array_copy_range_to_range(up->mc2p, mirror_geo_dg->mapc2p, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->mc2nu_pos, mirror_geo_dg->mc2nu_pos, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag, mirror_geo_dg->Bmag, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->g_ij, mirror_geo_dg->metric_covar, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->g_ij_neut, mirror_geo_dg->metric_covar_neut, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->dxdz, mirror_geo_dg->tang, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->dzdx, mirror_geo_dg->dual, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->dualmag, mirror_geo_dg->dualmag, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->normals, mirror_geo_dg->normals, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobgeo, mirror_geo_dg->Jc, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobgeo_inv, mirror_geo_dg->Jc_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gij, mirror_geo_dg->metric_contr, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gij_neut, mirror_geo_dg->metric_contr_neut, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->b_i, mirror_geo_dg->b_covar, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->bcart, mirror_geo_dg->b_cart, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->cmag, mirror_geo_dg->C, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobtot, mirror_geo_dg->JB, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobtot_inv, mirror_geo_dg->JB_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag_inv, mirror_geo_dg->Bmag_inv, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag_inv_sq, mirror_geo_dg->Bmag_inv_sq, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gxxj, mirror_geo_dg->gxxj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gxyj, mirror_geo_dg->gxyj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gyyj, mirror_geo_dg->gyyj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->gxzj, mirror_geo_dg->gxzj, &up->local, &global_sub_range);
    gkyl_array_copy_range_to_range(up->eps2, mirror_geo_dg->eps2, &up->local, &global_sub_range);
}