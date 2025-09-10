#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_mirror_geo_dg.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_geo_dg *
gkyl_mirror_geo_dg_inew(const struct gkyl_mirror_geo_dg_inp *inp)
{
  struct gkyl_mirror_geo_dg *geo = gkyl_malloc(sizeof *geo);

  struct gkyl_basis basis = inp->basis;

  // Create nodal range
  struct gkyl_range  nrange;
  int nodes[3];
  if (basis.poly_order == 1) {
    for (int d=0; d<3; ++d)
      nodes[d] = gkyl_range_shape(&inp->range, d) + 1;
  }
  if (basis.poly_order == 2) {
    for (int d=0; d<3; ++d)
      nodes[d] = 2*gkyl_range_shape(&inp->range, d) + 1;
  }
  gkyl_range_init_from_shape(&nrange, 3, nodes);

  // Allocate the arrays inside this geometry object
  geo->mapc2p = gkyl_array_new(GKYL_DOUBLE, basis.num_basis*inp->comp_grid->ndim, inp->range_ext.volume);
  geo->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, basis.num_basis*inp->comp_grid->ndim, inp->range_ext.volume);
  geo->tang = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, inp->range_ext.volume);
  geo->dual = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, inp->range_ext.volume);
  geo->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, inp->range_ext.volume);
  geo->normals = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, inp->range_ext.volume);
  geo->Jc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->Jc_inv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->JB = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->JB_inv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->metric_covar = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, inp->range_ext.volume);
  geo->metric_covar_neut = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, inp->range_ext.volume);
  geo->metric_contr = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, inp->range_ext.volume);
  geo->metric_contr_neut = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, inp->range_ext.volume);
  geo->gxxj = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->gxyj = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->gyyj = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->gxzj = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->Bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->Bmag_inv = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->Bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->b_covar = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, inp->range_ext.volume);
  geo->b_cart = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, inp->range_ext.volume);
  geo->C = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);
  geo->eps2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, inp->range_ext.volume);

  // First, we create nodal arrays for all the geometric quantities
  // This could be split into two functions. Determine nodal arrays, then n2m transform
  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, inp->comp_grid->ndim, inp->range_ext.volume);
  struct gkyl_array* mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, inp->comp_grid->ndim, inp->range_ext.volume);
  struct gkyl_array* tang_nodal = gkyl_array_new(GKYL_DOUBLE, 9, inp->range_ext.volume);
  struct gkyl_array* dual_nodal = gkyl_array_new(GKYL_DOUBLE, 9, inp->range_ext.volume);
  struct gkyl_array* dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, 3, inp->range_ext.volume);
  struct gkyl_array* normals_nodal = gkyl_array_new(GKYL_DOUBLE, 9, inp->range_ext.volume);

  struct gkyl_array* Jc_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* Jc_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* JB_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* JB_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);

  struct gkyl_array* metric_covar_nodal = gkyl_array_new(GKYL_DOUBLE, 6, inp->range_ext.volume);
  struct gkyl_array* metric_covar_neut_nodal = gkyl_array_new(GKYL_DOUBLE, 6, inp->range_ext.volume);
  struct gkyl_array* metric_contr_nodal = gkyl_array_new(GKYL_DOUBLE, 6, inp->range_ext.volume);
  struct gkyl_array* metric_contr_neut_nodal = gkyl_array_new(GKYL_DOUBLE, 6, inp->range_ext.volume);

  struct gkyl_array* gxxj_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* gxyj_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* gyyj_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* gxzj_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);

  struct gkyl_array* Bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* Bmag_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* Bmag_inv_sq_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* B_covar_nodal = gkyl_array_new(GKYL_DOUBLE, 3, inp->range_ext.volume);
  struct gkyl_array* B_cart_nodal = gkyl_array_new(GKYL_DOUBLE, 3, inp->range_ext.volume);

  struct gkyl_array* C_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);
  struct gkyl_array* eps2_nodal = gkyl_array_new(GKYL_DOUBLE, 1, inp->range_ext.volume);

  enum { PSI_IDX, AL_IDX, Z_IDX }; // arrangement of computational coordinates in the grid

  for(int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
    for(int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia) {
      for(int iz=nrange.lower[Z_IDX]; iz<=nrange.upper[Z_IDX]; ++iz) {

        int cidx[3] = { ip, ia, iz };
        long idx = gkyl_range_idx(&nrange, cidx);

        struct gkyl_mirror_geo_gen_geom *geom_n = gkyl_array_fetch(inp->mirror_geo->nodes_geom, idx);

        double *mapc2p_nodal_n = gkyl_array_fetch(mapc2p_nodal, idx);
        mapc2p_nodal_n[0] = geom_n->rza_coord[0];
        mapc2p_nodal_n[1] = geom_n->rza_coord[1];
        mapc2p_nodal_n[2] = -geom_n->rza_coord[2];

        double *mc2nu_pos_nodal_n = gkyl_array_fetch(mc2nu_pos_nodal, idx);
        mc2nu_pos_nodal_n[0] = geom_n->psi;
        mc2nu_pos_nodal_n[1] = geom_n->rza_coord[2];
        mc2nu_pos_nodal_n[2] = geom_n->rza_coord[1];

        // Fetch all array pointers at the beginning
        double *tang_nodal_n = gkyl_array_fetch(tang_nodal, idx);
        double *dual_nodal_n = gkyl_array_fetch(dual_nodal, idx);
        double *dualmag_nodal_n = gkyl_array_fetch(dualmag_nodal, idx);
        double *normals_nodal_n = gkyl_array_fetch(normals_nodal, idx);

        double *Jc_nodal_n = gkyl_array_fetch(Jc_nodal, idx);
        double *Jc_inv_nodal_n = gkyl_array_fetch(Jc_inv_nodal, idx);
        double *JB_nodal_n = gkyl_array_fetch(JB_nodal, idx);
        double *JB_inv_nodal_n = gkyl_array_fetch(JB_inv_nodal, idx);

        double *metric_covar_nodal_n = gkyl_array_fetch(metric_covar_nodal, idx);
        double *metric_covar_neut_nodal_n = gkyl_array_fetch(metric_covar_neut_nodal, idx);
        double *metric_contr_nodal_n = gkyl_array_fetch(metric_contr_nodal, idx);
        double *metric_contr_neut_nodal_n = gkyl_array_fetch(metric_contr_neut_nodal, idx);

        double *gxxj_nodal_n = gkyl_array_fetch(gxxj_nodal, idx);
        double *gxyj_nodal_n = gkyl_array_fetch(gxyj_nodal, idx);
        double *gxzj_nodal_n = gkyl_array_fetch(gxzj_nodal, idx);
        double *gyyj_nodal_n = gkyl_array_fetch(gyyj_nodal, idx);

        double *Bmag_nodal_n = gkyl_array_fetch(Bmag_nodal, idx);
        double *Bmag_inv_nodal_n = gkyl_array_fetch(Bmag_inv_nodal, idx);
        double *Bmag_inv_sq_nodal_n = gkyl_array_fetch(Bmag_inv_sq_nodal, idx);

        double *B_covar_nodal_n = gkyl_array_fetch(B_covar_nodal, idx);
        double *B_cart_nodal_n = gkyl_array_fetch(B_cart_nodal, idx);

        double *C_nodal_n = gkyl_array_fetch(C_nodal, idx);
        double *eps2_nodal_n = gkyl_array_fetch(eps2_nodal, idx);

        // Set all the scalar quantities
        Jc_nodal_n[0] = geom_n->Jc;
        Jc_inv_nodal_n[0] = geom_n->Jc_inv;
        JB_nodal_n[0] = geom_n->JB;
        JB_inv_nodal_n[0] = geom_n->JB_inv;
        Bmag_nodal_n[0] = geom_n->Bmag;
        Bmag_inv_nodal_n[0] = geom_n->Bmag_inv;
        Bmag_inv_sq_nodal_n[0] = geom_n->Bmag_inv_sq;

        gxxj_nodal_n[0] = geom_n->metric_contr[0] * geom_n->Jc;
        gxyj_nodal_n[0] = geom_n->metric_contr[1] * geom_n->Jc;
        gxzj_nodal_n[0] = geom_n->metric_contr[2] * geom_n->Jc;
        gyyj_nodal_n[0] = geom_n->metric_contr[3] * geom_n->Jc;

        C_nodal_n[0] = geom_n->C;
        eps2_nodal_n[0] = geom_n->eps2;

        // Loop over tensor quantities
        // This breaks conventions of functional programming, but I can't figure out a way to assign this data another way without it being disgustingly verbose.
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            tang_nodal_n[i*3 + j] = geom_n->tang[i].x[j];
            dual_nodal_n[i*3 + j] = geom_n->dual[i].x[j];
            normals_nodal_n[i*3 + j] = geom_n->normal[i].x[j];
          }
          dualmag_nodal_n[i] = geom_n->dualmag[i];
          B_covar_nodal_n[i] = geom_n->b_covar.x[i];
          B_cart_nodal_n[i] = geom_n->b_cart.x[i];
        }

        // Set metric tensors
        for (int i = 0; i < 6; ++i) {
          metric_covar_nodal_n[i] = geom_n->metric_covar[i];
          metric_covar_neut_nodal_n[i] = geom_n->metric_covar_neut[i];
          metric_contr_nodal_n[i] = geom_n->metric_contr[i];
          metric_contr_neut_nodal_n[i] = geom_n->metric_contr_neut[i];
        }
      }
    }
  }

  // Using the nodal arrays, we now convert them to the DG basis
  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&basis, inp->comp_grid, false);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 3, mapc2p_nodal, geo->mapc2p);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 3, mc2nu_pos_nodal, geo->mc2nu_pos);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 9, tang_nodal, geo->tang);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 9, dual_nodal, geo->dual);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 3, dualmag_nodal, geo->dualmag);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 9, normals_nodal, geo->normals);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, Jc_nodal, geo->Jc);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, Jc_inv_nodal, geo->Jc_inv);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, JB_nodal, geo->JB);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, JB_inv_nodal, geo->JB_inv);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 6, metric_covar_nodal, geo->metric_covar);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 6, metric_covar_neut_nodal, geo->metric_covar_neut);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 6, metric_contr_nodal, geo->metric_contr);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 6, metric_contr_neut_nodal, geo->metric_contr_neut);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, gxxj_nodal, geo->gxxj);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, gxyj_nodal, geo->gxyj);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, gyyj_nodal, geo->gyyj);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, gxzj_nodal, geo->gxzj);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, Bmag_nodal, geo->Bmag);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, Bmag_inv_nodal, geo->Bmag_inv);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, Bmag_inv_sq_nodal, geo->Bmag_inv_sq);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 3, B_covar_nodal, geo->b_covar);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 3, B_cart_nodal, geo->b_cart);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, C_nodal, geo->C);
  gkyl_nodal_ops_n2m(n2m, &basis, inp->comp_grid, &nrange, &inp->range, 1, eps2_nodal, geo->eps2);
  
  // Release memory for nodal operations
  gkyl_nodal_ops_release(n2m);

  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mc2nu_pos_nodal);
  gkyl_array_release(tang_nodal);
  gkyl_array_release(dual_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(normals_nodal);
  gkyl_array_release(Jc_nodal);
  gkyl_array_release(Jc_inv_nodal);
  gkyl_array_release(JB_nodal);
  gkyl_array_release(JB_inv_nodal);
  gkyl_array_release(metric_covar_nodal);
  gkyl_array_release(metric_covar_neut_nodal);
  gkyl_array_release(metric_contr_nodal);
  gkyl_array_release(metric_contr_neut_nodal);
  gkyl_array_release(gxxj_nodal);
  gkyl_array_release(gxyj_nodal);
  gkyl_array_release(gyyj_nodal);
  gkyl_array_release(gxzj_nodal);
  gkyl_array_release(Bmag_nodal);
  gkyl_array_release(Bmag_inv_nodal);
  gkyl_array_release(Bmag_inv_sq_nodal);
  gkyl_array_release(B_covar_nodal);
  gkyl_array_release(B_cart_nodal);
  gkyl_array_release(C_nodal);
  gkyl_array_release(eps2_nodal);

  return geo;
}

void 
gkyl_mirror_geo_dg_release(struct gkyl_mirror_geo_dg *geo) 
{
  gkyl_array_release(geo->mapc2p);
  gkyl_array_release(geo->mc2nu_pos);
  gkyl_array_release(geo->tang);
  gkyl_array_release(geo->dual);
  gkyl_array_release(geo->dualmag);
  gkyl_array_release(geo->normals);
  gkyl_array_release(geo->Jc);
  gkyl_array_release(geo->Jc_inv);
  gkyl_array_release(geo->JB);
  gkyl_array_release(geo->JB_inv);
  gkyl_array_release(geo->metric_covar);
  gkyl_array_release(geo->metric_covar_neut);
  gkyl_array_release(geo->metric_contr);
  gkyl_array_release(geo->metric_contr_neut);
  gkyl_array_release(geo->gxxj);
  gkyl_array_release(geo->gxyj);
  gkyl_array_release(geo->gyyj);
  gkyl_array_release(geo->gxzj);
  gkyl_array_release(geo->Bmag);
  gkyl_array_release(geo->Bmag_inv);
  gkyl_array_release(geo->Bmag_inv_sq);
  gkyl_array_release(geo->b_covar);
  gkyl_array_release(geo->b_cart);
  gkyl_array_release(geo->C);
  gkyl_array_release(geo->eps2);
  // Free the geometry object itself
  gkyl_free(geo);

}
