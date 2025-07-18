#include <gkyl_alloc.h>
#include <gkyl_math.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_geo_gen *
gkyl_mirror_geo_gen_inew(const struct gkyl_mirror_geo_gen_inp *inp)
{
  struct gkyl_mirror_geo_gen *geo = gkyl_malloc(sizeof *geo);

  // Create nodal ranges. Must extend to 3x for cartesian projection
  struct gkyl_range nodal_range_2x, nodal_range_3x;
  int nodes_2x[2];
  if (inp->basis.poly_order == 1) {
    nodes_2x[0] = gkyl_range_shape(&inp->range, 0) + 1;
    nodes_2x[1] = gkyl_range_shape(&inp->range, 2) + 1;
  }
  if (inp->basis.poly_order == 2) {
    nodes_2x[0] = 2*gkyl_range_shape(&inp->range, 0) + 1;
    nodes_2x[1] = 2*gkyl_range_shape(&inp->range, 2) + 1;
  }
  gkyl_range_init_from_shape(&nodal_range_2x, 2, nodes_2x);

  int nodes_3x[3];
  if (inp->basis.poly_order == 1) {
    for (int d=0; d<3; ++d)
      nodes_3x[d] = gkyl_range_shape(&inp->range, d) + 1;
  }
  if (inp->basis.poly_order == 2) {
    for (int d=0; d<3; ++d)
      nodes_3x[d] = 2*gkyl_range_shape(&inp->range, d) + 1;
  }
  gkyl_range_init_from_shape(&nodal_range_3x, 3, nodes_3x);

  enum { PSI_IDX, AL_IDX, Z_IDX }; // arrangement of computational coordinates in the grid
  double dalpha = inp->comp_grid->dx[AL_IDX];
  double alpha_lo = inp->comp_grid->lower[AL_IDX] + (inp->range.lower[AL_IDX] - inp->range.lower[AL_IDX]) * inp->comp_grid->dx[AL_IDX];

  // Compute geometry quantities at node locations
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_geo_gen_geom), nodal_range_3x.volume);

  // Loop over all node locations
  for (int ip=nodal_range_3x.lower[PSI_IDX]; ip<=nodal_range_3x.upper[PSI_IDX]; ++ip) {
    for (int ia=nodal_range_3x.lower[AL_IDX]; ia<=nodal_range_3x.upper[AL_IDX]; ++ia) {
      for (int iz=nodal_range_3x.lower[Z_IDX]; iz<=nodal_range_3x.upper[Z_IDX]; ++iz) {
        // Find our location in the array
        int idx_2x[2] = { ip, iz };
        int idx_3x[3] = { ip, ia, iz };
        long loc_2x = gkyl_range_idx(&nodal_range_2x, idx_2x);
        long loc_3x = gkyl_range_idx(&nodal_range_3x, idx_3x);

        // Access relevant array component
        struct gkyl_mirror_geo_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc_3x);
        const struct gkyl_mirror_grid_gen_geom *grid = gkyl_array_fetch(inp->mirror_grid->nodes_geom, loc_2x);
        const double *node_rz = gkyl_array_fetch(inp->mirror_grid->nodes_rz, loc_2x);
        double alpha_curr = alpha_lo + ia*dalpha;
        const double *node_psi = gkyl_array_fetch(inp->mirror_grid->nodes_psi, loc_2x);

        // Copy node locations
        for (int i=0; i<2; ++i)
          g->rza_coord[i] = node_rz[i];
        g->rza_coord[2] = alpha_curr; // z is negative in the mirror geometry
        g->psi = node_psi[0];
        
        // Dual and tangent vectors in cartesian coordinates
        for (int i=0; i<3; ++i) {
            g->tang[i] = gkyl_vec3_polar_con_to_cart(g->rza_coord[0], alpha_curr, grid->tang[i]);
            g->dual[i] = gkyl_vec3_polar_con_to_cart(g->rza_coord[0], alpha_curr, grid->dual[i]);
        }

        // Magnatude of dual vectors
        for (int i=0; i<3; ++i)
          g->dualmag[i] = gkyl_vec3_len(g->dual[i]);

        // Normal vectors
        for (int i=0; i<3; ++i) {
            g->normal[i] = gkyl_vec3_norm(gkyl_vec3_polar_con_to_cart(g->rza_coord[0], alpha_curr, grid->dual[i]));
            g->normal[i] = gkyl_vec3_norm(g->dual[i]);
        }

        // Determine metric tensor
        int count = 0;
        for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
            if (i > j)
              continue;
            g->metric_contr[count] = gkyl_vec3_dot(g->dual[i], g->dual[j]);
            g->metric_covar[count] = gkyl_vec3_dot(g->tang[i], g->tang[j]);
            // For neutrals, we assume the metric is the same for now
            g->metric_contr_neut[count] = g->metric_contr[count];
            g->metric_covar_neut[count] = g->metric_covar[count];
            count++;
          }
        }

          
        // Determine Cartesian components of magnetic field vector
        struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(g->rza_coord[0], alpha_curr, grid->B);

        g->Bmag = gkyl_vec3_len(B_cart);
        g->Bmag_inv = 1.0 / g->Bmag;
        g->Bmag_inv_sq = g->Bmag_inv * g->Bmag_inv;
        g->b_cart = gkyl_vec3_norm(B_cart);
        
        // Determine covariant magnetic field unit vector
        g->b_covar = gkyl_vec3_norm(gkyl_vec3_polar_con_to_cov(g->rza_coord[0], grid->B));

        // Determine scalar quantities combining J and B
        g->Jc = grid->Jc;
        g->Jc_inv = 1.0 / g->Jc;
        g->JB = g->Jc * g->Bmag;
        g->JB_inv = 1.0 / g->JB;

        // Compute cmag
        g->C = g->JB / sqrt(g->metric_covar[5]); // g_33 is the last element in the covariant metric tensor
        // I think eps2 changed because g^33 = g_33 = 1 unambiguously
        g->eps2 = g->Jc * (g->metric_contr[5] - 1 / g->metric_covar[5]);
      }
    }
  }
  return geo;
}

void 
gkyl_mirror_geo_gen_release(struct gkyl_mirror_geo_gen *geo) 
{
  gkyl_array_release(geo->nodes_geom);
  gkyl_free(geo);

}
