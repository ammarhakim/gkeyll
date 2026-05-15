#include <gkyl_fem_poisson_perp_kernels.h> 
 
void fem_poisson_perp_bias_line_lhs_3x_ser_p1_upx_periodicx_iny(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (edge[0] == -1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);
  }

}

void fem_poisson_perp_bias_line_lhs_3x_ser_p1_upx_nonperiodicx_iny(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (edge[0] == -1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);
  }

}

