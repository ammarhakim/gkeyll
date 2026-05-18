#include <gkyl_fem_parproj_kernels.h> 
 
void fem_parproj_bias_line_lhs_2x_ser_p1_iny_periodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  }

}

void fem_parproj_bias_line_lhs_2x_ser_p2_iny_periodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
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

void fem_parproj_bias_line_lhs_2x_ser_p1_upy_periodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  }

}

void fem_parproj_bias_line_lhs_2x_ser_p2_upy_periodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
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

void fem_parproj_bias_line_lhs_2x_ser_p1_upy_nonperiodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  }

}

void fem_parproj_bias_line_lhs_2x_ser_p2_upy_nonperiodicy(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
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

void fem_parproj_bias_line_lhs_3x_ser_p1_inz_periodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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

void fem_parproj_bias_line_lhs_3x_ser_p2_inz_periodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[19], 0.0);
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[12], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[15], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[17], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[14], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[16], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[19], 1.0);
  }

}

void fem_parproj_bias_line_lhs_3x_ser_p1_upz_periodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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

void fem_parproj_bias_line_lhs_3x_ser_p2_upz_periodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[19], 0.0);
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[12], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[15], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[17], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[14], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[16], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[19], 1.0);
  }

}

void fem_parproj_bias_line_lhs_3x_ser_p1_upz_nonperiodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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

void fem_parproj_bias_line_lhs_3x_ser_p2_upz_nonperiodicz(const int *edge, const int *perp_dir, const long *globalIdxs, gkyl_mat_triples *tri) 
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
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[19], 0.0);
  }

  if (edge[0] == -1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[12], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[12], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[15], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[15], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[17], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[17], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[7], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[7], globalIdxs[19], 0.0);
  }

  if (edge[0] == 1 && edge[1] == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[14], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[14], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[16], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[16], globalIdxs[19], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[8], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[9], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[10], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[11], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[12], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[13], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[14], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[15], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[16], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[17], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[18], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[19], globalIdxs[19], 1.0);
  }

}

