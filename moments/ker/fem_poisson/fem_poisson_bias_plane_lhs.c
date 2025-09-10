#include <gkyl_fem_poisson_kernels.h> 
 
void fem_poisson_bias_plane_lhs_1x_ser_p1_inx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_1x_ser_p2_inx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_1x_ser_p1_upx_periodicx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_1x_ser_p2_upx_periodicx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_1x_ser_p1_upx_nonperiodicx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_1x_ser_p2_upx_nonperiodicx(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_inx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_inx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_periodicx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_periodicx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_nonperiodicx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_nonperiodicx_iny(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_inx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_inx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_inx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_inx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_periodicx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_periodicx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_periodicx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_periodicx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_nonperiodicx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_nonperiodicx_upy_periodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
void fem_poisson_bias_plane_lhs_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
    }

  }

}
void fem_poisson_bias_plane_lhs_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(int edge, int dir, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  if (dir == 0) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[3], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[3], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[4], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[4], globalIdxs[7], 0.0);
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

  if (dir == 1) {
    if (edge == -1) {
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[0], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[0], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[1], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[1], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[2], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[2], globalIdxs[7], 0.0);
    }

    if (edge == 1) {
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[5], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[6], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[5], globalIdxs[7], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[0], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[1], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[2], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[3], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[4], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[5], 0.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[6], 1.0);
  gkyl_mat_triples_insert(tri, globalIdxs[6], globalIdxs[7], 0.0);
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

}
