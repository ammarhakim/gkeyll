#include <gkyl_fem_poisson_perp_kernels.h> 
 
void fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx(const double *epsilon, const double *kSq, const double *dx, const double *bcVals, const long *globalIdxs, gkyl_mat_triples *tri) 
{ 
  // epsilon: permittivity tensor.
  // kSq: wave number squared (factor multiplying phi).
  // dx: cell length in each direction.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // globalIdxs: global linear index of each basis function/node in current cell.
  // tri: triples object (i,j,val), i.e. contribute val to i,j element of the global matrix.

  const double *epsxx = &epsilon[0];

  double rdx00 = 4.0/(dx[0]*dx[0]);

  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[0], -(0.1443375672974064*epsxx[2]*rdx00)+0.16666666666666666*epsxx[0]*rdx00-0.16666666666666666*kSq[3]+0.19245008972987523*kSq[2]+0.19245008972987523*kSq[1]-0.2222222222222222*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[1], 0.1443375672974064*epsxx[2]*rdx00-0.16666666666666666*epsxx[0]*rdx00+0.09622504486493762*kSq[2]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[2], 0.08333333333333333*epsxx[0]*rdx00+0.09622504486493762*kSq[1]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[0], globalIdxs[3], -(0.08333333333333333*epsxx[0]*rdx00)-0.05555555555555555*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[0], 0.1443375672974064*epsxx[2]*rdx00-0.16666666666666666*epsxx[0]*rdx00+0.09622504486493762*kSq[2]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[1], -(0.1443375672974064*epsxx[2]*rdx00)+0.16666666666666666*epsxx[0]*rdx00+0.16666666666666666*kSq[3]+0.19245008972987523*kSq[2]-0.19245008972987523*kSq[1]-0.2222222222222222*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[2], -(0.08333333333333333*epsxx[0]*rdx00)-0.05555555555555555*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[1], globalIdxs[3], 0.08333333333333333*epsxx[0]*rdx00-0.09622504486493762*kSq[1]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[0], 0.08333333333333333*epsxx[0]*rdx00+0.09622504486493762*kSq[1]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[1], -(0.08333333333333333*epsxx[0]*rdx00)-0.05555555555555555*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[2], 0.1443375672974064*epsxx[2]*rdx00+0.16666666666666666*epsxx[0]*rdx00+0.16666666666666666*kSq[3]-0.19245008972987523*kSq[2]+0.19245008972987523*kSq[1]-0.2222222222222222*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[2], globalIdxs[3], -(0.1443375672974064*epsxx[2]*rdx00)-0.16666666666666666*epsxx[0]*rdx00-0.09622504486493762*kSq[2]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[0], -(0.08333333333333333*epsxx[0]*rdx00)-0.05555555555555555*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[1], 0.08333333333333333*epsxx[0]*rdx00-0.09622504486493762*kSq[1]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[2], -(0.1443375672974064*epsxx[2]*rdx00)-0.16666666666666666*epsxx[0]*rdx00-0.09622504486493762*kSq[2]-0.1111111111111111*kSq[0]);
  gkyl_mat_triples_accum(tri, globalIdxs[3], globalIdxs[3], 0.1443375672974064*epsxx[2]*rdx00+0.16666666666666666*epsxx[0]*rdx00-0.16666666666666666*kSq[3]-0.19245008972987523*kSq[2]-0.19245008972987523*kSq[1]-0.2222222222222222*kSq[0]);

}

