#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[0]],0.16666666666666666*rho[3]-0.28867513459481287*rho[2]-0.28867513459481287*rho[1]+0.5*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[0]] += 0.16666666666666666*rho[3]-0.28867513459481287*rho[2]-0.28867513459481287*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[1]],-(0.16666666666666666*rho[3])-0.28867513459481287*rho[2]+0.28867513459481287*rho[1]+0.5*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[1]] += -(0.16666666666666666*rho[3])-0.28867513459481287*rho[2]+0.28867513459481287*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[2]],-(0.16666666666666669*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[2]] += -(0.16666666666666669*rho[3])+0.2886751345948129*rho[2]-0.2886751345948129*rho[1]+0.5*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[3]],0.16666666666666669*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[3]] += 0.16666666666666669*rho[3]+0.2886751345948129*rho[2]+0.2886751345948129*rho[1]+0.5*rho[0];
  #endif

}

