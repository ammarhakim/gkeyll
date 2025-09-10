#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[1]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[perpOff+globalIdxs[1]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[3]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[perpOff+globalIdxs[3]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[5]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[perpOff+globalIdxs[5]] = bcVals[5];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[7]],__double_as_longlong(bcVals[5]));
  #else
  bsrc[perpOff+globalIdxs[7]] = bcVals[5];
  #endif

}

GKYL_CU_DH void fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs, double *bsrc) 
{ 
  // epsilon: permittivity tensor.
  // rho: right side source.
  // bcVals: values to impose as BCs, i.e. bcVals[off+0]*phi+bcVals[off+1]*epsilon^{ij}*d(phi)/dx_j=bcVals[off+2].
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  double rdx2[2]; 
  rdx2[0] = 2.0/dx[0]; 
  rdx2[1] = 2.0/dx[1]; 

  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[0]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[0]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[1]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[1]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[2]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[2]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]-0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[3]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[3]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]+0.1178511301977579*rho[4]-0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[4]],0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[4]] += 0.06804138174397718*rho[7]-0.1178511301977579*rho[6]-0.1178511301977579*rho[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[5]],(-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[5]] += (-0.06804138174397718*rho[7])-0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]-0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[6]],(-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[6]] += (-0.06804138174397718*rho[7])+0.1178511301977579*rho[6]-0.1178511301977579*rho[5]-0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]-0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif
  #ifdef __CUDA_ARCH__
  atomicAdd(&bsrc[perpOff+globalIdxs[7]],0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0]);
  #else
  bsrc[perpOff+globalIdxs[7]] += 0.06804138174397718*rho[7]+0.1178511301977579*rho[6]+0.1178511301977579*rho[5]+rdx2[0]*bcVals[5]+0.1178511301977579*rho[4]+0.2041241452319316*rho[3]+0.2041241452319316*rho[2]+0.2041241452319316*rho[1]+0.3535533905932738*rho[0];
  #endif

}

