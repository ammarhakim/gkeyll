#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_bias_line_src_3x_ser_p1_inx_iny(const int *edge, const int *perp_dirs, double val, long perpOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[6]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[3]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[perpOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[perpOff+globalIdxs[7]] = val;
  #endif
    }

}

