#include <gkyl_fem_poisson_kernels.h> 
 
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p1_inx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p2_inx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p1_upx_periodicx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p2_upx_periodicx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p1_upx_nonperiodicx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_1x_ser_p2_upx_nonperiodicx(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_inx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_inx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_periodicx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_periodicx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_nonperiodicx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_nonperiodicx_iny(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_inx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_inx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_inx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_inx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_periodicx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_periodicx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_periodicx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_periodicx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_nonperiodicx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_nonperiodicx_upy_periodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
    }

  }

}
GKYL_CU_DH void fem_poisson_bias_plane_src_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(int edge, int dir, double val, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell.
  // dir: direction perpendicular to the biased plane.
  // val: biasing value.
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (dir == 0) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

  if (dir == 1) {
    if (edge == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[2]] = val;
  #endif
    }

    if (edge == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[6]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[globalIdxs[7]] = val;
  #endif
    }

  }

}
