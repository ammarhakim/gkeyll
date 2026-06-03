#include <gkyl_fem_parproj_kernels.h> 
 
GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p1_iny_periodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p2_iny_periodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p1_upy_periodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p2_upy_periodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p1_upy_nonperiodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_2x_ser_p2_upy_nonperiodicy(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p1_inz_periodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[6]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p2_inz_periodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[12]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[12]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[15]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[15]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[17]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[17]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[14]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[14]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[16]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[16]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[19]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[19]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p1_upz_periodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[6]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p2_upz_periodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[12]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[12]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[15]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[15]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[17]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[17]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[14]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[14]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[16]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[16]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[19]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[19]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p1_upz_nonperiodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[6]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[6]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[1]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[1]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

}

GKYL_CU_DH void fem_parproj_bias_line_src_3x_ser_p2_upz_nonperiodicz(const int *edge, const int *perp_dirs, double val, long nodeOff, const long *globalIdxs, double *bsrc) 
{ 
  // edge: -1/+1 for lower or upper edge of the cell in each perp_dir.
  // perp_dirs: directions perpendicular to the biased line.
  // val: biasing value.
  // nodeOff: node offset (prob idx * global number of nodes).
  // globalIdxs: global linear index of each basis function/node in current cell.
  // bsrc: global right side source vector.

  if (edge[0] == -1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[0]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[0]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[3]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[3]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[5]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[5]] = val;
  #endif
    }

  if (edge[0] == -1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[12]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[12]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[15]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[15]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[17]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[17]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == -1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[2]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[2]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[4]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[4]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[7]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[7]] = val;
  #endif
    }

  if (edge[0] == 1 && edge[1] == 1) {
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[14]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[14]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[16]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[16]] = val;
  #endif
  #ifdef __CUDA_ARCH__
  atomicExch((unsigned long long int*) &bsrc[nodeOff+globalIdxs[19]],__double_as_longlong(val));
  #else
  bsrc[nodeOff+globalIdxs[19]] = val;
  #endif
    }

}

