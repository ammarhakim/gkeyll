#include <gkyl_fem_poisson_kernels.h> 
 
GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p1_inx_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = idx[0]+1; 

  }  else {
    globalIdxs[1] = idx[0]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p2_inx_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]; 

    globalIdxs[1] = 2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = 2*idx[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p1_inx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = idx[0]+1; 

  }  else {
    globalIdxs[1] = idx[0]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p2_inx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]; 

    globalIdxs[1] = 2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = 2*idx[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p1_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]; 

    globalIdxs[1] = 0; 


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p2_upx_periodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]; 

    globalIdxs[1] = 2*idx[0]+1; 

    globalIdxs[2] = 0; 


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p1_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]; 

    globalIdxs[1] = idx[0]+1; 


}

GKYL_CU_DH void fem_poisson_local_to_global_1x_ser_p2_upx_nonperiodicx(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx[0]; 

    globalIdxs[1] = 2*idx[0]+1; 

    globalIdxs[2] = 2*idx[0]+2; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+2; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+2; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+5*idx[1]; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+5*idx[1]+3; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+5*idx[1]+5; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+5*idx[1]+5; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+3; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+3; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

   if ((idx[0]+1==numCells[0]-1) && (idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+7; 

  }  else if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+7; 

  }  else if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

    globalIdxs[1] = idx[1]; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = idx[1]+1; 

  }  else {
    globalIdxs[3] = idx[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

    globalIdxs[2] = 3*idx[1]; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

    globalIdxs[4] = 3*idx[1]+2; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 3*idx[1]+3; 

  }  else {
    globalIdxs[7] = 3*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

    globalIdxs[1] = idx[1]; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = idx[1]+1; 

  }  else {
    globalIdxs[3] = idx[1]+1; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

    globalIdxs[2] = 3*idx[1]; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

    globalIdxs[4] = 3*idx[1]+2; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 3*idx[1]+3; 

  }  else {
    globalIdxs[7] = 3*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]; 

    globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+2; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+3; 

  }  else {
    globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+5*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+5*idx[1]+1; 

    globalIdxs[2] = 3*idx[0]*numCells[1]+5*idx[1]+2; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+5*idx[1]+3; 

    globalIdxs[4] = 3*idx[0]*numCells[1]+5*idx[1]+4; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+5*idx[1]+5; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+5*idx[1]+5; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+5*idx[1]+6; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+5*idx[1]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 3*idx[0]*numCells[1]+5*idx[1]+7; 

  }  else {
    globalIdxs[7] = 3*idx[0]*numCells[1]+5*idx[1]+7; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]+idx[0]; 

    globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+idx[0]+1; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 

  }  else {
    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+idx[0]+3; 

  }  else {
    globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+idx[0]+3; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+1; 

    globalIdxs[2] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+2; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+3; 

    globalIdxs[4] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+4; 

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[5] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[5] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+5; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[6] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+6; 

  }  else {
    globalIdxs[6] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+6; 

  }

   if ((idx[1]+1==numCells[1]-1)) {
    globalIdxs[7] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+7; 

  }  else {
    globalIdxs[7] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+7; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }

    globalIdxs[2] = idx[0]*numCells[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }

    globalIdxs[5] = 3*idx[0]*numCells[1]; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }

    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]; 

  }

    globalIdxs[2] = idx[0]*numCells[1]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+5*idx[1]; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+5*idx[1]+3; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2; 

  }

    globalIdxs[5] = 3*idx[0]*numCells[1]; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+1; 

  }  else {
    globalIdxs[1] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+1; 

  }

    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+2*idx[1]+idx[0]+3; 

  }  else {
    globalIdxs[3] = (idx[0]+1)*numCells[1]+idx[1]+idx[0]+2; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+2; 

  }

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+5; 

  }  else {
    globalIdxs[4] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+4; 

  }

    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

   if ((idx[0]+1==numCells[0]-1)) {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+5*idx[1]+2*idx[0]+7; 

  }  else {
    globalIdxs[7] = (3*idx[0]+3)*numCells[1]+3*idx[1]+2*idx[0]+5; 

  }


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]; 

    globalIdxs[1] = idx[1]; 

    globalIdxs[2] = idx[0]*numCells[1]; 

    globalIdxs[3] = 0; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+1; 

    globalIdxs[2] = 3*idx[1]; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2; 

    globalIdxs[4] = 3*idx[1]+2; 

    globalIdxs[5] = 3*idx[0]*numCells[1]; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+1; 

    globalIdxs[7] = 0; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+idx[1]+idx[0]; 

    globalIdxs[1] = idx[1]; 

    globalIdxs[2] = idx[0]*numCells[1]+idx[1]+idx[0]+1; 

    globalIdxs[3] = idx[1]+1; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_periodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+1; 

    globalIdxs[2] = 3*idx[1]; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+2; 

    globalIdxs[4] = 3*idx[1]+2; 

    globalIdxs[5] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+3; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+3*idx[1]+2*idx[0]+4; 

    globalIdxs[7] = 3*idx[1]+3; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]; 

    globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+1; 

    globalIdxs[2] = idx[0]*numCells[1]; 

    globalIdxs[3] = idx[0]*numCells[1]+1; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+5*idx[1]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+5*idx[1]+1; 

    globalIdxs[2] = 3*idx[0]*numCells[1]+5*idx[1]+2; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+5*idx[1]+3; 

    globalIdxs[4] = 3*idx[0]*numCells[1]+5*idx[1]+4; 

    globalIdxs[5] = 3*idx[0]*numCells[1]; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+1; 

    globalIdxs[7] = 3*idx[0]*numCells[1]+2; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx[0]*numCells[1]+2*idx[1]+idx[0]; 

    globalIdxs[1] = idx[0]*numCells[1]+2*idx[1]+idx[0]+1; 

    globalIdxs[2] = idx[0]*numCells[1]+2*idx[1]+idx[0]+2; 

    globalIdxs[3] = idx[0]*numCells[1]+2*idx[1]+idx[0]+3; 


}

GKYL_CU_DH void fem_poisson_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(const int *numCells, const int *idx, long *globalIdxs) 
{ 
  // numCells:   number of cells in each direction.
  // idx:        multi-dimensional index of current cell.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]; 

    globalIdxs[1] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+1; 

    globalIdxs[2] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+2; 

    globalIdxs[3] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+3; 

    globalIdxs[4] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+4; 

    globalIdxs[5] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+5; 

    globalIdxs[6] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+6; 

    globalIdxs[7] = 3*idx[0]*numCells[1]+5*idx[1]+2*idx[0]+7; 


}

