#include <gkyl_fem_parproj_couplex_kernels.h> 
 
GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_inx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = idx_curr[0]+1; 

  }  else {
    globalIdxs[1] = idx_curr[0]+1; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_inx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = 2*idx_curr[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx_curr[0]+2; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_inx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = idx_curr[0]+1; 

  }  else {
    globalIdxs[1] = idx_curr[0]+1; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_inx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = 2*idx_curr[0]+2; 

  }  else {
    globalIdxs[2] = 2*idx_curr[0]+2; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_upx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]; 

    globalIdxs[1] = 0; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_upx_periodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]+1; 

    globalIdxs[2] = 0; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p1_upx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]; 

    globalIdxs[1] = idx_curr[0]+1; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_1x_ser_p2_upx_nonperiodicx(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]+1; 

    globalIdxs[2] = 2*idx_curr[0]+2; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+idx_curr[1]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]; 

  }  else {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+idx_curr[1]+1; 

  }  else {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+idx_curr[1]+1; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+2; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+2; 

  }  else if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+1; 

  }  else {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+1; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]; 

  }  else {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]; 

  }

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+3; 

  }  else {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+3; 

  }  else {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+3; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+4; 

  }  else {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+4; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+5; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+5; 

  }  else if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+3; 

  }  else {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+3; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+idx_curr[1]+idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+idx_curr[0]+1; 

  }  else {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+idx_curr[0]+1; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+idx_curr[1]+idx_curr[0]+1; 

  }  else {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+idx_curr[1]+idx_curr[0]+1; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 

  }  else if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+idx_curr[0]+2; 

  }  else {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+idx_curr[0]+2; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+2; 

  }

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+3; 

  }  else {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+3; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 

  }  else if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+5; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+2*idx_curr[1]; 

    globalIdxs[1] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+1; 

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+2; 

  }  else {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+2; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+3; 

  }  else {
    globalIdxs[3] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+3; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+1; 

    globalIdxs[2] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2; 

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+3; 

    globalIdxs[4] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+4; 

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+5; 

  }  else {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+5; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+6; 

  }  else {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+6; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+7; 

  }  else {
    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+7; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]; 

    globalIdxs[1] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+1; 

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+2; 

  }  else {
    globalIdxs[2] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+2; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[3] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 

  }  else {
    globalIdxs[3] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+1; 

    globalIdxs[2] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+2; 

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+3; 

    globalIdxs[4] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+4; 

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+5; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+6; 

  }  else {
    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+6; 

  }

   if ((idx_curr[1]+1==num_cells[1]-1)) {
    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 

  }  else {
    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+idx_curr[1]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]; 

  }  else {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]; 

  }

    globalIdxs[2] = idx_curr[0]*num_cells[1]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]; 

  }  else {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]; 

  }  else {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]; 

  }

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+3; 

  }  else {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2; 

  }

    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]; 

    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]; 

  }  else {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+idx_curr[1]+idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+idx_curr[0]+1; 

  }  else {
    globalIdxs[1] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+idx_curr[0]+1; 

  }

    globalIdxs[2] = idx_curr[0]*num_cells[1]+idx_curr[1]+idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 

  }  else {
    globalIdxs[3] = (idx_curr[0]+1)*num_cells[1]+idx_curr[1]+idx_curr[0]+2; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[2] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+2; 

  }

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[4] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+4; 

  }

    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+3; 

    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 

  }  else {
    globalIdxs[7] = (3*idx_curr[0]+3)*num_cells[1]+3*idx_curr[1]+2*idx_curr[0]+5; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+2*idx_curr[1]; 

    globalIdxs[1] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+1; 

    globalIdxs[2] = idx_curr[0]*num_cells[1]; 

    globalIdxs[3] = idx_curr[0]*num_cells[1]+1; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+1; 

    globalIdxs[2] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2; 

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+3; 

    globalIdxs[4] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+4; 

    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]; 

    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+1; 

    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+2; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]; 

    globalIdxs[1] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+1; 

    globalIdxs[2] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+2; 

    globalIdxs[3] = idx_curr[0]*num_cells[1]+2*idx_curr[1]+idx_curr[0]+3; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]; 

    globalIdxs[1] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+1; 

    globalIdxs[2] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+2; 

    globalIdxs[3] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+3; 

    globalIdxs[4] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+4; 

    globalIdxs[5] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+5; 

    globalIdxs[6] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+6; 

    globalIdxs[7] = 3*idx_curr[0]*num_cells[1]+5*idx_curr[1]+2*idx_curr[0]+7; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }  else {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+8; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+8; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+5; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

  }

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+6; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }  else {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+12; 

  }  else {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+12; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+9; 

  }  else {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+9; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+20; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+20; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+14; 

  }  else {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+14; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }  else {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+22; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+22; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+15; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }  else {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+8; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+8; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+5; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

  }

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+6; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }  else {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+12; 

  }  else {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+12; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+9; 

  }  else {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+9; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+20; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+20; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+14; 

  }  else {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+14; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }  else {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

   if ((idx_curr[0]+1==num_cells[0]-1) && (idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+22; 

  }  else if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+22; 

  }  else if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+15; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+1; 

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

    globalIdxs[3] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+3; 

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+5; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+7; 

  }  else {
    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+7; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+1; 

    globalIdxs[2] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+2; 

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[4] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+4; 

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+6; 

    globalIdxs[7] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+7; 

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

    globalIdxs[9] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+9; 

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

    globalIdxs[11] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+11; 

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+12; 

  }  else {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+12; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+14; 

  }  else {
    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+14; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+16; 

  }  else {
    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+16; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+18; 

  }  else {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+18; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+19; 

  }  else {
    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+19; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+1; 

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

    globalIdxs[3] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+3; 

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+5; 

  }  else {
    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+5; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else {
    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+7; 

  }  else {
    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+7; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+1; 

    globalIdxs[2] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+2; 

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[4] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+4; 

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+6; 

    globalIdxs[7] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+7; 

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

    globalIdxs[9] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+9; 

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

    globalIdxs[11] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+11; 

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+12; 

  }  else {
    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+12; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+14; 

  }  else {
    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+14; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+16; 

  }  else {
    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+16; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else {
    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+18; 

  }  else {
    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+18; 

  }

   if ((idx_curr[2]+1==num_cells[2]-1)) {
    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+19; 

  }  else {
    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+19; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[0]+2; 

  }

    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[0]+3; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

  }

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+6; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]; 

    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+5; 

  }  else {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+5; 

  }

    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+7; 

  }

    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+3; 

    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+5*idx_curr[0]+8; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

  }  else {
    globalIdxs[1] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

  }

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

  }  else {
    globalIdxs[3] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

  }

    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

  }  else {
    globalIdxs[5] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+4; 

  }

    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+3; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+8; 

  }  else {
    globalIdxs[7] = (2*idx_curr[0]+2)*num_cells[2]+2*idx_curr[2]+2*idx_curr[0]+5; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+1; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

  }  else {
    globalIdxs[2] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

  }

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+2; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

  }  else {
    globalIdxs[4] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

  }

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+4; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

  }  else {
    globalIdxs[7] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

  }

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+5; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

  }  else {
    globalIdxs[9] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

  }

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+6; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

  }  else {
    globalIdxs[11] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

  }

    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+7; 

    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+8; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

  }  else {
    globalIdxs[14] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+12; 

  }

    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+9; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+20; 

  }  else {
    globalIdxs[16] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+14; 

  }

    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+10; 

    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+11; 

   if ((idx_curr[0]+1==num_cells[0]-1)) {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+22; 

  }  else {
    globalIdxs[19] = (7*idx_curr[0]+7)*num_cells[2]+7*idx_curr[2]+5*idx_curr[0]+15; 

  }


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+1; 

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

    globalIdxs[3] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+3; 

    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]; 

    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]+1; 

    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]+2; 

    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+2*idx_curr[0]+3; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_periodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+1; 

    globalIdxs[2] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+2; 

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[4] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+4; 

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+6; 

    globalIdxs[7] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+7; 

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

    globalIdxs[9] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+9; 

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

    globalIdxs[11] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+11; 

    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]; 

    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+1; 

    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+2; 

    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+3; 

    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+4; 

    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+5; 

    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+6; 

    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+5*idx_curr[0]+7; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]; 

    globalIdxs[1] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+1; 

    globalIdxs[2] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+2; 

    globalIdxs[3] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+3; 

    globalIdxs[4] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+4; 

    globalIdxs[5] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+5; 

    globalIdxs[6] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+6; 

    globalIdxs[7] = 2*idx_curr[0]*num_cells[2]+4*idx_curr[2]+2*idx_curr[0]+7; 


}

GKYL_CU_DH void fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_nonperiodicz(const int *num_cells, const int *idx_curr, long *globalIdxs) 
{ 
  // num_cells:  number of cells in each direction.
  // idx_curr:   index of current cell in each direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]; 

    globalIdxs[1] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+1; 

    globalIdxs[2] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+2; 

    globalIdxs[3] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+3; 

    globalIdxs[4] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+4; 

    globalIdxs[5] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+5; 

    globalIdxs[6] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+6; 

    globalIdxs[7] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+7; 

    globalIdxs[8] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+8; 

    globalIdxs[9] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+9; 

    globalIdxs[10] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+10; 

    globalIdxs[11] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+11; 

    globalIdxs[12] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+12; 

    globalIdxs[13] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+13; 

    globalIdxs[14] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+14; 

    globalIdxs[15] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+15; 

    globalIdxs[16] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+16; 

    globalIdxs[17] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+17; 

    globalIdxs[18] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+18; 

    globalIdxs[19] = 7*idx_curr[0]*num_cells[2]+12*idx_curr[2]+5*idx_curr[0]+19; 


}

