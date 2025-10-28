#include <gkyl_fem_parproj_kernels.h> 
 
GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_inx_periodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = parIdx; 

    globalIdxs[1] = parIdx+1; 



}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_inx_periodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 



}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_inx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = parIdx; 

    globalIdxs[1] = parIdx+1; 



}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_inx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 



}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_upx_periodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = parIdx; 

    globalIdxs[1] = 0; 


}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_upx_periodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 0; 


}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p1_upx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = parIdx; 

    globalIdxs[1] = parIdx+1; 


}

GKYL_CU_DH void fem_parproj_local_to_global_1x_ser_p2_upx_nonperiodicx(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 


}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_iny_periodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 


    globalIdxs[3] = 2*parIdx+3; 



}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_iny_periodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 5*parIdx; 

    globalIdxs[1] = 5*parIdx+1; 

    globalIdxs[2] = 5*parIdx+2; 

    globalIdxs[3] = 5*parIdx+3; 

    globalIdxs[4] = 5*parIdx+4; 

    globalIdxs[5] = 5*parIdx+5; 


    globalIdxs[6] = 5*parIdx+6; 


    globalIdxs[7] = 5*parIdx+7; 



}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_iny_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 


    globalIdxs[3] = 2*parIdx+3; 



}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_iny_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 5*parIdx; 

    globalIdxs[1] = 5*parIdx+1; 

    globalIdxs[2] = 5*parIdx+2; 

    globalIdxs[3] = 5*parIdx+3; 

    globalIdxs[4] = 5*parIdx+4; 

    globalIdxs[5] = 5*parIdx+5; 


    globalIdxs[6] = 5*parIdx+6; 


    globalIdxs[7] = 5*parIdx+7; 



}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_upy_periodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 0; 

    globalIdxs[3] = 1; 


}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_upy_periodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 5*parIdx; 

    globalIdxs[1] = 5*parIdx+1; 

    globalIdxs[2] = 5*parIdx+2; 

    globalIdxs[3] = 5*parIdx+3; 

    globalIdxs[4] = 5*parIdx+4; 

    globalIdxs[5] = 0; 

    globalIdxs[6] = 1; 

    globalIdxs[7] = 2; 


}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p1_upy_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 2*parIdx; 

    globalIdxs[1] = 2*parIdx+1; 

    globalIdxs[2] = 2*parIdx+2; 

    globalIdxs[3] = 2*parIdx+3; 


}

GKYL_CU_DH void fem_parproj_local_to_global_2x_ser_p2_upy_nonperiodicy(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 5*parIdx; 

    globalIdxs[1] = 5*parIdx+1; 

    globalIdxs[2] = 5*parIdx+2; 

    globalIdxs[3] = 5*parIdx+3; 

    globalIdxs[4] = 5*parIdx+4; 

    globalIdxs[5] = 5*parIdx+5; 

    globalIdxs[6] = 5*parIdx+6; 

    globalIdxs[7] = 5*parIdx+7; 


}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_inz_periodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 4*parIdx; 

    globalIdxs[1] = 4*parIdx+1; 

    globalIdxs[2] = 4*parIdx+2; 

    globalIdxs[3] = 4*parIdx+3; 

    globalIdxs[4] = 4*parIdx+4; 


    globalIdxs[5] = 4*parIdx+5; 


    globalIdxs[6] = 4*parIdx+6; 


    globalIdxs[7] = 4*parIdx+7; 



}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_inz_periodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 12*parIdx; 

    globalIdxs[1] = 12*parIdx+1; 

    globalIdxs[2] = 12*parIdx+2; 

    globalIdxs[3] = 12*parIdx+3; 

    globalIdxs[4] = 12*parIdx+4; 

    globalIdxs[5] = 12*parIdx+5; 

    globalIdxs[6] = 12*parIdx+6; 

    globalIdxs[7] = 12*parIdx+7; 

    globalIdxs[8] = 12*parIdx+8; 

    globalIdxs[9] = 12*parIdx+9; 

    globalIdxs[10] = 12*parIdx+10; 

    globalIdxs[11] = 12*parIdx+11; 

    globalIdxs[12] = 12*parIdx+12; 


    globalIdxs[13] = 12*parIdx+13; 


    globalIdxs[14] = 12*parIdx+14; 


    globalIdxs[15] = 12*parIdx+15; 


    globalIdxs[16] = 12*parIdx+16; 


    globalIdxs[17] = 12*parIdx+17; 


    globalIdxs[18] = 12*parIdx+18; 


    globalIdxs[19] = 12*parIdx+19; 



}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_inz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 4*parIdx; 

    globalIdxs[1] = 4*parIdx+1; 

    globalIdxs[2] = 4*parIdx+2; 

    globalIdxs[3] = 4*parIdx+3; 

    globalIdxs[4] = 4*parIdx+4; 


    globalIdxs[5] = 4*parIdx+5; 


    globalIdxs[6] = 4*parIdx+6; 


    globalIdxs[7] = 4*parIdx+7; 



}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_inz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 12*parIdx; 

    globalIdxs[1] = 12*parIdx+1; 

    globalIdxs[2] = 12*parIdx+2; 

    globalIdxs[3] = 12*parIdx+3; 

    globalIdxs[4] = 12*parIdx+4; 

    globalIdxs[5] = 12*parIdx+5; 

    globalIdxs[6] = 12*parIdx+6; 

    globalIdxs[7] = 12*parIdx+7; 

    globalIdxs[8] = 12*parIdx+8; 

    globalIdxs[9] = 12*parIdx+9; 

    globalIdxs[10] = 12*parIdx+10; 

    globalIdxs[11] = 12*parIdx+11; 

    globalIdxs[12] = 12*parIdx+12; 


    globalIdxs[13] = 12*parIdx+13; 


    globalIdxs[14] = 12*parIdx+14; 


    globalIdxs[15] = 12*parIdx+15; 


    globalIdxs[16] = 12*parIdx+16; 


    globalIdxs[17] = 12*parIdx+17; 


    globalIdxs[18] = 12*parIdx+18; 


    globalIdxs[19] = 12*parIdx+19; 



}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_upz_periodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 4*parIdx; 

    globalIdxs[1] = 4*parIdx+1; 

    globalIdxs[2] = 4*parIdx+2; 

    globalIdxs[3] = 4*parIdx+3; 

    globalIdxs[4] = 0; 

    globalIdxs[5] = 1; 

    globalIdxs[6] = 2; 

    globalIdxs[7] = 3; 


}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_upz_periodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 12*parIdx; 

    globalIdxs[1] = 12*parIdx+1; 

    globalIdxs[2] = 12*parIdx+2; 

    globalIdxs[3] = 12*parIdx+3; 

    globalIdxs[4] = 12*parIdx+4; 

    globalIdxs[5] = 12*parIdx+5; 

    globalIdxs[6] = 12*parIdx+6; 

    globalIdxs[7] = 12*parIdx+7; 

    globalIdxs[8] = 12*parIdx+8; 

    globalIdxs[9] = 12*parIdx+9; 

    globalIdxs[10] = 12*parIdx+10; 

    globalIdxs[11] = 12*parIdx+11; 

    globalIdxs[12] = 0; 

    globalIdxs[13] = 1; 

    globalIdxs[14] = 2; 

    globalIdxs[15] = 3; 

    globalIdxs[16] = 4; 

    globalIdxs[17] = 5; 

    globalIdxs[18] = 6; 

    globalIdxs[19] = 7; 


}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p1_upz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 4*parIdx; 

    globalIdxs[1] = 4*parIdx+1; 

    globalIdxs[2] = 4*parIdx+2; 

    globalIdxs[3] = 4*parIdx+3; 

    globalIdxs[4] = 4*parIdx+4; 

    globalIdxs[5] = 4*parIdx+5; 

    globalIdxs[6] = 4*parIdx+6; 

    globalIdxs[7] = 4*parIdx+7; 


}

GKYL_CU_DH void fem_parproj_local_to_global_3x_ser_p2_upz_nonperiodicz(int numCellsPar, int parIdx, long *globalIdxs) 
{ 
  // numCellsPar: number of cells in parallel direction.
  // parIdx:     index of current cell in parallel direction.
  // globalIdxs: global linear index of each basis function/node in current cell.

    globalIdxs[0] = 12*parIdx; 

    globalIdxs[1] = 12*parIdx+1; 

    globalIdxs[2] = 12*parIdx+2; 

    globalIdxs[3] = 12*parIdx+3; 

    globalIdxs[4] = 12*parIdx+4; 

    globalIdxs[5] = 12*parIdx+5; 

    globalIdxs[6] = 12*parIdx+6; 

    globalIdxs[7] = 12*parIdx+7; 

    globalIdxs[8] = 12*parIdx+8; 

    globalIdxs[9] = 12*parIdx+9; 

    globalIdxs[10] = 12*parIdx+10; 

    globalIdxs[11] = 12*parIdx+11; 

    globalIdxs[12] = 12*parIdx+12; 

    globalIdxs[13] = 12*parIdx+13; 

    globalIdxs[14] = 12*parIdx+14; 

    globalIdxs[15] = 12*parIdx+15; 

    globalIdxs[16] = 12*parIdx+16; 

    globalIdxs[17] = 12*parIdx+17; 

    globalIdxs[18] = 12*parIdx+18; 

    globalIdxs[19] = 12*parIdx+19; 


}

