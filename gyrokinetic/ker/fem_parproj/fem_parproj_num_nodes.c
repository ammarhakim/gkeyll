#include <gkyl_fem_parproj_kernels.h> 
 
long fem_parproj_num_nodes_global_1x_ser_p1_periodicx(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return numCellsPar;
}

long fem_parproj_num_nodes_global_1x_ser_p2_periodicx(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+2;
}

long fem_parproj_num_nodes_global_1x_ser_p1_nonperiodicx(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return numCellsPar+1;
}

long fem_parproj_num_nodes_global_1x_ser_p2_nonperiodicx(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+3;
}

long fem_parproj_num_nodes_global_2x_ser_p1_periodicy(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+2;
}

long fem_parproj_num_nodes_global_2x_ser_p2_periodicy(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 5*(numCellsPar-1)+5;
}

long fem_parproj_num_nodes_global_2x_ser_p1_nonperiodicy(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 2*(numCellsPar-1)+4;
}

long fem_parproj_num_nodes_global_2x_ser_p2_nonperiodicy(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 5*(numCellsPar-1)+8;
}

long fem_parproj_num_nodes_global_3x_ser_p1_periodicz(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 4*(numCellsPar-1)+4;
}

long fem_parproj_num_nodes_global_3x_ser_p2_periodicz(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 12*(numCellsPar-1)+12;
}

long fem_parproj_num_nodes_global_3x_ser_p1_nonperiodicz(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 4*(numCellsPar-1)+8;
}

long fem_parproj_num_nodes_global_3x_ser_p2_nonperiodicz(int numCellsPar) 
{ 
  // numCellsPar:  number of cells in parallel direction.

  return 12*(numCellsPar-1)+20;
}

