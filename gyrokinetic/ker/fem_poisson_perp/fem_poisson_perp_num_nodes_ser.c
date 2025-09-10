#include <gkyl_fem_poisson_perp_kernels.h> 
 
long fem_poisson_perp_num_nodes_global_2x_ser_p1_periodicx(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)+2;
}

long fem_poisson_perp_num_nodes_global_2x_ser_p1_nonperiodicx(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)+4;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+2*(numCells[0]-1)+2;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+4*(numCells[0]-1)+4;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+4*(numCells[1]-1)+2*(numCells[0]-1)+4;
}

long fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells: number of cells in each direction.

  return 2*(numCells[0]-1)*(numCells[1]-1)+4*(numCells[1]-1)+4*(numCells[0]-1)+8;
}

