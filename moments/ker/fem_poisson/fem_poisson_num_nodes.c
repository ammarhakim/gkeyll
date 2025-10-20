#include <gkyl_fem_poisson_kernels.h> 
 
long fem_poisson_num_nodes_global_1x_ser_p1_periodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[0];
}

long fem_poisson_num_nodes_global_1x_ser_p2_periodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 2*(numCells[0]-1)+2;
}

long fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[0]+1;
}

long fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 2*(numCells[0]-1)+3;
}

long fem_poisson_num_nodes_global_2x_ser_p1_periodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[1]+(numCells[0]-1)*(numCells[1]-1)+numCells[0]-1;
}

long fem_poisson_num_nodes_global_2x_ser_p2_periodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 3*(numCells[0]-1)*(numCells[1]-1)+3*(numCells[1]-1)+3*(numCells[0]-1)+3;
}

long fem_poisson_num_nodes_global_2x_ser_p1_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return numCells[1]+(numCells[0]-1)*(numCells[1]-1)+2*(numCells[0]-1)+1;
}

long fem_poisson_num_nodes_global_2x_ser_p2_periodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 3*(numCells[0]-1)*(numCells[1]-1)+3*(numCells[1]-1)+5*(numCells[0]-1)+5;
}

long fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return (numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+numCells[0]+1;
}

long fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_periodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 3*(numCells[0]-1)*(numCells[1]-1)+5*(numCells[1]-1)+3*(numCells[0]-1)+5;
}

long fem_poisson_num_nodes_global_2x_ser_p1_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return (numCells[0]-1)*(numCells[1]-1)+2*(numCells[1]-1)+2*(numCells[0]-1)+4;
}

long fem_poisson_num_nodes_global_2x_ser_p2_nonperiodicx_nonperiodicy(const int *numCells) 
{ 
  // numCells:       number of cells in each direction.

  return 3*(numCells[0]-1)*(numCells[1]-1)+5*(numCells[1]-1)+5*(numCells[0]-1)+8;
}

