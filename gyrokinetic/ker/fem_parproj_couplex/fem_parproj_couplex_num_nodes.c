#include <gkyl_fem_parproj_couplex_kernels.h> 
 
long fem_parproj_couplex_num_nodes_global_1x_ser_p1_periodicx(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return num_cells[0];
}

long fem_parproj_couplex_num_nodes_global_1x_ser_p2_periodicx(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 2*(num_cells[0]-1)+2;
}

long fem_parproj_couplex_num_nodes_global_1x_ser_p1_nonperiodicx(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return num_cells[0]+1;
}

long fem_parproj_couplex_num_nodes_global_1x_ser_p2_nonperiodicx(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 2*(num_cells[0]-1)+3;
}

long fem_parproj_couplex_num_nodes_global_2x_ser_p1_periodicy(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return (num_cells[0]-1)*(num_cells[1]-1)+2*(num_cells[1]-1)+num_cells[0]+1;
}

long fem_parproj_couplex_num_nodes_global_2x_ser_p2_periodicy(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 3*(num_cells[0]-1)*(num_cells[1]-1)+5*(num_cells[1]-1)+3*(num_cells[0]-1)+5;
}

long fem_parproj_couplex_num_nodes_global_2x_ser_p1_nonperiodicy(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return (num_cells[0]-1)*(num_cells[1]-1)+2*(num_cells[1]-1)+2*(num_cells[0]-1)+4;
}

long fem_parproj_couplex_num_nodes_global_2x_ser_p2_nonperiodicy(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 3*(num_cells[0]-1)*(num_cells[1]-1)+5*(num_cells[1]-1)+5*(num_cells[0]-1)+8;
}

long fem_parproj_couplex_num_nodes_global_3x_ser_p1_periodicz(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 2*(num_cells[0]-1)*(num_cells[2]-1)+4*(num_cells[2]-1)+2*(num_cells[0]-1)+4;
}

long fem_parproj_couplex_num_nodes_global_3x_ser_p2_periodicz(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 7*(num_cells[0]-1)*(num_cells[2]-1)+12*(num_cells[2]-1)+7*(num_cells[0]-1)+12;
}

long fem_parproj_couplex_num_nodes_global_3x_ser_p1_nonperiodicz(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 2*(num_cells[0]-1)*(num_cells[2]-1)+4*(num_cells[2]-1)+4*(num_cells[0]-1)+8;
}

long fem_parproj_couplex_num_nodes_global_3x_ser_p2_nonperiodicz(const int *num_cells) 
{ 
  // num_cells:  number of cells in each direction.

  return 7*(num_cells[0]-1)*(num_cells[2]-1)+12*(num_cells[2]-1)+12*(num_cells[0]-1)+20;
}

