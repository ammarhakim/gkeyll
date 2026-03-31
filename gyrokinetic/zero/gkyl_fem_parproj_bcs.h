#pragma once

// Boundary condition types.
enum gkyl_fem_parproj_bc_type {
  GKYL_FEM_PARPROJ_PERIODIC = 0,
  GKYL_FEM_PARPROJ_DIRICHLET, // sets the value.
  GKYL_FEM_PARPROJ_NONE,      // does not enforce a BC.
};

