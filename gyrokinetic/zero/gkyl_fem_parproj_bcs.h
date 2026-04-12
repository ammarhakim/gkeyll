#pragma once

// Boundary condition types.
enum gkyl_fem_parproj_bc_type {
  GKYL_FEM_PARPROJ_PERIODIC = 0,
  GKYL_FEM_PARPROJ_DIRICHLET_GHOST, // Solution = ghost evaluated at the boundary.
  GKYL_FEM_PARPROJ_DIRICHLET_SKIN, // Solution = skin evaluated at the boundary.
  GKYL_FEM_PARPROJ_NONE,      // does not enforce a BC.
};

