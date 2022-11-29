#pragma once

#include <stdbool.h>

// Update status
struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_REFLECT, // perfect reflector
  GKYL_SPECIES_ABSORB, // Absorbing BCs
  GKYL_SPECIES_NO_SLIP, // no-slip boundary conditions
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_SPECIES_FUNC, // Function boundary conditions
};

// Boundary conditions on fluids
enum gkyl_fluid_species_bc_type {
  GKYL_FLUID_SPECIES_COPY = 0, // copy BCs
  GKYL_FLUID_SPECIES_ABSORB, // Absorbing BCs
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_FIELD_FUNC, // Function boundary conditions
};
