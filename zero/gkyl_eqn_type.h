#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EQN_EULER,     // Euler equations
  GKYL_EQN_SR_EULER,  // SR Euler equations
  GKYL_EQN_ISO_EULER, // Isothermal Euler equations
  GKYL_EQN_COLDFLUID, // Cold fluid equations
  GKYL_EQN_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL, // Maxwell equations
  GKYL_EQN_MHD,  // Ideal MHD equations
  GKYL_EQN_BURGERS, // Burgers equations
  GKYL_EQN_ADVECTION, // Scalar advection equation
};

// Identifiers for specific gyrokinetic model types
enum gkyl_gkmodel_id {
  GKYL_GK_MODEL_GEN_GEO = 0, // General geometry GK. This is default
  GKYL_GK_MODEL_NO_BY = 1, // General geometry GK, but no toroidal field (by = 0)
};

// Identifiers for specific gyrokinetic field object types
enum gkyl_gkfield_id {
  GKYL_GK_FIELD_ES = 0, // Electrostatic GK. This is default
  GKYL_GK_FIELD_ADIABATIC = 1, // Adiabatic electrons GK field, phi = phi_sheath + (T_e/e)*ln(n_i/n_is)
  GKYL_GK_FIELD_ES_IWL = 2, // Inner-wall limited ES.
  GKYL_GK_FIELD_EM = 3, // Electromagnetic GK
};

// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_PHI = 1, // Poisson (only phi)
  GKYL_FIELD_PHI_A = 2, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL = 3, // no field is present
};

// Identifiers for subsidary models
// These are used to distinguish things like special relativistic from non-relativistic
enum gkyl_model_id {
  GKYL_MODEL_DEFAULT = 0, // No subsidiary model specified
  GKYL_MODEL_SR = 1,
  GKYL_MODEL_GEN_GEO = 2,
};

// Identifiers for specific collision object types
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default
  GKYL_BGK_COLLISIONS, // BGK Collision operator
  GKYL_LBO_COLLISIONS, // LBO Collision operator
  GKYL_FPO_COLLISIONS, // FPO Collision operator
};

// Identifiers for specific source object types
enum gkyl_source_id {
  GKYL_NO_SOURCE = 0, // No source. This is default
  GKYL_FUNC_SOURCE, // Function source
  GKYL_PROJ_SOURCE, // Source given by projection object determined by gkyl_projection_id
  GKYL_BFLUX_SOURCE // Source which scales to boundary fluxes
};

// Identifiers for specific projection object types
enum gkyl_projection_id {
  GKYL_PROJ_FUNC = 0, // Function projection. This is default
  GKYL_PROJ_MAXWELLIAN, // Maxwellian projection
  GKYL_PROJ_BIMAXWELLIAN, // Bi-Maxwellian projection
};

// Identifiers for specific radiation object types
enum gkyl_radiation_id {
  GKYL_NO_RADIATION = 0, // No radiation. This is default
  GKYL_GK_RADIATION, // Radiation in gyrokinetic equations.
};

// Identifiers for specific reaction object types
enum gkyl_react_id {
  GKYL_NO_REACT = 0, // No reactions. This is default
  GKYL_REACT_IZ, // Ionization.
  GKYL_REACT_CX, // Charge exchange.
  GKYL_REACT_RECOMB, // Recombination.
};

// Identifiers for different ion reaction types
enum gkyl_ion_type
{
  GKYL_ION_H = 0,  // Hydrogen ions
  GKYL_ION_HE = 1, // Helium ions
  GKYL_ION_LI = 2, // Lithium ions
  GKYL_ION_BE = 3, // Beryllium ions
  GKYL_ION_B = 4,  // Boron ions
  GKYL_ION_C = 5,  // Carbon ions
  GKYL_ION_N = 6,  // Nitrogen ions
  GKYL_ION_O = 7,  // Oxygen ions
  GKYL_ION_AR = 8,  // Argon ions
};

// Identifiers for different self in reaction
enum gkyl_react_self_type
{
  GKYL_SELF_ELC = 0, // Electron species in reaction
  GKYL_SELF_ION = 1, // Ion species in reaction 
  GKYL_SELF_DONOR = 2, // Donating species in reaction (giving up electron)
  GKYL_SELF_RECVR = 3, // Receiving species in reaction (receiving electron)
};

// Identifiers for specific geometry types
enum gkyl_geometry_id {
  GKYL_TOKAMAK, // Tokamak Geometry from Efit
  GKYL_MIRROR, // Mirror Geometry from Efit
  GKYL_MAPC2P, // General geometry from user provided mapc2p
  GKYL_GEOMETRY_FROMFILE, // Geometry from file
};

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD = 0, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD, // Gauss-Lobatto quadrature
};

/** Flags for indicating acting edge of velocity space */
enum gkyl_vel_edge { 
  GKYL_VX_LOWER, GKYL_VY_LOWER, GKYL_VZ_LOWER, 
  GKYL_VX_UPPER, GKYL_VY_UPPER, GKYL_VZ_UPPER 
};
