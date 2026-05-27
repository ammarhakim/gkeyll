#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EQN_EULER = 1,     // Euler equations
  GKYL_EQN_SR_EULER = 2,  // SR Euler equations
  GKYL_EQN_ISO_EULER = 3, // Isothermal Euler equations
  GKYL_EQN_COLDFLUID = 4, // Cold fluid equations
  GKYL_EQN_COLDFLUID_SR = 5, // Relativistic Cold fluid equations
  GKYL_EQN_TEN_MOMENT = 6, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL = 7, // Maxwell equations
  GKYL_EQN_MHD = 8,  // Ideal MHD equations
  GKYL_EQN_BURGERS = 9, // Burgers equations
  GKYL_EQN_ADVECTION = 10, // Scalar advection equation
  GKYL_EQN_GR_EULER = 11, // General relativistic Euler equations with ideal gas equation of state.
  GKYL_EQN_GR_EULER_TETRAD = 12, // General relativistic Euler equations in the tetrad basis with ideal gas equation of state.
  GKYL_EQN_GR_ULTRA_REL_EULER = 13, // General relativistic Euler equations with ultra-relativistic equation of state.
  GKYL_EQN_GR_ULTRA_REL_EULER_TETRAD = 14, // General relativistic Euler equations in the tetrad basis with ultra-relativistic equation of state.
  GKYL_EQN_GR_MAXWELL = 15, // General relativistic Maxwell equations.
  GKYL_EQN_GR_MAXWELL_TETRAD = 16, // General relativistic Maxwell equations in the tetrad basis.
  GKYL_EQN_GR_MEDIUM = 17, // Coupled fluid-Einstein equations in plane-symmetric spacetimes.
  GKYL_EQN_GR_TWOFLUID = 18, // General relativistic two-fluid equations.
  GKYL_EQN_VACUUM_EINSTEIN = 19, // Vacuum Einstein equations using the Bona-Masso formalism.
  GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL = 20, // Vacuum Einstein equations using the conformal Bona-Masso formalism.
  GKYL_EQN_GR_TWOFLUID_TETRAD = 21, // General relativistic two-fluid equations in the tetrad basis.
  GKYL_EQN_GR_MHD = 22, // General relativistic magnetohydrodynamics equations.
  GKYL_EQN_GR_MHD_TETRAD = 23, // General relativistic magnetohydrodynamics equations in the tetrad basis.
  GKYL_EQN_REACTIVE_EULER = 24, // Reactive Euler equations.
  GKYL_EQN_EULER_MIXTURE = 25, // Euler mixture equations.
  GKYL_EQN_ISO_EULER_MIXTURE = 26, // Isothermal Euler mixture equations.
  GKYL_EQN_EULER_RGFM = 27, // Euler Riemann ghost fluid equations.
  GKYL_EQN_CAN_PB_INCOMPRESS_EULER = 28, // Canonical Poisson Bracket form of incompressible Euler.
  GKYL_EQN_CAN_PB_HASEGAWA_MIMA = 29, // Canonical Poisson Bracket form of Hasegawa-Mima.
  GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI = 30, // Canonical Poisson Bracket form of Hasegawa-Wakatani.
};

// Identifiers for specific gyrokinetic field object types.
enum gkyl_gkfield_id {
  GKYL_GK_FIELD_ES = 0, // Electrostatic GK. This is default.
  GKYL_GK_FIELD_BOLTZMANN = 1, // GK Boltzmann, isothermal electrons, phi = phi_sheath + (T_e/e)*ln(n_i/n_is).
  GKYL_GK_FIELD_ADIABATIC = 2, // GK field with an adiabatic species.
  GKYL_GK_FIELD_ES_IWL = 3, // Inner-wall limited ES.
  GKYL_GK_FIELD_EM = 4, // Electromagnetic GK.
};

// Identifiers for specific field object types.
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default.
  GKYL_FIELD_PHI = 1, // Poisson (only phi).
  GKYL_FIELD_PHI_EXT_POTENTIALS = 2, // Poisson + external potentials (phi_ext, A_ext).
  GKYL_FIELD_PHI_EXT_FIELDS = 3, // Poisson + external fields (E_ext, B_ext).
  GKYL_FIELD_NULL = 4, // no field is present.
};

// Identifiers for subsidary models.
// These are used to distinguish things like special relativistic from non-relativistic.
enum gkyl_model_id {
  GKYL_MODEL_DEFAULT = 0, // No subsidiary model specified.
  GKYL_MODEL_SR = 1,
  GKYL_MODEL_CANONICAL_PB = 2,
  GKYL_MODEL_CANONICAL_PB_GR = 3,
};

// Identifiers for specific collision object types.
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default.
  GKYL_BGK_COLLISIONS = 1, // BGK Collision operator.
  GKYL_LBO_COLLISIONS = 2, // LBO Collision operator.
  GKYL_FPO_COLLISIONS = 3, // FPO Collision operator.
};

// Identifiers for specific source object types.
enum gkyl_source_id {
  GKYL_NO_SOURCE = 0, // No source. This is default.
  GKYL_FUNC_SOURCE = 1, // Function source.
  GKYL_PROJ_SOURCE = 2, // Source given by projection object determined by gkyl_projection_id.
  GKYL_BFLUX_SOURCE = 3 // Source which scales to boundary fluxes.
};

// Identifiers for specific bgk source object types.
enum gkyl_source_bgk_id {
  GKYL_NO_SOURCE_BGK = 0, // No source. This is default.
  GKYL_SOURCE_BGK_STATIC = 1, // BGK source with user specified rate and eq. function.
  GKYL_SOURCE_BGK_HEATING = 2, // Only energy source.
  GKYL_SOURCE_BGK_EXTERNAL = 3, // External particle, momentum, and energy source.
};

// Identifiers for specific projection object types
enum gkyl_projection_id {
  GKYL_PROJ_FUNC = 0, // Function projection. This is default
  GKYL_PROJ_MAXWELLIAN_PRIM = 1, // Maxwellian projection from primitive moments (n, u, T)
  GKYL_PROJ_MAXWELLIAN_LAB = 2, // Maxwellian projection from lab moments (M0, M1, M2)
  GKYL_PROJ_BIMAXWELLIAN = 3, // Bi-Maxwellian projection
  GKYL_PROJ_VLASOV_LTE = 4, // LTE (Local thermodynamic equilibrium) projection for Vlasov
                        // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
  GKYL_PROJ_MAXWELLIAN_GAUSSIAN = 5, // Maxwellian in velocity space + Gaussian in configuration space
};

// Identifiers for specific radiation object types
enum gkyl_radiation_id {
  GKYL_NO_RADIATION = 0, // No radiation. This is default
  GKYL_GK_RADIATION = 1, // Radiation in gyrokinetic equations.
  GKYL_VM_COMPTON_RADIATION = 2, // Vlasov simple Compton radiation model. 
};

// Identifiers for specific reaction object types
enum gkyl_react_id {
  GKYL_NO_REACT = 0, // No reactions. This is default
  GKYL_REACT_IZ = 1, // Ionization.
  GKYL_REACT_CX = 2, // Charge exchange.
  GKYL_REACT_RECOMB = 3, // Recombination.
};

enum gkyl_te_min_model {
  GKYL_VARY_TE_CONSERVATIVE = 0,  // Minimum temperature depends on V0, turns off at (relatively) high Te, so low chance of negative emissivity. This is default
  GKYL_VARY_TE_AGGRESSIVE = 1,  // Minimum temperature depends on V0, turns off at (relatively) low Te, so higher chance of negative emissivity
  GKYL_CONST_TE = 2,  // A constant minimum temperature, below which radiation is turned off
};

// Identifiers for different ion reaction types
enum gkyl_ion_type {
  GKYL_ION_H = 0,  // Hydrogen ions
  GKYL_ION_D = 1,  // Deuterium ions (for CX)
  GKYL_ION_HE = 2, // Helium ions
  GKYL_ION_LI = 3, // Lithium ions
  GKYL_ION_BE = 4, // Beryllium ions
  GKYL_ION_B = 5,  // Boron ions
  GKYL_ION_C = 6,  // Carbon ions
  GKYL_ION_N = 7,  // Nitrogen ions
  GKYL_ION_O = 8,  // Oxygen ions
  GKYL_ION_NE = 9, // Neon ions
  GKYL_ION_AR = 10,  // Argon ions
};

// Identifiers for different self in reaction
//  - For IZ: GKYL_SELF_ELC, GKYL_SELF_ION, GKYL_SELF_DONOR.
//  - For CX: GKYL_SELF_ION, GKYL_SELF_PARTNER.
//  - For RECOMB: GKYL_SELF_ELC, GKYL_SELF_ION, GKYL_SELF_RECVR.
enum gkyl_react_self_type {
  GKYL_SELF_ELC = 0, // Electron species in reaction
  GKYL_SELF_ION = 1, // Ion species in reaction 
  GKYL_SELF_DONOR = 2, // Donating species in reaction (giving up electron)
  GKYL_SELF_RECVR = 3, // Receiving species in reaction (receiving electron)
  GKYL_SELF_PARTNER = 4, // Neutral species in CX
};

// Identifiers for specific geometry types
enum gkyl_geometry_id {
  GKYL_GEOMETRY_NONE = 0, // No geometry, use Cartesian.
  GKYL_GEOMETRY_TOKAMAK = 1, // Tokamak Geometry from Efit.
  GKYL_GEOMETRY_MIRROR = 2, // Mirror Geometry from Efit.
  GKYL_GEOMETRY_MAPC2P = 3, // General geometry from user provided mapc2p.
  GKYL_GEOMETRY_FROMFILE = 4, // Geometry from file.
};

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD = 0, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD = 1, // Gauss-Lobatto quadrature
  GKYL_POSITIVITY_QUAD = 2 // Positivity quadrature nodes
};

/** Flags for indicating acting edge of velocity space */
enum gkyl_vel_edge { 
  GKYL_VX_LOWER = 0, 
  GKYL_VY_LOWER = 1, 
  GKYL_VZ_LOWER = 2, 
  GKYL_VX_UPPER = 3, 
  GKYL_VY_UPPER = 4, 
  GKYL_VZ_UPPER = 5
};

// Identifiers for FLR models (in gyrokinetics).
enum gkyl_gk_flr_type {
  GKYL_GK_FLR_NONE = 0, // No FLR effects.
  GKYL_GK_FLR_PADE_CONST = 1, // Pade-based approx. w/ const. rho_ts=sqrt(Tperp_s/m_s)
};

// Gyrokinetic anomaous diffusion models.
enum gkyl_gk_anomalous_diff_id {
  GKYL_GK_ANOMALOUS_DIFF_NONE = 0, // No anomalous diffusion.
  GKYL_GK_ANOMALOUS_DIFF_D = 1, // Specify D only.
  GKYL_GK_ANOMALOUS_DIFF_D_CHI = 2, // Specify D and chi.
};

enum gkyl_gk_collisionless_type {
  GKYL_GK_COLLISIONLESS_NONE = 0, // Turns off collisionless terms.
  GKYL_GK_COLLISIONLESS_ES = 1, // Electrostatic.
  GKYL_GK_COLLISIONLESS_ES_NO_BY = 2, // Electrostatic, w/o toroidal field (b_y=0).
  GKYL_GK_COLLISIONLESS_EM_BPERP = 3, // Electromagnetic with B_perp fluctuations.
  GKYL_GK_COLLISIONLESS_EM_BPAR = 4, // Electromagnetic with B_par fluctuations.
  GKYL_GK_COLLISIONLESS_EM = 5, // Electromagnetic with B_perp and B_par fluctuations.
  GKYL_GK_COLLISIONLESS_NEUTRAL = 6, // Neutral collisionless terms.
};
