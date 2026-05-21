#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

// GK-specific descriptions for each distribution moment.
// d^3v = J dvpar dmu is the GK velocity-space volume element, J = B_par*/m.
// After division by the configuration-space Jacobian J_c = sqrt(det(g_ij)), the fields live
// in physical configuration space.
static const char *gkyl_gk_distribution_moments_descriptions[] = {
  "M0 = int f d^3v: number density",
  "M1 = int v_par f d^3v: parallel momentum n*u_par",
  "M2 = int (v_par^2 + 2*mu*B/m) f d^3v: total kinetic energy density divided by mass",
  "M2par = int v_par^2 f d^3v: parallel kinetic energy density divided by mass",
  "M2perp = int (2*mu*B/m) f d^3v: perpendicular kinetic energy density divided by mass",
  "M2ij = int v_i*v_j f d^3v: kinetic energy tensor divided by mass",
  "M3 = int (v_par^2 + 2*mu*B/m)^(3/2) f d^3v: parallel heat flux divided by mass",
  "M3par = int v_par^3 f d^3v: parallel energy flux divided by mass",
  "M3perp = int (2*mu*B/m)^(3/2) f d^3v: perpendicular energy flux divided by mass",
  "M3ijk = int v_i*v_j*v_k f d^3v: heat flux tensor in lab frame divided by mass",
  "MaxwellianMoments: [n, u_par, v_ts^2] best-fit Maxwellian parameters",
  "BiMaxwellianMoments: [n, u_par, v_ts_par^2, v_ts_perp^2] best-fit bi-Maxwellian parameters",
  "LTEMoments: [n, u_par, v_ts^2] local thermodynamic equilibrium (Maxwellian) parameters",
  "M0M1M2: [M0, M1, M2] number density, parallel flux, total kinetic energy density divided by mass",
  "M0M1M2parM2perp: [M0, M1, M2par, M2perp] number density, parallel flux, parallel and perpendicular kinetic energy densities divided by mass",
  "HamiltonianMoments: [M0, M1, int H f d^3v] number density, parallel flux, Hamiltonian energy density divided by mass",
  "M1_from_H = int (dH/dv_par)/m_s f d^3v: parallel velocity moment from Hamiltonian",
  "EnergyMoment = int H f d^3v, H = (1/2)*m_s*v_par^2 + mu*B + q*phi: Hamiltonian energy density divided by mass",
  "M0EnergyM3: [M0, EnergyMoment, M3par] number density, Hamiltonian energy density divided by mass, parallel energy flux divided by mass",
  "Ni: [M0, M1_i] number density and parallel flux vector",
  "Tij = int v_i*v_j f d^3v: stress-energy tensor divided by mass",
};

// Same as gkyl_gk_distribution_moments_descriptions but for source terms S = df/dt|_source.
static const char *gkyl_gk_distribution_moments_source_descriptions[] = {
  "M0 = int S d^3v: source density flux",
  "M1 = int v_par S d^3v: source momentum flux n*u_par",
  "M2 = int (v_par^2 + 2*mu*B/m) S d^3v: total power of the source density divided by mass",
  "M2par = int v_par^2 S d^3v: parallel power of the source density divided by mass",
  "M2perp = int (2*mu*B/m) S d^3v: perpendicular power of the source density divided by mass",
  "M2ij = int v_i*v_j S d^3v: source momentum flux tensor divided by mass",
  "M3 = int (v_par^2 + 2*mu*B/m)^(3/2) S d^3v: source of isotropic heat flux divided by mass",
  "M3par = int v_par^3 S d^3v: source of parallel heat flux divided by mass",
  "M3perp = int (2*mu*B/m)^(3/2) S d^3v: source of perpendicular heat flux divided by mass",
  "M3ijk = int v_i*v_j*v_k S d^3v: source of 3rd-order velocity tensor divided by mass",
  "MaxwellianMoments: [n_s, u_par, v_ts^2] Maxwellian fit parameters for source S: source density flux, drift velocity, thermal speed squared",
  "BiMaxwellianMoments: [n_s, u_par, v_ts_par^2, v_ts_perp^2] bi-Maxwellian fit parameters for source S: source density flux, drift velocity, parallel and perpendicular thermal speeds squared",
  "LTEMoments: [n_s, u_par, v_ts^2] LTE fit parameters for source S: source density flux, drift velocity, thermal speed squared",
  "M0M1M2: [M0, M1, M2] source density flux, source momentum flux, total power of the source density divided by mass",
  "M0M1M2parM2perp: [M0, M1, M2par, M2perp] source density flux, source momentum flux, parallel and perpendicular power of the source density divided by mass",
  "HamiltonianMoments: [M0, M1, int H S d^3v] source density flux, source momentum flux, source Hamiltonian power density divided by mass",
  "M1_from_H = int (dH/dv_par)/m_s S d^3v: source momentum flux computed from the Hamiltonian gradient",
  "EnergyMoment = int H S d^3v, H = (1/2)*m_s*v_par^2 + mu*B + q*phi: source Hamiltonian power density divided by mass",
  "M0EnergyM3: [M0, EnergyMoment, M3par] source density flux, source Hamiltonian power density divided by mass, source of parallel heat flux divided by mass",
  "Ni: [M0, M1_i] source density flux and source momentum flux vector",
  "Tij = int v_i*v_j S d^3v: source stress-energy tensor divided by mass",
};

/**
 * Create new gyrokinetic moment type object. Valid 'mom' strings are "M0",
 * "M1", "M2", "M2par", "M2perp", "M3par", "M3perp", "ThreeMoments",
 * "FourMoments", "HamiltonianMoments".
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param conf_range Configuration-space range.
 * @param mass Mass of species.
 * @param charge Charge of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU.
 */
struct gkyl_mom_type* 
gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  double mass, double charge, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu);

/**
 * Create new integrated gyrokinetic moment type object.
 * Valid 'mom' strings are "ThreeMoments", "FourMoments", "HamiltonianMoments".
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param conf_range Configuration-space range.
 * @param mass Mass of species.
 * @param vel_map Velocity space mapping object.
 * @param gk_geom Geometry object.
 * @param phi Electrostatic potential (for Hamiltonian moment).
 * @param mom_type Name of moment to compute.
 * @param use_gpu bool to determine if on GPU.
 */
struct gkyl_mom_type* 
gkyl_int_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type, bool use_gpu);
