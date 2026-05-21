#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

// GK-specific descriptions for each distribution moment.
// d³v = J dvpar dμ is the GK velocity-space volume element, J = B_∥*/m.
// After division by the configuration-space Jacobian J_c = sqrt(det(g_ij)), the fields live
// in physical configuration space with the SI units listed in each entry.
static const char *gkyl_gk_distribution_moments_descriptions[] = {
  "M0 = ∫ f d³v: number density [SI: m⁻³]",
  "M1 = ∫ v_∥ f d³v: parallel momentum n·u_∥ [SI: m⁻² s⁻¹]",
  "M2 = ∫ (v_∥² + 2μB/m) f d³v: total kinetic energy density divided by mass [SI: m⁻¹ s⁻²]",
  "M2par = ∫ v_∥² f d³v: parallel kinetic energy density divided by mass [SI: m⁻¹ s⁻²]",
  "M2perp = ∫ (2μB/m) f d³v: perpendicular kinetic energy density divided by mass [SI: m⁻¹ s⁻²]",
  "M2ij = ∫ vᵢvⱼ f d³v: kinetic energy tensor divided by mass [SI: m⁻¹ s⁻²]",
  "M3 = ∫ (v_∥² + 2μB/m)³/² f d³v: parallel heat flux divided by mass [SI: s⁻³]",
  "M3par = ∫ v_∥³ f d³v: parallel energy flux divided by mass [SI: s⁻³]",
  "M3perp = ∫ (2μB/m)³/² f d³v: perpendicular energy flux divided by mass [SI: s⁻³]",
  "M3ijk = ∫ vᵢvⱼvₖ f d³v: heat flux tensor in lab frame divided by mass [SI: s⁻³]",
  "MaxwellianMoments: [n, u_∥, v_ts²] best-fit Maxwellian parameters [SI: m⁻³, m s⁻¹, m² s⁻²]",
  "BiMaxwellianMoments: [n, u_∥, v_ts∥², v_ts⊥²] best-fit bi-Maxwellian parameters [SI: m⁻³, m s⁻¹, m² s⁻², m² s⁻²]",
  "LTEMoments: [n, u_∥, v_ts²] local thermodynamic equilibrium (Maxwellian) parameters [SI: m⁻³, m s⁻¹, m² s⁻²]",
  "M0M1M2: [M0, M1, M2] number density, parallel flux, total kinetic energy density divided by mass [SI: m⁻³, m⁻² s⁻¹, m⁻¹ s⁻²]",
  "M0M1M2parM2perp: [M0, M1, M2par, M2perp] number density, parallel flux, parallel and perpendicular kinetic energy densities divided by mass [SI: m⁻³, m⁻² s⁻¹, m⁻¹ s⁻², m⁻¹ s⁻²]",
  "HamiltonianMoments: [M0, M1, ∫ H f d³v] number density, parallel flux, Hamiltonian energy density divided by mass [SI: m⁻³, m⁻² s⁻¹, m⁻¹ s⁻²]",
  "M1_from_H = ∫ (∂H/∂v_∥)/m_s f d³v: parallel velocity moment from Hamiltonian [SI: m⁻² s⁻¹]",
  "EnergyMoment = ∫ H f d³v, H = ½m_s v_∥² + μB + qφ: Hamiltonian energy density divided by mass [SI: m⁻¹ s⁻²]",
  "M0EnergyM3: [M0, EnergyMoment, M3par] number density, Hamiltonian energy density divided by mass, parallel energy flux divided by mass [SI: m⁻³, m⁻¹ s⁻², s⁻³]",
  "Ni: [M0, M1ᵢ] number density and parallel flux vector [SI: m⁻³, m⁻² s⁻¹]",
  "Tij = ∫ vᵢvⱼ f d³v: stress-energy tensor divided by mass [SI: m⁻¹ s⁻²]",
};

// Same as gkyl_gk_distribution_moments_descriptions but for source terms S = ∂f/∂t|_source.
// Units gain an extra s⁻¹ factor throughout since S has dimensions of f/time.
static const char *gkyl_gk_distribution_moments_source_descriptions[] = {
  "M0 = ∫ S d³v: source density flux [SI: m⁻³ s⁻¹]",
  "M1 = ∫ v_∥ S d³v: source momentum flux n·u_∥ [SI: m⁻² s⁻²]",
  "M2 = ∫ (v_∥² + 2μB/m) S d³v: total power of the source density divided by mass [SI: m⁻¹ s⁻³]",
  "M2par = ∫ v_∥² S d³v: parallel power of the source density divided by mass [SI: m⁻¹ s⁻³]",
  "M2perp = ∫ (2μB/m) S d³v: perpendicular power of the source density divided by mass [SI: m⁻¹ s⁻³]",
  "M2ij = ∫ vᵢvⱼ S d³v: source momentum flux tensor divided by mass [SI: m⁻¹ s⁻³]",
  "M3 = ∫ (v_∥² + 2μB/m)³/² S d³v: source of isotropic heat flux divided by mass [SI: s⁻⁴]",
  "M3par = ∫ v_∥³ S d³v: source of parallel heat flux divided by mass [SI: s⁻⁴]",
  "M3perp = ∫ (2μB/m)³/² S d³v: source of perpendicular heat flux divided by mass [SI: s⁻⁴]",
  "M3ijk = ∫ vᵢvⱼvₖ S d³v: source of 3rd-order velocity tensor divided by mass [SI: s⁻⁴]",
  "MaxwellianMoments: [n_s, u_∥, v_ts²] Maxwellian fit parameters for source S: source density flux, drift velocity, thermal speed squared [SI: m⁻³ s⁻¹, m s⁻¹, m² s⁻²]",
  "BiMaxwellianMoments: [n_s, u_∥, v_ts∥², v_ts⊥²] bi-Maxwellian fit parameters for source S: source density flux, drift velocity, parallel and perpendicular thermal speeds squared [SI: m⁻³ s⁻¹, m s⁻¹, m² s⁻², m² s⁻²]",
  "LTEMoments: [n_s, u_∥, v_ts²] LTE fit parameters for source S: source density flux, drift velocity, thermal speed squared [SI: m⁻³ s⁻¹, m s⁻¹, m² s⁻²]",
  "M0M1M2: [M0, M1, M2] source density flux, source momentum flux, total power of the source density divided by mass [SI: m⁻³ s⁻¹, m⁻² s⁻², m⁻¹ s⁻³]",
  "M0M1M2parM2perp: [M0, M1, M2par, M2perp] source density flux, source momentum flux, parallel and perpendicular power of the source density divided by mass [SI: m⁻³ s⁻¹, m⁻² s⁻², m⁻¹ s⁻³, m⁻¹ s⁻³]",
  "HamiltonianMoments: [M0, M1, ∫ H S d³v] source density flux, source momentum flux, source Hamiltonian power density divided by mass [SI: m⁻³ s⁻¹, m⁻² s⁻², m⁻¹ s⁻³]",
  "M1_from_H = ∫ (∂H/∂v_∥)/m_s S d³v: source momentum flux computed from the Hamiltonian gradient [SI: m⁻² s⁻²]",
  "EnergyMoment = ∫ H S d³v, H = ½m_s v_∥² + μB + qφ: source Hamiltonian power density divided by mass [SI: m⁻¹ s⁻³]",
  "M0EnergyM3: [M0, EnergyMoment, M3par] source density flux, source Hamiltonian power density divided by mass, source of parallel heat flux divided by mass [SI: m⁻³ s⁻¹, m⁻¹ s⁻³, s⁻⁴]",
  "Ni: [M0, M1ᵢ] source density flux and source momentum flux vector [SI: m⁻³ s⁻¹, m⁻² s⁻²]",
  "Tij = ∫ vᵢvⱼ S d³v: source stress-energy tensor divided by mass [SI: m⁻¹ s⁻³]",
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
