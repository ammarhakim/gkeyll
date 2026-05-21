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
