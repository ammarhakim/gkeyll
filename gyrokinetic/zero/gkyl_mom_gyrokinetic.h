#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

// GK-specific descriptions for each distribution moment.
static const char *gkyl_gk_distribution_moments_descriptions[] = {
  "Number density of the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel momentum density of the gyrokinetic species (number density times parallel drift velocity). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Total kinetic energy density of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel kinetic energy density of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Perpendicular kinetic energy density of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Second-order velocity-space moment tensor of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Isotropic heat flux density of the gyrokinetic species divided by mass (third-order velocity moment). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel energy flux density of the gyrokinetic species divided by mass (third-order parallel velocity moment). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Perpendicular energy flux density of the gyrokinetic species divided by mass (third-order perpendicular velocity moment). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Third-order heat flux tensor of the gyrokinetic species in the lab frame divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Maxwellian parameters for the gyrokinetic species: number density, parallel drift velocity, and thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Bi-Maxwellian parameters for the gyrokinetic species: number density, parallel drift velocity, parallel thermal speed squared, and perpendicular thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Local thermodynamic equilibrium (Maxwellian) parameters for the gyrokinetic species: number density, parallel drift velocity, and thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Number density, parallel momentum density, and total kinetic energy density divided by mass for the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Number density, parallel momentum density, and parallel and perpendicular kinetic energy densities divided by mass for the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Number density, parallel momentum density, and Hamiltonian energy density divided by mass for the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel momentum density of the gyrokinetic species computed from the Hamiltonian gradient. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Hamiltonian energy density of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Number density, Hamiltonian energy density divided by mass, and parallel energy flux divided by mass for the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Number density and parallel momentum flux vector for the gyrokinetic species. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Stress-energy tensor of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
};

// GK-specific descriptions for each distribution moment of the source term.
static const char *gkyl_gk_distribution_moments_source_descriptions[] = {
  "Source density flux of the gyrokinetic species source term (number density rate of change due to the source). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source momentum flux of the gyrokinetic species source term (parallel momentum density rate of change due to the source). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Total power of the gyrokinetic species source term divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel power of the gyrokinetic species source term divided by mass (parallel kinetic energy input rate). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Perpendicular power of the gyrokinetic species source term divided by mass (perpendicular kinetic energy input rate). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Second-order velocity-space moment tensor of the gyrokinetic species source term divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Isotropic heat flux source of the gyrokinetic species divided by mass (third-order velocity moment of the source term). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel heat flux source of the gyrokinetic species divided by mass (third-order parallel velocity moment of the source term). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Perpendicular heat flux source of the gyrokinetic species divided by mass (third-order perpendicular velocity moment of the source term). For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Third-order velocity tensor source of the gyrokinetic species divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Maxwellian parameters for the gyrokinetic species source term: source density flux, drift velocity, and thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Bi-Maxwellian parameters for the gyrokinetic species source term: source density flux, drift velocity, parallel thermal speed squared, and perpendicular thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Local thermodynamic equilibrium (Maxwellian) parameters for the gyrokinetic species source term: source density flux, drift velocity, and thermal speed squared. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source density flux, source momentum flux, and total power divided by mass for the gyrokinetic species source term. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source density flux, source momentum flux, and parallel and perpendicular powers divided by mass for the gyrokinetic species source term. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source density flux, source momentum flux, and Hamiltonian power density divided by mass for the gyrokinetic species source term. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Parallel momentum flux of the gyrokinetic species source term computed from the Hamiltonian gradient. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Hamiltonian power density of the gyrokinetic species source term divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source density flux, Hamiltonian power density divided by mass, and parallel heat flux source divided by mass for the gyrokinetic species source term. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Source density flux and source momentum flux vector for the gyrokinetic species source term. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
  "Stress-energy tensor of the gyrokinetic species source term divided by mass. For additional detail, documentation can be found at https://gkeyll.readthedocs.io/en/latest/",
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
