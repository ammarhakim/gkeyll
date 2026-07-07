#pragma once

#include <string.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_sources_implicit_priv.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_gr_spacetime.h>

struct gkyl_moment_em_coupling {
  struct gkyl_rect_grid grid; // Grid over which the equations are solved.
  int ndim; // Number of grid dimensions.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.

  bool is_charged_species; // Does the fluid species carry a charge?
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  bool static_field; // Is the plasma field static? If true, only J is updated to new time step. 
  bool ramp_app_E; // Use a linear ramp function for initializing external electric fields.
  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  bool ramp_app_curr; // Use a linear ramp function for initializing applied currents.
  double t_ramp_curr; // Ramp-up time for the linear ramp function for initializing applied currents.

  bool has_collision; // Run with collisions switched on.
  bool use_rel; // Assume special relativistic fluid species.

  // Matrix of scaling factors for collision frequencies. Should be symmetric (i.e. nu_base_sr = nu_base_rs).
  // These are defined such that nu_sr = nu_base_sr / rho_s, and nu_rs = nu_base_rs / rho_r.
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool use_explicit_em_coupling; // Use the explicit source-solver for handling moment-EM coupling (not operational yet).

  bool has_nT_sources; // Run with number density and temperature sources.

  // Source-family parameter bundles, copied verbatim from the input (see
  // gkyl_moment_em_coupling.h). Each family runs when its 'enabled' flag is
  // set; the solver reads e.g. mom_em->friction.Z.
  struct gkyl_moment_friction friction; // Frictional sources.
  struct gkyl_moment_volume_sources volume_sources; // Volume-based geometrical sources.
  struct gkyl_moment_reactivity reactivity; // Reactive sources.
  struct gkyl_moment_einstein_medium einstein_medium; // Coupled fluid-Einstein sources (plane-symmetric spacetimes).
  struct gkyl_moment_gr_sources gr_ultra_rel; // GR Euler sources, ultra-relativistic equation of state.
  struct gkyl_moment_gr_sources gr_euler; // GR Euler sources, general equation of state.
  struct gkyl_moment_gr_twofluid gr_twofluid; // GR two-fluid sources.
  struct gkyl_moment_vacuum_einstein vacuum_einstein; // Vacuum Einstein sources, Bona-Masso formalism.
  struct gkyl_moment_vacuum_einstein vacuum_einstein_conformal; // Vacuum Einstein sources, conformal Bona-Masso formalism.
  struct gkyl_moment_gr_sources gr_mhd; // GR MHD sources.

  bool has_gr_em_coupling; // Run with explicit (special-)relativistic multi-fluid + Maxwell coupling (SI units, separate fluid/EM inputs).
  struct gkyl_gr_euler_eos gr_em_eos[GKYL_MAX_SPECIES]; // Per-species EOS (IDEAL or APPROXIMATE_SYNGE/RCC) for relativistic-Euler primitive recovery in the GR-EM coupling.
};
