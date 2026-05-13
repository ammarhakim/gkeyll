#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_tetrad.h>  // share enum gkyl_wv_gr_euler_tetrad_rp with packed

// Auxfields handed to the modular tetrad-basis GR Euler equation object.
// Mirrors gkyl_wv_gr_euler_mod_auxfields — the spacetime products layout is
// identical for both mod variants (tetrad does not require stored basis
// vectors; the flat-frame factorization is purely a flux-computation
// strategy, see wv_gr_euler_tetrad_mod.c).
struct gkyl_wv_gr_euler_tetrad_mod_auxfields {
  const struct gkyl_array *prods; // spacetime products array
                                  // (layout: gkyl_moment_spacetime_products.h)
};

// Input context for the modular tetrad GR Euler constructor.
struct gkyl_wv_gr_euler_tetrad_mod_inp {
  double gas_gamma;                                // Adiabatic index.
  struct gkyl_range conf_range;                    // For indexing into prods.
  enum gkyl_wv_gr_euler_tetrad_rp rp_type;         // Lax / Roe / HLL.
  bool use_gpu;
};

// Construct a modular tetrad GR Euler equation object using the default HLL solver.
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_new(double gas_gamma,
  const struct gkyl_range *conf_range, bool use_gpu);

// Construct a modular tetrad GR Euler equation object from an input context.
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_mod_inew(
  const struct gkyl_wv_gr_euler_tetrad_mod_inp *inp);

// Attach the spacetime-products array to a modular tetrad GR Euler equation.
void
gkyl_gr_euler_tetrad_mod_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_tetrad_mod_auxfields auxin);

// Update the configuration-space range used to convert cell indices into
// offsets in the auxfields products array. Called once by the Lua flow after
// the app's range is built; tests that pass the right range at construction
// time normally do not need this setter.
void
gkyl_gr_euler_tetrad_mod_set_conf_range(const struct gkyl_wv_eqn *eqn,
  const struct gkyl_range *conf_range);

// Read back the adiabatic index from a modular tetrad GR Euler equation object.
double
gkyl_wv_gr_euler_tetrad_mod_gas_gamma(const struct gkyl_wv_eqn *eqn);
