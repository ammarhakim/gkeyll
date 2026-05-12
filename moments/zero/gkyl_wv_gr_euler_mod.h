#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler.h>  // share enum gkyl_wv_gr_euler_rp with the packed variant

// Auxfields handed to the modular GR Euler equation object. Following the
// vlasov auxfields pattern (see vlasov/zero/gkyl_dg_vlasov.h), this struct
// only carries pointers; the underlying array is owned by moment_spacetime.
struct gkyl_wv_gr_euler_mod_auxfields {
  const struct gkyl_array *prods; // spacetime products array
                                  // (layout: gkyl_moment_spacetime_products.h)
};

// Input context for the modular GR Euler constructor.
struct gkyl_wv_gr_euler_mod_inp {
  double gas_gamma;                       // Adiabatic index.
  struct gkyl_range conf_range;           // Configuration-space range used to
                                          // index into the products array.
  enum gkyl_wv_gr_euler_rp rp_type;       // Riemann-solver choice (Lax/Roe/HLL).
  bool use_gpu;
};

// Construct a modular GR Euler equation object using the default HLL solver.
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_mod_new(double gas_gamma, const struct gkyl_range *conf_range,
  bool use_gpu);

// Construct a modular GR Euler equation object from an input context.
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_mod_inew(const struct gkyl_wv_gr_euler_mod_inp *inp);

// Attach the spacetime-products array to a modular GR Euler equation. The
// array pointer must remain valid for the lifetime of the equation; in the
// standard app flow it is owned by struct moment_spacetime.
void
gkyl_gr_euler_mod_set_auxfields(const struct gkyl_wv_eqn *eqn,
  struct gkyl_wv_gr_euler_mod_auxfields auxin);

// Read back the adiabatic index from a modular GR Euler equation object.
double
gkyl_wv_gr_euler_mod_gas_gamma(const struct gkyl_wv_eqn *eqn);
