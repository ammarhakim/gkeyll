#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_gr_euler_prim_priv.h>  // struct gkyl_gr_euler_eos
#include <gkyl_wv_gr_euler_tetrad.h>  // share enum gkyl_wv_gr_euler_tetrad_rp with packed

// Auxfields handed to the modular tetrad-basis GR Euler equation object.
// Mirrors gkyl_wv_gr_euler_mod_auxfields — the spacetime products layout is
// identical for both mod variants (tetrad does not require stored basis
// vectors; the flat-frame factorization is purely a flux-computation
// strategy, see wv_gr_euler_tetrad_mod.c).
//
// Per-callsite instrumentation buckets:
//   prim_status_wave_prop:   updated by all wave-prop primitive recovery
//                            callsites (prim_vars hook + 4 Riemann solvers
//                            + max_speed_func).
//   repair_status_wave_prop: updated by the repair_state hook when
//                            cur_repair_ctx = 1 (wave_prop side).
//   repair_status_source:    updated by the repair_state hook when
//                            cur_repair_ctx = 0 (source-step side, called
//                            from gkyl_moment_spacetime_coupling_explicit_advance).
// The source-step primitive recovery has its OWN prim_status passed
// directly through explicit_advance — that one doesn't live here.
// All status pointers may be NULL (unit-test callsites etc.).
struct gkyl_wv_gr_euler_tetrad_mod_auxfields {
  const struct gkyl_array *prods;             // spacetime products array
                                              // (layout: gkyl_moment_spacetime_products.h)
  struct gkyl_gr_euler_prim_status   *prim_status_wave_prop;
  struct gkyl_gr_euler_repair_status *repair_status_wave_prop;
  struct gkyl_gr_euler_repair_status *repair_status_source;
};

// Input context for the modular tetrad GR Euler constructor.
//
// The struct embeds an EOS bundle (struct gkyl_gr_euler_eos) so callers
// can select IDEAL (with a gas_gamma) or APPROXIMATE_SYNGE (with the
// use_rcc closure selector). The legacy `gas_gamma` field is retained
// for source-compatibility — when it is set and eos.type is the default
// IDEAL, the constructor copies it into eos.gas_gamma. New call sites
// should populate the eos field directly.
//
// Roe Riemann solvers are only available with IDEAL gas: the eigenstructure
// in sr_roe_minkowski uses the Eulderink-Mellema ideal-gas Jacobian. The
// constructor asserts(rp_type != ROE || eos.type == IDEAL) on inew.
struct gkyl_wv_gr_euler_tetrad_mod_inp {
  double gas_gamma;                                // Adiabatic index (legacy IDEAL shortcut).
  struct gkyl_gr_euler_eos eos;                    // Full EOS specification.
  struct gkyl_range conf_range;                    // For indexing into prods.
  enum gkyl_wv_gr_euler_tetrad_rp rp_type;         // Lax / Roe / HLL / HLLC.
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

// Read back the full EOS bundle from a modular tetrad GR Euler equation object.
// For IDEAL the gas_gamma field is meaningful; for APPROXIMATE_SYNGE the
// use_rcc field selects the enthalpy closure (RCC vs TM).
struct gkyl_gr_euler_eos
gkyl_wv_gr_euler_tetrad_mod_eos(const struct gkyl_wv_eqn *eqn);
