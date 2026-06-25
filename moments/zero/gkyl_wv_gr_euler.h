#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_euler_rp {
  WV_GR_EULER_RP_HLL = 0, // Default (HLL fluxes).
  WV_GR_EULER_RP_ROE,
  WV_GR_EULER_RP_LAX,
  WV_GR_EULER_RP_HLL_LAX, // HLL in the bulk, Lax near the atmosphere/surface (low density switch).
};

struct gkyl_tov; // Frozen TOV discrete equilibrium (from the q[71] and q[72] slots), for optional well-balancing.

// Input context, packaged as a struct.
struct gkyl_wv_gr_euler_inp {
  double gas_gamma; // Adiabatic index.
  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.

  enum gkyl_wv_gr_euler_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).

  // Optional static-TOV well-balancing. Leave tov_eq = NULL (the default) for all other problems (WB inactive)
  const struct gkyl_tov *tov_eq; // Use the TOV slots for the well-balancing with discrete frozen equilibrium (NULL -> no WB and every other problem is untouched).

  double rho_atm; // Atmosphere rest-mass density (recovery rho floor).
  double p_atm; // Atmosphere pressure (recovery p floor).

  bool wb_family; // Use the Kappeli-Mishra equilibrium-family WB (false -> legacy frozen-discrete WB in slots 71,72).
  double equil_C; // Global equilibrium invariant C = h*sqrt(-g_tt) (relativistic Bernoulli constant).
  double equil_K_poly; // Cold polytrope constant K (p = K rho^gamma) defining the equilibrium EOS.
};

/**
* Create a new general relativistic Euler equations object with ideal gas equation of state.
*
* @param gas_gamma Adiabatic index.
* @param spacetime_gauge Spacetime gauge choice.
* @param reinit_freq Spacetime reinitialization frequency.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Euler equations object with ideal gas equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_new(double gas_gamma, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic Euler equations object with ideal gas equation of state, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Euler equations object with ideal gas equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_inew(const struct gkyl_wv_gr_euler_inp* inp);

/**
* Get adiabatic index.
*
* @param eqn General relativistic Euler equations object with ideal gas equation of state.
* @return Adiabatic index.
*/
double
gkyl_wv_gr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime gauge choice.
*
* @param eqn General relativistic Euler equations object with ideal gas equation of state.
* @return Spacetime gauge choice.
*/
enum gkyl_spacetime_gauge
gkyl_wv_gr_euler_spacetime_gauge(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime reinitialization frequency.
*
* @param eqn General relativistic Euler equations object with ideal gas equation of state.
* @return Spacetime reinitialization frequency.
*/
int
gkyl_wv_gr_euler_reinit_freq(const struct gkyl_wv_eqn* eqn);

/**
* Get base spacetime object.
*
* @param eqn General relativistic Euler equations object with ideal gas equation of state.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_spacetime(const struct gkyl_wv_eqn* eqn);