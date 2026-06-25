#pragma once

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include "tov_solver.h"

// Static TOV star spacetime, in Cartesian coordinates. Static because it uses Schwarzschild like areal coordinates - good for the static problem, but will need moving to Cartesian Kerr Schild for the collapse
// Metric read from the TOV table at the areal radius on a Cartesian grid. Covers interior and exterior, because of of gkyl_tov_eval_bl works.
// Schwarzschild-like (areal-radial) coordinates, metric:
// alpha = exp(Phi(r)), beta^i = 0, gamma_ij = delta_ij + (A(r) - 1) * x^i x^j / r^2 , A(r) = 1 / (1 - 2 m(r)/r) (regular at the origin)

// For now, this is the fixed background (use with GKYL_STATIC_GAUGE - metric is not evolved). 
// For now, not (horizon-penetratin) C Kerr-Schild: horizon-penetrating ((79),(154) from arXiv:2410.02549) - that will get upgraded for the collapse.

// Spacetime object.
struct gr_tov {
  struct gkyl_gr_spacetime spacetime;

  const struct gkyl_tov *tov; // From frozen TOV table (areal-radius lookup of Phi, m)

  double pos_x; // Position of the stellar center (x-direction)
  double pos_y; // Position of the stellar center (y-direction)
  double pos_z; // Position of the stellar center (z-direction)

  bool use_kerr_schild; // false: static Schwarzschild-like (areal) gauge (beta^i = 0, K_ij = 0)
                        // true: ingoing Cartesian Kerr-Schild gauge (nonzero beta^i, K_ij != 0), the horizon-penetrating, evolution-stable slicing.
};

struct gkyl_gr_tov_inp {
  bool use_gpu; // Whether the spacetime object is on the host (false) or the device (true).

  const struct gkyl_tov *tov; // Borrowed frozen TOV table (not owned)

  double pos_x; // Position of the stellar center (x-direction)
  double pos_y; // Position of the stellar center (y-direction)
  double pos_z; // Position of the stellar center (z-direction)

  bool use_kerr_schild; // Use the Cartesian Kerr-Schild gauge (see struct gr_tov).
};

/**
* Free TOV star spacetime object.
*
* @param ref Reference counter for the TOV spacetime.
*/
void
gkyl_gr_tov_spacetime_free(const struct gkyl_ref_count* ref);

/**
* Create a new static TOV neutron-star spacetime object.
*
* @param use_gpu Whether the spacetime object is on the host (false) or the device (true)
* @param tov From frozen TOV table (provides Phi(r), m(r) by areal radius)
* @param pos_x Stellar center (x-direction)
* @param pos_y Stellar center (y-direction)
* @param pos_z Stellar center (z-direction)
* @param use_kerr_schild Use the Cartesian Kerr-Schild gauge (false: static areal gauge).
* @return Pointer to the TOV spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_new(bool use_gpu, const struct gkyl_tov *tov, double pos_x, double pos_y, double pos_z, bool use_kerr_schild);

/**
* Create a new static TOV star spacetime object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the TOV spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_inew(const struct gkyl_gr_tov_inp* inp);
