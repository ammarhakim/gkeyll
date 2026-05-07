#pragma once

#include <gkyl_positivity_fdot_restrict.h>

struct gkyl_positivity_fdot_restrict {
  struct gkyl_basis basis;
  enum gkyl_positivity_fdot_restrict_type type;
  double safety_factor;
  bool use_gpu;
  struct gkyl_positivity_fdot_restrict *on_dev;

  void (*fdot_restriction_func)(struct gkyl_positivity_fdot_restrict *up,
    const double *f, double *dfdt, double dt);

  void (*advance_cu_func)(struct gkyl_positivity_fdot_restrict *up,
    const struct gkyl_range *range, const struct gkyl_array *f,
    struct gkyl_array *dfdt, double dt);
};

#ifdef GKYL_HAVE_CUDA

/**
 * Set CUDA kernel launcher function pointers based on restriction type.
 * Called at init time so advance_cu dispatches without runtime branching.
 */
void gkyl_positivity_fdot_restrict_cu_set_ptrs(struct gkyl_positivity_fdot_restrict *up);
#endif

/**
 * Restrict df/dt at cell average to maintain positivity.
 *
 * For cell-average mode, we restrict only the modal coefficient 0 (cell average)
 * to ensure the forward-Euler prediction stays non-negative.
 */
GKYL_CU_DH
static void
restrict_fdot_avg(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  double f_pred = f[0] + dt * dfdt[0]; // Euler step
  if (f_pred < 0.0) {
    dfdt[0] = (up->safety_factor - 1.0) * f[0] / dt; // Push towards zero with safety factor
  }
}

/**
 * Restrict df/dt at quadrature nodes to maintain positivity.
 *
 * For quadrature mode, we restrict dfdt at each quadrature node individually
 * to ensure predictions stay non-negative at those points, then project back
 * to modal space using L2-projection only if restrictions were needed.
 */
GKYL_CU_DH
static void
restrict_fdot_quad(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  const int nquad = up->basis.num_quad;
  const double safety_factor = up->safety_factor;
  bool reproject = false;

  double fquad[nquad];
  double dfdt_quad[nquad];

  // Evaluate dfdt at all quadrature nodes
  for (int k = 0; k < nquad; ++k) {
    up->basis.modal_to_quad_nodal(dfdt, dfdt_quad, k);
    up->basis.modal_to_quad_nodal(f, fquad, k);

    double f_pred_q = fquad[k] + dt * dfdt_quad[k]; // Euler step

    if (f_pred_q < 0.0) {
      dfdt_quad[k] = (safety_factor - 1.0) * fquad[k] / dt;
      reproject = true;
    }
  }

  // Project restricted dfdt only if needed
  if (reproject) {
    for (int k = 0; k < nquad; ++k) {
      up->basis.quad_nodal_to_modal(dfdt_quad, dfdt, k);
    }
  }
}

/**
 * Diode mode at cell average: set df/dt=0 if both f<0 and df/dt<0.
 *
 * This acts like a diode that blocks outflow when the distribution is already negative.
 */
GKYL_CU_DH
static void
restrict_fdot_diode_avg(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  // If both cell average f and dfdt are negative, set cell average dfdt to 0.
  if (f[0] < 0.0 && dfdt[0] < 0.0) {
    dfdt[0] = 0.0;
  }
}

/**
 * Diode mode at quadrature nodes: set df/dt=0 if both f<0 and df/dt<0.
 *
 * This acts like a diode at each quadrature point that blocks outflow
 * when the distribution is already negative.
 */
GKYL_CU_DH
static void
restrict_fdot_diode_quad(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  const int nquad = up->basis.num_quad;
  bool reproject = false;

  double dfdt_quad[nquad];
  double f_quad[nquad];

  for (int k = 0; k < nquad; ++k) {
    up->basis.modal_to_quad_nodal(f, f_quad, k);
    up->basis.modal_to_quad_nodal(dfdt, dfdt_quad, k);

    // If both f and dfdt are negative, set dfdt=0 (block outflow).
    if (f_quad[k] < 0.0 && dfdt_quad[k] < 0.0) {
      dfdt_quad[k] = 0.0;
      reproject = true;
    }
  }

  // Project modified dfdt back to modal space only if needed
  if (reproject) {
    for (int k = 0; k < nquad; ++k) {
      up->basis.quad_nodal_to_modal(dfdt_quad, dfdt, k);
    }
  }
}