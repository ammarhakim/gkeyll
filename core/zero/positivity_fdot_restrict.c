#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>

#include <gkyl_alloc.h>
#include <gkyl_range.h>

#include <float.h>
#include <math.h>

struct gkyl_positivity_fdot_restrict*
gkyl_positivity_fdot_restrict_new(struct gkyl_positivity_fdot_restrict_inp inp)
{
  struct gkyl_positivity_fdot_restrict *up = gkyl_malloc(sizeof(*up));
  
  up->basis = inp.basis;
  up->mode = inp.mode;
  up->safety_factor = inp.safety_factor;
  up->num_basis = inp.basis.num_basis;
  
  // Allocate temporary array for quadrature values.
  up->fquad = gkyl_malloc(sizeof(double[up->num_basis]));
  
  return up;
}

/**
 * Restrict df/dt at cell average to maintain positivity.
 * 
 * For cell-average mode, we restrict only the modal coefficient 0 (cell average)
 * to ensure the forward-Euler prediction stays non-negative.
 */
static void
restrict_fdot_avg(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  // Cell average is stored in modal coefficient 0.
  double f_avg = f[0];
  double dfdt_avg = dfdt[0];
  
  // Predict: f_pred = f_avg + dt * dfdt_avg
  double f_pred = f_avg + dt * dfdt_avg;
  
  // If prediction would be negative, clamp dfdt_avg to reach safety floor.
  if (f_pred < 0.0 && f_avg > 0.0) {
    // Target floor: safety_factor * f_avg
    // We want: f_avg + dt * dfdt_new = safety_factor * f_avg
    // So: dfdt_new = (safety_factor * f_avg - f_avg) / dt = (safety_factor - 1) * f_avg / dt
    double dfdt_new = (up->safety_factor - 1.0) * f_avg / dt;
    dfdt[0] = dfdt_new;
  }
}

/**
 * Restrict df/dt at quadrature nodes to maintain positivity.
 * 
 * For quadrature mode, we restrict dfdt at each quadrature node individually
 * to ensure predictions stay non-negative at those points, then project back
 * to modal space using L2-projection only if restrictions were needed.
 */
static void
restrict_fdot_quad(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt)
{
  const int nbasis = up->num_basis;
  const double safety_factor = up->safety_factor;
  bool reproject = false;
  
  // Temporary array to store dfdt values at quadrature nodes
  double dfdt_quad[nbasis];
  
  // Evaluate dfdt at all quadrature nodes
  for (int k = 0; k < nbasis; ++k) {
    up->basis.modal_to_quad_nodal(dfdt, dfdt_quad, k);
  }
  
  // Check and restrict dfdt at each quadrature node
  for (int k = 0; k < nbasis; ++k) {
    // Evaluate f at quadrature node k
    up->basis.modal_to_quad_nodal(f, up->fquad, k);
    double f_q = up->fquad[k];
    double dfdt_q = dfdt_quad[k];
    
    // Predict at this node: f_pred = f + dt * dfdt
    double f_pred_q = f_q + dt * dfdt_q;
    
    // If prediction would be negative, restrict dfdt at this node
    if (f_pred_q < 0.0 && f_q > 0.0) {
      double dfdt_new_q = (safety_factor - 1.0) * f_q / dt;
      dfdt_quad[k] = dfdt_new_q;
      reproject = true;
    }
  }
  
  // Project restricted dfdt back to modal space using L2-projection only if needed
  if (reproject) {
    for (int k = 0; k < nbasis; ++k) {
      up->basis.quad_nodal_to_modal(dfdt_quad, dfdt, k);
    }
  }
}

void
gkyl_positivity_fdot_restrict_advance(gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  
  if (up->mode == GKYL_POSITIVITY_FDOT_RESTRICT_QUAD) {
    // Quadrature mode: restrict at each quadrature node.
    while (gkyl_range_iter_next(&iter)) {
      long idx = gkyl_range_idx(range, iter.idx);
      const double *f_ptr = gkyl_array_cfetch(f, idx);
      double *dfdt_ptr = gkyl_array_fetch(dfdt, idx);
      
      restrict_fdot_quad(up, f_ptr, dfdt_ptr, dt);
    }
  } else {
    // Cell-average mode: restrict only the modal average.
    gkyl_range_iter_init(&iter, range);
    while (gkyl_range_iter_next(&iter)) {
      long idx = gkyl_range_idx(range, iter.idx);
      const double *f_ptr = gkyl_array_cfetch(f, idx);
      double *dfdt_ptr = gkyl_array_fetch(dfdt, idx);
      
      restrict_fdot_avg(up, f_ptr, dfdt_ptr, dt);
    }
  }
}

void
gkyl_positivity_fdot_restrict_release(gkyl_positivity_fdot_restrict* up)
{
  if (up) {
    gkyl_free(up->fquad);
    gkyl_free(up);
  }
}
