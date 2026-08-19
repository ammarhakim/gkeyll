#include <gkyl_gyrokinetic_priv.h>

// Take time-step using the SSP-RK3 method for the hyperbolic components
// Then, use the actual explicit timestep to update implicit BGK collisions
// and fluid-neutral diffusion.
struct gkyl_update_status
gyrokinetic_update_op_split(gkyl_gyrokinetic_app* app, double dt0)
{
  struct gkyl_update_status st = gyrokinetic_update_ssp_rk3(app, dt0);

  // Take the operator-split implicit timestep.
  gyrokinetic_update_implicit_coll(app, st.dt_actual);

  return st;
}
