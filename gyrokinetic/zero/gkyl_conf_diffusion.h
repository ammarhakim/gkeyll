#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

struct gkyl_conf_diffusion_auxfields {
  /** Full cdim-by-cdim tensor K^{ij}=J D g^{ij}, component-major. */
  const struct gkyl_array *diffusion_tensor;
  /** DG expansion of 1/J, used to estimate the CFL rate for N=J n. */
  const struct gkyl_array *jacobgeo_inv;
};

/**
 * Create a configuration-space diffusion equation object. The implementation
 * supports 1x-3x serendipity p1 bases on the host.
 */
struct gkyl_dg_eqn* gkyl_conf_diffusion_new(const struct gkyl_basis *basis,
  const struct gkyl_range *conf_range, bool use_gpu);

/** Set the mapped diffusion tensor and inverse Jacobian. */
void gkyl_conf_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn,
  struct gkyl_conf_diffusion_auxfields auxin);
