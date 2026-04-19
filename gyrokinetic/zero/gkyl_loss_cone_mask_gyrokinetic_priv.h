// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH
static inline void
nodal_coords(int ndim, int node, double *x)
{
  for (int d = 0; d < ndim; ++d) {
    int bit = ndim - 1 - d;
    x[d] = ((node >> bit) & 1) ? 1.0 : -1.0;
  }
}

struct gkyl_loss_cone_mask_gyrokinetic {
  int cdim; // Configuration-space dimension.
  int num_basis_conf; // Number of configuration-space basis functions.
  int num_nodes_conf; // Number of configuration-space cell nodes.
  bool use_gpu;

  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  const struct gkyl_array *bmag; // Magnetic field magnitude (cdim DG expansion).

  double mass; // Species mass.
  double charge; // Species charge.

  struct gkyl_array *basis_at_nodes_conf; // Basis functions at configuration-space nodes.
};
