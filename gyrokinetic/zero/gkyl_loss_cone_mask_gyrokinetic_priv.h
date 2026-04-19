// Private header: not for direct use
#pragma once

#include <float.h>
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
    int place_value = (int) pow(2.0, (double) (ndim - 1 - d));
    int digit = (node / place_value) % 2;
    x[d] = digit ? 1.0 : -1.0;
  }
}

struct gkyl_loss_cone_mask_gyrokinetic {
  int cdim; // Configuration-space dimension.
  int num_basis_conf; // Number of configuration-space basis functions.
  int num_nodes_conf; // Number of configuration-space cell nodes.
  bool use_gpu;

  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  double mass; // Species mass.
  double charge; // Species charge.

  struct gkyl_array *basis_at_nodes_conf; // Basis functions at configuration-space nodes.
};

GKYL_CU_DH
static inline int
conf_node_z_endpoint_index(int cdim, int conf_node, int zdim)
{
  double eta[GKYL_MAX_DIM] = { 0.0 };
  nodal_coords(cdim, conf_node, eta);
  return eta[zdim] > 0.0 ? 1 : 0;
}

GKYL_CU_DH
static inline int
conf_node_with_matching_perpendicular_coords(int cdim, int anchor_node, int zdim,
  int z_endpoint_index)
{
  int num_nodes_conf = (int)pow(2.0, (double)cdim);

  double eta_anchor[GKYL_MAX_DIM] = { 0.0 };
  nodal_coords(cdim, anchor_node, eta_anchor);

  for (int cand = 0; cand < num_nodes_conf; ++cand) {
    if (conf_node_z_endpoint_index(cdim, cand, zdim) != z_endpoint_index) {
      continue;
    }

    double eta_cand[GKYL_MAX_DIM] = { 0.0 };
    nodal_coords(cdim, cand, eta_cand);

    bool same_transverse = true;
    for (int d = 0; d < cdim; ++d) {
      if (d == zdim) {
        continue;
      }
      if (eta_cand[d] != eta_anchor[d]) {
        same_transverse = false;
        break;
      }
    }

    if (same_transverse) {
      return cand;
    }
  }

  return anchor_node;
}

GKYL_CU_DH
static inline double
field_node_val(const struct gkyl_array *arr, const struct gkyl_array *basis_at_nodes,
  int num_basis, long linidx, int node)
{
  const double *arr_d = (const double *)gkyl_array_cfetch(arr, linidx);
  const double *basis_d = (const double *)gkyl_array_cfetch(basis_at_nodes, node);

  double val = 0.0;
  for (int k = 0; k < num_basis; ++k) {
    val += arr_d[k] * basis_d[k];
  }

  return val;
}

GKYL_CU_DH
static inline void
escape_barriers(int cdim, int num_basis_conf, const struct gkyl_range *conf_range,
  const struct gkyl_array *basis_at_nodes_conf, const struct gkyl_array *phi,
  const struct gkyl_array *bmag, const int *base_idx, int z_cell,
  int anchor_conf_node, double mu, double charge,
  double *barrier_left, double *barrier_right)
{
  int zdim = cdim - 1;

  int z_endpoint_index = conf_node_z_endpoint_index(cdim, anchor_conf_node, zdim);
  int anchor_node = conf_node_with_matching_perpendicular_coords(cdim, anchor_conf_node, zdim,
    z_endpoint_index);
  int z_upper_node = conf_node_with_matching_perpendicular_coords(cdim, anchor_conf_node,
    zdim, 1);
  int z_lower_node = conf_node_with_matching_perpendicular_coords(cdim, anchor_conf_node,
    zdim, 0);

  int scan_idx[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    scan_idx[d] = base_idx[d];
  }

  *barrier_left = -DBL_MAX;
  *barrier_right = -DBL_MAX;

  for (int iz = conf_range->lower[zdim]; iz <= conf_range->upper[zdim]; ++iz) {
    scan_idx[zdim] = iz;
    long linidx = gkyl_range_idx(conf_range, scan_idx);

    int z_scan_node = anchor_node;
    if (iz < z_cell) {
      z_scan_node = z_upper_node;
    }
    else if (iz > z_cell) {
      z_scan_node = z_lower_node;
    }

    double phi_scan = field_node_val(phi, basis_at_nodes_conf, num_basis_conf,
      linidx, z_scan_node);
    double bmag_scan = field_node_val(bmag, basis_at_nodes_conf, num_basis_conf,
      linidx, z_scan_node);
    double u_scan = mu * bmag_scan + charge * phi_scan;

    if (iz <= z_cell && u_scan > *barrier_left) {
      *barrier_left = u_scan;
    }
    if (iz >= z_cell && u_scan > *barrier_right) {
      *barrier_right = u_scan;
    }
  }
}


#ifdef GKYL_HAVE_CUDA
void gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi, struct gkyl_array *mask_out);
#endif