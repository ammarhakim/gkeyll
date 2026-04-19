#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>

// The loss-cone mask is built from the escape barrier
// EB(z,mu) = min( max_{s in [z_L,z]} U(s,mu), max_{s in [z,z_R]} U(s,mu) )
// with U = mu*B + q*phi. A node is trapped if H < EB.

static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
    : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static int
init_node_values(int ndim, const struct gkyl_basis *basis,
  struct gkyl_array **basis_at_nodes, bool use_gpu)
{
  int num_nodes = 1 << ndim;
  struct gkyl_array *basis_at_nodes_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis,
    num_nodes);
  *basis_at_nodes = mkarr(basis->num_basis, num_nodes, use_gpu);

  double eta[GKYL_MAX_DIM] = { 0.0 };
  for (int node = 0; node < num_nodes; ++node) {
    nodal_coords(ndim, node, eta);
    basis->eval(eta, gkyl_array_fetch(basis_at_nodes_ho, node));
  }

  gkyl_array_copy(*basis_at_nodes, basis_at_nodes_ho);
  gkyl_array_release(basis_at_nodes_ho);

  return num_nodes;
}

static inline double
field_node_val(const struct gkyl_array *arr, const struct gkyl_array *basis_at_nodes,
  int num_basis, long linidx, int node)
{
  const double *arr_d = gkyl_array_cfetch(arr, linidx);
  const double *basis_d = gkyl_array_cfetch(basis_at_nodes, node);

  double val = 0.0;
  for (int k = 0; k < num_basis; ++k) {
    val += arr_d[k] * basis_d[k];
  }

  return val;
}

static inline void
escape_barriers(const gkyl_loss_cone_mask_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *bmag, const struct gkyl_range *conf_range, const int *base_idx,
  int z_cell, int target_node, int target_z_bit, double mu, double charge,
  double *barrier_left, double *barrier_right)
{
  int scan_idx[GKYL_MAX_DIM];
  for (int d = 0; d < conf_range->ndim; ++d) {
    scan_idx[d] = base_idx[d];
  }

  *barrier_left = -DBL_MAX;
  *barrier_right = -DBL_MAX;

  for (int iz = conf_range->lower[conf_range->ndim - 1]; iz <= conf_range->upper[conf_range->ndim - 1]; ++iz) {
    scan_idx[conf_range->ndim - 1] = iz;
    long linidx = gkyl_range_idx(conf_range, scan_idx);

    int left_node = target_node;
    int right_node = target_node;
    if (iz == z_cell) {
      if (target_z_bit) {
        left_node |= 1;
        right_node |= 1;
      }
      else {
        left_node &= ~1;
        right_node &= ~1;
      }
    }
    else if (iz < z_cell) {
      left_node |= 1;
      right_node &= ~1;
    }
    else {
      left_node |= 1;
      right_node &= ~1;
    }

    if (iz <= z_cell) {
      double phi_left = field_node_val(phi, up->basis_at_nodes_conf, up->num_basis_conf,
        linidx, left_node);
      double bmag_left = field_node_val(bmag, up->basis_at_nodes_conf, up->num_basis_conf,
        linidx, left_node);
      double u_left = mu * bmag_left + charge * phi_left;
      if (u_left > *barrier_left) {
        *barrier_left = u_left;
      }
    }
    if (iz >= z_cell) {
      double phi_right = field_node_val(phi, up->basis_at_nodes_conf, up->num_basis_conf,
        linidx, right_node);
      double bmag_right = field_node_val(bmag, up->basis_at_nodes_conf, up->num_basis_conf,
        linidx, right_node);
      double u_right = mu * bmag_right + charge * phi_right;
      if (u_right > *barrier_right) {
        *barrier_right = u_right;
      }
    }
  }
}

struct gkyl_loss_cone_mask_gyrokinetic*
gkyl_loss_cone_mask_gyrokinetic_inew(const struct gkyl_loss_cone_mask_gyrokinetic_inp *inp)
{
  gkyl_loss_cone_mask_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);
  up->mass = inp->mass;
  up->charge = inp->charge;
  up->use_gpu = inp->use_gpu;

  up->cdim = inp->conf_basis->ndim;
  up->num_basis_conf = inp->conf_basis->num_basis;
  up->num_nodes_conf = 1 << up->cdim;

  init_node_values(up->cdim, inp->conf_basis, &up->basis_at_nodes_conf, inp->use_gpu);

  return up;
}

void
gkyl_loss_cone_mask_gyrokinetic_advance(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *bmag, const struct gkyl_array *phi, struct gkyl_array *mask_out)
{
  int cdim = up->cdim;
  int pdim = phase_range->ndim;
  int vdim = pdim - cdim;
  int num_basis_conf = up->num_basis_conf;
  int num_phase_nodes = 1 << pdim;

  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    int conf_idx[GKYL_MAX_DIM] = { 0 };
    int vel_idx[GKYL_MAX_DIM] = { 0 };
    for (int d = 0; d < cdim; ++d) {
      conf_idx[d] = phase_iter.idx[d];
    }
    for (int d = 0; d < vdim; ++d) {
      vel_idx[d] = phase_iter.idx[cdim + d];
    }

    const struct gkyl_velocity_map *gvm = up->vel_map;
    long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_idx);
    const double *vmap_d = gkyl_array_cfetch(gvm->vmap, linidx_vel);

    bool cell_trapped = true;
    for (int node = 0; node < num_phase_nodes && cell_trapped; ++node) {
      int conf_node = node / (1 << vdim);
      int vel_node = node % (1 << vdim);

      double xmu[GKYL_MAX_DIM] = { 0.0 };
      for (int vd = 0; vd < vdim; ++vd) {
        double vel_eta[GKYL_MAX_DIM] = { 0.0 };
        nodal_coords(vdim, vel_node, vel_eta);
        double xcomp[1] = { vel_eta[vd] };
        xmu[cdim + vd] = gvm->vmap_basis->eval_expand(xcomp,
          vmap_d + vd * gvm->vmap_basis->num_basis);
      }

      double mu = xmu[cdim + 1];
      double vpar = xmu[cdim];

      long linidx_conf = gkyl_range_idx(conf_range, conf_idx);
      double bmag_curr = field_node_val(bmag, up->basis_at_nodes_conf, num_basis_conf,
        linidx_conf, conf_node);
      double phi_curr = field_node_val(phi, up->basis_at_nodes_conf, num_basis_conf,
        linidx_conf, conf_node);
      double h_curr = 0.5 * up->mass * vpar * vpar + mu * bmag_curr + up->charge * phi_curr;

      int zdim = cdim - 1;
      int target_z_bit = conf_node & 1;

      double barrier_left, barrier_right;
      escape_barriers(up, phi, bmag, conf_range, conf_idx, conf_idx[zdim], conf_node,
        target_z_bit, mu, up->charge, &barrier_left, &barrier_right);

      cell_trapped = h_curr < GKYL_MIN2(barrier_left, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(phase_range, phase_iter.idx);
    double *mask_d = gkyl_array_fetch(mask_out, linidx_phase);
    mask_d[0] = cell_trapped ? 1.0 : 0.0;
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_release(gkyl_loss_cone_mask_gyrokinetic *up)
{
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->basis_at_nodes_conf);
  gkyl_free(up);
}
