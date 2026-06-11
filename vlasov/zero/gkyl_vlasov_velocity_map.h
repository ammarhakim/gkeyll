#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

// Forward declarations to avoid pulling in comm/IO headers.
struct gkyl_comm;
struct gkyl_msgpack_data;

struct gkyl_vlasov_velocity_map_inp {
  evalf_t eval_vmap; // Velocity mapping in this velocity-space dimension used by dg_basis_ops to construct C^1 mapping. NULL => identity.
  void *ctx; // Context for function evaluation. Can be NULL.
};

// Representation of the stored velocity map, chosen by the velocity basis:
// tensor bases (p>1) use the C^1 cubic; Serendipity bases (and p=1, where
// tensor and Serendipity coincide) use the C^0 piecewise linear, stored as a
// degenerate cubic in the same vdim*4 layout (quadratic and cubic
// coefficients zero). The C^0 map's Jacobian is piecewise constant and
// discontinuous at cell interfaces: the Jacobian arrays hold that constant
// replicated at every quadrature point, and consumers needing the Jacobian
// at an interface must take the minimum of the two adjacent cells when
// estimating maximum frequencies.
enum gkyl_vlasov_vmap_rep {
  GKYL_VLASOV_VMAP_C1_CUBIC = 0,
  GKYL_VLASOV_VMAP_C0_LINEAR,
};

// Velocity-space mapping object for the Vlasov app. Packages the
// representation (C^1 cubic or C^0 linear, by velocity basis) of the
// (possibly non-uniform) computational-to-physical velocity mapping in each
// velocity dimension, along with the velocity-space Jacobian evaluated at the
// Gauss-Legendre quadrature points needed by volume, surface, and projection
// operations. This object is *always* created by the app: when no user
// mapping is supplied the identity map is used, and there is a single code
// path — the map is always populated, always written at frame 0, and the
// Jacobian is always divided/rescaled — so consumers handle uniform and
// non-uniform velocity grids transparently.
//
// This is a host-side, reference-counted container. There is deliberately no
// device-side mirror of the struct itself: when created with use_gpu=true the
// solver arrays (vmap, jacob_vel, jacob_vel_surf, jacob_vel_gauss) are
// device-resident gkyl_arrays, and consumers acquire the struct (for lifetime
// safety) and unpack raw device pointers (e.g. vvm->jacob_vel->on_dev) inside
// their _cu.cu initialization.
struct gkyl_vlasov_velocity_map {
  struct gkyl_rect_grid grid_vel; // Velocity-space grid.
  struct gkyl_range local_vel; // Velocity-space range the arrays are defined on.
  struct gkyl_basis basis_vel; // Velocity-space basis (b_type and poly_order drive kernel dispatch).
  enum gkyl_vlasov_vmap_rep rep; // Representation of the stored map.
  bool is_identity; // True if no user mapping was given in any direction.
  bool is_mapped; // Always true: the map is always created and consumed.
                  // Retained for consumers that set use_vmap from it; to be
                  // removed once those call sites are simplified.

  // Solver arrays; device-resident when created with use_gpu=true.
  struct gkyl_array *vmap; // C^1 cubic representation of mapping in each velocity dimension (vdim*4 components: vx, then vy, then vz).
  struct gkyl_array *jacob_vel; // Jacobian (derivative of vmap) in each velocity dimension at 1V Gauss-Legendre quadrature points (vdim*(p+1) components).
  struct gkyl_array *jacob_vel_surf; // Jacobian in each velocity dimension at the (higher order) 1V Gauss-Legendre quadrature points used by surface updates (vdim*(p+2) components).
  struct gkyl_array *jacob_vel_gauss; // Total velocity-space Jacobian at the full Gauss-Legendre quadrature points (tensor(vdim,p).num_basis components).

  // Host mirrors of the solver arrays (acquired aliases when !use_gpu).
  struct gkyl_array *vmap_host;
  struct gkyl_array *jacob_vel_host;
  struct gkyl_array *jacob_vel_surf_host;
  struct gkyl_array *jacob_vel_gauss_host;

  // I/O-only arrays; always host-resident.
  struct gkyl_array *vmap_pgkyl_host; // C^1 cubic representation of mapping for I/O (defined in the full 1V, 2V, or 3V).
  struct gkyl_array *vmap_avg_pgkyl_host; // Cell average of the mapping for I/O.

  uint32_t flags;
  struct gkyl_ref_count ref_count;
};

/**
 * Create a new velocity map based on the input mapping functions, using a
 * C^1 cubic representation for tensor velocity bases (p>1) and a C^0
 * piecewise linear representation for Serendipity velocity bases (and p=1).
 * Directions with eval_vmap==NULL get the identity map, so this constructor
 * is valid for every species and model. All arrays (including host mirrors
 * and I/O arrays) are allocated and populated here; when use_gpu=true the
 * solver arrays are copied onto device.
 *
 * @param vgrid Velocity-space grid object.
 * @param vrange Velocity-space range the arrays are defined on.
 * @param vel_basis Velocity-space basis (b_type and poly_order drive kernel dispatch).
 * @param inp_vmap Velocity mapping input (function and context) in each velocity-space dimension.
 * @param use_gpu Whether to place the solver arrays on device.
 * @return New velocity map object.
 */
struct gkyl_vlasov_velocity_map* gkyl_vlasov_velocity_map_new(
  const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange,
  const struct gkyl_basis *vel_basis,
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM],
  bool use_gpu);

/**
 * Acquire a pointer to the velocity map object, incrementing its reference
 * count. Release with gkyl_vlasov_velocity_map_release.
 *
 * @param vvm Velocity map object.
 * @return Acquired velocity map object.
 */
struct gkyl_vlasov_velocity_map* gkyl_vlasov_velocity_map_acquire(
  const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Release the velocity map object, decrementing its reference count and
 * freeing all contained arrays when the count reaches zero.
 *
 * @param vvm Velocity map object.
 */
void gkyl_vlasov_velocity_map_release(const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Check whether the solver arrays in this velocity map are device-resident.
 *
 * @param vvm Velocity map object.
 * @return True if the solver arrays live on device.
 */
bool gkyl_vlasov_velocity_map_is_cu_dev(const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Write the velocity map to %s-%s_vmap.gkyl and, if write_cell_avg, its cell
 * average to %s-%s_vmap_avg.gkyl, on rank 0 of comm. The map is static in
 * time so this is intended to be called once (e.g. at frame 0); it is always
 * written, uniform (identity) grids included.
 *
 * @param vvm Velocity map object.
 * @param comm Communicator of the species this map belongs to.
 * @param mt Msgpack metadata to attach to the files. Can be NULL.
 * @param app_name Name of the app.
 * @param species_name Name of the species.
 * @param write_cell_avg Whether to also write the cell-averaged map.
 */
void gkyl_vlasov_velocity_map_write(const struct gkyl_vlasov_velocity_map *vvm,
  struct gkyl_comm *comm, struct gkyl_msgpack_data *mt,
  const char *app_name, const char *species_name, bool write_cell_avg);

/**
 * Divide out the velocity-space Jacobian from Jf to obtain f for use in the
 * DG update. The division is performed at Gauss-Legendre quadrature points
 * because, by construction, we know Jv at these points and Jv is static in
 * time, so we obtain the weakly equivalent f from Jf in our DG expansion
 * (since we can always transform between the Gauss-Legendre nodal basis and
 * our modal basis). For the C^0 linear representation (Serendipity bases) the
 * Jacobian is constant in the cell so the division is exact, modally and
 * nodally. Runs on device when the map was created with use_gpu=true.
 *
 * @param vvm Velocity map object.
 * @param conf_basis Configuration-space basis.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param Jf Input array Jf.
 * @param f_no_J Output array f with velocity-space Jacobian divided out.
 */
void gkyl_vlasov_velocity_map_divide_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J);

/**
 * Multiply f by the velocity-space Jacobian to obtain Jf for use in the DG
 * update. Performed at Gauss-Legendre quadrature points; the inverse of
 * gkyl_vlasov_velocity_map_divide_jacobvel.
 *
 * @param vvm Velocity map object.
 * @param conf_basis Configuration-space basis.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param f_no_J Input array f with velocity-space Jacobian divided out.
 * @param Jf Output array Jf.
 */
void gkyl_vlasov_velocity_map_rescale_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf);
