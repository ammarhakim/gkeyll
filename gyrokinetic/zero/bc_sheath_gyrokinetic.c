#include "gkyl_bc_sheath_gyrokinetic.h"
#include "gkyl_bc_sheath_gyrokinetic_priv.h"
#include <gkyl_alloc.h>
#include <gkyl_dg_bin_ops.h>
#include <assert.h>

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void bc_gksheath_update_vcutsq_surrogate(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
   if (up->use_gpu) {
     bc_gksheath_update_vcutsq_surrogate_cu(up, phi, phi_wall, dens, temp, bmag, bimpact_angle, conf_r);
     return;
   }
#endif

  int cidx[GKYL_MAX_CDIM]; // Configuration space index.
  int vidx[GKYL_MAX_VDIM]; // Velocity space index.

  // Set the fixed parallel config index and vpar lower index.
  cidx[up->cdim-1] = up->edge == GKYL_LOWER_EDGE ? up->skin_r->lower[up->cdim-1] : up->skin_r->upper[up->cdim-1];
  vidx[0] = up->skin_r->lower[up->cdim];

  // Loop over perpendicular config cells only to fill the KANN input vectors.
  // build_input writes n_nodes_per_cell nodes into nn_inp_out using stride dim_in:
  //   feature f of local node n -> nn_inp_out[f + dim_in * n]
  // This matches gkyl_kn_vec layout: node g data at data[g*N], feature f at data[g*N + f].
  struct gkyl_range_iter iter_xy;
  gkyl_range_iter_init(&iter_xy, &up->perp_local);
  while (gkyl_range_iter_next(&iter_xy)) {
    for (int d = 0; d < up->cdim-1; d++) {
      cidx[d] = iter_xy.idx[d];
    }    
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *dens_p = (const double*) gkyl_array_cfetch(dens, conf_loc);
    const double *temp_p = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    const double *bimpact_angle_p = (const double*) gkyl_array_cfetch(bimpact_angle, conf_loc);

    // Contiguous block stride for the input parameter of the NN.
    long perp_k = gkyl_range_idx(&up->perp_local, iter_xy.idx);
    int stride = perp_k * up->perp_node_per_cell * up->kann_inp->N;

    up->kernels->vcutsq_input(phi_p, phi_wall_p, dens_p, temp_p, bmag_p, bimpact_angle_p,
      up->kann_inp->N, up->kann_inp->data + stride);
  }

  // Infer NN.
  gkyl_kann_net_apply(up->kann_net, up->kann_inp, up->kann_out);
  
  // Loop over perp config + mu cells to project NN output onto DG vcut_fact.
  // vcut_calc reads local node n's output block as nn_out + n*dim_out, matching
  // gkyl_kn_vec layout where node g data starts at data[g*N].
  struct gkyl_range_iter iter_xymu;
  gkyl_range_iter_init(&iter_xymu, &up->vcutsq_local);
  while (gkyl_range_iter_next(&iter_xymu)) {

    for (int d = 0; d < up->cdim-1; d++)
      cidx[d] = iter_xymu.idx[d];
    vidx[1] = iter_xymu.idx[up->cdim-1];

    long cperp_mu_loc = gkyl_range_idx(&up->vcutsq_local, iter_xymu.idx);
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    long vel_loc = gkyl_range_idx(&up->vel_map->local_vel, vidx);
    // 0-based sequential perp cell index from the perp part of the index.
    long perp_k = gkyl_range_idx(&up->perp_local, cidx);

    const double *vmap_p = (const double*) gkyl_array_cfetch(up->vel_map->vmap, vel_loc);
    const double *temp_p = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    double *vcutsq_p = (double*) gkyl_array_cfetch(up->vcutsq, cperp_mu_loc);

    int stride = perp_k * up->perp_node_per_cell * up->kann_out->N;
    // Compute the DG representation of the interpolated vcutsq values.
    up->kernels->vcutsq_surr(vmap_p, up->kann_out->data + stride, up->kann_out->N, temp_p, bmag_p, vcutsq_p);
  }
}

void bc_gksheath_update_vcutsq_const(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
   if (up->use_gpu) {
     bc_gksheath_update_vcutsq_const_cu(up, phi, phi_wall, dens, temp, bmag, bimpact_angle, conf_r);
     return;
   }
#endif

  int cidx[GKYL_MAX_CDIM]; // Configuration space index.

  // Set the fixed parallel config index and vpar lower index.
  cidx[up->cdim-1] = up->edge == GKYL_LOWER_EDGE ? up->skin_r->lower[up->cdim-1] : up->skin_r->upper[up->cdim-1];

  // Loop over perpendicular config cells to compute vcutsq using the conducting sheath model.
  struct gkyl_range_iter iter_perpmu;
  gkyl_range_iter_init(&iter_perpmu, &up->vcutsq_local);
  while (gkyl_range_iter_next(&iter_perpmu)) {
    
    long perpmu_loc = gkyl_range_idx(&up->vcutsq_local, iter_perpmu.idx);
    double *vcutsq_p = (double*) gkyl_array_cfetch(up->vcutsq, perpmu_loc);

    for (int d = 0; d < up->cdim-1; d++)
      cidx[d] = iter_perpmu.idx[d];
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);

    up->kernels->vcutsq_const(phi_p, phi_wall_p, up->q2Dm, vcutsq_p);
  }
}
struct gkyl_bc_sheath_gyrokinetic*
gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_velocity_map *vel_map,
  int cdim, double q2Dm, bool use_surrogate, const char *surrogate_model_path, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_sheath_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_surrogate = use_surrogate;
  up->use_gpu = use_gpu;
  up->q2Dm = q2Dm;
  up->basis = basis;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;
  up->vel_map = gkyl_velocity_map_acquire(vel_map);
  int vdim = skin_r->ndim - cdim;
  up->update_vcutsq = bc_gksheath_update_vcutsq_const;
  up->vcutsq_dim = cdim-1 + vdim-1;

  // Get polynomial order
  int poly_order;
  if (use_gpu) {
    struct gkyl_basis *basis_ho;
    basis_ho = gkyl_malloc(sizeof(struct gkyl_basis));
    gkyl_cu_memcpy(basis_ho, basis, sizeof(struct gkyl_basis), GKYL_CU_MEMCPY_D2H);
    poly_order = basis_ho->poly_order;
    gkyl_free(basis_ho);
  } else {
    poly_order = basis->poly_order;
  }
  gkyl_cart_modal_serendip(&up->vcutsq_basis, up->vcutsq_dim, poly_order);

  int lower[up->vcutsq_dim], upper[up->vcutsq_dim];
  for (int d=0; d < cdim-1; d++) {
    lower[d] = skin_r->lower[d];
    upper[d] = skin_r->upper[d];
  } 
  lower[up->vcutsq_dim-1] = skin_r->lower[skin_r->ndim-1];
  upper[up->vcutsq_dim-1] = skin_r->upper[skin_r->ndim-1];
  gkyl_range_init(&up->vcutsq_local, up->vcutsq_basis.ndim, lower, upper);

  int perp_lower[GKYL_MAX_CDIM], perp_upper[GKYL_MAX_CDIM];
  for (int d = 0; d < up->cdim-1; d++) {
    perp_lower[d] = up->skin_r->lower[d];
    perp_upper[d] = up->skin_r->upper[d];
  }
  gkyl_range_init(&up->perp_local, up->cdim-1, perp_lower, perp_upper);
  up->vcutsq = mkarr(use_gpu, up->vcutsq_basis.num_basis, up->vcutsq_local.volume);
  
  if (use_surrogate) {
    assert(vdim > 1);
    assert(q2Dm == -2*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS); // The surrogate is only valid for electrons.
    // Assert existence of surrogate model path.
    FILE *fp = fopen(surrogate_model_path, "r");
    if (fp == NULL) {
      gkyl_exit("surrogate_model_path does not point to a valid file");
    }
    up->perp_node_per_cell = 1 << (up->cdim-1);
    
    // gkyl KANN wrapper
    up->kann_net = gkyl_kann_net_load(surrogate_model_path, up->use_gpu);
    int dim_in = gkyl_kann_net_dim_in(up->kann_net);
    int dim_out = gkyl_kann_net_dim_out(up->kann_net);
    assert(dim_in == 3 && dim_out == SRGRZ_N_MU);
    int nperp_nodes = up->perp_node_per_cell * up->perp_local.volume;
    up->kann_inp = up->use_gpu ? gkyl_kn_vec_cu_dev_new(nperp_nodes, dim_in) : gkyl_kn_vec_new(nperp_nodes, dim_in);
    up->kann_out = up->use_gpu ? gkyl_kn_vec_cu_dev_new(nperp_nodes, dim_out) : gkyl_kn_vec_new(nperp_nodes, dim_out);
    up->kann_infer_xy_out = mkarr(up->use_gpu, dim_out*pow(2, up->cdim-1), up->skin_r->volume);
    up->nn_out = up->use_gpu ? gkyl_cu_malloc(dim_out*pow(2, up->cdim-1)*sizeof(double)) 
      : gkyl_malloc(dim_out*pow(2, up->cdim-1)*sizeof(double));

    up->update_vcutsq = bc_gksheath_update_vcutsq_surrogate;
  }
  
  // Choose the kernels that does the reflection/no reflection/partial reflection, 
  // and surrogate kernels if enabled.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
    gkyl_bc_gksheath_choose_reflectedf_kernel_cu(basis, edge, up->kernels_cu);
    if (use_surrogate) {
      gkyl_bc_gksheath_choose_surrogate_kernels_cu(basis, edge, up->kernels_cu);
    } else {
      gkyl_bc_gksheath_choose_const_kernels_cu(basis, edge, up->kernels_cu);
    }
  } else {
    up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
    assert(up->kernels->reflectedf);
    if (use_surrogate) {
      bc_gksheath_choose_vcutsq_surr_kernels(basis, edge, up->kernels);
    } else {
      bc_gksheath_choose_vcutsq_const_kernels(basis, edge, up->kernels);
    }
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
  assert(up->kernels->reflectedf);
  if (use_surrogate) {
    bc_gksheath_choose_vcutsq_surr_kernels(basis, edge, up->kernels);
  } else {
    bc_gksheath_choose_vcutsq_const_kernels(basis, edge, up->kernels);
  }
  up->kernels_cu = up->kernels;
#endif

  return up;
}

/* Modeled after gkyl_array_flip_copy_to_buffer_fn */
void
gkyl_bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_sheath_gyrokinetic_advance_cu(up, phi, phi_wall, distf, conf_r);
    return;
  }
#endif

  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int vidx[2];
  int vcutsq_idx[up->vcutsq_local.ndim];

  int pdim = up->skin_r->ndim; 
  int vpar_dir = up->cdim;
  int uplo = up->skin_r->upper[vpar_dir]+up->skin_r->lower[vpar_dir];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->skin_r);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_copy_int_arr(pdim, iter.idx, fidx);
    fidx[vpar_dir] = uplo - iter.idx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[up->dir] = up->ghost_r->lower[up->dir];

    long skin_loc = gkyl_range_idx(up->skin_r, iter.idx);
    long ghost_loc = gkyl_range_idx(up->ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);
    double *out = (double*) gkyl_array_fetch(distf, ghost_loc);

    for (int d=up->cdim; d<pdim; d++) vidx[d-up->cdim] = iter.idx[d]; 
    long conf_loc = gkyl_range_idx(conf_r, iter.idx);
    long vel_loc = gkyl_range_idx(&up->vel_map->local_vel, vidx);

    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *vmap_p = (const double*) gkyl_array_cfetch(up->vel_map->vmap, vel_loc);
    
    int vcutsq_dim = up->vcutsq_local.ndim;
    for (int d=0; d<vcutsq_dim-1; d++) vcutsq_idx[d] = iter.idx[d]; // config space perp directions.
    vcutsq_idx[vcutsq_dim-1] = iter.idx[pdim-1]; // mu direction.
    long vcutsq_loc = gkyl_range_idx(&up->vcutsq_local, vcutsq_idx);
    const double *vcutsq_p = (const double*) gkyl_array_cfetch(up->vcutsq, vcutsq_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->kernels->reflectedf(vmap_p, vcutsq_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void gkyl_bc_sheath_gyrokinetic_set_vcutsq(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *vcutsq)
{
  gkyl_array_copy_range(up->vcutsq, vcutsq, &up->vcutsq_local);
}

struct gkyl_array* gkyl_bc_sheath_gyrokinetic_acquire_vcutsq(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return gkyl_array_acquire(up->vcutsq);
}

struct gkyl_basis gkyl_bc_sheath_gyrokinetic_get_vcutsq_basis(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return up->vcutsq_basis;
}

struct gkyl_range* gkyl_bc_sheath_gyrokinetic_get_vcutsq_range(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return &up->vcutsq_local;
}

struct gkyl_kann_net *gkyl_bc_sheath_gyrokinetic_acquire_model(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return gkyl_kann_net_acquire(up->kann_net);
}

void gkyl_bc_sheath_gyrokinetic_update_vcutsq(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
  return up->update_vcutsq(up, phi, phi_wall, dens, temp, bmag, bimpact_angle, conf_r);
}

void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
  }
#endif
  gkyl_free(up->kernels);
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->vcutsq);

  if (up->use_surrogate) {
    gkyl_kann_net_release(up->kann_net);
    gkyl_kn_vec_release(up->kann_inp);
    gkyl_kn_vec_release(up->kann_out);
    gkyl_array_release(up->kann_infer_xy_out);
    if (up->use_gpu) {
      gkyl_cu_free(up->nn_out);
    } else {
      gkyl_free(up->nn_out);
    }
  }
  gkyl_free(up);
}