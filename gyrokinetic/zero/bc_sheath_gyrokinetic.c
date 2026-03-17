#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
#include <gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h>
#include <gkyl_alloc.h>
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

void bc_gksheath_update_vcut_fact_surrogate_enabled(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp, double q, double m,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
   // Not implemented on GPU yet.
#endif

  int cidx[GKYL_MAX_CDIM]; // Configuration space index.
  int vidx[GKYL_MAX_VDIM]; // Velocity space index.
  int vcut_idx[GKYL_MAX_CDIM]; // Index for vcut_fact array, which is in the space of perpendicular config space coords and mu.

  int pdim = up->skin_r->ndim; 
  int vpar_dir = up->cdim;

  double q2Dm = 2*q/m;
  
  // Set the independent indices.
  cidx[up->cdim-1] = up->edge == GKYL_LOWER_EDGE ? up->skin_r->lower[up->cdim-1] : up->skin_r->upper[up->cdim-1];
  vidx[0] = up->skin_r->lower[up->cdim];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->vcut_fact_local);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d < up->cdim-1; d++) {
      cidx[d] = iter.idx[d];
      vcut_idx[d] = iter.idx[d];
    }
    vidx[1] = iter.idx[up->cdim-1];
    vcut_idx[up->vcut_fact_dim-1] = iter.idx[up->vcut_fact_dim-1];

    long cperp_mu_loc = gkyl_range_idx(&up->vcut_fact_local, iter.idx);
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    long vel_loc = gkyl_range_idx(&up->vel_map->local_vel, vidx);

    const double *vmap_p = (const double*) gkyl_array_cfetch(up->vel_map->vmap, vel_loc);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *dens_p = (const double*) gkyl_array_cfetch(dens, conf_loc);
    const double *temp_p = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    const double *bimpact_angle_p = (const double*) gkyl_array_cfetch(bimpact_angle, conf_loc);
    double *vcut_fact_p = (double*) gkyl_array_cfetch(up->vcut_fact, cperp_mu_loc);

    up->kernels->surrogate(vmap_p, phi_p, phi_wall_p, dens_p, temp_p, q2Dm, bmag_p, bimpact_angle_p, vcut_fact_p);
  }
}

void bc_gksheath_update_vcut_fact_surrogate_disabled(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp, double q, double m,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
  return;
}

struct gkyl_bc_sheath_gyrokinetic*
gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_velocity_map *vel_map,
  int cdim, double q2Dm, bool use_surrogate, bool use_gpu)
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

  // Sheath physics surrogate setup.
  up->update_vcut_fact = bc_gksheath_update_vcut_fact_surrogate_disabled;
  // Setup vcut_fact to 1 by default.
  up->vcut_fact_dim = cdim-1 + vdim-1;

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
  gkyl_cart_modal_serendip(&up->vcut_fact_basis, up->vcut_fact_dim, poly_order);

  int lower[up->vcut_fact_dim], upper[up->vcut_fact_dim];
  for (int d=0; d < cdim-1; d++) {
    lower[d] = skin_r->lower[d];
    upper[d] = skin_r->upper[d];
  } 
  lower[up->vcut_fact_dim-1] = skin_r->lower[skin_r->ndim-1];
  upper[up->vcut_fact_dim-1] = skin_r->upper[skin_r->ndim-1];
  gkyl_range_init(&up->vcut_fact_local, up->vcut_fact_basis.ndim, lower, upper);
  up->vcut_fact = mkarr(use_gpu, up->vcut_fact_basis.num_basis, up->vcut_fact_local.volume);
  double dg_norm = pow(sqrt(2), up->vcut_fact_basis.ndim);
  gkyl_array_shiftc_range(up->vcut_fact, dg_norm, 0, &up->vcut_fact_local);
  if (use_surrogate) {
    assert(vdim > 1);
    assert(q2Dm == -2*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS); // The surrogate is only valid for electrons.
    up->update_vcut_fact = bc_gksheath_update_vcut_fact_surrogate_enabled;
  }
  
  // Choose the kernels that does the reflection/no reflection/partial reflection, 
  // and surrogate kernels if enabled.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
    gkyl_bc_gksheath_choose_reflectedf_kernel_cu(basis, edge, up->kernels_cu);
    if (use_surrogate)
      gkyl_bc_gksheath_choose_surrogate_kernel_cu(basis, edge, up->kernels_cu);
  } else {
    up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
    assert(up->kernels->reflectedf);
    if (use_surrogate) {
      up->kernels->surrogate = bc_gksheath_choose_surrogate_kernel(basis, edge);
      assert(up->kernels->surrogate);
    }
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
  assert(up->kernels->reflectedf);
  if (use_surrogate) {
    up->kernels->surrogate = bc_gksheath_choose_surrogate_kernel(basis, edge);
    assert(up->kernels->surrogate);
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
  int vcut_fact_idx[up->vcut_fact_local.ndim];

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
    
    for (int d=0; d<up->cdim-1; d++) vcut_fact_idx[d] = iter.idx[d]; // config space perp directions.
    vcut_fact_idx[up->cdim-1] = iter.idx[pdim-1]; // mu direction.
    long vcut_fact_loc = gkyl_range_idx(&up->vcut_fact_local, vcut_fact_idx);
    const double *vcut_fact_p = (const double*) gkyl_array_cfetch(up->vcut_fact, vcut_fact_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->kernels->reflectedf(vmap_p, up->q2Dm, phi_p, phi_wall_p, vcut_fact_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void gkyl_bc_sheath_gyrokinetic_set_vcut_fact(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *vcut_fact)
{
  gkyl_array_copy_range(up->vcut_fact, vcut_fact, &up->vcut_fact_local);
}

struct gkyl_array* gkyl_bc_sheath_gyrokinetic_get_vcut_fact(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return gkyl_array_acquire(up->vcut_fact);
}

struct gkyl_basis gkyl_bc_sheath_gyrokinetic_get_vcut_fact_basis(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return up->vcut_fact_basis;
}

struct gkyl_range* gkyl_bc_sheath_gyrokinetic_get_vcut_fact_range(struct gkyl_bc_sheath_gyrokinetic *up)
{
  return &up->vcut_fact_local;
}

void gkyl_bc_sheath_gyrokinetic_update_vcut_fact_surrogate(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp, double q, double m,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r)
{
  return up->update_vcut_fact(up, phi, phi_wall, dens, temp, q, m, bmag, bimpact_angle, conf_r);
}

void gkyl_bc_sheath_gyrokinetic_evaluate_vcut_fact_surrogate(const double *mu_new,  int n, double phi, double phi_wall,
    double dens_e, double temp_e, double q2Dm, double bmag, double bimpact_angle, double *out)
{
  return bc_sheath_gyrokinetic_srgrz_eval_physical_vcut_fact(mu_new, n, phi, phi_wall, dens_e, temp_e, q2Dm, bmag, bimpact_angle, out);
}

void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
  }
#endif
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->vcut_fact);
  gkyl_free(up->kernels);
  gkyl_free(up);
}