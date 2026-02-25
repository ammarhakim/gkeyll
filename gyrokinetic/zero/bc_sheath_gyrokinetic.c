#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

// To remove
#include <gkyl_proj_on_basis.h>

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

void
eval_func_alpha_mu(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[0], mu = xn[1];

  double Lmu = 3.0; // Characteristic scale length in mu direction.
  double alpha_mu_0 = 1.0;
  double dg_const = pow(sqrt(2.0), 2); // Normalization constant for gyrokinetic velocity space.

  double alpha_mu = alpha_mu_0 * (1 + 3*exp(-mu/Lmu));
  fout[0] =  pow(alpha_mu, 2)/dg_const;
}

struct gkyl_bc_sheath_gyrokinetic*
gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_velocity_map *vel_map,
  int cdim, double q2Dm, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_sheath_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->q2Dm = q2Dm;
  up->basis = basis;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;
  up->vel_map = gkyl_velocity_map_acquire(vel_map);

  gkyl_cart_modal_gkhybrid(&up->alpha_mu_basis, 0, vel_map->local_vel.ndim);
  up->alpha_mu = mkarr(use_gpu, up->alpha_mu_basis.num_basis, vel_map->local_vel.volume);
  // Set the alpha_mu array to 1 by default (constant vcut).
  double dg_norm = pow(sqrt(2), up->alpha_mu_basis.ndim);
  gkyl_array_shiftc_range(up->alpha_mu, dg_norm, 0, &vel_map->local_vel);
  // Test a function on alpha_mu
  gkyl_proj_on_basis *projalphamu = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vel_map->grid_vel,
      .basis = &up->alpha_mu_basis,
      .num_ret_vals = 1,
      .eval = eval_func_alpha_mu,
      .ctx = NULL,
    }
  );
  // gkyl_proj_on_basis_advance(projalphamu, 0.0, &vel_map->local_vel, up->alpha_mu);
  gkyl_proj_on_basis_release(projalphamu);

  // Choose the kernel that does the reflection/no reflection/partial
  // reflection.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
    gkyl_bc_gksheath_choose_reflectedf_kernel_cu(basis, edge, up->kernels_cu);
  } else {
    up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
    assert(up->kernels->reflectedf);
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
  assert(up->kernels->reflectedf);
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
    const double *alpha_mu_p = (const double*) gkyl_array_cfetch(up->alpha_mu, vel_loc);
    const double *vmap_p = (const double*) gkyl_array_cfetch(up->vel_map->vmap, vel_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->kernels->reflectedf(vmap_p, up->q2Dm, phi_p, phi_wall_p, alpha_mu_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void gkyl_bc_sheath_gyrokinetic_set_alpha_mu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *alpha_mu)
{
  gkyl_array_copy_range(up->alpha_mu, alpha_mu, &up->vel_map->local_vel);
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
  gkyl_array_release(up->alpha_mu);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
