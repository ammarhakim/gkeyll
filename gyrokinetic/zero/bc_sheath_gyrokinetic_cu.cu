/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
}

// CUDA kernel to set device pointers to kernel that computes the reflected f.
__global__ static void
gkyl_bc_gksheath_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type b_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->reflectedf = ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_bc_gksheath_choose_reflectedf_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  gkyl_bc_gksheath_set_cu_ker_ptrs<<<1,1>>>(basis, edge, kers);
}

// CUDA kernel to set device pointers to kernel that calls the surrogate.
__global__ static void
gkyl_bc_gksheath_surr_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, bool use_conv_check, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type b_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->surrogate = ser_sheath_surrogate_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_bc_gksheath_choose_surrogate_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, bool use_conv_check, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  gkyl_bc_gksheath_surr_set_cu_ker_ptrs<<<1,1>>>(basis, edge, use_conv_check, kers);
}

__global__ static void
gkyl_bc_sheath_gyrokinetic_advance_cu_ker(int cdim, int dir, const struct gkyl_range skin_r, 
  const struct gkyl_range ghost_r, const struct gkyl_range conf_r, const struct gkyl_range vcut_r, 
  const struct gkyl_range vel_r, const struct gkyl_basis *basis, const struct gkyl_array *vmap, 
  double q2Dm, const struct gkyl_array *phi, const struct gkyl_array *phi_wall, 
  const struct gkyl_array *vcut_fact, struct gkyl_bc_sheath_gyrokinetic_kernels *kers, 
  struct gkyl_array *distf)
{
  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int pidx[GKYL_MAX_DIM];
  int vidx[2];
  int vcut_fact_idx[GKYL_MAX_CDIM];

  int pdim = skin_r.ndim;
  int vpar_dir = cdim;
  int uplo = skin_r.upper[vpar_dir]+skin_r.lower[vpar_dir];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange

    gkyl_sub_range_inv_idx(&skin_r, linc, pidx);
    // gkyl_sub_range_inv_idx(&vcut_r, linc, vcut_fact_idx);

    gkyl_copy_int_arr(pdim, pidx, fidx);
    fidx[vpar_dir] = uplo - pidx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[dir] = ghost_r.lower[dir];

    long skin_loc = gkyl_range_idx(&skin_r, pidx);
    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);
    double *out = (double*) gkyl_array_fetch(distf, ghost_loc);

    for (int d=cdim; d<pdim; d++) vidx[d-cdim] = pidx[d]; 
    long conf_loc = gkyl_range_idx(&conf_r, pidx);
    long vel_loc = gkyl_range_idx(&vel_r, vidx);

    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *vmap_p = (const double*) gkyl_array_cfetch(vmap, vel_loc);

    // Get vcut factor.
    for (int d=0; d<cdim-1; d++) vcut_fact_idx[d] = pidx[d]; // config space perp directions.
    vcut_fact_idx[cdim-1] = pidx[pdim-1]; // mu direction.
    long vcut_fact_loc = gkyl_range_idx(&vcut_r, vcut_fact_idx);
    const double *vcut_fact_p = (const double*) gkyl_array_cfetch(vcut_fact, vcut_fact_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[112];  // MF 2022/08/24: hardcoded to number of DG coeffs in 3x2v p2 for now.
    kers->reflectedf(vmap_p, q2Dm, phi_p, phi_wall_p, vcut_fact_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(dir, basis, cdim, out, fhat);
  }
}

__global__ static void
bc_gksheath_update_vcut_fact_surrogate_cu_ker(enum gkyl_edge_loc edge, const struct gkyl_range conf_r,
  const struct gkyl_range skin_r,
  const struct gkyl_range vcut_fact_r, const struct gkyl_range vel_r, const struct gkyl_basis *basis, 
  const struct gkyl_array *vmap, double q2Dm, const struct gkyl_array *phi, const struct gkyl_array *phi_wall, 
  const struct gkyl_array *dens, const struct gkyl_array *temp, const struct gkyl_array *bmag, 
  const struct gkyl_array *bimpact_angle, struct gkyl_bc_sheath_gyrokinetic_kernels *kers, 
  struct gkyl_array *vcut_fact)
{
  int cidx[GKYL_MAX_CDIM];
  int vidx[2];
  int vcut_fact_idx[GKYL_MAX_CDIM];

  int cdim = conf_r.ndim;
  int vcut_dim = vcut_fact_r.ndim;

  // the last configuration index is the boundary along z. Upper or lower
  cidx[cdim-1] = (edge == GKYL_LOWER_EDGE) ? conf_r.lower[cdim-1] : conf_r.upper[cdim-1];
  // For mu we take vpar = 0
  vidx[0] = skin_r.lower[cdim]; // vpar = 0

  // We loop over the vcut_fact range, x,y,mu
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < vcut_fact_r.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&vcut_fact_r, linc, vcut_fact_idx);
    
    for (int d=0; d<cdim-1; d++) cidx[d] = vcut_fact_idx[d]; // x,y indices.
    vidx[1] = vcut_fact_idx[vcut_dim-1]; // mu index.

    long conf_loc = gkyl_range_idx(&conf_r, cidx);
    long vel_loc = gkyl_range_idx(&vel_r, vidx);
    long vcut_fact_loc = gkyl_range_idx(&vcut_fact_r, vcut_fact_idx);

    const double *vmap_p = (const double*) gkyl_array_cfetch(vmap, vel_loc);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *dens_p = (const double*) gkyl_array_cfetch(dens, conf_loc);
    const double *temp_p = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    const double *bimpact_angle_p = (const double*) gkyl_array_cfetch(bimpact_angle, conf_loc);
    double *vcut_fact_p = (double*) gkyl_array_cfetch(vcut_fact, vcut_fact_loc);

    kers->surrogate(vmap_p, phi_p, phi_wall_p, dens_p, temp_p, q2Dm, bmag_p, bimpact_angle_p, vcut_fact_p);
  }
}

void
gkyl_bc_sheath_gyrokinetic_advance_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r)
{
  if (up->skin_r->volume > 0) {
    int nblocks = up->skin_r->nblocks, nthreads = up->skin_r->nthreads;

    gkyl_bc_sheath_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, *up->skin_r, *up->ghost_r,
      *conf_r, up->vcut_fact_local, up->vel_map->local_vel, up->basis, up->vel_map->vmap->on_dev, up->q2Dm, phi->on_dev, phi_wall->on_dev, up->vcut_fact->on_dev,
      up->kernels_cu, distf->on_dev);
  }
}

void 
bc_gksheath_update_vcut_fact_surrogate_cu(const struct gkyl_bc_sheath_gyrokinetic *up, 
  const struct gkyl_array *phi, const struct gkyl_array *phi_wall, const struct gkyl_array *dens, 
  const struct gkyl_array *temp, const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, 
  const struct gkyl_range *conf_r)
{
  if (up->vcut_fact_local.volume > 0) {
    int nblocks = up->vcut_fact_local.nblocks, nthreads = up->vcut_fact_local.nthreads;
    bc_gksheath_update_vcut_fact_surrogate_cu_ker<<<nblocks, nthreads>>>(up->edge, *conf_r, *up->skin_r,
      up->vcut_fact_local, up->vel_map->local_vel, up->basis, up->vel_map->vmap->on_dev, up->q2Dm, 
      phi->on_dev, phi_wall->on_dev, dens->on_dev, temp->on_dev, bmag->on_dev, bimpact_angle->on_dev, 
      up->kernels_cu, up->vcut_fact->on_dev);
  }
}