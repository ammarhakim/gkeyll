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
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type b_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->vcut_calc = ser_sheath_vcut_calc_list[edge].dim_list[dim-2].kernels[poly_order-1];
      kers->build_input = ser_sheath_input_list[edge].dim_list[dim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_bc_gksheath_choose_surrogate_kernels_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  gkyl_bc_gksheath_surr_set_cu_ker_ptrs<<<1,1>>>(basis, edge, kers);
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
    int vcut_fact_dim = vcut_r.ndim;
    for (int d=0; d<vcut_fact_dim-1; d++) vcut_fact_idx[d] = pidx[d]; // config space perp directions.
    vcut_fact_idx[vcut_fact_dim-1] = pidx[pdim-1]; // mu direction.
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

// Phase 1: fill kann_inp->data for all perpendicular config cells.
// Each thread handles one perp cell; build_input writes n_nodes_per_cell nodes
// with stride dim_in, matching gkyl_kn_vec node-major layout.
__global__ static void
bc_gksheath_fill_kann_inp_cu_ker(enum gkyl_edge_loc edge,
  const struct gkyl_range conf_r, const struct gkyl_range perp_r,
  const struct gkyl_array *phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
  float *kann_inp_data, int n_nodes_per_cell, int dim_in,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int perp_idx[GKYL_MAX_CDIM];
  int cidx[GKYL_MAX_CDIM];
  int cdim = conf_r.ndim;

  cidx[cdim-1] = (edge == GKYL_LOWER_EDGE) ? conf_r.lower[cdim-1] : conf_r.upper[cdim-1];

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
       linc < perp_r.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&perp_r, linc, perp_idx);
    for (int d = 0; d < cdim-1; d++) cidx[d] = perp_idx[d];

    long conf_loc = gkyl_range_idx(&conf_r, cidx);
    // 0-based sequential perp cell index (gkyl_range_idx gives 0 at lower bound).
    long perp_k = gkyl_range_idx(&perp_r, perp_idx);

    const double *phi_p          = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p     = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);
    const double *dens_p         = (const double*) gkyl_array_cfetch(dens, conf_loc);
    const double *temp_p         = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p         = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    const double *bimpact_angle_p = (const double*) gkyl_array_cfetch(bimpact_angle, conf_loc);

    // Cell perp_k occupies nodes [perp_k*n_nodes_per_cell, (perp_k+1)*n_nodes_per_cell)
    // in kann_inp_data, each node taking dim_in floats.
    int kann_idx_shift = perp_k * n_nodes_per_cell * dim_in;
    kers->build_input(phi_p, phi_wall_p, dens_p, temp_p, bmag_p, bimpact_angle_p,
      dim_in, kann_inp_data + kann_idx_shift);
  }
}

// Phase 2 (after gkyl_kann_net_apply): read kann_out->data and project vcut_fact onto DG.
// Each thread handles one (perp, mu) cell.
__global__ static void
bc_gksheath_fill_vcut_fact_cu_ker(enum gkyl_edge_loc edge,
  const struct gkyl_range conf_r, const struct gkyl_range perp_r,
  const struct gkyl_range vcut_fact_r, const struct gkyl_range vel_r,
  const struct gkyl_array *vmap, int vpar_lower,
  const float *kann_out_data, int n_nodes_per_cell, int dim_out,
  const struct gkyl_array *temp, const struct gkyl_array *bmag,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers,
  struct gkyl_array *vcut_fact)
{
  int vcut_idx[GKYL_MAX_CDIM];
  int cidx[GKYL_MAX_CDIM];
  int vidx[2];

  int cdim = conf_r.ndim;
  int vcut_dim = vcut_fact_r.ndim;

  cidx[cdim-1] = (edge == GKYL_LOWER_EDGE) ? conf_r.lower[cdim-1] : conf_r.upper[cdim-1];
  vidx[0] = vpar_lower;

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
       linc < vcut_fact_r.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&vcut_fact_r, linc, vcut_idx);
    for (int d = 0; d < cdim-1; d++) cidx[d] = vcut_idx[d]; // perp config indices.
    vidx[1] = vcut_idx[vcut_dim-1]; // mu index.

    long conf_loc     = gkyl_range_idx(&conf_r, cidx);
    long vel_loc      = gkyl_range_idx(&vel_r, vidx);
    long vcut_fact_loc = gkyl_range_idx(&vcut_fact_r, vcut_idx);
    // cidx[0..cdim-2] are the perp indices; gkyl_range_idx on perp_r only reads those.
    long perp_k = gkyl_range_idx(&perp_r, cidx);

    const double *vmap_p      = (const double*) gkyl_array_cfetch(vmap, vel_loc);
    const double *temp_p      = (const double*) gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p      = (const double*) gkyl_array_cfetch(bmag, conf_loc);
    double *vcut_fact_p       = (double*) gkyl_array_fetch(vcut_fact, vcut_fact_loc);

    int kann_idx_shift = perp_k * n_nodes_per_cell * dim_out;
    kers->vcut_calc(vmap_p, kann_out_data + kann_idx_shift, dim_out, temp_p, bmag_p, vcut_fact_p);
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
  // Build a perp range (cdim-1 dimensional) from the perp dims of skin_r.
  // gkyl_range_idx on this range gives a sequential 0-based index for kann_inp/kann_out.
  struct gkyl_range perp_r;
  int perp_lower[GKYL_MAX_CDIM], perp_upper[GKYL_MAX_CDIM];
  for (int d = 0; d < up->cdim-1; d++) {
    perp_lower[d] = up->skin_r->lower[d];
    perp_upper[d] = up->skin_r->upper[d];
  }
  gkyl_range_init(&perp_r, up->cdim-1, perp_lower, perp_upper);

  int n_nodes_per_cell = 1 << (up->cdim-1);
  int dim_in  = up->kann_inp->N;
  int dim_out = up->kann_out->N;
  int vpar_lower = up->skin_r->lower[up->cdim];

  // Phase 1: fill kann_inp on GPU.
  if (perp_r.volume > 0) {
    int nblocks = perp_r.nblocks, nthreads = perp_r.nthreads;
    bc_gksheath_fill_kann_inp_cu_ker<<<nblocks, nthreads>>>(up->edge, *conf_r, perp_r,
      phi->on_dev, phi_wall->on_dev, dens->on_dev, temp->on_dev, bmag->on_dev, bimpact_angle->on_dev,
      up->kann_inp->data, n_nodes_per_cell, dim_in, up->kernels_cu);
  }

  // Phase 2: run NN inference (CUDA kernel dispatched inside gkyl_kann_net_apply).
  // CUDA launches are serialized in the default stream, so this runs after phase 1.
  gkyl_kann_net_apply(up->kann_net, up->kann_inp, up->kann_out);

  // Phase 3: project NN output onto DG vcut_fact.
  if (up->vcut_fact_local.volume > 0) {
    int nblocks = up->vcut_fact_local.nblocks, nthreads = up->vcut_fact_local.nthreads;
    bc_gksheath_fill_vcut_fact_cu_ker<<<nblocks, nthreads>>>(up->edge, *conf_r, perp_r,
      up->vcut_fact_local, up->vel_map->local_vel, up->vel_map->vmap->on_dev, vpar_lower,
      up->kann_out->data, n_nodes_per_cell, dim_out,
      temp->on_dev, bmag->on_dev, up->kernels_cu, up->vcut_fact->on_dev);
  }
}