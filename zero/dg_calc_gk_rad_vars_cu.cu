/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_dg_calc_gk_rad_vars_nu_advance_cu_kernel(struct gkyl_dg_calc_gk_rad_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
  double xc[GKYL_MAX_DIM] = {0.0};
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);

    double* vnu_surf_d = (double*) gkyl_array_fetch(vnu_surf, loc_phase);
    double* vnu_d = (double*) gkyl_array_fetch(vnu, loc_phase);
    double* vsqnu_surf_d = (double*) gkyl_array_fetch(vsqnu_surf, loc_phase);  
    double* vsqnu_d = (double*) gkyl_array_fetch(vsqnu, loc_phase);   

    up->rad_nu_vpar(xc, up->phase_grid.dx, up->charge, up->mass, 
      up->a, up->alpha, up->beta, up->gamma, up->v0, 
      bmag_d, vnu_surf_d, vnu_d);
    up->rad_nu_mu(xc, up->phase_grid.dx, up->charge, up->mass, 
      up->a, up->alpha, up->beta, up->gamma, up->v0, 
      bmag_d, vsqnu_surf_d, vsqnu_d);
  }  
}

// Host-side wrapper for radiation drag coefficient calculation
void 
gkyl_dg_calc_gk_rad_vars_nu_advance_cu(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_gk_rad_vars_nu_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, 
    vnu_surf->on_dev, vnu->on_dev, vsqnu_surf->on_dev, vsqnu->on_dev);
}

__global__ void
gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu_kernel(struct gkyl_dg_calc_gk_rad_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array* vnu_surf, const struct gkyl_array* vnu, 
  const struct gkyl_array* vsqnu_surf, const struct gkyl_array* vsqnu, 
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu)
{
  double xc[GKYL_MAX_DIM] = {0.0};
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double* vnu_surf_d = (const double*) gkyl_array_cfetch(vnu_surf, loc_phase);
    const double* vnu_d = (const double*) gkyl_array_cfetch(vnu, loc_phase);
    const double* vsqnu_surf_d = (const double*) gkyl_array_cfetch(vsqnu_surf, loc_phase);  
    const double* vsqnu_d = (const double*) gkyl_array_cfetch(vsqnu, loc_phase);   

    const double *nI_d = (const double*) gkyl_array_cfetch(nI, loc_conf);

    double* nvnu_surf_d = (double*) gkyl_array_fetch(nvnu_surf, loc_phase);
    double* nvnu_d = (double*) gkyl_array_fetch(nvnu, loc_phase);
    double* nvsqnu_surf_d = (double*) gkyl_array_fetch(nvsqnu_surf, loc_phase);  
    double* nvsqnu_d = (double*) gkyl_array_fetch(nvsqnu, loc_phase);   

    up->rad_nI_nu(vnu_surf_d, vnu_d, vsqnu_surf_d, vsqnu_d, nI_d, 
      nvnu_surf_d, nvnu_d, nvsqnu_surf_d, nvsqnu_d);
  }  
}

// Host-side wrapper for sum_s n_{i_s} nu_s(v) calculation for a given input n_{i_s} and nu_s(v)
void 
gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* vnu_surf, const struct gkyl_array* vnu, 
  const struct gkyl_array* vsqnu_surf, const struct gkyl_array* vsqnu, 
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, *phase_range, 
    vnu_surf->on_dev, vnu->on_dev, vsqnu_surf->on_dev, vsqnu->on_dev, 
    nI->on_dev, 
    nvnu_surf->on_dev, nvnu->on_dev, nvsqnu_surf->on_dev, nvsqnu->on_dev);
}

// CUDA kernel to set device pointers to gyrokinetic radiation vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_gk_rad_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_gk_rad_vars *up, 
  int cdim, int vdim, int poly_order)
{
  up->rad_nu_vpar = choose_rad_gyrokinetic_nu_vpar_kern(cdim, vdim, poly_order);
  up->rad_nu_mu = choose_rad_gyrokinetic_nu_mu_kern(cdim, vdim, poly_order);
  up->rad_nI_nu = choose_rad_gyrokinetic_nI_nu_kern(cdim, vdim, poly_order);
}

gkyl_dg_calc_gk_rad_vars*
gkyl_dg_calc_gk_rad_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass, const struct gk_geometry *gk_geom, 
  double a, double alpha, double beta, double gamma, double v0)
{
  struct gkyl_dg_calc_gk_rad_vars *up = (struct gkyl_dg_calc_gk_rad_vars*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  // Fitting parameters for a given collision type
  up->a = a;
  up->alpha = alpha;
  up->beta = beta;
  up->gamma = gamma;
  up->v0 = v0;

  // acquire pointer to geometry object
  struct gk_geometry *geom = gkyl_gk_geometry_acquire(gk_geom);
  up->gk_geom = geom->on_dev; // this is so the memcpy below has geometry on_dev

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_gk_rad_vars *up_cu = (struct gkyl_dg_calc_gk_rad_vars*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_gk_rad_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_gk_rad_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // updater should store host pointers
  up->gk_geom = geom; 
  
  return up;
}
