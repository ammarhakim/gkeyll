/* -*- c++ -*- */

extern "C" {
#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wave_prop_priv.h>
}

__global__ static void
gkyl_wave_prop_waves_qfluct_cu_ker(gkyl_wv_eqn *eqn, 
  int ndim, int dir, double cflm, double dtdx, 
  enum gkyl_wv_flux_type ftype, struct gkyl_range update_range,
  struct gkyl_range fluctuation_range, const struct gkyl_wave_geom *wg,
  const struct gkyl_array *phi, const struct gkyl_array *qin,
  struct gkyl_array *waves, struct gkyl_array *speeds, struct gkyl_array *amdq,
  struct gkyl_array *apdq, struct gkyl_array *cfla, struct gkyl_array *is_cfl_violated)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long tid_x = threadIdx.x + blockIdx.x*blockDim.x;
    tid_x < fluctuation_range.volume; tid_x += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&fluctuation_range, tid_x, idxl);
    gkyl_copy_int_arr(ndim, idxl, idxr);

    // perform 1D sweeps handled by second dimension of thread grid
    idxl[dir] = idxl[dir]-1; 

    long lidx = gkyl_range_idx(&update_range, idxl);
    long ridx = gkyl_range_idx(&update_range, idxr);

    // Fetch the geometry, waves, speeds, and A^-/+ Delta Q fluctuations arrays
    // Every cell owns its lower interface values, so we fetch using ridx since idxr[dir] = tid_y is
    // the current cell. 
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxr);
    double *waves_d = (double*) gkyl_array_fetch(waves, ridx);
    double *speeds_d = (double*) gkyl_array_fetch(speeds, ridx);
    double *amdq_d = (double*) gkyl_array_fetch(amdq, ridx);
    double *apdq_d = (double*) gkyl_array_fetch(apdq, ridx);

    const double *qinl = (const double*) gkyl_array_cfetch(qin, lidx);
    const double *qinr = (const double*) gkyl_array_cfetch(qin, ridx);

    const double *phil = (const double*) gkyl_array_cfetch(phi, lidx);
    const double *phir = (const double*) gkyl_array_cfetch(phi, ridx);

    double my_max_speed = gkyl_wv_eqn_fused_rotate_waves_qfluct(eqn, ftype, 
      cg->tau1[dir], cg->tau2[dir], cg->norm[dir], cg->lenr[dir], 
      qinl, qinr, phil[0], phir[0], waves_d, speeds_d, amdq_d, apdq_d); 

    double *cfla_d = (double*) gkyl_array_fetch(cfla, ridx);  
    cfla_d[0] = dtdx/cg->kappa*my_max_speed; 
    // check time-step before any updates are performed
    if (cfla_d[0] > cflm) {
      double *is_cfl_violated_d = (double*) gkyl_array_fetch(is_cfl_violated, ridx);  
      is_cfl_violated_d[0] = 1.0;
    }      
  }
}

__global__ static void
gkyl_wave_prop_redo_waves_qfluct_cu_ker(gkyl_wv_eqn *eqn, 
  int ndim, int dir, int meqn, 
  struct gkyl_range update_range, struct gkyl_range fluctuation_range, 
  const struct gkyl_wave_geom *wg, const struct gkyl_array *phi,
  const struct gkyl_array *qin, struct gkyl_array *waves,
  struct gkyl_array *speeds, struct gkyl_array *amdq, struct gkyl_array *apdq, 
  struct gkyl_array *flux2, struct gkyl_array *redo_fluct, const struct gkyl_array *qout)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long tid_x = threadIdx.x + blockIdx.x*blockDim.x;
    tid_x < fluctuation_range.volume; tid_x += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&fluctuation_range, tid_x, idxl);
    gkyl_copy_int_arr(ndim, idxl, idxr);

    // perform 1D sweeps handled by second dimension of thread grid
    idxl[dir] = idxl[dir]-1; 

    long lidx = gkyl_range_idx(&update_range, idxl);
    long ridx = gkyl_range_idx(&update_range, idxr);
    const double *qinl = (const double*) gkyl_array_cfetch(qin, lidx);
    const double *qinr = (const double*) gkyl_array_cfetch(qin, ridx);

    const double *qoutl = (const double*) gkyl_array_cfetch(qout, lidx);
    const double *qoutr = (const double*) gkyl_array_cfetch(qout, ridx);

    const double *phil = (const double*) gkyl_array_cfetch(phi, lidx);
    const double *phir = (const double*) gkyl_array_cfetch(phi, ridx);
    // Zero out redo fluctuation flag so no time steps are redone by default. 
    double *redo_fluct_c = (double*) gkyl_array_fetch(redo_fluct, ridx);
    redo_fluct_c[0] = 0.0; 
    // Check the solution on both sides of the interface
    if (!gkyl_wv_eqn_fuse_check_inv(eqn, qoutl, qoutr)) {
      // Set redo_fluct flag to be true for re-doing update with re-computed fluctuations
      redo_fluct_c[0] = 1.0; 
      // Fetch the geometry, waves, speeds, and A^-/+ Delta Q fluctuations arrays
      // Every cell owns its lower interface values, so we fetch using ridx since idxr[dir] = linc2 is
      // the current cell. 
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxr);
      double *waves_d = (double*) gkyl_array_fetch(waves, ridx);
      double *speeds_d = (double*) gkyl_array_fetch(speeds, ridx);
      double *amdq_d = (double*) gkyl_array_fetch(amdq, ridx);
      double *apdq_d = (double*) gkyl_array_fetch(apdq, ridx);
      double my_max_speed = gkyl_wv_eqn_fused_rotate_waves_qfluct(eqn,
        GKYL_WV_LOW_ORDER_FLUX, cg->tau1[dir], cg->tau2[dir], cg->norm[dir],
        cg->lenr[dir], qinl, qinr, phil[0], phir[0], waves_d, speeds_d, amdq_d, apdq_d); 

      // Reset second order flux correction to zero if invariant domain violated 
      // so that no second order fluxes are used if we have to redo any fluctuation
      // computation.
      double *flux2_d = (double*) gkyl_array_fetch(flux2, ridx);
      for (int i=0; i<meqn; ++i) {
        flux2_d[i] = 0.0; 
      }
    }
  }
}

__global__ static void
gkyl_wave_prop_second_order_flux_cu_ker(enum gkyl_wave_limiter limiter, double dtdx, 
  int ndim, int dir, int meqn, int mwaves, 
  struct gkyl_range update_range, struct gkyl_range second_order_range, 
  const struct gkyl_wave_geom *wg, struct gkyl_array *waves,
  struct gkyl_array *waves_scaled, struct gkyl_array *speeds, struct gkyl_array *flux2)
{
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long tid_x = threadIdx.x + blockIdx.x*blockDim.x;
    tid_x < second_order_range.volume; tid_x += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&second_order_range, tid_x, idxc);
    gkyl_copy_int_arr(ndim, idxc, idxl);
    gkyl_copy_int_arr(ndim, idxl, idxr);

    // perform 1D sweeps handled by second dimension of thread grid
    idxl[dir] = idxl[dir]-1; 
    idxr[dir] = idxr[dir]+1; 
    long linl = gkyl_range_idx(&update_range, idxl); 
    long linc = gkyl_range_idx(&update_range, idxc); 
    long linr = gkyl_range_idx(&update_range, idxr);
    const double *s = (const double*) gkyl_array_cfetch(speeds, linc);
    double *wc = (double*) gkyl_array_cfetch(waves, linc);
    double *ws = (double*) gkyl_array_cfetch(waves_scaled, linc);
    const double *wl = (const double*) gkyl_array_cfetch(waves, linl);
    const double *wr = (const double*) gkyl_array_cfetch(waves, linr);

    for (int mw=0; mw<mwaves; ++mw) {
      double dotl = wave_dot_prod(meqn, &wl[mw*meqn], &wc[mw*meqn]);
      double wnorm2 = wave_dot_prod(meqn, &wc[mw*meqn], &wc[mw*meqn]);
      double dotr = wave_dot_prod(meqn, &wc[mw*meqn], &wr[mw*meqn]);
      if (wnorm2 > 0) {
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, limiter);
        wave_rescale(meqn, theta, &wc[mw*meqn], &ws[mw*meqn]);
      }
    }
    // Each cell owns its flux at the lower interface, so fetch the 
    // lower and center geometry to get kappa on either side of the interface
    const struct gkyl_wave_cell_geom *cgl = gkyl_wave_geom_get(wg, idxl);
    const struct gkyl_wave_cell_geom *cgc = gkyl_wave_geom_get(wg, idxc);
    double kappal = cgl->kappa;
    double kappac = cgc->kappa;
    double *flux2_d = (double*) gkyl_array_fetch(flux2, linc);
    for (int mw=0; mw<mwaves; ++mw) {
      calc_second_order_qflux(meqn, dtdx/(0.5*(kappal+kappac)), s[mw], &ws[mw*meqn], flux2_d);
    }
  }
}

__global__ static void
gkyl_wave_prop_update_state_cu_kern(double dtdx, 
  int ndim, int dir, int meqn, struct gkyl_range update_range, 
  const struct gkyl_wave_geom *wg, const struct gkyl_array *redo_fluct, const struct gkyl_array *qin, 
  const struct gkyl_array *amdq, const struct gkyl_array *apdq, const struct gkyl_array *flux2, 
  struct gkyl_array *qout)
{
  int idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long tid_x = threadIdx.x + blockIdx.x*blockDim.x;
    tid_x < update_range.volume; tid_x += blockDim.x*gridDim.x) { 

    gkyl_sub_range_inv_idx(&update_range, tid_x, idxc);
    gkyl_copy_int_arr(ndim, idxc, idxr);

    idxr[dir] = idxr[dir]+1; 
    long linc = gkyl_range_idx(&update_range, idxc); 
    long linr = gkyl_range_idx(&update_range, idxr);
    const double *redo_fluct_c = (const double*) gkyl_array_cfetch(redo_fluct, linc);
    if (redo_fluct_c[0] > 0.0) {
      // Each cell owns its flux at the lower interface, so fetch the 
      // lower geometry to get kappa in the cell and the right-going fluctuations 
      // and second order fluxes on the lower interface and left-going fluctuations
      // and second order fluxes on the upper interface
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxc);
      const double *apdq_l = (const double*) gkyl_array_cfetch(apdq, linc);
      const double *amdq_r = (const double*) gkyl_array_cfetch(amdq, linr);
      const double *flux2_l = (const double*) gkyl_array_cfetch(flux2, linc);
      const double *flux2_r = (const double*) gkyl_array_cfetch(flux2, linr);
      const double *qin_c = (const double*) gkyl_array_cfetch(qin, linc);
      double *qout_c = (double*) gkyl_array_fetch(qout, linc);
      for (int i=0; i<meqn; ++i) {
        qout_c[i] = qin_c[i] 
          - dtdx/cg->kappa*(apdq_l[i] + amdq_r[i] + flux2_r[i] - flux2_l[i]);
      }
    }
  }
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance_cu(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout)
{
  wv->n_calls += 1;
  
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations;
  //  when forced to use Lax fluxes, we only have a single wave
  int mwaves = wv->force_low_order_flux ? 2 :  wv->equation->num_waves;

  double cfl = wv->cfl, cflm = 1.1*cfl;
  gkyl_array_clear(wv->cfla, 0.0); 
  gkyl_array_clear(wv->is_cfl_violated, 0.0);
  // Small arrays for reducing actual cfl value over the whole domain
  // and whether the cfl was violated anywhere in the domain. 
  double red_cfla[1];
  double red_is_cfl_violated[1];
  double red_cfl[2];
  double red_cfl_global[2];

  // Assume GPU simulation is genuinely dimensionally split and num_update_dirs = 1
  int dir = wv->update_dirs[0];
  double dtdx = dt/wv->grid.dx[dir];
  // Determine number of threads based on extended one-dimensional range 
  // and volume of perpendicular range    
  // upper/lower bounds in direction 'd'. These are edge indices
  int loidx = update_range->lower[dir]-1;
  int upidx = update_range->upper[dir]+2;

  struct gkyl_range perp_range;
  gkyl_range_shorten_from_above(&perp_range, update_range, dir, 1);
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks_e = (perp_range.volume*(upidx - loidx))/nthreads + 1;
  
  // Extend the local range for indexing fluctuations
  struct gkyl_range fluctuation_range;
  int extend_lo[] = {0, 0, 0};
  int extend_up[] = {0, 0, 0};
  extend_lo[dir] = 1;
  extend_up[dir] = 2;
  gkyl_range_extend(&fluctuation_range, update_range, extend_lo, extend_up);
  // Extend the local range for indexing second order fluxes
  struct gkyl_range second_order_range;
  extend_lo[dir] = 0;
  extend_up[dir] = 1;
  gkyl_range_extend(&second_order_range, update_range, extend_lo, extend_up);

  // Copy previous time step solution 
  gkyl_array_set(qout, 1.0, qin); 
  // Set the redo_fluct array so in the first sweep, we compute fluxes at every interface
  gkyl_array_clear(wv->redo_fluct, 1.0); 

  enum gkyl_wv_flux_type ftype = wv->force_low_order_flux ?
    GKYL_WV_LOW_ORDER_FLUX : GKYL_WV_HIGH_ORDER_FLUX;

  gkyl_wave_prop_waves_qfluct_cu_ker<<<nblocks_e, nthreads>>>(wv->equation->on_dev, 
    ndim, dir, cflm, dtdx, ftype, *update_range, fluctuation_range, 
    wv->geom->on_dev, phi->on_dev, qin->on_dev, wv->waves->on_dev,
    wv->speeds->on_dev, wv->amdq->on_dev, wv->apdq->on_dev, wv->cfla->on_dev,
    wv->is_cfl_violated->on_dev);
  
  // To avoid race conditions on the wave limiting, split the second order flux
  // computation into a separate kernel after waves and fluctuations computed
  int nblocks_f = (perp_range.volume*(upidx - loidx - 2))/nthreads + 1;
  gkyl_array_clear(wv->flux2, 0.0);
  gkyl_wave_prop_second_order_flux_cu_ker<<<nblocks_f, nthreads>>>(wv->limiter, dtdx, 
    ndim, dir, meqn, mwaves, *update_range, second_order_range, 
    wv->geom->on_dev, wv->waves->on_dev, wv->waves_scaled->on_dev,
    wv->speeds->on_dev, wv->flux2->on_dev);

  // Before updating the state, determine if the CFL was violated anywhere in the domain
  gkyl_array_reduce_range(wv->cfla_ptr, wv->cfla, GKYL_MAX, update_range);  
  gkyl_array_reduce_range(wv->is_cfl_violated_ptr, wv->is_cfl_violated, GKYL_MAX, update_range);
  gkyl_cu_memcpy(red_cfla, wv->cfla_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_memcpy(red_is_cfl_violated, wv->is_cfl_violated_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  red_cfl[0] = red_cfla[0]; 
  red_cfl[1] = red_is_cfl_violated[0]; 
  gkyl_comm_allreduce_host(wv->comm, GKYL_DOUBLE, GKYL_MAX, 2, &red_cfl, &red_cfl_global);
  double dt_suggested = dt*cfl/fmax(red_cfl_global[0], DBL_MIN);
  if (red_cfl_global[1] > 0.0) {
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_wave_prop_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
      .max_speed = 0.0, // Not currently used by GPUs
    };    
  }

  // Update the solution with both the first order update and the second order corrections. 
  gkyl_wave_prop_update_state_cu_kern<<<update_range->nblocks, update_range->nthreads>>>(dtdx, 
    ndim, dir, meqn, *update_range, 
    wv->geom->on_dev, wv->redo_fluct->on_dev, qin->on_dev, 
    wv->amdq->on_dev, wv->apdq->on_dev,  wv->flux2->on_dev, 
    qout->on_dev); 

  // Determine if we need to redo any flux computations 
  if (wv->check_inv_domain) {
    gkyl_wave_prop_redo_waves_qfluct_cu_ker<<<update_range->nblocks, update_range->nthreads>>>(wv->equation->on_dev, 
      ndim, dir, meqn, *update_range, fluctuation_range, 
      wv->geom->on_dev, phi->on_dev, qin->on_dev, wv->waves->on_dev,
      wv->speeds->on_dev, wv->amdq->on_dev, wv->apdq->on_dev, wv->flux2->on_dev,
      wv->redo_fluct->on_dev, qout->on_dev); 
    
    // Update state with re-computed first order fluxes. 
    // Note that the second order fluxes are zeroed out at re-done interfaces
    // so no second order corrections are included when re-doing the update. 
    gkyl_wave_prop_update_state_cu_kern<<<update_range->nblocks, update_range->nthreads>>>(dtdx, 
      ndim, dir, meqn, *update_range, 
      wv->geom->on_dev, wv->redo_fluct->on_dev, qin->on_dev, 
      wv->amdq->on_dev, wv->apdq->on_dev,  wv->flux2->on_dev, 
      qout->on_dev); 
  }

  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).
  return (struct gkyl_wave_prop_status) {
    .success = red_cfl_global[1] > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
    .max_speed = 0.0, // Not currently used by GPUs
  }; 
}

__global__ static void
gkyl_wave_prop_max_dt_cu_ker(gkyl_wv_eqn *eqn, double dx, double cfl,  
  struct gkyl_range update_range, const struct gkyl_wave_geom *wg, const struct gkyl_array *qin, 
  struct gkyl_array *cfla)
{
  int idxc[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
    tid < update_range.volume; tid += blockDim.x*gridDim.x) { 

    gkyl_sub_range_inv_idx(&update_range, tid, idxc);
    long linc = gkyl_range_idx(&update_range, idxc); 
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxc);
    const double *qc = (const double*) gkyl_array_cfetch(qin, linc);
    double *cfla_d = (double*) gkyl_array_fetch(cfla, linc);  
    double my_max_speed = gkyl_wv_eqn_max_speed(eqn, qc);
    cfla_d[0] = (cfl*dx*cg->kappa)/my_max_speed; 
  }  
}

// max dt method
double 
gkyl_wave_prop_max_dt_cu(const gkyl_wave_prop *wv,
  const struct gkyl_range *update_range, const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  double cfl = wv->cfl;
  gkyl_array_clear(wv->cfla, 0.0); 
  // Small arrays for reducing actual cfl value over the whole domain. 
  double red_cfla[1];
  double red_cfla_global[1];

  // Loop over update directions
  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];
    double dx = wv->grid.dx[dir];
    gkyl_wave_prop_max_dt_cu_ker<<<update_range->nblocks, update_range->nthreads>>>(wv->equation->on_dev, 
      dx, cfl, *update_range, wv->geom->on_dev, qin->on_dev, wv->cfla->on_dev);
    // Reduce over the domain to find maximum time step
    gkyl_array_reduce_range(wv->cfla_ptr, wv->cfla, GKYL_MIN, update_range);  
    gkyl_cu_memcpy(red_cfla, wv->cfla_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_comm_allreduce_host(wv->comm, GKYL_DOUBLE, GKYL_MIN, 1, &red_cfla, &red_cfla_global);
    max_dt = fmin(max_dt, red_cfla_global[0]);
  }  
  return max_dt;      
}
