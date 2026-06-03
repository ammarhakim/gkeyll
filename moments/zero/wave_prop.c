#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wave_prop_priv.h>

#include <gkyl_level_set.h>

gkyl_wave_prop*
gkyl_wave_prop_new(const struct gkyl_wave_prop_inp *winp)
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(*up));
  up->use_gpu = winp->use_gpu;

  up->grid = *(winp->grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = winp->num_up_dirs;
  for (int i=0; i<winp->num_up_dirs; ++i)
    up->update_dirs[i] = winp->update_dirs[i];

  up->limiter = winp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp->limiter;
  up->cfl = winp->cfl;
  up->equation = gkyl_wv_eqn_acquire(winp->equation);

  if (winp->comm)
    up->comm = gkyl_comm_acquire(winp->comm);
  else
    up->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  up->force_low_order_flux = winp->force_low_order_flux;
  up->check_inv_domain = winp->check_inv_domain;

  up->split_type = winp->split_type;

  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  long tot_vol = winp->update_range_ext->volume;

  // allocate memory to store 1D slices of waves, speeds and
  // second-order correction flux
  int meqn = winp->equation->num_equations, mwaves = winp->equation->num_waves;
  if (up->use_gpu) {
    up->waves = gkyl_array_cu_dev_new(GKYL_DOUBLE, meqn*mwaves, tot_vol);
    up->apdq = gkyl_array_cu_dev_new(GKYL_DOUBLE, meqn, tot_vol);
    up->amdq = gkyl_array_cu_dev_new(GKYL_DOUBLE, meqn, tot_vol);
    up->speeds = gkyl_array_cu_dev_new(GKYL_DOUBLE, mwaves, tot_vol);
    up->flux2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, meqn, tot_vol);
    up->redo_fluct = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_vol);
    up->cfla = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_vol);
    up->is_cfl_violated = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_vol);
    up->cfla_ptr = gkyl_cu_malloc(sizeof(double));
    up->is_cfl_violated_ptr = gkyl_cu_malloc(sizeof(double));
  }
  else {
    up->waves = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, tot_vol);
    up->waves_scaled = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, tot_vol);
    up->apdq = gkyl_array_new(GKYL_DOUBLE, meqn, tot_vol);
    up->amdq = gkyl_array_new(GKYL_DOUBLE, meqn, tot_vol);
    up->speeds = gkyl_array_new(GKYL_DOUBLE, mwaves, tot_vol);
    up->flux2 = gkyl_array_new(GKYL_DOUBLE, meqn, tot_vol);
    up->redo_fluct = gkyl_array_new(GKYL_DOUBLE, 1, tot_vol);
  }

  up->geom = gkyl_wave_geom_acquire(winp->geom);

  up->n_calls = up->n_bad_advance_calls = 0;
  up->n_bad_cells = up->n_max_bad_cells = 0;

  return up;
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout)
{
#ifdef GKYL_HAVE_CUDA  
  if (wv->use_gpu) {
    return gkyl_wave_prop_advance_cu(wv, tm, dt, update_range, qin, qout); 
  }
#endif
  wv->n_calls += 1;
  
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations;
  //  when forced to use Lax fluxes, we only have a single wave
  int mwaves = wv->force_low_order_flux ? 2 :  wv->equation->num_waves;

  double cfla[1] = { 0.0 }, cfl = wv->cfl, cflm = 1.1*cfl;
  double is_cfl_violated[1] = { 0.0 }; // delibrately a double
  
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  double max_speed[1] = { 0.0 };

  // Assume simulation is genuinely dimensionally split and num_update_dirs = 1
  int dir = wv->update_dirs[0];
  double dtdx = dt/wv->grid.dx[dir];

  struct gkyl_range perp_range;
  gkyl_range_shorten_from_above(&perp_range, update_range, dir, 1);

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

  // calculate waves and fluctuations
  wave_prop_waves_qfluct(wv->equation, ndim, dir, cflm, dtdx, ftype, update_range,
    fluctuation_range, wv->geom, phi, qin, wv->waves, wv->speeds, wv->amdq, wv->apdq,
    max_speed, cfla, is_cfl_violated);

  // calculate second order flux
  gkyl_array_clear(wv->flux2, 0.0);
  wave_prop_second_order_flux(wv->limiter, dtdx, ndim, dir, meqn, mwaves,
    update_range, second_order_range, wv->geom, wv->waves, wv->waves_scaled,
    wv->speeds, wv->flux2);
  
  // compute actual CFL, status & max-speed across all domains
  double red_vars[3] = { cfla[0], is_cfl_violated[0], max_speed[0] };
  double red_vars_global[3] = { 0.0, 0.0, 0.0 };
  gkyl_comm_allreduce(wv->comm, GKYL_DOUBLE, GKYL_MAX, 3, &red_vars, &red_vars_global);

  cfla[0] = red_vars_global[0];
  is_cfl_violated[0] = red_vars_global[1];
  max_speed[0] = red_vars_global[2];

  double dt_suggested = dt*cfl/fmax(cfla[0], DBL_MIN);
  if (is_cfl_violated[0] > 0.0) {
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_wave_prop_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
      .max_speed = max_speed[0],
    };
  }

  // first order update
  wave_prop_first_order(wv->equation, ndim, dir, dtdx, update_range, wv->geom,
    wv->amdq, wv->apdq, wv->flux2, qout);
  
  // Determine if we need to redo any flux computations
  if (wv->check_inv_domain) {
    wave_prop_check_inv(wv->equation, ndim, dir, cflm, dtdx, ftype, update_range,
      fluctuation_range, wv->geom, phi, qin, wv->waves, wv->speeds, wv->amdq, wv->apdq,
      wv->redo_fluct, max_speed, cfla, is_cfl_violated, qout);

    wave_prop_first_order_redo(wv->equation, ndim, dir, dtdx, update_range, wv->geom,
      wv->redo_fluct, wv->amdq, wv->apdq, qout);
  }
  return (struct gkyl_wave_prop_status) {
    .success = is_cfl_violated[0] > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
    .max_speed = max_speed[0],
  };
}

double
gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
#ifdef GKYL_HAVE_CUDA  
  if (wv->use_gpu) {
    return gkyl_wave_prop_max_dt_cu(wv, update_range, qin); 
  }
#endif
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];
      double dx = wv->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(wv->equation, q);
      max_dt = fmin(max_dt, wv->cfl*dx/maxs);
    }
    
  }

  return max_dt;
}

struct gkyl_wave_prop_stats
gkyl_wave_prop_stats(const gkyl_wave_prop *wv)
{
  return (struct gkyl_wave_prop_stats) {
    .n_calls = wv->n_calls,
    .n_bad_advance_calls = wv->n_bad_advance_calls,
    .n_bad_cells = wv->n_bad_cells,
    .n_max_bad_cells = wv->n_max_bad_cells
  };
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_array_release(up->waves);
  gkyl_array_release(up->apdq);
  gkyl_array_release(up->amdq);
  gkyl_array_release(up->speeds);
  gkyl_array_release(up->flux2);
  gkyl_array_release(up->redo_fluct);
  gkyl_comm_release(up->comm);
  if (up->use_gpu) {
    gkyl_array_release(up->cfla);
    gkyl_array_release(up->is_cfl_violated);   
    gkyl_cu_free(up->cfla_ptr); 
    gkyl_cu_free(up->is_cfl_violated_ptr); 
  }
  
  gkyl_wave_geom_release(up->geom);
  
  gkyl_free(up);
}
