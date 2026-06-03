#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>

struct gkyl_wave_prop {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_wave_limiter limiter; // limiter to use
  double cfl; // CFL number
  const struct gkyl_wv_eqn *equation; // equation object

  bool force_low_order_flux; // only use Lax flux
  bool check_inv_domain; // flag to indicate if invariant domains are checked

  enum gkyl_wave_split_type split_type; // type of splitting to use

  struct gkyl_wave_geom *geom; // geometry object
  struct gkyl_comm *comm; // communcator
  
  // data for 1D slice update
  struct gkyl_array *waves, *waves_scaled, *apdq, *amdq, *speeds, *flux2;
  // flags to indicate if fluctuations should be recomputed
  struct gkyl_array *redo_fluct;
  // CFL arrays
  struct gkyl_array *cfla, *is_cfl_violated; 
  double *cfla_ptr, *is_cfl_violated_ptr; 

  bool use_gpu; // are we using GPUs?

  // some stats
  long n_calls; // number of calls to updater
  long n_bad_advance_calls; // number of calls in which positivity had to be fixed
  long n_bad_cells; // number  of cells fixed
  long n_max_bad_cells; // maximum number of cells fixed in a call
};

static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
      theta = 1.0;
      break;
    
    // ** Fully formally-verified implementation of the minmod flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_minmod_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_minmod_tvd.rkt **
    case GKYL_MIN_MOD:
      theta = fmax(0.0, fmin(1.0, r));
      break;

    // ** Partially formally-verified implementation of the superbee flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): NOT PROVEN **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_superbee_tvd.rkt **
    case GKYL_SUPERBEE:
      theta = fmax3(0.0, fmin((2.0 * r), 1.0), fmin(r, 2.0));
      break;

    // ** Partially formally-verified implementation of the van Leer flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_van_leer_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): NOT PROVEN **
    case GKYL_VAN_LEER:
      theta = ((r + fabs(r)) / (1.0 + fabs(r)));
      break;

    // ** Fully formally-verified implementation of the monotonized-centered flux limiter **
    // ** Proof of symmetry (equivalent action on forward and backward gradients): ../proofs/finite_volume/proof_limiter_monotonized_centered_symmetry.rkt **
    // ** Proof of second-order TVD (total variation diminishing): ../proofs/finite_volume/proof_limiter_monotonized_centered_tvd.rkt **
    case GKYL_MONOTONIZED_CENTERED:
      theta = fmax(0.0, fmin3((2.0 * r), ((1.0 + r) / 2.0), 2.0));
      break;

    case GKYL_BEAM_WARMING:
      theta = r;
      break;

    case GKYL_ZERO:
      theta = 0;
      break;
  }
  return theta;
}

// some helper functions
static inline void
wave_prop_waves_qfluct(const gkyl_wv_eqn *eqn, int ndim, int dir, double cflm,
  double dtdx, enum gkyl_wv_flux_type ftype, const struct gkyl_range *update_range,
  struct gkyl_range fluctuation_range, const struct gkyl_wave_geom *wg,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *waves,
  struct gkyl_array *speeds, struct gkyl_array *amdq, struct gkyl_array *apdq,
  double *max_speed, double *cfla, double *is_cfl_violated)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &fluctuation_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxl);
    gkyl_copy_int_arr(ndim, iter.idx, idxr);

    idxl[dir] = idxl[dir]-1;

    long lidx = gkyl_range_idx(update_range, idxl);
    long ridx = gkyl_range_idx(update_range, idxr);

    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxr);
    double *waves_local = gkyl_array_fetch(waves, ridx);
    double *amdq_local = gkyl_array_fetch(amdq, ridx);
    double *apdq_local = gkyl_array_fetch(apdq, ridx);
    double *s = gkyl_array_fetch(speeds, ridx);
    
    const double *qinl = gkyl_array_cfetch(qin, lidx);
    const double *qinr = gkyl_array_cfetch(qin, ridx);

    const double *phil = gkyl_array_cfetch(phi, lidx);
    const double *phir = gkyl_array_cfetch(phi, ridx);

    double my_max_speed = gkyl_wv_eqn_fused_rotate_waves_qfluct(eqn,
      ftype, cg->tau1[dir], cg->tau2[dir], cg->norm[dir], cg->lenr[dir],
      qinl, qinr, phil[0], phir[0], waves_local, s, amdq_local, apdq_local);
    max_speed[0] = max_speed[0] > my_max_speed ? max_speed[0] : my_max_speed;
    
    // check time-step before any updates are performed
    cfla[0] = fmax(cfla[0], dtdx/cg->kappa*my_max_speed);
    if (cfla[0] > cflm)
      is_cfl_violated[0] = 1.0;
  }
}

static inline void
copy_wv_vec(int n, double * GKYL_RESTRICT out, const double * GKYL_RESTRICT inp)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT q, const double * GKYL_RESTRICT amdq_r, const double * GKYL_RESTRICT apdq_l, const double *fr, const double *fl)
{
  for (int i=0; i<meqn; ++i)
    q[i] = q[i] - dtdx*(apdq_l[i] + amdq_r[i] + fr[i] - fl[i]);
}

static inline void
wave_prop_first_order(const gkyl_wv_eqn *eqn, int ndim, int dir, double dtdx,
  const struct gkyl_range *update_range, const struct gkyl_wave_geom *wg,
  struct gkyl_array *amdq, struct gkyl_array *apdq, struct gkyl_array *flux2,
  struct gkyl_array *qout)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  int meqn = eqn->num_equations;

  // compute first-order update in each cell
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxl);
    gkyl_copy_int_arr(ndim, iter.idx, idxr);

    idxr[dir] = idxr[dir]+1;
    long lidx = gkyl_range_idx(update_range, idxl);
    long ridx = gkyl_range_idx(update_range, idxr);

    double *apdq_local = gkyl_array_fetch(apdq, lidx);
    double *amdq_local = gkyl_array_fetch(amdq, ridx);

    double *flux2l = gkyl_array_fetch(flux2, lidx);
    double *flux2r = gkyl_array_fetch(flux2, ridx);

    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxl);
    calc_first_order_update(meqn, dtdx/cg->kappa, gkyl_array_fetch(qout, lidx),
      amdq_local, apdq_local, flux2r, flux2l);
  }
}

static inline void
wave_prop_check_inv(const gkyl_wv_eqn *eqn, int ndim, int dir, double cflm,
  double dtdx, enum gkyl_wv_flux_type ftype, const struct gkyl_range *update_range,
  struct gkyl_range fluctuation_range, const struct gkyl_wave_geom *wg,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *waves,
  struct gkyl_array *speeds, struct gkyl_array *amdq, struct gkyl_array *apdq,
  struct gkyl_array *redo_fluct, double *max_speed, double *cfla,
  double *is_cfl_violated, const struct gkyl_array *qout)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &fluctuation_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxl);
    gkyl_copy_int_arr(ndim, iter.idx, idxr);

    idxl[dir] = idxl[dir]-1;

    long lidx = gkyl_range_idx(update_range, idxl);
    long ridx = gkyl_range_idx(update_range, idxr);

    double *redo_fluct_local = gkyl_array_fetch(redo_fluct, ridx);

    const double *qinl = gkyl_array_cfetch(qin, lidx);
    const double *qinr = gkyl_array_cfetch(qin, ridx);

    const double *qoutl = gkyl_array_cfetch(qout, lidx);
    const double *qoutr = gkyl_array_cfetch(qout, ridx);

    const double *phil = gkyl_array_cfetch(phi, lidx);
    const double *phir = gkyl_array_cfetch(phi, ridx);
  
    redo_fluct_local[0] = 0.0;

    if (!gkyl_wv_eqn_fuse_check_inv(eqn, qoutl, qoutr)) {
      redo_fluct_local[0] = 1.0;
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxr);
      double *waves_local = gkyl_array_fetch(waves, ridx);
      double *amdq_local = gkyl_array_fetch(amdq, ridx);
      double *apdq_local = gkyl_array_fetch(apdq, ridx);
      double *s = gkyl_array_fetch(speeds, ridx);
        
      double my_max_speed = gkyl_wv_eqn_fused_rotate_waves_qfluct(eqn,
        GKYL_WV_LOW_ORDER_FLUX, cg->tau1[dir], cg->tau2[dir], cg->norm[dir],
        cg->lenr[dir], qinl, qinr, phil[0], phir[0], waves_local, s, amdq_local,
        apdq_local);
      max_speed[0] = max_speed[0] > my_max_speed ? max_speed[0] : my_max_speed;
    }
  }
}

static inline void
wave_prop_first_order_redo(const gkyl_wv_eqn *eqn, int ndim, int dir, double dtdx,
  const struct gkyl_range *update_range, const struct gkyl_wave_geom *wg,
  struct gkyl_array *redo_fluct, struct gkyl_array *amdq, struct gkyl_array *apdq,
  struct gkyl_array *qout)
{
  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  int meqn = eqn->num_equations;

  double flux_zero[meqn]; // zero second order flux for redo
  for (int i=0; i<meqn; ++i)
    flux_zero[i] = 0.0;

  // compute first-order update in each cell
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxl);
    gkyl_copy_int_arr(ndim, iter.idx, idxr);

    idxr[dir] = idxr[dir]+1;
    long lidx = gkyl_range_idx(update_range, idxl);
    long ridx = gkyl_range_idx(update_range, idxr);

    double *amdq_local = gkyl_array_fetch(amdq, ridx);
    double *apdq_local = gkyl_array_fetch(apdq, lidx);

    const double *redo_fluct_c = gkyl_array_cfetch(redo_fluct, lidx);
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, idxl);
    if (redo_fluct_c[0] > 0.0) {
      calc_first_order_update(meqn, dtdx/cg->kappa, gkyl_array_fetch(qout, lidx),
        amdq_local, apdq_local, flux_zero, flux_zero);
    }
  }
}

static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

static inline void
wave_rescale(int meqn, double fact, double *w, double *ws)
{
  for (int i=0; i<meqn; ++i) ws[i] = w[i]*fact;
}

static inline void
calc_second_order_qflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

// this is the sign function for doubles
static inline int sign_double(double val) { return (0.0 < val) - (val < 0.0); }

static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

static inline void
wave_prop_second_order_flux(enum gkyl_wave_limiter limiter, double dtdx, 
  int ndim, int dir, int meqn, int mwaves, const struct gkyl_range *update_range,
  struct gkyl_range second_order_range, const struct gkyl_wave_geom *wg,
  struct gkyl_array *waves, struct gkyl_array *waves_scaled,
  struct gkyl_array *speeds, struct gkyl_array *flux2)
{
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &second_order_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_copy_int_arr(ndim, iter.idx, idxl);
    gkyl_copy_int_arr(ndim, iter.idx, idxr);

    idxl[dir] = idxc[dir]-1;
    idxr[dir] = idxc[dir]+1;
    long lidx = gkyl_range_idx(update_range, idxl);
    long cidx = gkyl_range_idx(update_range, idxc);
    long ridx = gkyl_range_idx(update_range, idxr);
    const double *s = gkyl_array_cfetch(speeds, cidx);
    double *wc = gkyl_array_fetch(waves, cidx);
    double *ws = gkyl_array_fetch(waves_scaled, cidx);
    const double *wl = gkyl_array_cfetch(waves, lidx);
    const double *wr = gkyl_array_cfetch(waves, ridx);

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
    
    double *flux2_d = gkyl_array_fetch(flux2, cidx);
    for (int mw=0; mw<mwaves; ++mw) {
      calc_second_order_qflux(meqn, dtdx/(0.5*(kappal+kappac)), s[mw], &ws[mw*meqn], flux2_d);
    }
  }
}

struct gkyl_wave_prop_status gkyl_wave_prop_advance_cu(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout);

double gkyl_wave_prop_max_dt_cu(const gkyl_wave_prop *wv,
  const struct gkyl_range *update_range, const struct gkyl_array *qin);
