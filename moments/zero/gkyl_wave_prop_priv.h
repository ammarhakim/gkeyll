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

GKYL_CU_DH
static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

GKYL_CU_DH
static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
GKYL_CU_DH
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
GKYL_CU_DH
static inline void
copy_wv_vec(int n, double * GKYL_RESTRICT out, const double * GKYL_RESTRICT inp)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

GKYL_CU_DH
static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

GKYL_CU_DH
static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT q, const double * GKYL_RESTRICT amdq_r, const double * GKYL_RESTRICT apdq_l, const double *fr, const double *fl)
{
  for (int i=0; i<meqn; ++i)
    q[i] = q[i] - dtdx*(apdq_l[i] + amdq_r[i] + fr[i] - fl[i]);
}

GKYL_CU_DH
static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

GKYL_CU_DH
static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

GKYL_CU_DH
static inline void
wave_rescale(int meqn, double fact, double *w, double *ws)
{
  for (int i=0; i<meqn; ++i) ws[i] = w[i]*fact;
}

GKYL_CU_DH
static inline void
calc_second_order_qflux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

// this is the sign function for doubles
GKYL_CU_DH
static inline int sign_double(double val) { return (0.0 < val) - (val < 0.0); }

GKYL_CU_DH
static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

#ifdef GKYL_HAVE_CUDA

/**
 * Compute wave-propagation update on the device. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param tm Current time
 * @param dt time-step
 * @param update_rng Range on which to compute.
 * @param phi Embedded mask
 * @param qin Input to updater
 * @param qout Solution at tm+dt
 */
struct gkyl_wave_prop_status gkyl_wave_prop_advance_cu(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *phi, const struct gkyl_array *qin, struct gkyl_array *qout);

/**
 * Compute an estimate of maximum stable time-step for given input
 * state 'qin' on device. 
 *
 * @param wv Updater object
 * @param qin Input to compute dt for
 * @return maximum stable time-step
 */
double gkyl_wave_prop_max_dt_cu(const gkyl_wave_prop *wv,
  const struct gkyl_range *update_range, const struct gkyl_array *qin);

#endif
