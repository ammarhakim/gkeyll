// Private header: not for direct use
#pragma once

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>

struct gkyl_vlasov_lte_correct
{
  int num_conf_basis; // Number of configuration-space basis functions
  int num_comp; // Number of components being corrected vdim+2 (n, V_drift, T/m)

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;
  int anderson_depth; // number of previous Anderson residuals to retain
  int anderson_hist_size; // number of Anderson residual/G(x) states
  int anderson_linsys_size; // Anderson constrained linear-system size
  int anderson_ncomp; // number of moment components in one Anderson state
  struct gkyl_array *aa_res_hist;
  struct gkyl_array *aa_g_hist;
  struct gkyl_array *aa_x_next;
  struct gkyl_array *aa_res_dot;
  struct gkyl_array *aa_moms_target_ho;

  struct gkyl_vlasov_lte_moments *moments_up;
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte;

  // error estimate (n, V_drift, T/m), 0 - success., num. correction iterations
  double *error; // absolute value of difference in cell averages between iteration and target
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  bool use_picard; // Use unaccelerated Picard iteration instead of Anderson acceleration?
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 

  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs 
  struct gkyl_array *abs_diff_moms;
};

// Status returned by one correction advance routine. The public status object
// is assembled from this after any Anderson fallback or backup path.
struct lte_correct_iter_status
{
  int niter; // Number of iterations completed by this advance routine
  bool ispositive; // true if density and temperature cell averages stayed positive
  double max_error; // maximum moment error at the end of this advance routine
};

// Weighted residual inner product used to build the Anderson Gram matrix.
// For residual histories r_i and r_j, compute
//
//   <r_i,r_j>_W = sum_x sum_d sum_k r_i,d,k(x) r_j,d,k(x) / s_d(x)^2,
//
// where d indexes the corrected moments, k indexes the configuration basis
// coefficients, and s_d is the local scaling used to make density and
// temperature relative errors while keeping small drift velocities absolute.
// The per-cell contributions are stored in aa_res_dot and summed with the
// existing array reduction machinery.
static double
lte_correct_res_dot(gkyl_vlasov_lte_correct *up,
  const struct gkyl_array *a, int aoff, const struct gkyl_array *b, int boff,
  const struct gkyl_array *moms_target, const struct gkyl_range *conf_range)
{
  int num_comp = up->num_comp;
  int nc = up->num_conf_basis;

  struct gkyl_range_iter biter;
  gkyl_range_iter_init(&biter, conf_range);
  while (gkyl_range_iter_next(&biter)) {
    long midx = gkyl_range_idx(conf_range, biter.idx);
    const double *a_d = (const double*) gkyl_array_cfetch(a, midx);
    const double *b_d = (const double*) gkyl_array_cfetch(b, midx);
    const double *moms_target_d = (const double*) gkyl_array_cfetch(moms_target, midx);
    double *res_dot_d = (double*) gkyl_array_fetch(up->aa_res_dot, midx);

    double dot = 0.0;
    for (int d=0; d<num_comp; ++d) {
      double scale = fabs(moms_target_d[d*nc]);
      if (d == 0 || d == num_comp-1)
        scale = scale > DBL_MIN ? scale : 1.0;
      else
        scale = scale < 1.0 ? 1.0 : scale;
      double w = 1.0/(scale*scale);
      for (int k=0; k<nc; ++k)
        dot += a_d[aoff+d*nc+k]*b_d[boff+d*nc+k]*w;
    }
    res_dot_d[0] = dot;
  }

  double dot[1] = { 0.0 };
  gkyl_array_reduce_range(dot, up->aa_res_dot, GKYL_SUM, conf_range);
  return dot[0];
}

// Compute Anderson coefficients for the fixed-point map x_{k+1}=G(x_k).
// With residuals r_i = G(x_i)-x_i, choose alpha to minimize
//
//   || sum_i alpha_i r_i ||_W^2  subject to  sum_i alpha_i = 1.
//
// This gives the constrained linear system
//
//   [ R^T W R   1 ][alpha ] = [0]
//   [   1^T     0 ][lambda]   [1],
//
// where R^T W R is the weighted residual Gram matrix computed above. The
// accelerated iterate is then x_{new} = sum_i alpha_i G(x_i).
static bool
lte_correct_anderson_coeffs(gkyl_vlasov_lte_correct *up,
  const struct gkyl_array *moms_target, const struct gkyl_range *conf_range,
  int hist_start, int hist_count, double alpha[8])
{
  for (int i=0; i<8; ++i)
    alpha[i] = 0.0;

  int linsys_size = hist_count+1;
  if (hist_count > up->anderson_hist_size ||
      linsys_size > up->anderson_linsys_size ||
      linsys_size > 8)
    return false;

  if (hist_count == 1) {
    alpha[0] = 1.0;
    return true;
  }

  double A_data[8*8] = { 0.0 };
  double rhs[8] = { 0.0 };
  long ipiv[8];
  double trace = 0.0;

  struct gkyl_mat A_mat;
  A_mat.nr = linsys_size;
  A_mat.nc = linsys_size;
  A_mat.data = A_data;

  struct gkyl_mat rhs_mat;
  rhs_mat.nr = linsys_size;
  rhs_mat.nc = 1;
  rhs_mat.data = rhs;

  for (int i=0; i<hist_count; ++i) {
    int si = (hist_start+i) % up->anderson_hist_size;
    int si_off = si*up->anderson_ncomp;
    for (int j=i; j<hist_count; ++j) {
      int sj = (hist_start+j) % up->anderson_hist_size;
      int sj_off = sj*up->anderson_ncomp;
      double dot = lte_correct_res_dot(up, up->aa_res_hist, si_off,
        up->aa_res_hist, sj_off, moms_target, conf_range);
      gkyl_mat_set(&A_mat, i, j, dot);
      gkyl_mat_set(&A_mat, j, i, dot);
    }
    trace += gkyl_mat_get(&A_mat, i, i);
  }

  double reg = fmax(1e-14*trace/(double) hist_count, 1e-30);
  for (int i=0; i<hist_count; ++i) {
    gkyl_mat_inc(&A_mat, i, i, reg);
    gkyl_mat_set(&A_mat, i, hist_count, 1.0);
    gkyl_mat_set(&A_mat, hist_count, i, 1.0);
  }
  gkyl_mat_set(&rhs_mat, hist_count, 0, 1.0);

  if (!gkyl_mat_linsolve_lu(&A_mat, &rhs_mat, ipiv))
    return false;

  double alpha_sum = 0.0;
  for (int i=0; i<hist_count; ++i) {
    double rhs_i = gkyl_mat_get(&rhs_mat, i, 0);
    if (!isfinite(rhs_i))
      return false;
    alpha[i] = rhs_i;
    alpha_sum += alpha[i];
  }

  if (fabs(alpha_sum) < DBL_MIN)
    return false;

  for (int i=0; i<hist_count; ++i)
    alpha[i] /= alpha_sum;

  return true;
}

static bool
lte_correct_update_error(gkyl_vlasov_lte_correct *up,
  const struct gkyl_array *moms_target, const struct gkyl_range *conf_local)
{
  int num_comp = up->num_comp;
  int nc = up->num_conf_basis;
  bool ispositive = true;

  if (up->use_gpu) {
    // We insure the reduction to find the maximum error is thread-safe on GPUs
    // by first calling a specialized kernel for computing the absolute value
    // of the difference of the cell averages, then calling reduce_range.
    gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(conf_local,
      num_comp, nc, moms_target, up->moms_iter, up->abs_diff_moms);
    gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_local);
    gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double)*num_comp, GKYL_CU_MEMCPY_D2H);
  }
  else {
    struct gkyl_range_iter biter;

    // Reset the maximum error.
    for (int i=0; i<num_comp; ++i) {
      up->error[i] = 0.0;
    }

    // Iterate over the input configuration-space range to find the maximum error.
    gkyl_range_iter_init(&biter, conf_local);
    while (gkyl_range_iter_next(&biter)) {
      long midx = gkyl_range_idx(conf_local, biter.idx);
      const double *moms_local = (const double*) gkyl_array_cfetch(up->moms_iter, midx);
      const double *moms_target_local = (const double*) gkyl_array_cfetch(moms_target, midx);

      // Check the error in the absolute value of the cell average.
      // Note: for density and temperature, this error is a relative error compared to the target moment value
      // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
      up->error[0] = fmax(fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc], fabs(up->error[0]));
      int T_idx = num_comp-1; // T/m is always the last component
      up->error[T_idx] = fmax(fabs(moms_local[T_idx*nc] - moms_target_local[T_idx*nc])/moms_target_local[T_idx*nc], fabs(up->error[T_idx]));

      // However, V_drift may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
      // absolute error if V_drift ~ O(1). Otherwise, we use relative error for V_drift.
      for (int d=1; d<num_comp-1; ++d) {
        if (fabs(moms_target_local[d*nc]) < 1.0) {
          up->error[d] = fmax(fabs(moms_local[d*nc] - moms_target_local[d*nc]), fabs(up->error[d]));
        }
        else {
          up->error[d] = fmax(fabs(moms_local[d*nc] - moms_target_local[d*nc])/moms_target_local[d*nc], fabs(up->error[d]));
        }
      }

      // Check if density and temperature are positive, if they aren't we will break out of the iteration.
      ispositive = (moms_local[0*nc] > 0.0) && ispositive;
      ispositive = (moms_local[T_idx*nc] > 0.0) && ispositive;
    }
  }

  return ispositive;
}

static struct lte_correct_iter_status
lte_correct_picard_advance(gkyl_vlasov_lte_correct *up,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  bool restart_from_target)
{
  int niter = 0;
  bool ispositive_f_lte = true;
  double max_error = 1.0;

  for (int i=0; i<up->num_comp; ++i) {
    up->error[i] = 1.0;
  }
  if (up->use_gpu) {
    gkyl_cu_memcpy(up->error_cu, up->error, sizeof(double)*up->num_comp, GKYL_CU_MEMCPY_H2D);
  }

  if (restart_from_target) {
    gkyl_vlasov_lte_proj_on_basis_advance(up->proj_lte,
      phase_local, conf_local, moms_target, f_lte);
  }

  // Clear the differences prior to iteration.
  gkyl_array_clear(up->d_moms, 0.0);
  gkyl_array_clear(up->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient for machine precision moments.
  while ((ispositive_f_lte) && ((niter < up->max_iter) && (max_error > up->eps))) {
    // 1. Calculate the LTE moments (n, V_drift, T) from the projected LTE distribution.
    gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new.
    // ddn = n_target - n.
    gkyl_array_set(up->dd_moms, -1.0, up->moms_iter);
    gkyl_array_accumulate(up->dd_moms, 1.0, moms_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0.
    // dm_new = dm_old + ddn.
    gkyl_array_accumulate(up->d_moms, 1.0, up->dd_moms);

    ispositive_f_lte = lte_correct_update_error(up, moms_target, conf_local);
    max_error = 0.0;
    for (int d=0; d<up->num_comp; ++d)
      max_error = fmax(max_error, up->error[d]);

    // c. Calculate  n^(k+1) = M^k + dM^(k+1).
    // n = n_target + dm_new.
    gkyl_array_set(up->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(up->moms_iter, 1.0, up->d_moms);

    // 2. Update the LTE distribution function using the corrected moments.
    // Projection routine also corrects the density before the next iteration.
    gkyl_vlasov_lte_proj_on_basis_advance(up->proj_lte,
      phase_local, conf_local, up->moms_iter, f_lte);

    niter += 1;
  }

  struct lte_correct_iter_status status;
  status.niter = niter;
  status.ispositive = ispositive_f_lte;
  status.max_error = max_error;
  return status;
}

static struct lte_correct_iter_status
lte_correct_anderson_advance(gkyl_vlasov_lte_correct *up,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  int niter = 0;
  bool ispositive_f_lte = true;
  double max_error = 1.0;

  for (int i=0; i<up->num_comp; ++i) {
    up->error[i] = 1.0;
  }
  if (up->use_gpu) {
    gkyl_cu_memcpy(up->error_cu, up->error, sizeof(double)*up->num_comp, GKYL_CU_MEMCPY_H2D);
  }

  // Initial corrected moments are the target moments. The iteration updates
  // this fixed-point variable, then projects f_lte from it.
  gkyl_array_set(up->d_moms, 1.0, moms_target);
  gkyl_array_clear(up->dd_moms, 0.0);

  const struct gkyl_array *moms_target_ho = moms_target;
  if (up->use_gpu) {
    gkyl_array_copy(up->aa_moms_target_ho, moms_target);
    moms_target_ho = up->aa_moms_target_ho;
  }

  // Anderson-accelerated fixed-point loop.
  while ((ispositive_f_lte) && ((niter < up->max_iter) && (max_error > up->eps))) {
    int iter = niter;

    // 1. Calculate the LTE moments (n, V_drift, T) from the projected LTE distribution.
    gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);

    // a. Calculate the fixed-point residual r^k = M_target - M(f_lte^k).
    gkyl_array_set(up->dd_moms, -1.0, up->moms_iter);
    gkyl_array_accumulate(up->dd_moms, 1.0, moms_target);

    int hist_slot = iter % up->anderson_hist_size;
    int hist_offset = hist_slot*up->anderson_ncomp;

    // Store the Anderson residual history r^k = G(x^k)-x^k.
    // Device arrays are copied through aa_x_next because the Anderson
    // history and coefficient solve are host-side.
    if (up->use_gpu) {
      gkyl_array_copy(up->aa_x_next, up->dd_moms);
      gkyl_array_set_offset(up->aa_res_hist, 1.0, up->aa_x_next, hist_offset);
    }
    else {
      gkyl_array_set_offset(up->aa_res_hist, 1.0, up->dd_moms, hist_offset);
    }

    ispositive_f_lte = lte_correct_update_error(up, moms_target, conf_local);
    max_error = 0.0;
    for (int d=0; d<up->num_comp; ++d)
      max_error = fmax(max_error, up->error[d]);

    niter += 1;

    if (!ispositive_f_lte)
      break;

    // b. Picard candidate G(x^k) = x^k + r^k.
    gkyl_array_set(up->moms_iter, 1.0, up->d_moms);
    gkyl_array_accumulate(up->moms_iter, 1.0, up->dd_moms);

    // Store the Picard image G(x^k) used to form the Anderson update
    // x_AA^{k+1} = sum_i alpha_i G(x_i).
    if (up->use_gpu) {
      gkyl_array_copy(up->aa_x_next, up->moms_iter);
      gkyl_array_set_offset(up->aa_g_hist, 1.0, up->aa_x_next, hist_offset);
    }
    else {
      gkyl_array_set_offset(up->aa_g_hist, 1.0, up->moms_iter, hist_offset);
    }

    int hist_count = iter+1 < up->anderson_hist_size ? iter+1 : up->anderson_hist_size;
    int hist_start = (iter+1 >= up->anderson_hist_size) ?
      (iter+1) % up->anderson_hist_size : 0;

    double alpha[8];
    bool aa_status = lte_correct_anderson_coeffs(up, moms_target_ho, conf_local,
      hist_start, hist_count, alpha);

    if (aa_status) {
      gkyl_array_clear(up->aa_x_next, 0.0);
      for (int i=0; i<hist_count; ++i) {
        int slot = (hist_start+i) % up->anderson_hist_size;
        gkyl_array_accumulate_offset(up->aa_x_next, alpha[i],
          up->aa_g_hist, slot*up->anderson_ncomp);
      }
    }
    else {
      gkyl_array_set_offset(up->aa_x_next, 1.0,
        up->aa_g_hist, hist_offset);
    }

    // Check density and temperature cell averages before accepting the
    // accelerated moments. If they fail positivity, use the unaccelerated
    // Picard image from the same history slot.
    bool aa_moms_positive = true;
    int T_idx = up->num_comp-1; // T/m is always the last component
    struct gkyl_range_iter biter;
    gkyl_range_iter_init(&biter, conf_local);
    while (gkyl_range_iter_next(&biter)) {
      long midx = gkyl_range_idx(conf_local, biter.idx);
      const double *moms_d = (const double*) gkyl_array_cfetch(up->aa_x_next, midx);
      aa_moms_positive = (moms_d[0] > 0.0) && aa_moms_positive;
      aa_moms_positive = (moms_d[T_idx*up->num_conf_basis] > 0.0) && aa_moms_positive;
    }

    if (!aa_moms_positive) {
      gkyl_array_set_offset(up->aa_x_next, 1.0, up->aa_g_hist, hist_offset);

      aa_moms_positive = true;
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)) {
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *moms_d = (const double*) gkyl_array_cfetch(up->aa_x_next, midx);
        aa_moms_positive = (moms_d[0] > 0.0) && aa_moms_positive;
        aa_moms_positive = (moms_d[T_idx*up->num_conf_basis] > 0.0) && aa_moms_positive;
      }
      if (!aa_moms_positive) {
        ispositive_f_lte = false;
        break;
      }
    }

    gkyl_array_copy(up->d_moms, up->aa_x_next);

    // 2. Update the LTE distribution function using the accelerated corrected moments.
    // Projection routine also corrects the density before the next iteration.
    gkyl_vlasov_lte_proj_on_basis_advance(up->proj_lte,
      phase_local, conf_local, up->d_moms, f_lte);
  }

  struct lte_correct_iter_status status;
  status.niter = niter;
  status.ispositive = ispositive_f_lte;
  status.max_error = max_error;
  return status;
}

static void
lte_correct_backup(gkyl_vlasov_lte_correct *up,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // Project the distribution function with the target moments. We correct
  // the density and then recompute moments/errors for this new projection.
  gkyl_vlasov_lte_proj_on_basis_advance(up->proj_lte,
    phase_local, conf_local, moms_target, f_lte);

  gkyl_vlasov_lte_moments_advance(up->moments_up, phase_local, conf_local, f_lte, up->moms_iter);
  lte_correct_update_error(up, moms_target, conf_local);
}
