#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_correct_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>

#include <assert.h>

struct gkyl_vlasov_lte_correct*
gkyl_vlasov_lte_correct_inew(const struct gkyl_vlasov_lte_correct_inp *inp)
{
  gkyl_vlasov_lte_correct *up = gkyl_malloc(sizeof(*up));

  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_gpu = inp->use_gpu;
  up->use_picard = inp->use_picard;
  up->use_last_converged = inp->use_last_converged;

  up->num_conf_basis = inp->conf_basis->num_basis;
  // (n, V_drift, T/m) being corrected
  // If the model is SR, V_drift is the spatial component of the four-velocity u_i = GammaV*V_drift
  int vdim = inp->phase_basis->ndim - inp->conf_basis->ndim;
  up->num_comp = vdim+2;
  up->anderson_depth = inp->anderson_depth > 0 ?
    inp->anderson_depth : 3;
  up->anderson_hist_size = up->anderson_depth+1;
  up->anderson_linsys_size = up->anderson_hist_size+1;
  up->anderson_ncomp = up->num_comp*up->num_conf_basis;
  if (!up->use_picard)
    assert(up->anderson_linsys_size <= 8);

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  if (up->use_gpu) {
    up->moms_iter = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
    // Two additional GPU-specific allocations for iterating over the grid to find the absolute value of 
    // the difference between the target and iterative moments, and the GPU-side array for performing the
    // thread-safe reduction to find the maximum error on the grid.
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_comp, conf_local_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[up->num_comp]));
  }
  else {
    up->moms_iter = gkyl_array_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
    up->d_moms = gkyl_array_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
    up->dd_moms = gkyl_array_new(GKYL_DOUBLE, up->num_comp*up->num_conf_basis, conf_local_ext_ncells);
  }
  if (!up->use_picard) {
    up->aa_res_hist = gkyl_array_new(GKYL_DOUBLE,
      up->anderson_hist_size*up->anderson_ncomp, conf_local_ext_ncells);
    up->aa_g_hist = gkyl_array_new(GKYL_DOUBLE,
      up->anderson_hist_size*up->anderson_ncomp, conf_local_ext_ncells);
    up->aa_x_next = gkyl_array_new(GKYL_DOUBLE, up->anderson_ncomp, conf_local_ext_ncells);
    up->aa_res_dot = gkyl_array_new(GKYL_DOUBLE, 1, conf_local_ext_ncells);
    up->aa_moms_target_ho = up->use_gpu ?
      gkyl_array_new(GKYL_DOUBLE, up->anderson_ncomp, conf_local_ext_ncells) : 0;
  }
  else {
    up->aa_res_hist = 0;
    up->aa_g_hist = 0;
    up->aa_x_next = 0;
    up->aa_res_dot = 0;
    up->aa_moms_target_ho = 0;
  }

  // Allocate host-side error for checking convergence and returning in the status object 
  up->error = gkyl_malloc(sizeof(double[up->num_comp]));

  // Moments structure 
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .gamma_inv = inp->gamma_inv,
    .h_ij = inp->h_ij,
    .h_ij_inv = inp->h_ij_inv,
    .det_h = inp->det_h,
    .use_extended_hamil_def = inp->use_extended_hamil_def,
    .effective_potential = inp->effective_potential,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_lte_moments_inew( &inp_mom );

  // Create a projection updater for projecting the LTE distribution function
  // Projection routine also corrects the density before returning 
  // the LTE distribution function.
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = inp->phase_grid,
    .vel_grid = inp->vel_grid,
    .conf_basis = inp->conf_basis,
    .vel_basis = inp->vel_basis,
    .phase_basis = inp->phase_basis,
    .conf_range = inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .phase_range = inp->phase_range,
    .vel_map = inp->vel_map,
    .hamil_range = inp->hamil_range,
    .hamil = inp->hamil,
    .model_id = inp->model_id,
    .gamma_inv = inp->gamma_inv,
    .quad_type = inp->quad_type,
    .h_ij = inp->h_ij,
    .h_ij_inv = inp->h_ij_inv,  
    .det_h = inp->det_h,
    .use_extended_hamil_def = inp->use_extended_hamil_def,
    .background_flows = inp->background_flows,
    .effective_potential = inp->effective_potential,
    .use_gpu = inp->use_gpu,
  };
  up->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  return up;
}

struct gkyl_vlasov_lte_correct_status
gkyl_vlasov_lte_correct_all_moments(gkyl_vlasov_lte_correct *up,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  struct lte_correct_iter_status iter_status = up->use_picard ?
    lte_correct_picard_advance(up, f_lte, moms_target, phase_local, conf_local, false) :
    lte_correct_anderson_advance(up, f_lte, moms_target, phase_local, conf_local);

  bool iter_succeeded = (iter_status.niter < up->max_iter) &&
    iter_status.ispositive && (iter_status.max_error < up->eps);
  if (!up->use_picard && !iter_succeeded && !up->use_last_converged) {
    int anderson_niter = iter_status.niter;
    iter_status = lte_correct_picard_advance(up, f_lte, moms_target,
      phase_local, conf_local, true);
    iter_succeeded = (iter_status.niter < up->max_iter) &&
      iter_status.ispositive && (iter_status.max_error < up->eps);
    iter_status.niter += anderson_niter;
  }

  bool corr_status = iter_succeeded ? 0 : 1;

  // If the algorithm fails to converge and we are *not* using the results of
  // the failed convergence, project the distribution function with the target
  // moments after correcting the density.
  if (corr_status == 1 && !up->use_last_converged) {
    lte_correct_backup(up, f_lte, moms_target, phase_local, conf_local);
  }

  struct gkyl_vlasov_lte_correct_status status;
  status.iter_converged = corr_status;
  status.num_iter = iter_status.niter;
  for (int i=0; i<up->num_comp; ++i) {
    status.error[i] = up->error[i];
  }
  return status;
}

void 
gkyl_vlasov_lte_correct_release(gkyl_vlasov_lte_correct *up)
{
  gkyl_array_release(up->moms_iter);
  gkyl_array_release(up->d_moms);
  gkyl_array_release(up->dd_moms);
  if (!up->use_picard) {
    gkyl_array_release(up->aa_res_hist);
    gkyl_array_release(up->aa_g_hist);
    gkyl_array_release(up->aa_x_next);
    gkyl_array_release(up->aa_res_dot);
    if (up->use_gpu)
      gkyl_array_release(up->aa_moms_target_ho);
  }
  if (up->use_gpu) {
    gkyl_array_release(up->abs_diff_moms);
    gkyl_cu_free(up->error_cu);
  }
  gkyl_free(up->error);

  gkyl_vlasov_lte_moments_release(up->moments_up);
  gkyl_vlasov_lte_proj_on_basis_release(up->proj_lte);

  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void 
gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  assert(false);
}

#endif
