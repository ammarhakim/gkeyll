/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_loss_cone_mask_gyrokinetic_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

// Kernel to compute Dbmag_quad = bmag_peak - bmag at quadrature nodes.
// bmag_peak is a per-field-line array (1D for 2x, scalar for 1x).
// For 1x: bmag_peak has a single value at index 0.
// For 2x: bmag_peak varies with psi (x-direction).
__global__ static void
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu_ker(int cdim, struct gkyl_range conf_range,
  struct gkyl_range bmag_peak_range, const struct gkyl_array* basis_at_ords_conf,
  const struct gkyl_array* bmag, const struct gkyl_array* bmag_peak,
  const struct gkyl_basis* bmag_peak_basis, struct gkyl_array* Dbmag_quad_out)
{    
  int num_basis_conf = basis_at_ords_conf->ncomp;
  int tot_quad_conf = basis_at_ords_conf->size;

  int cidx[GKYL_MAX_CDIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long linidx = gkyl_range_idx(&conf_range, cidx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(bmag, linidx);
    double *Dbmag_quad_d = (double*) gkyl_array_fetch(Dbmag_quad_out, linidx);

    // Get bmag_peak for this field line.
    // For 1x: single value (index 0).
    // For 2x: varies with psi, so use cidx[0].
    double bmag_peak_val;
    if (cdim == 1) {
      // 1x case: single value.
      const double *bmag_peak_d = (const double*) gkyl_array_cfetch(bmag_peak, 0);
      bmag_peak_val = bmag_peak_d[0]; // Just the constant coefficient.
    } else {
      // 2x case: evaluate bmag_peak at this psi cell.
      int psi_idx[1] = {cidx[0]};
      long psi_linidx = gkyl_range_idx(&bmag_peak_range, psi_idx);
      const double *bmag_peak_d = (const double*) gkyl_array_cfetch(bmag_peak, psi_linidx);
      // Evaluate at cell center (logical coord 0).
      double xc[1] = {0.0};
      bmag_peak_val = bmag_peak_basis->eval_expand(xc, bmag_peak_d);
    }

    // Sum over basis to get bmag at quadrature points, then compute difference.
    for (int n=0; n<tot_quad_conf; ++n) {
      const double *b_ord = (const double*) gkyl_array_cfetch(basis_at_ords_conf, n);

      double bmag_quad = 0.0;
      for (int k=0; k<num_basis_conf; ++k) {
        bmag_quad += bmag_d[k]*b_ord[k];
      }
      Dbmag_quad_d[n] = bmag_peak_val - bmag_quad;
    }
  }
}

void 
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag,
  struct gkyl_array *Dbmag_quad, const struct gkyl_array *bmag_peak)
{
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu_ker<<<nblocks, nthreads>>>(up->cdim, *conf_range,
    *up->bmag_max_range, up->basis_at_ords_conf->on_dev, bmag->on_dev, bmag_peak->on_dev,
    up->bmag_max_basis_on_dev, Dbmag_quad->on_dev);
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads 
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

// Kernel to compute qDphiDbmag_quad = charge*(phi-phi_m)/(bmag_max-bmag) at quadrature nodes.
// Supports per-field-line phi_m lookup for 2x mirrors.
__global__ static void
gkyl_loss_cone_mask_gyrokinetic_qDphiDbmag_quad_ker(int cdim, struct gkyl_range conf_range, 
  struct gkyl_range phi_m_range, const struct gkyl_array* basis_at_ords_conf, 
  const struct gkyl_basis* phi_m_basis, double charge, bool is_tandem,
  const struct gkyl_array* phi, const struct gkyl_array* phi_m, const struct gkyl_array* phi_tandem,
  const struct gkyl_array* Dbmag_quad, const struct gkyl_array* Dbmag_quad_wall, 
  const struct gkyl_array* Dbmag_quad_tandem,
  struct gkyl_array* qDphiDbmag_quad, struct gkyl_array* qDphiDbmag_quad_wall,
  struct gkyl_array* qDphiDbmag_quad_tandem)
{
  int num_basis_conf = basis_at_ords_conf->ncomp;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linc2 goes from 0 to tot_quad_conf= basis_at_ords_conf->size.
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linidx = gkyl_range_idx(&conf_range, cidx);

    const double *phi_d = (const double*) gkyl_array_cfetch(phi, linidx);
    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx);
    const double *Dbmag_quad_wall_d = (const double*) gkyl_array_cfetch(Dbmag_quad_wall, linidx);
    const double *Dbmag_quad_tandem_d = is_tandem ? 
      (const double*) gkyl_array_cfetch(Dbmag_quad_tandem, linidx) : Dbmag_quad_d;

    // Get phi_m value for this field line.
    // For 1x: single value (phi_m is a scalar stored as p=0 DG expansion).
    // For 2x: varies with psi, evaluate at this psi cell.
    double phi_m_val, phi_tandem_m_val;
    if (cdim == 1) {
      // 1x case: single scalar value stored as p=0 DG expansion.
      const double *phi_m_d = (const double*) gkyl_array_cfetch(phi_m, 0);
      phi_m_val = phi_m_d[0];
      if (is_tandem) {
        const double *phi_tandem_m_d = (const double*) gkyl_array_cfetch(phi_tandem, 0);
        phi_tandem_m_val = phi_tandem_m_d[0];
      }
    } else {
      // 2x case: evaluate phi_m at this psi cell center.
      int psi_idx[1] = {cidx[0]};
      long phi_m_linidx = gkyl_range_idx(&phi_m_range, psi_idx);
      const double *phi_m_d = (const double*) gkyl_array_cfetch(phi_m, phi_m_linidx);
      // Evaluate at cell center (logical coord 0).
      double xc[1] = {0.0};
      phi_m_val = phi_m_basis->eval_expand(xc, phi_m_d);
      if (is_tandem) {
        const double *phi_tandem_m_d = (const double*) gkyl_array_cfetch(phi_tandem, phi_m_linidx);
        phi_tandem_m_val = phi_m_basis->eval_expand(xc, phi_tandem_m_d);
      }
    }

    // Sum over basis at configuration-space quadrature points. 
    const double *b_ord = (const double*) gkyl_array_cfetch(basis_at_ords_conf, linc2);
    double phi_quad = 0;
    for (int k=0; k<num_basis_conf; ++k)
      phi_quad += phi_d[k]*b_ord[k];

    // Potential energy term at each quadrature point.
    double *qDphiDbmag_quad_d = (double*) gkyl_array_fetch(qDphiDbmag_quad, linidx);
    double *qDphiDbmag_quad_wall_d = (double*) gkyl_array_fetch(qDphiDbmag_quad_wall, linidx);
    
    if (Dbmag_quad_d[linc2] > 0.0)
      qDphiDbmag_quad_d[linc2] = charge*(phi_quad-phi_m_val)/Dbmag_quad_d[linc2];
    else
      qDphiDbmag_quad_d[linc2] = 0.0;

    if (Dbmag_quad_wall_d[linc2] > 0.0)
      qDphiDbmag_quad_wall_d[linc2] = charge*phi_quad/Dbmag_quad_wall_d[linc2];
    else
      qDphiDbmag_quad_wall_d[linc2] = 0.0;

    if (is_tandem) {
      double *qDphiDbmag_quad_tandem_d = (double*) gkyl_array_fetch(qDphiDbmag_quad_tandem, linidx);
      if (Dbmag_quad_tandem_d[linc2] > 0.0)
        qDphiDbmag_quad_tandem_d[linc2] = charge*(phi_quad-phi_tandem_m_val)/Dbmag_quad_tandem_d[linc2];
      else
        qDphiDbmag_quad_tandem_d[linc2] = 0.0;
    }
  }
}

// Cellwise kernel: determines if a cell is trapped or lost without quadrature.
// Supports tandem mirrors and per-field-line z-coordinate lookup.
__global__ static void
gkyl_loss_cone_mask_gyrokinetic_ker(int cdim, struct gkyl_rect_grid grid_phase,
  struct gkyl_range phase_range, struct gkyl_range conf_range, struct gkyl_range vel_range,
  struct gkyl_range bmag_max_range, const struct gkyl_basis* bmag_max_basis, bool is_tandem,
  double mass, const struct gkyl_array* phase_ordinates,
  const struct gkyl_array* bmag_max_z_coord, const struct gkyl_array* bmag_tandem_z_coord,
  const struct gkyl_array* qDphiDbmag_quad, const struct gkyl_array* qDphiDbmag_quad_wall,
  const struct gkyl_array* qDphiDbmag_quad_tandem,
  const struct gkyl_array* Dbmag_quad, const struct gkyl_array* Dbmag_quad_wall,
  const struct gkyl_array* Dbmag_quad_tandem,
  const int *p2c_qidx, struct gkyl_array* vmap, struct gkyl_basis* vmap_basis, struct gkyl_array* mask_out)
{
  int pdim = phase_range.ndim;
  int vdim = pdim-cdim;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  int tot_phase_quad = phase_ordinates->size;

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // Get configuration-space linear index.
    for (unsigned int k = 0; k < cdim; k++) cidx[k] = pidx[k];
    long linidx_conf = gkyl_range_idx(&conf_range, cidx);

    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx_conf);
    const double *Dbmag_quad_wall_d = (const double*) gkyl_array_cfetch(Dbmag_quad_wall, linidx_conf);
    const double *Dbmag_quad_tandem_d = is_tandem ? 
      (const double*) gkyl_array_cfetch(Dbmag_quad_tandem, linidx_conf) : Dbmag_quad_d;
    const double *qDphiDbmag_quad_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad, linidx_conf);
    const double *qDphiDbmag_quad_wall_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad_wall, linidx_conf);
    const double *qDphiDbmag_quad_tandem_d = is_tandem ?
      (const double*) gkyl_array_cfetch(qDphiDbmag_quad_tandem, linidx_conf) : qDphiDbmag_quad_d;

    // Get z-coordinates for field-line specific values.
    double bmag_max_z_val, bmag_tandem_z_val;
    if (cdim == 1) {
      const double *bmag_max_z_d = (const double*) gkyl_array_cfetch(bmag_max_z_coord, 0);
      bmag_max_z_val = bmag_max_z_d[0];
      if (is_tandem) {
        const double *bmag_tandem_z_d = (const double*) gkyl_array_cfetch(bmag_tandem_z_coord, 0);
        bmag_tandem_z_val = bmag_tandem_z_d[0];
      }
    } else {
      int psi_idx[1] = {cidx[0]};
      long psi_linidx = gkyl_range_idx(&bmag_max_range, psi_idx);
      const double *bmag_max_z_d = (const double*) gkyl_array_cfetch(bmag_max_z_coord, psi_linidx);
      double xc_log[1] = {0.0};
      bmag_max_z_val = bmag_max_basis->eval_expand(xc_log, bmag_max_z_d);
      if (is_tandem) {
        const double *bmag_tandem_z_d = (const double*) gkyl_array_cfetch(bmag_tandem_z_coord, psi_linidx);
        bmag_tandem_z_val = bmag_max_basis->eval_expand(xc_log, bmag_tandem_z_d);
      }
    }

    gkyl_rect_grid_cell_center(&grid_phase, pidx, xc);
    long linidx_phase = gkyl_range_idx(&phase_range, pidx);
    double *mask_d = (double*) gkyl_array_fetch(mask_out, linidx_phase);

    for (int d = cdim; d < pdim; d++) vidx[d-cdim] = pidx[d];
    long linidx_vel = gkyl_range_idx(&vel_range, vidx);
    const double *vmap_d = (const double*) gkyl_array_cfetch(vmap, linidx_vel);

    mask_d[0] = 1.0; // In this case the mask has ncomp=1.

    for (int n=0; n<tot_phase_quad; ++n) {
      int cqidx = p2c_qidx[n];

      const double *xcomp_d = (const double*) gkyl_array_cfetch(phase_ordinates, n);

      // Convert comp position coordinate to phys pos coord.
      log_to_comp(cdim, xcomp_d, grid_phase.dx, xc, xmu);
  
      // Convert comp velocity coordinate to phys velocity coord.
      double xcomp[1];
      for (int vd = 0; vd < vdim; vd++) {
        xcomp[0] = xcomp_d[cdim+vd];
        xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
      }
  
      // KEparDbmag = 0.5*mass*pow(vpar,2)/(bmag_peak-bmag).
      double KEparDbmag = 0.0, KEparDbmag_wall = 0.0, KEparDbmag_tandem = 0.0;
      if (Dbmag_quad_d[cqidx] > 0.0)
        KEparDbmag = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_d[cqidx];
  
      if (Dbmag_quad_wall_d[cqidx] > 0.0)
        KEparDbmag_wall = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_wall_d[cqidx];

      if (is_tandem && Dbmag_quad_tandem_d[cqidx] > 0.0)
        KEparDbmag_tandem = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_tandem_d[cqidx];

      double mu_bound = GKYL_MAX2(0.0, KEparDbmag+qDphiDbmag_quad_d[cqidx]);
      double mu_bound_wall = GKYL_MAX2(0.0, -(KEparDbmag_wall+qDphiDbmag_quad_wall_d[cqidx]));
      double mu_bound_tandem = is_tandem ? GKYL_MAX2(0.0, KEparDbmag_tandem+qDphiDbmag_quad_tandem_d[cqidx]) : 0.0;

      bool is_trapped;
      if (is_tandem) {
        // Tandem mirror trapping condition.
        bool in_outer_cell = fabs(xmu[cdim-1]) < fabs(bmag_max_z_val) &&
                             fabs(xmu[cdim-1]) > fabs(bmag_tandem_z_val);
        bool in_central_cell = fabs(xmu[cdim-1]) <= fabs(bmag_tandem_z_val);
        
        if (in_outer_cell) {
          is_trapped = mu_bound < xmu[cdim+1];
        } else if (in_central_cell) {
          double mu_bound_min = GKYL_MIN2(mu_bound, mu_bound_tandem);
          is_trapped = mu_bound_min < xmu[cdim+1];
        } else {
          is_trapped = mu_bound_wall > xmu[cdim+1] && fabs(xmu[cdim-1]) >= fabs(bmag_max_z_val);
        }
      } else {
        // Single mirror case.
        is_trapped = (mu_bound < xmu[cdim+1] && fabs(xmu[cdim-1]) < fabs(bmag_max_z_val)) ||
                     (mu_bound_wall > xmu[cdim+1] && fabs(xmu[cdim-1]) >= fabs(bmag_max_z_val));
      }
  
      if (!is_trapped) {
        mask_d[0] = 0.0;
        break;
      }
    }
  }
}

// Quadrature kernel: computes mask at phase-space quadrature nodes.
// Supports tandem mirrors and per-field-line z-coordinate lookup.
__global__ static void
gkyl_loss_cone_mask_gyrokinetic_quad_ker(int cdim, struct gkyl_rect_grid grid_phase,
  struct gkyl_range phase_range, struct gkyl_range conf_range, struct gkyl_range vel_range,
  struct gkyl_range bmag_max_range, const struct gkyl_basis* bmag_max_basis, bool is_tandem,
  double mass, double norm_fac, const struct gkyl_array* phase_ordinates, 
  const struct gkyl_array* bmag_max_z_coord, const struct gkyl_array* bmag_tandem_z_coord,
  const struct gkyl_array* qDphiDbmag_quad, const struct gkyl_array* qDphiDbmag_quad_wall,
  const struct gkyl_array* qDphiDbmag_quad_tandem,
  const struct gkyl_array* Dbmag_quad, const struct gkyl_array* Dbmag_quad_wall,
  const struct gkyl_array* Dbmag_quad_tandem,
  const int *p2c_qidx, struct gkyl_array* vmap, struct gkyl_basis* vmap_basis, struct gkyl_array* mask_out_quad)
{
  int pdim = phase_range.ndim;
  int vdim = pdim-cdim;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.0};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  // 2D thread grid
  // linc2 goes from 0 to tot_quad_phase
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // Get configuration-space linear index.
    for (unsigned int k = 0; k < cdim; k++) cidx[k] = pidx[k];

    long linidx_conf = gkyl_range_idx(&conf_range, cidx);

    const double *Dbmag_quad_d = (const double*) gkyl_array_cfetch(Dbmag_quad, linidx_conf);
    const double *Dbmag_quad_wall_d = (const double*) gkyl_array_cfetch(Dbmag_quad_wall, linidx_conf);
    const double *Dbmag_quad_tandem_d = is_tandem ?
      (const double*) gkyl_array_cfetch(Dbmag_quad_tandem, linidx_conf) : Dbmag_quad_d;
    const double *qDphiDbmag_quad_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad, linidx_conf);
    const double *qDphiDbmag_quad_wall_d = (const double*) gkyl_array_cfetch(qDphiDbmag_quad_wall, linidx_conf);
    const double *qDphiDbmag_quad_tandem_d = is_tandem ?
      (const double*) gkyl_array_cfetch(qDphiDbmag_quad_tandem, linidx_conf) : qDphiDbmag_quad_d;

    // Get z-coordinates for field-line specific values.
    double bmag_max_z_val, bmag_tandem_z_val;
    if (cdim == 1) {
      const double *bmag_max_z_d = (const double*) gkyl_array_cfetch(bmag_max_z_coord, 0);
      bmag_max_z_val = bmag_max_z_d[0];
      if (is_tandem) {
        const double *bmag_tandem_z_d = (const double*) gkyl_array_cfetch(bmag_tandem_z_coord, 0);
        bmag_tandem_z_val = bmag_tandem_z_d[0];
      }
    } else {
      int psi_idx[1] = {cidx[0]};
      long psi_linidx = gkyl_range_idx(&bmag_max_range, psi_idx);
      const double *bmag_max_z_d = (const double*) gkyl_array_cfetch(bmag_max_z_coord, psi_linidx);
      double xc_log[1] = {0.0};
      bmag_max_z_val = bmag_max_basis->eval_expand(xc_log, bmag_max_z_d);
      if (is_tandem) {
        const double *bmag_tandem_z_d = (const double*) gkyl_array_cfetch(bmag_tandem_z_coord, psi_linidx);
        bmag_tandem_z_val = bmag_max_basis->eval_expand(xc_log, bmag_tandem_z_d);
      }
    }

    gkyl_rect_grid_cell_center(&grid_phase, pidx, xc);
    long linidx_phase = gkyl_range_idx(&phase_range, pidx);

    int cqidx = p2c_qidx[linc2];
    for (int d = cdim; d < pdim; d++) vidx[d-cdim] = pidx[d];

    long linidx_vel = gkyl_range_idx(&vel_range, vidx);
    const double *vmap_d = (const double*) gkyl_array_cfetch(vmap, linidx_vel);
    const double *xcomp_d = (const double*) gkyl_array_cfetch(phase_ordinates, linc2);

    // Convert comp position coordinate to phys pos coord.
    log_to_comp(cdim, xcomp_d, grid_phase.dx, xc, xmu);

    // Convert comp velocity coordinate to phys velocity coord.
    double xcomp[1];
    for (int vd = 0; vd < vdim; vd++) {
      xcomp[0] = xcomp_d[cdim+vd];
      xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
    }

    // KEparDbmag = 0.5*mass*pow(vpar,2)/(bmag_peak-bmag).
    double KEparDbmag = 0.0, KEparDbmag_wall = 0.0, KEparDbmag_tandem = 0.0;
    if (Dbmag_quad_d[cqidx] > 0.0)
      KEparDbmag = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_d[cqidx];

    if (Dbmag_quad_wall_d[cqidx] > 0.0)
      KEparDbmag_wall = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_wall_d[cqidx];

    if (is_tandem && Dbmag_quad_tandem_d[cqidx] > 0.0)
      KEparDbmag_tandem = 0.5*mass*pow(xmu[cdim], 2.0)/Dbmag_quad_tandem_d[cqidx];

    double mu_bound = GKYL_MAX2(0.0, KEparDbmag+qDphiDbmag_quad_d[cqidx]);
    double mu_bound_wall = GKYL_MAX2(0.0, -(KEparDbmag_wall+qDphiDbmag_quad_wall_d[cqidx]));
    double mu_bound_tandem = is_tandem ? GKYL_MAX2(0.0, KEparDbmag_tandem+qDphiDbmag_quad_tandem_d[cqidx]) : 0.0;

    double *fq = (double*) gkyl_array_fetch(mask_out_quad, linidx_phase);
    
    if (is_tandem) {
      // Tandem mirror trapping condition.
      bool in_outer_cell = fabs(xmu[cdim-1]) < fabs(bmag_max_z_val) &&
                           fabs(xmu[cdim-1]) > fabs(bmag_tandem_z_val);
      bool in_central_cell = fabs(xmu[cdim-1]) <= fabs(bmag_tandem_z_val);
      
      if (in_outer_cell) {
        fq[linc2] = (mu_bound < xmu[cdim+1]) ? norm_fac : 0.0;
      } else if (in_central_cell) {
        double mu_bound_min = GKYL_MIN2(mu_bound, mu_bound_tandem);
        fq[linc2] = (mu_bound_min < xmu[cdim+1]) ? norm_fac : 0.0;
      } else {
        fq[linc2] = (mu_bound_wall > xmu[cdim+1] && fabs(xmu[cdim-1]) >= fabs(bmag_max_z_val)) ? norm_fac : 0.0;
      }
    } else {
      // Single mirror case.
      if (mu_bound < xmu[cdim+1] && fabs(xmu[cdim-1]) < fabs(bmag_max_z_val))
        fq[linc2] = norm_fac;
      else if (mu_bound_wall > xmu[cdim+1] && fabs(xmu[cdim-1]) >= fabs(bmag_max_z_val))
        fq[linc2] = norm_fac;
      else
        fq[linc2] = 0.0;
    }
  }
}

void
gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *phi, const struct gkyl_array *phi_m, const struct gkyl_array *phi_tandem,
  struct gkyl_array *mask_out)
{
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_quad_conf = up->basis_at_ords_conf->size;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_quad_conf);

  // Compute qDphiDbmag at quadrature points.
  gkyl_loss_cone_mask_gyrokinetic_qDphiDbmag_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(
    up->cdim, *conf_range, *up->bmag_max_range, 
    up->basis_at_ords_conf->on_dev, up->bmag_max_basis_on_dev, up->charge, up->is_tandem,
    phi->on_dev, phi_m->on_dev, phi_tandem->on_dev,
    up->Dbmag_quad->on_dev, up->Dbmag_quad_wall->on_dev, up->Dbmag_quad_tandem->on_dev,
    up->qDphiDbmag_quad->on_dev, up->qDphiDbmag_quad_wall->on_dev, up->qDphiDbmag_quad_tandem->on_dev);

  const struct gkyl_velocity_map *gvm = up->vel_map;

  if (up->cellwise_trap_loss) {
    // Don't do quadrature.
    int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
    gkyl_loss_cone_mask_gyrokinetic_ker<<<nblocks, nthreads>>>(up->cdim, *up->grid_phase, *phase_range, *conf_range,
      gvm->local_ext_vel, *up->bmag_max_range, up->bmag_max_basis_on_dev, up->is_tandem,
      up->mass, up->ordinates_phase->on_dev,
      up->bmag_max_z_coord->on_dev, up->bmag_tandem_z_coord->on_dev, 
      up->qDphiDbmag_quad->on_dev, up->qDphiDbmag_quad_wall->on_dev, up->qDphiDbmag_quad_tandem->on_dev,
      up->Dbmag_quad->on_dev, up->Dbmag_quad_wall->on_dev, up->Dbmag_quad_tandem->on_dev, 
      up->p2c_qidx, gvm->vmap->on_dev, gvm->vmap_basis, mask_out->on_dev);
  } else {
    // Use quadrature.
    dim3 dimGrid, dimBlock;
    int tot_quad_phase = up->basis_at_ords_phase->size;
    gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_quad_phase);

    gkyl_loss_cone_mask_gyrokinetic_quad_ker<<<dimGrid, dimBlock>>>(up->cdim, *up->grid_phase, *phase_range, *conf_range,
      gvm->local_ext_vel, *up->bmag_max_range, up->bmag_max_basis_on_dev, up->is_tandem,
      up->mass, up->norm_fac, up->ordinates_phase->on_dev,
      up->bmag_max_z_coord->on_dev, up->bmag_tandem_z_coord->on_dev,
      up->qDphiDbmag_quad->on_dev, up->qDphiDbmag_quad_wall->on_dev, up->qDphiDbmag_quad_tandem->on_dev,
      up->Dbmag_quad->on_dev, up->Dbmag_quad_wall->on_dev, up->Dbmag_quad_tandem->on_dev,
      up->p2c_qidx, gvm->vmap->on_dev, gvm->vmap_basis, up->mask_out_quad->on_dev);

    // Call cublas to do the matrix multiplication nodal to modal conversion.
    gkyl_mat_mm_array(up->phase_nodal_to_modal_mem, up->mask_out_quad, mask_out);
  }
}
