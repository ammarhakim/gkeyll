#pragma once

// Private header for mom_weighted_gyrokinetic updater, not for direct use in user code.

#include <gkyl_mom_weighted_gyrokinetic.h>
#include <gkyl_mom_weighted_gyrokinetic_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// Kernel signature.
typedef void (*wmomf_t)(const double *dxv, const double *vmap, double m_, double q_, const double *bmag,
  const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct {
  wmomf_t kernels[3]; // One for each poly_order.
} gkyl_gyrokinetic_mom_kern_list;

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0_wgtx_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p1, mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M0_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m1_wgtx_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p1, mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M1_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_wgtx_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p1, mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M2_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2par_wgtx_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p1, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M2PAR_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2perp_wgtx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M2PERP_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum and kinetic energy computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2_wgtx_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p1, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum, and parallel and perpendicular kinetic energy computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2parm2perp_wgtx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p1, mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p1, mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Hamiltonian moment kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_hamiltonian_nowgt_kernels[] = {
  // 1x kernels
  { NULL, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x1v_ser_p1, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x1v_ser_p2 }, // 0
  { NULL, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p1, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p2 }, // 1
  // 2x kernels                                                                                             
  { NULL, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p1, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0_wgtx_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p1, int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m1_wgtx_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p1, int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M1_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2_wgtx_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p1, int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2par_wgtx_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p1, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2PAR_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2perp_wgtx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M2PERP_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2_wgtx_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p1, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2PARM2PERP_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2parm2perp_wgtx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p1, int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p1, int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_3x2v_ser_p1, NULL }, // 3
};

// Integrated Hamiltonian kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_hamiltonian_nowgt_kernels[] = {
  // 1x kernels
  { NULL, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x1v_ser_p1, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x1v_ser_p2 }, // 0
  { NULL, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p1, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p1, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, int_mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_3x2v_ser_p1, NULL }, // 3
};

struct gkyl_mom_weighted_gyrokinetic_kernels {
  wmomf_t wmom; // Kernels that computes weighted moment.
};

// Object definition.
struct gkyl_mom_weighted_gyrokinetic {
  const struct gkyl_rect_grid *phase_grid;  // Phase-space grid.
  double mass;  // Species mass.
  double charge;  // Species charge.
  enum gkyl_mom_weight_type wgt_type; // Type of weight.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_mom_weighted_gyrokinetic_kernels *kernels;
  int num_mom; // Number of moments to compute.
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void
mom_weighted_choose_kernel_cu(struct gkyl_mom_weighted_gyrokinetic_kernels *kernels, int *num_mom,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type,
  enum gkyl_mom_weight_type wgt_type, bool is_integrated);

void
gkyl_mom_weighted_gyrokinetic_advance_cu(struct gkyl_mom_weighted_gyrokinetic *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng, const struct gkyl_range *wgt_rng,
  struct gkyl_array *phi, struct gkyl_array *wgt, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT mout);
#endif

GKYL_CU_D
static void mom_weighted_gk_choose_kernel(struct gkyl_mom_weighted_gyrokinetic_kernels *kernels, int *num_mom,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type,
  enum gkyl_mom_weight_type wgt_type, bool is_integrated)
{
  enum gkyl_basis_type cbasis_type = cbasis->b_type, pbasis_type = pbasis->b_type;
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = pbasis->poly_order;

  // Choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels,
    *m2par_kernels, *m2perp_kernels,
    *m0m1m2_kernels,
    *m0m1m2parm2perp_kernels,
    *hamiltonian_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (wgt_type == GKYL_F_MOMENT_WEIGHT_NONE) {
        if (is_integrated) {
          hamiltonian_kernels = ser_int_hamiltonian_nowgt_kernels;
        }
        else {
          hamiltonian_kernels = ser_hamiltonian_nowgt_kernels;
        }
      }
      else if (wgt_type == GKYL_F_MOMENT_WEIGHT_CONF) {
        if (is_integrated) {
          m0_kernels = ser_int_m0_wgtx_kernels;
          m1_kernels = ser_int_m1_wgtx_kernels;
          m2_kernels = ser_int_m2_wgtx_kernels;
          m2par_kernels = ser_int_m2par_wgtx_kernels;
          m2perp_kernels = ser_int_m2perp_wgtx_kernels;
          m0m1m2_kernels = ser_int_m0m1m2_wgtx_kernels;
          m0m1m2parm2perp_kernels = ser_int_m0m1m2parm2perp_wgtx_kernels;
        }
        else {
          m0_kernels = ser_m0_wgtx_kernels;
          m1_kernels = ser_m1_wgtx_kernels;
          m2_kernels = ser_m2_wgtx_kernels;
          m2par_kernels = ser_m2par_wgtx_kernels;
          m2perp_kernels = ser_m2perp_wgtx_kernels;
          m0m1m2_kernels = ser_m0m1m2_wgtx_kernels;
          m0m1m2parm2perp_kernels = ser_m0m1m2parm2perp_wgtx_kernels;
        }
      }
      else {
        fprintf(stderr, "mom_weighted_gyrokinetic: Weight type not yet supported.\n");
      }
      break;

    default:
      assert(false);
      break;
  }

  if (mom_type == GKYL_F_MOMENT_M0) { // Density.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M1) { // Parallel momentum.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2) { // Total kinetic energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PAR) { // Parallel energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M2PERP) { // Perpendicular energy.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2) { 
    // Density, parallel momentum, and total energy computed together.
    assert(cv_index[cdim].vdim[vdim] != -1);   
    assert(NULL != m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 3;
  }
  else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    // Density, parallel momentum, parallel and perpendicular
    // kinetic energy computed together.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = vdim+2;
  }
  else if (mom_type == GKYL_F_MOMENT_HAMILTONIAN) {
    // M0, mass*M0 and total particle energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    kernels->wmom = hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 3;
  }
  else {
    // string not recognized
    printf("Error: requested moment %d.\n", mom_type);
    gkyl_exit("gkyl_mom_type_gyrokinetic: Unrecognized moment requested!");
  }
}
