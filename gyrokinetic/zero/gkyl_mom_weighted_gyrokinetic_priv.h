#pragma once

// Private header for mom_weighted_gyrokinetic updater, not for direct use in user code.

#include <gkyl_mom_weighted_gyrokinetic.h>
#include <gkyl_mom_weighted_gyrokinetic_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

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
static const gkyl_gyrokinetic_mom_kern_list ser_m0_kernels[] = {
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
  struct gkyl_rect_grid *phase_grid;  // Phase-space grid.
  double mass;  // Species mass.
  double charge;  // Species charge.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_mom_weighted_gyrokinetic_kernels *kernels;
};
