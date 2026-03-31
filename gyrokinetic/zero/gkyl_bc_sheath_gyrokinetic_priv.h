#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_kernels.h>
#include <gkyl_bc_sheath_gyrokinetic_gyraze_surrogate.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*sheath_reflectedf_t)(const double *vmap, const double q2Dm,
  const double *phi, const double *phiWall, const double *vcut_fact, const double *f, double *fRefl);

typedef struct { sheath_reflectedf_t kernels[3]; } sheath_reflectedf_kern_list;  // For use in kernel tables.
typedef struct { sheath_reflectedf_kern_list list[4]; } edged_sheath_reflectedf_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_sheath_reflectedf_kern_list ser_sheath_reflect_list[] = {
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1, NULL },
          },
  },
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1, NULL },
          },
  },
};

// Function pointer type for surrogate sheath BC to determine vcut factor.
typedef void (*sheath_surrogate_t)(const double *vmap, const double *phi, const double *phi_wall, 
  const double *density, const double *temperature, double q2Dm, const double *bmag, const double *bimpact_angle, double *vcut_fact_out);

typedef struct { sheath_surrogate_t kernels[3]; } sheath_surrogate_kern_list;  // For use in kernel tables.
typedef struct { sheath_surrogate_kern_list list[4]; } edged_sheath_surrogate_kern_list;

// Function pointer for direct surrogate interface.
typedef void (*srgrz_eval_t) (const double *mu_new, int n, double phi, double phi_wall, double dens_e,
    double temp_e, double q2Dm, double bmag, double bimpact_angle, double *out);

// Serendipity surrogate kernels.
GKYL_CU_D
static const edged_sheath_surrogate_kern_list ser_sheath_surrogate_list[] = {
  { .list={
           { NULL, NULL },
           { bc_sheath_gyrokinetic_surrogate_lower_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_surrogate_lower_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_surrogate_lower_3x2v_ser_p1, NULL },
          },
  },
  { .list={
           { NULL, NULL },
           { bc_sheath_gyrokinetic_surrogate_upper_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_surrogate_upper_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_surrogate_upper_3x2v_ser_p1, NULL },
          },
  },
};

struct gkyl_bc_sheath_gyrokinetic_kernels {
  sheath_reflectedf_t reflectedf;  // reflectedf kernel.
  sheath_surrogate_t surrogate; // surrogate kernel to determine vcut factor.
};

// Primary struct in this updater.
struct gkyl_bc_sheath_gyrokinetic {
  int dir; // Direction perpendicular to the sheath boundary.
  int cdim; // Conf-space dimensionality.
  enum gkyl_edge_loc edge; // Lower or upper boundary.
  const struct gkyl_basis *basis; // Phase-space basis.
  bool use_gpu; // Whether to run on GPU.
  double q2Dm; // charge-to-mass ratio times 2.
  bool use_surrogate; // Whether to use surrogate sheath BC to determine vcut.
  bool use_surrogate_conv_check; // Whether to use the SVM classifier to check for convergence before using the surrogate.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels;  // reflectedf kernel.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r; // Skin and ghost ranges.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping.
  int vcut_fact_dim; // Dimensionality of vcut_fact array.
  struct gkyl_array *vcut_fact; // factor for mu dependent vcut in sheath BC.
  struct gkyl_basis vcut_fact_basis; // Basis for vcut_fact array (expansion in perpendicular config space and mu).
  struct gkyl_range vcut_fact_local; // Range for vcut_fact array.
  void (*update_vcut_fact) (const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
    const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
    const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r); // Function pointer to update vcut_fact array.
  srgrz_eval_t surrogate_eval; // Function pointer for direct surrogate interface.
};

void
gkyl_bc_gksheath_choose_reflectedf_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers);

GKYL_CU_D
static sheath_reflectedf_t
bc_gksheath_choose_reflectedf_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  sheath_reflectedf_t kern = NULL;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
      kern = ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kern = ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      kern = NULL;
      break;
  }
  assert(kern); // Reflection kernel must be non-null since we always use reflection BCs in the sheath updater.
  return kern;
}

GKYL_CU_D
static void
bc_gksheath_reflect(int dir, const struct gkyl_basis *basis, int cdim, double *out, const double *inp)
{
  basis->flip_odd_sign(dir, inp, out);
  basis->flip_odd_sign(cdim, out, out); // cdim is the vpar direction.
}

void 
gkyl_bc_gksheath_choose_surrogate_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, 
  bool use_conv_check, struct gkyl_bc_sheath_gyrokinetic_kernels *kers);

GKYL_CU_D
static sheath_surrogate_t
bc_gksheath_choose_surrogate_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, bool use_conv_check)
{
  
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  sheath_surrogate_t kern = NULL;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kern = ser_sheath_surrogate_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      kern = NULL;
      break;
  }
  assert(kern); // Surrogate kernel must be non-null if using surrogate.
  return kern;
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to apply the sheath BC.

 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_bc_sheath_gyrokinetic_advance_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r);

/**
 * CUDA device function to update the vcut_fact array.
 * 
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param dens Density array.
 * @param temp Temperature array.
 * @param bmag Magnetic field magnitude array.
 * @param bimpact_angle Magnetic field impact angle array.
 * @param conf_r Configuration space range (to index phi, dens, temp, bmag, and bimpact_angle).
 * 
 */
void bc_gksheath_update_vcut_fact_surrogate_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi, 
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle, const struct gkyl_range *conf_r);
#endif
