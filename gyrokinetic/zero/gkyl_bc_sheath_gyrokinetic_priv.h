#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include "gkyl_bc_sheath_gyrokinetic.h"
#include "gkyl_bc_sheath_gyrokinetic_kernels.h"
#include <assert.h>
#include <kann.h>
#include <gkyl_kann_net.h>

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

// Function pointer type for surrogate sheath BC to determine vcut factor from pre-computed NN output.
typedef void (*sheath_surrogate_t)(const double *vmap, const double *nn_out,
  const double *temperature, const double *bmag, double *vcut_fact_out);

// Function pointer type for the infer kernel: evaluates NN features and runs gkyl_kann_net_apply.
typedef void (*sheath_infer_t)(struct gkyl_kann_net *net, struct gkyl_kn_vec *inp_k, struct gkyl_kn_vec *out_k,
  const double *vmap, const double *phi, const double *phi_wall,
  const double *density, const double *temperature, const double *bmag,
  const double *bimpact_angle, double *out);

typedef struct { sheath_surrogate_t kernels[3]; } sheath_surrogate_kern_list;  // For use in kernel tables.
typedef struct { sheath_surrogate_kern_list dim_list[4]; } edged_sheath_surrogate_kern_list;

typedef struct { sheath_infer_t kernels[3]; } sheath_infer_kern_list;
typedef struct { sheath_infer_kern_list dim_list[4]; } edged_sheath_infer_kern_list;

// Serendipity surrogate kernels (vcut_fact projection using pre-computed nn_out).
GKYL_CU_D
static const edged_sheath_surrogate_kern_list ser_sheath_vcut_calc_list[] = {
  { .dim_list={
      { NULL, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_lower_1x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_lower_2x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_lower_3x2v_ser_p1, NULL },
    },
  },
  { .dim_list={
      { NULL, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_upper_1x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_upper_2x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_vparcut_calc_upper_3x2v_ser_p1, NULL },
    },
  },
};

// Serendipity infer kernels (NN feature evaluation + gkyl_kann_net_apply).
GKYL_CU_D
static const edged_sheath_infer_kern_list ser_sheath_infer_list[] = {
  { .dim_list={
      { NULL, NULL },
      { bc_sheath_gyrokinetic_infer_lower_1x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_infer_lower_2x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_infer_lower_3x2v_ser_p1, NULL },
    },
  },
  { .dim_list={
      { NULL, NULL },
      { bc_sheath_gyrokinetic_infer_upper_1x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_infer_upper_2x2v_ser_p1, NULL },
      { bc_sheath_gyrokinetic_infer_upper_3x2v_ser_p1, NULL },
    },
  },
};

struct gkyl_bc_sheath_gyrokinetic_kernels {
  sheath_reflectedf_t reflectedf; // reflectedf kernel.
  sheath_infer_t infer; // NN inference kernel (fills raw NN output).
  sheath_surrogate_t vcut_calc; // vcut_fact DG projection kernel (uses pre-computed NN output).
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
  kann_t *kann_model; // Loaded KANN model for the surrogate (CPU only; NULL on GPU or when not using surrogate).
  struct gkyl_kann_net *kann_net; // KANN sheath surrogate.
  struct gkyl_kn_vec *kann_inp, *kann_out; // Input and output vector structures.
  struct gkyl_array *kann_infer_xy_out; // Output array for KANN on a perpendicular plane.
  double *nn_out; // Buffer to hold raw NN outputs for all perp nodes (device memory when using GPU).
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
gkyl_bc_gksheath_choose_surrogate_kernels_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, 
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers);

GKYL_CU_D
static void
bc_gksheath_choose_surrogate_kernels(const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->vcut_calc = ser_sheath_vcut_calc_list[edge].dim_list[dim-2].kernels[poly_order-1];
      kers->infer = ser_sheath_infer_list[edge].dim_list[dim-2].kernels[poly_order-1];
      break;
    default:
      kers->vcut_calc = NULL;
      kers->infer = NULL;
      break;
  }
  assert(kers->vcut_calc);
  assert(kers->infer);
}

GKYL_CU_D
static sheath_infer_t
bc_gksheath_choose_infer_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  sheath_infer_t kern = NULL;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kern = ser_sheath_infer_list[edge].dim_list[dim-2].kernels[poly_order-1];
      break;
    default:
      kern = NULL;
      break;
  }
  assert(kern); // Infer kernel must be non-null if using surrogate.
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
