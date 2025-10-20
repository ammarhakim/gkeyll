#pragma once

// Private header for bc_basic_gyrokinetic updater, not for direct use in user code.

#include <gkyl_bc_basic_gyrokinetic.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>
#include <math.h>

// Primary struct in this updater.
struct gkyl_bc_basic_gyrokinetic {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  enum gkyl_gyrokinetic_bc_type bctype;
  const struct gkyl_range *skin_r, *ghost_r;
  struct gkyl_array_copy_func *array_copy_func;
  bool use_gpu;
};

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.

 * @param dir Direction in which to apply BC .
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param cdim Number of configuration space dimensions.
 * @param bctype Type of BC .
 * @param basis Basis in which to expand coefficients in array we apply BC to.
 * @param num_comp Number of components (DOFs) within a cell.
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods.
 */
struct gkyl_array_copy_func* gkyl_bc_basic_gyrokinetic_create_arr_copy_func_cu(int dir,
  enum gkyl_edge_loc edge, int cdim, enum gkyl_gyrokinetic_bc_type bctype,
  const struct gkyl_basis *basis, int num_comp);

#endif

// context for use in BCs
struct dg_bc_ctx {
  int dir; // direction for BCs.
  enum gkyl_edge_loc edge; // lower/upper edge boundary.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
copy_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  // Copy skin cell into ghost cell
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int num_comp = mc->ncomp;
  for (int c=0; c<num_comp; ++c) out[c] = inp[c];
}

GKYL_CU_D
static void
species_absorb_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  // Set ghost cell to zero.
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int num_comp = mc->ncomp;
  for (int c=0; c<num_comp; ++c) out[c] = 0.0;
}

GKYL_CU_D
static void
species_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  // Fill the ghost cell with a reflection of the distribution function
  // in the skin cell so that particle are reflected back into the domain.
  // (only works for the parallel direction, whose advection is given by the
  // vpar coordinate).
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;

  mc->basis->flip_odd_sign(dir, inp, out);
  mc->basis->flip_odd_sign(cdim, out, out);
}

GKYL_CU_D
static void
conf_boundary_value_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  // Fill the ghost cell with the skin cell evaluated at the boundary,
  // so it has no variation in the direction of the BC.
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int cdim = mc->cdim;
  enum gkyl_edge_loc edge = mc->edge;

  if (cdim == 1) {
    if (edge == GKYL_LOWER_EDGE) {
      out[0] = inp[0]-sqrt(3.0)*inp[1];
    }
    else {
      out[0] = inp[0]+sqrt(3.0)*inp[1];
    }
    out[1] = 0.0;
  }
  else if (cdim == 2) {
    if (dir == 0) {
      if (edge == GKYL_LOWER_EDGE) {
        out[0] = inp[0]-sqrt(3.0)*inp[1];
        out[2] = inp[2]-sqrt(3.0)*inp[3];
      }
      else {
        out[0] = inp[0]+sqrt(3.0)*inp[1];
        out[2] = inp[2]+sqrt(3.0)*inp[3];
      }
      out[1] = 0.0;
      out[3] = 0.0;
    }
    else if (dir == 1) {
      if (edge == GKYL_LOWER_EDGE) {
        out[0] = inp[0]-sqrt(3.0)*inp[2];
        out[1] = inp[1]-sqrt(3.0)*inp[3];
      }
      else {
        out[0] = inp[0]+sqrt(3.0)*inp[2];
        out[1] = inp[1]+sqrt(3.0)*inp[3];
      }
      out[2] = 0.0;
      out[3] = 0.0;
    }
  }
  else if (cdim == 3) {
    if (dir == 0) {
      if (edge == GKYL_LOWER_EDGE) {
        out[0] = inp[0]-sqrt(3.0)*inp[1];
        out[2] = inp[2]-sqrt(3.0)*inp[4];
        out[3] = inp[3]-sqrt(3.0)*inp[5];
        out[6] = inp[6]-sqrt(3.0)*inp[7];
      }
      else {
        out[0] = inp[0]+sqrt(3.0)*inp[1];
        out[2] = inp[2]+sqrt(3.0)*inp[4];
        out[3] = inp[3]+sqrt(3.0)*inp[5];
        out[6] = inp[6]+sqrt(3.0)*inp[7];
      }
      out[1] = 0.0;
      out[4] = 0.0;
      out[5] = 0.0;
      out[7] = 0.0;
    }
    else if (dir == 1) {
      if (edge == GKYL_LOWER_EDGE) {
        out[0] = inp[0]-sqrt(3.0)*inp[2];
        out[1] = inp[1]-sqrt(3.0)*inp[4];
        out[3] = inp[3]-sqrt(3.0)*inp[6];
        out[5] = inp[5]-sqrt(3.0)*inp[7];
      }
      else {
        out[0] = inp[0]+sqrt(3.0)*inp[2];
        out[1] = inp[1]+sqrt(3.0)*inp[4];
        out[3] = inp[3]+sqrt(3.0)*inp[6];
        out[5] = inp[5]+sqrt(3.0)*inp[7];
      }
      out[2] = 0.0;
      out[4] = 0.0;
      out[6] = 0.0;
      out[7] = 0.0;
    }
    else if (dir == 2) {
      if (edge == GKYL_LOWER_EDGE) {
        out[0] = inp[0]-sqrt(3.0)*inp[3];
        out[1] = inp[1]-sqrt(3.0)*inp[5];
        out[2] = inp[2]-sqrt(3.0)*inp[6];
        out[4] = inp[4]-sqrt(3.0)*inp[7];
      }
      else {
        out[0] = inp[0]+sqrt(3.0)*inp[3];
        out[1] = inp[1]+sqrt(3.0)*inp[5];
        out[2] = inp[2]+sqrt(3.0)*inp[6];
        out[4] = inp[4]+sqrt(3.0)*inp[7];
      }
      out[3] = 0.0;
      out[5] = 0.0;
      out[6] = 0.0;
      out[7] = 0.0;
    }
  }
}
