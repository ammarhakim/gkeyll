// Private header for fem_parproj_couplex updater.
#pragma once
#include <gkyl_fem_parproj_couplex.h>
#include <gkyl_fem_parproj_couplex_kernels.h>
#include <gkyl_util.h>
#include <assert.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

static long
gkyl_fem_parproj_couplex_global_num_nodes(const struct gkyl_basis *basis, bool isperiodic, const int *num_cells)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;
  enum gkyl_basis_type basis_type = basis->b_type;

  if (dim==1) {
    if (poly_order == 1) {
      return isperiodic? fem_parproj_couplex_num_nodes_global_1x_ser_p1_periodicx(num_cells)
                       : fem_parproj_couplex_num_nodes_global_1x_ser_p1_nonperiodicx(num_cells);
    } else if (poly_order == 2) {
      return isperiodic? fem_parproj_couplex_num_nodes_global_1x_ser_p2_periodicx(num_cells)
                       : fem_parproj_couplex_num_nodes_global_1x_ser_p2_nonperiodicx(num_cells);
    }
  } else if (dim==2) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_couplex_num_nodes_global_2x_ser_p1_periodicy(num_cells)
                         : fem_parproj_couplex_num_nodes_global_2x_ser_p1_nonperiodicy(num_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_couplex_num_nodes_global_2x_ser_p2_periodicy(num_cells)
                         : fem_parproj_couplex_num_nodes_global_2x_ser_p2_nonperiodicy(num_cells);
      }
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_couplex_num_nodes_global_3x_ser_p1_periodicz(num_cells)
                         : fem_parproj_couplex_num_nodes_global_3x_ser_p1_nonperiodicz(num_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_couplex_num_nodes_global_3x_ser_p2_periodicz(num_cells)
                         : fem_parproj_couplex_num_nodes_global_3x_ser_p2_nonperiodicz(num_cells);
      }
    }
  }
  assert(false);  // Other dimensionalities or basis not supported.
  return -1;
}

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int *num_cells, const int *idx_curr, long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[4]; } local2global_kern_loc_list;
typedef struct { local2global_kern_loc_list list[2]; } local2global_kern_bc_list;
typedef struct { local2global_kern_bc_list list[2]; } local2global_kern_list;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_list ser_loc2glob_list[] = {
  // 1x
  {.list =
    { 
      // periodicx
      {.list = 
        {
          {fem_parproj_couplex_local_to_global_1x_ser_p1_inx_periodicx, fem_parproj_couplex_local_to_global_1x_ser_p1_upx_periodicx,NULL,NULL,}, 
          {fem_parproj_couplex_local_to_global_1x_ser_p2_inx_periodicx, fem_parproj_couplex_local_to_global_1x_ser_p2_upx_periodicx,NULL,NULL,}, 
        },
      },
      // nonperiodicx
      {.list = 
        { 
          {fem_parproj_couplex_local_to_global_1x_ser_p1_inx_nonperiodicx, fem_parproj_couplex_local_to_global_1x_ser_p1_upx_nonperiodicx,NULL,NULL,}, 
          {fem_parproj_couplex_local_to_global_1x_ser_p2_inx_nonperiodicx, fem_parproj_couplex_local_to_global_1x_ser_p2_upx_nonperiodicx,NULL,NULL,}, 
        },
      },
    }
  },
  // 2x
  {.list = 
    { 
      // periodicy
      {.list = 
        {
          {fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_periodicy,}, 
          {fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_periodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_periodicy,}, 

        },
      },
      // nonperiodicy
      {.list = 
        { 
          {fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_iny_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_iny_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_inx_nonperiodicx_upy_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p1_upx_nonperiodicx_upy_nonperiodicy,}, 
          {fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_iny_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_iny_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_inx_nonperiodicx_upy_nonperiodicy, fem_parproj_couplex_local_to_global_2x_ser_p2_upx_nonperiodicx_upy_nonperiodicy,}, 
        },
      },
    }
  },
  // 3x
  {.list = 
    { 
      // periodicz
      {.list = 
        {
          {fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_periodicz,}, 
          {fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_periodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_periodicz,}, 
        },
      },
      // nonperiodicz
      {.list = 
        { 
          {fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_inz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_inz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_inx_nonperiodicx_upz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p1_upx_nonperiodicx_upz_nonperiodicz,}, 
          {fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_inz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_inz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_inx_nonperiodicx_upz_nonperiodicz, fem_parproj_couplex_local_to_global_3x_ser_p2_upx_nonperiodicx_upz_nonperiodicz,}, 
        },
      },
    }
  },
};

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { lhsstencil_t kernels[9]; } lhsstencil_kern_loc_list;
typedef struct { lhsstencil_kern_loc_list list[2]; } lhsstencil_kern_bc_list;
typedef struct { lhsstencil_kern_bc_list list[2]; } lhsstencil_kern_list;

// Serendipity unweighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_noweight[] = {
  // 1x
  {.list = 
    { 
      // nondirichletx
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
        },
      },
      // dirichletx
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_lox_dirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p1_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_lox_dirichletx, fem_parproj_couplex_lhs_stencil_noweight_1x_ser_p2_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
        },
      },
    }
  },
  // 2x
  {.list = 
    { 
      // nondirichlety
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_nondirichlety,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_nondirichlety,}, 
        },
      },
      // dirichlety
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_dirichlety,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_dirichlety,}, 
        },
      },
    }
  },
  // 3x
  {.list = 
    { 
      // nondirichletz
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_nondirichletz,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_nondirichletz,}, 
        },
      },
      // dirichletz
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_dirichletz,}, 
          {fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_dirichletz,}, 
        },
      },
    }
  },
};

// Serendipity weighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_weighted[] = {
  // 1x
  {.list = 
    { 
      // nondirichletx
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
        },
      },
      // dirichletx
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_lox_dirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p1_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_lox_dirichletx, fem_parproj_couplex_lhs_stencil_weighted_1x_ser_p2_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,}, 
        },
      },
    }
  },
  // 2x
  {.list = 
    { 
      // nondirichlety
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_nondirichlety,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_nondirichlety,}, 
        },
      },
      // dirichlety
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_dirichlety,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_lhs_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_dirichlety,}, 
        },
      },
    }
  },
  // 3x
  {.list = 
    { 
      // nondirichletz
      {.list = 
        {
          {fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_nondirichletz,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_nondirichletz,}, 
        },
      },
      // dirichletz
      {.list = 
        { 
          {fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_dirichletz,}, 
          {fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_lhs_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_dirichletz,}, 
        },
      },
    }
  },
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs,
  double *bsrc);

typedef struct { srcstencil_t kernels[9]; } srcstencil_kern_loc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_loc_list list[2]; } srcstencil_kern_bc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_bc_list list[2]; } srcstencil_kern_list;  // For use in kernel tables.

// Serendipity unweighted src kernels.
GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_noweight[] = {
  // 1x
  {.list = 
    {
      // nondirichletx
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
          {fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
        },
      },
      // dirichletx
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_lox_dirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p1_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
          {fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_lox_dirichletx, fem_parproj_couplex_src_stencil_noweight_1x_ser_p2_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
        },
      },
    }
  },
  // 2x
  {.list = 
    {
      // nondirichlety
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_nondirichlety,},
          {fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_nondirichlety,},
        },
      },
      // dirichlety
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p1_upx_nondirichletx_upy_dirichlety,},
          {fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_noweight_2x_ser_p2_upx_nondirichletx_upy_dirichlety,},
        },
      },
    }
  },
  // 3x
  {.list = 
    {
      // nondirichletz
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_nondirichletz,},
          {fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_nondirichletz,},
        },
      },
      // dirichletz
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p1_upx_nondirichletx_upz_dirichletz,},
          {fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_noweight_3x_ser_p2_upx_nondirichletx_upz_dirichletz,},
        },
      },
    }
  },
};

// Serendipity weighted src kernels.
GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_weighted[] = {
  // 1x
  {.list = 
    {
      // nondirichletx
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
          {fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_upx_nondirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
        },
      },
      // dirichletx
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_lox_dirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p1_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
          {fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_lox_dirichletx, fem_parproj_couplex_src_stencil_weighted_1x_ser_p2_upx_dirichletx,NULL,NULL,NULL,NULL,NULL,NULL,},
        },
      },
    }
  },
  // 2x
  {.list = 
    {
      // nondirichlety
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_nondirichlety,},
          {fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_nondirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_nondirichlety,},
        },
      },
      // dirichlety
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p1_upx_nondirichletx_upy_dirichlety,},
          {fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_iny_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_loy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_inx_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_lox_nondirichletx_upy_dirichlety, fem_parproj_couplex_src_stencil_weighted_2x_ser_p2_upx_nondirichletx_upy_dirichlety,},
        },
      },
    }
  },
  // 3x
  {.list = 
    {
      // nondirichletz
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_nondirichletz,},
          {fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_nondirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_nondirichletz,},
        },
      },
      // dirichletz
      {.list = 
        {
          {fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p1_upx_nondirichletx_upz_dirichletz,},
          {fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_inz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_loz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_inx_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_lox_nondirichletx_upz_dirichletz, fem_parproj_couplex_src_stencil_weighted_3x_ser_p2_upx_nondirichletx_upz_dirichletz,},
        },
      },
    }
  },
};

// Function pointer type for kernels that convert the solution from nodal to modal.
typedef void (*solstencil_t)(const double *sol_nodal_global, long nodeOff,
  const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;  // For use in kernel tables.

// Serendipity sol kernels.
GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { fem_parproj_couplex_sol_stencil_1x_ser_p1, fem_parproj_couplex_sol_stencil_1x_ser_p2 },
  { fem_parproj_couplex_sol_stencil_2x_ser_p1, fem_parproj_couplex_sol_stencil_2x_ser_p2 },
  { fem_parproj_couplex_sol_stencil_3x_ser_p1, fem_parproj_couplex_sol_stencil_3x_ser_p2 }
};

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_parproj_couplex_kernels {
  local2global_t l2g[4]; // Pointer to local-to-global kernel.
  lhsstencil_t lhsker[9]; // Weighted LHS kernel.
  srcstencil_t srcker[9]; // RHS source kernel.
  solstencil_t solker; // Kernel that takes the solution and converts it to modal.
};

struct gkyl_fem_parproj_couplex {
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // parallel (z) direction.
  int perpdir;   // perpendicular (y) direction.
  int num_cells[GKYL_MAX_CDIM]; // Number of cells in each direction.
  bool isperiodic; // =true if parallel direction is periodic.
  bool isdirichlet; // =true if parallel direction has periodic BCs.
  bool has_weight_rhs; // Whether there's a weight on the RHS.
  struct gkyl_array *weight_rhs; // The RHS weight.

  const struct gkyl_range *solve_range;
  struct gkyl_range perp_range; // 1D range of perpendicular cells (y).
  struct gkyl_range fem_range; // 2D range of cells spanned by FEM basis (x, z).
  struct gkyl_range_iter perp_iter;
  struct gkyl_range_iter fem_iter;

  int numnodes_local;
  long numnodes_global;

  long *globalidx;

  struct gkyl_superlu_prob* prob;
  struct gkyl_array *brhs;
  struct gkyl_fem_parproj_couplex_kernels *kernels;
#ifdef GKYL_HAVE_CUDA
  struct gkyl_culinsolver_prob* prob_cu;
  struct gkyl_array *brhs_cu;
  struct gkyl_fem_parproj_couplex_kernels *kernels_cu;
#endif

  bool use_gpu;  
};

GKYL_CU_DH
static inline int idx_to_inup_ker(const int dim, const int *range_upper, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between upper cells and
  // elsewhere.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == range_upper[d]) iout += (int)(pow(2,d)+0.5);
  }
  return iout;
}

GKYL_CU_DH
static inline int idx_to_inloup_ker(int dim, const int *range_lower, const int *range_upper, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == range_lower[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[d] == range_upper[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}

// "Choose Kernel" based on polyorder, stencil location and BCs.
#define CK(lst,dim,bc,poly_order,loc) lst[dim-1].list[bc].list[poly_order-1].kernels[loc]

GKYL_CU_D
static void
fem_parproj_couplex_choose_kernels(const struct gkyl_basis* basis, bool has_weight_lhs, bool has_weight_rhs,
  bool isperiodic, bool isdirichlet, struct gkyl_fem_parproj_couplex_kernels *kers)
{
  // Local-to-global kernel.
  int bckey[1] = {0};
  bckey[0] = isperiodic? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<4; k++)
        kers->l2g[k] = CK(ser_loc2glob_list, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }

  // LHS kernel.
  bckey[0] = isdirichlet? 1 : 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<9; k++)
        kers->lhsker[k] = has_weight_lhs? CK(ser_lhsstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                                        : CK(ser_lhsstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k); 
      break;
    default:
      assert(false);
      break;
  }

  // RHS source kernel.
  bckey[0] = isdirichlet? 1 : 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<9; k++)
        kers->srcker[k] = has_weight_rhs? CK(ser_srcstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                                        : CK(ser_srcstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }

  // Solution kernel.
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->solker = ser_solstencil_list[basis->ndim-1].kernels[basis->poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

#ifdef GKYL_HAVE_CUDA
/**
 * Assign the right-side vector with the discontinuous (DG) source field
 * on the NVIDIA GPU.
 */
void gkyl_fem_parproj_couplex_set_rhs_cu(struct gkyl_fem_parproj_couplex *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem
 * on the NVIDIA GPU.
 */
void gkyl_fem_parproj_couplex_solve_cu(struct gkyl_fem_parproj_couplex* up, struct gkyl_array *phiout);
#endif
