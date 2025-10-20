// Private header for fem_parproj updater.
#pragma once
#include <gkyl_fem_parproj_kernels.h>
#include <gkyl_basis.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

static long
gkyl_fem_parproj_global_num_nodes(const struct gkyl_basis *basis, bool isperiodic, int parnum_cells)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;
  enum gkyl_basis_type basis_type = basis->b_type;

  if (dim==1) {
    if (poly_order == 1) {
      return isperiodic? fem_parproj_num_nodes_global_1x_ser_p1_periodicx(parnum_cells)
                       : fem_parproj_num_nodes_global_1x_ser_p1_nonperiodicx(parnum_cells);
    } else if (poly_order == 2) {
      return isperiodic? fem_parproj_num_nodes_global_1x_ser_p2_periodicx(parnum_cells)
                       : fem_parproj_num_nodes_global_1x_ser_p2_nonperiodicx(parnum_cells);
    }
  } else if (dim==2) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_num_nodes_global_2x_ser_p1_periodicy(parnum_cells)
                         : fem_parproj_num_nodes_global_2x_ser_p1_nonperiodicy(parnum_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_num_nodes_global_2x_ser_p2_periodicy(parnum_cells)
                         : fem_parproj_num_nodes_global_2x_ser_p2_nonperiodicy(parnum_cells);
      }
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return isperiodic? fem_parproj_num_nodes_global_3x_ser_p1_periodicz(parnum_cells)
                         : fem_parproj_num_nodes_global_3x_ser_p1_nonperiodicz(parnum_cells);
      } else if (poly_order == 2) {
        return isperiodic? fem_parproj_num_nodes_global_3x_ser_p2_periodicz(parnum_cells)
                         : fem_parproj_num_nodes_global_3x_ser_p2_nonperiodicz(parnum_cells);
      }
    }
  }
  assert(false);  // Other dimensionalities or basis not supported.
  return -1;
}

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(int numCellsPar, int parIdx, long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[2]; } local2global_kern_loc_list;
typedef struct { local2global_kern_loc_list list[2]; } local2global_kern_bc_list;
typedef struct { local2global_kern_bc_list list[2]; } local2global_kern_list;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_list ser_loc2glob_list[] = {
  // 1x
  {.list={ 
          // periodicx
          {.list={
                  {fem_parproj_local_to_global_1x_ser_p1_inx_periodicx, fem_parproj_local_to_global_1x_ser_p1_upx_periodicx,}, 
                  {fem_parproj_local_to_global_1x_ser_p2_inx_periodicx, fem_parproj_local_to_global_1x_ser_p2_upx_periodicx,}, 
                 },
          },
          // nonperiodicx
          {.list={ 
                  {fem_parproj_local_to_global_1x_ser_p1_inx_nonperiodicx, fem_parproj_local_to_global_1x_ser_p1_upx_nonperiodicx,}, 
                  {fem_parproj_local_to_global_1x_ser_p2_inx_nonperiodicx, fem_parproj_local_to_global_1x_ser_p2_upx_nonperiodicx,}, 
                 },
          },
         }
  },
  // 2x
  {.list={ 
          // periodicy
          {.list={
                  {fem_parproj_local_to_global_2x_ser_p1_iny_periodicy, fem_parproj_local_to_global_2x_ser_p1_upy_periodicy,}, 
                  {fem_parproj_local_to_global_2x_ser_p2_iny_periodicy, fem_parproj_local_to_global_2x_ser_p2_upy_periodicy,}, 
                 },
          },
          // nonperiodicy
          {.list={ 
                  {fem_parproj_local_to_global_2x_ser_p1_iny_nonperiodicy, fem_parproj_local_to_global_2x_ser_p1_upy_nonperiodicy,}, 
                  {fem_parproj_local_to_global_2x_ser_p2_iny_nonperiodicy, fem_parproj_local_to_global_2x_ser_p2_upy_nonperiodicy,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // periodicz
          {.list={
                  {fem_parproj_local_to_global_3x_ser_p1_inz_periodicz, fem_parproj_local_to_global_3x_ser_p1_upz_periodicz,}, 
                  {fem_parproj_local_to_global_3x_ser_p2_inz_periodicz, fem_parproj_local_to_global_3x_ser_p2_upz_periodicz,}, 
                 },
          },
          // nonperiodicz
          {.list={ 
                  {fem_parproj_local_to_global_3x_ser_p1_inz_nonperiodicz, fem_parproj_local_to_global_3x_ser_p1_upz_nonperiodicz,}, 
                  {fem_parproj_local_to_global_3x_ser_p2_inz_nonperiodicz, fem_parproj_local_to_global_3x_ser_p2_upz_nonperiodicz,}, 
                 },
          },
         }
  },
};

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *weight, const long *globalIdxs, struct gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { lhsstencil_t kernels[3]; } lhsstencil_kern_loc_list;
typedef struct { lhsstencil_kern_loc_list list[2]; } lhsstencil_kern_bc_list;
typedef struct { lhsstencil_kern_bc_list list[2]; } lhsstencil_kern_list;

// Serendipity unweighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_noweight[] = {
  // 1x
  {.list={ 
          // nondirichletx
          {.list={
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_nondirichletx,}, 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_nondirichletx,}, 
                 },
          },
          // dirichletx
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_lox_dirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p1_upx_dirichletx,}, 
                  {fem_parproj_lhs_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_lox_dirichletx, fem_parproj_lhs_stencil_noweight_1x_ser_p2_upx_dirichletx,}, 
                 },
          },
         }
  },
  // 2x
  {.list={ 
          // nondirichlety
          {.list={
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_nondirichlety,}, 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_nondirichlety,}, 
                 },
          },
          // dirichlety
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_loy_dirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p1_upy_dirichlety,}, 
                  {fem_parproj_lhs_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_loy_dirichlety, fem_parproj_lhs_stencil_noweight_2x_ser_p2_upy_dirichlety,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // nondirichletz
          {.list={
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_nondirichletz,}, 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_nondirichletz,}, 
                 },
          },
          // dirichletz
          {.list={ 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_loz_dirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p1_upz_dirichletz,}, 
                  {fem_parproj_lhs_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_loz_dirichletz, fem_parproj_lhs_stencil_noweight_3x_ser_p2_upz_dirichletz,}, 
                 },
          },
         }
  },
};

// Serendipity weighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list_weighted[] = {
  // 1x
  {.list={ 
          // nondirichletx
          {.list={
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_nondirichletx,}, 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_nondirichletx,}, 
                 },
          },
          // dirichletx
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_lox_dirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p1_upx_dirichletx,}, 
                  {fem_parproj_lhs_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_lox_dirichletx, fem_parproj_lhs_stencil_weighted_1x_ser_p2_upx_dirichletx,}, 
                 },
          },
         }
  },
  // 2x
  {.list={ 
          // nondirichlety
          {.list={
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_nondirichlety,}, 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_nondirichlety,}, 
                 },
          },
          // dirichlety
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_loy_dirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p1_upy_dirichlety,}, 
                  {fem_parproj_lhs_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_loy_dirichlety, fem_parproj_lhs_stencil_weighted_2x_ser_p2_upy_dirichlety,}, 
                 },
          },
         }
  },
  // 3x
  {.list={ 
          // nondirichletz
          {.list={
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_nondirichletz,}, 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_nondirichletz,}, 
                 },
          },
          // dirichletz
          {.list={ 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_loz_dirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p1_upz_dirichletz,}, 
                  {fem_parproj_lhs_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_loz_dirichletz, fem_parproj_lhs_stencil_weighted_3x_ser_p2_upz_dirichletz,}, 
                 },
          },
         }
  },
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *weight, const double *rho, const double *phiBC, long nodeOff, const long *globalIdxs,
  double *bsrc);

typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_loc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_loc_list list[2]; } srcstencil_kern_bc_list;  // For use in kernel tables.
typedef struct { srcstencil_kern_bc_list list[2]; } srcstencil_kern_list;  // For use in kernel tables.

// Serendipity src kernels.
GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_noweight[] = {
  // 1x
  {.list={
          // nondirichletx
          {.list={
                  {fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_lox_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_upx_nondirichletx,},
                  {fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_lox_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_upx_nondirichletx,},
                 },
          },
          // dirichletx
          {.list={
                  {fem_parproj_src_stencil_noweight_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_lox_dirichletx, fem_parproj_src_stencil_noweight_1x_ser_p1_upx_dirichletx,},
                  {fem_parproj_src_stencil_noweight_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_lox_dirichletx, fem_parproj_src_stencil_noweight_1x_ser_p2_upx_dirichletx,},
                 },
          },
         }
  },
  // 2x
  {.list={
          // nondirichlety
          {.list={
                  {fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_loy_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_upy_nondirichlety,},
                  {fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_loy_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_upy_nondirichlety,},
                 },
          },
          // dirichlety
          {.list={
                  {fem_parproj_src_stencil_noweight_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_loy_dirichlety, fem_parproj_src_stencil_noweight_2x_ser_p1_upy_dirichlety,},
                  {fem_parproj_src_stencil_noweight_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_loy_dirichlety, fem_parproj_src_stencil_noweight_2x_ser_p2_upy_dirichlety,},
                 },
          },
         }
  },
  // 3x
  {.list={
          // nondirichletz
          {.list={
                  {fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_loz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_upz_nondirichletz,},
                  {fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_loz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_upz_nondirichletz,},
                 },
          },
          // dirichletz
          {.list={
                  {fem_parproj_src_stencil_noweight_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_loz_dirichletz, fem_parproj_src_stencil_noweight_3x_ser_p1_upz_dirichletz,},
                  {fem_parproj_src_stencil_noweight_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_loz_dirichletz, fem_parproj_src_stencil_noweight_3x_ser_p2_upz_dirichletz,},
                 },
          },
         }
  },
};

GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list_weighted[] = {
  // 1x
  {.list={
          // nondirichletx
          {.list={
                  {fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_lox_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_upx_nondirichletx,},
                  {fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_lox_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_upx_nondirichletx,},
                 },
          },
          // dirichletx
          {.list={
                  {fem_parproj_src_stencil_weighted_1x_ser_p1_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_lox_dirichletx, fem_parproj_src_stencil_weighted_1x_ser_p1_upx_dirichletx,},
                  {fem_parproj_src_stencil_weighted_1x_ser_p2_inx_nondirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_lox_dirichletx, fem_parproj_src_stencil_weighted_1x_ser_p2_upx_dirichletx,},
                 },
          },
         }
  },
  // 2x
  {.list={
          // nondirichlety
          {.list={
                  {fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_loy_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_upy_nondirichlety,},
                  {fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_loy_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_upy_nondirichlety,},
                 },
          },
          // dirichlety
          {.list={
                  {fem_parproj_src_stencil_weighted_2x_ser_p1_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_loy_dirichlety, fem_parproj_src_stencil_weighted_2x_ser_p1_upy_dirichlety,},
                  {fem_parproj_src_stencil_weighted_2x_ser_p2_iny_nondirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_loy_dirichlety, fem_parproj_src_stencil_weighted_2x_ser_p2_upy_dirichlety,},
                 },
          },
         }
  },
  // 3x
  {.list={
          // nondirichletz
          {.list={
                  {fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_loz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_upz_nondirichletz,},
                  {fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_loz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_upz_nondirichletz,},
                 },
          },
          // dirichletz
          {.list={
                  {fem_parproj_src_stencil_weighted_3x_ser_p1_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_loz_dirichletz, fem_parproj_src_stencil_weighted_3x_ser_p1_upz_dirichletz,},
                  {fem_parproj_src_stencil_weighted_3x_ser_p2_inz_nondirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_loz_dirichletz, fem_parproj_src_stencil_weighted_3x_ser_p2_upz_dirichletz,},
                 },
          },
         }
  },
};

// Function pointer type for kernels that convert the solution from nodal to
// modal.
typedef void (*solstencil_t)(const double *sol_nodal_global, long nodeOff,
  const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;  // For use in kernel tables.

// Serendipity sol kernels.
GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { fem_parproj_sol_stencil_1x_ser_p1, fem_parproj_sol_stencil_1x_ser_p2 },
  { fem_parproj_sol_stencil_2x_ser_p1, fem_parproj_sol_stencil_2x_ser_p2 },
  { fem_parproj_sol_stencil_3x_ser_p1, fem_parproj_sol_stencil_3x_ser_p2 }
};

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_parproj_kernels {
  local2global_t l2g[2];  // Pointer to local-to-global kernel.

  lhsstencil_t lhsker[3];  // Weighted LHS kernel.

  srcstencil_t srcker[3];  // RHS source kernel.

  solstencil_t solker;  // Kernel that takes the solution and converts it to modal.
};

struct gkyl_fem_parproj {
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // parallel (z) direction.
  int parnum_cells; // number of cells in parallel (z) direction.
  bool isperiodic; // =true if parallel direction is periodic.
  bool isdirichlet; // =true if parallel direction has periodic BCs.
  bool has_weight_rhs; // Whether there's a weight on the RHS.
  struct gkyl_array *weight_rhs; // The RHS weight.

  const struct gkyl_range *solve_range;
  struct gkyl_range perp_range2d; // 2D range of perpendicular cells.
  struct gkyl_range par_range1d; // 1D range of parallel cells.
  struct gkyl_range_iter perp_iter2d;
  struct gkyl_range_iter par_iter1d;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_array *brhs;

  struct gkyl_superlu_prob* prob;
#ifdef GKYL_HAVE_CUDA
  struct gkyl_culinsolver_prob* prob_cu;
  struct gkyl_array *brhs_cu;
#endif

  long *globalidx;

  struct gkyl_fem_parproj_kernels *kernels;
  struct gkyl_fem_parproj_kernels *kernels_cu;
  bool use_gpu;  
};

// "Choose Kernel" based on polyorder, stencil location and BCs.
#define CK(lst,dim,bc,poly_order,loc) lst[dim-1].list[bc].list[poly_order-1].kernels[loc]

void
fem_parproj_choose_kernels_cu(const struct gkyl_basis* basis, bool isweighted,
  bool isperiodic, bool isdirichlet, struct gkyl_fem_parproj_kernels *kers);

GKYL_CU_D
static void
fem_parproj_choose_local2global_kernel(const struct gkyl_basis *basis, bool isperiodic, local2global_t *l2gout)
{
  int bckey[1] = {-1};
  bckey[0] = isperiodic? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++)
        l2gout[k] = CK(ser_loc2glob_list, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void 
fem_parproj_choose_lhs_kernel(const struct gkyl_basis *basis, bool isdirichlet, bool isweighted, lhsstencil_t *lhsout)
{
  int bckey[1] = {-1};
  bckey[0] = isdirichlet? 1 : 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<3; k++)
        lhsout[k] = isweighted? CK(ser_lhsstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                              : CK(ser_lhsstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k); 
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_parproj_choose_srcstencil_kernel(const struct gkyl_basis *basis, bool isdirichlet, bool isweighted, srcstencil_t *srcout)
{
  int bckey[1] = {-1};
  bckey[0] = isdirichlet? 1 : 0;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<3; k++)
        srcout[k] = isweighted? CK(ser_srcstencil_list_weighted, basis->ndim, bckey[0], basis->poly_order, k)
                              : CK(ser_srcstencil_list_noweight, basis->ndim, bckey[0], basis->poly_order, k);
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static solstencil_t
fem_parproj_choose_solstencil_kernel(const struct gkyl_basis *basis)
{
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_solstencil_list[basis->ndim-1].kernels[basis->poly_order-1];
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_DH
static inline int idx_to_inloup_ker(int num_cells, int idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  if (idx == 1)
    iout = 1;
  else if (idx == num_cells)
    iout = 2;
  return iout;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Assign the right-side vector with the discontinuous (DG) source field
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_set_rhs_cu(struct gkyl_fem_parproj *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve_cu(struct gkyl_fem_parproj* up, struct gkyl_array *phiout);
#endif