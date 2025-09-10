// Private header for fem_poisson_perp updater.
#pragma once
#include <gkyl_fem_poisson_perp_kernels.h>
#include <gkyl_basis.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

#ifndef GKYL_IPOW
# define GKYL_IPOW(a,e) (int)(pow(a,e)+0.5)
#endif

#define PERP_DIM_MAX 2

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int *numCells, const int *idx,
  long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[2]; } local2global_kern_loc_list_2x;
typedef struct { local2global_kern_loc_list_2x list[3]; } local2global_kern_bcx_list_2x;

typedef struct { local2global_t kernels[4]; } local2global_kern_loc_list_3x;
typedef struct { local2global_kern_loc_list_3x list[3]; } local2global_kern_bcy_list_3x;
typedef struct { local2global_kern_bcy_list_3x list[2]; } local2global_kern_bcx_list_3x;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_bcx_list_2x ser_loc2glob_list_2x[] = {
  // periodicx
  { .list =
   {
     {NULL, NULL},
     {fem_poisson_perp_local_to_global_2x_ser_p1_inx_periodicx, fem_poisson_perp_local_to_global_2x_ser_p1_upx_periodicx},
     {NULL, NULL},
   },
  },
  // nonperiodicx
  { .list =
    { 
      {NULL, NULL},
      {fem_poisson_perp_local_to_global_2x_ser_p1_inx_nonperiodicx, fem_poisson_perp_local_to_global_2x_ser_p1_upx_nonperiodicx},
      {NULL, NULL},
    },
  },
};

GKYL_CU_D
static const local2global_kern_bcx_list_3x ser_loc2glob_list_3x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_periodicy},
               {NULL, NULL, NULL, NULL},
              },
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_iny_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_iny_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_inx_periodicx_upy_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_periodicx_upy_nonperiodicy},
               {NULL, NULL, NULL, NULL},
              }
    }}
  },
  // nonperiodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_periodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_periodicy},
               {NULL, NULL, NULL, NULL},
              }
    },
    // nonperiodicy
    { .list = {{NULL, NULL, NULL, NULL},
               {fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_iny_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_iny_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_inx_nonperiodicx_upy_nonperiodicy, fem_poisson_perp_local_to_global_3x_ser_p1_upx_nonperiodicx_upy_nonperiodicy},
               {NULL, NULL, NULL, NULL},
              }
    }}
  }
};

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *epsilon, const double *kSq, const double *dx, const double *bcVals,
  const long *globalIdxs, gkyl_mat_triples *tri);

// For use in kernel tables.
typedef struct { lhsstencil_t kernels[3]; } lhsstencil_kern_loc_list_2x;
typedef struct { lhsstencil_kern_loc_list_2x list[3]; } lhsstencil_kern_bcx_list_2x;

typedef struct { lhsstencil_t kernels[9]; } lhsstencil_kern_loc_list_3x;
typedef struct { lhsstencil_kern_loc_list_3x list[3]; } lhsstencil_kern_bcy_list_3x;
typedef struct { lhsstencil_kern_bcy_list_3x list[9]; } lhsstencil_kern_bcx_list_3x;

// Serendipity lhs kernels.
static const lhsstencil_kern_bcx_list_2x ser_lhsstencil_list_2x[] = {
  // periodicx
  { .list =
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_periodicx, fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_periodicx},
      {NULL, NULL, NULL},
    },
  },
  // dirichletx-dirichletx
  { .list =
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_dirichletx, fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_dirichletx},
      {NULL, NULL, NULL},
    },
  },
  // dirichletx-neumannx
  { .list = 
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_dirichletx, fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_neumannx},
      {NULL, NULL, NULL},
    },
  },
  // neumannx-dirichletx
  { .list =
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_lhs_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_lhs_stencil_2x_ser_p1_lox_neumannx, fem_poisson_perp_lhs_stencil_2x_ser_p1_upx_dirichletx},
      {NULL, NULL, NULL},
    }
  },
};

static const lhsstencil_kern_bcx_list_3x ser_lhsstencil_list_3x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_perp_lhs_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *epsilon, const double *dx, const double *rho, const double *bcVals, long perpOff, const long *globalIdxs,
  double *bsrc);

// For use in kernel tables.
typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_loc_list_2x;
typedef struct { srcstencil_kern_loc_list_2x list[3]; } srcstencil_kern_bcx_list_2x;

typedef struct { srcstencil_t kernels[9]; } srcstencil_kern_loc_list_3x;
typedef struct { srcstencil_kern_loc_list_3x list[3]; } srcstencil_kern_bcy_list_3x;
typedef struct { srcstencil_kern_bcy_list_3x list[9]; } srcstencil_kern_bcx_list_3x;

// Serendipity src kernels.
GKYL_CU_D
static const srcstencil_kern_bcx_list_2x ser_srcstencil_list_2x[] = {
  // periodicx
  { .list =
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_src_stencil_2x_ser_p1_lox_periodicx, fem_poisson_perp_src_stencil_2x_ser_p1_upx_periodicx},
      {NULL, NULL, NULL},
    },
  },
  // dirichletx-dirichletx
  { .list = 
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_src_stencil_2x_ser_p1_lox_dirichletx, fem_poisson_perp_src_stencil_2x_ser_p1_upx_dirichletx},
      {NULL, NULL, NULL},
    },
  },
  // dirichletx-neumannx
  { .list = 
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_src_stencil_2x_ser_p1_lox_dirichletx, fem_poisson_perp_src_stencil_2x_ser_p1_upx_neumannx},
      {NULL, NULL, NULL},
    },
  },
  // neumannx-dirichletx
  { .list = 
    {
      {NULL, NULL, NULL},
      {fem_poisson_perp_src_stencil_2x_ser_p1_inx_periodicx, fem_poisson_perp_src_stencil_2x_ser_p1_lox_neumannx, fem_poisson_perp_src_stencil_2x_ser_p1_upx_dirichletx},
      {NULL, NULL, NULL},
    },
  },
};

GKYL_CU_D
static const srcstencil_kern_bcx_list_3x ser_srcstencil_list_3x[] = {
  // periodicx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_periodicx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // dirichletx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // dirichletx-neumannx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_dirichletx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_neumannx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
  // neumannx-dirichletx
  { .list = {
    // periodicy
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_periodicy},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // dirichlety-neumanny
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_neumanny},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    // neumanny-dirichlety
    { .list = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
               {fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_iny_periodicy, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_inx_periodicx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_lox_neumannx_upy_dirichlety, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_loy_neumanny, fem_poisson_perp_src_stencil_3x_ser_p1_upx_dirichletx_upy_dirichlety},
               {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
              },
    },
    }
  },
};

// Function pointer type for sol kernels.
typedef void (*solstencil_t)(const double *sol_nodal_global, long perpOff, const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;

GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { NULL, NULL, NULL },
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, fem_poisson_perp_sol_stencil_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, fem_poisson_perp_sol_stencil_3x_ser_p1, NULL }, // 2
};

// "Choose Kernel" based on polyorder, stencil location and BCs.
#define CK2x(lst,poly_order,loc,bcx) lst[bcx].list[poly_order].kernels[loc]
#define CK3x(lst,poly_order,loc,bcx,bcy) lst[bcx].list[bcy].list[poly_order].kernels[loc]

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_poisson_perp_kernels { 
  // Pointer to local-to-global kernels. 2^3, 2 (interior and upper) in each direction.
  local2global_t l2g[8];

  // LHS kernels (one for each position in the domain, up to 3^3).
  lhsstencil_t lhsker[27];

  // RHS source kernels (one for each position in the domain, up to 3^3).
  srcstencil_t srcker[27];

  solstencil_t solker;
};

// Updater type
struct gkyl_fem_poisson_perp {
  void *ctx; // Evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // Grid's number of dimensions.
  int ndim_perp; // Grid's number of perpendicular dimensions.
  int num_basis; // Number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // Parallel (z) direction.
  struct gkyl_basis basis;
  int num_cells[GKYL_MAX_CDIM];
  double dx[GKYL_MAX_CDIM];
#ifdef GKYL_HAVE_CUDA
  double *dx_cu;
#endif
  bool isdirperiodic[GKYL_MAX_DIM]; // =true if direction is periodic.

  struct gkyl_array *epsilon; // Permittivity.
  bool ishelmholtz; // If solving Helmholtz equation (kSq is not zero/NULL).

  bool isdomperiodic; // =true if all directions are periodic.
  struct gkyl_array *rhs_cellavg;
  double *rhs_avg, mavgfac;
  double *rhs_avg_cu;

  double bcvals[PERP_DIM_MAX*2*3]; // BC values, bc[0]*phi+bc[1]*d(phi)/dx=phi[3] at each boundary.
  double *bcvals_cu; // BC values, bc[0]*phi+bc[1]*d(phi)/dx=phi[3] at each boundary.

  const struct gkyl_range *solve_range;
  struct gkyl_range_iter solve_iter;

  struct gkyl_range *perp_range;
  struct gkyl_range par_range;
  struct gkyl_range perp_range2d, par_range1d;
  struct gkyl_range_iter perp_iter2d, par_iter1d;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_superlu_prob *prob;
  struct gkyl_array *brhs;

#ifdef GKYL_HAVE_CUDA
  struct gkyl_culinsolver_prob *prob_cu;
  struct gkyl_array *brhs_cu;
#endif

  long *globalidx;

  struct gkyl_fem_poisson_perp_kernels *kernels;
  struct gkyl_fem_poisson_perp_kernels *kernels_cu;
  bool use_gpu;
};

void
fem_poisson_perp_choose_kernels_cu(const struct gkyl_basis* basis, const struct gkyl_poisson_bc* bcs, const bool *isdirperiodic, struct gkyl_fem_poisson_perp_kernels *kers);

static long
gkyl_fem_poisson_perp_global_num_nodes(int ndim, int poly_order, int basis_type, const int *num_cells, bool *isdirperiodic)
{
  if (ndim == 2) {
    if (poly_order == 1) {
      if (isdirperiodic[0]) {
        return fem_poisson_perp_num_nodes_global_2x_ser_p1_periodicx(num_cells);
      } else if (!isdirperiodic[0]) {
        return fem_poisson_perp_num_nodes_global_2x_ser_p1_nonperiodicx(num_cells);
      }
    }
  }
  else {
    if (poly_order == 1) {
      if (isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_periodicy(num_cells);
      } else if (!isdirperiodic[0] && isdirperiodic[1]) {
        return fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_periodicy(num_cells);
      } else if (isdirperiodic[0] && !isdirperiodic[1]) {
        return fem_poisson_perp_num_nodes_global_3x_ser_p1_periodicx_nonperiodicy(num_cells);
      } else {
        return fem_poisson_perp_num_nodes_global_3x_ser_p1_nonperiodicx_nonperiodicy(num_cells);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
  return -1;
}

GKYL_CU_D
static void
fem_poisson_perp_choose_local2global_kernels(const struct gkyl_basis* basis, const bool *isdirperiodic, local2global_t *l2gout)
{
  int poly_order = basis->poly_order;
  int ndim = basis->ndim;
  int ndim_perp = ndim-1;

  int bckey[GKYL_MAX_CDIM] = {-1};
  for (int d=0; d<ndim_perp; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<GKYL_IPOW(2,ndim_perp); k++)
        l2gout[k] = ndim == 2? CK2x(ser_loc2glob_list_2x, poly_order, k, bckey[0])
                             : CK3x(ser_loc2glob_list_3x, poly_order, k, bckey[0], bckey[1]);
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_poisson_perp_choose_lhs_kernels(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs, lhsstencil_t *lhsout)
{
  int poly_order = basis->poly_order;
  int ndim = basis->ndim;
  int ndim_perp = ndim-1;

  int bckey[GKYL_MAX_CDIM] = {-1, -1, -1};
  for (int d=0; d<ndim_perp; d++) {
         if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC  && bcs->up_type[d]==GKYL_POISSON_PERIODIC ) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN  ) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN   && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else { assert(false); }
  };

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<GKYL_IPOW(3,ndim_perp); k++)
        lhsout[k] = ndim == 2? CK2x(ser_lhsstencil_list_2x, poly_order, k, bckey[0])
                             : CK3x(ser_lhsstencil_list_3x, poly_order, k, bckey[0], bckey[1]);
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
fem_poisson_perp_choose_src_kernels(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs, srcstencil_t *srcout)
{
  int poly_order = basis->poly_order;
  int ndim = basis->ndim;
  int ndim_perp = ndim-1;

  int bckey[GKYL_MAX_CDIM] = {-1, -1, -1};
  for (int d=0; d<ndim_perp; d++) {
         if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC  && bcs->up_type[d]==GKYL_POISSON_PERIODIC ) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN  ) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN   && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else { assert(false); }
  };

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<GKYL_IPOW(3,ndim_perp); k++)
        srcout[k] = ndim == 2? CK2x(ser_srcstencil_list_2x, poly_order, k, bckey[0])
                             : CK3x(ser_srcstencil_list_3x, poly_order, k, bckey[0], bckey[1]);
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static solstencil_t
fem_poisson_perp_choose_sol_kernels(const struct gkyl_basis* basis)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_solstencil_list[dim].kernels[poly_order];
      
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Assign the right-side vector on the device.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_poisson_perp_set_rhs_cu(gkyl_fem_poisson_perp* up, struct gkyl_array *rhsin);

/**
 * Solve the linear problemon the device.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_perp_solve_cu(gkyl_fem_poisson_perp* up, struct gkyl_array *phiout);
#endif
