#pragma once

#include <gkyl_util.h>
#include <math.h>

// Boundary condition types.
enum gkyl_poisson_bc_type {
  GKYL_POISSON_PERIODIC = 0,
  GKYL_POISSON_DIRICHLET, // sets the value.
  GKYL_POISSON_NEUMANN,   // sets the slope normal to the boundary.
  GKYL_POISSON_ROBIN,  // a combination of dirichlet and neumann.  
  GKYL_POISSON_DIRICHLET_VARYING, // sets the value, spatially varying.
};

// Boundary condition values. Dirichlet and Neumann use only one value,
// Robin uses 3, and periodic ignores the value.
struct gkyl_poisson_bc_value { double v[3]; };

struct gkyl_poisson_bias_plane {
  int dir; // Direction perpendicular to the plane.
  double loc; // Location of the plane in the 'dir' dimension.
  double val; // Biasing value.
};

struct gkyl_poisson_bias_plane_list {
  int num_bias_plane; // Number of bias planes.
  struct gkyl_poisson_bias_plane *bp;
};

struct gkyl_poisson_bc {
  enum gkyl_poisson_bc_type lo_type[GKYL_MAX_CDIM], up_type[GKYL_MAX_CDIM];
  struct gkyl_poisson_bc_value lo_value[GKYL_MAX_CDIM], up_value[GKYL_MAX_CDIM];
  // Function specifying a spatially varying BC.
  void (*bc_value_func)(double t, const double *xn, double *phi_wall_up_out, void *ctx);
  void *bc_value_func_ctx;

  // Additional attributes to apply a bias plane at the extremal z values only.
  bool contains_lower_z_edge, contains_upper_z_edge;
};

GKYL_CU_DH
static inline int idx_to_inup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between upper cells and
  // elsewhere.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == num_cells[d]) iout += (int)(pow(2,d)+0.5);
  }
  return iout;
}

GKYL_CU_DH
static inline int idx_to_inloup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[d] == num_cells[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}
