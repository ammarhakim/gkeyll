// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h> 

typedef void (*dg_gr_maxwell_alpha_quad_conf_t)(const double *w, const double *dxv, const int geom_edge, const double *lapse, 
  const double *shift, const double *h_ij, const double *det_h,   
  const double *field_no_J_con_l, const double *field_no_J_con_r, 
  double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
  double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 

typedef double (*lax_flux_t)(const double *dxv, const double *det_h, const double *flux_l, 
  const double *flux_r, const double *alpha_quad, const double *field_no_J_con_l, 
  const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf); 

typedef double (*roe_flux_t)(const double *dxv, 
  const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
  double* GKYL_RESTRICT conf_flux_surf); 

typedef double (*conf_flux_surf_t)(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l, 
  const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf); 

// for use in kernel tables
typedef struct { dg_gr_maxwell_alpha_quad_conf_t kernels[4]; } gkyl_dg_gr_maxwell_alpha_quad_kern_list;     
typedef struct { lax_flux_t kernels[4]; } gkyl_lax_flux_kern_list;  
typedef struct { roe_flux_t kernels[4]; } gkyl_roe_flux_kern_list;  
typedef struct { conf_flux_surf_t kernels[4]; } gkyl_conf_flux_surf_kern_list;  

struct gkyl_dg_gr_maxwell_conf_flux_surf {
  struct gkyl_rect_grid conf_grid; // Conf-space grid. 
  int cdim; // Configuration-space dimensions.
  dg_gr_maxwell_alpha_quad_conf_t dg_gr_maxwell_alpha_quad[3]; // Contribution to alpha_c at quadrature points. 
  lax_flux_t lax_flux[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion. 
  roe_flux_t roe_flux[3]; // Convert nodal Roe-Friedrichs flux to modal surface expansion. 
  conf_flux_surf_t conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes. 
  int theta_pole_lo[GKYL_MAX_CDIM]; // (Lower side) Default zeros, but 1 if any directions use theta-pole BC's.
  int theta_pole_up[GKYL_MAX_CDIM]; // (Upper side) Default zeros, but 1 if any directions use theta-pole BC's.

  uint32_t flags;
  bool use_lax;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_conf_flux_surf *on_dev; // pointer to itself or device data.  
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static void 
no_dg_gr_maxwell_alpha_quad(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse, 
  const double *shift, const double *h_ij, const double *det_h,   
  const double *field_no_J_con_l, const double *field_no_J_con_r, 
  double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad)
{
}

GKYL_CU_DH
static void inline
conf_flux_surf_alpha_quad(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse, 
  const double *shift, const double *h_ij, const double *det_h,   
  const double *field_no_J_con_l, const double *field_no_J_con_r, 
  double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
  double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad)
{
  up->dg_gr_maxwell_alpha_quad[dir](w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad);  
}

// Configuration-space flux computation.
// Separated like this to correctly allocate stack-side arrays for quadrature evaluation.  
GKYL_CU_DH
static double 
conf_flux_surf_1x_p1(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 1 quadrature node on a 1x surface, 6 field components.
  double alpha_quad[1] = {0.0}; 
  double flux_l[6] = {0.0};
  double flux_r[6] = {0.0};
  double A_plus_dQ[6] = {0.0};
  double A_minus_dQ[6] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.  
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate; 
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x_p2(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 1 quadrature node on a 1x surface, 6 field components.
  double alpha_quad[1] = {0.0}; 
  double flux_l[6] = {0.0};
  double flux_r[6] = {0.0};
  double A_plus_dQ[6] = {0.0};
  double A_minus_dQ[6] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.    
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate; 
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x_p1(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 2 quadrature nodes on a 2x surface, 6 field components.
  double alpha_quad[2] = {0.0}; 
  double flux_l[12] = {0.0};
  double flux_r[12] = {0.0};
  double A_plus_dQ[12] = {0.0};
  double A_minus_dQ[12] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.    
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate; 
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x_p2(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 3 quadrature nodes on a 2x p2 surface, 6 field components.
  double alpha_quad[3] = {0.0}; 
  double flux_l[18] = {0.0};
  double flux_r[18] = {0.0};
  double A_plus_dQ[18] = {0.0};
  double A_minus_dQ[18] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.    
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate;  
} 

GKYL_CU_DH
static double 
conf_flux_surf_3x_p1(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 4 quadrature nodes on a 3x surface, 6 field components.
  double alpha_quad[4] = {0.0}; 
  double flux_l[24] = {0.0};
  double flux_r[24] = {0.0};
  double A_plus_dQ[24] = {0.0};
  double A_minus_dQ[24] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.    
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate; 
} 

GKYL_CU_DH
static double 
conf_flux_surf_3x_p2(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_no_J_con_l,
   const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf)
{
  // 4 quadrature nodes on a 3x surface, 6 field components.
  double alpha_quad[9] = {0.0}; 
  double flux_l[54] = {0.0};
  double flux_r[54] = {0.0};
  double A_plus_dQ[54] = {0.0};
  double A_minus_dQ[54] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.    
  conf_flux_surf_alpha_quad(up, dir, w, dxv, geom_edge, lapse, shift, 
    h_ij, det_h, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, det_h, flux_l, flux_r, alpha_quad,
      field_no_J_con_l, field_no_J_con_r, conf_flux_surf); 
  }
  else {
    // Compute nodal Roe flux 
    cflrate = up->roe_flux[dir](dxv, flux_l, A_minus_dQ, alpha_quad, conf_flux_surf); 
  }

  // Return the cfl rate
  return cflrate; 
} 

// Configuration-space flux computation assembly function. 
GKYL_CU_D
static const gkyl_conf_flux_surf_kern_list conf_flux_surf_kernels[] = {
  // 1x kernels
  { NULL, conf_flux_surf_1x_p1, conf_flux_surf_1x_p2, NULL }, // 0
  // 2x kernels
  { NULL, conf_flux_surf_2x_p1, conf_flux_surf_2x_p2, NULL }, // 1
  // 3x kernels
  { NULL, conf_flux_surf_3x_p1, conf_flux_surf_3x_p2, NULL }, // 2
};


// Nodal Lax-Friedrichs flux (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_kern_list ser_lax_flux_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_x_1x_ser_p1, lax_flux_x_1x_ser_p2, NULL }, // 0
  // 2x kernels
  { NULL, lax_flux_x_2x_ser_p1, lax_flux_x_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_x_3x_ser_p1, lax_flux_x_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_lax_flux_kern_list ser_lax_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, lax_flux_y_2x_ser_p1, lax_flux_y_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_y_3x_ser_p1, lax_flux_y_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_lax_flux_kern_list ser_lax_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_z_3x_ser_p1, lax_flux_z_3x_ser_p2, NULL }, // 2
};

// Nodal Lax-Friedrichs flux (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_kern_list ten_lax_flux_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_x_1x_tensor_p1, lax_flux_x_1x_tensor_p2, NULL }, // 0
  // 2x kernels
  { NULL, lax_flux_x_2x_tensor_p1, lax_flux_x_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_x_3x_tensor_p1, lax_flux_x_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_lax_flux_kern_list ten_lax_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, lax_flux_y_2x_tensor_p1, lax_flux_y_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_y_3x_tensor_p1, lax_flux_y_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_lax_flux_kern_list ten_lax_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, lax_flux_z_3x_tensor_p1, lax_flux_z_3x_tensor_p2, NULL }, // 2
};

// Nodal roe-Friedrichs flux (Serendipity basis). 
GKYL_CU_D
static const gkyl_roe_flux_kern_list ser_roe_flux_x_kernels[] = {
  // 1x kernels
  { NULL, roe_flux_x_1x_ser_p1, roe_flux_x_1x_ser_p2, NULL }, // 0
  // 2x kernels
  { NULL, roe_flux_x_2x_ser_p1, roe_flux_x_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_x_3x_ser_p1, roe_flux_x_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_roe_flux_kern_list ser_roe_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, roe_flux_y_2x_ser_p1, roe_flux_y_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_y_3x_ser_p1, roe_flux_y_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_roe_flux_kern_list ser_roe_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_z_3x_ser_p1, roe_flux_z_3x_ser_p2, NULL }, // 2
};

// Nodal roe-Friedrichs flux (Tensor basis). 
GKYL_CU_D
static const gkyl_roe_flux_kern_list ten_roe_flux_x_kernels[] = {
  // 1x kernels
  { NULL, roe_flux_x_1x_tensor_p1, roe_flux_x_1x_tensor_p2, NULL }, // 0
  // 2x kernels
  { NULL, roe_flux_x_2x_tensor_p1, roe_flux_x_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_x_3x_tensor_p1, roe_flux_x_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_roe_flux_kern_list ten_roe_flux_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, roe_flux_y_2x_tensor_p1, roe_flux_y_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_y_3x_tensor_p1, roe_flux_y_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_roe_flux_kern_list ten_roe_flux_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, roe_flux_z_3x_tensor_p1, roe_flux_z_3x_tensor_p2, NULL }, // 2
};


// alpha_c evaluated at quadrature points. 
GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ser_dg_gr_maxwell_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_1x_ser_p1, dg_gr_maxwell_alpha_quad_x_1x_ser_p2, NULL }, // 0
  // 2x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_2x_ser_p1, dg_gr_maxwell_alpha_quad_x_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_3x_ser_p1, dg_gr_maxwell_alpha_quad_x_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ser_dg_gr_maxwell_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, dg_gr_maxwell_alpha_quad_y_2x_ser_p1, dg_gr_maxwell_alpha_quad_y_2x_ser_p2, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_y_3x_ser_p1, dg_gr_maxwell_alpha_quad_y_3x_ser_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ser_dg_gr_maxwell_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_z_3x_ser_p1, dg_gr_maxwell_alpha_quad_z_3x_ser_p2, NULL }, // 2
};

// alpha_c evaluated at quadrature points. 
GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ten_dg_gr_maxwell_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_1x_tensor_p1, dg_gr_maxwell_alpha_quad_x_1x_tensor_p2, NULL }, // 0
  // 2x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_2x_tensor_p1, dg_gr_maxwell_alpha_quad_x_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_x_3x_tensor_p1, dg_gr_maxwell_alpha_quad_x_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ten_dg_gr_maxwell_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, dg_gr_maxwell_alpha_quad_y_2x_tensor_p1, dg_gr_maxwell_alpha_quad_y_2x_tensor_p2, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_y_3x_tensor_p1, dg_gr_maxwell_alpha_quad_y_3x_tensor_p2, NULL }, // 2
};

GKYL_CU_D
static const gkyl_dg_gr_maxwell_alpha_quad_kern_list ten_dg_gr_maxwell_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, dg_gr_maxwell_alpha_quad_z_3x_tensor_p1, dg_gr_maxwell_alpha_quad_z_3x_tensor_p2, NULL }, // 2
};
