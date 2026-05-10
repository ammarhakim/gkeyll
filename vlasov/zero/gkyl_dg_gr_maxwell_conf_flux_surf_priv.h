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
  const double *field_con_l, const double *field_con_r, 
  const double *field_no_J_con_l, const double *field_no_J_con_r, 
  double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
  double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad); 

typedef double (*lax_flux_t)(const double *dxv, const int geom_edge, const double *det_h, const double *flux_l,
  const double *flux_r, const double *alpha_quad, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf);

// Curved-norm LLF flux: takes the face nodal h_ij. The dissipation matrix
// is the spatial 3-metric h_ij (block-diagonal between D and B sub-vectors),
// which is positive-definite on the spatial slice (only the spacetime metric
// becomes indefinite in the ergoregion via g_{tt}). Provides positive-
// definite damping in the spatial 3-norm. Tested with J_c*alpha*h_ij at N=24
// — that variant gave 40% lower max Re(lambda) at N=24 but lost the machine-
// zero result at N=12, with the same overall N^3.8 scaling, so the simpler
// h_ij-only form was kept.
typedef double (*lax_flux_curved_t)(const double *dxv, const int geom_edge,
  const double *h_ij_nodal, const double *det_h,
  const double *flux_l, const double *flux_r, const double *alpha_quad,
  const double *field_con_l, const double *field_con_r,
  double* GKYL_RESTRICT conf_flux_surf);

// Tetrad-frame Roe flux: per-quad-node tetrad transform + flat-Maxwell wave
// decomposition (paper eq. 57) + back-transform. Self-contained -- computes
// its own wave speeds, so no alpha_quad input. Currently implemented for
// dir=0 (radial face) only; theta-direction kernels would need separate
// pole-face handling.
typedef double (*tetrad_roe_flux_t)(const double *dxv, const int geom_edge,
  const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *det_h,
  const double *flux_l, const double *flux_r,
  const double *field_con_l, const double *field_con_r,
  double* GKYL_RESTRICT conf_flux_surf);

typedef double (*roe_flux_t)(const double *dxv, 
  const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
  double* GKYL_RESTRICT conf_flux_surf); 

typedef double (*conf_flux_surf_t)(struct gkyl_dg_gr_maxwell_conf_flux_surf *up,
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h, const double *field_con_l,
  const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf);

// One-sided outflow flux kernel: takes only the SKIN cell field (no ghost) and the face
// nodal metric, builds A^outgoing from per-quad-node eigenvalues (positive at upper edge,
// negative at lower edge), applies it to U_skin via A_dQ_<dir>_calc, and writes the
// resulting outgoing-only numerical flux into the conf_flux_surf slot for this face.
// `edge` is GKYL_LOWER_EDGE (=0) or GKYL_UPPER_EDGE (=1).
typedef double (*outflow_flux_t)(const double *dxv, int edge,
  const double *lapse, const double *shift, const double *h_ij, const double *det_h,
  const double *field_con_skin, double* GKYL_RESTRICT conf_flux_surf);

// for use in kernel tables
typedef struct { dg_gr_maxwell_alpha_quad_conf_t kernels[4]; } gkyl_dg_gr_maxwell_alpha_quad_kern_list;
typedef struct { lax_flux_t kernels[4]; } gkyl_lax_flux_kern_list;
typedef struct { lax_flux_curved_t kernels[4]; } gkyl_lax_flux_curved_kern_list;
typedef struct { tetrad_roe_flux_t kernels[4]; } gkyl_tetrad_roe_flux_kern_list;
typedef struct { roe_flux_t kernels[4]; } gkyl_roe_flux_kern_list;
typedef struct { conf_flux_surf_t kernels[4]; } gkyl_conf_flux_surf_kern_list;
typedef struct { outflow_flux_t kernels[4]; } gkyl_outflow_flux_kern_list;

struct gkyl_dg_gr_maxwell_conf_flux_surf {
  struct gkyl_rect_grid conf_grid; // Conf-space grid.
  int cdim; // Configuration-space dimensions.
  dg_gr_maxwell_alpha_quad_conf_t dg_gr_maxwell_alpha_quad[3]; // Contribution to alpha_c at quadrature points.
  lax_flux_t lax_flux[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion.
  lax_flux_curved_t lax_flux_curved[3]; // Curved-norm LLF flux (positive-definite spatial-3-metric dissipation).
  tetrad_roe_flux_t tetrad_roe_flux[3]; // Tetrad-frame Roe flux (per-quad-node, paper eq. 57). Currently dir=0 only.
  roe_flux_t roe_flux[3]; // Convert nodal Roe-Friedrichs flux to modal surface expansion.
  outflow_flux_t outflow_flux[3]; // One-sided outflow flux for characteristic-based BC (per direction).
  outflow_flux_t outflow_flux_curved[3]; // Curved-norm variant: A^outgoing*U_skin + h_ij-weighted dissipation.
  outflow_flux_t outflow_flux_tetrad_roe[3]; // Tetrad-frame variant: paper eq. 57 wave decomposition in tetrad,
                                              // outgoing-only projection. Currently only x-direction implemented.
  conf_flux_surf_t conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes.
  int theta_pole_lo[GKYL_MAX_CDIM]; // (Lower side) Default zeros, but 1 if any directions use theta-pole BC's.
  int theta_pole_up[GKYL_MAX_CDIM]; // (Upper side) Default zeros, but 1 if any directions use theta-pole BC's.
  int outflow_lo[GKYL_MAX_CDIM]; // (Lower side) Default zeros, 1 if direction uses GKYL_FIELD_OUTFLOW.
  int outflow_up[GKYL_MAX_CDIM]; // (Upper side) Default zeros, 1 if direction uses GKYL_FIELD_OUTFLOW.

  uint32_t flags;
  bool use_lax;
  bool use_curved_norm; // when use_lax is also true: dispatch lax_flux_curved instead of lax_flux
  bool use_tetrad_flux; // when use_lax is also true and dir == 0: override alpha_quad with the tetrad
                        // wave speed (alpha/sqrt(gamma^xx) + |beta^x|/gamma^xx) and dispatch to the
                        // standard Euclidean lax_flux kernel. See note in gkyl_dg_gr_maxwell_conf_flux_surf.h.
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_conf_flux_surf *on_dev; // pointer to itself or device data.
};

GKYL_CU_DH
static void inline
conf_flux_surf_alpha_quad(struct gkyl_dg_gr_maxwell_conf_flux_surf *up,
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h,
  const double *field_con_l, const double *field_con_r,
  const double *field_no_J_con_l, const double *field_no_J_con_r,
  double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
  double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT alpha_quad)
{
  up->dg_gr_maxwell_alpha_quad[dir](w, dxv, geom_edge, lapse, shift,
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad);
}

// Note: an earlier prototype here implemented "tetrad LLF" by overriding
// alpha_quad with a tetrad-frame max wave speed and dispatching to the
// standard (Euclidean) lax_flux kernel. The N=12 spectrum was strictly
// worse than the curved-LLF baseline (4 unstable modes appeared at max
// Re ~ 5e-2 vs the baseline's machine zero), because the tetrad-LLF
// dissipation in coord basis collapses to a SCALAR * unweighted jump,
// which is weaker than the h_ij-weighted dissipation in the curved-LLF
// kernel at our (gamma^xx ~= 0.7) parameters. The natural follow-up
// (per-quadrature-point tetrad-Roe) does not have this collapse: the
// flat-Maxwell eigenvectors back-transform into non-trivial coord-basis
// directions, so the kernel actually needs to do the tetrad math. That
// is what tetrad_roe_flux[] dispatches to below.

// Configuration-space flux computation.
// Separated like this to correctly allocate stack-side arrays for quadrature evaluation.  
GKYL_CU_DH
static double 
conf_flux_surf_1x_p1(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int geom_edge, const double *lapse,
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, geom_edge,det_h, flux_l, flux_r, alpha_quad,
      field_con_l, field_con_r, conf_flux_surf); 
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
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, geom_edge,det_h, flux_l, flux_r, alpha_quad,
      field_con_l, field_con_r, conf_flux_surf); 
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
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad);

  double cflrate;
  if (up->use_lax) {
    if (up->use_tetrad_flux && up->tetrad_roe_flux[dir] != 0) {
      // Tetrad-frame Roe: per-quad-node tetrad transform, flat-Maxwell wave
      // decomposition, and back-transform. Self-contained kernel; computes its
      // own wave speeds internally. Currently implemented for dir=0 (radial,
      // x-direction flat-Maxwell eigenvectors) and dir=1 (theta, y-direction
      // eigenvectors with theta-pole reflective fallback). Falls through if
      // the dir's kernel pointer is NULL (e.g., dir=2 in 2D code is always
      // unused).
      cflrate = up->tetrad_roe_flux[dir](dxv, geom_edge, lapse, shift,
        h_ij, det_h, flux_l, flux_r, field_con_l, field_con_r, conf_flux_surf);
    } else if (up->use_curved_norm && up->lax_flux_curved[dir] != 0) {
      // Curved-norm LLF: dissipation matrix is the spatial 3-metric h_ij,
      // positive-definite on the spatial slice. Damps in the spatial 3-norm.
      cflrate = up->lax_flux_curved[dir](dxv, geom_edge, h_ij, det_h,
        flux_l, flux_r, alpha_quad, field_con_l, field_con_r, conf_flux_surf);
    } else {
      // Standard Lax-Friedrichs: per-component diagonal dissipation in flat
      // L^2 norm.
      cflrate = up->lax_flux[dir](dxv, geom_edge, det_h, flux_l, flux_r, alpha_quad,
        field_con_l, field_con_r, conf_flux_surf);
    }
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
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, geom_edge,det_h, flux_l, flux_r, alpha_quad,
      field_con_l, field_con_r, conf_flux_surf); 
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
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, geom_edge,det_h, flux_l, flux_r, alpha_quad,
      field_con_l, field_con_r, conf_flux_surf); 
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
  const double *shift, const double *h_ij, const double *det_h,const double *field_con_l,
   const double *field_con_r, const double *field_no_J_con_l,
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
    h_ij, det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r,
     A_plus_dQ, A_minus_dQ, flux_l, flux_r, alpha_quad); 

  double cflrate;   
  if (up->use_lax) {
    // Compute nodal Lax-Friedrichs flux 
    cflrate = up->lax_flux[dir](dxv, geom_edge,det_h, flux_l, flux_r, alpha_quad,
      field_con_l, field_con_r, conf_flux_surf); 
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


// Curved-norm LLF flux kernels: only 2x ser p1 currently implemented (the
// prototype). Other entries are NULL — when use_curved_norm is set but the
// kernel is NULL for the active basis/order, conf_flux_surf_2x_p1 falls
// back to standard lax_flux.
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ser_lax_flux_curved_x_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 1x — not yet implemented
  { NULL, lax_flux_x_2x_ser_p1_curved, NULL, NULL }, // 2x ser p1 only
  { NULL, NULL, NULL, NULL }, // 3x — not yet implemented
};
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ser_lax_flux_curved_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, lax_flux_y_2x_ser_p1_curved, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ser_lax_flux_curved_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ten_lax_flux_curved_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ten_lax_flux_curved_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_lax_flux_curved_kern_list ten_lax_flux_curved_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

// Tetrad-frame Roe flux kernels (paper eq. 57 wave decomposition per
// quadrature node). Currently only the radial face (x-direction) at 2x
// ser p=1 is implemented. y/z directions and other basis/order combinations
// are NULL and the wrapper falls back to use_curved_norm / standard LLF.
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ser_tetrad_roe_flux_x_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 1x — not yet implemented
  { NULL, lax_flux_x_2x_ser_p1_tetrad_roe, NULL, NULL }, // 2x ser p1 only
  { NULL, NULL, NULL, NULL }, // 3x — not yet implemented
};
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ser_tetrad_roe_flux_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, lax_flux_y_2x_ser_p1_tetrad_roe, NULL, NULL }, // theta-pole reflective fallback inside kernel
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ser_tetrad_roe_flux_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ten_tetrad_roe_flux_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ten_tetrad_roe_flux_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_tetrad_roe_flux_kern_list ten_tetrad_roe_flux_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

// One-sided outflow flux kernels (characteristic-based BC).
// NULL entries are unimplemented variants; outflow currently only ships for
// 2x serendipity p=1. Adding p=2, tensor and 3D follows the same recipe:
// see outflow_flux_x_2x_ser_p1.c for the template.
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_x_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0: 1x — not yet implemented
  { NULL, outflow_flux_x_2x_ser_p1, NULL, NULL }, // 1: 2x ser p1 only
  { NULL, NULL, NULL, NULL }, // 2: 3x — not yet implemented
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, outflow_flux_y_2x_ser_p1, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
};

// Curved-norm outflow kernels (h_ij-weighted dissipation added on top of
// the standard A^outgoing*U_skin). Only 2x ser p1 implemented as prototype.
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_curved_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, outflow_flux_x_2x_ser_p1_curved, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_curved_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, outflow_flux_y_2x_ser_p1_curved, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_curved_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_curved_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_curved_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_curved_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2: 3x — not yet implemented
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // 2x tensor — not yet implemented
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};

// Tetrad-Roe outflow kernels (paper eq. 57 wave decomposition in tetrad,
// outgoing-only projection). Only x-direction 2x ser p1 implemented; the
// y-direction analog would need the same theta-pole handling as the LLF
// kernel, deferred.
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_tetrad_roe_x_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 1x — not yet implemented
  { NULL, outflow_flux_x_2x_ser_p1_tetrad_roe, NULL, NULL }, // 2x ser p1
  { NULL, NULL, NULL, NULL }, // 3x — not yet implemented
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_tetrad_roe_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL }, // y-direction tetrad-Roe outflow not yet implemented
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ser_outflow_flux_tetrad_roe_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_tetrad_roe_x_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_tetrad_roe_y_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
};
GKYL_CU_D
static const gkyl_outflow_flux_kern_list ten_outflow_flux_tetrad_roe_z_kernels[] = {
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
  { NULL, NULL, NULL, NULL },
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
