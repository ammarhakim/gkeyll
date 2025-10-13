// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h> 

typedef void (*hamil_alpha_quad_conf_t)(const double *w, const double *dxv, const int hamil_pt_edge, 
  const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad); 

typedef double (*lax_flux_nodal_to_modal_t)(const double *dxv, const double *jacob_vel_surf, const double *alpha_quad,
  const double *f_l, const double *f_c, double *lax_nodal_quad, double* GKYL_RESTRICT conf_flux_surf); 

typedef double (*conf_flux_surf_t)(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf); 

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { hamil_alpha_quad_conf_t kernels[4]; } gkyl_hamil_alpha_quad_kern_list;     
typedef struct { lax_flux_nodal_to_modal_t kernels[4]; } gkyl_lax_flux_nodal_to_modal_kern_list;  
typedef struct { conf_flux_surf_t kernels[4]; } gkyl_conf_flux_surf_kern_list;  

struct gkyl_dg_vlasov_conf_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid. 
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index. 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in velocity-space Jacobian. 
  hamil_alpha_quad_conf_t hamil_alpha_quad[3]; // Hamiltonian contribution to alpha_c at quadrature points. 
  lax_flux_nodal_to_modal_t lax_flux_nodal_to_modal[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion. 
  conf_flux_surf_t conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes. 

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_vlasov_conf_flux_surf *on_dev; // pointer to itself or device data.  
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static void 
no_hamil_alpha_quad(const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad)
{
}

GKYL_CU_DH
static void inline
conf_flux_surf_alpha_quad(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad)
{
  up->hamil_alpha_quad[dir](w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad);  
}

// Configuration-space flux computation.
// Separated like this to correctly allocate stack-side arrays for quadrature evaluation.  
GKYL_CU_DH
static double 
conf_flux_surf_1x1v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, so 2^1 = 2 quadrature points. 
  double alpha_quad[2] = {0.0};
  double lax_nodal_quad[2] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x1v_p2(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, so 4^1 = 4 quadrature points. 
  double alpha_quad[4] = {0.0};
  double lax_nodal_quad[4] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x2v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^1*2^1 = 4 quadrature points. 
  double alpha_quad[4] = {0.0};
  double lax_nodal_quad[4] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  } 
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x2v_p2(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^1*4^1 = 16 quadrature points. 
  double alpha_quad[16] = {0.0};
  double lax_nodal_quad[16] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x3v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^1*2^2 = 8 quadrature points. 
  double alpha_quad[8] = {0.0};
  double lax_nodal_quad[8] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_1x3v_p2(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^1*4^2 = 36 quadrature points. 
  double alpha_quad[64] = {0.0};
  double lax_nodal_quad[64] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x2v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^2*2^1 = 8 quadrature points. 
  double alpha_quad[8] = {0.0};
  double lax_nodal_quad[8] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x2v_p2(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 3^3 = 27 quadrature points. 
  double alpha_quad[64] = {0.0};
  double lax_nodal_quad[64] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  } 
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x3v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic in velocity space so 2^2*2^2 = 16 quadrature points. 
  double alpha_quad[16] = {0.0};
  double lax_nodal_quad[16] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

GKYL_CU_DH
static double 
conf_flux_surf_2x3v_p2(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^2*4^2 = 256 quadrature points. 
  double alpha_quad[256] = {0.0};
  double lax_nodal_quad[256] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  } 
} 

GKYL_CU_DH
static double 
conf_flux_surf_3x3v_p1(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge, const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_c, double* GKYL_RESTRICT conf_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic in velocity space so 2^3*2^2 = 32 quadrature points. 
  double alpha_quad[32] = {0.0};
  double lax_nodal_quad[32] = {0.0};

  // Accumulate forces to construct total alpha_c at each quadrature point. 
  conf_flux_surf_alpha_quad(up, dir, w, dxv, hamil_pt_edge, poisson_tensor_conf, hamil, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, 0, alpha_quad,
    f_l, f_c, lax_nodal_quad, conf_flux_surf); 

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  else {
    return cflrate; 
  }
} 

// Velocity-space flux computation assembly function. 
GKYL_CU_D
static const gkyl_conf_flux_surf_kern_list conf_flux_surf_kernels[] = {
  // 1x kernels
  { NULL, conf_flux_surf_1x1v_p1, conf_flux_surf_1x1v_p2, NULL }, // 0
  { NULL, conf_flux_surf_1x2v_p1, conf_flux_surf_1x2v_p2, NULL }, // 1
  { NULL, conf_flux_surf_1x3v_p1, conf_flux_surf_1x3v_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, conf_flux_surf_2x2v_p1, conf_flux_surf_2x2v_p2, NULL }, // 4
  { NULL, conf_flux_surf_2x3v_p1, conf_flux_surf_2x3v_p2, NULL }, // 5
  // 3x kernels
  { NULL, conf_flux_surf_3x3v_p1, NULL, NULL }, // 6
};


// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_to_modal_x_1x1v_ser_p1, lax_flux_nodal_to_modal_x_1x1v_ser_p2, NULL }, // 0
  { NULL, lax_flux_nodal_to_modal_x_1x2v_ser_p1, lax_flux_nodal_to_modal_x_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_to_modal_x_1x3v_ser_p1, lax_flux_nodal_to_modal_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_x_2x2v_ser_p1, lax_flux_nodal_to_modal_x_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_x_2x3v_ser_p1, lax_flux_nodal_to_modal_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_y_2x2v_ser_p1, lax_flux_nodal_to_modal_y_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_y_2x3v_ser_p1, lax_flux_nodal_to_modal_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_z_3x3v_ser_p1, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_ho_lax_flux_nodal_to_modal_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_to_modal_x_1x1v_ser_p1, ho_lax_flux_nodal_to_modal_x_1x1v_ser_p2, NULL }, // 0
  { NULL, lax_flux_nodal_to_modal_x_1x2v_ser_p1, ho_lax_flux_nodal_to_modal_x_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_to_modal_x_1x3v_ser_p1, ho_lax_flux_nodal_to_modal_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_x_2x2v_ser_p1, ho_lax_flux_nodal_to_modal_x_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_x_2x3v_ser_p1, ho_lax_flux_nodal_to_modal_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_ho_lax_flux_nodal_to_modal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_y_2x2v_ser_p1, ho_lax_flux_nodal_to_modal_y_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_y_2x3v_ser_p1, ho_lax_flux_nodal_to_modal_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_ho_lax_flux_nodal_to_modal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_z_3x3v_ser_p1, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_ho_lax_flux_nodal_to_modal_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_ho_lax_flux_nodal_to_modal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_ho_lax_flux_nodal_to_modal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_alpha_quad_x_1x1v_ser_p1, hamil_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_alpha_quad_x_1x2v_ser_p1, hamil_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_alpha_quad_x_1x3v_ser_p1, hamil_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_alpha_quad_x_2x2v_ser_p1, hamil_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_x_2x3v_ser_p1, hamil_alpha_quad_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_alpha_quad_y_2x2v_ser_p1, hamil_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_y_2x3v_ser_p1, hamil_alpha_quad_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_alpha_quad_x_1x1v_ser_p1, hamil_ho_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_alpha_quad_x_1x2v_ser_p1, hamil_ho_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_alpha_quad_x_1x3v_ser_p1, hamil_ho_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_alpha_quad_x_2x2v_ser_p1, hamil_ho_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_x_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_alpha_quad_y_2x2v_ser_p1, hamil_ho_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_y_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};
