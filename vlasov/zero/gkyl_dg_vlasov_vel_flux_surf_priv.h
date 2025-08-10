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

typedef void (*hamil_alpha_quad_t)(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad); 

typedef void (*E_alpha_quad_t)(const double *dxv, 
  const double *qmem, double* GKYL_RESTRICT alpha_quad);  

typedef void (*phi_alpha_quad_t)(const double *dxv, 
  const double *phi_tot, double* GKYL_RESTRICT alpha_quad); 

typedef void (*B_alpha_quad_t)(const double *dxv, 
  const double *jacob_vel, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad); 

typedef void (*rad_alpha_quad_t)(const double *dxv, 
  const double *rad, double* GKYL_RESTRICT alpha_quad); 

typedef double (*lax_flux_nodal_to_modal_t)(const double *dxv, 
  const double *jacob_vel, const double *alpha_quad, 
  const double *f_l, const double *f_c, 
  double *lax_nodal_quad, double* GKYL_RESTRICT vel_flux_surf); 

typedef double (*vel_flux_surf_t)(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf); 

typedef double (*vel_flux_surf_edge_t)(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_c, double* GKYL_RESTRICT vel_flux_surf); 

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
typedef struct { hamil_alpha_quad_t kernels[4]; } gkyl_hamil_alpha_quad_kern_list;  
typedef struct { E_alpha_quad_t kernels[4]; } gkyl_E_alpha_quad_kern_list;  
typedef struct { phi_alpha_quad_t kernels[4]; } gkyl_phi_alpha_quad_kern_list;  
typedef struct { B_alpha_quad_t kernels[4]; } gkyl_B_alpha_quad_kern_list;  
typedef struct { rad_alpha_quad_t kernels[4]; } gkyl_rad_alpha_quad_kern_list;    
typedef struct { lax_flux_nodal_to_modal_t kernels[4]; } gkyl_lax_flux_nodal_to_modal_kern_list;  
typedef struct { vel_flux_surf_t kernels[4]; } gkyl_vel_flux_surf_kern_list;  

struct gkyl_dg_vlasov_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid. 
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index. 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in velocity-space Jacobian. 
  hamil_alpha_quad_t hamil_alpha_quad[3]; // Hamiltonian contribution to alpha_v at quadrature points. 
  E_alpha_quad_t E_alpha_quad[3]; // Lorentz force contribution from electric field to alpha_v at quadrature points. 
  phi_alpha_quad_t phi_alpha_quad[3]; // Scalar potential, -grad(phi), force contribution to alpha_v at quadrature points. 
  B_alpha_quad_t B_alpha_quad[3]; // Lorentz force contribution from magnetic field to alpha_v at quadrature points. 
  rad_alpha_quad_t rad_alpha_quad[3]; // Radiation drag force contribution to alpha_v at quadrature points. 
  lax_flux_nodal_to_modal_t lax_flux_nodal_to_modal[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion. 
  vel_flux_surf_t vel_flux_surf; // Assembly function for computing modal surface expansion of velocity-space fluxes. 
  vel_flux_surf_edge_t vel_flux_surf_edge; // Edge-of-velocity-space Assembly function for computing velocity-space fluxes. 

  uint32_t flags;
  struct gkyl_dg_vlasov_vel_flux_surf *on_dev; // pointer to itself or device data.  
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static void 
no_hamil_alpha_quad(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_E_alpha_quad(const double *dxv, 
  const double *qmem, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_phi_alpha_quad(const double *dxv, 
  const double *phi, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_B_alpha_quad(const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_rad_alpha_quad(const double *dxv, 
  const double *rad, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static double
no_vel_flux_surf_edge(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Zero-flux boundary conditions, so immediately return. 
  return 0.0; 
}

GKYL_CU_DH
static void inline
vel_flux_surf_alpha_quad(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  double* GKYL_RESTRICT alpha_quad)
{
  up->hamil_alpha_quad[dir](dxv, hamil, alpha_quad); 
  up->E_alpha_quad[dir](dxv, qmem, alpha_quad); 
  up->phi_alpha_quad[dir](dxv, pot_tot, alpha_quad); 
  up->B_alpha_quad[dir](dxv, jacob_vel, hamil, qmem, alpha_quad); 
  up->rad_alpha_quad[dir](dxv, rad, alpha_quad); 
}

// Velocity-space flux computation.
// Separated like this to correctly allocate stack-side arrays for quadrature evaluation.  
GKYL_CU_DH
static double 
vel_flux_surf_1x1v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, so 2^1 = 2 quadrature points. 
  double alpha_quad[2] = {0.0};
  double lax_nodal_quad[2] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_1x1v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, so 4^1 = 4 quadrature points. 
  double alpha_quad[4] = {0.0};
  double lax_nodal_quad[4] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_1x2v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^1*2^1 = 4 quadrature points. 
  double alpha_quad[4] = {0.0};
  double lax_nodal_quad[4] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_1x2v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^1*4^1 = 16 quadrature points. 
  double alpha_quad[16] = {0.0};
  double lax_nodal_quad[16] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_1x3v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^1*2^2 = 8 quadrature points. 
  double alpha_quad[8] = {0.0};
  double lax_nodal_quad[8] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_1x3v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^1*4^2 = 36 quadrature points. 
  double alpha_quad[64] = {0.0};
  double lax_nodal_quad[64] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x1v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, so 2^2 = 4 quadrature points. 
  double alpha_quad[4] = {0.0};
  double lax_nodal_quad[4] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x1v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, so 4^2 = 16 quadrature points. 
  double alpha_quad[16] = {0.0};
  double lax_nodal_quad[16] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x2v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic/cubic in velocity space so 2^2*2^1 = 8 quadrature points. 
  double alpha_quad[8] = {0.0};
  double lax_nodal_quad[8] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x2v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^2*4 = 64 quadrature points. 
  double alpha_quad[64] = {0.0};
  double lax_nodal_quad[64] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x3v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic in velocity space so 2^2*2^2 = 16 quadrature points. 
  double alpha_quad[16] = {0.0};
  double lax_nodal_quad[16] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_2x3v_p2(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a sixth order polynomial in configuration space, and a
  // sixth order polynomial in velocity space, so 4^2*4^2 = 256 quadrature points. 
  double alpha_quad[256] = {0.0};
  double lax_nodal_quad[256] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
vel_flux_surf_3x3v_p1(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Needs to integrate a cubic in configuration space, and a 
  // quadratic in velocity space so 2^3*2^2 = 32 quadrature points. 
  double alpha_quad[32] = {0.0};
  double lax_nodal_quad[32] = {0.0};

  // Accumulate forces to construct total alpha_v at each quadrature point. 
  vel_flux_surf_alpha_quad(up, dir, dxv, jacob_vel, hamil, qmem, pot_tot, rad, alpha_quad); 

  // Compute nodal Lax-Friedrichs flux and convert back to modal expansion of flux.
  double cflrate = up->lax_flux_nodal_to_modal[dir](dxv, jacob_vel, alpha_quad, 
    f_l, f_c, lax_nodal_quad, vel_flux_surf); 

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
static const gkyl_vel_flux_surf_kern_list vel_flux_surf_kernels[] = {
  // 1x kernels
  { NULL, vel_flux_surf_1x1v_p1, vel_flux_surf_1x1v_p2, NULL }, // 0
  { NULL, vel_flux_surf_1x2v_p1, vel_flux_surf_1x2v_p2, NULL }, // 1
  { NULL, vel_flux_surf_1x3v_p1, vel_flux_surf_1x3v_p2, NULL }, // 2
  // 2x kernels
  { NULL, vel_flux_surf_2x1v_p1, vel_flux_surf_2x1v_p2, NULL }, // 3
  { NULL, vel_flux_surf_2x2v_p1, vel_flux_surf_2x2v_p2, NULL }, // 4
  { NULL, vel_flux_surf_2x3v_p1, vel_flux_surf_2x3v_p2, NULL }, // 5
  // 3x kernels
  { NULL, vel_flux_surf_3x3v_p1, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_vx_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_to_modal_vx_1x1v_ser_p1, lax_flux_nodal_to_modal_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, lax_flux_nodal_to_modal_vx_1x2v_ser_p1, lax_flux_nodal_to_modal_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_to_modal_vx_1x3v_ser_p1, lax_flux_nodal_to_modal_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_to_modal_vx_2x1v_ser_p1, lax_flux_nodal_to_modal_vx_2x1v_ser_p2, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_vx_2x2v_ser_p1, lax_flux_nodal_to_modal_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_vx_2x3v_ser_p1, lax_flux_nodal_to_modal_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_to_modal_vy_1x2v_ser_p1, lax_flux_nodal_to_modal_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_to_modal_vy_1x3v_ser_p1, lax_flux_nodal_to_modal_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_to_modal_vy_2x2v_ser_p1, lax_flux_nodal_to_modal_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_vy_2x3v_ser_p1, lax_flux_nodal_to_modal_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list ser_lax_flux_nodal_to_modal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_to_modal_vz_1x3v_ser_p1, lax_flux_nodal_to_modal_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_to_modal_vz_2x3v_ser_p1, lax_flux_nodal_to_modal_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_to_modal_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, lax_flux_nodal_to_modal_vx_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, lax_flux_nodal_to_modal_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, lax_flux_nodal_to_modal_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, lax_flux_nodal_to_modal_vx_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, lax_flux_nodal_to_modal_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, lax_flux_nodal_to_modal_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, lax_flux_nodal_to_modal_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, lax_flux_nodal_to_modal_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_to_modal_kern_list tensor_lax_flux_nodal_to_modal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, lax_flux_nodal_to_modal_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, hamil_alpha_quad_vx_1x1v_ser_p1, hamil_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_alpha_quad_vx_1x2v_ser_p1, hamil_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_alpha_quad_vx_1x3v_ser_p1, hamil_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, hamil_alpha_quad_vx_2x1v_ser_p1, hamil_alpha_quad_vx_2x1v_ser_p2, NULL }, // 3
  { NULL, hamil_alpha_quad_vx_2x2v_ser_p1, hamil_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_vx_2x3v_ser_p1, hamil_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_alpha_quad_vy_2x2v_ser_p1, hamil_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_alpha_quad_vy_2x3v_ser_p1, hamil_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, hamil_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Tensor basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, hamil_alpha_quad_vx_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, hamil_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, hamil_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, hamil_alpha_quad_vx_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, hamil_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, hamil_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Serendipity basis). 
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_ser_p1, E_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, E_alpha_quad_vx_1x2v_ser_p1, E_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_ser_p1, E_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_ser_p1, E_alpha_quad_vx_2x1v_ser_p2, NULL }, // 3
  { NULL, E_alpha_quad_vx_2x2v_ser_p1, E_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_ser_p1, E_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_ser_p1, E_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_ser_p1, E_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_ser_p1, E_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_ser_p1, E_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_ser_p1, E_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_ser_p1, E_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Tensor basis). 
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, E_alpha_quad_vx_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, E_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, E_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, E_alpha_quad_vx_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, E_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, E_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, E_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, E_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, E_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_ser_p1, phi_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_ser_p1, phi_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_ser_p1, phi_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_ser_p1, phi_alpha_quad_vx_2x1v_ser_p2, NULL }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_ser_p1, phi_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_ser_p1, phi_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_ser_p1, phi_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_ser_p1, phi_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Tensor basis). 
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, phi_alpha_quad_vx_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, phi_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, phi_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, phi_alpha_quad_vx_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, phi_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, phi_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, NULL }, // 0
  { NULL, B_alpha_quad_vx_1x2v_ser_p1, B_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_alpha_quad_vx_1x3v_ser_p1, B_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, NULL }, // 3
  { NULL, B_alpha_quad_vx_2x2v_ser_p1, B_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_alpha_quad_vx_2x3v_ser_p1, B_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_alpha_quad_vy_1x2v_ser_p1, B_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_alpha_quad_vy_1x3v_ser_p1, B_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_alpha_quad_vy_2x2v_ser_p1, B_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_alpha_quad_vy_2x3v_ser_p1, B_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_alpha_quad_vz_1x3v_ser_p1, B_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_alpha_quad_vz_2x3v_ser_p1, B_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Tensor basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, no_B_alpha_quad, NULL }, // 0
  { NULL, NULL, B_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, B_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, no_B_alpha_quad, NULL }, // 3
  { NULL, NULL, B_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, B_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, B_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, B_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, B_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};


// alpha_v evaluated at quadrature points for the radiation drag force (Serendipity basis). 
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_ser_p1, rad_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_ser_p1, rad_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_ser_p1, rad_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_ser_p1, rad_alpha_quad_vx_2x1v_ser_p2, NULL }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_ser_p1, rad_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_ser_p1, rad_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_ser_p1, rad_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_ser_p1, rad_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_ser_p1, rad_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_ser_p1, rad_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_ser_p1, rad_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_ser_p1, rad_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Tensor basis). 
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, rad_alpha_quad_vx_1x1v_tensor_p2, NULL }, // 0
  { NULL, NULL, rad_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, rad_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, rad_alpha_quad_vx_2x1v_tensor_p2, NULL }, // 3
  { NULL, NULL, rad_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, rad_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, NULL, rad_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, rad_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, rad_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};