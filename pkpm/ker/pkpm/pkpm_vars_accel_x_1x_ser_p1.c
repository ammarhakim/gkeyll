#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p1(const double *dxv, 
  const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
  const double *prim_c, const double *bvar_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // prim_surf_l/c/r: Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // prim_c:          Input volume expansion of primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp] in center cell.
  // bvar_c:          Input volume expansion of magnetic field unit vector and tensor in center cell.
  // nu_c:            Input volume expansion of collisionality in center cell.
  // pkpm_accel:      Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[0]; 
  const double *ux_c = &prim_c[0]; 
  const double *uy_c = &prim_c[2]; 
  const double *uz_c = &prim_c[4]; 

  const double *bxbx = &bvar_c[6]; 
  const double *bxby = &bvar_c[8]; 
  const double *bxbz = &bvar_c[10]; 
  const double *byby = &bvar_c[12]; 
  const double *bybz = &bvar_c[14]; 
  const double *bzbz = &bvar_c[16]; 

  const double *ux_surf_lr = &prim_surf_l[1]; 
  const double *uy_surf_lr = &prim_surf_l[3]; 
  const double *uz_surf_lr = &prim_surf_l[5]; 

  const double *ux_surf_cl = &prim_surf_c[0]; 
  const double *uy_surf_cl = &prim_surf_c[2]; 
  const double *uz_surf_cl = &prim_surf_c[4]; 

  const double *ux_surf_cr = &prim_surf_c[1]; 
  const double *uy_surf_cr = &prim_surf_c[3]; 
  const double *uz_surf_cr = &prim_surf_c[5]; 

  const double *ux_surf_rl = &prim_surf_r[0]; 
  const double *uy_surf_rl = &prim_surf_r[2]; 
  const double *uz_surf_rl = &prim_surf_r[4]; 

  double *bb_grad_u = &pkpm_accel[2]; 
  double *p_perp_source = &pkpm_accel[6]; 


  double grad_u_x[2] = {0.0}; 
  double grad_u_y[2] = {0.0}; 
  double grad_u_z[2] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 

  double bb_grad_u_comp[2] = {0.0}; 
  bb_grad_u_comp[0] = 0.7071067811865475*bxbz[1]*grad_u_z[1]+0.7071067811865475*bxby[1]*grad_u_y[1]+0.7071067811865475*bxbx[1]*grad_u_x[1]+0.7071067811865475*bxbz[0]*grad_u_z[0]+0.7071067811865475*bxby[0]*grad_u_y[0]+0.7071067811865475*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.7071067811865475*bxbz[0]*grad_u_z[1]+0.7071067811865475*bxby[0]*grad_u_y[1]+0.7071067811865475*bxbx[0]*grad_u_x[1]+0.7071067811865475*grad_u_z[0]*bxbz[1]+0.7071067811865475*grad_u_y[0]*bxby[1]+0.7071067811865475*grad_u_x[0]*bxbx[1]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 

  p_perp_source[0] += (-2.0*nu_c[0])-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-2.0*nu_c[1])-1.0*grad_u_x[1]+bb_grad_u_comp[1]; 

} 
