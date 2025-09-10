#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p1(const double *dxv, 
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

  const double dx1 = 2.0/dxv[1]; 
  const double *ux_c = &prim_c[0]; 
  const double *uy_c = &prim_c[8]; 
  const double *uz_c = &prim_c[16]; 

  const double *bxbx = &bvar_c[24]; 
  const double *bxby = &bvar_c[32]; 
  const double *bxbz = &bvar_c[40]; 
  const double *byby = &bvar_c[48]; 
  const double *bybz = &bvar_c[56]; 
  const double *bzbz = &bvar_c[64]; 

  const double *ux_surf_lr = &prim_surf_l[36]; 
  const double *uy_surf_lr = &prim_surf_l[44]; 
  const double *uz_surf_lr = &prim_surf_l[52]; 

  const double *ux_surf_cl = &prim_surf_c[32]; 
  const double *uy_surf_cl = &prim_surf_c[40]; 
  const double *uz_surf_cl = &prim_surf_c[48]; 

  const double *ux_surf_cr = &prim_surf_c[36]; 
  const double *uy_surf_cr = &prim_surf_c[44]; 
  const double *uz_surf_cr = &prim_surf_c[52]; 

  const double *ux_surf_rl = &prim_surf_r[32]; 
  const double *uy_surf_rl = &prim_surf_r[40]; 
  const double *uz_surf_rl = &prim_surf_r[48]; 

  double *bb_grad_u = &pkpm_accel[8]; 
  double *p_perp_source = &pkpm_accel[24]; 


  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.3535533905932737*ux_surf_rl[1]-0.3535533905932737*ux_surf_lr[1]+0.3535533905932737*ux_surf_cr[1]-0.3535533905932737*ux_surf_cl[1])*dx1; 
  grad_u_x[2] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[3] = (0.3535533905932737*ux_surf_rl[2]-0.3535533905932737*ux_surf_lr[2]+0.3535533905932737*ux_surf_cr[2]-0.3535533905932737*ux_surf_cl[2])*dx1; 
  grad_u_x[4] = (0.6123724356957944*(ux_surf_rl[1]+ux_surf_lr[1]+ux_surf_cr[1]+ux_surf_cl[1])-1.732050807568877*ux_c[1])*dx1; 
  grad_u_x[5] = (0.3535533905932737*ux_surf_rl[3]-0.3535533905932737*ux_surf_lr[3]+0.3535533905932737*ux_surf_cr[3]-0.3535533905932737*ux_surf_cl[3])*dx1; 
  grad_u_x[6] = (0.6123724356957944*(ux_surf_rl[2]+ux_surf_lr[2]+ux_surf_cr[2]+ux_surf_cl[2])-1.732050807568877*ux_c[3])*dx1; 
  grad_u_x[7] = (0.6123724356957944*(ux_surf_rl[3]+ux_surf_lr[3]+ux_surf_cr[3]+ux_surf_cl[3])-1.732050807568877*ux_c[5])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.3535533905932737*uy_surf_rl[1]-0.3535533905932737*uy_surf_lr[1]+0.3535533905932737*uy_surf_cr[1]-0.3535533905932737*uy_surf_cl[1])*dx1; 
  grad_u_y[2] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[3] = (0.3535533905932737*uy_surf_rl[2]-0.3535533905932737*uy_surf_lr[2]+0.3535533905932737*uy_surf_cr[2]-0.3535533905932737*uy_surf_cl[2])*dx1; 
  grad_u_y[4] = (0.6123724356957944*(uy_surf_rl[1]+uy_surf_lr[1]+uy_surf_cr[1]+uy_surf_cl[1])-1.732050807568877*uy_c[1])*dx1; 
  grad_u_y[5] = (0.3535533905932737*uy_surf_rl[3]-0.3535533905932737*uy_surf_lr[3]+0.3535533905932737*uy_surf_cr[3]-0.3535533905932737*uy_surf_cl[3])*dx1; 
  grad_u_y[6] = (0.6123724356957944*(uy_surf_rl[2]+uy_surf_lr[2]+uy_surf_cr[2]+uy_surf_cl[2])-1.732050807568877*uy_c[3])*dx1; 
  grad_u_y[7] = (0.6123724356957944*(uy_surf_rl[3]+uy_surf_lr[3]+uy_surf_cr[3]+uy_surf_cl[3])-1.732050807568877*uy_c[5])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.3535533905932737*uz_surf_rl[1]-0.3535533905932737*uz_surf_lr[1]+0.3535533905932737*uz_surf_cr[1]-0.3535533905932737*uz_surf_cl[1])*dx1; 
  grad_u_z[2] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[3] = (0.3535533905932737*uz_surf_rl[2]-0.3535533905932737*uz_surf_lr[2]+0.3535533905932737*uz_surf_cr[2]-0.3535533905932737*uz_surf_cl[2])*dx1; 
  grad_u_z[4] = (0.6123724356957944*(uz_surf_rl[1]+uz_surf_lr[1]+uz_surf_cr[1]+uz_surf_cl[1])-1.732050807568877*uz_c[1])*dx1; 
  grad_u_z[5] = (0.3535533905932737*uz_surf_rl[3]-0.3535533905932737*uz_surf_lr[3]+0.3535533905932737*uz_surf_cr[3]-0.3535533905932737*uz_surf_cl[3])*dx1; 
  grad_u_z[6] = (0.6123724356957944*(uz_surf_rl[2]+uz_surf_lr[2]+uz_surf_cr[2]+uz_surf_cl[2])-1.732050807568877*uz_c[3])*dx1; 
  grad_u_z[7] = (0.6123724356957944*(uz_surf_rl[3]+uz_surf_lr[3]+uz_surf_cr[3]+uz_surf_cl[3])-1.732050807568877*uz_c[5])*dx1; 

  double bb_grad_u_comp[8] = {0.0}; 
  bb_grad_u_comp[0] = 0.3535533905932737*bybz[7]*grad_u_z[7]+0.3535533905932737*byby[7]*grad_u_y[7]+0.3535533905932737*bxby[7]*grad_u_x[7]+0.3535533905932737*bybz[6]*grad_u_z[6]+0.3535533905932737*byby[6]*grad_u_y[6]+0.3535533905932737*bxby[6]*grad_u_x[6]+0.3535533905932737*bybz[5]*grad_u_z[5]+0.3535533905932737*byby[5]*grad_u_y[5]+0.3535533905932737*bxby[5]*grad_u_x[5]+0.3535533905932737*bybz[4]*grad_u_z[4]+0.3535533905932737*byby[4]*grad_u_y[4]+0.3535533905932737*bxby[4]*grad_u_x[4]+0.3535533905932737*bybz[3]*grad_u_z[3]+0.3535533905932737*byby[3]*grad_u_y[3]+0.3535533905932737*bxby[3]*grad_u_x[3]+0.3535533905932737*bybz[2]*grad_u_z[2]+0.3535533905932737*byby[2]*grad_u_y[2]+0.3535533905932737*bxby[2]*grad_u_x[2]+0.3535533905932737*bybz[1]*grad_u_z[1]+0.3535533905932737*byby[1]*grad_u_y[1]+0.3535533905932737*bxby[1]*grad_u_x[1]+0.3535533905932737*bybz[0]*grad_u_z[0]+0.3535533905932737*byby[0]*grad_u_y[0]+0.3535533905932737*bxby[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.3535533905932737*bybz[6]*grad_u_z[7]+0.3535533905932737*byby[6]*grad_u_y[7]+0.3535533905932737*bxby[6]*grad_u_x[7]+0.3535533905932737*grad_u_z[6]*bybz[7]+0.3535533905932737*grad_u_y[6]*byby[7]+0.3535533905932737*grad_u_x[6]*bxby[7]+0.3535533905932737*bybz[3]*grad_u_z[5]+0.3535533905932737*byby[3]*grad_u_y[5]+0.3535533905932737*bxby[3]*grad_u_x[5]+0.3535533905932737*grad_u_z[3]*bybz[5]+0.3535533905932737*grad_u_y[3]*byby[5]+0.3535533905932737*grad_u_x[3]*bxby[5]+0.3535533905932737*bybz[2]*grad_u_z[4]+0.3535533905932737*byby[2]*grad_u_y[4]+0.3535533905932737*bxby[2]*grad_u_x[4]+0.3535533905932737*grad_u_z[2]*bybz[4]+0.3535533905932737*grad_u_y[2]*byby[4]+0.3535533905932737*grad_u_x[2]*bxby[4]+0.3535533905932737*bybz[0]*grad_u_z[1]+0.3535533905932737*byby[0]*grad_u_y[1]+0.3535533905932737*bxby[0]*grad_u_x[1]+0.3535533905932737*grad_u_z[0]*bybz[1]+0.3535533905932737*grad_u_y[0]*byby[1]+0.3535533905932737*grad_u_x[0]*bxby[1]; 
  bb_grad_u_comp[2] = 0.3535533905932737*bybz[5]*grad_u_z[7]+0.3535533905932737*byby[5]*grad_u_y[7]+0.3535533905932737*bxby[5]*grad_u_x[7]+0.3535533905932737*grad_u_z[5]*bybz[7]+0.3535533905932737*grad_u_y[5]*byby[7]+0.3535533905932737*grad_u_x[5]*bxby[7]+0.3535533905932737*bybz[3]*grad_u_z[6]+0.3535533905932737*byby[3]*grad_u_y[6]+0.3535533905932737*bxby[3]*grad_u_x[6]+0.3535533905932737*grad_u_z[3]*bybz[6]+0.3535533905932737*grad_u_y[3]*byby[6]+0.3535533905932737*grad_u_x[3]*bxby[6]+0.3535533905932737*bybz[1]*grad_u_z[4]+0.3535533905932737*byby[1]*grad_u_y[4]+0.3535533905932737*bxby[1]*grad_u_x[4]+0.3535533905932737*grad_u_z[1]*bybz[4]+0.3535533905932737*grad_u_y[1]*byby[4]+0.3535533905932737*grad_u_x[1]*bxby[4]+0.3535533905932737*bybz[0]*grad_u_z[2]+0.3535533905932737*byby[0]*grad_u_y[2]+0.3535533905932737*bxby[0]*grad_u_x[2]+0.3535533905932737*grad_u_z[0]*bybz[2]+0.3535533905932737*grad_u_y[0]*byby[2]+0.3535533905932737*grad_u_x[0]*bxby[2]; 
  bb_grad_u_comp[3] = 0.3535533905932737*bybz[4]*grad_u_z[7]+0.3535533905932737*byby[4]*grad_u_y[7]+0.3535533905932737*bxby[4]*grad_u_x[7]+0.3535533905932737*grad_u_z[4]*bybz[7]+0.3535533905932737*grad_u_y[4]*byby[7]+0.3535533905932737*grad_u_x[4]*bxby[7]+0.3535533905932737*bybz[2]*grad_u_z[6]+0.3535533905932737*byby[2]*grad_u_y[6]+0.3535533905932737*bxby[2]*grad_u_x[6]+0.3535533905932737*grad_u_z[2]*bybz[6]+0.3535533905932737*grad_u_y[2]*byby[6]+0.3535533905932737*grad_u_x[2]*bxby[6]+0.3535533905932737*bybz[1]*grad_u_z[5]+0.3535533905932737*byby[1]*grad_u_y[5]+0.3535533905932737*bxby[1]*grad_u_x[5]+0.3535533905932737*grad_u_z[1]*bybz[5]+0.3535533905932737*grad_u_y[1]*byby[5]+0.3535533905932737*grad_u_x[1]*bxby[5]+0.3535533905932737*bybz[0]*grad_u_z[3]+0.3535533905932737*byby[0]*grad_u_y[3]+0.3535533905932737*bxby[0]*grad_u_x[3]+0.3535533905932737*grad_u_z[0]*bybz[3]+0.3535533905932737*grad_u_y[0]*byby[3]+0.3535533905932737*grad_u_x[0]*bxby[3]; 
  bb_grad_u_comp[4] = 0.3535533905932737*bybz[3]*grad_u_z[7]+0.3535533905932737*byby[3]*grad_u_y[7]+0.3535533905932737*bxby[3]*grad_u_x[7]+0.3535533905932737*grad_u_z[3]*bybz[7]+0.3535533905932737*grad_u_y[3]*byby[7]+0.3535533905932737*grad_u_x[3]*bxby[7]+0.3535533905932737*bybz[5]*grad_u_z[6]+0.3535533905932737*byby[5]*grad_u_y[6]+0.3535533905932737*bxby[5]*grad_u_x[6]+0.3535533905932737*grad_u_z[5]*bybz[6]+0.3535533905932737*grad_u_y[5]*byby[6]+0.3535533905932737*grad_u_x[5]*bxby[6]+0.3535533905932737*bybz[0]*grad_u_z[4]+0.3535533905932737*byby[0]*grad_u_y[4]+0.3535533905932737*bxby[0]*grad_u_x[4]+0.3535533905932737*grad_u_z[0]*bybz[4]+0.3535533905932737*grad_u_y[0]*byby[4]+0.3535533905932737*grad_u_x[0]*bxby[4]+0.3535533905932737*bybz[1]*grad_u_z[2]+0.3535533905932737*byby[1]*grad_u_y[2]+0.3535533905932737*bxby[1]*grad_u_x[2]+0.3535533905932737*grad_u_z[1]*bybz[2]+0.3535533905932737*grad_u_y[1]*byby[2]+0.3535533905932737*grad_u_x[1]*bxby[2]; 
  bb_grad_u_comp[5] = 0.3535533905932737*bybz[2]*grad_u_z[7]+0.3535533905932737*byby[2]*grad_u_y[7]+0.3535533905932737*bxby[2]*grad_u_x[7]+0.3535533905932737*grad_u_z[2]*bybz[7]+0.3535533905932737*grad_u_y[2]*byby[7]+0.3535533905932737*grad_u_x[2]*bxby[7]+0.3535533905932737*bybz[4]*grad_u_z[6]+0.3535533905932737*byby[4]*grad_u_y[6]+0.3535533905932737*bxby[4]*grad_u_x[6]+0.3535533905932737*grad_u_z[4]*bybz[6]+0.3535533905932737*grad_u_y[4]*byby[6]+0.3535533905932737*grad_u_x[4]*bxby[6]+0.3535533905932737*bybz[0]*grad_u_z[5]+0.3535533905932737*byby[0]*grad_u_y[5]+0.3535533905932737*bxby[0]*grad_u_x[5]+0.3535533905932737*grad_u_z[0]*bybz[5]+0.3535533905932737*grad_u_y[0]*byby[5]+0.3535533905932737*grad_u_x[0]*bxby[5]+0.3535533905932737*bybz[1]*grad_u_z[3]+0.3535533905932737*byby[1]*grad_u_y[3]+0.3535533905932737*bxby[1]*grad_u_x[3]+0.3535533905932737*grad_u_z[1]*bybz[3]+0.3535533905932737*grad_u_y[1]*byby[3]+0.3535533905932737*grad_u_x[1]*bxby[3]; 
  bb_grad_u_comp[6] = 0.3535533905932737*bybz[1]*grad_u_z[7]+0.3535533905932737*byby[1]*grad_u_y[7]+0.3535533905932737*bxby[1]*grad_u_x[7]+0.3535533905932737*grad_u_z[1]*bybz[7]+0.3535533905932737*grad_u_y[1]*byby[7]+0.3535533905932737*grad_u_x[1]*bxby[7]+0.3535533905932737*bybz[0]*grad_u_z[6]+0.3535533905932737*byby[0]*grad_u_y[6]+0.3535533905932737*bxby[0]*grad_u_x[6]+0.3535533905932737*grad_u_z[0]*bybz[6]+0.3535533905932737*grad_u_y[0]*byby[6]+0.3535533905932737*grad_u_x[0]*bxby[6]+0.3535533905932737*bybz[4]*grad_u_z[5]+0.3535533905932737*byby[4]*grad_u_y[5]+0.3535533905932737*bxby[4]*grad_u_x[5]+0.3535533905932737*grad_u_z[4]*bybz[5]+0.3535533905932737*grad_u_y[4]*byby[5]+0.3535533905932737*grad_u_x[4]*bxby[5]+0.3535533905932737*bybz[2]*grad_u_z[3]+0.3535533905932737*byby[2]*grad_u_y[3]+0.3535533905932737*bxby[2]*grad_u_x[3]+0.3535533905932737*grad_u_z[2]*bybz[3]+0.3535533905932737*grad_u_y[2]*byby[3]+0.3535533905932737*grad_u_x[2]*bxby[3]; 
  bb_grad_u_comp[7] = 0.3535533905932737*bybz[0]*grad_u_z[7]+0.3535533905932737*byby[0]*grad_u_y[7]+0.3535533905932737*bxby[0]*grad_u_x[7]+0.3535533905932737*grad_u_z[0]*bybz[7]+0.3535533905932737*grad_u_y[0]*byby[7]+0.3535533905932737*grad_u_x[0]*bxby[7]+0.3535533905932737*bybz[1]*grad_u_z[6]+0.3535533905932737*byby[1]*grad_u_y[6]+0.3535533905932737*bxby[1]*grad_u_x[6]+0.3535533905932737*grad_u_z[1]*bybz[6]+0.3535533905932737*grad_u_y[1]*byby[6]+0.3535533905932737*grad_u_x[1]*bxby[6]+0.3535533905932737*bybz[2]*grad_u_z[5]+0.3535533905932737*byby[2]*grad_u_y[5]+0.3535533905932737*bxby[2]*grad_u_x[5]+0.3535533905932737*grad_u_z[2]*bybz[5]+0.3535533905932737*grad_u_y[2]*byby[5]+0.3535533905932737*grad_u_x[2]*bxby[5]+0.3535533905932737*bybz[3]*grad_u_z[4]+0.3535533905932737*byby[3]*grad_u_y[4]+0.3535533905932737*bxby[3]*grad_u_x[4]+0.3535533905932737*grad_u_z[3]*bybz[4]+0.3535533905932737*grad_u_y[3]*byby[4]+0.3535533905932737*grad_u_x[3]*bxby[4]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 

  p_perp_source[0] += (-0.6666666666666666*nu_c[0])-1.0*grad_u_y[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-0.6666666666666666*nu_c[1])-1.0*grad_u_y[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-0.6666666666666666*nu_c[2])-1.0*grad_u_y[2]+bb_grad_u_comp[2]; 
  p_perp_source[3] += (-0.6666666666666666*nu_c[3])-1.0*grad_u_y[3]+bb_grad_u_comp[3]; 
  p_perp_source[4] += (-0.6666666666666666*nu_c[4])-1.0*grad_u_y[4]+bb_grad_u_comp[4]; 
  p_perp_source[5] += (-0.6666666666666666*nu_c[5])-1.0*grad_u_y[5]+bb_grad_u_comp[5]; 
  p_perp_source[6] += (-0.6666666666666666*nu_c[6])-1.0*grad_u_y[6]+bb_grad_u_comp[6]; 
  p_perp_source[7] += (-0.6666666666666666*nu_c[7])-1.0*grad_u_y[7]+bb_grad_u_comp[7]; 

} 
