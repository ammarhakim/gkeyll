#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_x_2x_tensor_p2(const double *dxv, 
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
  const double *uy_c = &prim_c[9]; 
  const double *uz_c = &prim_c[18]; 

  const double *bxbx = &bvar_c[27]; 
  const double *bxby = &bvar_c[36]; 
  const double *bxbz = &bvar_c[45]; 
  const double *byby = &bvar_c[54]; 
  const double *bybz = &bvar_c[63]; 
  const double *bzbz = &bvar_c[72]; 

  const double *ux_surf_lr = &prim_surf_l[3]; 
  const double *uy_surf_lr = &prim_surf_l[9]; 
  const double *uz_surf_lr = &prim_surf_l[15]; 

  const double *ux_surf_cl = &prim_surf_c[0]; 
  const double *uy_surf_cl = &prim_surf_c[6]; 
  const double *uz_surf_cl = &prim_surf_c[12]; 

  const double *ux_surf_cr = &prim_surf_c[3]; 
  const double *uy_surf_cr = &prim_surf_c[9]; 
  const double *uz_surf_cr = &prim_surf_c[15]; 

  const double *ux_surf_rl = &prim_surf_r[0]; 
  const double *uy_surf_rl = &prim_surf_r[6]; 
  const double *uz_surf_rl = &prim_surf_r[12]; 

  double *bb_grad_u = &pkpm_accel[9]; 
  double *p_perp_source = &pkpm_accel[27]; 


  double grad_u_x[9] = {0.0}; 
  double grad_u_y[9] = {0.0}; 
  double grad_u_z[9] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[2] = (0.3535533905932737*ux_surf_rl[1]-0.3535533905932737*ux_surf_lr[1]+0.3535533905932737*ux_surf_cr[1]-0.3535533905932737*ux_surf_cl[1])*dx1; 
  grad_u_x[3] = (0.6123724356957944*(ux_surf_rl[1]+ux_surf_lr[1]+ux_surf_cr[1]+ux_surf_cl[1])-1.732050807568877*ux_c[2])*dx1; 
  grad_u_x[4] = ((-3.872983346207417*ux_c[1])+0.7905694150420947*ux_surf_rl[0]-0.7905694150420947*ux_surf_lr[0]+0.7905694150420947*ux_surf_cr[0]-0.7905694150420947*ux_surf_cl[0])*dx1; 
  grad_u_x[5] = (0.3535533905932737*ux_surf_rl[2]-0.3535533905932737*ux_surf_lr[2]+0.3535533905932737*ux_surf_cr[2]-0.3535533905932737*ux_surf_cl[2])*dx1; 
  grad_u_x[6] = ((-3.872983346207417*ux_c[3])+0.7905694150420948*ux_surf_rl[1]-0.7905694150420948*ux_surf_lr[1]+0.7905694150420948*ux_surf_cr[1]-0.7905694150420948*ux_surf_cl[1])*dx1; 
  grad_u_x[7] = (0.6123724356957944*(ux_surf_rl[2]+ux_surf_lr[2]+ux_surf_cr[2]+ux_surf_cl[2])-1.732050807568877*ux_c[5])*dx1; 
  grad_u_x[8] = ((-3.872983346207417*ux_c[7])+0.7905694150420947*ux_surf_rl[2]-0.7905694150420947*ux_surf_lr[2]+0.7905694150420947*ux_surf_cr[2]-0.7905694150420947*ux_surf_cl[2])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[2] = (0.3535533905932737*uy_surf_rl[1]-0.3535533905932737*uy_surf_lr[1]+0.3535533905932737*uy_surf_cr[1]-0.3535533905932737*uy_surf_cl[1])*dx1; 
  grad_u_y[3] = (0.6123724356957944*(uy_surf_rl[1]+uy_surf_lr[1]+uy_surf_cr[1]+uy_surf_cl[1])-1.732050807568877*uy_c[2])*dx1; 
  grad_u_y[4] = ((-3.872983346207417*uy_c[1])+0.7905694150420947*uy_surf_rl[0]-0.7905694150420947*uy_surf_lr[0]+0.7905694150420947*uy_surf_cr[0]-0.7905694150420947*uy_surf_cl[0])*dx1; 
  grad_u_y[5] = (0.3535533905932737*uy_surf_rl[2]-0.3535533905932737*uy_surf_lr[2]+0.3535533905932737*uy_surf_cr[2]-0.3535533905932737*uy_surf_cl[2])*dx1; 
  grad_u_y[6] = ((-3.872983346207417*uy_c[3])+0.7905694150420948*uy_surf_rl[1]-0.7905694150420948*uy_surf_lr[1]+0.7905694150420948*uy_surf_cr[1]-0.7905694150420948*uy_surf_cl[1])*dx1; 
  grad_u_y[7] = (0.6123724356957944*(uy_surf_rl[2]+uy_surf_lr[2]+uy_surf_cr[2]+uy_surf_cl[2])-1.732050807568877*uy_c[5])*dx1; 
  grad_u_y[8] = ((-3.872983346207417*uy_c[7])+0.7905694150420947*uy_surf_rl[2]-0.7905694150420947*uy_surf_lr[2]+0.7905694150420947*uy_surf_cr[2]-0.7905694150420947*uy_surf_cl[2])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[2] = (0.3535533905932737*uz_surf_rl[1]-0.3535533905932737*uz_surf_lr[1]+0.3535533905932737*uz_surf_cr[1]-0.3535533905932737*uz_surf_cl[1])*dx1; 
  grad_u_z[3] = (0.6123724356957944*(uz_surf_rl[1]+uz_surf_lr[1]+uz_surf_cr[1]+uz_surf_cl[1])-1.732050807568877*uz_c[2])*dx1; 
  grad_u_z[4] = ((-3.872983346207417*uz_c[1])+0.7905694150420947*uz_surf_rl[0]-0.7905694150420947*uz_surf_lr[0]+0.7905694150420947*uz_surf_cr[0]-0.7905694150420947*uz_surf_cl[0])*dx1; 
  grad_u_z[5] = (0.3535533905932737*uz_surf_rl[2]-0.3535533905932737*uz_surf_lr[2]+0.3535533905932737*uz_surf_cr[2]-0.3535533905932737*uz_surf_cl[2])*dx1; 
  grad_u_z[6] = ((-3.872983346207417*uz_c[3])+0.7905694150420948*uz_surf_rl[1]-0.7905694150420948*uz_surf_lr[1]+0.7905694150420948*uz_surf_cr[1]-0.7905694150420948*uz_surf_cl[1])*dx1; 
  grad_u_z[7] = (0.6123724356957944*(uz_surf_rl[2]+uz_surf_lr[2]+uz_surf_cr[2]+uz_surf_cl[2])-1.732050807568877*uz_c[5])*dx1; 
  grad_u_z[8] = ((-3.872983346207417*uz_c[7])+0.7905694150420947*uz_surf_rl[2]-0.7905694150420947*uz_surf_lr[2]+0.7905694150420947*uz_surf_cr[2]-0.7905694150420947*uz_surf_cl[2])*dx1; 

  double bb_grad_u_comp[9] = {0.0}; 
  bb_grad_u_comp[0] = 0.5*bxbz[8]*grad_u_z[8]+0.5*bxby[8]*grad_u_y[8]+0.5*bxbx[8]*grad_u_x[8]+0.5*bxbz[7]*grad_u_z[7]+0.5*bxby[7]*grad_u_y[7]+0.5*bxbx[7]*grad_u_x[7]+0.5*bxbz[6]*grad_u_z[6]+0.5*bxby[6]*grad_u_y[6]+0.5*bxbx[6]*grad_u_x[6]+0.5*bxbz[5]*grad_u_z[5]+0.5*bxby[5]*grad_u_y[5]+0.5*bxbx[5]*grad_u_x[5]+0.5*bxbz[4]*grad_u_z[4]+0.5*bxby[4]*grad_u_y[4]+0.5*bxbx[4]*grad_u_x[4]+0.5*bxbz[3]*grad_u_z[3]+0.5*bxby[3]*grad_u_y[3]+0.5*bxbx[3]*grad_u_x[3]+0.5*bxbz[2]*grad_u_z[2]+0.5*bxby[2]*grad_u_y[2]+0.5*bxbx[2]*grad_u_x[2]+0.5*bxbz[1]*grad_u_z[1]+0.5*bxby[1]*grad_u_y[1]+0.5*bxbx[1]*grad_u_x[1]+0.5*bxbz[0]*grad_u_z[0]+0.5*bxby[0]*grad_u_y[0]+0.5*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.447213595499958*bxbz[7]*grad_u_z[8]+0.447213595499958*bxby[7]*grad_u_y[8]+0.447213595499958*bxbx[7]*grad_u_x[8]+0.447213595499958*grad_u_z[7]*bxbz[8]+0.447213595499958*grad_u_y[7]*bxby[8]+0.447213595499958*grad_u_x[7]*bxbx[8]+0.5000000000000001*bxbz[5]*grad_u_z[7]+0.5000000000000001*bxby[5]*grad_u_y[7]+0.5000000000000001*bxbx[5]*grad_u_x[7]+0.5000000000000001*grad_u_z[5]*bxbz[7]+0.5000000000000001*grad_u_y[5]*bxby[7]+0.5000000000000001*grad_u_x[5]*bxbx[7]+0.447213595499958*bxbz[3]*grad_u_z[6]+0.447213595499958*bxby[3]*grad_u_y[6]+0.447213595499958*bxbx[3]*grad_u_x[6]+0.447213595499958*grad_u_z[3]*bxbz[6]+0.447213595499958*grad_u_y[3]*bxby[6]+0.447213595499958*grad_u_x[3]*bxbx[6]+0.4472135954999579*bxbz[1]*grad_u_z[4]+0.4472135954999579*bxby[1]*grad_u_y[4]+0.4472135954999579*bxbx[1]*grad_u_x[4]+0.4472135954999579*grad_u_z[1]*bxbz[4]+0.4472135954999579*grad_u_y[1]*bxby[4]+0.4472135954999579*grad_u_x[1]*bxbx[4]+0.5*bxbz[2]*grad_u_z[3]+0.5*bxby[2]*grad_u_y[3]+0.5*bxbx[2]*grad_u_x[3]+0.5*grad_u_z[2]*bxbz[3]+0.5*grad_u_y[2]*bxby[3]+0.5*grad_u_x[2]*bxbx[3]+0.5*bxbz[0]*grad_u_z[1]+0.5*bxby[0]*grad_u_y[1]+0.5*bxbx[0]*grad_u_x[1]+0.5*grad_u_z[0]*bxbz[1]+0.5*grad_u_y[0]*bxby[1]+0.5*grad_u_x[0]*bxbx[1]; 
  bb_grad_u_comp[2] = 0.447213595499958*bxbz[6]*grad_u_z[8]+0.447213595499958*bxby[6]*grad_u_y[8]+0.447213595499958*bxbx[6]*grad_u_x[8]+0.447213595499958*grad_u_z[6]*bxbz[8]+0.447213595499958*grad_u_y[6]*bxby[8]+0.447213595499958*grad_u_x[6]*bxbx[8]+0.447213595499958*bxbz[3]*grad_u_z[7]+0.447213595499958*bxby[3]*grad_u_y[7]+0.447213595499958*bxbx[3]*grad_u_x[7]+0.447213595499958*grad_u_z[3]*bxbz[7]+0.447213595499958*grad_u_y[3]*bxby[7]+0.447213595499958*grad_u_x[3]*bxbx[7]+0.5000000000000001*bxbz[4]*grad_u_z[6]+0.5000000000000001*bxby[4]*grad_u_y[6]+0.5000000000000001*bxbx[4]*grad_u_x[6]+0.5000000000000001*grad_u_z[4]*bxbz[6]+0.5000000000000001*grad_u_y[4]*bxby[6]+0.5000000000000001*grad_u_x[4]*bxbx[6]+0.4472135954999579*bxbz[2]*grad_u_z[5]+0.4472135954999579*bxby[2]*grad_u_y[5]+0.4472135954999579*bxbx[2]*grad_u_x[5]+0.4472135954999579*grad_u_z[2]*bxbz[5]+0.4472135954999579*grad_u_y[2]*bxby[5]+0.4472135954999579*grad_u_x[2]*bxbx[5]+0.5*bxbz[1]*grad_u_z[3]+0.5*bxby[1]*grad_u_y[3]+0.5*bxbx[1]*grad_u_x[3]+0.5*grad_u_z[1]*bxbz[3]+0.5*grad_u_y[1]*bxby[3]+0.5*grad_u_x[1]*bxbx[3]+0.5*bxbz[0]*grad_u_z[2]+0.5*bxby[0]*grad_u_y[2]+0.5*bxbx[0]*grad_u_x[2]+0.5*grad_u_z[0]*bxbz[2]+0.5*grad_u_y[0]*bxby[2]+0.5*grad_u_x[0]*bxbx[2]; 
  bb_grad_u_comp[3] = 0.4*bxbz[3]*grad_u_z[8]+0.4*bxby[3]*grad_u_y[8]+0.4*bxbx[3]*grad_u_x[8]+0.4*grad_u_z[3]*bxbz[8]+0.4*grad_u_y[3]*bxby[8]+0.4*grad_u_x[3]*bxbx[8]+0.4*bxbz[6]*grad_u_z[7]+0.447213595499958*bxbz[2]*grad_u_z[7]+0.4*bxby[6]*grad_u_y[7]+0.447213595499958*bxby[2]*grad_u_y[7]+0.4*bxbx[6]*grad_u_x[7]+0.447213595499958*bxbx[2]*grad_u_x[7]+0.4*grad_u_z[6]*bxbz[7]+0.447213595499958*grad_u_z[2]*bxbz[7]+0.4*grad_u_y[6]*bxby[7]+0.447213595499958*grad_u_y[2]*bxby[7]+0.4*grad_u_x[6]*bxbx[7]+0.447213595499958*grad_u_x[2]*bxbx[7]+0.447213595499958*bxbz[1]*grad_u_z[6]+0.447213595499958*bxby[1]*grad_u_y[6]+0.447213595499958*bxbx[1]*grad_u_x[6]+0.447213595499958*grad_u_z[1]*bxbz[6]+0.447213595499958*grad_u_y[1]*bxby[6]+0.447213595499958*grad_u_x[1]*bxbx[6]+0.4472135954999579*bxbz[3]*grad_u_z[5]+0.4472135954999579*bxby[3]*grad_u_y[5]+0.4472135954999579*bxbx[3]*grad_u_x[5]+0.4472135954999579*grad_u_z[3]*bxbz[5]+0.4472135954999579*grad_u_y[3]*bxby[5]+0.4472135954999579*grad_u_x[3]*bxbx[5]+0.4472135954999579*bxbz[3]*grad_u_z[4]+0.4472135954999579*bxby[3]*grad_u_y[4]+0.4472135954999579*bxbx[3]*grad_u_x[4]+0.4472135954999579*grad_u_z[3]*bxbz[4]+0.4472135954999579*grad_u_y[3]*bxby[4]+0.4472135954999579*grad_u_x[3]*bxbx[4]+0.5*bxbz[0]*grad_u_z[3]+0.5*bxby[0]*grad_u_y[3]+0.5*bxbx[0]*grad_u_x[3]+0.5*grad_u_z[0]*bxbz[3]+0.5*grad_u_y[0]*bxby[3]+0.5*grad_u_x[0]*bxbx[3]+0.5*bxbz[1]*grad_u_z[2]+0.5*bxby[1]*grad_u_y[2]+0.5*bxbx[1]*grad_u_x[2]+0.5*grad_u_z[1]*bxbz[2]+0.5*grad_u_y[1]*bxby[2]+0.5*grad_u_x[1]*bxbx[2]; 
  bb_grad_u_comp[4] = 0.31943828249997*bxbz[8]*grad_u_z[8]+0.5*bxbz[5]*grad_u_z[8]+0.31943828249997*bxby[8]*grad_u_y[8]+0.5*bxby[5]*grad_u_y[8]+0.31943828249997*bxbx[8]*grad_u_x[8]+0.5*bxbx[5]*grad_u_x[8]+0.5*grad_u_z[5]*bxbz[8]+0.5*grad_u_y[5]*bxby[8]+0.5*grad_u_x[5]*bxbx[8]+0.4472135954999579*bxbz[7]*grad_u_z[7]+0.4472135954999579*bxby[7]*grad_u_y[7]+0.4472135954999579*bxbx[7]*grad_u_x[7]+0.31943828249997*bxbz[6]*grad_u_z[6]+0.5000000000000001*bxbz[2]*grad_u_z[6]+0.31943828249997*bxby[6]*grad_u_y[6]+0.5000000000000001*bxby[2]*grad_u_y[6]+0.31943828249997*bxbx[6]*grad_u_x[6]+0.5000000000000001*bxbx[2]*grad_u_x[6]+0.5000000000000001*grad_u_z[2]*bxbz[6]+0.5000000000000001*grad_u_y[2]*bxby[6]+0.5000000000000001*grad_u_x[2]*bxbx[6]+0.31943828249997*bxbz[4]*grad_u_z[4]+0.5*bxbz[0]*grad_u_z[4]+0.31943828249997*bxby[4]*grad_u_y[4]+0.5*bxby[0]*grad_u_y[4]+0.31943828249997*bxbx[4]*grad_u_x[4]+0.5*bxbx[0]*grad_u_x[4]+0.5*grad_u_z[0]*bxbz[4]+0.5*grad_u_y[0]*bxby[4]+0.5*grad_u_x[0]*bxbx[4]+0.4472135954999579*bxbz[3]*grad_u_z[3]+0.4472135954999579*bxby[3]*grad_u_y[3]+0.4472135954999579*bxbx[3]*grad_u_x[3]+0.4472135954999579*bxbz[1]*grad_u_z[1]+0.4472135954999579*bxby[1]*grad_u_y[1]+0.4472135954999579*bxbx[1]*grad_u_x[1]; 
  bb_grad_u_comp[5] = 0.31943828249997*bxbz[8]*grad_u_z[8]+0.5*bxbz[4]*grad_u_z[8]+0.31943828249997*bxby[8]*grad_u_y[8]+0.5*bxby[4]*grad_u_y[8]+0.31943828249997*bxbx[8]*grad_u_x[8]+0.5*bxbx[4]*grad_u_x[8]+0.5*grad_u_z[4]*bxbz[8]+0.5*grad_u_y[4]*bxby[8]+0.5*grad_u_x[4]*bxbx[8]+0.31943828249997*bxbz[7]*grad_u_z[7]+0.5000000000000001*bxbz[1]*grad_u_z[7]+0.31943828249997*bxby[7]*grad_u_y[7]+0.5000000000000001*bxby[1]*grad_u_y[7]+0.31943828249997*bxbx[7]*grad_u_x[7]+0.5000000000000001*bxbx[1]*grad_u_x[7]+0.5000000000000001*grad_u_z[1]*bxbz[7]+0.5000000000000001*grad_u_y[1]*bxby[7]+0.5000000000000001*grad_u_x[1]*bxbx[7]+0.4472135954999579*bxbz[6]*grad_u_z[6]+0.4472135954999579*bxby[6]*grad_u_y[6]+0.4472135954999579*bxbx[6]*grad_u_x[6]+0.31943828249997*bxbz[5]*grad_u_z[5]+0.5*bxbz[0]*grad_u_z[5]+0.31943828249997*bxby[5]*grad_u_y[5]+0.5*bxby[0]*grad_u_y[5]+0.31943828249997*bxbx[5]*grad_u_x[5]+0.5*bxbx[0]*grad_u_x[5]+0.5*grad_u_z[0]*bxbz[5]+0.5*grad_u_y[0]*bxby[5]+0.5*grad_u_x[0]*bxbx[5]+0.4472135954999579*bxbz[3]*grad_u_z[3]+0.4472135954999579*bxby[3]*grad_u_y[3]+0.4472135954999579*bxbx[3]*grad_u_x[3]+0.4472135954999579*bxbz[2]*grad_u_z[2]+0.4472135954999579*bxby[2]*grad_u_y[2]+0.4472135954999579*bxbx[2]*grad_u_x[2]; 
  bb_grad_u_comp[6] = 0.2857142857142857*bxbz[6]*grad_u_z[8]+0.447213595499958*bxbz[2]*grad_u_z[8]+0.2857142857142857*bxby[6]*grad_u_y[8]+0.447213595499958*bxby[2]*grad_u_y[8]+0.2857142857142857*bxbx[6]*grad_u_x[8]+0.447213595499958*bxbx[2]*grad_u_x[8]+0.2857142857142857*grad_u_z[6]*bxbz[8]+0.447213595499958*grad_u_z[2]*bxbz[8]+0.2857142857142857*grad_u_y[6]*bxby[8]+0.447213595499958*grad_u_y[2]*bxby[8]+0.2857142857142857*grad_u_x[6]*bxbx[8]+0.447213595499958*grad_u_x[2]*bxbx[8]+0.4*bxbz[3]*grad_u_z[7]+0.4*bxby[3]*grad_u_y[7]+0.4*bxbx[3]*grad_u_x[7]+0.4*grad_u_z[3]*bxbz[7]+0.4*grad_u_y[3]*bxby[7]+0.4*grad_u_x[3]*bxbx[7]+0.4472135954999579*bxbz[5]*grad_u_z[6]+0.31943828249997*bxbz[4]*grad_u_z[6]+0.5*bxbz[0]*grad_u_z[6]+0.4472135954999579*bxby[5]*grad_u_y[6]+0.31943828249997*bxby[4]*grad_u_y[6]+0.5*bxby[0]*grad_u_y[6]+0.4472135954999579*bxbx[5]*grad_u_x[6]+0.31943828249997*bxbx[4]*grad_u_x[6]+0.5*bxbx[0]*grad_u_x[6]+0.4472135954999579*grad_u_z[5]*bxbz[6]+0.31943828249997*grad_u_z[4]*bxbz[6]+0.5*grad_u_z[0]*bxbz[6]+0.4472135954999579*grad_u_y[5]*bxby[6]+0.31943828249997*grad_u_y[4]*bxby[6]+0.5*grad_u_y[0]*bxby[6]+0.4472135954999579*grad_u_x[5]*bxbx[6]+0.31943828249997*grad_u_x[4]*bxbx[6]+0.5*grad_u_x[0]*bxbx[6]+0.5000000000000001*bxbz[2]*grad_u_z[4]+0.5000000000000001*bxby[2]*grad_u_y[4]+0.5000000000000001*bxbx[2]*grad_u_x[4]+0.5000000000000001*grad_u_z[2]*bxbz[4]+0.5000000000000001*grad_u_y[2]*bxby[4]+0.5000000000000001*grad_u_x[2]*bxbx[4]+0.447213595499958*bxbz[1]*grad_u_z[3]+0.447213595499958*bxby[1]*grad_u_y[3]+0.447213595499958*bxbx[1]*grad_u_x[3]+0.447213595499958*grad_u_z[1]*bxbz[3]+0.447213595499958*grad_u_y[1]*bxby[3]+0.447213595499958*grad_u_x[1]*bxbx[3]; 
  bb_grad_u_comp[7] = 0.2857142857142857*bxbz[7]*grad_u_z[8]+0.447213595499958*bxbz[1]*grad_u_z[8]+0.2857142857142857*bxby[7]*grad_u_y[8]+0.447213595499958*bxby[1]*grad_u_y[8]+0.2857142857142857*bxbx[7]*grad_u_x[8]+0.447213595499958*bxbx[1]*grad_u_x[8]+0.2857142857142857*grad_u_z[7]*bxbz[8]+0.447213595499958*grad_u_z[1]*bxbz[8]+0.2857142857142857*grad_u_y[7]*bxby[8]+0.447213595499958*grad_u_y[1]*bxby[8]+0.2857142857142857*grad_u_x[7]*bxbx[8]+0.447213595499958*grad_u_x[1]*bxbx[8]+0.31943828249997*bxbz[5]*grad_u_z[7]+0.4472135954999579*bxbz[4]*grad_u_z[7]+0.5*bxbz[0]*grad_u_z[7]+0.31943828249997*bxby[5]*grad_u_y[7]+0.4472135954999579*bxby[4]*grad_u_y[7]+0.5*bxby[0]*grad_u_y[7]+0.31943828249997*bxbx[5]*grad_u_x[7]+0.4472135954999579*bxbx[4]*grad_u_x[7]+0.5*bxbx[0]*grad_u_x[7]+0.31943828249997*grad_u_z[5]*bxbz[7]+0.4472135954999579*grad_u_z[4]*bxbz[7]+0.5*grad_u_z[0]*bxbz[7]+0.31943828249997*grad_u_y[5]*bxby[7]+0.4472135954999579*grad_u_y[4]*bxby[7]+0.5*grad_u_y[0]*bxby[7]+0.31943828249997*grad_u_x[5]*bxbx[7]+0.4472135954999579*grad_u_x[4]*bxbx[7]+0.5*grad_u_x[0]*bxbx[7]+0.4*bxbz[3]*grad_u_z[6]+0.4*bxby[3]*grad_u_y[6]+0.4*bxbx[3]*grad_u_x[6]+0.4*grad_u_z[3]*bxbz[6]+0.4*grad_u_y[3]*bxby[6]+0.4*grad_u_x[3]*bxbx[6]+0.5000000000000001*bxbz[1]*grad_u_z[5]+0.5000000000000001*bxby[1]*grad_u_y[5]+0.5000000000000001*bxbx[1]*grad_u_x[5]+0.5000000000000001*grad_u_z[1]*bxbz[5]+0.5000000000000001*grad_u_y[1]*bxby[5]+0.5000000000000001*grad_u_x[1]*bxbx[5]+0.447213595499958*bxbz[2]*grad_u_z[3]+0.447213595499958*bxby[2]*grad_u_y[3]+0.447213595499958*bxbx[2]*grad_u_x[3]+0.447213595499958*grad_u_z[2]*bxbz[3]+0.447213595499958*grad_u_y[2]*bxby[3]+0.447213595499958*grad_u_x[2]*bxbx[3]; 
  bb_grad_u_comp[8] = 0.2040816326530612*bxbz[8]*grad_u_z[8]+0.31943828249997*bxbz[5]*grad_u_z[8]+0.31943828249997*bxbz[4]*grad_u_z[8]+0.5*bxbz[0]*grad_u_z[8]+0.2040816326530612*bxby[8]*grad_u_y[8]+0.31943828249997*bxby[5]*grad_u_y[8]+0.31943828249997*bxby[4]*grad_u_y[8]+0.5*bxby[0]*grad_u_y[8]+0.2040816326530612*bxbx[8]*grad_u_x[8]+0.31943828249997*bxbx[5]*grad_u_x[8]+0.31943828249997*bxbx[4]*grad_u_x[8]+0.5*bxbx[0]*grad_u_x[8]+0.31943828249997*grad_u_z[5]*bxbz[8]+0.31943828249997*grad_u_z[4]*bxbz[8]+0.5*grad_u_z[0]*bxbz[8]+0.31943828249997*grad_u_y[5]*bxby[8]+0.31943828249997*grad_u_y[4]*bxby[8]+0.5*grad_u_y[0]*bxby[8]+0.31943828249997*grad_u_x[5]*bxbx[8]+0.31943828249997*grad_u_x[4]*bxbx[8]+0.5*grad_u_x[0]*bxbx[8]+0.2857142857142857*bxbz[7]*grad_u_z[7]+0.447213595499958*bxbz[1]*grad_u_z[7]+0.2857142857142857*bxby[7]*grad_u_y[7]+0.447213595499958*bxby[1]*grad_u_y[7]+0.2857142857142857*bxbx[7]*grad_u_x[7]+0.447213595499958*bxbx[1]*grad_u_x[7]+0.447213595499958*grad_u_z[1]*bxbz[7]+0.447213595499958*grad_u_y[1]*bxby[7]+0.447213595499958*grad_u_x[1]*bxbx[7]+0.2857142857142857*bxbz[6]*grad_u_z[6]+0.447213595499958*bxbz[2]*grad_u_z[6]+0.2857142857142857*bxby[6]*grad_u_y[6]+0.447213595499958*bxby[2]*grad_u_y[6]+0.2857142857142857*bxbx[6]*grad_u_x[6]+0.447213595499958*bxbx[2]*grad_u_x[6]+0.447213595499958*grad_u_z[2]*bxbz[6]+0.447213595499958*grad_u_y[2]*bxby[6]+0.447213595499958*grad_u_x[2]*bxbx[6]+0.5*bxbz[4]*grad_u_z[5]+0.5*bxby[4]*grad_u_y[5]+0.5*bxbx[4]*grad_u_x[5]+0.5*grad_u_z[4]*bxbz[5]+0.5*grad_u_y[4]*bxby[5]+0.5*grad_u_x[4]*bxbx[5]+0.4*bxbz[3]*grad_u_z[3]+0.4*bxby[3]*grad_u_y[3]+0.4*bxbx[3]*grad_u_x[3]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 
  bb_grad_u[8] += bb_grad_u_comp[8]; 

  p_perp_source[0] += bb_grad_u_comp[0]-1.0*(nu_c[0]+grad_u_x[0]); 
  p_perp_source[1] += bb_grad_u_comp[1]-1.0*(nu_c[1]+grad_u_x[1]); 
  p_perp_source[2] += bb_grad_u_comp[2]-1.0*(nu_c[2]+grad_u_x[2]); 
  p_perp_source[3] += bb_grad_u_comp[3]-1.0*(nu_c[3]+grad_u_x[3]); 
  p_perp_source[4] += bb_grad_u_comp[4]-1.0*(nu_c[4]+grad_u_x[4]); 
  p_perp_source[5] += bb_grad_u_comp[5]-1.0*(nu_c[5]+grad_u_x[5]); 
  p_perp_source[6] += bb_grad_u_comp[6]-1.0*(nu_c[6]+grad_u_x[6]); 
  p_perp_source[7] += bb_grad_u_comp[7]-1.0*(nu_c[7]+grad_u_x[7]); 
  p_perp_source[8] += bb_grad_u_comp[8]-1.0*(nu_c[8]+grad_u_x[8]); 

} 
