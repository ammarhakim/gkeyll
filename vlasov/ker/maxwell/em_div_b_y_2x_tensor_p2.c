#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void em_div_b_y_2x_tensor_p2(const double *dxv, 
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector in left/center/right cells in each direction. 
  // bvar_c:          Input volume expansion of magnetic field unit vector in center cell. 
  // max_b:           Output surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|. 
  // div_b:           Output volume expansion of div(b).

  const double dx1 = 2.0/dxv[1]; 
  const double *b_c = &bvar_c[9]; 
  const double *b_surf_lr = &bvar_surf_l[27]; 
  const double *b_surf_cl = &bvar_surf_c[24]; 
  const double *b_surf_cr = &bvar_surf_c[27]; 
  const double *b_surf_rl = &bvar_surf_r[24]; 

  double *max_b_l = &max_b[6]; 
  double *max_b_r = &max_b[9]; 

  double bl_r = 0.0; 
  double bc_l = 0.0; 
  double bc_r = 0.0; 
  double br_l = 0.0; 
  double max_b_quad_l[3] = {0.0}; 
  double max_b_quad_r[3] = {0.0}; 

  bl_r = 0.6324555320336759*b_surf_lr[2]-0.9486832980505137*b_surf_lr[1]+0.7071067811865475*b_surf_lr[0]; 
  bc_l = 0.6324555320336759*b_surf_cl[2]-0.9486832980505137*b_surf_cl[1]+0.7071067811865475*b_surf_cl[0]; 
  bc_r = 0.6324555320336759*b_surf_cr[2]-0.9486832980505137*b_surf_cr[1]+0.7071067811865475*b_surf_cr[0]; 
  br_l = 0.6324555320336759*b_surf_rl[2]-0.9486832980505137*b_surf_rl[1]+0.7071067811865475*b_surf_rl[0]; 
  max_b_quad_l[0] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[0] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.7071067811865475*b_surf_lr[0]-0.7905694150420947*b_surf_lr[2]; 
  bc_l = 0.7071067811865475*b_surf_cl[0]-0.7905694150420947*b_surf_cl[2]; 
  bc_r = 0.7071067811865475*b_surf_cr[0]-0.7905694150420947*b_surf_cr[2]; 
  br_l = 0.7071067811865475*b_surf_rl[0]-0.7905694150420947*b_surf_rl[2]; 
  max_b_quad_l[1] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[1] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.6324555320336759*b_surf_lr[2]+0.9486832980505137*b_surf_lr[1]+0.7071067811865475*b_surf_lr[0]; 
  bc_l = 0.6324555320336759*b_surf_cl[2]+0.9486832980505137*b_surf_cl[1]+0.7071067811865475*b_surf_cl[0]; 
  bc_r = 0.6324555320336759*b_surf_cr[2]+0.9486832980505137*b_surf_cr[1]+0.7071067811865475*b_surf_cr[0]; 
  br_l = 0.6324555320336759*b_surf_rl[2]+0.9486832980505137*b_surf_rl[1]+0.7071067811865475*b_surf_rl[0]; 
  max_b_quad_l[2] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[2] = fmax(fabs(bc_r), fabs(br_l)); 

  tensor_2x_p2_upwind_quad_to_modal(max_b_quad_l, max_b_l); 
  tensor_2x_p2_upwind_quad_to_modal(max_b_quad_r, max_b_r); 

  div_b[0] += (0.3535533905932737*b_surf_rl[0]-0.3535533905932737*b_surf_lr[0]+0.3535533905932737*b_surf_cr[0]-0.3535533905932737*b_surf_cl[0])*dx1; 
  div_b[1] += (0.3535533905932737*b_surf_rl[1]-0.3535533905932737*b_surf_lr[1]+0.3535533905932737*b_surf_cr[1]-0.3535533905932737*b_surf_cl[1])*dx1; 
  div_b[2] += (0.6123724356957944*(b_surf_rl[0]+b_surf_lr[0]+b_surf_cr[0]+b_surf_cl[0])-1.732050807568877*b_c[0])*dx1; 
  div_b[3] += (0.6123724356957944*(b_surf_rl[1]+b_surf_lr[1]+b_surf_cr[1]+b_surf_cl[1])-1.732050807568877*b_c[1])*dx1; 
  div_b[4] += (0.3535533905932737*b_surf_rl[2]-0.3535533905932737*b_surf_lr[2]+0.3535533905932737*b_surf_cr[2]-0.3535533905932737*b_surf_cl[2])*dx1; 
  div_b[5] += ((-3.872983346207417*b_c[2])+0.7905694150420947*b_surf_rl[0]-0.7905694150420947*b_surf_lr[0]+0.7905694150420947*b_surf_cr[0]-0.7905694150420947*b_surf_cl[0])*dx1; 
  div_b[6] += (0.6123724356957944*(b_surf_rl[2]+b_surf_lr[2]+b_surf_cr[2]+b_surf_cl[2])-1.732050807568877*b_c[4])*dx1; 
  div_b[7] += ((-3.872983346207417*b_c[3])+0.7905694150420948*b_surf_rl[1]-0.7905694150420948*b_surf_lr[1]+0.7905694150420948*b_surf_cr[1]-0.7905694150420948*b_surf_cl[1])*dx1; 
  div_b[8] += ((-3.872983346207417*b_c[6])+0.7905694150420947*b_surf_rl[2]-0.7905694150420947*b_surf_lr[2]+0.7905694150420947*b_surf_cr[2]-0.7905694150420947*b_surf_cl[2])*dx1; 

} 
