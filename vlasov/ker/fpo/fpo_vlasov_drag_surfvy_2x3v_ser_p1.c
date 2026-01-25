#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 


GKYL_CU_DH double fpo_vlasov_drag_surfvy_2x3v_ser_p1(const double* dxv, 
  const double *alpha_surf_L, const double *alpha_surf_R,
  const double *sgn_alpha_surf_L, const double *sgn_alpha_surf_R,
  const int *const_sgn_alpha_L, const int *const_sgn_alpha_R,
  const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: Cell spacing in each direction. 
  // alpha_surf_L,R: Surface expansion of drag coefficient on left,right boundary of center cell. 
  // sgn_alpha_L,R: sign(alpha_surf_l,r) at quadrature points. 
  // const_sgn_alpha_L,R: Boolean array true if sign(alpha_surf_l,r) is only one sign. 
  // fL, fC, fR: Distribution function in left, center, and right cells. 
  // out: Incremented output. 


  // Index into drag coefficient surface expansion arrays 
  const double *drag_coeff_surf_L = &alpha_surf_L[16]; 
  const double *drag_coeff_surf_R = &alpha_surf_R[16]; 
  const double *sgn_drag_coeff_surf_L = &sgn_alpha_surf_L[16]; 
  const double *sgn_drag_coeff_surf_R = &sgn_alpha_surf_R[16]; 
  const int *const_sgn_drag_coeff_L = &const_sgn_alpha_L[1]; 
  const int *const_sgn_drag_coeff_R = &const_sgn_alpha_R[1]; 
  double dv_inv = 2.0/dxv[3]; 


  double fUp_L[16] = {0.0}; 
  if (const_sgn_drag_coeff_L[0] == 1) { 
    if (sgn_drag_coeff_surf_L[0] == 1.0) { 
  fUp_L[0] = 1.5811388300841895*fL[48]+1.224744871391589*fL[4]+0.7071067811865475*fL[0]; 
  fUp_L[1] = 1.5811388300841898*fL[49]+1.224744871391589*fL[9]+0.7071067811865475*fL[1]; 
  fUp_L[2] = 1.5811388300841898*fL[50]+1.224744871391589*fL[10]+0.7071067811865475*fL[2]; 
  fUp_L[3] = 1.5811388300841898*fL[51]+1.224744871391589*fL[11]+0.7071067811865475*fL[3]; 
  fUp_L[4] = 1.5811388300841898*fL[52]+1.224744871391589*fL[15]+0.7071067811865475*fL[5]; 
  fUp_L[5] = 1.5811388300841895*fL[53]+1.224744871391589*fL[17]+0.7071067811865475*fL[6]; 
  fUp_L[6] = 1.5811388300841895*fL[54]+1.224744871391589*fL[18]+0.7071067811865475*fL[7]; 
  fUp_L[7] = 1.5811388300841895*fL[55]+1.224744871391589*fL[19]+0.7071067811865475*fL[8]; 
  fUp_L[8] = 1.5811388300841895*fL[56]+1.224744871391589*fL[23]+0.7071067811865475*fL[12]; 
  fUp_L[9] = 1.5811388300841895*fL[57]+1.224744871391589*fL[24]+0.7071067811865475*fL[13]; 
  fUp_L[10] = 1.5811388300841895*fL[58]+1.224744871391589*fL[25]+0.7071067811865475*fL[14]; 
  fUp_L[11] = 1.5811388300841898*fL[59]+1.224744871391589*fL[26]+0.7071067811865475*fL[16]; 
  fUp_L[12] = 1.5811388300841898*fL[60]+1.224744871391589*fL[28]+0.7071067811865475*fL[20]; 
  fUp_L[13] = 1.5811388300841898*fL[61]+1.224744871391589*fL[29]+0.7071067811865475*fL[21]; 
  fUp_L[14] = 1.5811388300841898*fL[62]+1.224744871391589*fL[30]+0.7071067811865475*fL[22]; 
  fUp_L[15] = 1.5811388300841895*fL[63]+1.224744871391589*fL[31]+0.7071067811865475*fL[27]; 
    } else { 
  fUp_L[0] = 1.5811388300841895*fC[48]-1.224744871391589*fC[4]+0.7071067811865475*fC[0]; 
  fUp_L[1] = 1.5811388300841898*fC[49]-1.224744871391589*fC[9]+0.7071067811865475*fC[1]; 
  fUp_L[2] = 1.5811388300841898*fC[50]-1.224744871391589*fC[10]+0.7071067811865475*fC[2]; 
  fUp_L[3] = 1.5811388300841898*fC[51]-1.224744871391589*fC[11]+0.7071067811865475*fC[3]; 
  fUp_L[4] = 1.5811388300841898*fC[52]-1.224744871391589*fC[15]+0.7071067811865475*fC[5]; 
  fUp_L[5] = 1.5811388300841895*fC[53]-1.224744871391589*fC[17]+0.7071067811865475*fC[6]; 
  fUp_L[6] = 1.5811388300841895*fC[54]-1.224744871391589*fC[18]+0.7071067811865475*fC[7]; 
  fUp_L[7] = 1.5811388300841895*fC[55]-1.224744871391589*fC[19]+0.7071067811865475*fC[8]; 
  fUp_L[8] = 1.5811388300841895*fC[56]-1.224744871391589*fC[23]+0.7071067811865475*fC[12]; 
  fUp_L[9] = 1.5811388300841895*fC[57]-1.224744871391589*fC[24]+0.7071067811865475*fC[13]; 
  fUp_L[10] = 1.5811388300841895*fC[58]-1.224744871391589*fC[25]+0.7071067811865475*fC[14]; 
  fUp_L[11] = 1.5811388300841898*fC[59]-1.224744871391589*fC[26]+0.7071067811865475*fC[16]; 
  fUp_L[12] = 1.5811388300841898*fC[60]-1.224744871391589*fC[28]+0.7071067811865475*fC[20]; 
  fUp_L[13] = 1.5811388300841898*fC[61]-1.224744871391589*fC[29]+0.7071067811865475*fC[21]; 
  fUp_L[14] = 1.5811388300841898*fC[62]-1.224744871391589*fC[30]+0.7071067811865475*fC[22]; 
  fUp_L[15] = 1.5811388300841895*fC[63]-1.224744871391589*fC[31]+0.7071067811865475*fC[27]; 
   } 
  } else { 
  double fL_r[16] = {0.0}; 
  double fC_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_L[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_L, sgn_drag_coeff_Up_L); 

  fL_r[0] = 1.5811388300841895*fL[48]+1.224744871391589*fL[4]+0.7071067811865475*fL[0]; 
  fL_r[1] = 1.5811388300841898*fL[49]+1.224744871391589*fL[9]+0.7071067811865475*fL[1]; 
  fL_r[2] = 1.5811388300841898*fL[50]+1.224744871391589*fL[10]+0.7071067811865475*fL[2]; 
  fL_r[3] = 1.5811388300841898*fL[51]+1.224744871391589*fL[11]+0.7071067811865475*fL[3]; 
  fL_r[4] = 1.5811388300841898*fL[52]+1.224744871391589*fL[15]+0.7071067811865475*fL[5]; 
  fL_r[5] = 1.5811388300841895*fL[53]+1.224744871391589*fL[17]+0.7071067811865475*fL[6]; 
  fL_r[6] = 1.5811388300841895*fL[54]+1.224744871391589*fL[18]+0.7071067811865475*fL[7]; 
  fL_r[7] = 1.5811388300841895*fL[55]+1.224744871391589*fL[19]+0.7071067811865475*fL[8]; 
  fL_r[8] = 1.5811388300841895*fL[56]+1.224744871391589*fL[23]+0.7071067811865475*fL[12]; 
  fL_r[9] = 1.5811388300841895*fL[57]+1.224744871391589*fL[24]+0.7071067811865475*fL[13]; 
  fL_r[10] = 1.5811388300841895*fL[58]+1.224744871391589*fL[25]+0.7071067811865475*fL[14]; 
  fL_r[11] = 1.5811388300841898*fL[59]+1.224744871391589*fL[26]+0.7071067811865475*fL[16]; 
  fL_r[12] = 1.5811388300841898*fL[60]+1.224744871391589*fL[28]+0.7071067811865475*fL[20]; 
  fL_r[13] = 1.5811388300841898*fL[61]+1.224744871391589*fL[29]+0.7071067811865475*fL[21]; 
  fL_r[14] = 1.5811388300841898*fL[62]+1.224744871391589*fL[30]+0.7071067811865475*fL[22]; 
  fL_r[15] = 1.5811388300841895*fL[63]+1.224744871391589*fL[31]+0.7071067811865475*fL[27]; 

  fC_l[0] = 1.5811388300841895*fC[48]-1.224744871391589*fC[4]+0.7071067811865475*fC[0]; 
  fC_l[1] = 1.5811388300841898*fC[49]-1.224744871391589*fC[9]+0.7071067811865475*fC[1]; 
  fC_l[2] = 1.5811388300841898*fC[50]-1.224744871391589*fC[10]+0.7071067811865475*fC[2]; 
  fC_l[3] = 1.5811388300841898*fC[51]-1.224744871391589*fC[11]+0.7071067811865475*fC[3]; 
  fC_l[4] = 1.5811388300841898*fC[52]-1.224744871391589*fC[15]+0.7071067811865475*fC[5]; 
  fC_l[5] = 1.5811388300841895*fC[53]-1.224744871391589*fC[17]+0.7071067811865475*fC[6]; 
  fC_l[6] = 1.5811388300841895*fC[54]-1.224744871391589*fC[18]+0.7071067811865475*fC[7]; 
  fC_l[7] = 1.5811388300841895*fC[55]-1.224744871391589*fC[19]+0.7071067811865475*fC[8]; 
  fC_l[8] = 1.5811388300841895*fC[56]-1.224744871391589*fC[23]+0.7071067811865475*fC[12]; 
  fC_l[9] = 1.5811388300841895*fC[57]-1.224744871391589*fC[24]+0.7071067811865475*fC[13]; 
  fC_l[10] = 1.5811388300841895*fC[58]-1.224744871391589*fC[25]+0.7071067811865475*fC[14]; 
  fC_l[11] = 1.5811388300841898*fC[59]-1.224744871391589*fC[26]+0.7071067811865475*fC[16]; 
  fC_l[12] = 1.5811388300841898*fC[60]-1.224744871391589*fC[28]+0.7071067811865475*fC[20]; 
  fC_l[13] = 1.5811388300841898*fC[61]-1.224744871391589*fC[29]+0.7071067811865475*fC[21]; 
  fC_l[14] = 1.5811388300841898*fC[62]-1.224744871391589*fC[30]+0.7071067811865475*fC[22]; 
  fC_l[15] = 1.5811388300841895*fC[63]-1.224744871391589*fC[31]+0.7071067811865475*fC[27]; 

  fUp_L[0] = 0.125*(fL_r[15]*sgn_drag_coeff_Up_L[15]+fL_r[14]*sgn_drag_coeff_Up_L[14]+fL_r[13]*sgn_drag_coeff_Up_L[13]+fL_r[12]*sgn_drag_coeff_Up_L[12]+fL_r[11]*sgn_drag_coeff_Up_L[11]+fL_r[10]*sgn_drag_coeff_Up_L[10]+fL_r[9]*sgn_drag_coeff_Up_L[9]+fL_r[8]*sgn_drag_coeff_Up_L[8]+fL_r[7]*sgn_drag_coeff_Up_L[7]+fL_r[6]*sgn_drag_coeff_Up_L[6]+fL_r[5]*sgn_drag_coeff_Up_L[5]+fL_r[4]*sgn_drag_coeff_Up_L[4]+fL_r[3]*sgn_drag_coeff_Up_L[3]+fL_r[2]*sgn_drag_coeff_Up_L[2]+fL_r[1]*sgn_drag_coeff_Up_L[1]+fL_r[0]*sgn_drag_coeff_Up_L[0])-0.125*(fC_l[15]*sgn_drag_coeff_Up_L[15]+fC_l[14]*sgn_drag_coeff_Up_L[14]+fC_l[13]*sgn_drag_coeff_Up_L[13]+fC_l[12]*sgn_drag_coeff_Up_L[12]+fC_l[11]*sgn_drag_coeff_Up_L[11]+fC_l[10]*sgn_drag_coeff_Up_L[10]+fC_l[9]*sgn_drag_coeff_Up_L[9]+fC_l[8]*sgn_drag_coeff_Up_L[8]+fC_l[7]*sgn_drag_coeff_Up_L[7]+fC_l[6]*sgn_drag_coeff_Up_L[6]+fC_l[5]*sgn_drag_coeff_Up_L[5]+fC_l[4]*sgn_drag_coeff_Up_L[4]+fC_l[3]*sgn_drag_coeff_Up_L[3]+fC_l[2]*sgn_drag_coeff_Up_L[2]+fC_l[1]*sgn_drag_coeff_Up_L[1]+fC_l[0]*sgn_drag_coeff_Up_L[0])+0.5*(fL_r[0]+fC_l[0]); 
  fUp_L[1] = 0.125*(fL_r[14]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[14]*fL_r[15]+fL_r[10]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[10]*fL_r[13]+fL_r[9]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[9]*fL_r[12]+fL_r[7]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[7]*fL_r[11]+fL_r[4]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[4]*fL_r[8]+fL_r[3]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[3]*fL_r[6]+fL_r[2]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[2]*fL_r[5]+fL_r[0]*sgn_drag_coeff_Up_L[1]+sgn_drag_coeff_Up_L[0]*fL_r[1])-0.125*(fC_l[14]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[14]*fC_l[15]+fC_l[10]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[10]*fC_l[13]+fC_l[9]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[9]*fC_l[12]+fC_l[7]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[7]*fC_l[11]+fC_l[4]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[4]*fC_l[8]+fC_l[3]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[3]*fC_l[6]+fC_l[2]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[2]*fC_l[5]+fC_l[0]*sgn_drag_coeff_Up_L[1]+sgn_drag_coeff_Up_L[0]*fC_l[1])+0.5*(fL_r[1]+fC_l[1]); 
  fUp_L[2] = 0.125*(fL_r[13]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[13]*fL_r[15]+fL_r[10]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[10]*fL_r[14]+fL_r[8]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[8]*fL_r[12]+fL_r[6]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[6]*fL_r[11]+fL_r[4]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[4]*fL_r[9]+fL_r[3]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[3]*fL_r[7]+fL_r[1]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[1]*fL_r[5]+fL_r[0]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[0]*fL_r[2])-0.125*(fC_l[13]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[13]*fC_l[15]+fC_l[10]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[10]*fC_l[14]+fC_l[8]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[8]*fC_l[12]+fC_l[6]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[6]*fC_l[11]+fC_l[4]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[4]*fC_l[9]+fC_l[3]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[3]*fC_l[7]+fC_l[1]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[1]*fC_l[5]+fC_l[0]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[0]*fC_l[2])+0.5*(fL_r[2]+fC_l[2]); 
  fUp_L[3] = 0.125*(fL_r[12]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[12]*fL_r[15]+fL_r[9]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[9]*fL_r[14]+fL_r[8]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[8]*fL_r[13]+fL_r[5]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[5]*fL_r[11]+fL_r[4]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[4]*fL_r[10]+fL_r[2]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[2]*fL_r[7]+fL_r[1]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[1]*fL_r[6]+fL_r[0]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[0]*fL_r[3])-0.125*(fC_l[12]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[12]*fC_l[15]+fC_l[9]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[9]*fC_l[14]+fC_l[8]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[8]*fC_l[13]+fC_l[5]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[5]*fC_l[11]+fC_l[4]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[4]*fC_l[10]+fC_l[2]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[2]*fC_l[7]+fC_l[1]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[1]*fC_l[6]+fC_l[0]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[0]*fC_l[3])+0.5*(fL_r[3]+fC_l[3]); 
  fUp_L[4] = 0.125*(fL_r[11]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[11]*fL_r[15]+fL_r[7]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[7]*fL_r[14]+fL_r[6]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[6]*fL_r[13]+fL_r[5]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[5]*fL_r[12]+fL_r[3]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[3]*fL_r[10]+fL_r[2]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[2]*fL_r[9]+fL_r[1]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[1]*fL_r[8]+fL_r[0]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[0]*fL_r[4])-0.125*(fC_l[11]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[11]*fC_l[15]+fC_l[7]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[7]*fC_l[14]+fC_l[6]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[6]*fC_l[13]+fC_l[5]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[5]*fC_l[12]+fC_l[3]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[3]*fC_l[10]+fC_l[2]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[2]*fC_l[9]+fC_l[1]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[1]*fC_l[8]+fC_l[0]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[0]*fC_l[4])+0.5*(fL_r[4]+fC_l[4]); 
  fUp_L[5] = 0.125*(fL_r[10]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[10]*fL_r[15]+fL_r[13]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[13]*fL_r[14]+fL_r[4]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[4]*fL_r[12]+fL_r[3]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[3]*fL_r[11]+fL_r[8]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[8]*fL_r[9]+fL_r[6]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[6]*fL_r[7]+fL_r[0]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[0]*fL_r[5]+fL_r[1]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[1]*fL_r[2])-0.125*(fC_l[10]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[10]*fC_l[15]+fC_l[13]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[13]*fC_l[14]+fC_l[4]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[4]*fC_l[12]+fC_l[3]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[3]*fC_l[11]+fC_l[8]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[8]*fC_l[9]+fC_l[6]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[6]*fC_l[7]+fC_l[0]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[0]*fC_l[5]+fC_l[1]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[1]*fC_l[2])+0.5*(fL_r[5]+fC_l[5]); 
  fUp_L[6] = 0.125*(fL_r[9]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[9]*fL_r[15]+fL_r[12]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[12]*fL_r[14]+fL_r[4]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[4]*fL_r[13]+fL_r[2]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[2]*fL_r[11]+fL_r[8]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[8]*fL_r[10]+fL_r[5]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[5]*fL_r[7]+fL_r[0]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[0]*fL_r[6]+fL_r[1]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[1]*fL_r[3])-0.125*(fC_l[9]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[9]*fC_l[15]+fC_l[12]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[12]*fC_l[14]+fC_l[4]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[4]*fC_l[13]+fC_l[2]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[2]*fC_l[11]+fC_l[8]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[8]*fC_l[10]+fC_l[5]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[5]*fC_l[7]+fC_l[0]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[0]*fC_l[6]+fC_l[1]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[1]*fC_l[3])+0.5*(fL_r[6]+fC_l[6]); 
  fUp_L[7] = 0.125*(fL_r[8]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[8]*fL_r[15]+fL_r[4]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[4]*fL_r[14]+fL_r[12]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[12]*fL_r[13]+fL_r[1]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[1]*fL_r[11]+fL_r[9]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[9]*fL_r[10]+fL_r[0]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[0]*fL_r[7]+fL_r[5]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[5]*fL_r[6]+fL_r[2]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[2]*fL_r[3])-0.125*(fC_l[8]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[8]*fC_l[15]+fC_l[4]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[4]*fC_l[14]+fC_l[12]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[12]*fC_l[13]+fC_l[1]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[1]*fC_l[11]+fC_l[9]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[9]*fC_l[10]+fC_l[0]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[0]*fC_l[7]+fC_l[5]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[5]*fC_l[6]+fC_l[2]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[2]*fC_l[3])+0.5*(fL_r[7]+fC_l[7]); 
  fUp_L[8] = 0.125*(fL_r[7]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[7]*fL_r[15]+fL_r[11]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[11]*fL_r[14]+fL_r[3]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[3]*fL_r[13]+fL_r[2]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[2]*fL_r[12]+fL_r[6]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[6]*fL_r[10]+fL_r[5]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[5]*fL_r[9]+fL_r[0]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[0]*fL_r[8]+fL_r[1]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[1]*fL_r[4])-0.125*(fC_l[7]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[7]*fC_l[15]+fC_l[11]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[11]*fC_l[14]+fC_l[3]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[3]*fC_l[13]+fC_l[2]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[2]*fC_l[12]+fC_l[6]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[6]*fC_l[10]+fC_l[5]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[5]*fC_l[9]+fC_l[0]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[0]*fC_l[8]+fC_l[1]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[1]*fC_l[4])+0.5*(fL_r[8]+fC_l[8]); 
  fUp_L[9] = 0.125*(fL_r[6]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[6]*fL_r[15]+fL_r[3]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[3]*fL_r[14]+fL_r[11]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[11]*fL_r[13]+fL_r[1]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[1]*fL_r[12]+fL_r[7]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[7]*fL_r[10]+fL_r[0]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[0]*fL_r[9]+fL_r[5]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[5]*fL_r[8]+fL_r[2]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[2]*fL_r[4])-0.125*(fC_l[6]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[6]*fC_l[15]+fC_l[3]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[3]*fC_l[14]+fC_l[11]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[11]*fC_l[13]+fC_l[1]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[1]*fC_l[12]+fC_l[7]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[7]*fC_l[10]+fC_l[0]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[0]*fC_l[9]+fC_l[5]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[5]*fC_l[8]+fC_l[2]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[2]*fC_l[4])+0.5*(fL_r[9]+fC_l[9]); 
  fUp_L[10] = 0.125*(fL_r[5]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[5]*fL_r[15]+fL_r[2]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[2]*fL_r[14]+fL_r[1]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[1]*fL_r[13]+fL_r[11]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[11]*fL_r[12]+fL_r[0]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[0]*fL_r[10]+fL_r[7]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[7]*fL_r[9]+fL_r[6]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[6]*fL_r[8]+fL_r[3]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[3]*fL_r[4])-0.125*(fC_l[5]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[5]*fC_l[15]+fC_l[2]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[2]*fC_l[14]+fC_l[1]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[1]*fC_l[13]+fC_l[11]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[11]*fC_l[12]+fC_l[0]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[0]*fC_l[10]+fC_l[7]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[7]*fC_l[9]+fC_l[6]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[6]*fC_l[8]+fC_l[3]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[3]*fC_l[4])+0.5*(fL_r[10]+fC_l[10]); 
  fUp_L[11] = 0.125*(fL_r[4]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[4]*fL_r[15]+fL_r[8]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[8]*fL_r[14]+fL_r[9]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[9]*fL_r[13]+fL_r[10]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[10]*fL_r[12]+fL_r[0]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[0]*fL_r[11]+fL_r[1]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[1]*fL_r[7]+fL_r[2]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[2]*fL_r[6]+fL_r[3]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[3]*fL_r[5])-0.125*(fC_l[4]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[4]*fC_l[15]+fC_l[8]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[8]*fC_l[14]+fC_l[9]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[9]*fC_l[13]+fC_l[10]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[10]*fC_l[12]+fC_l[0]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[0]*fC_l[11]+fC_l[1]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[1]*fC_l[7]+fC_l[2]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[2]*fC_l[6]+fC_l[3]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[3]*fC_l[5])+0.5*(fL_r[11]+fC_l[11]); 
  fUp_L[12] = 0.125*(fL_r[3]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[3]*fL_r[15]+fL_r[6]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[6]*fL_r[14]+fL_r[7]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[7]*fL_r[13]+fL_r[0]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[0]*fL_r[12]+fL_r[10]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[10]*fL_r[11]+fL_r[1]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[1]*fL_r[9]+fL_r[2]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[2]*fL_r[8]+fL_r[4]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[4]*fL_r[5])-0.125*(fC_l[3]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[3]*fC_l[15]+fC_l[6]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[6]*fC_l[14]+fC_l[7]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[7]*fC_l[13]+fC_l[0]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[0]*fC_l[12]+fC_l[10]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[10]*fC_l[11]+fC_l[1]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[1]*fC_l[9]+fC_l[2]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[2]*fC_l[8]+fC_l[4]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[4]*fC_l[5])+0.5*(fL_r[12]+fC_l[12]); 
  fUp_L[13] = 0.125*(fL_r[2]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[2]*fL_r[15]+fL_r[5]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[5]*fL_r[14]+fL_r[0]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[0]*fL_r[13]+fL_r[7]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[7]*fL_r[12]+fL_r[9]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[9]*fL_r[11]+fL_r[1]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[1]*fL_r[10]+fL_r[3]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[3]*fL_r[8]+fL_r[4]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[4]*fL_r[6])-0.125*(fC_l[2]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[2]*fC_l[15]+fC_l[5]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[5]*fC_l[14]+fC_l[0]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[0]*fC_l[13]+fC_l[7]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[7]*fC_l[12]+fC_l[9]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[9]*fC_l[11]+fC_l[1]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[1]*fC_l[10]+fC_l[3]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[3]*fC_l[8]+fC_l[4]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[4]*fC_l[6])+0.5*(fL_r[13]+fC_l[13]); 
  fUp_L[14] = 0.125*(fL_r[1]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[1]*fL_r[15]+fL_r[0]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[0]*fL_r[14]+fL_r[5]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[5]*fL_r[13]+fL_r[6]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[6]*fL_r[12]+fL_r[8]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[8]*fL_r[11]+fL_r[2]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[2]*fL_r[10]+fL_r[3]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[3]*fL_r[9]+fL_r[4]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[4]*fL_r[7])-0.125*(fC_l[1]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[1]*fC_l[15]+fC_l[0]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[0]*fC_l[14]+fC_l[5]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[5]*fC_l[13]+fC_l[6]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[6]*fC_l[12]+fC_l[8]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[8]*fC_l[11]+fC_l[2]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[2]*fC_l[10]+fC_l[3]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[3]*fC_l[9]+fC_l[4]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[4]*fC_l[7])+0.5*(fL_r[14]+fC_l[14]); 
  fUp_L[15] = 0.125*(fL_r[0]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[0]*fL_r[15]+fL_r[1]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[1]*fL_r[14]+fL_r[2]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[2]*fL_r[13]+fL_r[3]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[3]*fL_r[12]+fL_r[4]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[4]*fL_r[11]+fL_r[5]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[5]*fL_r[10]+fL_r[6]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[6]*fL_r[9]+fL_r[7]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[7]*fL_r[8])-0.125*(fC_l[0]*sgn_drag_coeff_Up_L[15]+sgn_drag_coeff_Up_L[0]*fC_l[15]+fC_l[1]*sgn_drag_coeff_Up_L[14]+sgn_drag_coeff_Up_L[1]*fC_l[14]+fC_l[2]*sgn_drag_coeff_Up_L[13]+sgn_drag_coeff_Up_L[2]*fC_l[13]+fC_l[3]*sgn_drag_coeff_Up_L[12]+sgn_drag_coeff_Up_L[3]*fC_l[12]+fC_l[4]*sgn_drag_coeff_Up_L[11]+sgn_drag_coeff_Up_L[4]*fC_l[11]+fC_l[5]*sgn_drag_coeff_Up_L[10]+sgn_drag_coeff_Up_L[5]*fC_l[10]+fC_l[6]*sgn_drag_coeff_Up_L[9]+sgn_drag_coeff_Up_L[6]*fC_l[9]+fC_l[7]*sgn_drag_coeff_Up_L[8]+sgn_drag_coeff_Up_L[7]*fC_l[8])+0.5*(fL_r[15]+fC_l[15]); 

  } 
  double fUp_R[16] = {0.0}; 
  if (const_sgn_drag_coeff_R[0] == 1) { 
    if (sgn_drag_coeff_surf_R[0] == 1.0) { 
  fUp_R[0] = 1.5811388300841895*fC[48]+1.224744871391589*fC[4]+0.7071067811865475*fC[0]; 
  fUp_R[1] = 1.5811388300841898*fC[49]+1.224744871391589*fC[9]+0.7071067811865475*fC[1]; 
  fUp_R[2] = 1.5811388300841898*fC[50]+1.224744871391589*fC[10]+0.7071067811865475*fC[2]; 
  fUp_R[3] = 1.5811388300841898*fC[51]+1.224744871391589*fC[11]+0.7071067811865475*fC[3]; 
  fUp_R[4] = 1.5811388300841898*fC[52]+1.224744871391589*fC[15]+0.7071067811865475*fC[5]; 
  fUp_R[5] = 1.5811388300841895*fC[53]+1.224744871391589*fC[17]+0.7071067811865475*fC[6]; 
  fUp_R[6] = 1.5811388300841895*fC[54]+1.224744871391589*fC[18]+0.7071067811865475*fC[7]; 
  fUp_R[7] = 1.5811388300841895*fC[55]+1.224744871391589*fC[19]+0.7071067811865475*fC[8]; 
  fUp_R[8] = 1.5811388300841895*fC[56]+1.224744871391589*fC[23]+0.7071067811865475*fC[12]; 
  fUp_R[9] = 1.5811388300841895*fC[57]+1.224744871391589*fC[24]+0.7071067811865475*fC[13]; 
  fUp_R[10] = 1.5811388300841895*fC[58]+1.224744871391589*fC[25]+0.7071067811865475*fC[14]; 
  fUp_R[11] = 1.5811388300841898*fC[59]+1.224744871391589*fC[26]+0.7071067811865475*fC[16]; 
  fUp_R[12] = 1.5811388300841898*fC[60]+1.224744871391589*fC[28]+0.7071067811865475*fC[20]; 
  fUp_R[13] = 1.5811388300841898*fC[61]+1.224744871391589*fC[29]+0.7071067811865475*fC[21]; 
  fUp_R[14] = 1.5811388300841898*fC[62]+1.224744871391589*fC[30]+0.7071067811865475*fC[22]; 
  fUp_R[15] = 1.5811388300841895*fC[63]+1.224744871391589*fC[31]+0.7071067811865475*fC[27]; 
    } else { 
  fUp_R[0] = 1.5811388300841895*fR[48]-1.224744871391589*fR[4]+0.7071067811865475*fR[0]; 
  fUp_R[1] = 1.5811388300841898*fR[49]-1.224744871391589*fR[9]+0.7071067811865475*fR[1]; 
  fUp_R[2] = 1.5811388300841898*fR[50]-1.224744871391589*fR[10]+0.7071067811865475*fR[2]; 
  fUp_R[3] = 1.5811388300841898*fR[51]-1.224744871391589*fR[11]+0.7071067811865475*fR[3]; 
  fUp_R[4] = 1.5811388300841898*fR[52]-1.224744871391589*fR[15]+0.7071067811865475*fR[5]; 
  fUp_R[5] = 1.5811388300841895*fR[53]-1.224744871391589*fR[17]+0.7071067811865475*fR[6]; 
  fUp_R[6] = 1.5811388300841895*fR[54]-1.224744871391589*fR[18]+0.7071067811865475*fR[7]; 
  fUp_R[7] = 1.5811388300841895*fR[55]-1.224744871391589*fR[19]+0.7071067811865475*fR[8]; 
  fUp_R[8] = 1.5811388300841895*fR[56]-1.224744871391589*fR[23]+0.7071067811865475*fR[12]; 
  fUp_R[9] = 1.5811388300841895*fR[57]-1.224744871391589*fR[24]+0.7071067811865475*fR[13]; 
  fUp_R[10] = 1.5811388300841895*fR[58]-1.224744871391589*fR[25]+0.7071067811865475*fR[14]; 
  fUp_R[11] = 1.5811388300841898*fR[59]-1.224744871391589*fR[26]+0.7071067811865475*fR[16]; 
  fUp_R[12] = 1.5811388300841898*fR[60]-1.224744871391589*fR[28]+0.7071067811865475*fR[20]; 
  fUp_R[13] = 1.5811388300841898*fR[61]-1.224744871391589*fR[29]+0.7071067811865475*fR[21]; 
  fUp_R[14] = 1.5811388300841898*fR[62]-1.224744871391589*fR[30]+0.7071067811865475*fR[22]; 
  fUp_R[15] = 1.5811388300841895*fR[63]-1.224744871391589*fR[31]+0.7071067811865475*fR[27]; 
   } 
  } else { 
  double fC_r[16] = {0.0}; 
  double fR_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_R[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_R, sgn_drag_coeff_Up_R); 

  fC_r[0] = 1.5811388300841895*fC[48]+1.224744871391589*fC[4]+0.7071067811865475*fC[0]; 
  fC_r[1] = 1.5811388300841898*fC[49]+1.224744871391589*fC[9]+0.7071067811865475*fC[1]; 
  fC_r[2] = 1.5811388300841898*fC[50]+1.224744871391589*fC[10]+0.7071067811865475*fC[2]; 
  fC_r[3] = 1.5811388300841898*fC[51]+1.224744871391589*fC[11]+0.7071067811865475*fC[3]; 
  fC_r[4] = 1.5811388300841898*fC[52]+1.224744871391589*fC[15]+0.7071067811865475*fC[5]; 
  fC_r[5] = 1.5811388300841895*fC[53]+1.224744871391589*fC[17]+0.7071067811865475*fC[6]; 
  fC_r[6] = 1.5811388300841895*fC[54]+1.224744871391589*fC[18]+0.7071067811865475*fC[7]; 
  fC_r[7] = 1.5811388300841895*fC[55]+1.224744871391589*fC[19]+0.7071067811865475*fC[8]; 
  fC_r[8] = 1.5811388300841895*fC[56]+1.224744871391589*fC[23]+0.7071067811865475*fC[12]; 
  fC_r[9] = 1.5811388300841895*fC[57]+1.224744871391589*fC[24]+0.7071067811865475*fC[13]; 
  fC_r[10] = 1.5811388300841895*fC[58]+1.224744871391589*fC[25]+0.7071067811865475*fC[14]; 
  fC_r[11] = 1.5811388300841898*fC[59]+1.224744871391589*fC[26]+0.7071067811865475*fC[16]; 
  fC_r[12] = 1.5811388300841898*fC[60]+1.224744871391589*fC[28]+0.7071067811865475*fC[20]; 
  fC_r[13] = 1.5811388300841898*fC[61]+1.224744871391589*fC[29]+0.7071067811865475*fC[21]; 
  fC_r[14] = 1.5811388300841898*fC[62]+1.224744871391589*fC[30]+0.7071067811865475*fC[22]; 
  fC_r[15] = 1.5811388300841895*fC[63]+1.224744871391589*fC[31]+0.7071067811865475*fC[27]; 

  fR_l[0] = 1.5811388300841895*fR[48]-1.224744871391589*fR[4]+0.7071067811865475*fR[0]; 
  fR_l[1] = 1.5811388300841898*fR[49]-1.224744871391589*fR[9]+0.7071067811865475*fR[1]; 
  fR_l[2] = 1.5811388300841898*fR[50]-1.224744871391589*fR[10]+0.7071067811865475*fR[2]; 
  fR_l[3] = 1.5811388300841898*fR[51]-1.224744871391589*fR[11]+0.7071067811865475*fR[3]; 
  fR_l[4] = 1.5811388300841898*fR[52]-1.224744871391589*fR[15]+0.7071067811865475*fR[5]; 
  fR_l[5] = 1.5811388300841895*fR[53]-1.224744871391589*fR[17]+0.7071067811865475*fR[6]; 
  fR_l[6] = 1.5811388300841895*fR[54]-1.224744871391589*fR[18]+0.7071067811865475*fR[7]; 
  fR_l[7] = 1.5811388300841895*fR[55]-1.224744871391589*fR[19]+0.7071067811865475*fR[8]; 
  fR_l[8] = 1.5811388300841895*fR[56]-1.224744871391589*fR[23]+0.7071067811865475*fR[12]; 
  fR_l[9] = 1.5811388300841895*fR[57]-1.224744871391589*fR[24]+0.7071067811865475*fR[13]; 
  fR_l[10] = 1.5811388300841895*fR[58]-1.224744871391589*fR[25]+0.7071067811865475*fR[14]; 
  fR_l[11] = 1.5811388300841898*fR[59]-1.224744871391589*fR[26]+0.7071067811865475*fR[16]; 
  fR_l[12] = 1.5811388300841898*fR[60]-1.224744871391589*fR[28]+0.7071067811865475*fR[20]; 
  fR_l[13] = 1.5811388300841898*fR[61]-1.224744871391589*fR[29]+0.7071067811865475*fR[21]; 
  fR_l[14] = 1.5811388300841898*fR[62]-1.224744871391589*fR[30]+0.7071067811865475*fR[22]; 
  fR_l[15] = 1.5811388300841895*fR[63]-1.224744871391589*fR[31]+0.7071067811865475*fR[27]; 

  fUp_R[0] = -(0.125*(fR_l[15]*sgn_drag_coeff_Up_R[15]+fR_l[14]*sgn_drag_coeff_Up_R[14]+fR_l[13]*sgn_drag_coeff_Up_R[13]+fR_l[12]*sgn_drag_coeff_Up_R[12]+fR_l[11]*sgn_drag_coeff_Up_R[11]+fR_l[10]*sgn_drag_coeff_Up_R[10]+fR_l[9]*sgn_drag_coeff_Up_R[9]+fR_l[8]*sgn_drag_coeff_Up_R[8]+fR_l[7]*sgn_drag_coeff_Up_R[7]+fR_l[6]*sgn_drag_coeff_Up_R[6]+fR_l[5]*sgn_drag_coeff_Up_R[5]+fR_l[4]*sgn_drag_coeff_Up_R[4]+fR_l[3]*sgn_drag_coeff_Up_R[3]+fR_l[2]*sgn_drag_coeff_Up_R[2]+fR_l[1]*sgn_drag_coeff_Up_R[1]+fR_l[0]*sgn_drag_coeff_Up_R[0]))+0.125*(fC_r[15]*sgn_drag_coeff_Up_R[15]+fC_r[14]*sgn_drag_coeff_Up_R[14]+fC_r[13]*sgn_drag_coeff_Up_R[13]+fC_r[12]*sgn_drag_coeff_Up_R[12]+fC_r[11]*sgn_drag_coeff_Up_R[11]+fC_r[10]*sgn_drag_coeff_Up_R[10]+fC_r[9]*sgn_drag_coeff_Up_R[9]+fC_r[8]*sgn_drag_coeff_Up_R[8]+fC_r[7]*sgn_drag_coeff_Up_R[7]+fC_r[6]*sgn_drag_coeff_Up_R[6]+fC_r[5]*sgn_drag_coeff_Up_R[5]+fC_r[4]*sgn_drag_coeff_Up_R[4]+fC_r[3]*sgn_drag_coeff_Up_R[3]+fC_r[2]*sgn_drag_coeff_Up_R[2]+fC_r[1]*sgn_drag_coeff_Up_R[1]+fC_r[0]*sgn_drag_coeff_Up_R[0])+0.5*(fR_l[0]+fC_r[0]); 
  fUp_R[1] = -(0.125*(fR_l[14]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[14]*fR_l[15]+fR_l[10]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[10]*fR_l[13]+fR_l[9]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[9]*fR_l[12]+fR_l[7]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[7]*fR_l[11]+fR_l[4]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[4]*fR_l[8]+fR_l[3]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[3]*fR_l[6]+fR_l[2]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[2]*fR_l[5]+fR_l[0]*sgn_drag_coeff_Up_R[1]+sgn_drag_coeff_Up_R[0]*fR_l[1]))+0.125*(fC_r[14]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[14]*fC_r[15]+fC_r[10]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[10]*fC_r[13]+fC_r[9]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[9]*fC_r[12]+fC_r[7]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[7]*fC_r[11]+fC_r[4]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[4]*fC_r[8]+fC_r[3]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[3]*fC_r[6]+fC_r[2]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[2]*fC_r[5]+fC_r[0]*sgn_drag_coeff_Up_R[1]+sgn_drag_coeff_Up_R[0]*fC_r[1])+0.5*(fR_l[1]+fC_r[1]); 
  fUp_R[2] = -(0.125*(fR_l[13]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[13]*fR_l[15]+fR_l[10]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[10]*fR_l[14]+fR_l[8]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[8]*fR_l[12]+fR_l[6]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[6]*fR_l[11]+fR_l[4]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[4]*fR_l[9]+fR_l[3]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[3]*fR_l[7]+fR_l[1]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[1]*fR_l[5]+fR_l[0]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[0]*fR_l[2]))+0.125*(fC_r[13]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[13]*fC_r[15]+fC_r[10]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[10]*fC_r[14]+fC_r[8]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[8]*fC_r[12]+fC_r[6]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[6]*fC_r[11]+fC_r[4]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[4]*fC_r[9]+fC_r[3]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[3]*fC_r[7]+fC_r[1]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[1]*fC_r[5]+fC_r[0]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[0]*fC_r[2])+0.5*(fR_l[2]+fC_r[2]); 
  fUp_R[3] = -(0.125*(fR_l[12]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[12]*fR_l[15]+fR_l[9]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[9]*fR_l[14]+fR_l[8]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[8]*fR_l[13]+fR_l[5]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[5]*fR_l[11]+fR_l[4]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[4]*fR_l[10]+fR_l[2]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[2]*fR_l[7]+fR_l[1]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[1]*fR_l[6]+fR_l[0]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[0]*fR_l[3]))+0.125*(fC_r[12]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[12]*fC_r[15]+fC_r[9]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[9]*fC_r[14]+fC_r[8]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[8]*fC_r[13]+fC_r[5]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[5]*fC_r[11]+fC_r[4]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[4]*fC_r[10]+fC_r[2]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[2]*fC_r[7]+fC_r[1]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[1]*fC_r[6]+fC_r[0]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[0]*fC_r[3])+0.5*(fR_l[3]+fC_r[3]); 
  fUp_R[4] = -(0.125*(fR_l[11]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[11]*fR_l[15]+fR_l[7]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[7]*fR_l[14]+fR_l[6]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[6]*fR_l[13]+fR_l[5]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[5]*fR_l[12]+fR_l[3]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[3]*fR_l[10]+fR_l[2]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[2]*fR_l[9]+fR_l[1]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[1]*fR_l[8]+fR_l[0]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[0]*fR_l[4]))+0.125*(fC_r[11]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[11]*fC_r[15]+fC_r[7]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[7]*fC_r[14]+fC_r[6]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[6]*fC_r[13]+fC_r[5]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[5]*fC_r[12]+fC_r[3]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[3]*fC_r[10]+fC_r[2]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[2]*fC_r[9]+fC_r[1]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[1]*fC_r[8]+fC_r[0]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[0]*fC_r[4])+0.5*(fR_l[4]+fC_r[4]); 
  fUp_R[5] = -(0.125*(fR_l[10]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[10]*fR_l[15]+fR_l[13]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[13]*fR_l[14]+fR_l[4]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[4]*fR_l[12]+fR_l[3]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[3]*fR_l[11]+fR_l[8]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[8]*fR_l[9]+fR_l[6]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[6]*fR_l[7]+fR_l[0]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[0]*fR_l[5]+fR_l[1]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[1]*fR_l[2]))+0.125*(fC_r[10]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[10]*fC_r[15]+fC_r[13]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[13]*fC_r[14]+fC_r[4]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[4]*fC_r[12]+fC_r[3]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[3]*fC_r[11]+fC_r[8]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[8]*fC_r[9]+fC_r[6]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[6]*fC_r[7]+fC_r[0]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[0]*fC_r[5]+fC_r[1]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[1]*fC_r[2])+0.5*(fR_l[5]+fC_r[5]); 
  fUp_R[6] = -(0.125*(fR_l[9]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[9]*fR_l[15]+fR_l[12]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[12]*fR_l[14]+fR_l[4]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[4]*fR_l[13]+fR_l[2]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[2]*fR_l[11]+fR_l[8]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[8]*fR_l[10]+fR_l[5]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[5]*fR_l[7]+fR_l[0]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[0]*fR_l[6]+fR_l[1]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[1]*fR_l[3]))+0.125*(fC_r[9]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[9]*fC_r[15]+fC_r[12]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[12]*fC_r[14]+fC_r[4]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[4]*fC_r[13]+fC_r[2]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[2]*fC_r[11]+fC_r[8]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[8]*fC_r[10]+fC_r[5]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[5]*fC_r[7]+fC_r[0]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[0]*fC_r[6]+fC_r[1]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[1]*fC_r[3])+0.5*(fR_l[6]+fC_r[6]); 
  fUp_R[7] = -(0.125*(fR_l[8]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[8]*fR_l[15]+fR_l[4]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[4]*fR_l[14]+fR_l[12]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[12]*fR_l[13]+fR_l[1]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[1]*fR_l[11]+fR_l[9]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[9]*fR_l[10]+fR_l[0]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[0]*fR_l[7]+fR_l[5]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[5]*fR_l[6]+fR_l[2]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[2]*fR_l[3]))+0.125*(fC_r[8]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[8]*fC_r[15]+fC_r[4]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[4]*fC_r[14]+fC_r[12]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[12]*fC_r[13]+fC_r[1]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[1]*fC_r[11]+fC_r[9]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[9]*fC_r[10]+fC_r[0]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[0]*fC_r[7]+fC_r[5]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[5]*fC_r[6]+fC_r[2]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[2]*fC_r[3])+0.5*(fR_l[7]+fC_r[7]); 
  fUp_R[8] = -(0.125*(fR_l[7]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[7]*fR_l[15]+fR_l[11]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[11]*fR_l[14]+fR_l[3]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[3]*fR_l[13]+fR_l[2]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[2]*fR_l[12]+fR_l[6]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[6]*fR_l[10]+fR_l[5]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[5]*fR_l[9]+fR_l[0]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[0]*fR_l[8]+fR_l[1]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[1]*fR_l[4]))+0.125*(fC_r[7]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[7]*fC_r[15]+fC_r[11]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[11]*fC_r[14]+fC_r[3]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[3]*fC_r[13]+fC_r[2]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[2]*fC_r[12]+fC_r[6]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[6]*fC_r[10]+fC_r[5]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[5]*fC_r[9]+fC_r[0]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[0]*fC_r[8]+fC_r[1]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[1]*fC_r[4])+0.5*(fR_l[8]+fC_r[8]); 
  fUp_R[9] = -(0.125*(fR_l[6]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[6]*fR_l[15]+fR_l[3]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[3]*fR_l[14]+fR_l[11]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[11]*fR_l[13]+fR_l[1]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[1]*fR_l[12]+fR_l[7]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[7]*fR_l[10]+fR_l[0]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[0]*fR_l[9]+fR_l[5]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[5]*fR_l[8]+fR_l[2]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[2]*fR_l[4]))+0.125*(fC_r[6]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[6]*fC_r[15]+fC_r[3]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[3]*fC_r[14]+fC_r[11]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[11]*fC_r[13]+fC_r[1]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[1]*fC_r[12]+fC_r[7]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[7]*fC_r[10]+fC_r[0]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[0]*fC_r[9]+fC_r[5]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[5]*fC_r[8]+fC_r[2]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[2]*fC_r[4])+0.5*(fR_l[9]+fC_r[9]); 
  fUp_R[10] = -(0.125*(fR_l[5]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[5]*fR_l[15]+fR_l[2]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[2]*fR_l[14]+fR_l[1]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[1]*fR_l[13]+fR_l[11]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[11]*fR_l[12]+fR_l[0]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[0]*fR_l[10]+fR_l[7]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[7]*fR_l[9]+fR_l[6]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[6]*fR_l[8]+fR_l[3]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[3]*fR_l[4]))+0.125*(fC_r[5]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[5]*fC_r[15]+fC_r[2]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[2]*fC_r[14]+fC_r[1]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[1]*fC_r[13]+fC_r[11]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[11]*fC_r[12]+fC_r[0]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[0]*fC_r[10]+fC_r[7]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[7]*fC_r[9]+fC_r[6]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[6]*fC_r[8]+fC_r[3]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[3]*fC_r[4])+0.5*(fR_l[10]+fC_r[10]); 
  fUp_R[11] = -(0.125*(fR_l[4]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[4]*fR_l[15]+fR_l[8]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[8]*fR_l[14]+fR_l[9]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[9]*fR_l[13]+fR_l[10]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[10]*fR_l[12]+fR_l[0]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[0]*fR_l[11]+fR_l[1]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[1]*fR_l[7]+fR_l[2]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[2]*fR_l[6]+fR_l[3]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[3]*fR_l[5]))+0.125*(fC_r[4]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[4]*fC_r[15]+fC_r[8]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[8]*fC_r[14]+fC_r[9]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[9]*fC_r[13]+fC_r[10]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[10]*fC_r[12]+fC_r[0]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[0]*fC_r[11]+fC_r[1]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[1]*fC_r[7]+fC_r[2]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[2]*fC_r[6]+fC_r[3]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[3]*fC_r[5])+0.5*(fR_l[11]+fC_r[11]); 
  fUp_R[12] = -(0.125*(fR_l[3]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[3]*fR_l[15]+fR_l[6]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[6]*fR_l[14]+fR_l[7]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[7]*fR_l[13]+fR_l[0]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[0]*fR_l[12]+fR_l[10]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[10]*fR_l[11]+fR_l[1]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[1]*fR_l[9]+fR_l[2]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[2]*fR_l[8]+fR_l[4]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[4]*fR_l[5]))+0.125*(fC_r[3]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[3]*fC_r[15]+fC_r[6]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[6]*fC_r[14]+fC_r[7]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[7]*fC_r[13]+fC_r[0]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[0]*fC_r[12]+fC_r[10]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[10]*fC_r[11]+fC_r[1]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[1]*fC_r[9]+fC_r[2]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[2]*fC_r[8]+fC_r[4]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[4]*fC_r[5])+0.5*(fR_l[12]+fC_r[12]); 
  fUp_R[13] = -(0.125*(fR_l[2]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[2]*fR_l[15]+fR_l[5]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[5]*fR_l[14]+fR_l[0]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[0]*fR_l[13]+fR_l[7]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[7]*fR_l[12]+fR_l[9]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[9]*fR_l[11]+fR_l[1]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[1]*fR_l[10]+fR_l[3]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[3]*fR_l[8]+fR_l[4]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[4]*fR_l[6]))+0.125*(fC_r[2]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[2]*fC_r[15]+fC_r[5]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[5]*fC_r[14]+fC_r[0]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[0]*fC_r[13]+fC_r[7]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[7]*fC_r[12]+fC_r[9]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[9]*fC_r[11]+fC_r[1]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[1]*fC_r[10]+fC_r[3]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[3]*fC_r[8]+fC_r[4]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[4]*fC_r[6])+0.5*(fR_l[13]+fC_r[13]); 
  fUp_R[14] = -(0.125*(fR_l[1]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[1]*fR_l[15]+fR_l[0]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[0]*fR_l[14]+fR_l[5]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[5]*fR_l[13]+fR_l[6]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[6]*fR_l[12]+fR_l[8]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[8]*fR_l[11]+fR_l[2]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[2]*fR_l[10]+fR_l[3]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[3]*fR_l[9]+fR_l[4]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[4]*fR_l[7]))+0.125*(fC_r[1]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[1]*fC_r[15]+fC_r[0]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[0]*fC_r[14]+fC_r[5]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[5]*fC_r[13]+fC_r[6]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[6]*fC_r[12]+fC_r[8]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[8]*fC_r[11]+fC_r[2]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[2]*fC_r[10]+fC_r[3]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[3]*fC_r[9]+fC_r[4]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[4]*fC_r[7])+0.5*(fR_l[14]+fC_r[14]); 
  fUp_R[15] = -(0.125*(fR_l[0]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[0]*fR_l[15]+fR_l[1]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[1]*fR_l[14]+fR_l[2]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[2]*fR_l[13]+fR_l[3]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[3]*fR_l[12]+fR_l[4]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[4]*fR_l[11]+fR_l[5]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[5]*fR_l[10]+fR_l[6]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[6]*fR_l[9]+fR_l[7]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[7]*fR_l[8]))+0.125*(fC_r[0]*sgn_drag_coeff_Up_R[15]+sgn_drag_coeff_Up_R[0]*fC_r[15]+fC_r[1]*sgn_drag_coeff_Up_R[14]+sgn_drag_coeff_Up_R[1]*fC_r[14]+fC_r[2]*sgn_drag_coeff_Up_R[13]+sgn_drag_coeff_Up_R[2]*fC_r[13]+fC_r[3]*sgn_drag_coeff_Up_R[12]+sgn_drag_coeff_Up_R[3]*fC_r[12]+fC_r[4]*sgn_drag_coeff_Up_R[11]+sgn_drag_coeff_Up_R[4]*fC_r[11]+fC_r[5]*sgn_drag_coeff_Up_R[10]+sgn_drag_coeff_Up_R[5]*fC_r[10]+fC_r[6]*sgn_drag_coeff_Up_R[9]+sgn_drag_coeff_Up_R[6]*fC_r[9]+fC_r[7]*sgn_drag_coeff_Up_R[8]+sgn_drag_coeff_Up_R[7]*fC_r[8])+0.5*(fR_l[15]+fC_r[15]); 

  } 
  double GhatL[16] = {0.0}; 
  double GhatR[16] = {0.0}; 


  out[0] += (0.35355339059327373*GhatL[0]-0.35355339059327373*GhatR[0])*dv_inv; 
  out[1] += (0.35355339059327373*GhatL[1]-0.35355339059327373*GhatR[1])*dv_inv; 
  out[2] += (0.35355339059327373*GhatL[2]-0.35355339059327373*GhatR[2])*dv_inv; 
  out[3] += (0.35355339059327373*GhatL[3]-0.35355339059327373*GhatR[3])*dv_inv; 
  out[4] += -(0.6123724356957945*(GhatR[0]+GhatL[0])*dv_inv); 
  out[5] += (0.35355339059327373*GhatL[4]-0.35355339059327373*GhatR[4])*dv_inv; 
  out[6] += (0.35355339059327373*GhatL[5]-0.35355339059327373*GhatR[5])*dv_inv; 
  out[7] += (0.35355339059327373*GhatL[6]-0.35355339059327373*GhatR[6])*dv_inv; 
  out[8] += (0.35355339059327373*GhatL[7]-0.35355339059327373*GhatR[7])*dv_inv; 
  out[9] += -(0.6123724356957945*(GhatR[1]+GhatL[1])*dv_inv); 
  out[10] += -(0.6123724356957945*(GhatR[2]+GhatL[2])*dv_inv); 
  out[11] += -(0.6123724356957945*(GhatR[3]+GhatL[3])*dv_inv); 
  out[12] += (0.35355339059327373*GhatL[8]-0.35355339059327373*GhatR[8])*dv_inv; 
  out[13] += (0.35355339059327373*GhatL[9]-0.35355339059327373*GhatR[9])*dv_inv; 
  out[14] += (0.35355339059327373*GhatL[10]-0.35355339059327373*GhatR[10])*dv_inv; 
  out[15] += -(0.6123724356957945*(GhatR[4]+GhatL[4])*dv_inv); 
  out[16] += (0.35355339059327373*GhatL[11]-0.35355339059327373*GhatR[11])*dv_inv; 
  out[17] += -(0.6123724356957945*(GhatR[5]+GhatL[5])*dv_inv); 
  out[18] += -(0.6123724356957945*(GhatR[6]+GhatL[6])*dv_inv); 
  out[19] += -(0.6123724356957945*(GhatR[7]+GhatL[7])*dv_inv); 
  out[20] += (0.35355339059327373*GhatL[12]-0.35355339059327373*GhatR[12])*dv_inv; 
  out[21] += (0.35355339059327373*GhatL[13]-0.35355339059327373*GhatR[13])*dv_inv; 
  out[22] += (0.35355339059327373*GhatL[14]-0.35355339059327373*GhatR[14])*dv_inv; 
  out[23] += -(0.6123724356957945*(GhatR[8]+GhatL[8])*dv_inv); 
  out[24] += -(0.6123724356957945*(GhatR[9]+GhatL[9])*dv_inv); 
  out[25] += -(0.6123724356957945*(GhatR[10]+GhatL[10])*dv_inv); 
  out[26] += -(0.6123724356957945*(GhatR[11]+GhatL[11])*dv_inv); 
  out[27] += (0.35355339059327373*GhatL[15]-0.35355339059327373*GhatR[15])*dv_inv; 
  out[28] += -(0.6123724356957945*(GhatR[12]+GhatL[12])*dv_inv); 
  out[29] += -(0.6123724356957945*(GhatR[13]+GhatL[13])*dv_inv); 
  out[30] += -(0.6123724356957945*(GhatR[14]+GhatL[14])*dv_inv); 
  out[31] += -(0.6123724356957945*(GhatR[15]+GhatL[15])*dv_inv); 
  out[48] += (0.7905694150420948*GhatL[0]-0.7905694150420948*GhatR[0])*dv_inv; 
  out[49] += (0.7905694150420949*GhatL[1]-0.7905694150420949*GhatR[1])*dv_inv; 
  out[50] += (0.7905694150420949*GhatL[2]-0.7905694150420949*GhatR[2])*dv_inv; 
  out[51] += (0.7905694150420949*GhatL[3]-0.7905694150420949*GhatR[3])*dv_inv; 
  out[52] += (0.7905694150420949*GhatL[4]-0.7905694150420949*GhatR[4])*dv_inv; 
  out[53] += (0.7905694150420948*GhatL[5]-0.7905694150420948*GhatR[5])*dv_inv; 
  out[54] += (0.7905694150420948*GhatL[6]-0.7905694150420948*GhatR[6])*dv_inv; 
  out[55] += (0.7905694150420948*GhatL[7]-0.7905694150420948*GhatR[7])*dv_inv; 
  out[56] += (0.7905694150420948*GhatL[8]-0.7905694150420948*GhatR[8])*dv_inv; 
  out[57] += (0.7905694150420948*GhatL[9]-0.7905694150420948*GhatR[9])*dv_inv; 
  out[58] += (0.7905694150420948*GhatL[10]-0.7905694150420948*GhatR[10])*dv_inv; 
  out[59] += (0.7905694150420949*GhatL[11]-0.7905694150420949*GhatR[11])*dv_inv; 
  out[60] += (0.7905694150420949*GhatL[12]-0.7905694150420949*GhatR[12])*dv_inv; 
  out[61] += (0.7905694150420949*GhatL[13]-0.7905694150420949*GhatR[13])*dv_inv; 
  out[62] += (0.7905694150420949*GhatL[14]-0.7905694150420949*GhatR[14])*dv_inv; 
  out[63] += (0.7905694150420948*GhatL[15]-0.7905694150420948*GhatR[15])*dv_inv; 
  double cflFreq = fmax(fabs(drag_coeff_surf_L[0]), fabs(drag_coeff_surf_R[0])); 

  return 0.375*dv_inv*cflFreq; 
}
