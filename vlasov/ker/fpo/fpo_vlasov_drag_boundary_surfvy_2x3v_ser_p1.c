#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 


GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvy_2x3v_ser_p1(const double* dxv,
  const double *alpha_surf_Edge, const double *alpha_surf_Skin,
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin,
  const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin,
  const int edge, const double *fEdge, const double *fSkin,
  double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // alpha_surf_Edge,Skin: Surface expansion of drag coefficient on lower boundary of Skin/Edge cell. 
  // sgn_alpha_Edge,Skin: sign(alpha_surf_l,r) at quadrature points. 
  // const_sgn_alpha_Edge,Skin: Boolean array true if sign(alpha_surf_l,r) is only one sign. 
  // fEdge, fSkin: Distribution function in left, center, and right cells. 
  // out: Incremented output. 


  // Index into drag coefficient surface expansion arrays 
  const double *drag_coeff_surf_Edge = &alpha_surf_Edge[16]; 
  const double *drag_coeff_surf_Skin = &alpha_surf_Skin[16]; 
  const double *sgn_drag_coeff_surf_Edge = &sgn_alpha_surf_Edge[16]; 
  const double *sgn_drag_coeff_surf_Skin = &sgn_alpha_surf_Skin[16]; 
  const int *const_sgn_drag_coeff_Edge = &const_sgn_alpha_Edge[1]; 
  const int *const_sgn_drag_coeff_Skin = &const_sgn_alpha_Skin[1]; 
  double dv_inv = 2.0/dxv[3]; 

  double cflFreq = 0.0; 

  if (edge == -1) { 
  const double *drag_coeff_surf_R = drag_coeff_surf_Edge; 
  const double *sgn_drag_coeff_surf_R = sgn_drag_coeff_surf_Edge; 
  const int *const_sgn_drag_coeff_R = const_sgn_drag_coeff_Edge; 
  double fUp_R[16] = {0.0}; 
  if (const_sgn_drag_coeff_R[0] == 1) { 
    if (sgn_drag_coeff_surf_R[0] == 1.0) { 
  fUp_R[0] = 1.5811388300841895*fSkin[48]+1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[0]; 
  fUp_R[1] = 1.5811388300841898*fSkin[49]+1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[1]; 
  fUp_R[2] = 1.5811388300841898*fSkin[50]+1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[2]; 
  fUp_R[3] = 1.5811388300841898*fSkin[51]+1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[3]; 
  fUp_R[4] = 1.5811388300841898*fSkin[52]+1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[5]; 
  fUp_R[5] = 1.5811388300841895*fSkin[53]+1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[6]; 
  fUp_R[6] = 1.5811388300841895*fSkin[54]+1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[7]; 
  fUp_R[7] = 1.5811388300841895*fSkin[55]+1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[8]; 
  fUp_R[8] = 1.5811388300841895*fSkin[56]+1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[12]; 
  fUp_R[9] = 1.5811388300841895*fSkin[57]+1.224744871391589*fSkin[24]+0.7071067811865475*fSkin[13]; 
  fUp_R[10] = 1.5811388300841895*fSkin[58]+1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[14]; 
  fUp_R[11] = 1.5811388300841898*fSkin[59]+1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[16]; 
  fUp_R[12] = 1.5811388300841898*fSkin[60]+1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[20]; 
  fUp_R[13] = 1.5811388300841898*fSkin[61]+1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[21]; 
  fUp_R[14] = 1.5811388300841898*fSkin[62]+1.224744871391589*fSkin[30]+0.7071067811865475*fSkin[22]; 
  fUp_R[15] = 1.5811388300841895*fSkin[63]+1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[27]; 
    } else { 
  fUp_R[0] = 1.5811388300841895*fEdge[48]-1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[0]; 
  fUp_R[1] = 1.5811388300841898*fEdge[49]-1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[1]; 
  fUp_R[2] = 1.5811388300841898*fEdge[50]-1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[2]; 
  fUp_R[3] = 1.5811388300841898*fEdge[51]-1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[3]; 
  fUp_R[4] = 1.5811388300841898*fEdge[52]-1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[5]; 
  fUp_R[5] = 1.5811388300841895*fEdge[53]-1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[6]; 
  fUp_R[6] = 1.5811388300841895*fEdge[54]-1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[7]; 
  fUp_R[7] = 1.5811388300841895*fEdge[55]-1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[8]; 
  fUp_R[8] = 1.5811388300841895*fEdge[56]-1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[12]; 
  fUp_R[9] = 1.5811388300841895*fEdge[57]-1.224744871391589*fEdge[24]+0.7071067811865475*fEdge[13]; 
  fUp_R[10] = 1.5811388300841895*fEdge[58]-1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[14]; 
  fUp_R[11] = 1.5811388300841898*fEdge[59]-1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[16]; 
  fUp_R[12] = 1.5811388300841898*fEdge[60]-1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[20]; 
  fUp_R[13] = 1.5811388300841898*fEdge[61]-1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[21]; 
  fUp_R[14] = 1.5811388300841898*fEdge[62]-1.224744871391589*fEdge[30]+0.7071067811865475*fEdge[22]; 
  fUp_R[15] = 1.5811388300841895*fEdge[63]-1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[27]; 
   } 
  } else { 
  double fC_r[16] = {0.0}; 
  double fR_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_R[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_R, sgn_drag_coeff_Up_R); 

  fC_r[0] = 1.5811388300841895*fSkin[48]+1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[0]; 
  fC_r[1] = 1.5811388300841898*fSkin[49]+1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[1]; 
  fC_r[2] = 1.5811388300841898*fSkin[50]+1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[2]; 
  fC_r[3] = 1.5811388300841898*fSkin[51]+1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[3]; 
  fC_r[4] = 1.5811388300841898*fSkin[52]+1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[5]; 
  fC_r[5] = 1.5811388300841895*fSkin[53]+1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[6]; 
  fC_r[6] = 1.5811388300841895*fSkin[54]+1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[7]; 
  fC_r[7] = 1.5811388300841895*fSkin[55]+1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[8]; 
  fC_r[8] = 1.5811388300841895*fSkin[56]+1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[12]; 
  fC_r[9] = 1.5811388300841895*fSkin[57]+1.224744871391589*fSkin[24]+0.7071067811865475*fSkin[13]; 
  fC_r[10] = 1.5811388300841895*fSkin[58]+1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[14]; 
  fC_r[11] = 1.5811388300841898*fSkin[59]+1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[16]; 
  fC_r[12] = 1.5811388300841898*fSkin[60]+1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[20]; 
  fC_r[13] = 1.5811388300841898*fSkin[61]+1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[21]; 
  fC_r[14] = 1.5811388300841898*fSkin[62]+1.224744871391589*fSkin[30]+0.7071067811865475*fSkin[22]; 
  fC_r[15] = 1.5811388300841895*fSkin[63]+1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[27]; 

  fR_l[0] = 1.5811388300841895*fEdge[48]-1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[0]; 
  fR_l[1] = 1.5811388300841898*fEdge[49]-1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[1]; 
  fR_l[2] = 1.5811388300841898*fEdge[50]-1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[2]; 
  fR_l[3] = 1.5811388300841898*fEdge[51]-1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[3]; 
  fR_l[4] = 1.5811388300841898*fEdge[52]-1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[5]; 
  fR_l[5] = 1.5811388300841895*fEdge[53]-1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[6]; 
  fR_l[6] = 1.5811388300841895*fEdge[54]-1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[7]; 
  fR_l[7] = 1.5811388300841895*fEdge[55]-1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[8]; 
  fR_l[8] = 1.5811388300841895*fEdge[56]-1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[12]; 
  fR_l[9] = 1.5811388300841895*fEdge[57]-1.224744871391589*fEdge[24]+0.7071067811865475*fEdge[13]; 
  fR_l[10] = 1.5811388300841895*fEdge[58]-1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[14]; 
  fR_l[11] = 1.5811388300841898*fEdge[59]-1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[16]; 
  fR_l[12] = 1.5811388300841898*fEdge[60]-1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[20]; 
  fR_l[13] = 1.5811388300841898*fEdge[61]-1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[21]; 
  fR_l[14] = 1.5811388300841898*fEdge[62]-1.224744871391589*fEdge[30]+0.7071067811865475*fEdge[22]; 
  fR_l[15] = 1.5811388300841895*fEdge[63]-1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[27]; 

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
  double GhatR[16] = {0.0}; 
  GhatR[0] = 0.25*(drag_coeff_surf_Edge[15]*fUp_R[15]+drag_coeff_surf_Edge[14]*fUp_R[14]+drag_coeff_surf_Edge[13]*fUp_R[13]+drag_coeff_surf_Edge[12]*fUp_R[12]+drag_coeff_surf_Edge[11]*fUp_R[11]+drag_coeff_surf_Edge[10]*fUp_R[10]+drag_coeff_surf_Edge[9]*fUp_R[9]+drag_coeff_surf_Edge[8]*fUp_R[8]+drag_coeff_surf_Edge[7]*fUp_R[7]+drag_coeff_surf_Edge[6]*fUp_R[6]+drag_coeff_surf_Edge[5]*fUp_R[5]+drag_coeff_surf_Edge[4]*fUp_R[4]+drag_coeff_surf_Edge[3]*fUp_R[3]+drag_coeff_surf_Edge[2]*fUp_R[2]+drag_coeff_surf_Edge[1]*fUp_R[1]+drag_coeff_surf_Edge[0]*fUp_R[0]); 
  GhatR[1] = 0.25*(drag_coeff_surf_Edge[14]*fUp_R[15]+fUp_R[14]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[10]*fUp_R[13]+fUp_R[10]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[9]*fUp_R[12]+fUp_R[9]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[7]*fUp_R[11]+fUp_R[7]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[4]*fUp_R[8]+fUp_R[4]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[3]*fUp_R[6]+fUp_R[3]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[2]*fUp_R[5]+fUp_R[2]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[0]*fUp_R[1]+fUp_R[0]*drag_coeff_surf_Edge[1]); 
  GhatR[2] = 0.25*(drag_coeff_surf_Edge[13]*fUp_R[15]+fUp_R[13]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[10]*fUp_R[14]+fUp_R[10]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[8]*fUp_R[12]+fUp_R[8]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[6]*fUp_R[11]+fUp_R[6]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[4]*fUp_R[9]+fUp_R[4]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[3]*fUp_R[7]+fUp_R[3]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[1]*fUp_R[5]+fUp_R[1]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[0]*fUp_R[2]+fUp_R[0]*drag_coeff_surf_Edge[2]); 
  GhatR[3] = 0.25*(drag_coeff_surf_Edge[12]*fUp_R[15]+fUp_R[12]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[9]*fUp_R[14]+fUp_R[9]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[8]*fUp_R[13]+fUp_R[8]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[5]*fUp_R[11]+fUp_R[5]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[4]*fUp_R[10]+fUp_R[4]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[2]*fUp_R[7]+fUp_R[2]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[1]*fUp_R[6]+fUp_R[1]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[0]*fUp_R[3]+fUp_R[0]*drag_coeff_surf_Edge[3]); 
  GhatR[4] = 0.25*(drag_coeff_surf_Edge[11]*fUp_R[15]+fUp_R[11]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[7]*fUp_R[14]+fUp_R[7]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[6]*fUp_R[13]+fUp_R[6]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[5]*fUp_R[12]+fUp_R[5]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[3]*fUp_R[10]+fUp_R[3]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[2]*fUp_R[9]+fUp_R[2]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[1]*fUp_R[8]+fUp_R[1]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[0]*fUp_R[4]+fUp_R[0]*drag_coeff_surf_Edge[4]); 
  GhatR[5] = 0.25*(drag_coeff_surf_Edge[10]*fUp_R[15]+fUp_R[10]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[13]*fUp_R[14]+fUp_R[13]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[4]*fUp_R[12]+fUp_R[4]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[3]*fUp_R[11]+fUp_R[3]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[8]*fUp_R[9]+fUp_R[8]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[6]*fUp_R[7]+fUp_R[6]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[0]*fUp_R[5]+fUp_R[0]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[1]*fUp_R[2]+fUp_R[1]*drag_coeff_surf_Edge[2]); 
  GhatR[6] = 0.25*(drag_coeff_surf_Edge[9]*fUp_R[15]+fUp_R[9]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[12]*fUp_R[14]+fUp_R[12]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[4]*fUp_R[13]+fUp_R[4]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[2]*fUp_R[11]+fUp_R[2]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[8]*fUp_R[10]+fUp_R[8]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[5]*fUp_R[7]+fUp_R[5]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[0]*fUp_R[6]+fUp_R[0]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[1]*fUp_R[3]+fUp_R[1]*drag_coeff_surf_Edge[3]); 
  GhatR[7] = 0.25*(drag_coeff_surf_Edge[8]*fUp_R[15]+fUp_R[8]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[4]*fUp_R[14]+fUp_R[4]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[12]*fUp_R[13]+fUp_R[12]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[1]*fUp_R[11]+fUp_R[1]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[9]*fUp_R[10]+fUp_R[9]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[0]*fUp_R[7]+fUp_R[0]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[5]*fUp_R[6]+fUp_R[5]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[2]*fUp_R[3]+fUp_R[2]*drag_coeff_surf_Edge[3]); 
  GhatR[8] = 0.25*(drag_coeff_surf_Edge[7]*fUp_R[15]+fUp_R[7]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[11]*fUp_R[14]+fUp_R[11]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[3]*fUp_R[13]+fUp_R[3]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[2]*fUp_R[12]+fUp_R[2]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[6]*fUp_R[10]+fUp_R[6]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[5]*fUp_R[9]+fUp_R[5]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[0]*fUp_R[8]+fUp_R[0]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[1]*fUp_R[4]+fUp_R[1]*drag_coeff_surf_Edge[4]); 
  GhatR[9] = 0.25*(drag_coeff_surf_Edge[6]*fUp_R[15]+fUp_R[6]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[3]*fUp_R[14]+fUp_R[3]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[11]*fUp_R[13]+fUp_R[11]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[1]*fUp_R[12]+fUp_R[1]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[7]*fUp_R[10]+fUp_R[7]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[0]*fUp_R[9]+fUp_R[0]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[5]*fUp_R[8]+fUp_R[5]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[2]*fUp_R[4]+fUp_R[2]*drag_coeff_surf_Edge[4]); 
  GhatR[10] = 0.25*(drag_coeff_surf_Edge[5]*fUp_R[15]+fUp_R[5]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[2]*fUp_R[14]+fUp_R[2]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[1]*fUp_R[13]+fUp_R[1]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[11]*fUp_R[12]+fUp_R[11]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[0]*fUp_R[10]+fUp_R[0]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[7]*fUp_R[9]+fUp_R[7]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[6]*fUp_R[8]+fUp_R[6]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[3]*fUp_R[4]+fUp_R[3]*drag_coeff_surf_Edge[4]); 
  GhatR[11] = 0.25*(drag_coeff_surf_Edge[4]*fUp_R[15]+fUp_R[4]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[8]*fUp_R[14]+fUp_R[8]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[9]*fUp_R[13]+fUp_R[9]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[10]*fUp_R[12]+fUp_R[10]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[0]*fUp_R[11]+fUp_R[0]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[1]*fUp_R[7]+fUp_R[1]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[2]*fUp_R[6]+fUp_R[2]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[3]*fUp_R[5]+fUp_R[3]*drag_coeff_surf_Edge[5]); 
  GhatR[12] = 0.25*(drag_coeff_surf_Edge[3]*fUp_R[15]+fUp_R[3]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[6]*fUp_R[14]+fUp_R[6]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[7]*fUp_R[13]+fUp_R[7]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[0]*fUp_R[12]+fUp_R[0]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[10]*fUp_R[11]+fUp_R[10]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[1]*fUp_R[9]+fUp_R[1]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[2]*fUp_R[8]+fUp_R[2]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[4]*fUp_R[5]+fUp_R[4]*drag_coeff_surf_Edge[5]); 
  GhatR[13] = 0.25*(drag_coeff_surf_Edge[2]*fUp_R[15]+fUp_R[2]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[5]*fUp_R[14]+fUp_R[5]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[0]*fUp_R[13]+fUp_R[0]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[7]*fUp_R[12]+fUp_R[7]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[9]*fUp_R[11]+fUp_R[9]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[1]*fUp_R[10]+fUp_R[1]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[3]*fUp_R[8]+fUp_R[3]*drag_coeff_surf_Edge[8]+drag_coeff_surf_Edge[4]*fUp_R[6]+fUp_R[4]*drag_coeff_surf_Edge[6]); 
  GhatR[14] = 0.25*(drag_coeff_surf_Edge[1]*fUp_R[15]+fUp_R[1]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[0]*fUp_R[14]+fUp_R[0]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[5]*fUp_R[13]+fUp_R[5]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[6]*fUp_R[12]+fUp_R[6]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[8]*fUp_R[11]+fUp_R[8]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[2]*fUp_R[10]+fUp_R[2]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[3]*fUp_R[9]+fUp_R[3]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[4]*fUp_R[7]+fUp_R[4]*drag_coeff_surf_Edge[7]); 
  GhatR[15] = 0.25*(drag_coeff_surf_Edge[0]*fUp_R[15]+fUp_R[0]*drag_coeff_surf_Edge[15]+drag_coeff_surf_Edge[1]*fUp_R[14]+fUp_R[1]*drag_coeff_surf_Edge[14]+drag_coeff_surf_Edge[2]*fUp_R[13]+fUp_R[2]*drag_coeff_surf_Edge[13]+drag_coeff_surf_Edge[3]*fUp_R[12]+fUp_R[3]*drag_coeff_surf_Edge[12]+drag_coeff_surf_Edge[4]*fUp_R[11]+fUp_R[4]*drag_coeff_surf_Edge[11]+drag_coeff_surf_Edge[5]*fUp_R[10]+fUp_R[5]*drag_coeff_surf_Edge[10]+drag_coeff_surf_Edge[6]*fUp_R[9]+fUp_R[6]*drag_coeff_surf_Edge[9]+drag_coeff_surf_Edge[7]*fUp_R[8]+fUp_R[7]*drag_coeff_surf_Edge[8]); 

  out[0] += -(0.35355339059327373*GhatR[0]*dv_inv); 
  out[1] += -(0.35355339059327373*GhatR[1]*dv_inv); 
  out[2] += -(0.35355339059327373*GhatR[2]*dv_inv); 
  out[3] += -(0.35355339059327373*GhatR[3]*dv_inv); 
  out[4] += -(0.6123724356957945*GhatR[0]*dv_inv); 
  out[5] += -(0.35355339059327373*GhatR[4]*dv_inv); 
  out[6] += -(0.35355339059327373*GhatR[5]*dv_inv); 
  out[7] += -(0.35355339059327373*GhatR[6]*dv_inv); 
  out[8] += -(0.35355339059327373*GhatR[7]*dv_inv); 
  out[9] += -(0.6123724356957945*GhatR[1]*dv_inv); 
  out[10] += -(0.6123724356957945*GhatR[2]*dv_inv); 
  out[11] += -(0.6123724356957945*GhatR[3]*dv_inv); 
  out[12] += -(0.35355339059327373*GhatR[8]*dv_inv); 
  out[13] += -(0.35355339059327373*GhatR[9]*dv_inv); 
  out[14] += -(0.35355339059327373*GhatR[10]*dv_inv); 
  out[15] += -(0.6123724356957945*GhatR[4]*dv_inv); 
  out[16] += -(0.35355339059327373*GhatR[11]*dv_inv); 
  out[17] += -(0.6123724356957945*GhatR[5]*dv_inv); 
  out[18] += -(0.6123724356957945*GhatR[6]*dv_inv); 
  out[19] += -(0.6123724356957945*GhatR[7]*dv_inv); 
  out[20] += -(0.35355339059327373*GhatR[12]*dv_inv); 
  out[21] += -(0.35355339059327373*GhatR[13]*dv_inv); 
  out[22] += -(0.35355339059327373*GhatR[14]*dv_inv); 
  out[23] += -(0.6123724356957945*GhatR[8]*dv_inv); 
  out[24] += -(0.6123724356957945*GhatR[9]*dv_inv); 
  out[25] += -(0.6123724356957945*GhatR[10]*dv_inv); 
  out[26] += -(0.6123724356957945*GhatR[11]*dv_inv); 
  out[27] += -(0.35355339059327373*GhatR[15]*dv_inv); 
  out[28] += -(0.6123724356957945*GhatR[12]*dv_inv); 
  out[29] += -(0.6123724356957945*GhatR[13]*dv_inv); 
  out[30] += -(0.6123724356957945*GhatR[14]*dv_inv); 
  out[31] += -(0.6123724356957945*GhatR[15]*dv_inv); 
  out[48] += -(0.7905694150420948*GhatR[0]*dv_inv); 
  out[49] += -(0.7905694150420949*GhatR[1]*dv_inv); 
  out[50] += -(0.7905694150420949*GhatR[2]*dv_inv); 
  out[51] += -(0.7905694150420949*GhatR[3]*dv_inv); 
  out[52] += -(0.7905694150420949*GhatR[4]*dv_inv); 
  out[53] += -(0.7905694150420948*GhatR[5]*dv_inv); 
  out[54] += -(0.7905694150420948*GhatR[6]*dv_inv); 
  out[55] += -(0.7905694150420948*GhatR[7]*dv_inv); 
  out[56] += -(0.7905694150420948*GhatR[8]*dv_inv); 
  out[57] += -(0.7905694150420948*GhatR[9]*dv_inv); 
  out[58] += -(0.7905694150420948*GhatR[10]*dv_inv); 
  out[59] += -(0.7905694150420949*GhatR[11]*dv_inv); 
  out[60] += -(0.7905694150420949*GhatR[12]*dv_inv); 
  out[61] += -(0.7905694150420949*GhatR[13]*dv_inv); 
  out[62] += -(0.7905694150420949*GhatR[14]*dv_inv); 
  out[63] += -(0.7905694150420948*GhatR[15]*dv_inv); 
  cflFreq = fabs(drag_coeff_surf_Edge[0]);

  } else { 
  const double *drag_coeff_surf_L = drag_coeff_surf_Skin; 
  const double *sgn_drag_coeff_surf_L = sgn_drag_coeff_surf_Skin; 
  const int *const_sgn_drag_coeff_L = const_sgn_drag_coeff_Skin; 
  double fUp_L[16] = {0.0}; 
  if (const_sgn_drag_coeff_L[0] == 1) { 
    if (sgn_drag_coeff_surf_L[0] == 1.0) { 
  fUp_L[0] = 1.5811388300841895*fEdge[48]+1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[0]; 
  fUp_L[1] = 1.5811388300841898*fEdge[49]+1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[1]; 
  fUp_L[2] = 1.5811388300841898*fEdge[50]+1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[2]; 
  fUp_L[3] = 1.5811388300841898*fEdge[51]+1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[3]; 
  fUp_L[4] = 1.5811388300841898*fEdge[52]+1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[5]; 
  fUp_L[5] = 1.5811388300841895*fEdge[53]+1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[6]; 
  fUp_L[6] = 1.5811388300841895*fEdge[54]+1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[7]; 
  fUp_L[7] = 1.5811388300841895*fEdge[55]+1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[8]; 
  fUp_L[8] = 1.5811388300841895*fEdge[56]+1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[12]; 
  fUp_L[9] = 1.5811388300841895*fEdge[57]+1.224744871391589*fEdge[24]+0.7071067811865475*fEdge[13]; 
  fUp_L[10] = 1.5811388300841895*fEdge[58]+1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[14]; 
  fUp_L[11] = 1.5811388300841898*fEdge[59]+1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[16]; 
  fUp_L[12] = 1.5811388300841898*fEdge[60]+1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[20]; 
  fUp_L[13] = 1.5811388300841898*fEdge[61]+1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[21]; 
  fUp_L[14] = 1.5811388300841898*fEdge[62]+1.224744871391589*fEdge[30]+0.7071067811865475*fEdge[22]; 
  fUp_L[15] = 1.5811388300841895*fEdge[63]+1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[27]; 
    } else { 
  fUp_L[0] = 1.5811388300841895*fSkin[48]-1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[0]; 
  fUp_L[1] = 1.5811388300841898*fSkin[49]-1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[1]; 
  fUp_L[2] = 1.5811388300841898*fSkin[50]-1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[2]; 
  fUp_L[3] = 1.5811388300841898*fSkin[51]-1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[3]; 
  fUp_L[4] = 1.5811388300841898*fSkin[52]-1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[5]; 
  fUp_L[5] = 1.5811388300841895*fSkin[53]-1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[6]; 
  fUp_L[6] = 1.5811388300841895*fSkin[54]-1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[7]; 
  fUp_L[7] = 1.5811388300841895*fSkin[55]-1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[8]; 
  fUp_L[8] = 1.5811388300841895*fSkin[56]-1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[12]; 
  fUp_L[9] = 1.5811388300841895*fSkin[57]-1.224744871391589*fSkin[24]+0.7071067811865475*fSkin[13]; 
  fUp_L[10] = 1.5811388300841895*fSkin[58]-1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[14]; 
  fUp_L[11] = 1.5811388300841898*fSkin[59]-1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[16]; 
  fUp_L[12] = 1.5811388300841898*fSkin[60]-1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[20]; 
  fUp_L[13] = 1.5811388300841898*fSkin[61]-1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[21]; 
  fUp_L[14] = 1.5811388300841898*fSkin[62]-1.224744871391589*fSkin[30]+0.7071067811865475*fSkin[22]; 
  fUp_L[15] = 1.5811388300841895*fSkin[63]-1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[27]; 
   } 
  } else { 
  double fL_r[16] = {0.0}; 
  double fC_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_L[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_L, sgn_drag_coeff_Up_L); 

  fL_r[0] = 1.5811388300841895*fEdge[48]+1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[0]; 
  fL_r[1] = 1.5811388300841898*fEdge[49]+1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[1]; 
  fL_r[2] = 1.5811388300841898*fEdge[50]+1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[2]; 
  fL_r[3] = 1.5811388300841898*fEdge[51]+1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[3]; 
  fL_r[4] = 1.5811388300841898*fEdge[52]+1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[5]; 
  fL_r[5] = 1.5811388300841895*fEdge[53]+1.224744871391589*fEdge[17]+0.7071067811865475*fEdge[6]; 
  fL_r[6] = 1.5811388300841895*fEdge[54]+1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[7]; 
  fL_r[7] = 1.5811388300841895*fEdge[55]+1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[8]; 
  fL_r[8] = 1.5811388300841895*fEdge[56]+1.224744871391589*fEdge[23]+0.7071067811865475*fEdge[12]; 
  fL_r[9] = 1.5811388300841895*fEdge[57]+1.224744871391589*fEdge[24]+0.7071067811865475*fEdge[13]; 
  fL_r[10] = 1.5811388300841895*fEdge[58]+1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[14]; 
  fL_r[11] = 1.5811388300841898*fEdge[59]+1.224744871391589*fEdge[26]+0.7071067811865475*fEdge[16]; 
  fL_r[12] = 1.5811388300841898*fEdge[60]+1.224744871391589*fEdge[28]+0.7071067811865475*fEdge[20]; 
  fL_r[13] = 1.5811388300841898*fEdge[61]+1.224744871391589*fEdge[29]+0.7071067811865475*fEdge[21]; 
  fL_r[14] = 1.5811388300841898*fEdge[62]+1.224744871391589*fEdge[30]+0.7071067811865475*fEdge[22]; 
  fL_r[15] = 1.5811388300841895*fEdge[63]+1.224744871391589*fEdge[31]+0.7071067811865475*fEdge[27]; 

  fC_l[0] = 1.5811388300841895*fSkin[48]-1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[0]; 
  fC_l[1] = 1.5811388300841898*fSkin[49]-1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[1]; 
  fC_l[2] = 1.5811388300841898*fSkin[50]-1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[2]; 
  fC_l[3] = 1.5811388300841898*fSkin[51]-1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[3]; 
  fC_l[4] = 1.5811388300841898*fSkin[52]-1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[5]; 
  fC_l[5] = 1.5811388300841895*fSkin[53]-1.224744871391589*fSkin[17]+0.7071067811865475*fSkin[6]; 
  fC_l[6] = 1.5811388300841895*fSkin[54]-1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[7]; 
  fC_l[7] = 1.5811388300841895*fSkin[55]-1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[8]; 
  fC_l[8] = 1.5811388300841895*fSkin[56]-1.224744871391589*fSkin[23]+0.7071067811865475*fSkin[12]; 
  fC_l[9] = 1.5811388300841895*fSkin[57]-1.224744871391589*fSkin[24]+0.7071067811865475*fSkin[13]; 
  fC_l[10] = 1.5811388300841895*fSkin[58]-1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[14]; 
  fC_l[11] = 1.5811388300841898*fSkin[59]-1.224744871391589*fSkin[26]+0.7071067811865475*fSkin[16]; 
  fC_l[12] = 1.5811388300841898*fSkin[60]-1.224744871391589*fSkin[28]+0.7071067811865475*fSkin[20]; 
  fC_l[13] = 1.5811388300841898*fSkin[61]-1.224744871391589*fSkin[29]+0.7071067811865475*fSkin[21]; 
  fC_l[14] = 1.5811388300841898*fSkin[62]-1.224744871391589*fSkin[30]+0.7071067811865475*fSkin[22]; 
  fC_l[15] = 1.5811388300841895*fSkin[63]-1.224744871391589*fSkin[31]+0.7071067811865475*fSkin[27]; 

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
  double GhatL[16] = {0.0}; 
  GhatL[0] = 0.25*(drag_coeff_surf_Skin[15]*fUp_L[15]+drag_coeff_surf_Skin[14]*fUp_L[14]+drag_coeff_surf_Skin[13]*fUp_L[13]+drag_coeff_surf_Skin[12]*fUp_L[12]+drag_coeff_surf_Skin[11]*fUp_L[11]+drag_coeff_surf_Skin[10]*fUp_L[10]+drag_coeff_surf_Skin[9]*fUp_L[9]+drag_coeff_surf_Skin[8]*fUp_L[8]+drag_coeff_surf_Skin[7]*fUp_L[7]+drag_coeff_surf_Skin[6]*fUp_L[6]+drag_coeff_surf_Skin[5]*fUp_L[5]+drag_coeff_surf_Skin[4]*fUp_L[4]+drag_coeff_surf_Skin[3]*fUp_L[3]+drag_coeff_surf_Skin[2]*fUp_L[2]+drag_coeff_surf_Skin[1]*fUp_L[1]+drag_coeff_surf_Skin[0]*fUp_L[0]); 
  GhatL[1] = 0.25*(drag_coeff_surf_Skin[14]*fUp_L[15]+fUp_L[14]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[10]*fUp_L[13]+fUp_L[10]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[9]*fUp_L[12]+fUp_L[9]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[7]*fUp_L[11]+fUp_L[7]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[4]*fUp_L[8]+fUp_L[4]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[3]*fUp_L[6]+fUp_L[3]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[2]*fUp_L[5]+fUp_L[2]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[0]*fUp_L[1]+fUp_L[0]*drag_coeff_surf_Skin[1]); 
  GhatL[2] = 0.25*(drag_coeff_surf_Skin[13]*fUp_L[15]+fUp_L[13]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[10]*fUp_L[14]+fUp_L[10]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[8]*fUp_L[12]+fUp_L[8]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[6]*fUp_L[11]+fUp_L[6]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[4]*fUp_L[9]+fUp_L[4]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[3]*fUp_L[7]+fUp_L[3]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[1]*fUp_L[5]+fUp_L[1]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[0]*fUp_L[2]+fUp_L[0]*drag_coeff_surf_Skin[2]); 
  GhatL[3] = 0.25*(drag_coeff_surf_Skin[12]*fUp_L[15]+fUp_L[12]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[9]*fUp_L[14]+fUp_L[9]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[8]*fUp_L[13]+fUp_L[8]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[5]*fUp_L[11]+fUp_L[5]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[4]*fUp_L[10]+fUp_L[4]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[2]*fUp_L[7]+fUp_L[2]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[1]*fUp_L[6]+fUp_L[1]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[0]*fUp_L[3]+fUp_L[0]*drag_coeff_surf_Skin[3]); 
  GhatL[4] = 0.25*(drag_coeff_surf_Skin[11]*fUp_L[15]+fUp_L[11]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[7]*fUp_L[14]+fUp_L[7]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[6]*fUp_L[13]+fUp_L[6]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[5]*fUp_L[12]+fUp_L[5]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[3]*fUp_L[10]+fUp_L[3]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[2]*fUp_L[9]+fUp_L[2]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[1]*fUp_L[8]+fUp_L[1]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[0]*fUp_L[4]+fUp_L[0]*drag_coeff_surf_Skin[4]); 
  GhatL[5] = 0.25*(drag_coeff_surf_Skin[10]*fUp_L[15]+fUp_L[10]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[13]*fUp_L[14]+fUp_L[13]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[4]*fUp_L[12]+fUp_L[4]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[3]*fUp_L[11]+fUp_L[3]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[8]*fUp_L[9]+fUp_L[8]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[6]*fUp_L[7]+fUp_L[6]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[0]*fUp_L[5]+fUp_L[0]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[1]*fUp_L[2]+fUp_L[1]*drag_coeff_surf_Skin[2]); 
  GhatL[6] = 0.25*(drag_coeff_surf_Skin[9]*fUp_L[15]+fUp_L[9]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[12]*fUp_L[14]+fUp_L[12]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[4]*fUp_L[13]+fUp_L[4]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[2]*fUp_L[11]+fUp_L[2]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[8]*fUp_L[10]+fUp_L[8]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[5]*fUp_L[7]+fUp_L[5]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[0]*fUp_L[6]+fUp_L[0]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[1]*fUp_L[3]+fUp_L[1]*drag_coeff_surf_Skin[3]); 
  GhatL[7] = 0.25*(drag_coeff_surf_Skin[8]*fUp_L[15]+fUp_L[8]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[4]*fUp_L[14]+fUp_L[4]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[12]*fUp_L[13]+fUp_L[12]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[1]*fUp_L[11]+fUp_L[1]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[9]*fUp_L[10]+fUp_L[9]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[0]*fUp_L[7]+fUp_L[0]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[5]*fUp_L[6]+fUp_L[5]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[2]*fUp_L[3]+fUp_L[2]*drag_coeff_surf_Skin[3]); 
  GhatL[8] = 0.25*(drag_coeff_surf_Skin[7]*fUp_L[15]+fUp_L[7]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[11]*fUp_L[14]+fUp_L[11]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[3]*fUp_L[13]+fUp_L[3]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[2]*fUp_L[12]+fUp_L[2]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[6]*fUp_L[10]+fUp_L[6]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[5]*fUp_L[9]+fUp_L[5]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[0]*fUp_L[8]+fUp_L[0]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[1]*fUp_L[4]+fUp_L[1]*drag_coeff_surf_Skin[4]); 
  GhatL[9] = 0.25*(drag_coeff_surf_Skin[6]*fUp_L[15]+fUp_L[6]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[3]*fUp_L[14]+fUp_L[3]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[11]*fUp_L[13]+fUp_L[11]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[1]*fUp_L[12]+fUp_L[1]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[7]*fUp_L[10]+fUp_L[7]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[0]*fUp_L[9]+fUp_L[0]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[5]*fUp_L[8]+fUp_L[5]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[2]*fUp_L[4]+fUp_L[2]*drag_coeff_surf_Skin[4]); 
  GhatL[10] = 0.25*(drag_coeff_surf_Skin[5]*fUp_L[15]+fUp_L[5]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[2]*fUp_L[14]+fUp_L[2]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[1]*fUp_L[13]+fUp_L[1]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[11]*fUp_L[12]+fUp_L[11]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[0]*fUp_L[10]+fUp_L[0]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[7]*fUp_L[9]+fUp_L[7]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[6]*fUp_L[8]+fUp_L[6]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[3]*fUp_L[4]+fUp_L[3]*drag_coeff_surf_Skin[4]); 
  GhatL[11] = 0.25*(drag_coeff_surf_Skin[4]*fUp_L[15]+fUp_L[4]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[8]*fUp_L[14]+fUp_L[8]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[9]*fUp_L[13]+fUp_L[9]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[10]*fUp_L[12]+fUp_L[10]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[0]*fUp_L[11]+fUp_L[0]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[1]*fUp_L[7]+fUp_L[1]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[2]*fUp_L[6]+fUp_L[2]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[3]*fUp_L[5]+fUp_L[3]*drag_coeff_surf_Skin[5]); 
  GhatL[12] = 0.25*(drag_coeff_surf_Skin[3]*fUp_L[15]+fUp_L[3]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[6]*fUp_L[14]+fUp_L[6]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[7]*fUp_L[13]+fUp_L[7]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[0]*fUp_L[12]+fUp_L[0]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[10]*fUp_L[11]+fUp_L[10]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[1]*fUp_L[9]+fUp_L[1]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[2]*fUp_L[8]+fUp_L[2]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[4]*fUp_L[5]+fUp_L[4]*drag_coeff_surf_Skin[5]); 
  GhatL[13] = 0.25*(drag_coeff_surf_Skin[2]*fUp_L[15]+fUp_L[2]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[5]*fUp_L[14]+fUp_L[5]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[0]*fUp_L[13]+fUp_L[0]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[7]*fUp_L[12]+fUp_L[7]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[9]*fUp_L[11]+fUp_L[9]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[1]*fUp_L[10]+fUp_L[1]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[3]*fUp_L[8]+fUp_L[3]*drag_coeff_surf_Skin[8]+drag_coeff_surf_Skin[4]*fUp_L[6]+fUp_L[4]*drag_coeff_surf_Skin[6]); 
  GhatL[14] = 0.25*(drag_coeff_surf_Skin[1]*fUp_L[15]+fUp_L[1]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[0]*fUp_L[14]+fUp_L[0]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[5]*fUp_L[13]+fUp_L[5]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[6]*fUp_L[12]+fUp_L[6]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[8]*fUp_L[11]+fUp_L[8]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[2]*fUp_L[10]+fUp_L[2]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[3]*fUp_L[9]+fUp_L[3]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[4]*fUp_L[7]+fUp_L[4]*drag_coeff_surf_Skin[7]); 
  GhatL[15] = 0.25*(drag_coeff_surf_Skin[0]*fUp_L[15]+fUp_L[0]*drag_coeff_surf_Skin[15]+drag_coeff_surf_Skin[1]*fUp_L[14]+fUp_L[1]*drag_coeff_surf_Skin[14]+drag_coeff_surf_Skin[2]*fUp_L[13]+fUp_L[2]*drag_coeff_surf_Skin[13]+drag_coeff_surf_Skin[3]*fUp_L[12]+fUp_L[3]*drag_coeff_surf_Skin[12]+drag_coeff_surf_Skin[4]*fUp_L[11]+fUp_L[4]*drag_coeff_surf_Skin[11]+drag_coeff_surf_Skin[5]*fUp_L[10]+fUp_L[5]*drag_coeff_surf_Skin[10]+drag_coeff_surf_Skin[6]*fUp_L[9]+fUp_L[6]*drag_coeff_surf_Skin[9]+drag_coeff_surf_Skin[7]*fUp_L[8]+fUp_L[7]*drag_coeff_surf_Skin[8]); 

  out[0] += 0.35355339059327373*GhatL[0]*dv_inv; 
  out[1] += 0.35355339059327373*GhatL[1]*dv_inv; 
  out[2] += 0.35355339059327373*GhatL[2]*dv_inv; 
  out[3] += 0.35355339059327373*GhatL[3]*dv_inv; 
  out[4] += -(0.6123724356957945*GhatL[0]*dv_inv); 
  out[5] += 0.35355339059327373*GhatL[4]*dv_inv; 
  out[6] += 0.35355339059327373*GhatL[5]*dv_inv; 
  out[7] += 0.35355339059327373*GhatL[6]*dv_inv; 
  out[8] += 0.35355339059327373*GhatL[7]*dv_inv; 
  out[9] += -(0.6123724356957945*GhatL[1]*dv_inv); 
  out[10] += -(0.6123724356957945*GhatL[2]*dv_inv); 
  out[11] += -(0.6123724356957945*GhatL[3]*dv_inv); 
  out[12] += 0.35355339059327373*GhatL[8]*dv_inv; 
  out[13] += 0.35355339059327373*GhatL[9]*dv_inv; 
  out[14] += 0.35355339059327373*GhatL[10]*dv_inv; 
  out[15] += -(0.6123724356957945*GhatL[4]*dv_inv); 
  out[16] += 0.35355339059327373*GhatL[11]*dv_inv; 
  out[17] += -(0.6123724356957945*GhatL[5]*dv_inv); 
  out[18] += -(0.6123724356957945*GhatL[6]*dv_inv); 
  out[19] += -(0.6123724356957945*GhatL[7]*dv_inv); 
  out[20] += 0.35355339059327373*GhatL[12]*dv_inv; 
  out[21] += 0.35355339059327373*GhatL[13]*dv_inv; 
  out[22] += 0.35355339059327373*GhatL[14]*dv_inv; 
  out[23] += -(0.6123724356957945*GhatL[8]*dv_inv); 
  out[24] += -(0.6123724356957945*GhatL[9]*dv_inv); 
  out[25] += -(0.6123724356957945*GhatL[10]*dv_inv); 
  out[26] += -(0.6123724356957945*GhatL[11]*dv_inv); 
  out[27] += 0.35355339059327373*GhatL[15]*dv_inv; 
  out[28] += -(0.6123724356957945*GhatL[12]*dv_inv); 
  out[29] += -(0.6123724356957945*GhatL[13]*dv_inv); 
  out[30] += -(0.6123724356957945*GhatL[14]*dv_inv); 
  out[31] += -(0.6123724356957945*GhatL[15]*dv_inv); 
  out[48] += 0.7905694150420948*GhatL[0]*dv_inv; 
  out[49] += 0.7905694150420949*GhatL[1]*dv_inv; 
  out[50] += 0.7905694150420949*GhatL[2]*dv_inv; 
  out[51] += 0.7905694150420949*GhatL[3]*dv_inv; 
  out[52] += 0.7905694150420949*GhatL[4]*dv_inv; 
  out[53] += 0.7905694150420948*GhatL[5]*dv_inv; 
  out[54] += 0.7905694150420948*GhatL[6]*dv_inv; 
  out[55] += 0.7905694150420948*GhatL[7]*dv_inv; 
  out[56] += 0.7905694150420948*GhatL[8]*dv_inv; 
  out[57] += 0.7905694150420948*GhatL[9]*dv_inv; 
  out[58] += 0.7905694150420948*GhatL[10]*dv_inv; 
  out[59] += 0.7905694150420949*GhatL[11]*dv_inv; 
  out[60] += 0.7905694150420949*GhatL[12]*dv_inv; 
  out[61] += 0.7905694150420949*GhatL[13]*dv_inv; 
  out[62] += 0.7905694150420949*GhatL[14]*dv_inv; 
  out[63] += 0.7905694150420948*GhatL[15]*dv_inv; 
  cflFreq = fabs(drag_coeff_surf_Skin[0]);

  } 

  return 0.375*dv_inv*cflFreq; 
}
