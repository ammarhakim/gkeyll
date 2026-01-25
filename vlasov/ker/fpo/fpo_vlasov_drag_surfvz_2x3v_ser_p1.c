#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 


GKYL_CU_DH double fpo_vlasov_drag_surfvz_2x3v_ser_p1(const double* dxv, 
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
  const double *drag_coeff_surf_L = &alpha_surf_L[32]; 
  const double *drag_coeff_surf_R = &alpha_surf_R[32]; 
  const double *sgn_drag_coeff_surf_L = &sgn_alpha_surf_L[32]; 
  const double *sgn_drag_coeff_surf_R = &sgn_alpha_surf_R[32]; 
  const int *const_sgn_drag_coeff_L = &const_sgn_alpha_L[2]; 
  const int *const_sgn_drag_coeff_R = &const_sgn_alpha_R[2]; 
  double dv_inv = 2.0/dxv[4]; 


  double fUp_L[16] = {0.0}; 
  if (const_sgn_drag_coeff_L[0] == 1) { 
    if (sgn_drag_coeff_surf_L[0] == 1.0) { 
  fUp_L[0] = 1.5811388300841895*fL[64]+1.224744871391589*fL[5]+0.7071067811865475*fL[0]; 
  fUp_L[1] = 1.5811388300841898*fL[65]+1.224744871391589*fL[12]+0.7071067811865475*fL[1]; 
  fUp_L[2] = 1.5811388300841898*fL[66]+1.224744871391589*fL[13]+0.7071067811865475*fL[2]; 
  fUp_L[3] = 1.5811388300841898*fL[67]+1.224744871391589*fL[14]+0.7071067811865475*fL[3]; 
  fUp_L[4] = 1.5811388300841898*fL[68]+1.224744871391589*fL[15]+0.7071067811865475*fL[4]; 
  fUp_L[5] = 1.5811388300841895*fL[69]+1.224744871391589*fL[20]+0.7071067811865475*fL[6]; 
  fUp_L[6] = 1.5811388300841895*fL[70]+1.224744871391589*fL[21]+0.7071067811865475*fL[7]; 
  fUp_L[7] = 1.5811388300841895*fL[71]+1.224744871391589*fL[22]+0.7071067811865475*fL[8]; 
  fUp_L[8] = 1.5811388300841895*fL[72]+1.224744871391589*fL[23]+0.7071067811865475*fL[9]; 
  fUp_L[9] = 1.5811388300841895*fL[73]+1.224744871391589*fL[24]+0.7071067811865475*fL[10]; 
  fUp_L[10] = 1.5811388300841895*fL[74]+1.224744871391589*fL[25]+0.7071067811865475*fL[11]; 
  fUp_L[11] = 1.5811388300841898*fL[75]+1.224744871391589*fL[27]+0.7071067811865475*fL[16]; 
  fUp_L[12] = 1.5811388300841898*fL[76]+1.224744871391589*fL[28]+0.7071067811865475*fL[17]; 
  fUp_L[13] = 1.5811388300841898*fL[77]+1.224744871391589*fL[29]+0.7071067811865475*fL[18]; 
  fUp_L[14] = 1.5811388300841898*fL[78]+1.224744871391589*fL[30]+0.7071067811865475*fL[19]; 
  fUp_L[15] = 1.5811388300841895*fL[79]+1.224744871391589*fL[31]+0.7071067811865475*fL[26]; 
    } else { 
  fUp_L[0] = 1.5811388300841895*fC[64]-1.224744871391589*fC[5]+0.7071067811865475*fC[0]; 
  fUp_L[1] = 1.5811388300841898*fC[65]-1.224744871391589*fC[12]+0.7071067811865475*fC[1]; 
  fUp_L[2] = 1.5811388300841898*fC[66]-1.224744871391589*fC[13]+0.7071067811865475*fC[2]; 
  fUp_L[3] = 1.5811388300841898*fC[67]-1.224744871391589*fC[14]+0.7071067811865475*fC[3]; 
  fUp_L[4] = 1.5811388300841898*fC[68]-1.224744871391589*fC[15]+0.7071067811865475*fC[4]; 
  fUp_L[5] = 1.5811388300841895*fC[69]-1.224744871391589*fC[20]+0.7071067811865475*fC[6]; 
  fUp_L[6] = 1.5811388300841895*fC[70]-1.224744871391589*fC[21]+0.7071067811865475*fC[7]; 
  fUp_L[7] = 1.5811388300841895*fC[71]-1.224744871391589*fC[22]+0.7071067811865475*fC[8]; 
  fUp_L[8] = 1.5811388300841895*fC[72]-1.224744871391589*fC[23]+0.7071067811865475*fC[9]; 
  fUp_L[9] = 1.5811388300841895*fC[73]-1.224744871391589*fC[24]+0.7071067811865475*fC[10]; 
  fUp_L[10] = 1.5811388300841895*fC[74]-1.224744871391589*fC[25]+0.7071067811865475*fC[11]; 
  fUp_L[11] = 1.5811388300841898*fC[75]-1.224744871391589*fC[27]+0.7071067811865475*fC[16]; 
  fUp_L[12] = 1.5811388300841898*fC[76]-1.224744871391589*fC[28]+0.7071067811865475*fC[17]; 
  fUp_L[13] = 1.5811388300841898*fC[77]-1.224744871391589*fC[29]+0.7071067811865475*fC[18]; 
  fUp_L[14] = 1.5811388300841898*fC[78]-1.224744871391589*fC[30]+0.7071067811865475*fC[19]; 
  fUp_L[15] = 1.5811388300841895*fC[79]-1.224744871391589*fC[31]+0.7071067811865475*fC[26]; 
   } 
  } else { 
  double fL_r[16] = {0.0}; 
  double fC_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_L[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_L, sgn_drag_coeff_Up_L); 

  fL_r[0] = 1.5811388300841895*fL[64]+1.224744871391589*fL[5]+0.7071067811865475*fL[0]; 
  fL_r[1] = 1.5811388300841898*fL[65]+1.224744871391589*fL[12]+0.7071067811865475*fL[1]; 
  fL_r[2] = 1.5811388300841898*fL[66]+1.224744871391589*fL[13]+0.7071067811865475*fL[2]; 
  fL_r[3] = 1.5811388300841898*fL[67]+1.224744871391589*fL[14]+0.7071067811865475*fL[3]; 
  fL_r[4] = 1.5811388300841898*fL[68]+1.224744871391589*fL[15]+0.7071067811865475*fL[4]; 
  fL_r[5] = 1.5811388300841895*fL[69]+1.224744871391589*fL[20]+0.7071067811865475*fL[6]; 
  fL_r[6] = 1.5811388300841895*fL[70]+1.224744871391589*fL[21]+0.7071067811865475*fL[7]; 
  fL_r[7] = 1.5811388300841895*fL[71]+1.224744871391589*fL[22]+0.7071067811865475*fL[8]; 
  fL_r[8] = 1.5811388300841895*fL[72]+1.224744871391589*fL[23]+0.7071067811865475*fL[9]; 
  fL_r[9] = 1.5811388300841895*fL[73]+1.224744871391589*fL[24]+0.7071067811865475*fL[10]; 
  fL_r[10] = 1.5811388300841895*fL[74]+1.224744871391589*fL[25]+0.7071067811865475*fL[11]; 
  fL_r[11] = 1.5811388300841898*fL[75]+1.224744871391589*fL[27]+0.7071067811865475*fL[16]; 
  fL_r[12] = 1.5811388300841898*fL[76]+1.224744871391589*fL[28]+0.7071067811865475*fL[17]; 
  fL_r[13] = 1.5811388300841898*fL[77]+1.224744871391589*fL[29]+0.7071067811865475*fL[18]; 
  fL_r[14] = 1.5811388300841898*fL[78]+1.224744871391589*fL[30]+0.7071067811865475*fL[19]; 
  fL_r[15] = 1.5811388300841895*fL[79]+1.224744871391589*fL[31]+0.7071067811865475*fL[26]; 

  fC_l[0] = 1.5811388300841895*fC[64]-1.224744871391589*fC[5]+0.7071067811865475*fC[0]; 
  fC_l[1] = 1.5811388300841898*fC[65]-1.224744871391589*fC[12]+0.7071067811865475*fC[1]; 
  fC_l[2] = 1.5811388300841898*fC[66]-1.224744871391589*fC[13]+0.7071067811865475*fC[2]; 
  fC_l[3] = 1.5811388300841898*fC[67]-1.224744871391589*fC[14]+0.7071067811865475*fC[3]; 
  fC_l[4] = 1.5811388300841898*fC[68]-1.224744871391589*fC[15]+0.7071067811865475*fC[4]; 
  fC_l[5] = 1.5811388300841895*fC[69]-1.224744871391589*fC[20]+0.7071067811865475*fC[6]; 
  fC_l[6] = 1.5811388300841895*fC[70]-1.224744871391589*fC[21]+0.7071067811865475*fC[7]; 
  fC_l[7] = 1.5811388300841895*fC[71]-1.224744871391589*fC[22]+0.7071067811865475*fC[8]; 
  fC_l[8] = 1.5811388300841895*fC[72]-1.224744871391589*fC[23]+0.7071067811865475*fC[9]; 
  fC_l[9] = 1.5811388300841895*fC[73]-1.224744871391589*fC[24]+0.7071067811865475*fC[10]; 
  fC_l[10] = 1.5811388300841895*fC[74]-1.224744871391589*fC[25]+0.7071067811865475*fC[11]; 
  fC_l[11] = 1.5811388300841898*fC[75]-1.224744871391589*fC[27]+0.7071067811865475*fC[16]; 
  fC_l[12] = 1.5811388300841898*fC[76]-1.224744871391589*fC[28]+0.7071067811865475*fC[17]; 
  fC_l[13] = 1.5811388300841898*fC[77]-1.224744871391589*fC[29]+0.7071067811865475*fC[18]; 
  fC_l[14] = 1.5811388300841898*fC[78]-1.224744871391589*fC[30]+0.7071067811865475*fC[19]; 
  fC_l[15] = 1.5811388300841895*fC[79]-1.224744871391589*fC[31]+0.7071067811865475*fC[26]; 

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
  fUp_R[0] = 1.5811388300841895*fC[64]+1.224744871391589*fC[5]+0.7071067811865475*fC[0]; 
  fUp_R[1] = 1.5811388300841898*fC[65]+1.224744871391589*fC[12]+0.7071067811865475*fC[1]; 
  fUp_R[2] = 1.5811388300841898*fC[66]+1.224744871391589*fC[13]+0.7071067811865475*fC[2]; 
  fUp_R[3] = 1.5811388300841898*fC[67]+1.224744871391589*fC[14]+0.7071067811865475*fC[3]; 
  fUp_R[4] = 1.5811388300841898*fC[68]+1.224744871391589*fC[15]+0.7071067811865475*fC[4]; 
  fUp_R[5] = 1.5811388300841895*fC[69]+1.224744871391589*fC[20]+0.7071067811865475*fC[6]; 
  fUp_R[6] = 1.5811388300841895*fC[70]+1.224744871391589*fC[21]+0.7071067811865475*fC[7]; 
  fUp_R[7] = 1.5811388300841895*fC[71]+1.224744871391589*fC[22]+0.7071067811865475*fC[8]; 
  fUp_R[8] = 1.5811388300841895*fC[72]+1.224744871391589*fC[23]+0.7071067811865475*fC[9]; 
  fUp_R[9] = 1.5811388300841895*fC[73]+1.224744871391589*fC[24]+0.7071067811865475*fC[10]; 
  fUp_R[10] = 1.5811388300841895*fC[74]+1.224744871391589*fC[25]+0.7071067811865475*fC[11]; 
  fUp_R[11] = 1.5811388300841898*fC[75]+1.224744871391589*fC[27]+0.7071067811865475*fC[16]; 
  fUp_R[12] = 1.5811388300841898*fC[76]+1.224744871391589*fC[28]+0.7071067811865475*fC[17]; 
  fUp_R[13] = 1.5811388300841898*fC[77]+1.224744871391589*fC[29]+0.7071067811865475*fC[18]; 
  fUp_R[14] = 1.5811388300841898*fC[78]+1.224744871391589*fC[30]+0.7071067811865475*fC[19]; 
  fUp_R[15] = 1.5811388300841895*fC[79]+1.224744871391589*fC[31]+0.7071067811865475*fC[26]; 
    } else { 
  fUp_R[0] = 1.5811388300841895*fR[64]-1.224744871391589*fR[5]+0.7071067811865475*fR[0]; 
  fUp_R[1] = 1.5811388300841898*fR[65]-1.224744871391589*fR[12]+0.7071067811865475*fR[1]; 
  fUp_R[2] = 1.5811388300841898*fR[66]-1.224744871391589*fR[13]+0.7071067811865475*fR[2]; 
  fUp_R[3] = 1.5811388300841898*fR[67]-1.224744871391589*fR[14]+0.7071067811865475*fR[3]; 
  fUp_R[4] = 1.5811388300841898*fR[68]-1.224744871391589*fR[15]+0.7071067811865475*fR[4]; 
  fUp_R[5] = 1.5811388300841895*fR[69]-1.224744871391589*fR[20]+0.7071067811865475*fR[6]; 
  fUp_R[6] = 1.5811388300841895*fR[70]-1.224744871391589*fR[21]+0.7071067811865475*fR[7]; 
  fUp_R[7] = 1.5811388300841895*fR[71]-1.224744871391589*fR[22]+0.7071067811865475*fR[8]; 
  fUp_R[8] = 1.5811388300841895*fR[72]-1.224744871391589*fR[23]+0.7071067811865475*fR[9]; 
  fUp_R[9] = 1.5811388300841895*fR[73]-1.224744871391589*fR[24]+0.7071067811865475*fR[10]; 
  fUp_R[10] = 1.5811388300841895*fR[74]-1.224744871391589*fR[25]+0.7071067811865475*fR[11]; 
  fUp_R[11] = 1.5811388300841898*fR[75]-1.224744871391589*fR[27]+0.7071067811865475*fR[16]; 
  fUp_R[12] = 1.5811388300841898*fR[76]-1.224744871391589*fR[28]+0.7071067811865475*fR[17]; 
  fUp_R[13] = 1.5811388300841898*fR[77]-1.224744871391589*fR[29]+0.7071067811865475*fR[18]; 
  fUp_R[14] = 1.5811388300841898*fR[78]-1.224744871391589*fR[30]+0.7071067811865475*fR[19]; 
  fUp_R[15] = 1.5811388300841895*fR[79]-1.224744871391589*fR[31]+0.7071067811865475*fR[26]; 
   } 
  } else { 
  double fC_r[16] = {0.0}; 
  double fR_l[16] = {0.0}; 
  double sgn_drag_coeff_Up_R[16] = {0.0}; 
  ser_5x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_R, sgn_drag_coeff_Up_R); 

  fC_r[0] = 1.5811388300841895*fC[64]+1.224744871391589*fC[5]+0.7071067811865475*fC[0]; 
  fC_r[1] = 1.5811388300841898*fC[65]+1.224744871391589*fC[12]+0.7071067811865475*fC[1]; 
  fC_r[2] = 1.5811388300841898*fC[66]+1.224744871391589*fC[13]+0.7071067811865475*fC[2]; 
  fC_r[3] = 1.5811388300841898*fC[67]+1.224744871391589*fC[14]+0.7071067811865475*fC[3]; 
  fC_r[4] = 1.5811388300841898*fC[68]+1.224744871391589*fC[15]+0.7071067811865475*fC[4]; 
  fC_r[5] = 1.5811388300841895*fC[69]+1.224744871391589*fC[20]+0.7071067811865475*fC[6]; 
  fC_r[6] = 1.5811388300841895*fC[70]+1.224744871391589*fC[21]+0.7071067811865475*fC[7]; 
  fC_r[7] = 1.5811388300841895*fC[71]+1.224744871391589*fC[22]+0.7071067811865475*fC[8]; 
  fC_r[8] = 1.5811388300841895*fC[72]+1.224744871391589*fC[23]+0.7071067811865475*fC[9]; 
  fC_r[9] = 1.5811388300841895*fC[73]+1.224744871391589*fC[24]+0.7071067811865475*fC[10]; 
  fC_r[10] = 1.5811388300841895*fC[74]+1.224744871391589*fC[25]+0.7071067811865475*fC[11]; 
  fC_r[11] = 1.5811388300841898*fC[75]+1.224744871391589*fC[27]+0.7071067811865475*fC[16]; 
  fC_r[12] = 1.5811388300841898*fC[76]+1.224744871391589*fC[28]+0.7071067811865475*fC[17]; 
  fC_r[13] = 1.5811388300841898*fC[77]+1.224744871391589*fC[29]+0.7071067811865475*fC[18]; 
  fC_r[14] = 1.5811388300841898*fC[78]+1.224744871391589*fC[30]+0.7071067811865475*fC[19]; 
  fC_r[15] = 1.5811388300841895*fC[79]+1.224744871391589*fC[31]+0.7071067811865475*fC[26]; 

  fR_l[0] = 1.5811388300841895*fR[64]-1.224744871391589*fR[5]+0.7071067811865475*fR[0]; 
  fR_l[1] = 1.5811388300841898*fR[65]-1.224744871391589*fR[12]+0.7071067811865475*fR[1]; 
  fR_l[2] = 1.5811388300841898*fR[66]-1.224744871391589*fR[13]+0.7071067811865475*fR[2]; 
  fR_l[3] = 1.5811388300841898*fR[67]-1.224744871391589*fR[14]+0.7071067811865475*fR[3]; 
  fR_l[4] = 1.5811388300841898*fR[68]-1.224744871391589*fR[15]+0.7071067811865475*fR[4]; 
  fR_l[5] = 1.5811388300841895*fR[69]-1.224744871391589*fR[20]+0.7071067811865475*fR[6]; 
  fR_l[6] = 1.5811388300841895*fR[70]-1.224744871391589*fR[21]+0.7071067811865475*fR[7]; 
  fR_l[7] = 1.5811388300841895*fR[71]-1.224744871391589*fR[22]+0.7071067811865475*fR[8]; 
  fR_l[8] = 1.5811388300841895*fR[72]-1.224744871391589*fR[23]+0.7071067811865475*fR[9]; 
  fR_l[9] = 1.5811388300841895*fR[73]-1.224744871391589*fR[24]+0.7071067811865475*fR[10]; 
  fR_l[10] = 1.5811388300841895*fR[74]-1.224744871391589*fR[25]+0.7071067811865475*fR[11]; 
  fR_l[11] = 1.5811388300841898*fR[75]-1.224744871391589*fR[27]+0.7071067811865475*fR[16]; 
  fR_l[12] = 1.5811388300841898*fR[76]-1.224744871391589*fR[28]+0.7071067811865475*fR[17]; 
  fR_l[13] = 1.5811388300841898*fR[77]-1.224744871391589*fR[29]+0.7071067811865475*fR[18]; 
  fR_l[14] = 1.5811388300841898*fR[78]-1.224744871391589*fR[30]+0.7071067811865475*fR[19]; 
  fR_l[15] = 1.5811388300841895*fR[79]-1.224744871391589*fR[31]+0.7071067811865475*fR[26]; 

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
  out[4] += (0.35355339059327373*GhatL[4]-0.35355339059327373*GhatR[4])*dv_inv; 
  out[5] += -(0.6123724356957945*(GhatR[0]+GhatL[0])*dv_inv); 
  out[6] += (0.35355339059327373*GhatL[5]-0.35355339059327373*GhatR[5])*dv_inv; 
  out[7] += (0.35355339059327373*GhatL[6]-0.35355339059327373*GhatR[6])*dv_inv; 
  out[8] += (0.35355339059327373*GhatL[7]-0.35355339059327373*GhatR[7])*dv_inv; 
  out[9] += (0.35355339059327373*GhatL[8]-0.35355339059327373*GhatR[8])*dv_inv; 
  out[10] += (0.35355339059327373*GhatL[9]-0.35355339059327373*GhatR[9])*dv_inv; 
  out[11] += (0.35355339059327373*GhatL[10]-0.35355339059327373*GhatR[10])*dv_inv; 
  out[12] += -(0.6123724356957945*(GhatR[1]+GhatL[1])*dv_inv); 
  out[13] += -(0.6123724356957945*(GhatR[2]+GhatL[2])*dv_inv); 
  out[14] += -(0.6123724356957945*(GhatR[3]+GhatL[3])*dv_inv); 
  out[15] += -(0.6123724356957945*(GhatR[4]+GhatL[4])*dv_inv); 
  out[16] += (0.35355339059327373*GhatL[11]-0.35355339059327373*GhatR[11])*dv_inv; 
  out[17] += (0.35355339059327373*GhatL[12]-0.35355339059327373*GhatR[12])*dv_inv; 
  out[18] += (0.35355339059327373*GhatL[13]-0.35355339059327373*GhatR[13])*dv_inv; 
  out[19] += (0.35355339059327373*GhatL[14]-0.35355339059327373*GhatR[14])*dv_inv; 
  out[20] += -(0.6123724356957945*(GhatR[5]+GhatL[5])*dv_inv); 
  out[21] += -(0.6123724356957945*(GhatR[6]+GhatL[6])*dv_inv); 
  out[22] += -(0.6123724356957945*(GhatR[7]+GhatL[7])*dv_inv); 
  out[23] += -(0.6123724356957945*(GhatR[8]+GhatL[8])*dv_inv); 
  out[24] += -(0.6123724356957945*(GhatR[9]+GhatL[9])*dv_inv); 
  out[25] += -(0.6123724356957945*(GhatR[10]+GhatL[10])*dv_inv); 
  out[26] += (0.35355339059327373*GhatL[15]-0.35355339059327373*GhatR[15])*dv_inv; 
  out[27] += -(0.6123724356957945*(GhatR[11]+GhatL[11])*dv_inv); 
  out[28] += -(0.6123724356957945*(GhatR[12]+GhatL[12])*dv_inv); 
  out[29] += -(0.6123724356957945*(GhatR[13]+GhatL[13])*dv_inv); 
  out[30] += -(0.6123724356957945*(GhatR[14]+GhatL[14])*dv_inv); 
  out[31] += -(0.6123724356957945*(GhatR[15]+GhatL[15])*dv_inv); 
  out[64] += (0.7905694150420948*GhatL[0]-0.7905694150420948*GhatR[0])*dv_inv; 
  out[65] += (0.7905694150420949*GhatL[1]-0.7905694150420949*GhatR[1])*dv_inv; 
  out[66] += (0.7905694150420949*GhatL[2]-0.7905694150420949*GhatR[2])*dv_inv; 
  out[67] += (0.7905694150420949*GhatL[3]-0.7905694150420949*GhatR[3])*dv_inv; 
  out[68] += (0.7905694150420949*GhatL[4]-0.7905694150420949*GhatR[4])*dv_inv; 
  out[69] += (0.7905694150420948*GhatL[5]-0.7905694150420948*GhatR[5])*dv_inv; 
  out[70] += (0.7905694150420948*GhatL[6]-0.7905694150420948*GhatR[6])*dv_inv; 
  out[71] += (0.7905694150420948*GhatL[7]-0.7905694150420948*GhatR[7])*dv_inv; 
  out[72] += (0.7905694150420948*GhatL[8]-0.7905694150420948*GhatR[8])*dv_inv; 
  out[73] += (0.7905694150420948*GhatL[9]-0.7905694150420948*GhatR[9])*dv_inv; 
  out[74] += (0.7905694150420948*GhatL[10]-0.7905694150420948*GhatR[10])*dv_inv; 
  out[75] += (0.7905694150420949*GhatL[11]-0.7905694150420949*GhatR[11])*dv_inv; 
  out[76] += (0.7905694150420949*GhatL[12]-0.7905694150420949*GhatR[12])*dv_inv; 
  out[77] += (0.7905694150420949*GhatL[13]-0.7905694150420949*GhatR[13])*dv_inv; 
  out[78] += (0.7905694150420949*GhatL[14]-0.7905694150420949*GhatR[14])*dv_inv; 
  out[79] += (0.7905694150420948*GhatL[15]-0.7905694150420948*GhatR[15])*dv_inv; 
  double cflFreq = fmax(fabs(drag_coeff_surf_L[0]), fabs(drag_coeff_surf_R[0])); 

  return 0.375*dv_inv*cflFreq; 
}
