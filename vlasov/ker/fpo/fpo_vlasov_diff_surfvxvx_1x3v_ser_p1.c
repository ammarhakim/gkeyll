#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[1]/dxv[1]; 

  double f_rec_lo[8] = {0.0}; 
  double f_rec_up[8] = {0.0}; 
  double df_rec_lo[8] = {0.0}; 
  double df_rec_up[8] = {0.0}; 
  double surft1_lo[8] = {0.0}; 
  double surft1_up[8] = {0.0}; 
  double surft2_lo[8] = {0.0}; 
  double surft2_up[8] = {0.0}; 
  double vol[40] = {0.0}; 

  // Index into D and f stencils. 
  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 
  const double* fR = f_stencil[2]; 

  const double* DC = &diff_coeff_C[0]; 
  const double* D_surf_C = &diff_coeff_surf_stencil[1][0]; 
  const double* D_surf_R = &diff_coeff_surf_stencil[2][0]; 

  f_rec_lo[0] = 0.408248290463863*fL[2]-0.408248290463863*fC[2]+0.3535533905932737*(fL[0]+fC[0]); 
  f_rec_lo[1] = 0.408248290463863*fL[5]-0.408248290463863*fC[5]+0.3535533905932737*(fL[1]+fC[1]); 
  f_rec_lo[2] = 0.408248290463863*fL[7]-0.408248290463863*fC[7]+0.3535533905932737*(fL[3]+fC[3]); 
  f_rec_lo[3] = 0.408248290463863*fL[9]-0.408248290463863*fC[9]+0.3535533905932737*(fL[4]+fC[4]); 
  f_rec_lo[4] = 0.408248290463863*fL[11]-0.408248290463863*fC[11]+0.3535533905932737*(fL[6]+fC[6]); 
  f_rec_lo[5] = 0.408248290463863*fL[12]-0.408248290463863*fC[12]+0.3535533905932737*(fL[8]+fC[8]); 
  f_rec_lo[6] = 0.408248290463863*fL[14]-0.408248290463863*fC[14]+0.3535533905932737*(fL[10]+fC[10]); 
  f_rec_lo[7] = 0.408248290463863*fL[15]-0.408248290463863*fC[15]+0.3535533905932737*(fL[13]+fC[13]); 
  f_rec_up[0] = -(0.408248290463863*fR[2])+0.408248290463863*fC[2]+0.3535533905932737*(fR[0]+fC[0]); 
  f_rec_up[1] = -(0.408248290463863*fR[5])+0.408248290463863*fC[5]+0.3535533905932737*(fR[1]+fC[1]); 
  f_rec_up[2] = -(0.408248290463863*fR[7])+0.408248290463863*fC[7]+0.3535533905932737*(fR[3]+fC[3]); 
  f_rec_up[3] = -(0.408248290463863*fR[9])+0.408248290463863*fC[9]+0.3535533905932737*(fR[4]+fC[4]); 
  f_rec_up[4] = -(0.408248290463863*fR[11])+0.408248290463863*fC[11]+0.3535533905932737*(fR[6]+fC[6]); 
  f_rec_up[5] = -(0.408248290463863*fR[12])+0.408248290463863*fC[12]+0.3535533905932737*(fR[8]+fC[8]); 
  f_rec_up[6] = -(0.408248290463863*fR[14])+0.408248290463863*fC[14]+0.3535533905932737*(fR[10]+fC[10]); 
  f_rec_up[7] = -(0.408248290463863*fR[15])+0.408248290463863*fC[15]+0.3535533905932737*(fR[13]+fC[13]); 

  df_rec_lo[0] = -(0.7654655446197428*(fL[2]+fC[2]))-0.7954951288348656*fL[0]+0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*(fL[5]+fC[5]))-0.7954951288348656*fL[1]+0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*(fL[7]+fC[7]))-0.7954951288348656*fL[3]+0.7954951288348656*fC[3]; 
  df_rec_lo[3] = -(0.7654655446197428*(fL[9]+fC[9]))-0.7954951288348656*fL[4]+0.7954951288348656*fC[4]; 
  df_rec_lo[4] = -(0.7654655446197428*(fL[11]+fC[11]))-0.7954951288348656*fL[6]+0.7954951288348656*fC[6]; 
  df_rec_lo[5] = -(0.7654655446197428*(fL[12]+fC[12]))-0.7954951288348656*fL[8]+0.7954951288348656*fC[8]; 
  df_rec_lo[6] = -(0.7654655446197428*(fL[14]+fC[14]))-0.7954951288348656*fL[10]+0.7954951288348656*fC[10]; 
  df_rec_lo[7] = -(0.7654655446197428*(fL[15]+fC[15]))-0.7954951288348656*fL[13]+0.7954951288348656*fC[13]; 
  df_rec_up[0] = -(0.7654655446197428*(fR[2]+fC[2]))+0.7954951288348656*fR[0]-0.7954951288348656*fC[0]; 
  df_rec_up[1] = -(0.7654655446197428*(fR[5]+fC[5]))+0.7954951288348656*fR[1]-0.7954951288348656*fC[1]; 
  df_rec_up[2] = -(0.7654655446197428*(fR[7]+fC[7]))+0.7954951288348656*fR[3]-0.7954951288348656*fC[3]; 
  df_rec_up[3] = -(0.7654655446197428*(fR[9]+fC[9]))+0.7954951288348656*fR[4]-0.7954951288348656*fC[4]; 
  df_rec_up[4] = -(0.7654655446197428*(fR[11]+fC[11]))+0.7954951288348656*fR[6]-0.7954951288348656*fC[6]; 
  df_rec_up[5] = -(0.7654655446197428*(fR[12]+fC[12]))+0.7954951288348656*fR[8]-0.7954951288348656*fC[8]; 
  df_rec_up[6] = -(0.7654655446197428*(fR[14]+fC[14]))+0.7954951288348656*fR[10]-0.7954951288348656*fC[10]; 
  df_rec_up[7] = -(0.7654655446197428*(fR[15]+fC[15]))+0.7954951288348656*fR[13]-0.7954951288348656*fC[13]; 

  surft1_lo[0] = 0.3535533905932737*(D_surf_C[7]*df_rec_lo[7]+D_surf_C[6]*df_rec_lo[6]+D_surf_C[5]*df_rec_lo[5]+D_surf_C[4]*df_rec_lo[4]+D_surf_C[3]*df_rec_lo[3]+D_surf_C[2]*df_rec_lo[2]+D_surf_C[1]*df_rec_lo[1]+D_surf_C[0]*df_rec_lo[0]); 
  surft1_lo[1] = 0.3535533905932737*(D_surf_C[6]*df_rec_lo[7]+df_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[5]+df_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[4]+df_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[1]+df_rec_lo[0]*D_surf_C[1]); 
  surft1_lo[2] = 0.3535533905932737*(D_surf_C[5]*df_rec_lo[7]+df_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*df_rec_lo[6]+df_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[4]+df_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*df_rec_lo[2]+df_rec_lo[0]*D_surf_C[2]); 
  surft1_lo[3] = 0.3535533905932737*(D_surf_C[4]*df_rec_lo[7]+df_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*df_rec_lo[6]+df_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*df_rec_lo[5]+df_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*df_rec_lo[3]+df_rec_lo[0]*D_surf_C[3]); 
  surft1_lo[4] = 0.3535533905932737*(D_surf_C[3]*df_rec_lo[7]+df_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*df_rec_lo[6]+df_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[4]+df_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*df_rec_lo[2]+df_rec_lo[1]*D_surf_C[2]); 
  surft1_lo[5] = 0.3535533905932737*(D_surf_C[2]*df_rec_lo[7]+df_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*df_rec_lo[6]+df_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*df_rec_lo[5]+df_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*df_rec_lo[3]+df_rec_lo[1]*D_surf_C[3]); 
  surft1_lo[6] = 0.3535533905932737*(D_surf_C[1]*df_rec_lo[7]+df_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*df_rec_lo[6]+df_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*df_rec_lo[5]+df_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*df_rec_lo[3]+df_rec_lo[2]*D_surf_C[3]); 
  surft1_lo[7] = 0.3535533905932737*(D_surf_C[0]*df_rec_lo[7]+df_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*df_rec_lo[6]+df_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*df_rec_lo[5]+df_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*df_rec_lo[4]+df_rec_lo[3]*D_surf_C[4]); 
  surft1_up[0] = 0.3535533905932737*(D_surf_R[7]*df_rec_up[7]+D_surf_R[6]*df_rec_up[6]+D_surf_R[5]*df_rec_up[5]+D_surf_R[4]*df_rec_up[4]+D_surf_R[3]*df_rec_up[3]+D_surf_R[2]*df_rec_up[2]+D_surf_R[1]*df_rec_up[1]+D_surf_R[0]*df_rec_up[0]); 
  surft1_up[1] = 0.3535533905932737*(D_surf_R[6]*df_rec_up[7]+df_rec_up[6]*D_surf_R[7]+D_surf_R[3]*df_rec_up[5]+df_rec_up[3]*D_surf_R[5]+D_surf_R[2]*df_rec_up[4]+df_rec_up[2]*D_surf_R[4]+D_surf_R[0]*df_rec_up[1]+df_rec_up[0]*D_surf_R[1]); 
  surft1_up[2] = 0.3535533905932737*(D_surf_R[5]*df_rec_up[7]+df_rec_up[5]*D_surf_R[7]+D_surf_R[3]*df_rec_up[6]+df_rec_up[3]*D_surf_R[6]+D_surf_R[1]*df_rec_up[4]+df_rec_up[1]*D_surf_R[4]+D_surf_R[0]*df_rec_up[2]+df_rec_up[0]*D_surf_R[2]); 
  surft1_up[3] = 0.3535533905932737*(D_surf_R[4]*df_rec_up[7]+df_rec_up[4]*D_surf_R[7]+D_surf_R[2]*df_rec_up[6]+df_rec_up[2]*D_surf_R[6]+D_surf_R[1]*df_rec_up[5]+df_rec_up[1]*D_surf_R[5]+D_surf_R[0]*df_rec_up[3]+df_rec_up[0]*D_surf_R[3]); 
  surft1_up[4] = 0.3535533905932737*(D_surf_R[3]*df_rec_up[7]+df_rec_up[3]*D_surf_R[7]+D_surf_R[5]*df_rec_up[6]+df_rec_up[5]*D_surf_R[6]+D_surf_R[0]*df_rec_up[4]+df_rec_up[0]*D_surf_R[4]+D_surf_R[1]*df_rec_up[2]+df_rec_up[1]*D_surf_R[2]); 
  surft1_up[5] = 0.3535533905932737*(D_surf_R[2]*df_rec_up[7]+df_rec_up[2]*D_surf_R[7]+D_surf_R[4]*df_rec_up[6]+df_rec_up[4]*D_surf_R[6]+D_surf_R[0]*df_rec_up[5]+df_rec_up[0]*D_surf_R[5]+D_surf_R[1]*df_rec_up[3]+df_rec_up[1]*D_surf_R[3]); 
  surft1_up[6] = 0.3535533905932737*(D_surf_R[1]*df_rec_up[7]+df_rec_up[1]*D_surf_R[7]+D_surf_R[0]*df_rec_up[6]+df_rec_up[0]*D_surf_R[6]+D_surf_R[4]*df_rec_up[5]+df_rec_up[4]*D_surf_R[5]+D_surf_R[2]*df_rec_up[3]+df_rec_up[2]*D_surf_R[3]); 
  surft1_up[7] = 0.3535533905932737*(D_surf_R[0]*df_rec_up[7]+df_rec_up[0]*D_surf_R[7]+D_surf_R[1]*df_rec_up[6]+df_rec_up[1]*D_surf_R[6]+D_surf_R[2]*df_rec_up[5]+df_rec_up[2]*D_surf_R[5]+D_surf_R[3]*df_rec_up[4]+df_rec_up[3]*D_surf_R[4]); 

  surft2_lo[0] = 0.3535533905932737*(D_surf_C[7]*f_rec_lo[7]+D_surf_C[6]*f_rec_lo[6]+D_surf_C[5]*f_rec_lo[5]+D_surf_C[4]*f_rec_lo[4]+D_surf_C[3]*f_rec_lo[3]+D_surf_C[2]*f_rec_lo[2]+D_surf_C[1]*f_rec_lo[1]+D_surf_C[0]*f_rec_lo[0]); 
  surft2_lo[1] = 0.3535533905932737*(D_surf_C[6]*f_rec_lo[7]+f_rec_lo[6]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[5]+f_rec_lo[3]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[4]+f_rec_lo[2]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[1]+f_rec_lo[0]*D_surf_C[1]); 
  surft2_lo[2] = 0.3535533905932737*(D_surf_C[5]*f_rec_lo[7]+f_rec_lo[5]*D_surf_C[7]+D_surf_C[3]*f_rec_lo[6]+f_rec_lo[3]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[4]+f_rec_lo[1]*D_surf_C[4]+D_surf_C[0]*f_rec_lo[2]+f_rec_lo[0]*D_surf_C[2]); 
  surft2_lo[3] = 0.3535533905932737*(D_surf_C[4]*f_rec_lo[7]+f_rec_lo[4]*D_surf_C[7]+D_surf_C[2]*f_rec_lo[6]+f_rec_lo[2]*D_surf_C[6]+D_surf_C[1]*f_rec_lo[5]+f_rec_lo[1]*D_surf_C[5]+D_surf_C[0]*f_rec_lo[3]+f_rec_lo[0]*D_surf_C[3]); 
  surft2_lo[4] = 0.3535533905932737*(D_surf_C[3]*f_rec_lo[7]+f_rec_lo[3]*D_surf_C[7]+D_surf_C[5]*f_rec_lo[6]+f_rec_lo[5]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[4]+f_rec_lo[0]*D_surf_C[4]+D_surf_C[1]*f_rec_lo[2]+f_rec_lo[1]*D_surf_C[2]); 
  surft2_lo[5] = 0.3535533905932737*(D_surf_C[2]*f_rec_lo[7]+f_rec_lo[2]*D_surf_C[7]+D_surf_C[4]*f_rec_lo[6]+f_rec_lo[4]*D_surf_C[6]+D_surf_C[0]*f_rec_lo[5]+f_rec_lo[0]*D_surf_C[5]+D_surf_C[1]*f_rec_lo[3]+f_rec_lo[1]*D_surf_C[3]); 
  surft2_lo[6] = 0.3535533905932737*(D_surf_C[1]*f_rec_lo[7]+f_rec_lo[1]*D_surf_C[7]+D_surf_C[0]*f_rec_lo[6]+f_rec_lo[0]*D_surf_C[6]+D_surf_C[4]*f_rec_lo[5]+f_rec_lo[4]*D_surf_C[5]+D_surf_C[2]*f_rec_lo[3]+f_rec_lo[2]*D_surf_C[3]); 
  surft2_lo[7] = 0.3535533905932737*(D_surf_C[0]*f_rec_lo[7]+f_rec_lo[0]*D_surf_C[7]+D_surf_C[1]*f_rec_lo[6]+f_rec_lo[1]*D_surf_C[6]+D_surf_C[2]*f_rec_lo[5]+f_rec_lo[2]*D_surf_C[5]+D_surf_C[3]*f_rec_lo[4]+f_rec_lo[3]*D_surf_C[4]); 
  surft2_up[0] = 0.3535533905932737*(D_surf_R[7]*f_rec_up[7]+D_surf_R[6]*f_rec_up[6]+D_surf_R[5]*f_rec_up[5]+D_surf_R[4]*f_rec_up[4]+D_surf_R[3]*f_rec_up[3]+D_surf_R[2]*f_rec_up[2]+D_surf_R[1]*f_rec_up[1]+D_surf_R[0]*f_rec_up[0]); 
  surft2_up[1] = 0.3535533905932737*(D_surf_R[6]*f_rec_up[7]+f_rec_up[6]*D_surf_R[7]+D_surf_R[3]*f_rec_up[5]+f_rec_up[3]*D_surf_R[5]+D_surf_R[2]*f_rec_up[4]+f_rec_up[2]*D_surf_R[4]+D_surf_R[0]*f_rec_up[1]+f_rec_up[0]*D_surf_R[1]); 
  surft2_up[2] = 0.3535533905932737*(D_surf_R[5]*f_rec_up[7]+f_rec_up[5]*D_surf_R[7]+D_surf_R[3]*f_rec_up[6]+f_rec_up[3]*D_surf_R[6]+D_surf_R[1]*f_rec_up[4]+f_rec_up[1]*D_surf_R[4]+D_surf_R[0]*f_rec_up[2]+f_rec_up[0]*D_surf_R[2]); 
  surft2_up[3] = 0.3535533905932737*(D_surf_R[4]*f_rec_up[7]+f_rec_up[4]*D_surf_R[7]+D_surf_R[2]*f_rec_up[6]+f_rec_up[2]*D_surf_R[6]+D_surf_R[1]*f_rec_up[5]+f_rec_up[1]*D_surf_R[5]+D_surf_R[0]*f_rec_up[3]+f_rec_up[0]*D_surf_R[3]); 
  surft2_up[4] = 0.3535533905932737*(D_surf_R[3]*f_rec_up[7]+f_rec_up[3]*D_surf_R[7]+D_surf_R[5]*f_rec_up[6]+f_rec_up[5]*D_surf_R[6]+D_surf_R[0]*f_rec_up[4]+f_rec_up[0]*D_surf_R[4]+D_surf_R[1]*f_rec_up[2]+f_rec_up[1]*D_surf_R[2]); 
  surft2_up[5] = 0.3535533905932737*(D_surf_R[2]*f_rec_up[7]+f_rec_up[2]*D_surf_R[7]+D_surf_R[4]*f_rec_up[6]+f_rec_up[4]*D_surf_R[6]+D_surf_R[0]*f_rec_up[5]+f_rec_up[0]*D_surf_R[5]+D_surf_R[1]*f_rec_up[3]+f_rec_up[1]*D_surf_R[3]); 
  surft2_up[6] = 0.3535533905932737*(D_surf_R[1]*f_rec_up[7]+f_rec_up[1]*D_surf_R[7]+D_surf_R[0]*f_rec_up[6]+f_rec_up[0]*D_surf_R[6]+D_surf_R[4]*f_rec_up[5]+f_rec_up[4]*D_surf_R[5]+D_surf_R[2]*f_rec_up[3]+f_rec_up[2]*D_surf_R[3]); 
  surft2_up[7] = 0.3535533905932737*(D_surf_R[0]*f_rec_up[7]+f_rec_up[0]*D_surf_R[7]+D_surf_R[1]*f_rec_up[6]+f_rec_up[1]*D_surf_R[6]+D_surf_R[2]*f_rec_up[5]+f_rec_up[2]*D_surf_R[5]+D_surf_R[3]*f_rec_up[4]+f_rec_up[3]*D_surf_R[4]); 

  vol[2] = 0.75*(fC[37]*DC[39]+fC[35]*DC[38]+fC[33]*DC[36]+fC[32]*DC[34]+fC[29]*DC[31]+fC[27]*DC[30]+fC[25]*DC[28]+fC[24]*DC[26]+fC[13]*DC[15]+fC[10]*DC[14]+fC[8]*DC[12]+fC[6]*DC[11]+fC[4]*DC[9]+fC[3]*DC[7]+fC[1]*DC[5]+fC[0]*DC[2])+1.6770509831248421*(fC[15]*DC[23]+fC[9]*DC[19]+fC[7]*DC[18]+fC[5]*DC[17])+1.6770509831248424*(fC[14]*DC[22]+fC[12]*DC[21]+fC[11]*DC[20]+fC[2]*DC[16]); 
  vol[5] = 0.75*(fC[35]*DC[39]+fC[37]*DC[38]+fC[32]*DC[36]+fC[33]*DC[34]+fC[27]*DC[31]+fC[29]*DC[30]+fC[24]*DC[28]+fC[25]*DC[26]+fC[10]*DC[15]+fC[13]*DC[14]+fC[4]*DC[12]+fC[3]*DC[11]+fC[8]*DC[9]+fC[6]*DC[7]+fC[0]*DC[5]+fC[1]*DC[2])+1.6770509831248421*(fC[14]*DC[23]+fC[12]*DC[19]+fC[11]*DC[18]+fC[2]*DC[17])+1.6770509831248424*(fC[15]*DC[22]+fC[9]*DC[21]+fC[7]*DC[20]+fC[5]*DC[16]); 
  vol[7] = 0.75*(fC[33]*DC[39]+fC[32]*DC[38]+DC[36]*fC[37]+DC[34]*fC[35]+fC[8]*DC[15]+fC[4]*DC[14]+DC[12]*fC[13]+fC[1]*DC[11]+DC[9]*fC[10]+fC[0]*DC[7]+DC[5]*fC[6]+DC[2]*fC[3])+1.5*(DC[23]*fC[31]+DC[22]*fC[30]+DC[20]*fC[28]+DC[18]*fC[26])+0.6708203932499369*(fC[13]*DC[31]+fC[10]*DC[30]+DC[15]*fC[29]+fC[6]*DC[28]+DC[14]*fC[27]+fC[3]*DC[26]+DC[11]*fC[25]+DC[7]*fC[24])+1.6770509831248421*(fC[12]*DC[23]+fC[14]*DC[19]+fC[2]*DC[18]+fC[11]*DC[17])+1.6770509831248424*(fC[9]*DC[22]+fC[15]*DC[21]+fC[5]*DC[20]+fC[7]*DC[16]); 
  vol[9] = 1.5*(DC[23]*fC[39]+DC[22]*fC[38]+DC[21]*fC[36]+DC[19]*fC[34])+0.6708203932499369*(fC[13]*DC[39]+fC[10]*DC[38]+DC[15]*fC[37]+fC[8]*DC[36]+DC[14]*fC[35]+fC[4]*DC[34]+DC[12]*fC[33]+DC[9]*fC[32])+0.75*(fC[25]*DC[31]+fC[24]*DC[30]+DC[28]*fC[29]+DC[26]*fC[27]+fC[6]*DC[15]+fC[3]*DC[14]+DC[11]*fC[13]+fC[1]*DC[12]+DC[7]*fC[10]+fC[0]*DC[9]+DC[5]*fC[8]+DC[2]*fC[4])+1.6770509831248421*(fC[11]*DC[23]+fC[2]*DC[19]+fC[14]*DC[18]+fC[12]*DC[17])+1.6770509831248424*(fC[7]*DC[22]+fC[5]*DC[21]+fC[15]*DC[20]+fC[9]*DC[16]); 
  vol[11] = 0.75*(fC[32]*DC[39]+fC[33]*DC[38]+DC[34]*fC[37]+fC[35]*DC[36]+fC[4]*DC[15]+fC[8]*DC[14]+DC[9]*fC[13]+fC[10]*DC[12]+fC[0]*DC[11]+fC[1]*DC[7]+DC[2]*fC[6]+fC[3]*DC[5])+1.5*(DC[22]*fC[31]+DC[23]*fC[30]+DC[18]*fC[28]+DC[20]*fC[26])+0.6708203932499369*(fC[10]*DC[31]+fC[13]*DC[30]+DC[14]*fC[29]+fC[3]*DC[28]+DC[15]*fC[27]+fC[6]*DC[26]+DC[7]*fC[25]+DC[11]*fC[24])+1.6770509831248421*(fC[9]*DC[23]+fC[15]*DC[19]+fC[5]*DC[18]+fC[7]*DC[17])+1.6770509831248424*(fC[12]*DC[22]+fC[14]*DC[21]+fC[2]*DC[20]+fC[11]*DC[16]); 
  vol[12] = 1.5*(DC[22]*fC[39]+DC[23]*fC[38]+DC[19]*fC[36]+DC[21]*fC[34])+0.6708203932499369*(fC[10]*DC[39]+fC[13]*DC[38]+DC[14]*fC[37]+fC[4]*DC[36]+DC[15]*fC[35]+fC[8]*DC[34]+DC[9]*fC[33]+DC[12]*fC[32])+0.75*(fC[24]*DC[31]+fC[25]*DC[30]+DC[26]*fC[29]+fC[27]*DC[28]+fC[3]*DC[15]+fC[6]*DC[14]+DC[7]*fC[13]+fC[0]*DC[12]+fC[10]*DC[11]+fC[1]*DC[9]+DC[2]*fC[8]+fC[4]*DC[5])+1.6770509831248421*(fC[7]*DC[23]+fC[5]*DC[19]+fC[15]*DC[18]+fC[9]*DC[17])+1.6770509831248424*(fC[11]*DC[22]+fC[2]*DC[21]+fC[14]*DC[20]+fC[12]*DC[16]); 
  vol[14] = 1.5*(DC[21]*fC[39]+DC[19]*fC[38]+DC[23]*fC[36]+DC[22]*fC[34]+DC[20]*fC[31]+DC[18]*fC[30]+DC[23]*fC[28]+DC[22]*fC[26])+0.5999999999999999*(fC[29]*DC[39]+fC[27]*DC[38]+DC[31]*fC[37]+DC[30]*fC[35])+0.6708203932499369*(fC[8]*DC[39]+fC[4]*DC[38]+DC[12]*fC[37]+fC[13]*DC[36]+DC[9]*fC[35]+fC[10]*DC[34]+DC[15]*fC[33]+DC[14]*fC[32]+fC[6]*DC[31]+fC[3]*DC[30]+DC[11]*fC[29]+fC[13]*DC[28]+DC[7]*fC[27]+fC[10]*DC[26]+DC[15]*fC[25]+DC[14]*fC[24])+1.6770509831248421*(fC[5]*DC[23]+fC[7]*DC[19]+fC[9]*DC[18]+fC[15]*DC[17])+1.6770509831248424*(fC[2]*DC[22]+fC[11]*DC[21]+fC[12]*DC[20]+fC[14]*DC[16])+0.75*(fC[1]*DC[15]+fC[0]*DC[14]+DC[5]*fC[13]+fC[6]*DC[12]+fC[8]*DC[11]+DC[2]*fC[10]+fC[3]*DC[9]+fC[4]*DC[7]); 
  vol[15] = 1.5*(DC[19]*fC[39]+DC[21]*fC[38]+DC[22]*fC[36]+DC[23]*fC[34]+DC[18]*fC[31]+DC[20]*fC[30]+DC[22]*fC[28]+DC[23]*fC[26])+0.6*(fC[27]*DC[39]+fC[29]*DC[38]+DC[30]*fC[37]+DC[31]*fC[35])+0.6708203932499369*(fC[4]*DC[39]+fC[8]*DC[38]+DC[9]*fC[37]+fC[10]*DC[36]+DC[12]*fC[35]+fC[13]*DC[34]+DC[14]*fC[33]+DC[15]*fC[32]+fC[3]*DC[31]+fC[6]*DC[30]+DC[7]*fC[29]+fC[10]*DC[28]+DC[11]*fC[27]+fC[13]*DC[26]+DC[14]*fC[25]+DC[15]*fC[24])+1.6770509831248421*(fC[2]*DC[23]+fC[11]*DC[19]+fC[12]*DC[18]+fC[14]*DC[17])+1.6770509831248424*(fC[5]*DC[22]+fC[7]*DC[21]+fC[9]*DC[20]+fC[15]*DC[16])+0.75*(fC[0]*DC[15]+fC[1]*DC[14]+DC[2]*fC[13]+fC[3]*DC[12]+fC[4]*DC[11]+DC[5]*fC[10]+fC[6]*DC[9]+DC[7]*fC[8]); 
  vol[16] = 3.3541019662496847*(DC[39]*fC[39]+DC[38]*fC[38]+DC[36]*fC[36]+DC[34]*fC[34]+DC[31]*fC[31]+DC[30]*fC[30]+DC[28]*fC[28]+DC[26]*fC[26]+DC[15]*fC[15]+DC[14]*fC[14]+DC[12]*fC[12]+DC[11]*fC[11]+DC[9]*fC[9]+DC[7]*fC[7]+DC[5]*fC[5]+DC[2]*fC[2])+1.6770509831248424*(DC[37]*fC[37]+DC[35]*fC[35]+DC[33]*fC[33]+DC[32]*fC[32]+DC[29]*fC[29]+DC[27]*fC[27]+DC[25]*fC[25]+DC[24]*fC[24]+DC[13]*fC[13]+DC[10]*fC[10]+DC[8]*fC[8]+DC[6]*fC[6]+DC[4]*fC[4]+DC[3]*fC[3]+DC[1]*fC[1]+DC[0]*fC[0])+5.031152949374527*(DC[23]*fC[23]+DC[22]*fC[22]+DC[21]*fC[21]+DC[20]*fC[20]+DC[19]*fC[19]+DC[18]*fC[18]+DC[17]*fC[17]+DC[16]*fC[16])+3.7500000000000004*(fC[13]*DC[23]+fC[4]*DC[19]+fC[3]*DC[18]+fC[1]*DC[17])+3.75*(fC[10]*DC[22]+fC[8]*DC[21]+fC[6]*DC[20]+fC[0]*DC[16]); 
  vol[17] = 3.3541019662496847*(DC[38]*fC[39]+fC[38]*DC[39]+DC[34]*fC[36]+fC[34]*DC[36]+DC[30]*fC[31]+fC[30]*DC[31]+DC[26]*fC[28]+fC[26]*DC[28])+1.6770509831248424*(DC[35]*fC[37]+fC[35]*DC[37]+DC[32]*fC[33]+fC[32]*DC[33]+DC[27]*fC[29]+fC[27]*DC[29]+DC[24]*fC[25]+fC[24]*DC[25])+5.031152949374527*(DC[22]*fC[23]+fC[22]*DC[23]+DC[19]*fC[21]+fC[19]*DC[21]+DC[18]*fC[20]+fC[18]*DC[20]+DC[16]*fC[17]+fC[16]*DC[17])+3.75*(fC[10]*DC[23]+fC[8]*DC[19]+fC[6]*DC[18]+fC[0]*DC[17])+3.7500000000000004*(fC[13]*DC[22]+fC[4]*DC[21]+fC[3]*DC[20]+fC[1]*DC[16])+3.3541019662496843*(DC[14]*fC[15]+fC[14]*DC[15]+DC[9]*fC[12]+fC[9]*DC[12]+DC[7]*fC[11]+fC[7]*DC[11]+DC[2]*fC[5]+fC[2]*DC[5])+1.6770509831248421*(DC[10]*fC[13]+fC[10]*DC[13]+DC[4]*fC[8]+fC[4]*DC[8]+DC[3]*fC[6]+fC[3]*DC[6]+DC[0]*fC[1]+fC[0]*DC[1]); 
  vol[18] = 3.3541019662496847*(DC[36]*fC[39]+fC[36]*DC[39]+DC[34]*fC[38]+fC[34]*DC[38]+DC[23]*fC[29]+DC[22]*fC[27]+DC[20]*fC[25]+DC[18]*fC[24])+1.6770509831248424*(DC[33]*fC[37]+fC[33]*DC[37]+DC[32]*fC[35]+fC[32]*DC[35])+3.0*(DC[15]*fC[31]+fC[15]*DC[31]+DC[14]*fC[30]+fC[14]*DC[30]+DC[11]*fC[28]+fC[11]*DC[28]+DC[7]*fC[26]+fC[7]*DC[26])+1.5*(DC[13]*fC[29]+fC[13]*DC[29]+DC[10]*fC[27]+fC[10]*DC[27]+DC[6]*fC[25]+fC[6]*DC[25]+DC[3]*fC[24]+fC[3]*DC[24])+5.031152949374527*(DC[21]*fC[23]+fC[21]*DC[23]+DC[19]*fC[22]+fC[19]*DC[22]+DC[17]*fC[20]+fC[17]*DC[20]+DC[16]*fC[18]+fC[16]*DC[18])+3.75*(fC[8]*DC[23]+fC[10]*DC[19]+fC[0]*DC[18]+fC[6]*DC[17])+3.7500000000000004*(fC[4]*DC[22]+fC[13]*DC[21]+fC[1]*DC[20]+fC[3]*DC[16])+3.3541019662496843*(DC[12]*fC[15]+fC[12]*DC[15]+DC[9]*fC[14]+fC[9]*DC[14]+DC[5]*fC[11]+fC[5]*DC[11]+DC[2]*fC[7]+fC[2]*DC[7])+1.6770509831248421*(DC[8]*fC[13]+fC[8]*DC[13]+DC[4]*fC[10]+fC[4]*DC[10]+DC[1]*fC[6]+fC[1]*DC[6]+DC[0]*fC[3]+fC[0]*DC[3]); 
  vol[19] = 3.0*(DC[15]*fC[39]+fC[15]*DC[39]+DC[14]*fC[38]+fC[14]*DC[38]+DC[12]*fC[36]+fC[12]*DC[36]+DC[9]*fC[34]+fC[9]*DC[34])+3.3541019662496847*(DC[23]*fC[37]+DC[22]*fC[35]+DC[21]*fC[33]+DC[19]*fC[32]+DC[28]*fC[31]+fC[28]*DC[31]+DC[26]*fC[30]+fC[26]*DC[30])+1.5*(DC[13]*fC[37]+fC[13]*DC[37]+DC[10]*fC[35]+fC[10]*DC[35]+DC[8]*fC[33]+fC[8]*DC[33]+DC[4]*fC[32]+fC[4]*DC[32])+1.6770509831248424*(DC[25]*fC[29]+fC[25]*DC[29]+DC[24]*fC[27]+fC[24]*DC[27])+5.031152949374527*(DC[20]*fC[23]+fC[20]*DC[23]+DC[18]*fC[22]+fC[18]*DC[22]+DC[17]*fC[21]+fC[17]*DC[21]+DC[16]*fC[19]+fC[16]*DC[19])+3.75*(fC[6]*DC[23]+fC[0]*DC[19]+fC[10]*DC[18]+fC[8]*DC[17])+3.7500000000000004*(fC[3]*DC[22]+fC[1]*DC[21]+fC[13]*DC[20]+fC[4]*DC[16])+3.3541019662496843*(DC[11]*fC[15]+fC[11]*DC[15]+DC[7]*fC[14]+fC[7]*DC[14]+DC[5]*fC[12]+fC[5]*DC[12]+DC[2]*fC[9]+fC[2]*DC[9])+1.6770509831248421*(DC[6]*fC[13]+fC[6]*DC[13]+DC[3]*fC[10]+fC[3]*DC[10]+DC[1]*fC[8]+fC[1]*DC[8]+DC[0]*fC[4]+fC[0]*DC[4]); 
  vol[20] = 3.3541019662496847*(DC[34]*fC[39]+fC[34]*DC[39]+DC[36]*fC[38]+fC[36]*DC[38]+DC[22]*fC[29]+DC[23]*fC[27]+DC[18]*fC[25]+DC[20]*fC[24]+DC[9]*fC[15]+fC[9]*DC[15]+DC[12]*fC[14]+fC[12]*DC[14]+DC[2]*fC[11]+fC[2]*DC[11]+DC[5]*fC[7]+fC[5]*DC[7])+1.6770509831248424*(DC[32]*fC[37]+fC[32]*DC[37]+DC[33]*fC[35]+fC[33]*DC[35]+DC[4]*fC[13]+fC[4]*DC[13]+DC[8]*fC[10]+fC[8]*DC[10]+DC[0]*fC[6]+fC[0]*DC[6]+DC[1]*fC[3]+fC[1]*DC[3])+3.0*(DC[14]*fC[31]+fC[14]*DC[31]+DC[15]*fC[30]+fC[15]*DC[30]+DC[7]*fC[28]+fC[7]*DC[28]+DC[11]*fC[26]+fC[11]*DC[26])+1.5*(DC[10]*fC[29]+fC[10]*DC[29]+DC[13]*fC[27]+fC[13]*DC[27]+DC[3]*fC[25]+fC[3]*DC[25]+DC[6]*fC[24]+fC[6]*DC[24])+5.031152949374527*(DC[19]*fC[23]+fC[19]*DC[23]+DC[21]*fC[22]+fC[21]*DC[22]+DC[16]*fC[20]+fC[16]*DC[20]+DC[17]*fC[18]+fC[17]*DC[18])+3.7500000000000004*(fC[4]*DC[23]+fC[13]*DC[19]+fC[1]*DC[18]+fC[3]*DC[17])+3.75*(fC[8]*DC[22]+fC[10]*DC[21]+fC[0]*DC[20]+fC[6]*DC[16]); 
  vol[21] = 3.0*(DC[14]*fC[39]+fC[14]*DC[39]+DC[15]*fC[38]+fC[15]*DC[38]+DC[9]*fC[36]+fC[9]*DC[36]+DC[12]*fC[34]+fC[12]*DC[34])+3.3541019662496847*(DC[22]*fC[37]+DC[23]*fC[35]+DC[19]*fC[33]+DC[21]*fC[32]+DC[26]*fC[31]+fC[26]*DC[31]+DC[28]*fC[30]+fC[28]*DC[30]+DC[7]*fC[15]+fC[7]*DC[15]+DC[11]*fC[14]+fC[11]*DC[14]+DC[2]*fC[12]+fC[2]*DC[12]+DC[5]*fC[9]+fC[5]*DC[9])+1.5*(DC[10]*fC[37]+fC[10]*DC[37]+DC[13]*fC[35]+fC[13]*DC[35]+DC[4]*fC[33]+fC[4]*DC[33]+DC[8]*fC[32]+fC[8]*DC[32])+1.6770509831248424*(DC[24]*fC[29]+fC[24]*DC[29]+DC[25]*fC[27]+fC[25]*DC[27]+DC[3]*fC[13]+fC[3]*DC[13]+DC[6]*fC[10]+fC[6]*DC[10]+DC[0]*fC[8]+fC[0]*DC[8]+DC[1]*fC[4]+fC[1]*DC[4])+5.031152949374527*(DC[18]*fC[23]+fC[18]*DC[23]+DC[20]*fC[22]+fC[20]*DC[22]+DC[16]*fC[21]+fC[16]*DC[21]+DC[17]*fC[19]+fC[17]*DC[19])+3.7500000000000004*(fC[3]*DC[23]+fC[1]*DC[19]+fC[13]*DC[18]+fC[4]*DC[17])+3.75*(fC[6]*DC[22]+fC[0]*DC[21]+fC[10]*DC[20]+fC[8]*DC[16]); 
  vol[22] = 2.6832815729997477*(DC[31]*fC[39]+fC[31]*DC[39]+DC[30]*fC[38]+fC[30]*DC[38])+3.0*(DC[12]*fC[39]+fC[12]*DC[39]+DC[9]*fC[38]+fC[9]*DC[38]+DC[15]*fC[36]+fC[15]*DC[36]+DC[14]*fC[34]+fC[14]*DC[34]+DC[11]*fC[31]+fC[11]*DC[31]+DC[7]*fC[30]+fC[7]*DC[30]+DC[15]*fC[28]+fC[15]*DC[28]+DC[14]*fC[26]+fC[14]*DC[26])+1.3416407864998738*(DC[29]*fC[37]+fC[29]*DC[37]+DC[27]*fC[35]+fC[27]*DC[35])+3.3541019662496847*(DC[21]*fC[37]+DC[19]*fC[35]+DC[23]*fC[33]+DC[22]*fC[32]+DC[20]*fC[29]+DC[18]*fC[27]+DC[23]*fC[25]+DC[22]*fC[24]+DC[5]*fC[15]+fC[5]*DC[15]+DC[2]*fC[14]+fC[2]*DC[14]+DC[11]*fC[12]+fC[11]*DC[12]+DC[7]*fC[9]+fC[7]*DC[9])+1.5*(DC[8]*fC[37]+fC[8]*DC[37]+DC[4]*fC[35]+fC[4]*DC[35]+DC[13]*fC[33]+fC[13]*DC[33]+DC[10]*fC[32]+fC[10]*DC[32]+DC[6]*fC[29]+fC[6]*DC[29]+DC[3]*fC[27]+fC[3]*DC[27]+DC[13]*fC[25]+fC[13]*DC[25]+DC[10]*fC[24]+fC[10]*DC[24])+5.031152949374527*(DC[17]*fC[23]+fC[17]*DC[23]+DC[16]*fC[22]+fC[16]*DC[22]+DC[20]*fC[21]+fC[20]*DC[21]+DC[18]*fC[19]+fC[18]*DC[19])+3.7500000000000004*(fC[1]*DC[23]+fC[3]*DC[19]+fC[4]*DC[18]+fC[13]*DC[17])+3.75*(fC[0]*DC[22]+fC[6]*DC[21]+fC[8]*DC[20]+fC[10]*DC[16])+1.6770509831248424*(DC[1]*fC[13]+fC[1]*DC[13]+DC[0]*fC[10]+fC[0]*DC[10]+DC[6]*fC[8]+fC[6]*DC[8]+DC[3]*fC[4]+fC[3]*DC[4]); 
  vol[23] = 2.6832815729997477*(DC[30]*fC[39]+fC[30]*DC[39]+DC[31]*fC[38]+fC[31]*DC[38])+3.0*(DC[9]*fC[39]+fC[9]*DC[39]+DC[12]*fC[38]+fC[12]*DC[38]+DC[14]*fC[36]+fC[14]*DC[36]+DC[15]*fC[34]+fC[15]*DC[34]+DC[7]*fC[31]+fC[7]*DC[31]+DC[11]*fC[30]+fC[11]*DC[30]+DC[14]*fC[28]+fC[14]*DC[28]+DC[15]*fC[26]+fC[15]*DC[26])+1.3416407864998738*(DC[27]*fC[37]+fC[27]*DC[37]+DC[29]*fC[35]+fC[29]*DC[35])+3.3541019662496847*(DC[19]*fC[37]+DC[21]*fC[35]+DC[22]*fC[33]+DC[23]*fC[32]+DC[18]*fC[29]+DC[20]*fC[27]+DC[22]*fC[25]+DC[23]*fC[24])+1.5*(DC[4]*fC[37]+fC[4]*DC[37]+DC[8]*fC[35]+fC[8]*DC[35]+DC[10]*fC[33]+fC[10]*DC[33]+DC[13]*fC[32]+fC[13]*DC[32]+DC[3]*fC[29]+fC[3]*DC[29]+DC[6]*fC[27]+fC[6]*DC[27]+DC[10]*fC[25]+fC[10]*DC[25]+DC[13]*fC[24]+fC[13]*DC[24])+5.031152949374527*(DC[16]*fC[23]+fC[16]*DC[23]+DC[17]*fC[22]+fC[17]*DC[22]+DC[18]*fC[21]+fC[18]*DC[21]+DC[19]*fC[20]+fC[19]*DC[20])+3.75*(fC[0]*DC[23]+fC[6]*DC[19]+fC[8]*DC[18]+fC[10]*DC[17])+3.7500000000000004*(fC[1]*DC[22]+fC[3]*DC[21]+fC[4]*DC[20]+fC[13]*DC[16])+3.3541019662496843*(DC[2]*fC[15]+fC[2]*DC[15]+DC[5]*fC[14]+fC[5]*DC[14]+DC[7]*fC[12]+fC[7]*DC[12]+DC[9]*fC[11]+fC[9]*DC[11])+1.6770509831248421*(DC[0]*fC[13]+fC[0]*DC[13]+DC[1]*fC[10]+fC[1]*DC[10]+DC[3]*fC[8]+fC[3]*DC[8]+DC[4]*fC[6]+fC[4]*DC[6]); 
  vol[26] = 0.6708203932499369*(fC[37]*DC[39]+fC[35]*DC[38]+fC[13]*DC[15]+fC[10]*DC[14]+fC[6]*DC[11]+fC[3]*DC[7])+1.6770509831248424*(DC[21]*fC[31]+DC[19]*fC[30]+DC[17]*fC[28]+DC[16]*fC[26])+0.47915742374995496*(fC[29]*DC[31]+fC[27]*DC[30]+fC[25]*DC[28]+fC[24]*DC[26])+0.75*(fC[8]*DC[31]+fC[4]*DC[30]+DC[12]*fC[29]+fC[1]*DC[28]+DC[9]*fC[27]+fC[0]*DC[26]+DC[5]*fC[25]+DC[2]*fC[24])+1.5*(fC[15]*DC[23]+fC[14]*DC[22]+fC[11]*DC[20]+fC[7]*DC[18]); 
  vol[28] = 0.6708203932499369*(fC[35]*DC[39]+fC[37]*DC[38]+fC[10]*DC[15]+fC[13]*DC[14]+fC[3]*DC[11]+fC[6]*DC[7])+1.6770509831248424*(DC[19]*fC[31]+DC[21]*fC[30]+DC[16]*fC[28]+DC[17]*fC[26])+0.47915742374995496*(fC[27]*DC[31]+fC[29]*DC[30]+fC[24]*DC[28]+fC[25]*DC[26])+0.75*(fC[4]*DC[31]+fC[8]*DC[30]+DC[9]*fC[29]+fC[0]*DC[28]+DC[12]*fC[27]+fC[1]*DC[26]+DC[2]*fC[25]+DC[5]*fC[24])+1.5*(fC[14]*DC[23]+fC[15]*DC[22]+fC[7]*DC[20]+fC[11]*DC[18]); 
  vol[30] = 1.3416407864998738*(DC[23]*fC[39]+DC[22]*fC[38])+0.5999999999999999*(fC[13]*DC[39]+DC[14]*fC[35])+0.6*(fC[10]*DC[38]+DC[15]*fC[37])+0.6708203932499369*(fC[29]*DC[36]+fC[27]*DC[34]+DC[31]*fC[33]+DC[30]*fC[32]+fC[6]*DC[15]+fC[3]*DC[14]+DC[11]*fC[13]+DC[7]*fC[10])+1.6770509831248424*(DC[17]*fC[31]+DC[16]*fC[30]+DC[21]*fC[28]+DC[19]*fC[26])+0.47915742374995496*(fC[25]*DC[31]+fC[24]*DC[30]+DC[28]*fC[29]+DC[26]*fC[27])+0.75*(fC[1]*DC[31]+fC[0]*DC[30]+DC[5]*fC[29]+fC[8]*DC[28]+DC[2]*fC[27]+fC[4]*DC[26]+DC[12]*fC[25]+DC[9]*fC[24])+1.5*(fC[11]*DC[23]+fC[7]*DC[22]+fC[15]*DC[20]+fC[14]*DC[18]); 
  vol[31] = 1.3416407864998738*(DC[22]*fC[39]+DC[23]*fC[38])+0.6*(fC[10]*DC[39]+DC[15]*fC[35])+0.5999999999999999*(fC[13]*DC[38]+DC[14]*fC[37])+0.6708203932499369*(fC[27]*DC[36]+fC[29]*DC[34]+DC[30]*fC[33]+DC[31]*fC[32]+fC[3]*DC[15]+fC[6]*DC[14]+DC[7]*fC[13]+fC[10]*DC[11])+1.6770509831248424*(DC[16]*fC[31]+DC[17]*fC[30]+DC[19]*fC[28]+DC[21]*fC[26])+0.47915742374995496*(fC[24]*DC[31]+fC[25]*DC[30]+DC[26]*fC[29]+fC[27]*DC[28])+0.75*(fC[0]*DC[31]+fC[1]*DC[30]+DC[2]*fC[29]+fC[4]*DC[28]+DC[5]*fC[27]+fC[8]*DC[26]+DC[9]*fC[25]+DC[12]*fC[24])+1.5*(fC[7]*DC[23]+fC[11]*DC[22]+fC[14]*DC[20]+fC[15]*DC[18]); 
  vol[34] = 1.6770509831248424*(DC[20]*fC[39]+DC[18]*fC[38]+DC[17]*fC[36]+DC[16]*fC[34])+0.47915742374995496*(fC[37]*DC[39]+fC[35]*DC[38]+fC[33]*DC[36]+fC[32]*DC[34])+0.75*(fC[6]*DC[39]+fC[3]*DC[38]+DC[11]*fC[37]+fC[1]*DC[36]+DC[7]*fC[35]+fC[0]*DC[34]+DC[5]*fC[33]+DC[2]*fC[32])+0.6708203932499369*(fC[29]*DC[31]+fC[27]*DC[30]+fC[13]*DC[15]+fC[10]*DC[14]+fC[8]*DC[12]+fC[4]*DC[9])+1.5*(fC[15]*DC[23]+fC[14]*DC[22]+fC[12]*DC[21]+fC[9]*DC[19]); 
  vol[36] = 1.6770509831248424*(DC[18]*fC[39]+DC[20]*fC[38]+DC[16]*fC[36]+DC[17]*fC[34])+0.47915742374995496*(fC[35]*DC[39]+fC[37]*DC[38]+fC[32]*DC[36]+fC[33]*DC[34])+0.75*(fC[3]*DC[39]+fC[6]*DC[38]+DC[7]*fC[37]+fC[0]*DC[36]+DC[11]*fC[35]+fC[1]*DC[34]+DC[2]*fC[33]+DC[5]*fC[32])+0.6708203932499369*(fC[27]*DC[31]+fC[29]*DC[30]+fC[10]*DC[15]+fC[13]*DC[14]+fC[4]*DC[12]+fC[8]*DC[9])+1.5*(fC[14]*DC[23]+fC[15]*DC[22]+fC[9]*DC[21]+fC[12]*DC[19]); 
  vol[38] = 1.6770509831248424*(DC[17]*fC[39]+DC[16]*fC[38]+DC[20]*fC[36]+DC[18]*fC[34])+0.47915742374995496*(fC[33]*DC[39]+fC[32]*DC[38]+DC[36]*fC[37]+DC[34]*fC[35])+0.6708203932499369*(fC[25]*DC[39]+fC[24]*DC[38]+DC[28]*fC[37]+DC[26]*fC[35]+fC[8]*DC[15]+fC[4]*DC[14]+DC[12]*fC[13]+DC[9]*fC[10])+0.75*(fC[1]*DC[39]+fC[0]*DC[38]+DC[5]*fC[37]+fC[6]*DC[36]+DC[2]*fC[35]+fC[3]*DC[34]+DC[11]*fC[33]+DC[7]*fC[32])+1.3416407864998738*(DC[23]*fC[31]+DC[22]*fC[30])+0.5999999999999999*(fC[13]*DC[31]+DC[14]*fC[27])+0.6*(fC[10]*DC[30]+DC[15]*fC[29])+1.5*(fC[12]*DC[23]+fC[9]*DC[22]+fC[15]*DC[21]+fC[14]*DC[19]); 
  vol[39] = 1.6770509831248424*(DC[16]*fC[39]+DC[17]*fC[38]+DC[18]*fC[36]+DC[20]*fC[34])+0.47915742374995496*(fC[32]*DC[39]+fC[33]*DC[38]+DC[34]*fC[37]+fC[35]*DC[36])+0.6708203932499369*(fC[24]*DC[39]+fC[25]*DC[38]+DC[26]*fC[37]+DC[28]*fC[35]+fC[4]*DC[15]+fC[8]*DC[14]+DC[9]*fC[13]+fC[10]*DC[12])+0.75*(fC[0]*DC[39]+fC[1]*DC[38]+DC[2]*fC[37]+fC[3]*DC[36]+DC[5]*fC[35]+fC[6]*DC[34]+DC[7]*fC[33]+DC[11]*fC[32])+1.3416407864998738*(DC[22]*fC[31]+DC[23]*fC[30])+0.6*(fC[10]*DC[31]+DC[15]*fC[27])+0.5999999999999999*(fC[13]*DC[30]+DC[14]*fC[29])+1.5*(fC[9]*DC[23]+fC[12]*DC[22]+fC[14]*DC[21]+fC[15]*DC[19]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[4] += (0.5*vol[4]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[5] += (0.5*vol[5]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[7] += (0.5*vol[7]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[8] += (0.5*vol[8]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[9] += (0.5*vol[9]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[11] += (0.5*vol[11]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[12] += (0.5*vol[12]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[14] += (0.5*vol[14]-0.6123724356957945*surft2_up[6]+0.6123724356957945*(surft2_lo[6]+surft1_up[6]+surft1_lo[6]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]-0.6123724356957945*surft2_up[7]+0.6123724356957945*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv_inv_sq; 
  out[16] += (0.5*vol[16]-2.3717082451262845*(surft2_up[0]+surft2_lo[0])+0.7905694150420948*surft1_up[0]-0.7905694150420948*surft1_lo[0])*dv_inv_sq; 
  out[17] += (0.5*vol[17]-2.3717082451262845*(surft2_up[1]+surft2_lo[1])+0.7905694150420949*surft1_up[1]-0.7905694150420949*surft1_lo[1])*dv_inv_sq; 
  out[18] += (0.5*vol[18]-2.3717082451262845*(surft2_up[2]+surft2_lo[2])+0.7905694150420949*surft1_up[2]-0.7905694150420949*surft1_lo[2])*dv_inv_sq; 
  out[19] += (0.5*vol[19]-2.3717082451262845*(surft2_up[3]+surft2_lo[3])+0.7905694150420949*surft1_up[3]-0.7905694150420949*surft1_lo[3])*dv_inv_sq; 
  out[20] += (0.5*vol[20]-2.3717082451262845*(surft2_up[4]+surft2_lo[4])+0.7905694150420948*surft1_up[4]-0.7905694150420948*surft1_lo[4])*dv_inv_sq; 
  out[21] += (0.5*vol[21]-2.3717082451262845*(surft2_up[5]+surft2_lo[5])+0.7905694150420948*surft1_up[5]-0.7905694150420948*surft1_lo[5])*dv_inv_sq; 
  out[22] += (0.5*vol[22]-2.3717082451262845*(surft2_up[6]+surft2_lo[6])+0.7905694150420948*surft1_up[6]-0.7905694150420948*surft1_lo[6])*dv_inv_sq; 
  out[23] += (0.5*vol[23]-2.3717082451262845*(surft2_up[7]+surft2_lo[7])+0.7905694150420949*surft1_up[7]-0.7905694150420949*surft1_lo[7])*dv_inv_sq; 
  out[24] += 0.5*vol[24]*dv_inv_sq; 
  out[25] += 0.5*vol[25]*dv_inv_sq; 
  out[26] += 0.5*vol[26]*dv_inv_sq; 
  out[27] += 0.5*vol[27]*dv_inv_sq; 
  out[28] += 0.5*vol[28]*dv_inv_sq; 
  out[29] += 0.5*vol[29]*dv_inv_sq; 
  out[30] += 0.5*vol[30]*dv_inv_sq; 
  out[31] += 0.5*vol[31]*dv_inv_sq; 
  out[32] += 0.5*vol[32]*dv_inv_sq; 
  out[33] += 0.5*vol[33]*dv_inv_sq; 
  out[34] += 0.5*vol[34]*dv_inv_sq; 
  out[35] += 0.5*vol[35]*dv_inv_sq; 
  out[36] += 0.5*vol[36]*dv_inv_sq; 
  out[37] += 0.5*vol[37]*dv_inv_sq; 
  out[38] += 0.5*vol[38]*dv_inv_sq; 
  out[39] += 0.5*vol[39]*dv_inv_sq; 

  double cflFreq = fmax(fabs(D_surf_C[0]), fabs(D_surf_R[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
