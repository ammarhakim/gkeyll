#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* dxv, const double* diff_coeff_C, 
  const double* diff_coeff_surf_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_C: Diffusion tensor in center cell. 
  // diff_coeff_surf_stencil[9]: 9-cell stencil of surface expansion of recovered diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 

  double dv_inv_sq = 4.0/dxv[3]/dxv[3]; 

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

  const double* DC = &diff_coeff_C[320]; 
  const double* D_surf_C = &diff_coeff_surf_stencil[1][128]; 
  const double* D_surf_R = &diff_coeff_surf_stencil[2][128]; 

  f_rec_lo[0] = 0.408248290463863*fL[4]-0.408248290463863*fC[4]+0.3535533905932737*(fL[0]+fC[0]); 
  f_rec_lo[1] = 0.408248290463863*fL[8]-0.408248290463863*fC[8]+0.3535533905932737*(fL[1]+fC[1]); 
  f_rec_lo[2] = 0.408248290463863*fL[9]-0.408248290463863*fC[9]+0.3535533905932737*(fL[2]+fC[2]); 
  f_rec_lo[3] = 0.408248290463863*fL[10]-0.408248290463863*fC[10]+0.3535533905932737*(fL[3]+fC[3]); 
  f_rec_lo[4] = 0.408248290463863*fL[12]-0.408248290463863*fC[12]+0.3535533905932737*(fL[5]+fC[5]); 
  f_rec_lo[5] = 0.408248290463863*fL[13]-0.408248290463863*fC[13]+0.3535533905932737*(fL[6]+fC[6]); 
  f_rec_lo[6] = 0.408248290463863*fL[14]-0.408248290463863*fC[14]+0.3535533905932737*(fL[7]+fC[7]); 
  f_rec_lo[7] = 0.408248290463863*fL[15]-0.408248290463863*fC[15]+0.3535533905932737*(fL[11]+fC[11]); 
  f_rec_up[0] = -(0.408248290463863*fR[4])+0.408248290463863*fC[4]+0.3535533905932737*(fR[0]+fC[0]); 
  f_rec_up[1] = -(0.408248290463863*fR[8])+0.408248290463863*fC[8]+0.3535533905932737*(fR[1]+fC[1]); 
  f_rec_up[2] = -(0.408248290463863*fR[9])+0.408248290463863*fC[9]+0.3535533905932737*(fR[2]+fC[2]); 
  f_rec_up[3] = -(0.408248290463863*fR[10])+0.408248290463863*fC[10]+0.3535533905932737*(fR[3]+fC[3]); 
  f_rec_up[4] = -(0.408248290463863*fR[12])+0.408248290463863*fC[12]+0.3535533905932737*(fR[5]+fC[5]); 
  f_rec_up[5] = -(0.408248290463863*fR[13])+0.408248290463863*fC[13]+0.3535533905932737*(fR[6]+fC[6]); 
  f_rec_up[6] = -(0.408248290463863*fR[14])+0.408248290463863*fC[14]+0.3535533905932737*(fR[7]+fC[7]); 
  f_rec_up[7] = -(0.408248290463863*fR[15])+0.408248290463863*fC[15]+0.3535533905932737*(fR[11]+fC[11]); 

  df_rec_lo[0] = -(0.7654655446197428*(fL[4]+fC[4]))-0.7954951288348656*fL[0]+0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*(fL[8]+fC[8]))-0.7954951288348656*fL[1]+0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*(fL[9]+fC[9]))-0.7954951288348656*fL[2]+0.7954951288348656*fC[2]; 
  df_rec_lo[3] = -(0.7654655446197428*(fL[10]+fC[10]))-0.7954951288348656*fL[3]+0.7954951288348656*fC[3]; 
  df_rec_lo[4] = -(0.7654655446197428*(fL[12]+fC[12]))-0.7954951288348656*fL[5]+0.7954951288348656*fC[5]; 
  df_rec_lo[5] = -(0.7654655446197428*(fL[13]+fC[13]))-0.7954951288348656*fL[6]+0.7954951288348656*fC[6]; 
  df_rec_lo[6] = -(0.7654655446197428*(fL[14]+fC[14]))-0.7954951288348656*fL[7]+0.7954951288348656*fC[7]; 
  df_rec_lo[7] = -(0.7654655446197428*(fL[15]+fC[15]))-0.7954951288348656*fL[11]+0.7954951288348656*fC[11]; 
  df_rec_up[0] = -(0.7654655446197428*(fR[4]+fC[4]))+0.7954951288348656*fR[0]-0.7954951288348656*fC[0]; 
  df_rec_up[1] = -(0.7654655446197428*(fR[8]+fC[8]))+0.7954951288348656*fR[1]-0.7954951288348656*fC[1]; 
  df_rec_up[2] = -(0.7654655446197428*(fR[9]+fC[9]))+0.7954951288348656*fR[2]-0.7954951288348656*fC[2]; 
  df_rec_up[3] = -(0.7654655446197428*(fR[10]+fC[10]))+0.7954951288348656*fR[3]-0.7954951288348656*fC[3]; 
  df_rec_up[4] = -(0.7654655446197428*(fR[12]+fC[12]))+0.7954951288348656*fR[5]-0.7954951288348656*fC[5]; 
  df_rec_up[5] = -(0.7654655446197428*(fR[13]+fC[13]))+0.7954951288348656*fR[6]-0.7954951288348656*fC[6]; 
  df_rec_up[6] = -(0.7654655446197428*(fR[14]+fC[14]))+0.7954951288348656*fR[7]-0.7954951288348656*fC[7]; 
  df_rec_up[7] = -(0.7654655446197428*(fR[15]+fC[15]))+0.7954951288348656*fR[11]-0.7954951288348656*fC[11]; 

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

  vol[4] = 1.6770509831248421*(fC[15]*DC[39]+fC[10]*DC[35]+fC[9]*DC[34]+fC[8]*DC[33])+1.6770509831248424*(fC[14]*DC[38]+fC[13]*DC[37]+fC[12]*DC[36]+fC[4]*DC[32])+0.75*(fC[28]*DC[31]+fC[26]*DC[30]+fC[25]*DC[29]+fC[24]*DC[27]+fC[20]*DC[23]+fC[18]*DC[22]+fC[17]*DC[21]+fC[16]*DC[19]+fC[11]*DC[15]+fC[7]*DC[14]+fC[6]*DC[13]+fC[5]*DC[12]+fC[3]*DC[10]+fC[2]*DC[9]+fC[1]*DC[8]+fC[0]*DC[4]); 
  vol[8] = 1.6770509831248421*(fC[14]*DC[39]+fC[13]*DC[35]+fC[12]*DC[34]+fC[4]*DC[33])+1.6770509831248424*(fC[15]*DC[38]+fC[10]*DC[37]+fC[9]*DC[36]+fC[8]*DC[32])+0.75*(fC[26]*DC[31]+fC[28]*DC[30]+fC[24]*DC[29]+fC[25]*DC[27]+fC[18]*DC[23]+fC[20]*DC[22]+fC[16]*DC[21]+fC[17]*DC[19]+fC[7]*DC[15]+fC[11]*DC[14]+fC[3]*DC[13]+fC[2]*DC[12]+fC[6]*DC[10]+fC[5]*DC[9]+fC[0]*DC[8]+fC[1]*DC[4]); 
  vol[9] = 1.5*(fC[23]*DC[39]+fC[22]*DC[38]+fC[21]*DC[36]+fC[19]*DC[34])+1.6770509831248421*(fC[13]*DC[39]+fC[14]*DC[35]+fC[4]*DC[34]+fC[12]*DC[33])+1.6770509831248424*(fC[10]*DC[38]+fC[15]*DC[37]+fC[8]*DC[36]+fC[9]*DC[32])+0.75*(fC[25]*DC[31]+fC[24]*DC[30]+fC[28]*DC[29]+fC[26]*DC[27]+fC[6]*DC[15]+fC[3]*DC[14]+fC[11]*DC[13]+fC[1]*DC[12]+fC[7]*DC[10]+fC[0]*DC[9]+fC[5]*DC[8]+fC[2]*DC[4])+0.6708203932499369*(fC[11]*DC[23]+fC[7]*DC[22]+fC[5]*DC[21]+DC[15]*fC[20]+fC[2]*DC[19]+DC[14]*fC[18]+DC[12]*fC[17]+DC[9]*fC[16]); 
  vol[10] = 1.5*(fC[31]*DC[39]+fC[30]*DC[38]+fC[29]*DC[37]+fC[27]*DC[35])+1.6770509831248421*(fC[12]*DC[39]+fC[4]*DC[35]+fC[14]*DC[34]+fC[13]*DC[33])+1.6770509831248424*(fC[9]*DC[38]+fC[8]*DC[37]+fC[15]*DC[36]+fC[10]*DC[32])+0.6708203932499369*(fC[11]*DC[31]+fC[7]*DC[30]+fC[6]*DC[29]+DC[15]*fC[28]+fC[3]*DC[27]+DC[14]*fC[26]+DC[13]*fC[25]+DC[10]*fC[24])+0.75*(fC[17]*DC[23]+fC[16]*DC[22]+fC[20]*DC[21]+fC[18]*DC[19]+fC[5]*DC[15]+fC[2]*DC[14]+fC[1]*DC[13]+fC[11]*DC[12]+fC[0]*DC[10]+fC[7]*DC[9]+fC[6]*DC[8]+fC[3]*DC[4]); 
  vol[12] = 1.5*(fC[22]*DC[39]+fC[23]*DC[38]+fC[19]*DC[36]+fC[21]*DC[34])+1.6770509831248421*(fC[10]*DC[39]+fC[15]*DC[35]+fC[8]*DC[34]+fC[9]*DC[33])+1.6770509831248424*(fC[13]*DC[38]+fC[14]*DC[37]+fC[4]*DC[36]+fC[12]*DC[32])+0.75*(fC[24]*DC[31]+fC[25]*DC[30]+fC[26]*DC[29]+DC[27]*fC[28]+fC[3]*DC[15]+fC[6]*DC[14]+fC[7]*DC[13]+fC[0]*DC[12]+DC[10]*fC[11]+fC[1]*DC[9]+fC[2]*DC[8]+DC[4]*fC[5])+0.6708203932499369*(fC[7]*DC[23]+fC[11]*DC[22]+fC[2]*DC[21]+DC[14]*fC[20]+fC[5]*DC[19]+DC[15]*fC[18]+DC[9]*fC[17]+DC[12]*fC[16]); 
  vol[13] = 1.5*(fC[30]*DC[39]+fC[31]*DC[38]+fC[27]*DC[37]+fC[29]*DC[35])+1.6770509831248421*(fC[9]*DC[39]+fC[8]*DC[35]+fC[15]*DC[34]+fC[10]*DC[33])+1.6770509831248424*(fC[12]*DC[38]+fC[4]*DC[37]+fC[14]*DC[36]+fC[13]*DC[32])+0.6708203932499369*(fC[7]*DC[31]+fC[11]*DC[30]+fC[3]*DC[29]+DC[14]*fC[28]+fC[6]*DC[27]+DC[15]*fC[26]+DC[10]*fC[25]+DC[13]*fC[24])+0.75*(fC[16]*DC[23]+fC[17]*DC[22]+fC[18]*DC[21]+DC[19]*fC[20]+fC[2]*DC[15]+fC[5]*DC[14]+fC[0]*DC[13]+fC[7]*DC[12]+DC[9]*fC[11]+fC[1]*DC[10]+fC[3]*DC[8]+DC[4]*fC[6]); 
  vol[14] = 1.5*(fC[29]*DC[39]+fC[21]*DC[39]+fC[27]*DC[38]+fC[19]*DC[38]+fC[31]*DC[37]+fC[23]*DC[36]+fC[30]*DC[35]+fC[22]*DC[34])+1.6770509831248421*(fC[8]*DC[39]+fC[9]*DC[35]+fC[10]*DC[34]+fC[15]*DC[33])+1.6770509831248424*(fC[4]*DC[38]+fC[12]*DC[37]+fC[13]*DC[36]+fC[14]*DC[32])+0.5999999999999999*(fC[20]*DC[31]+fC[18]*DC[30]+DC[23]*fC[28]+DC[22]*fC[26])+0.6708203932499369*(fC[6]*DC[31]+fC[3]*DC[30]+fC[11]*DC[29]+DC[13]*fC[28]+fC[7]*DC[27]+DC[10]*fC[26]+DC[15]*fC[25]+DC[14]*fC[24]+fC[5]*DC[23]+fC[2]*DC[22]+fC[11]*DC[21]+DC[12]*fC[20]+fC[7]*DC[19]+DC[9]*fC[18]+DC[15]*fC[17]+DC[14]*fC[16])+0.75*(fC[1]*DC[15]+fC[0]*DC[14]+fC[5]*DC[13]+fC[6]*DC[12]+DC[8]*fC[11]+fC[2]*DC[10]+fC[3]*DC[9]+DC[4]*fC[7]); 
  vol[15] = 1.5*(fC[27]*DC[39]+fC[19]*DC[39]+fC[29]*DC[38]+fC[21]*DC[38]+fC[30]*DC[37]+fC[22]*DC[36]+fC[31]*DC[35]+fC[23]*DC[34])+1.6770509831248421*(fC[4]*DC[39]+fC[12]*DC[35]+fC[13]*DC[34]+fC[14]*DC[33])+1.6770509831248424*(fC[8]*DC[38]+fC[9]*DC[37]+fC[10]*DC[36]+fC[15]*DC[32])+0.6*(fC[18]*DC[31]+fC[20]*DC[30]+DC[22]*fC[28]+DC[23]*fC[26])+0.6708203932499369*(fC[3]*DC[31]+fC[6]*DC[30]+fC[7]*DC[29]+DC[10]*fC[28]+fC[11]*DC[27]+DC[13]*fC[26]+DC[14]*fC[25]+DC[15]*fC[24]+fC[2]*DC[23]+fC[5]*DC[22]+fC[7]*DC[21]+DC[9]*fC[20]+fC[11]*DC[19]+DC[12]*fC[18]+DC[14]*fC[17]+DC[15]*fC[16])+0.75*(fC[0]*DC[15]+fC[1]*DC[14]+fC[2]*DC[13]+fC[3]*DC[12]+DC[4]*fC[11]+fC[5]*DC[10]+fC[6]*DC[9]+fC[7]*DC[8]); 
  vol[19] = 1.5*(fC[15]*DC[39]+fC[14]*DC[38]+fC[12]*DC[36]+fC[9]*DC[34])+1.6770509831248424*(fC[23]*DC[37]+fC[22]*DC[35]+fC[21]*DC[33]+fC[19]*DC[32])+0.6708203932499369*(fC[28]*DC[31]+fC[26]*DC[30]+fC[11]*DC[15]+fC[7]*DC[14]+fC[5]*DC[12]+fC[2]*DC[9])+0.47915742374995496*(fC[20]*DC[23]+fC[18]*DC[22]+fC[17]*DC[21]+fC[16]*DC[19])+0.75*(fC[6]*DC[23]+fC[3]*DC[22]+fC[1]*DC[21]+DC[13]*fC[20]+fC[0]*DC[19]+DC[10]*fC[18]+DC[8]*fC[17]+DC[4]*fC[16]); 
  vol[21] = 1.5*(fC[14]*DC[39]+fC[15]*DC[38]+fC[9]*DC[36]+fC[12]*DC[34])+1.6770509831248424*(fC[22]*DC[37]+fC[23]*DC[35]+fC[19]*DC[33]+fC[21]*DC[32])+0.6708203932499369*(fC[26]*DC[31]+fC[28]*DC[30]+fC[7]*DC[15]+fC[11]*DC[14]+fC[2]*DC[12]+fC[5]*DC[9])+0.47915742374995496*(fC[18]*DC[23]+fC[20]*DC[22]+fC[16]*DC[21]+fC[17]*DC[19])+0.75*(fC[3]*DC[23]+fC[6]*DC[22]+fC[0]*DC[21]+DC[10]*fC[20]+fC[1]*DC[19]+DC[13]*fC[18]+DC[4]*fC[17]+DC[8]*fC[16]); 
  vol[22] = 1.3416407864998738*(fC[31]*DC[39]+fC[30]*DC[38])+1.5*(fC[12]*DC[39]+fC[9]*DC[38]+fC[15]*DC[36]+fC[14]*DC[34])+1.6770509831248424*(fC[21]*DC[37]+fC[19]*DC[35]+fC[23]*DC[33]+fC[22]*DC[32])+0.5999999999999999*(fC[11]*DC[31]+DC[14]*fC[26])+0.6*(fC[7]*DC[30]+DC[15]*fC[28])+0.6708203932499369*(fC[20]*DC[29]+fC[18]*DC[27]+DC[23]*fC[25]+DC[22]*fC[24]+fC[5]*DC[15]+fC[2]*DC[14]+fC[11]*DC[12]+fC[7]*DC[9])+0.47915742374995496*(fC[17]*DC[23]+fC[16]*DC[22]+fC[20]*DC[21]+fC[18]*DC[19])+0.75*(fC[1]*DC[23]+fC[0]*DC[22]+fC[6]*DC[21]+DC[8]*fC[20]+fC[3]*DC[19]+DC[4]*fC[18]+DC[13]*fC[17]+DC[10]*fC[16]); 
  vol[23] = 1.3416407864998738*(fC[30]*DC[39]+fC[31]*DC[38])+1.5*(fC[9]*DC[39]+fC[12]*DC[38]+fC[14]*DC[36]+fC[15]*DC[34])+1.6770509831248424*(fC[19]*DC[37]+fC[21]*DC[35]+fC[22]*DC[33]+fC[23]*DC[32])+0.6*(fC[7]*DC[31]+DC[15]*fC[26])+0.5999999999999999*(fC[11]*DC[30]+DC[14]*fC[28])+0.6708203932499369*(fC[18]*DC[29]+fC[20]*DC[27]+DC[22]*fC[25]+DC[23]*fC[24]+fC[2]*DC[15]+fC[5]*DC[14]+fC[7]*DC[12]+DC[9]*fC[11])+0.47915742374995496*(fC[16]*DC[23]+fC[17]*DC[22]+fC[18]*DC[21]+DC[19]*fC[20])+0.75*(fC[0]*DC[23]+fC[1]*DC[22]+fC[3]*DC[21]+DC[4]*fC[20]+fC[6]*DC[19]+DC[8]*fC[18]+DC[10]*fC[17]+DC[13]*fC[16]); 
  vol[27] = 1.5*(fC[15]*DC[39]+fC[14]*DC[38]+fC[13]*DC[37]+fC[10]*DC[35])+1.6770509831248424*(fC[31]*DC[36]+fC[30]*DC[34]+fC[29]*DC[33]+fC[27]*DC[32])+0.47915742374995496*(fC[28]*DC[31]+fC[26]*DC[30]+fC[25]*DC[29]+fC[24]*DC[27])+0.75*(fC[5]*DC[31]+fC[2]*DC[30]+fC[1]*DC[29]+DC[12]*fC[28]+fC[0]*DC[27]+DC[9]*fC[26]+DC[8]*fC[25]+DC[4]*fC[24])+0.6708203932499369*(fC[20]*DC[23]+fC[18]*DC[22]+fC[11]*DC[15]+fC[7]*DC[14]+fC[6]*DC[13]+fC[3]*DC[10]); 
  vol[29] = 1.5*(fC[14]*DC[39]+fC[15]*DC[38]+fC[10]*DC[37]+fC[13]*DC[35])+1.6770509831248424*(fC[30]*DC[36]+fC[31]*DC[34]+fC[27]*DC[33]+fC[29]*DC[32])+0.47915742374995496*(fC[26]*DC[31]+fC[28]*DC[30]+fC[24]*DC[29]+fC[25]*DC[27])+0.75*(fC[2]*DC[31]+fC[5]*DC[30]+fC[0]*DC[29]+DC[9]*fC[28]+fC[1]*DC[27]+DC[12]*fC[26]+DC[4]*fC[25]+DC[8]*fC[24])+0.6708203932499369*(fC[18]*DC[23]+fC[20]*DC[22]+fC[7]*DC[15]+fC[11]*DC[14]+fC[3]*DC[13]+fC[6]*DC[10]); 
  vol[30] = 1.3416407864998738*(fC[23]*DC[39]+fC[22]*DC[38])+1.5*(fC[13]*DC[39]+fC[10]*DC[38]+fC[15]*DC[37]+fC[14]*DC[35])+1.6770509831248424*(fC[29]*DC[36]+fC[27]*DC[34]+fC[31]*DC[33]+fC[30]*DC[32])+0.47915742374995496*(fC[25]*DC[31]+fC[24]*DC[30]+fC[28]*DC[29]+fC[26]*DC[27])+0.6708203932499369*(fC[17]*DC[31]+fC[16]*DC[30]+DC[21]*fC[28]+DC[19]*fC[26]+fC[6]*DC[15]+fC[3]*DC[14]+fC[11]*DC[13]+fC[7]*DC[10])+0.75*(fC[1]*DC[31]+fC[0]*DC[30]+fC[5]*DC[29]+DC[8]*fC[28]+fC[2]*DC[27]+DC[4]*fC[26]+DC[12]*fC[25]+DC[9]*fC[24])+0.5999999999999999*(fC[11]*DC[23]+DC[14]*fC[18])+0.6*(fC[7]*DC[22]+DC[15]*fC[20]); 
  vol[31] = 1.3416407864998738*(fC[22]*DC[39]+fC[23]*DC[38])+1.5*(fC[10]*DC[39]+fC[13]*DC[38]+fC[14]*DC[37]+fC[15]*DC[35])+1.6770509831248424*(fC[27]*DC[36]+fC[29]*DC[34]+fC[30]*DC[33]+fC[31]*DC[32])+0.47915742374995496*(fC[24]*DC[31]+fC[25]*DC[30]+fC[26]*DC[29]+DC[27]*fC[28])+0.6708203932499369*(fC[16]*DC[31]+fC[17]*DC[30]+DC[19]*fC[28]+DC[21]*fC[26]+fC[3]*DC[15]+fC[6]*DC[14]+fC[7]*DC[13]+DC[10]*fC[11])+0.75*(fC[0]*DC[31]+fC[1]*DC[30]+fC[2]*DC[29]+DC[4]*fC[28]+fC[5]*DC[27]+DC[8]*fC[26]+DC[9]*fC[25]+DC[12]*fC[24])+0.6*(fC[7]*DC[23]+DC[15]*fC[18])+0.5999999999999999*(fC[11]*DC[22]+DC[14]*fC[20]); 
  vol[32] = 5.031152949374527*(DC[39]*fC[39]+DC[38]*fC[38]+DC[37]*fC[37]+DC[36]*fC[36]+DC[35]*fC[35]+DC[34]*fC[34]+DC[33]*fC[33]+DC[32]*fC[32])+3.7500000000000004*(fC[11]*DC[39]+fC[3]*DC[35]+fC[2]*DC[34]+fC[1]*DC[33])+3.75*(fC[7]*DC[38]+fC[6]*DC[37]+fC[5]*DC[36]+fC[0]*DC[32])+3.3541019662496847*(DC[31]*fC[31]+DC[30]*fC[30]+DC[29]*fC[29]+DC[27]*fC[27]+DC[23]*fC[23]+DC[22]*fC[22]+DC[21]*fC[21]+DC[19]*fC[19]+DC[15]*fC[15]+DC[14]*fC[14]+DC[13]*fC[13]+DC[12]*fC[12]+DC[10]*fC[10]+DC[9]*fC[9]+DC[8]*fC[8]+DC[4]*fC[4])+1.6770509831248424*(DC[28]*fC[28]+DC[26]*fC[26]+DC[25]*fC[25]+DC[24]*fC[24]+DC[20]*fC[20]+DC[18]*fC[18]+DC[17]*fC[17]+DC[16]*fC[16]+DC[11]*fC[11]+DC[7]*fC[7]+DC[6]*fC[6]+DC[5]*fC[5]+DC[3]*fC[3]+DC[2]*fC[2]+DC[1]*fC[1]+DC[0]*fC[0]); 
  vol[33] = 5.031152949374527*(DC[38]*fC[39]+fC[38]*DC[39]+DC[35]*fC[37]+fC[35]*DC[37]+DC[34]*fC[36]+fC[34]*DC[36]+DC[32]*fC[33]+fC[32]*DC[33])+3.75*(fC[7]*DC[39]+fC[6]*DC[35]+fC[5]*DC[34]+fC[0]*DC[33])+3.7500000000000004*(fC[11]*DC[38]+fC[3]*DC[37]+fC[2]*DC[36]+fC[1]*DC[32])+3.3541019662496847*(DC[30]*fC[31]+fC[30]*DC[31]+DC[27]*fC[29]+fC[27]*DC[29]+DC[22]*fC[23]+fC[22]*DC[23]+DC[19]*fC[21]+fC[19]*DC[21])+1.6770509831248424*(DC[26]*fC[28]+fC[26]*DC[28]+DC[24]*fC[25]+fC[24]*DC[25]+DC[18]*fC[20]+fC[18]*DC[20]+DC[16]*fC[17]+fC[16]*DC[17])+3.3541019662496843*(DC[14]*fC[15]+fC[14]*DC[15]+DC[10]*fC[13]+fC[10]*DC[13]+DC[9]*fC[12]+fC[9]*DC[12]+DC[4]*fC[8]+fC[4]*DC[8])+1.6770509831248421*(DC[7]*fC[11]+fC[7]*DC[11]+DC[3]*fC[6]+fC[3]*DC[6]+DC[2]*fC[5]+fC[2]*DC[5]+DC[0]*fC[1]+fC[0]*DC[1]); 
  vol[34] = 5.031152949374527*(DC[37]*fC[39]+fC[37]*DC[39]+DC[35]*fC[38]+fC[35]*DC[38]+DC[33]*fC[36]+fC[33]*DC[36]+DC[32]*fC[34]+fC[32]*DC[34])+3.3541019662496847*(fC[20]*DC[39]+fC[18]*DC[38]+fC[17]*DC[36]+fC[16]*DC[34]+DC[29]*fC[31]+fC[29]*DC[31]+DC[27]*fC[30]+fC[27]*DC[30])+3.75*(fC[6]*DC[39]+fC[7]*DC[35]+fC[0]*DC[34]+fC[5]*DC[33])+3.7500000000000004*(fC[3]*DC[38]+fC[11]*DC[37]+fC[1]*DC[36]+fC[2]*DC[32])+1.6770509831248424*(DC[25]*fC[28]+fC[25]*DC[28]+DC[24]*fC[26]+fC[24]*DC[26])+3.0*(DC[15]*fC[23]+fC[15]*DC[23]+DC[14]*fC[22]+fC[14]*DC[22]+DC[12]*fC[21]+fC[12]*DC[21]+DC[9]*fC[19]+fC[9]*DC[19])+1.5*(DC[11]*fC[20]+fC[11]*DC[20]+DC[7]*fC[18]+fC[7]*DC[18]+DC[5]*fC[17]+fC[5]*DC[17]+DC[2]*fC[16]+fC[2]*DC[16])+3.3541019662496843*(DC[13]*fC[15]+fC[13]*DC[15]+DC[10]*fC[14]+fC[10]*DC[14]+DC[8]*fC[12]+fC[8]*DC[12]+DC[4]*fC[9]+fC[4]*DC[9])+1.6770509831248421*(DC[6]*fC[11]+fC[6]*DC[11]+DC[3]*fC[7]+fC[3]*DC[7]+DC[1]*fC[5]+fC[1]*DC[5]+DC[0]*fC[2]+fC[0]*DC[2]); 
  vol[35] = 5.031152949374527*(DC[36]*fC[39]+fC[36]*DC[39]+DC[34]*fC[38]+fC[34]*DC[38]+DC[33]*fC[37]+fC[33]*DC[37]+DC[32]*fC[35]+fC[32]*DC[35])+3.3541019662496847*(fC[28]*DC[39]+fC[26]*DC[38]+fC[25]*DC[37]+fC[24]*DC[35]+DC[21]*fC[23]+fC[21]*DC[23]+DC[19]*fC[22]+fC[19]*DC[22])+3.75*(fC[5]*DC[39]+fC[0]*DC[35]+fC[7]*DC[34]+fC[6]*DC[33])+3.7500000000000004*(fC[2]*DC[38]+fC[1]*DC[37]+fC[11]*DC[36]+fC[3]*DC[32])+3.0*(DC[15]*fC[31]+fC[15]*DC[31]+DC[14]*fC[30]+fC[14]*DC[30]+DC[13]*fC[29]+fC[13]*DC[29]+DC[10]*fC[27]+fC[10]*DC[27])+1.5*(DC[11]*fC[28]+fC[11]*DC[28]+DC[7]*fC[26]+fC[7]*DC[26]+DC[6]*fC[25]+fC[6]*DC[25]+DC[3]*fC[24]+fC[3]*DC[24])+1.6770509831248424*(DC[17]*fC[20]+fC[17]*DC[20]+DC[16]*fC[18]+fC[16]*DC[18])+3.3541019662496843*(DC[12]*fC[15]+fC[12]*DC[15]+DC[9]*fC[14]+fC[9]*DC[14]+DC[8]*fC[13]+fC[8]*DC[13]+DC[4]*fC[10]+fC[4]*DC[10])+1.6770509831248421*(DC[5]*fC[11]+fC[5]*DC[11]+DC[2]*fC[7]+fC[2]*DC[7]+DC[1]*fC[6]+fC[1]*DC[6]+DC[0]*fC[3]+fC[0]*DC[3]); 
  vol[36] = 5.031152949374527*(DC[35]*fC[39]+fC[35]*DC[39]+DC[37]*fC[38]+fC[37]*DC[38]+DC[32]*fC[36]+fC[32]*DC[36]+DC[33]*fC[34]+fC[33]*DC[34])+3.3541019662496847*(fC[18]*DC[39]+fC[20]*DC[38]+fC[16]*DC[36]+fC[17]*DC[34]+DC[27]*fC[31]+fC[27]*DC[31]+DC[29]*fC[30]+fC[29]*DC[30]+DC[10]*fC[15]+fC[10]*DC[15]+DC[13]*fC[14]+fC[13]*DC[14]+DC[4]*fC[12]+fC[4]*DC[12]+DC[8]*fC[9]+fC[8]*DC[9])+3.7500000000000004*(fC[3]*DC[39]+fC[11]*DC[35]+fC[1]*DC[34]+fC[2]*DC[33])+3.75*(fC[6]*DC[38]+fC[7]*DC[37]+fC[0]*DC[36]+fC[5]*DC[32])+1.6770509831248424*(DC[24]*fC[28]+fC[24]*DC[28]+DC[25]*fC[26]+fC[25]*DC[26]+DC[3]*fC[11]+fC[3]*DC[11]+DC[6]*fC[7]+fC[6]*DC[7]+DC[0]*fC[5]+fC[0]*DC[5]+DC[1]*fC[2]+fC[1]*DC[2])+3.0*(DC[14]*fC[23]+fC[14]*DC[23]+DC[15]*fC[22]+fC[15]*DC[22]+DC[9]*fC[21]+fC[9]*DC[21]+DC[12]*fC[19]+fC[12]*DC[19])+1.5*(DC[7]*fC[20]+fC[7]*DC[20]+DC[11]*fC[18]+fC[11]*DC[18]+DC[2]*fC[17]+fC[2]*DC[17]+DC[5]*fC[16]+fC[5]*DC[16]); 
  vol[37] = 5.031152949374527*(DC[34]*fC[39]+fC[34]*DC[39]+DC[36]*fC[38]+fC[36]*DC[38]+DC[32]*fC[37]+fC[32]*DC[37]+DC[33]*fC[35]+fC[33]*DC[35])+3.3541019662496847*(fC[26]*DC[39]+fC[28]*DC[38]+fC[24]*DC[37]+fC[25]*DC[35]+DC[19]*fC[23]+fC[19]*DC[23]+DC[21]*fC[22]+fC[21]*DC[22]+DC[9]*fC[15]+fC[9]*DC[15]+DC[12]*fC[14]+fC[12]*DC[14]+DC[4]*fC[13]+fC[4]*DC[13]+DC[8]*fC[10]+fC[8]*DC[10])+3.7500000000000004*(fC[2]*DC[39]+fC[1]*DC[35]+fC[11]*DC[34]+fC[3]*DC[33])+3.75*(fC[5]*DC[38]+fC[0]*DC[37]+fC[7]*DC[36]+fC[6]*DC[32])+3.0*(DC[14]*fC[31]+fC[14]*DC[31]+DC[15]*fC[30]+fC[15]*DC[30]+DC[10]*fC[29]+fC[10]*DC[29]+DC[13]*fC[27]+fC[13]*DC[27])+1.5*(DC[7]*fC[28]+fC[7]*DC[28]+DC[11]*fC[26]+fC[11]*DC[26]+DC[3]*fC[25]+fC[3]*DC[25]+DC[6]*fC[24]+fC[6]*DC[24])+1.6770509831248424*(DC[16]*fC[20]+fC[16]*DC[20]+DC[17]*fC[18]+fC[17]*DC[18]+DC[2]*fC[11]+fC[2]*DC[11]+DC[5]*fC[7]+fC[5]*DC[7]+DC[0]*fC[6]+fC[0]*DC[6]+DC[1]*fC[3]+fC[1]*DC[3]); 
  vol[38] = 5.031152949374527*(DC[33]*fC[39]+fC[33]*DC[39]+DC[32]*fC[38]+fC[32]*DC[38]+DC[36]*fC[37]+fC[36]*DC[37]+DC[34]*fC[35]+fC[34]*DC[35])+3.3541019662496847*(fC[25]*DC[39]+fC[17]*DC[39]+fC[24]*DC[38]+fC[16]*DC[38]+fC[28]*DC[37]+fC[20]*DC[36]+fC[26]*DC[35]+fC[18]*DC[34]+DC[8]*fC[15]+fC[8]*DC[15]+DC[4]*fC[14]+fC[4]*DC[14]+DC[12]*fC[13]+fC[12]*DC[13]+DC[9]*fC[10]+fC[9]*DC[10])+3.7500000000000004*(fC[1]*DC[39]+fC[2]*DC[35]+fC[3]*DC[34]+fC[11]*DC[33])+3.75*(fC[0]*DC[38]+fC[5]*DC[37]+fC[6]*DC[36]+fC[7]*DC[32])+2.6832815729997477*(DC[23]*fC[31]+fC[23]*DC[31]+DC[22]*fC[30]+fC[22]*DC[30])+3.0*(DC[13]*fC[31]+fC[13]*DC[31]+DC[10]*fC[30]+fC[10]*DC[30]+DC[15]*fC[29]+fC[15]*DC[29]+DC[14]*fC[27]+fC[14]*DC[27]+DC[12]*fC[23]+fC[12]*DC[23]+DC[9]*fC[22]+fC[9]*DC[22]+DC[15]*fC[21]+fC[15]*DC[21]+DC[14]*fC[19]+fC[14]*DC[19])+1.3416407864998738*(DC[20]*fC[28]+fC[20]*DC[28]+DC[18]*fC[26]+fC[18]*DC[26])+1.5*(DC[6]*fC[28]+fC[6]*DC[28]+DC[3]*fC[26]+fC[3]*DC[26]+DC[11]*fC[25]+fC[11]*DC[25]+DC[7]*fC[24]+fC[7]*DC[24]+DC[5]*fC[20]+fC[5]*DC[20]+DC[2]*fC[18]+fC[2]*DC[18]+DC[11]*fC[17]+fC[11]*DC[17]+DC[7]*fC[16]+fC[7]*DC[16])+1.6770509831248424*(DC[1]*fC[11]+fC[1]*DC[11]+DC[0]*fC[7]+fC[0]*DC[7]+DC[5]*fC[6]+fC[5]*DC[6]+DC[2]*fC[3]+fC[2]*DC[3]); 
  vol[39] = 5.031152949374527*(DC[32]*fC[39]+fC[32]*DC[39]+DC[33]*fC[38]+fC[33]*DC[38]+DC[34]*fC[37]+fC[34]*DC[37]+DC[35]*fC[36]+fC[35]*DC[36])+3.3541019662496847*(fC[24]*DC[39]+fC[16]*DC[39]+fC[25]*DC[38]+fC[17]*DC[38]+fC[26]*DC[37]+fC[18]*DC[36]+fC[28]*DC[35]+fC[20]*DC[34])+3.75*(fC[0]*DC[39]+fC[5]*DC[35]+fC[6]*DC[34]+fC[7]*DC[33])+3.7500000000000004*(fC[1]*DC[38]+fC[2]*DC[37]+fC[3]*DC[36]+fC[11]*DC[32])+2.6832815729997477*(DC[22]*fC[31]+fC[22]*DC[31]+DC[23]*fC[30]+fC[23]*DC[30])+3.0*(DC[10]*fC[31]+fC[10]*DC[31]+DC[13]*fC[30]+fC[13]*DC[30]+DC[14]*fC[29]+fC[14]*DC[29]+DC[15]*fC[27]+fC[15]*DC[27]+DC[9]*fC[23]+fC[9]*DC[23]+DC[12]*fC[22]+fC[12]*DC[22]+DC[14]*fC[21]+fC[14]*DC[21]+DC[15]*fC[19]+fC[15]*DC[19])+1.3416407864998738*(DC[18]*fC[28]+fC[18]*DC[28]+DC[20]*fC[26]+fC[20]*DC[26])+1.5*(DC[3]*fC[28]+fC[3]*DC[28]+DC[6]*fC[26]+fC[6]*DC[26]+DC[7]*fC[25]+fC[7]*DC[25]+DC[11]*fC[24]+fC[11]*DC[24]+DC[2]*fC[20]+fC[2]*DC[20]+DC[5]*fC[18]+fC[5]*DC[18]+DC[7]*fC[17]+fC[7]*DC[17]+DC[11]*fC[16]+fC[11]*DC[16])+3.3541019662496843*(DC[4]*fC[15]+fC[4]*DC[15]+DC[8]*fC[14]+fC[8]*DC[14]+DC[9]*fC[13]+fC[9]*DC[13]+DC[10]*fC[12]+fC[10]*DC[12])+1.6770509831248421*(DC[0]*fC[11]+fC[0]*DC[11]+DC[1]*fC[7]+fC[1]*DC[7]+DC[2]*fC[6]+fC[2]*DC[6]+DC[3]*fC[5]+fC[3]*DC[5]); 

  out[0] += (0.5*vol[0]+0.35355339059327373*surft1_up[0]-0.35355339059327373*surft1_lo[0])*dv_inv_sq; 
  out[1] += (0.5*vol[1]+0.35355339059327373*surft1_up[1]-0.35355339059327373*surft1_lo[1])*dv_inv_sq; 
  out[2] += (0.5*vol[2]+0.35355339059327373*surft1_up[2]-0.35355339059327373*surft1_lo[2])*dv_inv_sq; 
  out[3] += (0.5*vol[3]+0.35355339059327373*surft1_up[3]-0.35355339059327373*surft1_lo[3])*dv_inv_sq; 
  out[4] += (0.5*vol[4]-0.6123724356957945*surft2_up[0]+0.6123724356957945*(surft2_lo[0]+surft1_up[0]+surft1_lo[0]))*dv_inv_sq; 
  out[5] += (0.5*vol[5]+0.35355339059327373*surft1_up[4]-0.35355339059327373*surft1_lo[4])*dv_inv_sq; 
  out[6] += (0.5*vol[6]+0.35355339059327373*surft1_up[5]-0.35355339059327373*surft1_lo[5])*dv_inv_sq; 
  out[7] += (0.5*vol[7]+0.35355339059327373*surft1_up[6]-0.35355339059327373*surft1_lo[6])*dv_inv_sq; 
  out[8] += (0.5*vol[8]-0.6123724356957945*surft2_up[1]+0.6123724356957945*(surft2_lo[1]+surft1_up[1]+surft1_lo[1]))*dv_inv_sq; 
  out[9] += (0.5*vol[9]-0.6123724356957945*surft2_up[2]+0.6123724356957945*(surft2_lo[2]+surft1_up[2]+surft1_lo[2]))*dv_inv_sq; 
  out[10] += (0.5*vol[10]-0.6123724356957945*surft2_up[3]+0.6123724356957945*(surft2_lo[3]+surft1_up[3]+surft1_lo[3]))*dv_inv_sq; 
  out[11] += (0.5*vol[11]+0.35355339059327373*surft1_up[7]-0.35355339059327373*surft1_lo[7])*dv_inv_sq; 
  out[12] += (0.5*vol[12]-0.6123724356957945*surft2_up[4]+0.6123724356957945*(surft2_lo[4]+surft1_up[4]+surft1_lo[4]))*dv_inv_sq; 
  out[13] += (0.5*vol[13]-0.6123724356957945*surft2_up[5]+0.6123724356957945*(surft2_lo[5]+surft1_up[5]+surft1_lo[5]))*dv_inv_sq; 
  out[14] += (0.5*vol[14]-0.6123724356957945*surft2_up[6]+0.6123724356957945*(surft2_lo[6]+surft1_up[6]+surft1_lo[6]))*dv_inv_sq; 
  out[15] += (0.5*vol[15]-0.6123724356957945*surft2_up[7]+0.6123724356957945*(surft2_lo[7]+surft1_up[7]+surft1_lo[7]))*dv_inv_sq; 
  out[16] += 0.5*vol[16]*dv_inv_sq; 
  out[17] += 0.5*vol[17]*dv_inv_sq; 
  out[18] += 0.5*vol[18]*dv_inv_sq; 
  out[19] += 0.5*vol[19]*dv_inv_sq; 
  out[20] += 0.5*vol[20]*dv_inv_sq; 
  out[21] += 0.5*vol[21]*dv_inv_sq; 
  out[22] += 0.5*vol[22]*dv_inv_sq; 
  out[23] += 0.5*vol[23]*dv_inv_sq; 
  out[24] += 0.5*vol[24]*dv_inv_sq; 
  out[25] += 0.5*vol[25]*dv_inv_sq; 
  out[26] += 0.5*vol[26]*dv_inv_sq; 
  out[27] += 0.5*vol[27]*dv_inv_sq; 
  out[28] += 0.5*vol[28]*dv_inv_sq; 
  out[29] += 0.5*vol[29]*dv_inv_sq; 
  out[30] += 0.5*vol[30]*dv_inv_sq; 
  out[31] += 0.5*vol[31]*dv_inv_sq; 
  out[32] += (0.5*vol[32]-2.3717082451262845*(surft2_up[0]+surft2_lo[0])+0.7905694150420948*surft1_up[0]-0.7905694150420948*surft1_lo[0])*dv_inv_sq; 
  out[33] += (0.5*vol[33]-2.3717082451262845*(surft2_up[1]+surft2_lo[1])+0.7905694150420949*surft1_up[1]-0.7905694150420949*surft1_lo[1])*dv_inv_sq; 
  out[34] += (0.5*vol[34]-2.3717082451262845*(surft2_up[2]+surft2_lo[2])+0.7905694150420949*surft1_up[2]-0.7905694150420949*surft1_lo[2])*dv_inv_sq; 
  out[35] += (0.5*vol[35]-2.3717082451262845*(surft2_up[3]+surft2_lo[3])+0.7905694150420949*surft1_up[3]-0.7905694150420949*surft1_lo[3])*dv_inv_sq; 
  out[36] += (0.5*vol[36]-2.3717082451262845*(surft2_up[4]+surft2_lo[4])+0.7905694150420948*surft1_up[4]-0.7905694150420948*surft1_lo[4])*dv_inv_sq; 
  out[37] += (0.5*vol[37]-2.3717082451262845*(surft2_up[5]+surft2_lo[5])+0.7905694150420948*surft1_up[5]-0.7905694150420948*surft1_lo[5])*dv_inv_sq; 
  out[38] += (0.5*vol[38]-2.3717082451262845*(surft2_up[6]+surft2_lo[6])+0.7905694150420948*surft1_up[6]-0.7905694150420948*surft1_lo[6])*dv_inv_sq; 
  out[39] += (0.5*vol[39]-2.3717082451262845*(surft2_up[7]+surft2_lo[7])+0.7905694150420949*surft1_up[7]-0.7905694150420949*surft1_lo[7])*dv_inv_sq; 

  double cflFreq = fmax(fabs(D_surf_C[0]), fabs(D_surf_R[0])); 

  return 1.4142135623730947*dv_inv_sq*cflFreq; 
} 
