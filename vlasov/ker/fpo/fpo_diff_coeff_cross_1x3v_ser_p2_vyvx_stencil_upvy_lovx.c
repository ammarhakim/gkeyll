#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_upvy_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[2]; 
  double dv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[2]/dxv[1]; 
 
  double surft1_lo[20] = {0.0}; 
  double surft1_up[20] = {0.0}; 
  double surft2_lo[20] = {0.0}; 
  double surft2_up[20] = {0.0}; 
  double vol[48] = {0.0}; 
  double *out = &diff_coeff[144]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vy = &fpo_g_surf_stencil[4][20]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvy_surf_CC_vx = &fpo_dgdv_surf[0]; 

  surft1_lo[0] = dGdvy_surf_CC_vx[0]/dv1; 
  surft1_lo[1] = dGdvy_surf_CC_vx[1]/dv1; 
  surft1_lo[2] = dGdvy_surf_CC_vx[2]/dv1; 
  surft1_lo[3] = dGdvy_surf_CC_vx[3]/dv1; 
  surft1_lo[4] = dGdvy_surf_CC_vx[4]/dv1; 
  surft1_lo[5] = dGdvy_surf_CC_vx[5]/dv1; 
  surft1_lo[6] = dGdvy_surf_CC_vx[6]/dv1; 
  surft1_lo[7] = dGdvy_surf_CC_vx[7]/dv1; 
  surft1_lo[8] = dGdvy_surf_CC_vx[8]/dv1; 
  surft1_lo[9] = dGdvy_surf_CC_vx[9]/dv1; 
  surft1_lo[10] = dGdvy_surf_CC_vx[10]/dv1; 
  surft1_lo[11] = dGdvy_surf_CC_vx[11]/dv1; 
  surft1_lo[12] = dGdvy_surf_CC_vx[12]/dv1; 
  surft1_lo[13] = dGdvy_surf_CC_vx[13]/dv1; 
  surft1_lo[14] = dGdvy_surf_CC_vx[14]/dv1; 
  surft1_lo[15] = dGdvy_surf_CC_vx[15]/dv1; 
  surft1_lo[16] = dGdvy_surf_CC_vx[16]/dv1; 
  surft1_lo[17] = dGdvy_surf_CC_vx[17]/dv1; 
  surft1_lo[18] = dGdvy_surf_CC_vx[18]/dv1; 
  surft1_lo[19] = dGdvy_surf_CC_vx[19]/dv1; 
  surft1_up[0] = 0.599071547271275*GCR[22]+0.599071547271275*GCC[22]-0.8617863895711042*GCR[7]+0.8617863895711042*GCC[7]+0.6123724356957944*GCR[3]+0.6123724356957944*GCC[3]; 
  surft1_up[1] = 0.599071547271275*GCR[33]+0.599071547271275*GCC[33]-0.8617863895711042*GCR[15]+0.8617863895711042*GCC[15]+0.6123724356957944*GCR[6]+0.6123724356957944*GCC[6]; 
  surft1_up[2] = -(1.9270129491651047*GCR[24])+1.9270129491651047*GCC[24]+1.369306393762915*GCR[13]+1.369306393762915*GCC[13]; 
  surft1_up[3] = 0.599071547271275*GCR[38]+0.599071547271275*GCC[38]-0.8617863895711042*GCR[18]+0.8617863895711042*GCC[18]+0.6123724356957944*GCR[10]+0.6123724356957944*GCC[10]; 
  surft1_up[4] = -(1.927012949165105*GCR[34])+1.927012949165105*GCC[34]+1.369306393762915*GCR[23]+1.369306393762915*GCC[23]; 
  surft1_up[5] = 0.599071547271275*GCR[45]+0.599071547271275*GCC[45]-0.8617863895711042*GCR[31]+0.8617863895711042*GCC[31]+0.6123724356957944*GCR[17]+0.6123724356957944*GCC[17]; 
  surft1_up[6] = -(1.927012949165105*GCR[40])+1.927012949165105*GCC[40]+1.369306393762915*GCR[27]+1.369306393762915*GCC[27]; 
  surft1_up[7] = -(0.8617863895711042*GCR[32])+0.8617863895711042*GCC[32]+0.6123724356957944*GCR[21]+0.6123724356957944*GCC[21]; 
  surft1_up[9] = -(0.8617863895711042*GCR[43])+0.8617863895711042*GCC[43]+0.6123724356957944*GCR[30]+0.6123724356957944*GCC[30]; 
  surft1_up[10] = -(1.9270129491651047*GCR[46])+1.9270129491651047*GCC[46]+1.369306393762915*GCR[39]+1.369306393762915*GCC[39]; 
  surft1_up[13] = -(0.8617863895711042*GCR[44])+0.8617863895711042*GCC[44]+0.6123724356957944*GCR[37]+0.6123724356957944*GCC[37]; 
  surft1_up[15] = -(0.8617863895711042*GCR[47])+0.8617863895711042*GCC[47]+0.6123724356957944*GCR[42]+0.6123724356957944*GCC[42]; 

  surft2_lo[0] = 0.34587411908091625*GCC[13]+0.34587411908091625*GBC[13]-0.49755260400283263*GCC[3]+0.49755260400283263*GBC[3]+0.3535533905932737*GCC[0]+0.3535533905932737*GBC[0]; 
  surft2_lo[1] = 0.34587411908091625*GCC[23]+0.34587411908091625*GBC[23]-0.49755260400283263*GCC[6]+0.49755260400283263*GBC[6]+0.3535533905932737*GCC[1]+0.3535533905932737*GBC[1]; 
  surft2_lo[2] = 0.34587411908091625*GCC[24]+0.34587411908091625*GBC[24]-0.49755260400283263*GCC[7]+0.49755260400283263*GBC[7]+0.3535533905932737*GCC[2]+0.3535533905932737*GBC[2]; 
  surft2_lo[3] = 0.34587411908091625*GCC[27]+0.34587411908091625*GBC[27]-0.49755260400283263*GCC[10]+0.49755260400283263*GBC[10]+0.3535533905932737*GCC[4]+0.3535533905932737*GBC[4]; 
  surft2_lo[4] = 0.34587411908091625*GCC[34]+0.34587411908091625*GBC[34]-0.49755260400283263*GCC[15]+0.49755260400283263*GBC[15]+0.3535533905932737*GCC[5]+0.3535533905932737*GBC[5]; 
  surft2_lo[5] = 0.34587411908091625*GCC[39]+0.34587411908091625*GBC[39]-0.49755260400283263*GCC[17]+0.49755260400283263*GBC[17]+0.3535533905932737*GCC[8]+0.3535533905932737*GBC[8]; 
  surft2_lo[6] = 0.34587411908091625*GCC[40]+0.34587411908091625*GBC[40]-0.49755260400283263*GCC[18]+0.49755260400283263*GBC[18]+0.3535533905932737*GCC[9]+0.3535533905932737*GBC[9]; 
  surft2_lo[7] = -(0.49755260400283263*GCC[21])+0.49755260400283263*GBC[21]+0.3535533905932737*GCC[11]+0.3535533905932737*GBC[11]; 
  surft2_lo[8] = -(0.49755260400283263*GCC[22])+0.49755260400283263*GBC[22]+0.3535533905932737*GCC[12]+0.3535533905932737*GBC[12]; 
  surft2_lo[9] = -(0.49755260400283263*GCC[30])+0.49755260400283263*GBC[30]+0.3535533905932737*GCC[14]+0.3535533905932737*GBC[14]; 
  surft2_lo[10] = 0.34587411908091625*GCC[46]+0.34587411908091625*GBC[46]-0.49755260400283263*GCC[31]+0.49755260400283263*GBC[31]+0.3535533905932737*GCC[16]+0.3535533905932737*GBC[16]; 
  surft2_lo[11] = -(0.49755260400283263*GCC[32])+0.49755260400283263*GBC[32]+0.3535533905932737*GCC[19]+0.3535533905932737*GBC[19]; 
  surft2_lo[12] = -(0.49755260400283263*GCC[33])+0.49755260400283263*GBC[33]+0.3535533905932737*GCC[20]+0.3535533905932737*GBC[20]; 
  surft2_lo[13] = -(0.49755260400283263*GCC[37])+0.49755260400283263*GBC[37]+0.3535533905932737*GCC[25]+0.3535533905932737*GBC[25]; 
  surft2_lo[14] = -(0.49755260400283263*GCC[38])+0.49755260400283263*GBC[38]+0.3535533905932737*GCC[26]+0.3535533905932737*GBC[26]; 
  surft2_lo[15] = -(0.49755260400283263*GCC[42])+0.49755260400283263*GBC[42]+0.3535533905932737*GCC[28]+0.3535533905932737*GBC[28]; 
  surft2_lo[16] = -(0.49755260400283263*GCC[43])+0.49755260400283263*GBC[43]+0.3535533905932737*GCC[29]+0.3535533905932737*GBC[29]; 
  surft2_lo[17] = -(0.49755260400283263*GCC[44])+0.49755260400283263*GBC[44]+0.3535533905932737*GCC[35]+0.3535533905932737*GBC[35]; 
  surft2_lo[18] = -(0.49755260400283263*GCC[45])+0.49755260400283263*GBC[45]+0.3535533905932737*GCC[36]+0.3535533905932737*GBC[36]; 
  surft2_lo[19] = -(0.49755260400283263*GCC[47])+0.49755260400283263*GBC[47]+0.3535533905932737*GCC[41]+0.3535533905932737*GBC[41]; 
  surft2_up[0] = G_surf_CC_vy[0]; 
  surft2_up[1] = G_surf_CC_vy[1]; 
  surft2_up[2] = G_surf_CC_vy[2]; 
  surft2_up[3] = G_surf_CC_vy[3]; 
  surft2_up[4] = G_surf_CC_vy[4]; 
  surft2_up[5] = G_surf_CC_vy[5]; 
  surft2_up[6] = G_surf_CC_vy[6]; 
  surft2_up[7] = G_surf_CC_vy[7]; 
  surft2_up[8] = G_surf_CC_vy[8]; 
  surft2_up[9] = G_surf_CC_vy[9]; 
  surft2_up[10] = G_surf_CC_vy[10]; 
  surft2_up[11] = G_surf_CC_vy[11]; 
  surft2_up[12] = G_surf_CC_vy[12]; 
  surft2_up[13] = G_surf_CC_vy[13]; 
  surft2_up[14] = G_surf_CC_vy[14]; 
  surft2_up[15] = G_surf_CC_vy[15]; 
  surft2_up[16] = G_surf_CC_vy[16]; 
  surft2_up[17] = G_surf_CC_vy[17]; 
  surft2_up[18] = G_surf_CC_vy[18]; 
  surft2_up[19] = G_surf_CC_vy[19]; 

  vol[7] = 3.0*GCC[0]; 
  vol[15] = 3.0*GCC[1]; 
  vol[18] = 3.0*GCC[4]; 
  vol[22] = 6.7082039324993685*GCC[2]; 
  vol[24] = 6.7082039324993685*GCC[3]; 
  vol[31] = 3.0*GCC[8]; 
  vol[32] = 3.0*GCC[11]; 
  vol[33] = 6.708203932499369*GCC[5]; 
  vol[34] = 6.708203932499369*GCC[6]; 
  vol[38] = 6.708203932499369*GCC[9]; 
  vol[40] = 6.708203932499369*GCC[10]; 
  vol[43] = 3.0*GCC[14]; 
  vol[44] = 3.0*GCC[25]; 
  vol[45] = 6.7082039324993685*GCC[16]; 
  vol[46] = 6.7082039324993685*GCC[17]; 
  vol[47] = 3.0*GCC[28]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]-1.224744871391589*surft2_up[0]+1.224744871391589*surft2_lo[0]+1.224744871391589*surft1_up[0]+1.224744871391589*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]-1.224744871391589*surft2_up[1]+1.224744871391589*surft2_lo[1]+1.224744871391589*surft1_up[1]+1.224744871391589*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]+1.224744871391589*surft1_up[2]+1.224744871391589*surft1_lo[2]-2.1213203435596424*surft2_up[0]-2.1213203435596424*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]-1.224744871391589*surft2_up[3]+1.224744871391589*surft2_lo[3]+1.224744871391589*surft1_up[3]+1.224744871391589*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*surft1_up[4]+1.224744871391589*surft1_lo[4]-2.1213203435596424*surft2_up[1]-2.1213203435596424*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]-1.224744871391589*surft2_up[5]+1.224744871391589*surft2_lo[5]+1.224744871391589*surft1_up[5]+1.224744871391589*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*surft1_up[6]+1.224744871391589*surft1_lo[6]-2.1213203435596424*surft2_up[3]-2.1213203435596424*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]-1.224744871391589*surft2_up[7]+1.224744871391589*surft2_lo[7]+1.224744871391589*surft1_up[7]+1.224744871391589*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]-4.743416490252569*surft2_up[2]-4.743416490252569*surft2_lo[2]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+1.224744871391589*surft1_up[8]+1.224744871391589*surft1_lo[8]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]-1.224744871391589*surft2_up[9]+1.224744871391589*surft2_lo[9]+1.224744871391589*surft1_up[9]+1.224744871391589*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+0.7071067811865475*surft1_up[16]-0.7071067811865475*surft1_lo[16])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*surft1_up[10]+1.224744871391589*surft1_lo[10]-2.1213203435596424*surft2_up[5]-2.1213203435596424*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[32] = (vol[32]+1.224744871391589*surft1_up[11]+1.224744871391589*surft1_lo[11]-2.1213203435596424*surft2_up[7]-2.1213203435596424*surft2_lo[7])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]-4.743416490252569*surft2_up[4]-4.743416490252569*surft2_lo[4]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]+1.224744871391589*surft1_up[12]+1.224744871391589*surft1_lo[12]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]-1.224744871391589*surft2_up[13]+1.224744871391589*surft2_lo[13]+1.224744871391589*surft1_up[13]+1.224744871391589*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]-2.7386127875258306*surft2_up[10]+2.7386127875258306*surft2_lo[10]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]+0.7071067811865475*surft1_up[17]-0.7071067811865475*surft1_lo[17])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]+0.7071067811865475*surft1_up[18]-0.7071067811865475*surft1_lo[18])*dv1_sq*gamma_avg; 
  out[40] = (vol[40]+1.224744871391589*surft1_up[14]+1.224744871391589*surft1_lo[14]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[41] = (vol[41]-1.224744871391589*surft2_up[15]+1.224744871391589*surft2_lo[15]+1.224744871391589*surft1_up[15]+1.224744871391589*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[42] = (vol[42]+0.7071067811865475*surft1_up[19]-0.7071067811865475*surft1_lo[19])*dv1_sq*gamma_avg; 
  out[43] = (vol[43]+1.224744871391589*surft1_up[16]+1.224744871391589*surft1_lo[16]-2.1213203435596424*surft2_up[9]-2.1213203435596424*surft2_lo[9])*dv1_sq*gamma_avg; 
  out[44] = (vol[44]+1.224744871391589*surft1_up[17]+1.224744871391589*surft1_lo[17]-2.1213203435596424*surft2_up[13]-2.1213203435596424*surft2_lo[13])*dv1_sq*gamma_avg; 
  out[45] = (vol[45]-4.743416490252569*surft2_up[10]-4.743416490252569*surft2_lo[10]+1.5811388300841898*surft1_up[10]-1.5811388300841898*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[46] = (vol[46]+1.224744871391589*surft1_up[18]+1.224744871391589*surft1_lo[18]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[47] = (vol[47]+1.224744871391589*surft1_up[19]+1.224744871391589*surft1_lo[19]-2.1213203435596424*surft2_up[15]-2.1213203435596424*surft2_lo[15])*dv1_sq*gamma_avg; 
} 
