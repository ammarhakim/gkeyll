#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_lovz_lovy(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[3]/dxv[2]; 
 
  double surft1_lo[20] = {0.0}; 
  double surft1_up[20] = {0.0}; 
  double surft2_lo[20] = {0.0}; 
  double surft2_up[20] = {0.0}; 
  double vol[48] = {0.0}; 
  double *out = &diff_coeff[336]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vz = &fpo_g_surf_stencil[4][40]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvz_surf_CC_vy = &fpo_dgdv_surf[60]; 

  surft1_lo[0] = dGdvz_surf_CC_vy[0]/dv1; 
  surft1_lo[1] = dGdvz_surf_CC_vy[1]/dv1; 
  surft1_lo[2] = dGdvz_surf_CC_vy[2]/dv1; 
  surft1_lo[3] = dGdvz_surf_CC_vy[3]/dv1; 
  surft1_lo[4] = dGdvz_surf_CC_vy[4]/dv1; 
  surft1_lo[5] = dGdvz_surf_CC_vy[5]/dv1; 
  surft1_lo[6] = dGdvz_surf_CC_vy[6]/dv1; 
  surft1_lo[7] = dGdvz_surf_CC_vy[7]/dv1; 
  surft1_lo[8] = dGdvz_surf_CC_vy[8]/dv1; 
  surft1_lo[9] = dGdvz_surf_CC_vy[9]/dv1; 
  surft1_lo[10] = dGdvz_surf_CC_vy[10]/dv1; 
  surft1_lo[11] = dGdvz_surf_CC_vy[11]/dv1; 
  surft1_lo[12] = dGdvz_surf_CC_vy[12]/dv1; 
  surft1_lo[13] = dGdvz_surf_CC_vy[13]/dv1; 
  surft1_lo[14] = dGdvz_surf_CC_vy[14]/dv1; 
  surft1_lo[15] = dGdvz_surf_CC_vy[15]/dv1; 
  surft1_lo[16] = dGdvz_surf_CC_vy[16]/dv1; 
  surft1_lo[17] = dGdvz_surf_CC_vy[17]/dv1; 
  surft1_lo[18] = dGdvz_surf_CC_vy[18]/dv1; 
  surft1_lo[19] = dGdvz_surf_CC_vy[19]/dv1; 
  surft1_up[0] = 0.599071547271275*GCR[27]+0.599071547271275*GCC[27]-0.8617863895711042*GCR[10]+0.8617863895711042*GCC[10]+0.6123724356957944*GCR[4]+0.6123724356957944*GCC[4]; 
  surft1_up[1] = 0.599071547271275*GCR[39]+0.599071547271275*GCC[39]-0.8617863895711042*GCR[17]+0.8617863895711042*GCC[17]+0.6123724356957944*GCR[8]+0.6123724356957944*GCC[8]; 
  surft1_up[2] = 0.599071547271275*GCR[40]+0.599071547271275*GCC[40]-0.8617863895711042*GCR[18]+0.8617863895711042*GCC[18]+0.6123724356957944*GCR[9]+0.6123724356957944*GCC[9]; 
  surft1_up[3] = -(1.9270129491651047*GCR[30])+1.9270129491651047*GCC[30]+1.369306393762915*GCR[14]+1.369306393762915*GCC[14]; 
  surft1_up[4] = 0.599071547271275*GCR[46]+0.599071547271275*GCC[46]-0.8617863895711042*GCR[31]+0.8617863895711042*GCC[31]+0.6123724356957944*GCR[16]+0.6123724356957944*GCC[16]; 
  surft1_up[5] = -(1.927012949165105*GCR[42])+1.927012949165105*GCC[42]+1.369306393762915*GCR[28]+1.369306393762915*GCC[28]; 
  surft1_up[6] = -(1.927012949165105*GCR[43])+1.927012949165105*GCC[43]+1.369306393762915*GCR[29]+1.369306393762915*GCC[29]; 
  surft1_up[7] = -(0.8617863895711042*GCR[37])+0.8617863895711042*GCC[37]+0.6123724356957944*GCR[25]+0.6123724356957944*GCC[25]; 
  surft1_up[8] = -(0.8617863895711042*GCR[38])+0.8617863895711042*GCC[38]+0.6123724356957944*GCR[26]+0.6123724356957944*GCC[26]; 
  surft1_up[10] = -(1.9270129491651047*GCR[47])+1.9270129491651047*GCC[47]+1.369306393762915*GCR[41]+1.369306393762915*GCC[41]; 
  surft1_up[11] = -(0.8617863895711042*GCR[44])+0.8617863895711042*GCC[44]+0.6123724356957944*GCR[35]+0.6123724356957944*GCC[35]; 
  surft1_up[12] = -(0.8617863895711042*GCR[45])+0.8617863895711042*GCC[45]+0.6123724356957944*GCR[36]+0.6123724356957944*GCC[36]; 

  surft2_lo[0] = G_surf_CC_vz[0]; 
  surft2_lo[1] = G_surf_CC_vz[1]; 
  surft2_lo[2] = G_surf_CC_vz[2]; 
  surft2_lo[3] = G_surf_CC_vz[3]; 
  surft2_lo[4] = G_surf_CC_vz[4]; 
  surft2_lo[5] = G_surf_CC_vz[5]; 
  surft2_lo[6] = G_surf_CC_vz[6]; 
  surft2_lo[7] = G_surf_CC_vz[7]; 
  surft2_lo[8] = G_surf_CC_vz[8]; 
  surft2_lo[9] = G_surf_CC_vz[9]; 
  surft2_lo[10] = G_surf_CC_vz[10]; 
  surft2_lo[11] = G_surf_CC_vz[11]; 
  surft2_lo[12] = G_surf_CC_vz[12]; 
  surft2_lo[13] = G_surf_CC_vz[13]; 
  surft2_lo[14] = G_surf_CC_vz[14]; 
  surft2_lo[15] = G_surf_CC_vz[15]; 
  surft2_lo[16] = G_surf_CC_vz[16]; 
  surft2_lo[17] = G_surf_CC_vz[17]; 
  surft2_lo[18] = G_surf_CC_vz[18]; 
  surft2_lo[19] = G_surf_CC_vz[19]; 
  surft2_up[0] = 0.34587411908091625*GTC[14]+0.34587411908091625*GCC[14]-0.49755260400283263*GTC[4]+0.49755260400283263*GCC[4]+0.3535533905932737*GTC[0]+0.3535533905932737*GCC[0]; 
  surft2_up[1] = 0.34587411908091625*GTC[28]+0.34587411908091625*GCC[28]-0.49755260400283263*GTC[8]+0.49755260400283263*GCC[8]+0.3535533905932737*GTC[1]+0.3535533905932737*GCC[1]; 
  surft2_up[2] = 0.34587411908091625*GTC[29]+0.34587411908091625*GCC[29]-0.49755260400283263*GTC[9]+0.49755260400283263*GCC[9]+0.3535533905932737*GTC[2]+0.3535533905932737*GCC[2]; 
  surft2_up[3] = 0.34587411908091625*GTC[30]+0.34587411908091625*GCC[30]-0.49755260400283263*GTC[10]+0.49755260400283263*GCC[10]+0.3535533905932737*GTC[3]+0.3535533905932737*GCC[3]; 
  surft2_up[4] = 0.34587411908091625*GTC[41]+0.34587411908091625*GCC[41]-0.49755260400283263*GTC[16]+0.49755260400283263*GCC[16]+0.3535533905932737*GTC[5]+0.3535533905932737*GCC[5]; 
  surft2_up[5] = 0.34587411908091625*GTC[42]+0.34587411908091625*GCC[42]-0.49755260400283263*GTC[17]+0.49755260400283263*GCC[17]+0.3535533905932737*GTC[6]+0.3535533905932737*GCC[6]; 
  surft2_up[6] = 0.34587411908091625*GTC[43]+0.34587411908091625*GCC[43]-0.49755260400283263*GTC[18]+0.49755260400283263*GCC[18]+0.3535533905932737*GTC[7]+0.3535533905932737*GCC[7]; 
  surft2_up[7] = -(0.49755260400283263*GTC[25])+0.49755260400283263*GCC[25]+0.3535533905932737*GTC[11]+0.3535533905932737*GCC[11]; 
  surft2_up[8] = -(0.49755260400283263*GTC[26])+0.49755260400283263*GCC[26]+0.3535533905932737*GTC[12]+0.3535533905932737*GCC[12]; 
  surft2_up[9] = -(0.49755260400283263*GTC[27])+0.49755260400283263*GCC[27]+0.3535533905932737*GTC[13]+0.3535533905932737*GCC[13]; 
  surft2_up[10] = 0.34587411908091625*GTC[47]+0.34587411908091625*GCC[47]-0.49755260400283263*GTC[31]+0.49755260400283263*GCC[31]+0.3535533905932737*GTC[15]+0.3535533905932737*GCC[15]; 
  surft2_up[11] = -(0.49755260400283263*GTC[35])+0.49755260400283263*GCC[35]+0.3535533905932737*GTC[19]+0.3535533905932737*GCC[19]; 
  surft2_up[12] = -(0.49755260400283263*GTC[36])+0.49755260400283263*GCC[36]+0.3535533905932737*GTC[20]+0.3535533905932737*GCC[20]; 
  surft2_up[13] = -(0.49755260400283263*GTC[37])+0.49755260400283263*GCC[37]+0.3535533905932737*GTC[21]+0.3535533905932737*GCC[21]; 
  surft2_up[14] = -(0.49755260400283263*GTC[38])+0.49755260400283263*GCC[38]+0.3535533905932737*GTC[22]+0.3535533905932737*GCC[22]; 
  surft2_up[15] = -(0.49755260400283263*GTC[39])+0.49755260400283263*GCC[39]+0.3535533905932737*GTC[23]+0.3535533905932737*GCC[23]; 
  surft2_up[16] = -(0.49755260400283263*GTC[40])+0.49755260400283263*GCC[40]+0.3535533905932737*GTC[24]+0.3535533905932737*GCC[24]; 
  surft2_up[17] = -(0.49755260400283263*GTC[44])+0.49755260400283263*GCC[44]+0.3535533905932737*GTC[32]+0.3535533905932737*GCC[32]; 
  surft2_up[18] = -(0.49755260400283263*GTC[45])+0.49755260400283263*GCC[45]+0.3535533905932737*GTC[33]+0.3535533905932737*GCC[33]; 
  surft2_up[19] = -(0.49755260400283263*GTC[46])+0.49755260400283263*GCC[46]+0.3535533905932737*GTC[34]+0.3535533905932737*GCC[34]; 

  vol[10] = 3.0*GCC[0]; 
  vol[17] = 3.0*GCC[1]; 
  vol[18] = 3.0*GCC[2]; 
  vol[27] = 6.7082039324993685*GCC[3]; 
  vol[30] = 6.7082039324993685*GCC[4]; 
  vol[31] = 3.0*GCC[5]; 
  vol[37] = 3.0*GCC[11]; 
  vol[38] = 3.0*GCC[12]; 
  vol[39] = 6.708203932499369*GCC[6]; 
  vol[40] = 6.708203932499369*GCC[7]; 
  vol[42] = 6.708203932499369*GCC[8]; 
  vol[43] = 6.708203932499369*GCC[9]; 
  vol[44] = 3.0*GCC[19]; 
  vol[45] = 3.0*GCC[20]; 
  vol[46] = 6.7082039324993685*GCC[15]; 
  vol[47] = 6.7082039324993685*GCC[16]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]-1.224744871391589*surft2_up[0]+1.224744871391589*surft2_lo[0]+1.224744871391589*surft1_up[0]+1.224744871391589*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]-1.224744871391589*surft2_up[1]+1.224744871391589*surft2_lo[1]+1.224744871391589*surft1_up[1]+1.224744871391589*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]-1.224744871391589*surft2_up[2]+1.224744871391589*surft2_lo[2]+1.224744871391589*surft1_up[2]+1.224744871391589*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+1.224744871391589*surft1_up[3]+1.224744871391589*surft1_lo[3]-2.1213203435596424*surft2_up[0]-2.1213203435596424*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[15] = (vol[15]-1.224744871391589*surft2_up[4]+1.224744871391589*surft2_lo[4]+1.224744871391589*surft1_up[4]+1.224744871391589*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+1.224744871391589*surft1_up[5]+1.224744871391589*surft1_lo[5]-2.1213203435596424*surft2_up[1]-2.1213203435596424*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*surft1_up[6]+1.224744871391589*surft1_lo[6]-2.1213203435596424*surft2_up[2]-2.1213203435596424*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]-1.224744871391589*surft2_up[7]+1.224744871391589*surft2_lo[7]+1.224744871391589*surft1_up[7]+1.224744871391589*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]-1.224744871391589*surft2_up[8]+1.224744871391589*surft2_lo[8]+1.224744871391589*surft1_up[8]+1.224744871391589*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-4.743416490252569*surft2_up[3]-4.743416490252569*surft2_lo[3]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]+0.7071067811865475*surft1_up[16]-0.7071067811865475*surft1_lo[16])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]+1.224744871391589*surft1_up[9]+1.224744871391589*surft1_lo[9]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*surft1_up[10]+1.224744871391589*surft1_lo[10]-2.1213203435596424*surft2_up[4]-2.1213203435596424*surft2_lo[4])*dv1_sq*gamma_avg; 
  out[32] = (vol[32]-1.224744871391589*surft2_up[11]+1.224744871391589*surft2_lo[11]+1.224744871391589*surft1_up[11]+1.224744871391589*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]-1.224744871391589*surft2_up[12]+1.224744871391589*surft2_lo[12]+1.224744871391589*surft1_up[12]+1.224744871391589*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]-2.7386127875258306*surft2_up[10]+2.7386127875258306*surft2_lo[10]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]+0.7071067811865475*surft1_up[17]-0.7071067811865475*surft1_lo[17])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]+0.7071067811865475*surft1_up[18]-0.7071067811865475*surft1_lo[18])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]+1.224744871391589*surft1_up[13]+1.224744871391589*surft1_lo[13]-2.1213203435596424*surft2_up[7]-2.1213203435596424*surft2_lo[7])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]+1.224744871391589*surft1_up[14]+1.224744871391589*surft1_lo[14]-2.1213203435596424*surft2_up[8]-2.1213203435596424*surft2_lo[8])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]-4.743416490252569*surft2_up[5]-4.743416490252569*surft2_lo[5]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[40] = (vol[40]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[41] = (vol[41]+0.7071067811865475*surft1_up[19]-0.7071067811865475*surft1_lo[19])*dv1_sq*gamma_avg; 
  out[42] = (vol[42]+1.224744871391589*surft1_up[15]+1.224744871391589*surft1_lo[15]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[43] = (vol[43]+1.224744871391589*surft1_up[16]+1.224744871391589*surft1_lo[16]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[44] = (vol[44]+1.224744871391589*surft1_up[17]+1.224744871391589*surft1_lo[17]-2.1213203435596424*surft2_up[11]-2.1213203435596424*surft2_lo[11])*dv1_sq*gamma_avg; 
  out[45] = (vol[45]+1.224744871391589*surft1_up[18]+1.224744871391589*surft1_lo[18]-2.1213203435596424*surft2_up[12]-2.1213203435596424*surft2_lo[12])*dv1_sq*gamma_avg; 
  out[46] = (vol[46]-4.743416490252569*surft2_up[10]-4.743416490252569*surft2_lo[10]+1.5811388300841898*surft1_up[10]-1.5811388300841898*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[47] = (vol[47]+1.224744871391589*surft1_up[19]+1.224744871391589*surft1_lo[19]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4])*dv1_sq*gamma_avg; 
} 
