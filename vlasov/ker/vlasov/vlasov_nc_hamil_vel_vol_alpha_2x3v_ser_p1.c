#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_vol_alpha_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv2 = dxv[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double wx3 = w[4]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_0[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_0[1]*dv10; 
  out[2] += 1.732050807568877*hamil[1]*pt_conf_0[2]*dv10; 
  out[4] += 1.732050807568877*pt_conf_0[0]*hamil[4]*dv10; 
  out[5] += 1.732050807568877*pt_conf_0[0]*hamil[5]*dv10; 
  out[6] += 1.732050807568877*hamil[1]*pt_conf_0[3]*dv10; 
  out[9] += 1.732050807568877*pt_conf_0[1]*hamil[4]*dv10; 
  out[10] += 1.732050807568877*pt_conf_0[2]*hamil[4]*dv10; 
  out[12] += 1.732050807568877*pt_conf_0[1]*hamil[5]*dv10; 
  out[13] += 1.732050807568877*pt_conf_0[2]*hamil[5]*dv10; 
  out[15] += 1.732050807568877*pt_conf_0[0]*hamil[7]*dv10; 
  out[17] += 1.732050807568877*pt_conf_0[3]*hamil[4]*dv10; 
  out[20] += 1.732050807568877*pt_conf_0[3]*hamil[5]*dv10; 
  out[23] += 1.732050807568877*pt_conf_0[1]*hamil[7]*dv10; 
  out[24] += 1.732050807568877*pt_conf_0[2]*hamil[7]*dv10; 
  out[28] += 1.732050807568877*pt_conf_0[3]*hamil[7]*dv10; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[4]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_1[0]*hamil[2]*dv11; 
  out[1] += 1.732050807568877*pt_conf_1[1]*hamil[2]*dv11; 
  out[2] += 1.732050807568877*hamil[2]*pt_conf_1[2]*dv11; 
  out[3] += 1.732050807568877*pt_conf_1[0]*hamil[4]*dv11; 
  out[5] += 1.732050807568877*pt_conf_1[0]*hamil[6]*dv11; 
  out[6] += 1.732050807568877*hamil[2]*pt_conf_1[3]*dv11; 
  out[7] += 1.732050807568877*pt_conf_1[1]*hamil[4]*dv11; 
  out[8] += 1.732050807568877*pt_conf_1[2]*hamil[4]*dv11; 
  out[12] += 1.732050807568877*pt_conf_1[1]*hamil[6]*dv11; 
  out[13] += 1.732050807568877*pt_conf_1[2]*hamil[6]*dv11; 
  out[14] += 1.732050807568877*pt_conf_1[0]*hamil[7]*dv11; 
  out[16] += 1.732050807568877*pt_conf_1[3]*hamil[4]*dv11; 
  out[20] += 1.732050807568877*pt_conf_1[3]*hamil[6]*dv11; 
  out[21] += 1.732050807568877*pt_conf_1[1]*hamil[7]*dv11; 
  out[22] += 1.732050807568877*pt_conf_1[2]*hamil[7]*dv11; 
  out[27] += 1.732050807568877*pt_conf_1[3]*hamil[7]*dv11; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[8]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_2[0]*hamil[3]*dv12; 
  out[1] += 1.732050807568877*pt_conf_2[1]*hamil[3]*dv12; 
  out[2] += 1.732050807568877*pt_conf_2[2]*hamil[3]*dv12; 
  out[3] += 1.732050807568877*pt_conf_2[0]*hamil[5]*dv12; 
  out[4] += 1.732050807568877*pt_conf_2[0]*hamil[6]*dv12; 
  out[6] += 1.732050807568877*hamil[3]*pt_conf_2[3]*dv12; 
  out[7] += 1.732050807568877*pt_conf_2[1]*hamil[5]*dv12; 
  out[8] += 1.732050807568877*pt_conf_2[2]*hamil[5]*dv12; 
  out[9] += 1.732050807568877*pt_conf_2[1]*hamil[6]*dv12; 
  out[10] += 1.732050807568877*pt_conf_2[2]*hamil[6]*dv12; 
  out[11] += 1.732050807568877*pt_conf_2[0]*hamil[7]*dv12; 
  out[16] += 1.732050807568877*pt_conf_2[3]*hamil[5]*dv12; 
  out[17] += 1.732050807568877*pt_conf_2[3]*hamil[6]*dv12; 
  out[18] += 1.732050807568877*pt_conf_2[1]*hamil[7]*dv12; 
  out[19] += 1.732050807568877*pt_conf_2[2]*hamil[7]*dv12; 
  out[26] += 1.732050807568877*pt_conf_2[3]*hamil[7]*dv12; 
  } 
  if (dir == 4) { 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[12]; 
  if (dir == 1) { 
  out[0] += 1.732050807568877*pt_conf_3[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_3[1]*dv10; 
  out[2] += 1.732050807568877*hamil[1]*pt_conf_3[2]*dv10; 
  out[4] += 1.732050807568877*pt_conf_3[0]*hamil[4]*dv10; 
  out[5] += 1.732050807568877*pt_conf_3[0]*hamil[5]*dv10; 
  out[6] += 1.732050807568877*hamil[1]*pt_conf_3[3]*dv10; 
  out[9] += 1.732050807568877*pt_conf_3[1]*hamil[4]*dv10; 
  out[10] += 1.732050807568877*pt_conf_3[2]*hamil[4]*dv10; 
  out[12] += 1.732050807568877*pt_conf_3[1]*hamil[5]*dv10; 
  out[13] += 1.732050807568877*pt_conf_3[2]*hamil[5]*dv10; 
  out[15] += 1.732050807568877*pt_conf_3[0]*hamil[7]*dv10; 
  out[17] += 1.732050807568877*pt_conf_3[3]*hamil[4]*dv10; 
  out[20] += 1.732050807568877*pt_conf_3[3]*hamil[5]*dv10; 
  out[23] += 1.732050807568877*pt_conf_3[1]*hamil[7]*dv10; 
  out[24] += 1.732050807568877*pt_conf_3[2]*hamil[7]*dv10; 
  out[28] += 1.732050807568877*pt_conf_3[3]*hamil[7]*dv10; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[16]; 
  if (dir == 1) { 
  out[0] += 1.732050807568877*pt_conf_4[0]*hamil[2]*dv11; 
  out[1] += 1.732050807568877*pt_conf_4[1]*hamil[2]*dv11; 
  out[2] += 1.732050807568877*hamil[2]*pt_conf_4[2]*dv11; 
  out[3] += 1.732050807568877*pt_conf_4[0]*hamil[4]*dv11; 
  out[5] += 1.732050807568877*pt_conf_4[0]*hamil[6]*dv11; 
  out[6] += 1.732050807568877*hamil[2]*pt_conf_4[3]*dv11; 
  out[7] += 1.732050807568877*pt_conf_4[1]*hamil[4]*dv11; 
  out[8] += 1.732050807568877*pt_conf_4[2]*hamil[4]*dv11; 
  out[12] += 1.732050807568877*pt_conf_4[1]*hamil[6]*dv11; 
  out[13] += 1.732050807568877*pt_conf_4[2]*hamil[6]*dv11; 
  out[14] += 1.732050807568877*pt_conf_4[0]*hamil[7]*dv11; 
  out[16] += 1.732050807568877*pt_conf_4[3]*hamil[4]*dv11; 
  out[20] += 1.732050807568877*pt_conf_4[3]*hamil[6]*dv11; 
  out[21] += 1.732050807568877*pt_conf_4[1]*hamil[7]*dv11; 
  out[22] += 1.732050807568877*pt_conf_4[2]*hamil[7]*dv11; 
  out[27] += 1.732050807568877*pt_conf_4[3]*hamil[7]*dv11; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[20]; 
  if (dir == 1) { 
  out[0] += 1.732050807568877*pt_conf_5[0]*hamil[3]*dv12; 
  out[1] += 1.732050807568877*pt_conf_5[1]*hamil[3]*dv12; 
  out[2] += 1.732050807568877*pt_conf_5[2]*hamil[3]*dv12; 
  out[3] += 1.732050807568877*pt_conf_5[0]*hamil[5]*dv12; 
  out[4] += 1.732050807568877*pt_conf_5[0]*hamil[6]*dv12; 
  out[6] += 1.732050807568877*hamil[3]*pt_conf_5[3]*dv12; 
  out[7] += 1.732050807568877*pt_conf_5[1]*hamil[5]*dv12; 
  out[8] += 1.732050807568877*pt_conf_5[2]*hamil[5]*dv12; 
  out[9] += 1.732050807568877*pt_conf_5[1]*hamil[6]*dv12; 
  out[10] += 1.732050807568877*pt_conf_5[2]*hamil[6]*dv12; 
  out[11] += 1.732050807568877*pt_conf_5[0]*hamil[7]*dv12; 
  out[16] += 1.732050807568877*pt_conf_5[3]*hamil[5]*dv12; 
  out[17] += 1.732050807568877*pt_conf_5[3]*hamil[6]*dv12; 
  out[18] += 1.732050807568877*pt_conf_5[1]*hamil[7]*dv12; 
  out[19] += 1.732050807568877*pt_conf_5[2]*hamil[7]*dv12; 
  out[26] += 1.732050807568877*pt_conf_5[3]*hamil[7]*dv12; 
  } 
  if (dir == 4) { 
  } 

  const double *pt_conf_9 = &poisson_tensor_conf[36]; 
  if (dir == 2) { 
  out[0] += dv11*(1.732050807568877*pt_conf_9[0]*hamil[2]*wx1+0.5*pt_conf_9[0]*hamil[4]*dv0); 
  out[1] += dv11*(1.732050807568877*pt_conf_9[1]*hamil[2]*wx1+0.5*pt_conf_9[1]*hamil[4]*dv0); 
  out[2] += dv11*(1.732050807568877*hamil[2]*pt_conf_9[2]*wx1+0.5*pt_conf_9[2]*hamil[4]*dv0); 
  out[3] += dv11*(1.732050807568877*pt_conf_9[0]*hamil[4]*wx1+0.5*pt_conf_9[0]*hamil[2]*dv0); 
  out[5] += dv11*(1.732050807568877*pt_conf_9[0]*hamil[6]*wx1+0.5*pt_conf_9[0]*hamil[7]*dv0); 
  out[6] += dv11*(1.732050807568877*hamil[2]*pt_conf_9[3]*wx1+0.5*pt_conf_9[3]*hamil[4]*dv0); 
  out[7] += dv11*(1.732050807568877*pt_conf_9[1]*hamil[4]*wx1+0.5*pt_conf_9[1]*hamil[2]*dv0); 
  out[8] += dv11*(1.732050807568877*pt_conf_9[2]*hamil[4]*wx1+0.5*hamil[2]*pt_conf_9[2]*dv0); 
  out[12] += dv11*(1.732050807568877*pt_conf_9[1]*hamil[6]*wx1+0.5*pt_conf_9[1]*hamil[7]*dv0); 
  out[13] += dv11*(1.732050807568877*pt_conf_9[2]*hamil[6]*wx1+0.5*pt_conf_9[2]*hamil[7]*dv0); 
  out[14] += dv11*(1.732050807568877*pt_conf_9[0]*hamil[7]*wx1+0.5*pt_conf_9[0]*hamil[6]*dv0); 
  out[16] += dv11*(1.732050807568877*pt_conf_9[3]*hamil[4]*wx1+0.5*hamil[2]*pt_conf_9[3]*dv0); 
  out[20] += dv11*(1.732050807568877*pt_conf_9[3]*hamil[6]*wx1+0.5*pt_conf_9[3]*hamil[7]*dv0); 
  out[21] += dv11*(1.732050807568877*pt_conf_9[1]*hamil[7]*wx1+0.5*pt_conf_9[1]*hamil[6]*dv0); 
  out[22] += dv11*(1.732050807568877*pt_conf_9[2]*hamil[7]*wx1+0.5*pt_conf_9[2]*hamil[6]*dv0); 
  out[27] += dv11*(1.732050807568877*pt_conf_9[3]*hamil[7]*wx1+0.5*pt_conf_9[3]*hamil[6]*dv0); 
  } 
  if (dir == 3) { 
  out[0] += -1.732050807568877*pt_conf_9[0]*hamil[1]*dv11*wx1; 
  out[1] += -1.732050807568877*hamil[1]*pt_conf_9[1]*dv11*wx1; 
  out[2] += -1.732050807568877*hamil[1]*pt_conf_9[2]*dv11*wx1; 
  out[3] += -0.5*pt_conf_9[0]*hamil[1]*dv0*dv11; 
  out[4] += -1.732050807568877*pt_conf_9[0]*hamil[4]*dv11*wx1; 
  out[5] += -1.732050807568877*pt_conf_9[0]*hamil[5]*dv11*wx1; 
  out[6] += -1.732050807568877*hamil[1]*pt_conf_9[3]*dv11*wx1; 
  out[7] += -0.5*hamil[1]*pt_conf_9[1]*dv0*dv11; 
  out[8] += -0.5*hamil[1]*pt_conf_9[2]*dv0*dv11; 
  out[9] += -1.732050807568877*pt_conf_9[1]*hamil[4]*dv11*wx1; 
  out[10] += -1.732050807568877*pt_conf_9[2]*hamil[4]*dv11*wx1; 
  out[11] += -0.5*pt_conf_9[0]*hamil[4]*dv0*dv11; 
  out[12] += -1.732050807568877*pt_conf_9[1]*hamil[5]*dv11*wx1; 
  out[13] += -1.732050807568877*pt_conf_9[2]*hamil[5]*dv11*wx1; 
  out[14] += -0.5*pt_conf_9[0]*hamil[5]*dv0*dv11; 
  out[15] += -1.732050807568877*pt_conf_9[0]*hamil[7]*dv11*wx1; 
  out[16] += -0.5*hamil[1]*pt_conf_9[3]*dv0*dv11; 
  out[17] += -1.732050807568877*pt_conf_9[3]*hamil[4]*dv11*wx1; 
  out[18] += -0.5*pt_conf_9[1]*hamil[4]*dv0*dv11; 
  out[19] += -0.5*pt_conf_9[2]*hamil[4]*dv0*dv11; 
  out[20] += -1.732050807568877*pt_conf_9[3]*hamil[5]*dv11*wx1; 
  out[21] += -0.5*pt_conf_9[1]*hamil[5]*dv0*dv11; 
  out[22] += -0.5*pt_conf_9[2]*hamil[5]*dv0*dv11; 
  out[23] += -1.732050807568877*pt_conf_9[1]*hamil[7]*dv11*wx1; 
  out[24] += -1.732050807568877*pt_conf_9[2]*hamil[7]*dv11*wx1; 
  out[25] += -0.5*pt_conf_9[0]*hamil[7]*dv0*dv11; 
  out[26] += -0.5*pt_conf_9[3]*hamil[4]*dv0*dv11; 
  out[27] += -0.5*pt_conf_9[3]*hamil[5]*dv0*dv11; 
  out[28] += -1.732050807568877*pt_conf_9[3]*hamil[7]*dv11*wx1; 
  out[29] += -0.5*pt_conf_9[1]*hamil[7]*dv0*dv11; 
  out[30] += -0.5*pt_conf_9[2]*hamil[7]*dv0*dv11; 
  out[31] += -0.5*pt_conf_9[3]*hamil[7]*dv0*dv11; 
  } 

  const double *pt_conf_10 = &poisson_tensor_conf[40]; 
  if (dir == 2) { 
  out[0] += dv12*(1.732050807568877*pt_conf_10[0]*hamil[3]*wx1+0.5*pt_conf_10[0]*hamil[5]*dv0); 
  out[1] += dv12*(1.732050807568877*pt_conf_10[1]*hamil[3]*wx1+0.5*pt_conf_10[1]*hamil[5]*dv0); 
  out[2] += dv12*(1.732050807568877*pt_conf_10[2]*hamil[3]*wx1+0.5*pt_conf_10[2]*hamil[5]*dv0); 
  out[3] += dv12*(1.732050807568877*pt_conf_10[0]*hamil[5]*wx1+0.5*pt_conf_10[0]*hamil[3]*dv0); 
  out[4] += dv12*(1.732050807568877*pt_conf_10[0]*hamil[6]*wx1+0.5*pt_conf_10[0]*hamil[7]*dv0); 
  out[6] += dv12*(1.732050807568877*hamil[3]*pt_conf_10[3]*wx1+0.5*pt_conf_10[3]*hamil[5]*dv0); 
  out[7] += dv12*(1.732050807568877*pt_conf_10[1]*hamil[5]*wx1+0.5*pt_conf_10[1]*hamil[3]*dv0); 
  out[8] += dv12*(1.732050807568877*pt_conf_10[2]*hamil[5]*wx1+0.5*pt_conf_10[2]*hamil[3]*dv0); 
  out[9] += dv12*(1.732050807568877*pt_conf_10[1]*hamil[6]*wx1+0.5*pt_conf_10[1]*hamil[7]*dv0); 
  out[10] += dv12*(1.732050807568877*pt_conf_10[2]*hamil[6]*wx1+0.5*pt_conf_10[2]*hamil[7]*dv0); 
  out[11] += dv12*(1.732050807568877*pt_conf_10[0]*hamil[7]*wx1+0.5*pt_conf_10[0]*hamil[6]*dv0); 
  out[16] += dv12*(1.732050807568877*pt_conf_10[3]*hamil[5]*wx1+0.5*hamil[3]*pt_conf_10[3]*dv0); 
  out[17] += dv12*(1.732050807568877*pt_conf_10[3]*hamil[6]*wx1+0.5*pt_conf_10[3]*hamil[7]*dv0); 
  out[18] += dv12*(1.732050807568877*pt_conf_10[1]*hamil[7]*wx1+0.5*pt_conf_10[1]*hamil[6]*dv0); 
  out[19] += dv12*(1.732050807568877*pt_conf_10[2]*hamil[7]*wx1+0.5*pt_conf_10[2]*hamil[6]*dv0); 
  out[26] += dv12*(1.732050807568877*pt_conf_10[3]*hamil[7]*wx1+0.5*pt_conf_10[3]*hamil[6]*dv0); 
  } 
  if (dir == 4) { 
  out[0] += -1.732050807568877*pt_conf_10[0]*hamil[1]*dv12*wx1; 
  out[1] += -1.732050807568877*hamil[1]*pt_conf_10[1]*dv12*wx1; 
  out[2] += -1.732050807568877*hamil[1]*pt_conf_10[2]*dv12*wx1; 
  out[3] += -0.5*pt_conf_10[0]*hamil[1]*dv0*dv12; 
  out[4] += -1.732050807568877*pt_conf_10[0]*hamil[4]*dv12*wx1; 
  out[5] += -1.732050807568877*pt_conf_10[0]*hamil[5]*dv12*wx1; 
  out[6] += -1.732050807568877*hamil[1]*pt_conf_10[3]*dv12*wx1; 
  out[7] += -0.5*hamil[1]*pt_conf_10[1]*dv0*dv12; 
  out[8] += -0.5*hamil[1]*pt_conf_10[2]*dv0*dv12; 
  out[9] += -1.732050807568877*pt_conf_10[1]*hamil[4]*dv12*wx1; 
  out[10] += -1.732050807568877*pt_conf_10[2]*hamil[4]*dv12*wx1; 
  out[11] += -0.5*pt_conf_10[0]*hamil[4]*dv0*dv12; 
  out[12] += -1.732050807568877*pt_conf_10[1]*hamil[5]*dv12*wx1; 
  out[13] += -1.732050807568877*pt_conf_10[2]*hamil[5]*dv12*wx1; 
  out[14] += -0.5*pt_conf_10[0]*hamil[5]*dv0*dv12; 
  out[15] += -1.732050807568877*pt_conf_10[0]*hamil[7]*dv12*wx1; 
  out[16] += -0.5*hamil[1]*pt_conf_10[3]*dv0*dv12; 
  out[17] += -1.732050807568877*pt_conf_10[3]*hamil[4]*dv12*wx1; 
  out[18] += -0.5*pt_conf_10[1]*hamil[4]*dv0*dv12; 
  out[19] += -0.5*pt_conf_10[2]*hamil[4]*dv0*dv12; 
  out[20] += -1.732050807568877*pt_conf_10[3]*hamil[5]*dv12*wx1; 
  out[21] += -0.5*pt_conf_10[1]*hamil[5]*dv0*dv12; 
  out[22] += -0.5*pt_conf_10[2]*hamil[5]*dv0*dv12; 
  out[23] += -1.732050807568877*pt_conf_10[1]*hamil[7]*dv12*wx1; 
  out[24] += -1.732050807568877*pt_conf_10[2]*hamil[7]*dv12*wx1; 
  out[25] += -0.5*pt_conf_10[0]*hamil[7]*dv0*dv12; 
  out[26] += -0.5*pt_conf_10[3]*hamil[4]*dv0*dv12; 
  out[27] += -0.5*pt_conf_10[3]*hamil[5]*dv0*dv12; 
  out[28] += -1.732050807568877*pt_conf_10[3]*hamil[7]*dv12*wx1; 
  out[29] += -0.5*pt_conf_10[1]*hamil[7]*dv0*dv12; 
  out[30] += -0.5*pt_conf_10[2]*hamil[7]*dv0*dv12; 
  out[31] += -0.5*pt_conf_10[3]*hamil[7]*dv0*dv12; 
  } 

  const double *pt_conf_11 = &poisson_tensor_conf[44]; 
  if (dir == 3) { 
  out[0] += dv12*(1.732050807568877*pt_conf_11[0]*hamil[3]*wx1+0.5*pt_conf_11[0]*hamil[5]*dv0); 
  out[1] += dv12*(1.732050807568877*pt_conf_11[1]*hamil[3]*wx1+0.5*pt_conf_11[1]*hamil[5]*dv0); 
  out[2] += dv12*(1.732050807568877*pt_conf_11[2]*hamil[3]*wx1+0.5*pt_conf_11[2]*hamil[5]*dv0); 
  out[3] += dv12*(1.732050807568877*pt_conf_11[0]*hamil[5]*wx1+0.5*pt_conf_11[0]*hamil[3]*dv0); 
  out[4] += dv12*(1.732050807568877*pt_conf_11[0]*hamil[6]*wx1+0.5*pt_conf_11[0]*hamil[7]*dv0); 
  out[6] += dv12*(1.732050807568877*hamil[3]*pt_conf_11[3]*wx1+0.5*pt_conf_11[3]*hamil[5]*dv0); 
  out[7] += dv12*(1.732050807568877*pt_conf_11[1]*hamil[5]*wx1+0.5*pt_conf_11[1]*hamil[3]*dv0); 
  out[8] += dv12*(1.732050807568877*pt_conf_11[2]*hamil[5]*wx1+0.5*pt_conf_11[2]*hamil[3]*dv0); 
  out[9] += dv12*(1.732050807568877*pt_conf_11[1]*hamil[6]*wx1+0.5*pt_conf_11[1]*hamil[7]*dv0); 
  out[10] += dv12*(1.732050807568877*pt_conf_11[2]*hamil[6]*wx1+0.5*pt_conf_11[2]*hamil[7]*dv0); 
  out[11] += dv12*(1.732050807568877*pt_conf_11[0]*hamil[7]*wx1+0.5*pt_conf_11[0]*hamil[6]*dv0); 
  out[16] += dv12*(1.732050807568877*pt_conf_11[3]*hamil[5]*wx1+0.5*hamil[3]*pt_conf_11[3]*dv0); 
  out[17] += dv12*(1.732050807568877*pt_conf_11[3]*hamil[6]*wx1+0.5*pt_conf_11[3]*hamil[7]*dv0); 
  out[18] += dv12*(1.732050807568877*pt_conf_11[1]*hamil[7]*wx1+0.5*pt_conf_11[1]*hamil[6]*dv0); 
  out[19] += dv12*(1.732050807568877*pt_conf_11[2]*hamil[7]*wx1+0.5*pt_conf_11[2]*hamil[6]*dv0); 
  out[26] += dv12*(1.732050807568877*pt_conf_11[3]*hamil[7]*wx1+0.5*pt_conf_11[3]*hamil[6]*dv0); 
  } 
  if (dir == 4) { 
  out[0] += dv12*((-1.732050807568877*pt_conf_11[0]*hamil[2]*wx1)-0.5*pt_conf_11[0]*hamil[4]*dv0); 
  out[1] += dv12*((-1.732050807568877*pt_conf_11[1]*hamil[2]*wx1)-0.5*pt_conf_11[1]*hamil[4]*dv0); 
  out[2] += dv12*((-1.732050807568877*hamil[2]*pt_conf_11[2]*wx1)-0.5*pt_conf_11[2]*hamil[4]*dv0); 
  out[3] += dv12*((-1.732050807568877*pt_conf_11[0]*hamil[4]*wx1)-0.5*pt_conf_11[0]*hamil[2]*dv0); 
  out[5] += dv12*((-1.732050807568877*pt_conf_11[0]*hamil[6]*wx1)-0.5*pt_conf_11[0]*hamil[7]*dv0); 
  out[6] += dv12*((-1.732050807568877*hamil[2]*pt_conf_11[3]*wx1)-0.5*pt_conf_11[3]*hamil[4]*dv0); 
  out[7] += dv12*((-1.732050807568877*pt_conf_11[1]*hamil[4]*wx1)-0.5*pt_conf_11[1]*hamil[2]*dv0); 
  out[8] += dv12*((-1.732050807568877*pt_conf_11[2]*hamil[4]*wx1)-0.5*hamil[2]*pt_conf_11[2]*dv0); 
  out[12] += dv12*((-1.732050807568877*pt_conf_11[1]*hamil[6]*wx1)-0.5*pt_conf_11[1]*hamil[7]*dv0); 
  out[13] += dv12*((-1.732050807568877*pt_conf_11[2]*hamil[6]*wx1)-0.5*pt_conf_11[2]*hamil[7]*dv0); 
  out[14] += dv12*((-1.732050807568877*pt_conf_11[0]*hamil[7]*wx1)-0.5*pt_conf_11[0]*hamil[6]*dv0); 
  out[16] += dv12*((-1.732050807568877*pt_conf_11[3]*hamil[4]*wx1)-0.5*hamil[2]*pt_conf_11[3]*dv0); 
  out[20] += dv12*((-1.732050807568877*pt_conf_11[3]*hamil[6]*wx1)-0.5*pt_conf_11[3]*hamil[7]*dv0); 
  out[21] += dv12*((-1.732050807568877*pt_conf_11[1]*hamil[7]*wx1)-0.5*pt_conf_11[1]*hamil[6]*dv0); 
  out[22] += dv12*((-1.732050807568877*pt_conf_11[2]*hamil[7]*wx1)-0.5*pt_conf_11[2]*hamil[6]*dv0); 
  out[27] += dv12*((-1.732050807568877*pt_conf_11[3]*hamil[7]*wx1)-0.5*pt_conf_11[3]*hamil[6]*dv0); 
  } 

  const double *pt_conf_12 = &poisson_tensor_conf[48]; 
  if (dir == 2) { 
  out[0] += 1.732050807568877*pt_conf_12[0]*hamil[2]*dv11*wx2; 
  out[1] += 1.732050807568877*pt_conf_12[1]*hamil[2]*dv11*wx2; 
  out[2] += 1.732050807568877*hamil[2]*pt_conf_12[2]*dv11*wx2; 
  out[3] += 1.732050807568877*pt_conf_12[0]*hamil[4]*dv11*wx2; 
  out[4] += 0.5*pt_conf_12[0]*hamil[2]*dv1*dv11; 
  out[5] += 1.732050807568877*pt_conf_12[0]*hamil[6]*dv11*wx2; 
  out[6] += 1.732050807568877*hamil[2]*pt_conf_12[3]*dv11*wx2; 
  out[7] += 1.732050807568877*pt_conf_12[1]*hamil[4]*dv11*wx2; 
  out[8] += 1.732050807568877*pt_conf_12[2]*hamil[4]*dv11*wx2; 
  out[9] += 0.5*pt_conf_12[1]*hamil[2]*dv1*dv11; 
  out[10] += 0.5*hamil[2]*pt_conf_12[2]*dv1*dv11; 
  out[11] += 0.5*pt_conf_12[0]*hamil[4]*dv1*dv11; 
  out[12] += 1.732050807568877*pt_conf_12[1]*hamil[6]*dv11*wx2; 
  out[13] += 1.732050807568877*pt_conf_12[2]*hamil[6]*dv11*wx2; 
  out[14] += 1.732050807568877*pt_conf_12[0]*hamil[7]*dv11*wx2; 
  out[15] += 0.5*pt_conf_12[0]*hamil[6]*dv1*dv11; 
  out[16] += 1.732050807568877*pt_conf_12[3]*hamil[4]*dv11*wx2; 
  out[17] += 0.5*hamil[2]*pt_conf_12[3]*dv1*dv11; 
  out[18] += 0.5*pt_conf_12[1]*hamil[4]*dv1*dv11; 
  out[19] += 0.5*pt_conf_12[2]*hamil[4]*dv1*dv11; 
  out[20] += 1.732050807568877*pt_conf_12[3]*hamil[6]*dv11*wx2; 
  out[21] += 1.732050807568877*pt_conf_12[1]*hamil[7]*dv11*wx2; 
  out[22] += 1.732050807568877*pt_conf_12[2]*hamil[7]*dv11*wx2; 
  out[23] += 0.5*pt_conf_12[1]*hamil[6]*dv1*dv11; 
  out[24] += 0.5*pt_conf_12[2]*hamil[6]*dv1*dv11; 
  out[25] += 0.5*pt_conf_12[0]*hamil[7]*dv1*dv11; 
  out[26] += 0.5*pt_conf_12[3]*hamil[4]*dv1*dv11; 
  out[27] += 1.732050807568877*pt_conf_12[3]*hamil[7]*dv11*wx2; 
  out[28] += 0.5*pt_conf_12[3]*hamil[6]*dv1*dv11; 
  out[29] += 0.5*pt_conf_12[1]*hamil[7]*dv1*dv11; 
  out[30] += 0.5*pt_conf_12[2]*hamil[7]*dv1*dv11; 
  out[31] += 0.5*pt_conf_12[3]*hamil[7]*dv1*dv11; 
  } 
  if (dir == 3) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_12[0]*hamil[1]*wx2)-0.5*pt_conf_12[0]*hamil[4]*dv1); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_12[1]*wx2)-0.5*pt_conf_12[1]*hamil[4]*dv1); 
  out[2] += dv11*((-1.732050807568877*hamil[1]*pt_conf_12[2]*wx2)-0.5*pt_conf_12[2]*hamil[4]*dv1); 
  out[4] += dv11*((-1.732050807568877*pt_conf_12[0]*hamil[4]*wx2)-0.5*pt_conf_12[0]*hamil[1]*dv1); 
  out[5] += dv11*((-1.732050807568877*pt_conf_12[0]*hamil[5]*wx2)-0.5*pt_conf_12[0]*hamil[7]*dv1); 
  out[6] += dv11*((-1.732050807568877*hamil[1]*pt_conf_12[3]*wx2)-0.5*pt_conf_12[3]*hamil[4]*dv1); 
  out[9] += dv11*((-1.732050807568877*pt_conf_12[1]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_12[1]*dv1); 
  out[10] += dv11*((-1.732050807568877*pt_conf_12[2]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_12[2]*dv1); 
  out[12] += dv11*((-1.732050807568877*pt_conf_12[1]*hamil[5]*wx2)-0.5*pt_conf_12[1]*hamil[7]*dv1); 
  out[13] += dv11*((-1.732050807568877*pt_conf_12[2]*hamil[5]*wx2)-0.5*pt_conf_12[2]*hamil[7]*dv1); 
  out[15] += dv11*((-1.732050807568877*pt_conf_12[0]*hamil[7]*wx2)-0.5*pt_conf_12[0]*hamil[5]*dv1); 
  out[17] += dv11*((-1.732050807568877*pt_conf_12[3]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_12[3]*dv1); 
  out[20] += dv11*((-1.732050807568877*pt_conf_12[3]*hamil[5]*wx2)-0.5*pt_conf_12[3]*hamil[7]*dv1); 
  out[23] += dv11*((-1.732050807568877*pt_conf_12[1]*hamil[7]*wx2)-0.5*pt_conf_12[1]*hamil[5]*dv1); 
  out[24] += dv11*((-1.732050807568877*pt_conf_12[2]*hamil[7]*wx2)-0.5*pt_conf_12[2]*hamil[5]*dv1); 
  out[28] += dv11*((-1.732050807568877*pt_conf_12[3]*hamil[7]*wx2)-0.5*pt_conf_12[3]*hamil[5]*dv1); 
  } 

  const double *pt_conf_13 = &poisson_tensor_conf[52]; 
  if (dir == 2) { 
  out[0] += dv12*(1.732050807568877*pt_conf_13[0]*hamil[3]*wx2+0.5*pt_conf_13[0]*hamil[6]*dv1); 
  out[1] += dv12*(1.732050807568877*pt_conf_13[1]*hamil[3]*wx2+0.5*pt_conf_13[1]*hamil[6]*dv1); 
  out[2] += dv12*(1.732050807568877*pt_conf_13[2]*hamil[3]*wx2+0.5*pt_conf_13[2]*hamil[6]*dv1); 
  out[3] += dv12*(1.732050807568877*pt_conf_13[0]*hamil[5]*wx2+0.5*pt_conf_13[0]*hamil[7]*dv1); 
  out[4] += dv12*(1.732050807568877*pt_conf_13[0]*hamil[6]*wx2+0.5*pt_conf_13[0]*hamil[3]*dv1); 
  out[6] += dv12*(1.732050807568877*hamil[3]*pt_conf_13[3]*wx2+0.5*pt_conf_13[3]*hamil[6]*dv1); 
  out[7] += dv12*(1.732050807568877*pt_conf_13[1]*hamil[5]*wx2+0.5*pt_conf_13[1]*hamil[7]*dv1); 
  out[8] += dv12*(1.732050807568877*pt_conf_13[2]*hamil[5]*wx2+0.5*pt_conf_13[2]*hamil[7]*dv1); 
  out[9] += dv12*(1.732050807568877*pt_conf_13[1]*hamil[6]*wx2+0.5*pt_conf_13[1]*hamil[3]*dv1); 
  out[10] += dv12*(1.732050807568877*pt_conf_13[2]*hamil[6]*wx2+0.5*pt_conf_13[2]*hamil[3]*dv1); 
  out[11] += dv12*(1.732050807568877*pt_conf_13[0]*hamil[7]*wx2+0.5*pt_conf_13[0]*hamil[5]*dv1); 
  out[16] += dv12*(1.732050807568877*pt_conf_13[3]*hamil[5]*wx2+0.5*pt_conf_13[3]*hamil[7]*dv1); 
  out[17] += dv12*(1.732050807568877*pt_conf_13[3]*hamil[6]*wx2+0.5*hamil[3]*pt_conf_13[3]*dv1); 
  out[18] += dv12*(1.732050807568877*pt_conf_13[1]*hamil[7]*wx2+0.5*pt_conf_13[1]*hamil[5]*dv1); 
  out[19] += dv12*(1.732050807568877*pt_conf_13[2]*hamil[7]*wx2+0.5*pt_conf_13[2]*hamil[5]*dv1); 
  out[26] += dv12*(1.732050807568877*pt_conf_13[3]*hamil[7]*wx2+0.5*pt_conf_13[3]*hamil[5]*dv1); 
  } 
  if (dir == 4) { 
  out[0] += dv12*((-1.732050807568877*pt_conf_13[0]*hamil[1]*wx2)-0.5*pt_conf_13[0]*hamil[4]*dv1); 
  out[1] += dv12*((-1.732050807568877*hamil[1]*pt_conf_13[1]*wx2)-0.5*pt_conf_13[1]*hamil[4]*dv1); 
  out[2] += dv12*((-1.732050807568877*hamil[1]*pt_conf_13[2]*wx2)-0.5*pt_conf_13[2]*hamil[4]*dv1); 
  out[4] += dv12*((-1.732050807568877*pt_conf_13[0]*hamil[4]*wx2)-0.5*pt_conf_13[0]*hamil[1]*dv1); 
  out[5] += dv12*((-1.732050807568877*pt_conf_13[0]*hamil[5]*wx2)-0.5*pt_conf_13[0]*hamil[7]*dv1); 
  out[6] += dv12*((-1.732050807568877*hamil[1]*pt_conf_13[3]*wx2)-0.5*pt_conf_13[3]*hamil[4]*dv1); 
  out[9] += dv12*((-1.732050807568877*pt_conf_13[1]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_13[1]*dv1); 
  out[10] += dv12*((-1.732050807568877*pt_conf_13[2]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_13[2]*dv1); 
  out[12] += dv12*((-1.732050807568877*pt_conf_13[1]*hamil[5]*wx2)-0.5*pt_conf_13[1]*hamil[7]*dv1); 
  out[13] += dv12*((-1.732050807568877*pt_conf_13[2]*hamil[5]*wx2)-0.5*pt_conf_13[2]*hamil[7]*dv1); 
  out[15] += dv12*((-1.732050807568877*pt_conf_13[0]*hamil[7]*wx2)-0.5*pt_conf_13[0]*hamil[5]*dv1); 
  out[17] += dv12*((-1.732050807568877*pt_conf_13[3]*hamil[4]*wx2)-0.5*hamil[1]*pt_conf_13[3]*dv1); 
  out[20] += dv12*((-1.732050807568877*pt_conf_13[3]*hamil[5]*wx2)-0.5*pt_conf_13[3]*hamil[7]*dv1); 
  out[23] += dv12*((-1.732050807568877*pt_conf_13[1]*hamil[7]*wx2)-0.5*pt_conf_13[1]*hamil[5]*dv1); 
  out[24] += dv12*((-1.732050807568877*pt_conf_13[2]*hamil[7]*wx2)-0.5*pt_conf_13[2]*hamil[5]*dv1); 
  out[28] += dv12*((-1.732050807568877*pt_conf_13[3]*hamil[7]*wx2)-0.5*pt_conf_13[3]*hamil[5]*dv1); 
  } 

  const double *pt_conf_14 = &poisson_tensor_conf[56]; 
  if (dir == 3) { 
  out[0] += dv12*(1.732050807568877*pt_conf_14[0]*hamil[3]*wx2+0.5*pt_conf_14[0]*hamil[6]*dv1); 
  out[1] += dv12*(1.732050807568877*pt_conf_14[1]*hamil[3]*wx2+0.5*pt_conf_14[1]*hamil[6]*dv1); 
  out[2] += dv12*(1.732050807568877*pt_conf_14[2]*hamil[3]*wx2+0.5*pt_conf_14[2]*hamil[6]*dv1); 
  out[3] += dv12*(1.732050807568877*pt_conf_14[0]*hamil[5]*wx2+0.5*pt_conf_14[0]*hamil[7]*dv1); 
  out[4] += dv12*(1.732050807568877*pt_conf_14[0]*hamil[6]*wx2+0.5*pt_conf_14[0]*hamil[3]*dv1); 
  out[6] += dv12*(1.732050807568877*hamil[3]*pt_conf_14[3]*wx2+0.5*pt_conf_14[3]*hamil[6]*dv1); 
  out[7] += dv12*(1.732050807568877*pt_conf_14[1]*hamil[5]*wx2+0.5*pt_conf_14[1]*hamil[7]*dv1); 
  out[8] += dv12*(1.732050807568877*pt_conf_14[2]*hamil[5]*wx2+0.5*pt_conf_14[2]*hamil[7]*dv1); 
  out[9] += dv12*(1.732050807568877*pt_conf_14[1]*hamil[6]*wx2+0.5*pt_conf_14[1]*hamil[3]*dv1); 
  out[10] += dv12*(1.732050807568877*pt_conf_14[2]*hamil[6]*wx2+0.5*pt_conf_14[2]*hamil[3]*dv1); 
  out[11] += dv12*(1.732050807568877*pt_conf_14[0]*hamil[7]*wx2+0.5*pt_conf_14[0]*hamil[5]*dv1); 
  out[16] += dv12*(1.732050807568877*pt_conf_14[3]*hamil[5]*wx2+0.5*pt_conf_14[3]*hamil[7]*dv1); 
  out[17] += dv12*(1.732050807568877*pt_conf_14[3]*hamil[6]*wx2+0.5*hamil[3]*pt_conf_14[3]*dv1); 
  out[18] += dv12*(1.732050807568877*pt_conf_14[1]*hamil[7]*wx2+0.5*pt_conf_14[1]*hamil[5]*dv1); 
  out[19] += dv12*(1.732050807568877*pt_conf_14[2]*hamil[7]*wx2+0.5*pt_conf_14[2]*hamil[5]*dv1); 
  out[26] += dv12*(1.732050807568877*pt_conf_14[3]*hamil[7]*wx2+0.5*pt_conf_14[3]*hamil[5]*dv1); 
  } 
  if (dir == 4) { 
  out[0] += -1.732050807568877*pt_conf_14[0]*hamil[2]*dv12*wx2; 
  out[1] += -1.732050807568877*pt_conf_14[1]*hamil[2]*dv12*wx2; 
  out[2] += -1.732050807568877*hamil[2]*pt_conf_14[2]*dv12*wx2; 
  out[3] += -1.732050807568877*pt_conf_14[0]*hamil[4]*dv12*wx2; 
  out[4] += -0.5*pt_conf_14[0]*hamil[2]*dv1*dv12; 
  out[5] += -1.732050807568877*pt_conf_14[0]*hamil[6]*dv12*wx2; 
  out[6] += -1.732050807568877*hamil[2]*pt_conf_14[3]*dv12*wx2; 
  out[7] += -1.732050807568877*pt_conf_14[1]*hamil[4]*dv12*wx2; 
  out[8] += -1.732050807568877*pt_conf_14[2]*hamil[4]*dv12*wx2; 
  out[9] += -0.5*pt_conf_14[1]*hamil[2]*dv1*dv12; 
  out[10] += -0.5*hamil[2]*pt_conf_14[2]*dv1*dv12; 
  out[11] += -0.5*pt_conf_14[0]*hamil[4]*dv1*dv12; 
  out[12] += -1.732050807568877*pt_conf_14[1]*hamil[6]*dv12*wx2; 
  out[13] += -1.732050807568877*pt_conf_14[2]*hamil[6]*dv12*wx2; 
  out[14] += -1.732050807568877*pt_conf_14[0]*hamil[7]*dv12*wx2; 
  out[15] += -0.5*pt_conf_14[0]*hamil[6]*dv1*dv12; 
  out[16] += -1.732050807568877*pt_conf_14[3]*hamil[4]*dv12*wx2; 
  out[17] += -0.5*hamil[2]*pt_conf_14[3]*dv1*dv12; 
  out[18] += -0.5*pt_conf_14[1]*hamil[4]*dv1*dv12; 
  out[19] += -0.5*pt_conf_14[2]*hamil[4]*dv1*dv12; 
  out[20] += -1.732050807568877*pt_conf_14[3]*hamil[6]*dv12*wx2; 
  out[21] += -1.732050807568877*pt_conf_14[1]*hamil[7]*dv12*wx2; 
  out[22] += -1.732050807568877*pt_conf_14[2]*hamil[7]*dv12*wx2; 
  out[23] += -0.5*pt_conf_14[1]*hamil[6]*dv1*dv12; 
  out[24] += -0.5*pt_conf_14[2]*hamil[6]*dv1*dv12; 
  out[25] += -0.5*pt_conf_14[0]*hamil[7]*dv1*dv12; 
  out[26] += -0.5*pt_conf_14[3]*hamil[4]*dv1*dv12; 
  out[27] += -1.732050807568877*pt_conf_14[3]*hamil[7]*dv12*wx2; 
  out[28] += -0.5*pt_conf_14[3]*hamil[6]*dv1*dv12; 
  out[29] += -0.5*pt_conf_14[1]*hamil[7]*dv1*dv12; 
  out[30] += -0.5*pt_conf_14[2]*hamil[7]*dv1*dv12; 
  out[31] += -0.5*pt_conf_14[3]*hamil[7]*dv1*dv12; 
  } 

  const double *pt_conf_15 = &poisson_tensor_conf[60]; 
  if (dir == 2) { 
  out[0] += dv11*(1.732050807568877*pt_conf_15[0]*hamil[2]*wx3+0.5*pt_conf_15[0]*hamil[6]*dv2); 
  out[1] += dv11*(1.732050807568877*pt_conf_15[1]*hamil[2]*wx3+0.5*pt_conf_15[1]*hamil[6]*dv2); 
  out[2] += dv11*(1.732050807568877*hamil[2]*pt_conf_15[2]*wx3+0.5*pt_conf_15[2]*hamil[6]*dv2); 
  out[3] += dv11*(1.732050807568877*pt_conf_15[0]*hamil[4]*wx3+0.5*pt_conf_15[0]*hamil[7]*dv2); 
  out[5] += dv11*(1.732050807568877*pt_conf_15[0]*hamil[6]*wx3+0.5*pt_conf_15[0]*hamil[2]*dv2); 
  out[6] += dv11*(1.732050807568877*hamil[2]*pt_conf_15[3]*wx3+0.5*pt_conf_15[3]*hamil[6]*dv2); 
  out[7] += dv11*(1.732050807568877*pt_conf_15[1]*hamil[4]*wx3+0.5*pt_conf_15[1]*hamil[7]*dv2); 
  out[8] += dv11*(1.732050807568877*pt_conf_15[2]*hamil[4]*wx3+0.5*pt_conf_15[2]*hamil[7]*dv2); 
  out[12] += dv11*(1.732050807568877*pt_conf_15[1]*hamil[6]*wx3+0.5*pt_conf_15[1]*hamil[2]*dv2); 
  out[13] += dv11*(1.732050807568877*pt_conf_15[2]*hamil[6]*wx3+0.5*hamil[2]*pt_conf_15[2]*dv2); 
  out[14] += dv11*(1.732050807568877*pt_conf_15[0]*hamil[7]*wx3+0.5*pt_conf_15[0]*hamil[4]*dv2); 
  out[16] += dv11*(1.732050807568877*pt_conf_15[3]*hamil[4]*wx3+0.5*pt_conf_15[3]*hamil[7]*dv2); 
  out[20] += dv11*(1.732050807568877*pt_conf_15[3]*hamil[6]*wx3+0.5*hamil[2]*pt_conf_15[3]*dv2); 
  out[21] += dv11*(1.732050807568877*pt_conf_15[1]*hamil[7]*wx3+0.5*pt_conf_15[1]*hamil[4]*dv2); 
  out[22] += dv11*(1.732050807568877*pt_conf_15[2]*hamil[7]*wx3+0.5*pt_conf_15[2]*hamil[4]*dv2); 
  out[27] += dv11*(1.732050807568877*pt_conf_15[3]*hamil[7]*wx3+0.5*pt_conf_15[3]*hamil[4]*dv2); 
  } 
  if (dir == 3) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_15[0]*hamil[1]*wx3)-0.5*pt_conf_15[0]*hamil[5]*dv2); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_15[1]*wx3)-0.5*pt_conf_15[1]*hamil[5]*dv2); 
  out[2] += dv11*((-1.732050807568877*hamil[1]*pt_conf_15[2]*wx3)-0.5*pt_conf_15[2]*hamil[5]*dv2); 
  out[4] += dv11*((-1.732050807568877*pt_conf_15[0]*hamil[4]*wx3)-0.5*pt_conf_15[0]*hamil[7]*dv2); 
  out[5] += dv11*((-1.732050807568877*pt_conf_15[0]*hamil[5]*wx3)-0.5*pt_conf_15[0]*hamil[1]*dv2); 
  out[6] += dv11*((-1.732050807568877*hamil[1]*pt_conf_15[3]*wx3)-0.5*pt_conf_15[3]*hamil[5]*dv2); 
  out[9] += dv11*((-1.732050807568877*pt_conf_15[1]*hamil[4]*wx3)-0.5*pt_conf_15[1]*hamil[7]*dv2); 
  out[10] += dv11*((-1.732050807568877*pt_conf_15[2]*hamil[4]*wx3)-0.5*pt_conf_15[2]*hamil[7]*dv2); 
  out[12] += dv11*((-1.732050807568877*pt_conf_15[1]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_15[1]*dv2); 
  out[13] += dv11*((-1.732050807568877*pt_conf_15[2]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_15[2]*dv2); 
  out[15] += dv11*((-1.732050807568877*pt_conf_15[0]*hamil[7]*wx3)-0.5*pt_conf_15[0]*hamil[4]*dv2); 
  out[17] += dv11*((-1.732050807568877*pt_conf_15[3]*hamil[4]*wx3)-0.5*pt_conf_15[3]*hamil[7]*dv2); 
  out[20] += dv11*((-1.732050807568877*pt_conf_15[3]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_15[3]*dv2); 
  out[23] += dv11*((-1.732050807568877*pt_conf_15[1]*hamil[7]*wx3)-0.5*pt_conf_15[1]*hamil[4]*dv2); 
  out[24] += dv11*((-1.732050807568877*pt_conf_15[2]*hamil[7]*wx3)-0.5*pt_conf_15[2]*hamil[4]*dv2); 
  out[28] += dv11*((-1.732050807568877*pt_conf_15[3]*hamil[7]*wx3)-0.5*pt_conf_15[3]*hamil[4]*dv2); 
  } 

  const double *pt_conf_16 = &poisson_tensor_conf[64]; 
  if (dir == 2) { 
  out[0] += 1.732050807568877*pt_conf_16[0]*hamil[3]*dv12*wx3; 
  out[1] += 1.732050807568877*pt_conf_16[1]*hamil[3]*dv12*wx3; 
  out[2] += 1.732050807568877*pt_conf_16[2]*hamil[3]*dv12*wx3; 
  out[3] += 1.732050807568877*pt_conf_16[0]*hamil[5]*dv12*wx3; 
  out[4] += 1.732050807568877*pt_conf_16[0]*hamil[6]*dv12*wx3; 
  out[5] += 0.5*pt_conf_16[0]*hamil[3]*dv12*dv2; 
  out[6] += 1.732050807568877*hamil[3]*pt_conf_16[3]*dv12*wx3; 
  out[7] += 1.732050807568877*pt_conf_16[1]*hamil[5]*dv12*wx3; 
  out[8] += 1.732050807568877*pt_conf_16[2]*hamil[5]*dv12*wx3; 
  out[9] += 1.732050807568877*pt_conf_16[1]*hamil[6]*dv12*wx3; 
  out[10] += 1.732050807568877*pt_conf_16[2]*hamil[6]*dv12*wx3; 
  out[11] += 1.732050807568877*pt_conf_16[0]*hamil[7]*dv12*wx3; 
  out[12] += 0.5*pt_conf_16[1]*hamil[3]*dv12*dv2; 
  out[13] += 0.5*pt_conf_16[2]*hamil[3]*dv12*dv2; 
  out[14] += 0.5*pt_conf_16[0]*hamil[5]*dv12*dv2; 
  out[15] += 0.5*pt_conf_16[0]*hamil[6]*dv12*dv2; 
  out[16] += 1.732050807568877*pt_conf_16[3]*hamil[5]*dv12*wx3; 
  out[17] += 1.732050807568877*pt_conf_16[3]*hamil[6]*dv12*wx3; 
  out[18] += 1.732050807568877*pt_conf_16[1]*hamil[7]*dv12*wx3; 
  out[19] += 1.732050807568877*pt_conf_16[2]*hamil[7]*dv12*wx3; 
  out[20] += 0.5*hamil[3]*pt_conf_16[3]*dv12*dv2; 
  out[21] += 0.5*pt_conf_16[1]*hamil[5]*dv12*dv2; 
  out[22] += 0.5*pt_conf_16[2]*hamil[5]*dv12*dv2; 
  out[23] += 0.5*pt_conf_16[1]*hamil[6]*dv12*dv2; 
  out[24] += 0.5*pt_conf_16[2]*hamil[6]*dv12*dv2; 
  out[25] += 0.5*pt_conf_16[0]*hamil[7]*dv12*dv2; 
  out[26] += 1.732050807568877*pt_conf_16[3]*hamil[7]*dv12*wx3; 
  out[27] += 0.5*pt_conf_16[3]*hamil[5]*dv12*dv2; 
  out[28] += 0.5*pt_conf_16[3]*hamil[6]*dv12*dv2; 
  out[29] += 0.5*pt_conf_16[1]*hamil[7]*dv12*dv2; 
  out[30] += 0.5*pt_conf_16[2]*hamil[7]*dv12*dv2; 
  out[31] += 0.5*pt_conf_16[3]*hamil[7]*dv12*dv2; 
  } 
  if (dir == 4) { 
  out[0] += dv12*((-1.732050807568877*pt_conf_16[0]*hamil[1]*wx3)-0.5*pt_conf_16[0]*hamil[5]*dv2); 
  out[1] += dv12*((-1.732050807568877*hamil[1]*pt_conf_16[1]*wx3)-0.5*pt_conf_16[1]*hamil[5]*dv2); 
  out[2] += dv12*((-1.732050807568877*hamil[1]*pt_conf_16[2]*wx3)-0.5*pt_conf_16[2]*hamil[5]*dv2); 
  out[4] += dv12*((-1.732050807568877*pt_conf_16[0]*hamil[4]*wx3)-0.5*pt_conf_16[0]*hamil[7]*dv2); 
  out[5] += dv12*((-1.732050807568877*pt_conf_16[0]*hamil[5]*wx3)-0.5*pt_conf_16[0]*hamil[1]*dv2); 
  out[6] += dv12*((-1.732050807568877*hamil[1]*pt_conf_16[3]*wx3)-0.5*pt_conf_16[3]*hamil[5]*dv2); 
  out[9] += dv12*((-1.732050807568877*pt_conf_16[1]*hamil[4]*wx3)-0.5*pt_conf_16[1]*hamil[7]*dv2); 
  out[10] += dv12*((-1.732050807568877*pt_conf_16[2]*hamil[4]*wx3)-0.5*pt_conf_16[2]*hamil[7]*dv2); 
  out[12] += dv12*((-1.732050807568877*pt_conf_16[1]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_16[1]*dv2); 
  out[13] += dv12*((-1.732050807568877*pt_conf_16[2]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_16[2]*dv2); 
  out[15] += dv12*((-1.732050807568877*pt_conf_16[0]*hamil[7]*wx3)-0.5*pt_conf_16[0]*hamil[4]*dv2); 
  out[17] += dv12*((-1.732050807568877*pt_conf_16[3]*hamil[4]*wx3)-0.5*pt_conf_16[3]*hamil[7]*dv2); 
  out[20] += dv12*((-1.732050807568877*pt_conf_16[3]*hamil[5]*wx3)-0.5*hamil[1]*pt_conf_16[3]*dv2); 
  out[23] += dv12*((-1.732050807568877*pt_conf_16[1]*hamil[7]*wx3)-0.5*pt_conf_16[1]*hamil[4]*dv2); 
  out[24] += dv12*((-1.732050807568877*pt_conf_16[2]*hamil[7]*wx3)-0.5*pt_conf_16[2]*hamil[4]*dv2); 
  out[28] += dv12*((-1.732050807568877*pt_conf_16[3]*hamil[7]*wx3)-0.5*pt_conf_16[3]*hamil[4]*dv2); 
  } 

  const double *pt_conf_17 = &poisson_tensor_conf[68]; 
  if (dir == 3) { 
  out[0] += 1.732050807568877*pt_conf_17[0]*hamil[3]*dv12*wx3; 
  out[1] += 1.732050807568877*pt_conf_17[1]*hamil[3]*dv12*wx3; 
  out[2] += 1.732050807568877*pt_conf_17[2]*hamil[3]*dv12*wx3; 
  out[3] += 1.732050807568877*pt_conf_17[0]*hamil[5]*dv12*wx3; 
  out[4] += 1.732050807568877*pt_conf_17[0]*hamil[6]*dv12*wx3; 
  out[5] += 0.5*pt_conf_17[0]*hamil[3]*dv12*dv2; 
  out[6] += 1.732050807568877*hamil[3]*pt_conf_17[3]*dv12*wx3; 
  out[7] += 1.732050807568877*pt_conf_17[1]*hamil[5]*dv12*wx3; 
  out[8] += 1.732050807568877*pt_conf_17[2]*hamil[5]*dv12*wx3; 
  out[9] += 1.732050807568877*pt_conf_17[1]*hamil[6]*dv12*wx3; 
  out[10] += 1.732050807568877*pt_conf_17[2]*hamil[6]*dv12*wx3; 
  out[11] += 1.732050807568877*pt_conf_17[0]*hamil[7]*dv12*wx3; 
  out[12] += 0.5*pt_conf_17[1]*hamil[3]*dv12*dv2; 
  out[13] += 0.5*pt_conf_17[2]*hamil[3]*dv12*dv2; 
  out[14] += 0.5*pt_conf_17[0]*hamil[5]*dv12*dv2; 
  out[15] += 0.5*pt_conf_17[0]*hamil[6]*dv12*dv2; 
  out[16] += 1.732050807568877*pt_conf_17[3]*hamil[5]*dv12*wx3; 
  out[17] += 1.732050807568877*pt_conf_17[3]*hamil[6]*dv12*wx3; 
  out[18] += 1.732050807568877*pt_conf_17[1]*hamil[7]*dv12*wx3; 
  out[19] += 1.732050807568877*pt_conf_17[2]*hamil[7]*dv12*wx3; 
  out[20] += 0.5*hamil[3]*pt_conf_17[3]*dv12*dv2; 
  out[21] += 0.5*pt_conf_17[1]*hamil[5]*dv12*dv2; 
  out[22] += 0.5*pt_conf_17[2]*hamil[5]*dv12*dv2; 
  out[23] += 0.5*pt_conf_17[1]*hamil[6]*dv12*dv2; 
  out[24] += 0.5*pt_conf_17[2]*hamil[6]*dv12*dv2; 
  out[25] += 0.5*pt_conf_17[0]*hamil[7]*dv12*dv2; 
  out[26] += 1.732050807568877*pt_conf_17[3]*hamil[7]*dv12*wx3; 
  out[27] += 0.5*pt_conf_17[3]*hamil[5]*dv12*dv2; 
  out[28] += 0.5*pt_conf_17[3]*hamil[6]*dv12*dv2; 
  out[29] += 0.5*pt_conf_17[1]*hamil[7]*dv12*dv2; 
  out[30] += 0.5*pt_conf_17[2]*hamil[7]*dv12*dv2; 
  out[31] += 0.5*pt_conf_17[3]*hamil[7]*dv12*dv2; 
  } 
  if (dir == 4) { 
  out[0] += dv12*((-1.732050807568877*pt_conf_17[0]*hamil[2]*wx3)-0.5*pt_conf_17[0]*hamil[6]*dv2); 
  out[1] += dv12*((-1.732050807568877*pt_conf_17[1]*hamil[2]*wx3)-0.5*pt_conf_17[1]*hamil[6]*dv2); 
  out[2] += dv12*((-1.732050807568877*hamil[2]*pt_conf_17[2]*wx3)-0.5*pt_conf_17[2]*hamil[6]*dv2); 
  out[3] += dv12*((-1.732050807568877*pt_conf_17[0]*hamil[4]*wx3)-0.5*pt_conf_17[0]*hamil[7]*dv2); 
  out[5] += dv12*((-1.732050807568877*pt_conf_17[0]*hamil[6]*wx3)-0.5*pt_conf_17[0]*hamil[2]*dv2); 
  out[6] += dv12*((-1.732050807568877*hamil[2]*pt_conf_17[3]*wx3)-0.5*pt_conf_17[3]*hamil[6]*dv2); 
  out[7] += dv12*((-1.732050807568877*pt_conf_17[1]*hamil[4]*wx3)-0.5*pt_conf_17[1]*hamil[7]*dv2); 
  out[8] += dv12*((-1.732050807568877*pt_conf_17[2]*hamil[4]*wx3)-0.5*pt_conf_17[2]*hamil[7]*dv2); 
  out[12] += dv12*((-1.732050807568877*pt_conf_17[1]*hamil[6]*wx3)-0.5*pt_conf_17[1]*hamil[2]*dv2); 
  out[13] += dv12*((-1.732050807568877*pt_conf_17[2]*hamil[6]*wx3)-0.5*hamil[2]*pt_conf_17[2]*dv2); 
  out[14] += dv12*((-1.732050807568877*pt_conf_17[0]*hamil[7]*wx3)-0.5*pt_conf_17[0]*hamil[4]*dv2); 
  out[16] += dv12*((-1.732050807568877*pt_conf_17[3]*hamil[4]*wx3)-0.5*pt_conf_17[3]*hamil[7]*dv2); 
  out[20] += dv12*((-1.732050807568877*pt_conf_17[3]*hamil[6]*wx3)-0.5*hamil[2]*pt_conf_17[3]*dv2); 
  out[21] += dv12*((-1.732050807568877*pt_conf_17[1]*hamil[7]*wx3)-0.5*pt_conf_17[1]*hamil[4]*dv2); 
  out[22] += dv12*((-1.732050807568877*pt_conf_17[2]*hamil[7]*wx3)-0.5*pt_conf_17[2]*hamil[4]*dv2); 
  out[27] += dv12*((-1.732050807568877*pt_conf_17[3]*hamil[7]*wx3)-0.5*pt_conf_17[3]*hamil[4]*dv2); 
  } 

} 
