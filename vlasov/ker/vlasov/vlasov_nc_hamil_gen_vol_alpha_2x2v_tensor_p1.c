#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_alpha_2x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv0 = dxv[2]; 
  const double dv1 = dxv[3]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_0[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_0[1]*dv10; 
  out[2] += 1.732050807568877*hamil[1]*pt_conf_0[2]*dv10; 
  out[3] += 3.872983346207417*pt_conf_0[0]*hamil[4]*dv10; 
  out[4] += 1.732050807568877*pt_conf_0[0]*hamil[3]*dv10; 
  out[5] += 1.732050807568877*hamil[1]*pt_conf_0[3]*dv10; 
  out[6] += 3.872983346207417*pt_conf_0[1]*hamil[4]*dv10; 
  out[7] += 3.872983346207417*pt_conf_0[2]*hamil[4]*dv10; 
  out[8] += 1.732050807568877*pt_conf_0[1]*hamil[3]*dv10; 
  out[9] += 1.732050807568877*pt_conf_0[2]*hamil[3]*dv10; 
  out[10] += 3.872983346207417*pt_conf_0[0]*hamil[6]*dv10; 
  out[12] += 1.732050807568877*pt_conf_0[0]*hamil[7]*dv10; 
  out[13] += 3.872983346207417*pt_conf_0[3]*hamil[4]*dv10; 
  out[14] += 1.732050807568877*hamil[3]*pt_conf_0[3]*dv10; 
  out[15] += 3.872983346207417*pt_conf_0[1]*hamil[6]*dv10; 
  out[16] += 3.872983346207417*pt_conf_0[2]*hamil[6]*dv10; 
  out[20] += 1.732050807568877*pt_conf_0[1]*hamil[7]*dv10; 
  out[21] += 1.732050807568877*pt_conf_0[2]*hamil[7]*dv10; 
  out[22] += 3.872983346207417*pt_conf_0[0]*hamil[8]*dv10; 
  out[23] += 3.872983346207417*pt_conf_0[3]*hamil[6]*dv10; 
  out[27] += 1.732050807568877*pt_conf_0[3]*hamil[7]*dv10; 
  out[28] += 3.872983346207417*pt_conf_0[1]*hamil[8]*dv10; 
  out[29] += 3.872983346207417*pt_conf_0[2]*hamil[8]*dv10; 
  out[32] += 3.872983346207417*pt_conf_0[3]*hamil[8]*dv10; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[4]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_1[0]*hamil[2]*dv11; 
  out[1] += 1.732050807568877*pt_conf_1[1]*hamil[2]*dv11; 
  out[2] += 1.732050807568877*hamil[2]*pt_conf_1[2]*dv11; 
  out[3] += 1.732050807568877*pt_conf_1[0]*hamil[3]*dv11; 
  out[4] += 3.872983346207417*pt_conf_1[0]*hamil[5]*dv11; 
  out[5] += 1.732050807568877*hamil[2]*pt_conf_1[3]*dv11; 
  out[6] += 1.732050807568877*pt_conf_1[1]*hamil[3]*dv11; 
  out[7] += 1.732050807568877*pt_conf_1[2]*hamil[3]*dv11; 
  out[8] += 3.872983346207417*pt_conf_1[1]*hamil[5]*dv11; 
  out[9] += 3.872983346207417*pt_conf_1[2]*hamil[5]*dv11; 
  out[10] += 3.872983346207417*pt_conf_1[0]*hamil[7]*dv11; 
  out[11] += 1.732050807568877*pt_conf_1[0]*hamil[6]*dv11; 
  out[13] += 1.732050807568877*hamil[3]*pt_conf_1[3]*dv11; 
  out[14] += 3.872983346207417*pt_conf_1[3]*hamil[5]*dv11; 
  out[15] += 3.872983346207417*pt_conf_1[1]*hamil[7]*dv11; 
  out[16] += 3.872983346207417*pt_conf_1[2]*hamil[7]*dv11; 
  out[17] += 1.732050807568877*pt_conf_1[1]*hamil[6]*dv11; 
  out[18] += 1.732050807568877*pt_conf_1[2]*hamil[6]*dv11; 
  out[19] += 3.872983346207417*pt_conf_1[0]*hamil[8]*dv11; 
  out[23] += 3.872983346207417*pt_conf_1[3]*hamil[7]*dv11; 
  out[24] += 1.732050807568877*pt_conf_1[3]*hamil[6]*dv11; 
  out[25] += 3.872983346207417*pt_conf_1[1]*hamil[8]*dv11; 
  out[26] += 3.872983346207417*pt_conf_1[2]*hamil[8]*dv11; 
  out[31] += 3.872983346207417*pt_conf_1[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_2 = &poisson_tensor_conf[8]; 
  if (dir == 1) { 
  out[0] += 1.732050807568877*pt_conf_2[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_2[1]*dv10; 
  out[2] += 1.732050807568877*hamil[1]*pt_conf_2[2]*dv10; 
  out[3] += 3.872983346207417*pt_conf_2[0]*hamil[4]*dv10; 
  out[4] += 1.732050807568877*pt_conf_2[0]*hamil[3]*dv10; 
  out[5] += 1.732050807568877*hamil[1]*pt_conf_2[3]*dv10; 
  out[6] += 3.872983346207417*pt_conf_2[1]*hamil[4]*dv10; 
  out[7] += 3.872983346207417*pt_conf_2[2]*hamil[4]*dv10; 
  out[8] += 1.732050807568877*pt_conf_2[1]*hamil[3]*dv10; 
  out[9] += 1.732050807568877*pt_conf_2[2]*hamil[3]*dv10; 
  out[10] += 3.872983346207417*pt_conf_2[0]*hamil[6]*dv10; 
  out[12] += 1.732050807568877*pt_conf_2[0]*hamil[7]*dv10; 
  out[13] += 3.872983346207417*pt_conf_2[3]*hamil[4]*dv10; 
  out[14] += 1.732050807568877*hamil[3]*pt_conf_2[3]*dv10; 
  out[15] += 3.872983346207417*pt_conf_2[1]*hamil[6]*dv10; 
  out[16] += 3.872983346207417*pt_conf_2[2]*hamil[6]*dv10; 
  out[20] += 1.732050807568877*pt_conf_2[1]*hamil[7]*dv10; 
  out[21] += 1.732050807568877*pt_conf_2[2]*hamil[7]*dv10; 
  out[22] += 3.872983346207417*pt_conf_2[0]*hamil[8]*dv10; 
  out[23] += 3.872983346207417*pt_conf_2[3]*hamil[6]*dv10; 
  out[27] += 1.732050807568877*pt_conf_2[3]*hamil[7]*dv10; 
  out[28] += 3.872983346207417*pt_conf_2[1]*hamil[8]*dv10; 
  out[29] += 3.872983346207417*pt_conf_2[2]*hamil[8]*dv10; 
  out[32] += 3.872983346207417*pt_conf_2[3]*hamil[8]*dv10; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_3 = &poisson_tensor_conf[12]; 
  if (dir == 1) { 
  out[0] += 1.732050807568877*pt_conf_3[0]*hamil[2]*dv11; 
  out[1] += 1.732050807568877*pt_conf_3[1]*hamil[2]*dv11; 
  out[2] += 1.732050807568877*hamil[2]*pt_conf_3[2]*dv11; 
  out[3] += 1.732050807568877*pt_conf_3[0]*hamil[3]*dv11; 
  out[4] += 3.872983346207417*pt_conf_3[0]*hamil[5]*dv11; 
  out[5] += 1.732050807568877*hamil[2]*pt_conf_3[3]*dv11; 
  out[6] += 1.732050807568877*pt_conf_3[1]*hamil[3]*dv11; 
  out[7] += 1.732050807568877*pt_conf_3[2]*hamil[3]*dv11; 
  out[8] += 3.872983346207417*pt_conf_3[1]*hamil[5]*dv11; 
  out[9] += 3.872983346207417*pt_conf_3[2]*hamil[5]*dv11; 
  out[10] += 3.872983346207417*pt_conf_3[0]*hamil[7]*dv11; 
  out[11] += 1.732050807568877*pt_conf_3[0]*hamil[6]*dv11; 
  out[13] += 1.732050807568877*hamil[3]*pt_conf_3[3]*dv11; 
  out[14] += 3.872983346207417*pt_conf_3[3]*hamil[5]*dv11; 
  out[15] += 3.872983346207417*pt_conf_3[1]*hamil[7]*dv11; 
  out[16] += 3.872983346207417*pt_conf_3[2]*hamil[7]*dv11; 
  out[17] += 1.732050807568877*pt_conf_3[1]*hamil[6]*dv11; 
  out[18] += 1.732050807568877*pt_conf_3[2]*hamil[6]*dv11; 
  out[19] += 3.872983346207417*pt_conf_3[0]*hamil[8]*dv11; 
  out[23] += 3.872983346207417*pt_conf_3[3]*hamil[7]*dv11; 
  out[24] += 1.732050807568877*pt_conf_3[3]*hamil[6]*dv11; 
  out[25] += 3.872983346207417*pt_conf_3[1]*hamil[8]*dv11; 
  out[26] += 3.872983346207417*pt_conf_3[2]*hamil[8]*dv11; 
  out[31] += 3.872983346207417*pt_conf_3[3]*hamil[8]*dv11; 
  } 
  if (dir == 3) { 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[16]; 
  if (dir == 2) { 
  out[0] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[2]*wx1+0.5*pt_conf_4[0]*hamil[3]*dv0); 
  out[1] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[2]*wx1+0.5*pt_conf_4[1]*hamil[3]*dv0); 
  out[2] += dv11*(1.732050807568877*hamil[2]*pt_conf_4[2]*wx1+0.5*pt_conf_4[2]*hamil[3]*dv0); 
  out[3] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[3]*wx1+0.447213595499958*pt_conf_4[0]*hamil[6]*dv0+0.5*pt_conf_4[0]*hamil[2]*dv0); 
  out[4] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[5]*wx1+1.118033988749895*pt_conf_4[0]*hamil[7]*dv0); 
  out[5] += dv11*(1.732050807568877*hamil[2]*pt_conf_4[3]*wx1+0.5*hamil[3]*pt_conf_4[3]*dv0); 
  out[6] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[3]*wx1+0.447213595499958*pt_conf_4[1]*hamil[6]*dv0+0.5*pt_conf_4[1]*hamil[2]*dv0); 
  out[7] += dv11*(1.732050807568877*pt_conf_4[2]*hamil[3]*wx1+0.447213595499958*pt_conf_4[2]*hamil[6]*dv0+0.5*hamil[2]*pt_conf_4[2]*dv0); 
  out[8] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[5]*wx1+1.118033988749895*pt_conf_4[1]*hamil[7]*dv0); 
  out[9] += dv11*(3.872983346207417*pt_conf_4[2]*hamil[5]*wx1+1.118033988749895*pt_conf_4[2]*hamil[7]*dv0); 
  out[10] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[7]*wx1+pt_conf_4[0]*hamil[8]*dv0+1.118033988749895*pt_conf_4[0]*hamil[5]*dv0); 
  out[11] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[6]*wx1+0.4472135954999579*pt_conf_4[0]*hamil[3]*dv0); 
  out[13] += dv11*(1.732050807568877*hamil[3]*pt_conf_4[3]*wx1+0.447213595499958*pt_conf_4[3]*hamil[6]*dv0+0.5*hamil[2]*pt_conf_4[3]*dv0); 
  out[14] += dv11*(3.872983346207417*pt_conf_4[3]*hamil[5]*wx1+1.118033988749895*pt_conf_4[3]*hamil[7]*dv0); 
  out[15] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[7]*wx1+pt_conf_4[1]*hamil[8]*dv0+1.118033988749895*pt_conf_4[1]*hamil[5]*dv0); 
  out[16] += dv11*(3.872983346207417*pt_conf_4[2]*hamil[7]*wx1+pt_conf_4[2]*hamil[8]*dv0+1.118033988749895*pt_conf_4[2]*hamil[5]*dv0); 
  out[17] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[6]*wx1+0.447213595499958*pt_conf_4[1]*hamil[3]*dv0); 
  out[18] += dv11*(1.732050807568877*pt_conf_4[2]*hamil[6]*wx1+0.447213595499958*pt_conf_4[2]*hamil[3]*dv0); 
  out[19] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[8]*wx1+pt_conf_4[0]*hamil[7]*dv0); 
  out[23] += dv11*(3.872983346207417*pt_conf_4[3]*hamil[7]*wx1+pt_conf_4[3]*hamil[8]*dv0+1.118033988749895*pt_conf_4[3]*hamil[5]*dv0); 
  out[24] += dv11*(1.732050807568877*pt_conf_4[3]*hamil[6]*wx1+0.4472135954999579*hamil[3]*pt_conf_4[3]*dv0); 
  out[25] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[8]*wx1+1.0*pt_conf_4[1]*hamil[7]*dv0); 
  out[26] += dv11*(3.872983346207417*pt_conf_4[2]*hamil[8]*wx1+1.0*pt_conf_4[2]*hamil[7]*dv0); 
  out[31] += dv11*(3.872983346207417*pt_conf_4[3]*hamil[8]*wx1+pt_conf_4[3]*hamil[7]*dv0); 
  } 
  if (dir == 3) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[1]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[4]*dv0); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_4[1]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[4]*dv0); 
  out[2] += dv11*((-1.732050807568877*hamil[1]*pt_conf_4[2]*wx1)-1.118033988749895*pt_conf_4[2]*hamil[4]*dv0); 
  out[3] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[4]*wx1)-0.5*pt_conf_4[0]*hamil[1]*dv0); 
  out[4] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[3]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[6]*dv0); 
  out[5] += dv11*((-1.732050807568877*hamil[1]*pt_conf_4[3]*wx1)-1.118033988749895*pt_conf_4[3]*hamil[4]*dv0); 
  out[6] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[4]*wx1)-0.5*hamil[1]*pt_conf_4[1]*dv0); 
  out[7] += dv11*((-3.872983346207417*pt_conf_4[2]*hamil[4]*wx1)-0.5*hamil[1]*pt_conf_4[2]*dv0); 
  out[8] += dv11*((-1.732050807568877*pt_conf_4[1]*hamil[3]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[6]*dv0); 
  out[9] += dv11*((-1.732050807568877*pt_conf_4[2]*hamil[3]*wx1)-1.118033988749895*pt_conf_4[2]*hamil[6]*dv0); 
  out[10] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[6]*wx1)-0.5*pt_conf_4[0]*hamil[3]*dv0); 
  out[11] += -1.0*pt_conf_4[0]*hamil[4]*dv0*dv11; 
  out[12] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[8]*dv0); 
  out[13] += dv11*((-3.872983346207417*pt_conf_4[3]*hamil[4]*wx1)-0.5*hamil[1]*pt_conf_4[3]*dv0); 
  out[14] += dv11*((-1.732050807568877*hamil[3]*pt_conf_4[3]*wx1)-1.118033988749895*pt_conf_4[3]*hamil[6]*dv0); 
  out[15] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[6]*wx1)-0.5*pt_conf_4[1]*hamil[3]*dv0); 
  out[16] += dv11*((-3.872983346207417*pt_conf_4[2]*hamil[6]*wx1)-0.5*pt_conf_4[2]*hamil[3]*dv0); 
  out[17] += -1.0*pt_conf_4[1]*hamil[4]*dv0*dv11; 
  out[18] += -1.0*pt_conf_4[2]*hamil[4]*dv0*dv11; 
  out[19] += -1.0*pt_conf_4[0]*hamil[6]*dv0*dv11; 
  out[20] += dv11*((-1.732050807568877*pt_conf_4[1]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[8]*dv0); 
  out[21] += dv11*((-1.732050807568877*pt_conf_4[2]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[2]*hamil[8]*dv0); 
  out[22] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[8]*wx1)-0.5*pt_conf_4[0]*hamil[7]*dv0); 
  out[23] += dv11*((-3.872983346207417*pt_conf_4[3]*hamil[6]*wx1)-0.5*hamil[3]*pt_conf_4[3]*dv0); 
  out[24] += -1.0*pt_conf_4[3]*hamil[4]*dv0*dv11; 
  out[25] += -1.0*pt_conf_4[1]*hamil[6]*dv0*dv11; 
  out[26] += -1.0*pt_conf_4[2]*hamil[6]*dv0*dv11; 
  out[27] += dv11*((-1.732050807568877*pt_conf_4[3]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[3]*hamil[8]*dv0); 
  out[28] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[8]*wx1)-0.5000000000000001*pt_conf_4[1]*hamil[7]*dv0); 
  out[29] += dv11*((-3.872983346207417*pt_conf_4[2]*hamil[8]*wx1)-0.5000000000000001*pt_conf_4[2]*hamil[7]*dv0); 
  out[30] += -1.0*pt_conf_4[0]*hamil[8]*dv0*dv11; 
  out[31] += -1.0*pt_conf_4[3]*hamil[6]*dv0*dv11; 
  out[32] += dv11*((-3.872983346207417*pt_conf_4[3]*hamil[8]*wx1)-0.5*pt_conf_4[3]*hamil[7]*dv0); 
  out[33] += -1.0*pt_conf_4[1]*hamil[8]*dv0*dv11; 
  out[34] += -1.0*pt_conf_4[2]*hamil[8]*dv0*dv11; 
  out[35] += -1.0*pt_conf_4[3]*hamil[8]*dv0*dv11; 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[20]; 
  if (dir == 2) { 
  out[0] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[2]*wx2+1.118033988749895*pt_conf_5[0]*hamil[5]*dv1); 
  out[1] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[2]*wx2+1.118033988749895*pt_conf_5[1]*hamil[5]*dv1); 
  out[2] += dv11*(1.732050807568877*hamil[2]*pt_conf_5[2]*wx2+1.118033988749895*pt_conf_5[2]*hamil[5]*dv1); 
  out[3] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[3]*wx2+1.118033988749895*pt_conf_5[0]*hamil[7]*dv1); 
  out[4] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[5]*wx2+0.5*pt_conf_5[0]*hamil[2]*dv1); 
  out[5] += dv11*(1.732050807568877*hamil[2]*pt_conf_5[3]*wx2+1.118033988749895*pt_conf_5[3]*hamil[5]*dv1); 
  out[6] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[3]*wx2+1.118033988749895*pt_conf_5[1]*hamil[7]*dv1); 
  out[7] += dv11*(1.732050807568877*pt_conf_5[2]*hamil[3]*wx2+1.118033988749895*pt_conf_5[2]*hamil[7]*dv1); 
  out[8] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[5]*wx2+0.5*pt_conf_5[1]*hamil[2]*dv1); 
  out[9] += dv11*(3.872983346207417*pt_conf_5[2]*hamil[5]*wx2+0.5*hamil[2]*pt_conf_5[2]*dv1); 
  out[10] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[7]*wx2+0.5*pt_conf_5[0]*hamil[3]*dv1); 
  out[11] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[6]*wx2+1.118033988749895*pt_conf_5[0]*hamil[8]*dv1); 
  out[12] += pt_conf_5[0]*hamil[5]*dv1*dv11; 
  out[13] += dv11*(1.732050807568877*hamil[3]*pt_conf_5[3]*wx2+1.118033988749895*pt_conf_5[3]*hamil[7]*dv1); 
  out[14] += dv11*(3.872983346207417*pt_conf_5[3]*hamil[5]*wx2+0.5*hamil[2]*pt_conf_5[3]*dv1); 
  out[15] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[7]*wx2+0.5*pt_conf_5[1]*hamil[3]*dv1); 
  out[16] += dv11*(3.872983346207417*pt_conf_5[2]*hamil[7]*wx2+0.5*pt_conf_5[2]*hamil[3]*dv1); 
  out[17] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[6]*wx2+1.118033988749895*pt_conf_5[1]*hamil[8]*dv1); 
  out[18] += dv11*(1.732050807568877*pt_conf_5[2]*hamil[6]*wx2+1.118033988749895*pt_conf_5[2]*hamil[8]*dv1); 
  out[19] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[8]*wx2+0.5*pt_conf_5[0]*hamil[6]*dv1); 
  out[20] += 1.0*pt_conf_5[1]*hamil[5]*dv1*dv11; 
  out[21] += 1.0*pt_conf_5[2]*hamil[5]*dv1*dv11; 
  out[22] += pt_conf_5[0]*hamil[7]*dv1*dv11; 
  out[23] += dv11*(3.872983346207417*pt_conf_5[3]*hamil[7]*wx2+0.5*hamil[3]*pt_conf_5[3]*dv1); 
  out[24] += dv11*(1.732050807568877*pt_conf_5[3]*hamil[6]*wx2+1.118033988749895*pt_conf_5[3]*hamil[8]*dv1); 
  out[25] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[8]*wx2+0.5000000000000001*pt_conf_5[1]*hamil[6]*dv1); 
  out[26] += dv11*(3.872983346207417*pt_conf_5[2]*hamil[8]*wx2+0.5000000000000001*pt_conf_5[2]*hamil[6]*dv1); 
  out[27] += pt_conf_5[3]*hamil[5]*dv1*dv11; 
  out[28] += 1.0*pt_conf_5[1]*hamil[7]*dv1*dv11; 
  out[29] += 1.0*pt_conf_5[2]*hamil[7]*dv1*dv11; 
  out[30] += pt_conf_5[0]*hamil[8]*dv1*dv11; 
  out[31] += dv11*(3.872983346207417*pt_conf_5[3]*hamil[8]*wx2+0.5*pt_conf_5[3]*hamil[6]*dv1); 
  out[32] += pt_conf_5[3]*hamil[7]*dv1*dv11; 
  out[33] += pt_conf_5[1]*hamil[8]*dv1*dv11; 
  out[34] += pt_conf_5[2]*hamil[8]*dv1*dv11; 
  out[35] += pt_conf_5[3]*hamil[8]*dv1*dv11; 
  } 
  if (dir == 3) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[1]*wx2)-0.5*pt_conf_5[0]*hamil[3]*dv1); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_5[1]*wx2)-0.5*pt_conf_5[1]*hamil[3]*dv1); 
  out[2] += dv11*((-1.732050807568877*hamil[1]*pt_conf_5[2]*wx2)-0.5*pt_conf_5[2]*hamil[3]*dv1); 
  out[3] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[0]*hamil[6]*dv1); 
  out[4] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[3]*wx2)-0.447213595499958*pt_conf_5[0]*hamil[7]*dv1-0.5*pt_conf_5[0]*hamil[1]*dv1); 
  out[5] += dv11*((-1.732050807568877*hamil[1]*pt_conf_5[3]*wx2)-0.5*hamil[3]*pt_conf_5[3]*dv1); 
  out[6] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[1]*hamil[6]*dv1); 
  out[7] += dv11*((-3.872983346207417*pt_conf_5[2]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[2]*hamil[6]*dv1); 
  out[8] += dv11*((-1.732050807568877*pt_conf_5[1]*hamil[3]*wx2)-0.447213595499958*pt_conf_5[1]*hamil[7]*dv1-0.5*hamil[1]*pt_conf_5[1]*dv1); 
  out[9] += dv11*((-1.732050807568877*pt_conf_5[2]*hamil[3]*wx2)-0.447213595499958*pt_conf_5[2]*hamil[7]*dv1-0.5*hamil[1]*pt_conf_5[2]*dv1); 
  out[10] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[6]*wx2)-1.0*pt_conf_5[0]*hamil[8]*dv1-1.118033988749895*pt_conf_5[0]*hamil[4]*dv1); 
  out[12] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[7]*wx2)-0.4472135954999579*pt_conf_5[0]*hamil[3]*dv1); 
  out[13] += dv11*((-3.872983346207417*pt_conf_5[3]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[3]*hamil[6]*dv1); 
  out[14] += dv11*((-1.732050807568877*hamil[3]*pt_conf_5[3]*wx2)-0.447213595499958*pt_conf_5[3]*hamil[7]*dv1-0.5*hamil[1]*pt_conf_5[3]*dv1); 
  out[15] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[6]*wx2)-1.0*pt_conf_5[1]*hamil[8]*dv1-1.118033988749895*pt_conf_5[1]*hamil[4]*dv1); 
  out[16] += dv11*((-3.872983346207417*pt_conf_5[2]*hamil[6]*wx2)-1.0*pt_conf_5[2]*hamil[8]*dv1-1.118033988749895*pt_conf_5[2]*hamil[4]*dv1); 
  out[20] += dv11*((-1.732050807568877*pt_conf_5[1]*hamil[7]*wx2)-0.447213595499958*pt_conf_5[1]*hamil[3]*dv1); 
  out[21] += dv11*((-1.732050807568877*pt_conf_5[2]*hamil[7]*wx2)-0.447213595499958*pt_conf_5[2]*hamil[3]*dv1); 
  out[22] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[8]*wx2)-1.0*pt_conf_5[0]*hamil[6]*dv1); 
  out[23] += dv11*((-3.872983346207417*pt_conf_5[3]*hamil[6]*wx2)-1.0*pt_conf_5[3]*hamil[8]*dv1-1.118033988749895*pt_conf_5[3]*hamil[4]*dv1); 
  out[27] += dv11*((-1.732050807568877*pt_conf_5[3]*hamil[7]*wx2)-0.4472135954999579*hamil[3]*pt_conf_5[3]*dv1); 
  out[28] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[8]*wx2)-1.0*pt_conf_5[1]*hamil[6]*dv1); 
  out[29] += dv11*((-3.872983346207417*pt_conf_5[2]*hamil[8]*wx2)-1.0*pt_conf_5[2]*hamil[6]*dv1); 
  out[32] += dv11*((-3.872983346207417*pt_conf_5[3]*hamil[8]*wx2)-1.0*pt_conf_5[3]*hamil[6]*dv1); 
  } 

} 
