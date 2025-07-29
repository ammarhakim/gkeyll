#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[2]; 
  const double wv = w[2]; 

  double Ghat[16]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[23]+0.4330127018922193*hamil[7]*fskin[15]+0.75*hamil[5]*fskin[12]+0.75*hamil[4]*fskin[9]+0.4330127018922193*fskin[5]*hamil[5]+0.4330127018922193*fskin[4]*hamil[4]+0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fskin[28]+0.4330127018922193*hamil[7]*fskin[24]+0.75*hamil[5]*fskin[20]+0.75*hamil[4]*fskin[17]+0.4330127018922193*hamil[5]*fskin[13]+0.4330127018922193*hamil[4]*fskin[10]+0.75*hamil[1]*fskin[6]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = 0.75*hamil[7]*fskin[29]+0.4330127018922193*hamil[7]*fskin[25]+0.75*hamil[5]*fskin[21]+0.75*hamil[4]*fskin[18]+0.4330127018922193*hamil[5]*fskin[14]+0.4330127018922193*hamil[4]*fskin[11]+0.75*hamil[1]*fskin[7]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = 0.75*hamil[5]*fskin[23]+0.4330127018922193*hamil[5]*fskin[15]+0.75*hamil[7]*fskin[12]+0.75*hamil[1]*fskin[9]+0.4330127018922193*fskin[5]*hamil[7]+0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = 0.75*hamil[4]*fskin[23]+0.4330127018922193*hamil[4]*fskin[15]+0.75*hamil[1]*fskin[12]+0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[4]*hamil[7]+0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[1]*fskin[5]; 
  Ghat[5] = 0.75*hamil[7]*fskin[31]+0.4330127018922193*hamil[7]*fskin[30]+0.75*hamil[5]*fskin[27]+0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[5]*fskin[22]+0.4330127018922193*hamil[4]*fskin[19]+0.75*hamil[1]*fskin[16]+0.4330127018922193*hamil[1]*fskin[8]; 
  Ghat[6] = 0.75*hamil[5]*fskin[28]+0.4330127018922193*hamil[5]*fskin[24]+0.75*hamil[7]*fskin[20]+0.75*hamil[1]*fskin[17]+0.4330127018922193*hamil[7]*fskin[13]+0.4330127018922193*hamil[1]*fskin[10]+0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[7] = 0.75*hamil[5]*fskin[29]+0.4330127018922193*hamil[5]*fskin[25]+0.75*hamil[7]*fskin[21]+0.75*hamil[1]*fskin[18]+0.4330127018922193*hamil[7]*fskin[14]+0.4330127018922193*hamil[1]*fskin[11]+0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[8] = 0.75*hamil[4]*fskin[28]+0.4330127018922193*hamil[4]*fskin[24]+0.75*hamil[1]*fskin[20]+0.75*hamil[7]*fskin[17]+0.4330127018922193*hamil[1]*fskin[13]+0.4330127018922193*hamil[7]*fskin[10]+0.75*hamil[5]*fskin[6]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[9] = 0.75*hamil[4]*fskin[29]+0.4330127018922193*hamil[4]*fskin[25]+0.75*hamil[1]*fskin[21]+0.75*hamil[7]*fskin[18]+0.4330127018922193*hamil[1]*fskin[14]+0.4330127018922193*hamil[7]*fskin[11]+0.75*hamil[5]*fskin[7]+0.4330127018922193*fskin[3]*hamil[5]; 
  Ghat[10] = 0.75*hamil[1]*fskin[23]+0.4330127018922193*hamil[1]*fskin[15]+0.75*hamil[4]*fskin[12]+0.75*hamil[5]*fskin[9]+0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[5]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[11] = 0.75*hamil[5]*fskin[31]+0.4330127018922193*hamil[5]*fskin[30]+0.75*hamil[7]*fskin[27]+0.75*hamil[1]*fskin[26]+0.4330127018922193*hamil[7]*fskin[22]+0.4330127018922193*hamil[1]*fskin[19]+0.75*hamil[4]*fskin[16]+0.4330127018922193*hamil[4]*fskin[8]; 
  Ghat[12] = 0.75*hamil[4]*fskin[31]+0.4330127018922193*hamil[4]*fskin[30]+0.75*hamil[1]*fskin[27]+0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[1]*fskin[22]+0.4330127018922193*hamil[7]*fskin[19]+0.75*hamil[5]*fskin[16]+0.4330127018922193*hamil[5]*fskin[8]; 
  Ghat[13] = 0.75*hamil[1]*fskin[28]+0.4330127018922193*hamil[1]*fskin[24]+0.75*hamil[4]*fskin[20]+0.75*hamil[5]*fskin[17]+0.4330127018922193*hamil[4]*fskin[13]+0.4330127018922193*hamil[5]*fskin[10]+0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[14] = 0.75*hamil[1]*fskin[29]+0.4330127018922193*hamil[1]*fskin[25]+0.75*hamil[4]*fskin[21]+0.75*hamil[5]*fskin[18]+0.4330127018922193*hamil[4]*fskin[14]+0.4330127018922193*hamil[5]*fskin[11]+0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[15] = 0.75*hamil[1]*fskin[31]+0.4330127018922193*hamil[1]*fskin[30]+0.75*hamil[4]*fskin[27]+0.75*hamil[5]*fskin[26]+0.4330127018922193*hamil[4]*fskin[22]+0.4330127018922193*hamil[5]*fskin[19]+0.75*hamil[7]*fskin[16]+0.4330127018922193*hamil[7]*fskin[8]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[23])+0.4330127018922193*hamil[7]*fedge[15]-0.75*hamil[5]*fedge[12]-0.75*hamil[4]*fedge[9]+0.4330127018922193*fedge[5]*hamil[5]+0.4330127018922193*fedge[4]*hamil[4]-0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[28])+0.4330127018922193*hamil[7]*fedge[24]-0.75*hamil[5]*fedge[20]-0.75*hamil[4]*fedge[17]+0.4330127018922193*hamil[5]*fedge[13]+0.4330127018922193*hamil[4]*fedge[10]-0.75*hamil[1]*fedge[6]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = -(0.75*hamil[7]*fedge[29])+0.4330127018922193*hamil[7]*fedge[25]-0.75*hamil[5]*fedge[21]-0.75*hamil[4]*fedge[18]+0.4330127018922193*hamil[5]*fedge[14]+0.4330127018922193*hamil[4]*fedge[11]-0.75*hamil[1]*fedge[7]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = -(0.75*hamil[5]*fedge[23])+0.4330127018922193*hamil[5]*fedge[15]-0.75*hamil[7]*fedge[12]-0.75*hamil[1]*fedge[9]+0.4330127018922193*fedge[5]*hamil[7]-0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[4]*fedge[23])+0.4330127018922193*hamil[4]*fedge[15]-0.75*hamil[1]*fedge[12]-0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[4]*hamil[7]-0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[1]*fedge[5]; 
  Ghat[5] = -(0.75*hamil[7]*fedge[31])+0.4330127018922193*hamil[7]*fedge[30]-0.75*hamil[5]*fedge[27]-0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[5]*fedge[22]+0.4330127018922193*hamil[4]*fedge[19]-0.75*hamil[1]*fedge[16]+0.4330127018922193*hamil[1]*fedge[8]; 
  Ghat[6] = -(0.75*hamil[5]*fedge[28])+0.4330127018922193*hamil[5]*fedge[24]-0.75*hamil[7]*fedge[20]-0.75*hamil[1]*fedge[17]+0.4330127018922193*hamil[7]*fedge[13]+0.4330127018922193*hamil[1]*fedge[10]-0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[7] = -(0.75*hamil[5]*fedge[29])+0.4330127018922193*hamil[5]*fedge[25]-0.75*hamil[7]*fedge[21]-0.75*hamil[1]*fedge[18]+0.4330127018922193*hamil[7]*fedge[14]+0.4330127018922193*hamil[1]*fedge[11]-0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[4]*fedge[28])+0.4330127018922193*hamil[4]*fedge[24]-0.75*hamil[1]*fedge[20]-0.75*hamil[7]*fedge[17]+0.4330127018922193*hamil[1]*fedge[13]+0.4330127018922193*hamil[7]*fedge[10]-0.75*hamil[5]*fedge[6]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[9] = -(0.75*hamil[4]*fedge[29])+0.4330127018922193*hamil[4]*fedge[25]-0.75*hamil[1]*fedge[21]-0.75*hamil[7]*fedge[18]+0.4330127018922193*hamil[1]*fedge[14]+0.4330127018922193*hamil[7]*fedge[11]-0.75*hamil[5]*fedge[7]+0.4330127018922193*fedge[3]*hamil[5]; 
  Ghat[10] = -(0.75*hamil[1]*fedge[23])+0.4330127018922193*hamil[1]*fedge[15]-0.75*hamil[4]*fedge[12]-0.75*hamil[5]*fedge[9]-0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[5]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[11] = -(0.75*hamil[5]*fedge[31])+0.4330127018922193*hamil[5]*fedge[30]-0.75*hamil[7]*fedge[27]-0.75*hamil[1]*fedge[26]+0.4330127018922193*hamil[7]*fedge[22]+0.4330127018922193*hamil[1]*fedge[19]-0.75*hamil[4]*fedge[16]+0.4330127018922193*hamil[4]*fedge[8]; 
  Ghat[12] = -(0.75*hamil[4]*fedge[31])+0.4330127018922193*hamil[4]*fedge[30]-0.75*hamil[1]*fedge[27]-0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[1]*fedge[22]+0.4330127018922193*hamil[7]*fedge[19]-0.75*hamil[5]*fedge[16]+0.4330127018922193*hamil[5]*fedge[8]; 
  Ghat[13] = -(0.75*hamil[1]*fedge[28])+0.4330127018922193*hamil[1]*fedge[24]-0.75*hamil[4]*fedge[20]-0.75*hamil[5]*fedge[17]+0.4330127018922193*hamil[4]*fedge[13]+0.4330127018922193*hamil[5]*fedge[10]-0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[14] = -(0.75*hamil[1]*fedge[29])+0.4330127018922193*hamil[1]*fedge[25]-0.75*hamil[4]*fedge[21]-0.75*hamil[5]*fedge[18]+0.4330127018922193*hamil[4]*fedge[14]+0.4330127018922193*hamil[5]*fedge[11]-0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[15] = -(0.75*hamil[1]*fedge[31])+0.4330127018922193*hamil[1]*fedge[30]-0.75*hamil[4]*fedge[27]-0.75*hamil[5]*fedge[26]+0.4330127018922193*hamil[4]*fedge[22]+0.4330127018922193*hamil[5]*fedge[19]-0.75*hamil[7]*fedge[16]+0.4330127018922193*hamil[7]*fedge[8]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[8] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[9] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[13] += -(0.7071067811865475*Ghat[8]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv10*dx10); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv10*dx10); 
  out[16] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[18] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[19] += -(0.7071067811865475*Ghat[11]*dv10*dx10); 
  out[20] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[21] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[22] += -(0.7071067811865475*Ghat[12]*dv10*dx10); 
  out[23] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[24] += -(0.7071067811865475*Ghat[13]*dv10*dx10); 
  out[25] += -(0.7071067811865475*Ghat[14]*dv10*dx10); 
  out[26] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[27] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[28] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[29] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[30] += -(0.7071067811865475*Ghat[15]*dv10*dx10); 
  out[31] += -(1.224744871391589*Ghat[15]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[23]+0.4330127018922193*hamil[7]*fedge[15]+0.75*hamil[5]*fedge[12]+0.75*hamil[4]*fedge[9]+0.4330127018922193*fedge[5]*hamil[5]+0.4330127018922193*fedge[4]*hamil[4]+0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fedge[28]+0.4330127018922193*hamil[7]*fedge[24]+0.75*hamil[5]*fedge[20]+0.75*hamil[4]*fedge[17]+0.4330127018922193*hamil[5]*fedge[13]+0.4330127018922193*hamil[4]*fedge[10]+0.75*hamil[1]*fedge[6]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = 0.75*hamil[7]*fedge[29]+0.4330127018922193*hamil[7]*fedge[25]+0.75*hamil[5]*fedge[21]+0.75*hamil[4]*fedge[18]+0.4330127018922193*hamil[5]*fedge[14]+0.4330127018922193*hamil[4]*fedge[11]+0.75*hamil[1]*fedge[7]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = 0.75*hamil[5]*fedge[23]+0.4330127018922193*hamil[5]*fedge[15]+0.75*hamil[7]*fedge[12]+0.75*hamil[1]*fedge[9]+0.4330127018922193*fedge[5]*hamil[7]+0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = 0.75*hamil[4]*fedge[23]+0.4330127018922193*hamil[4]*fedge[15]+0.75*hamil[1]*fedge[12]+0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[4]*hamil[7]+0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[1]*fedge[5]; 
  Ghat[5] = 0.75*hamil[7]*fedge[31]+0.4330127018922193*hamil[7]*fedge[30]+0.75*hamil[5]*fedge[27]+0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[5]*fedge[22]+0.4330127018922193*hamil[4]*fedge[19]+0.75*hamil[1]*fedge[16]+0.4330127018922193*hamil[1]*fedge[8]; 
  Ghat[6] = 0.75*hamil[5]*fedge[28]+0.4330127018922193*hamil[5]*fedge[24]+0.75*hamil[7]*fedge[20]+0.75*hamil[1]*fedge[17]+0.4330127018922193*hamil[7]*fedge[13]+0.4330127018922193*hamil[1]*fedge[10]+0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[7] = 0.75*hamil[5]*fedge[29]+0.4330127018922193*hamil[5]*fedge[25]+0.75*hamil[7]*fedge[21]+0.75*hamil[1]*fedge[18]+0.4330127018922193*hamil[7]*fedge[14]+0.4330127018922193*hamil[1]*fedge[11]+0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[8] = 0.75*hamil[4]*fedge[28]+0.4330127018922193*hamil[4]*fedge[24]+0.75*hamil[1]*fedge[20]+0.75*hamil[7]*fedge[17]+0.4330127018922193*hamil[1]*fedge[13]+0.4330127018922193*hamil[7]*fedge[10]+0.75*hamil[5]*fedge[6]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[9] = 0.75*hamil[4]*fedge[29]+0.4330127018922193*hamil[4]*fedge[25]+0.75*hamil[1]*fedge[21]+0.75*hamil[7]*fedge[18]+0.4330127018922193*hamil[1]*fedge[14]+0.4330127018922193*hamil[7]*fedge[11]+0.75*hamil[5]*fedge[7]+0.4330127018922193*fedge[3]*hamil[5]; 
  Ghat[10] = 0.75*hamil[1]*fedge[23]+0.4330127018922193*hamil[1]*fedge[15]+0.75*hamil[4]*fedge[12]+0.75*hamil[5]*fedge[9]+0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[5]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[11] = 0.75*hamil[5]*fedge[31]+0.4330127018922193*hamil[5]*fedge[30]+0.75*hamil[7]*fedge[27]+0.75*hamil[1]*fedge[26]+0.4330127018922193*hamil[7]*fedge[22]+0.4330127018922193*hamil[1]*fedge[19]+0.75*hamil[4]*fedge[16]+0.4330127018922193*hamil[4]*fedge[8]; 
  Ghat[12] = 0.75*hamil[4]*fedge[31]+0.4330127018922193*hamil[4]*fedge[30]+0.75*hamil[1]*fedge[27]+0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[1]*fedge[22]+0.4330127018922193*hamil[7]*fedge[19]+0.75*hamil[5]*fedge[16]+0.4330127018922193*hamil[5]*fedge[8]; 
  Ghat[13] = 0.75*hamil[1]*fedge[28]+0.4330127018922193*hamil[1]*fedge[24]+0.75*hamil[4]*fedge[20]+0.75*hamil[5]*fedge[17]+0.4330127018922193*hamil[4]*fedge[13]+0.4330127018922193*hamil[5]*fedge[10]+0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[14] = 0.75*hamil[1]*fedge[29]+0.4330127018922193*hamil[1]*fedge[25]+0.75*hamil[4]*fedge[21]+0.75*hamil[5]*fedge[18]+0.4330127018922193*hamil[4]*fedge[14]+0.4330127018922193*hamil[5]*fedge[11]+0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[15] = 0.75*hamil[1]*fedge[31]+0.4330127018922193*hamil[1]*fedge[30]+0.75*hamil[4]*fedge[27]+0.75*hamil[5]*fedge[26]+0.4330127018922193*hamil[4]*fedge[22]+0.4330127018922193*hamil[5]*fedge[19]+0.75*hamil[7]*fedge[16]+0.4330127018922193*hamil[7]*fedge[8]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[23])+0.4330127018922193*hamil[7]*fskin[15]-0.75*hamil[5]*fskin[12]-0.75*hamil[4]*fskin[9]+0.4330127018922193*fskin[5]*hamil[5]+0.4330127018922193*fskin[4]*hamil[4]-0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[28])+0.4330127018922193*hamil[7]*fskin[24]-0.75*hamil[5]*fskin[20]-0.75*hamil[4]*fskin[17]+0.4330127018922193*hamil[5]*fskin[13]+0.4330127018922193*hamil[4]*fskin[10]-0.75*hamil[1]*fskin[6]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = -(0.75*hamil[7]*fskin[29])+0.4330127018922193*hamil[7]*fskin[25]-0.75*hamil[5]*fskin[21]-0.75*hamil[4]*fskin[18]+0.4330127018922193*hamil[5]*fskin[14]+0.4330127018922193*hamil[4]*fskin[11]-0.75*hamil[1]*fskin[7]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = -(0.75*hamil[5]*fskin[23])+0.4330127018922193*hamil[5]*fskin[15]-0.75*hamil[7]*fskin[12]-0.75*hamil[1]*fskin[9]+0.4330127018922193*fskin[5]*hamil[7]-0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[4]*fskin[23])+0.4330127018922193*hamil[4]*fskin[15]-0.75*hamil[1]*fskin[12]-0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[4]*hamil[7]-0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[1]*fskin[5]; 
  Ghat[5] = -(0.75*hamil[7]*fskin[31])+0.4330127018922193*hamil[7]*fskin[30]-0.75*hamil[5]*fskin[27]-0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[5]*fskin[22]+0.4330127018922193*hamil[4]*fskin[19]-0.75*hamil[1]*fskin[16]+0.4330127018922193*hamil[1]*fskin[8]; 
  Ghat[6] = -(0.75*hamil[5]*fskin[28])+0.4330127018922193*hamil[5]*fskin[24]-0.75*hamil[7]*fskin[20]-0.75*hamil[1]*fskin[17]+0.4330127018922193*hamil[7]*fskin[13]+0.4330127018922193*hamil[1]*fskin[10]-0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[7] = -(0.75*hamil[5]*fskin[29])+0.4330127018922193*hamil[5]*fskin[25]-0.75*hamil[7]*fskin[21]-0.75*hamil[1]*fskin[18]+0.4330127018922193*hamil[7]*fskin[14]+0.4330127018922193*hamil[1]*fskin[11]-0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[4]*fskin[28])+0.4330127018922193*hamil[4]*fskin[24]-0.75*hamil[1]*fskin[20]-0.75*hamil[7]*fskin[17]+0.4330127018922193*hamil[1]*fskin[13]+0.4330127018922193*hamil[7]*fskin[10]-0.75*hamil[5]*fskin[6]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[9] = -(0.75*hamil[4]*fskin[29])+0.4330127018922193*hamil[4]*fskin[25]-0.75*hamil[1]*fskin[21]-0.75*hamil[7]*fskin[18]+0.4330127018922193*hamil[1]*fskin[14]+0.4330127018922193*hamil[7]*fskin[11]-0.75*hamil[5]*fskin[7]+0.4330127018922193*fskin[3]*hamil[5]; 
  Ghat[10] = -(0.75*hamil[1]*fskin[23])+0.4330127018922193*hamil[1]*fskin[15]-0.75*hamil[4]*fskin[12]-0.75*hamil[5]*fskin[9]-0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[5]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[11] = -(0.75*hamil[5]*fskin[31])+0.4330127018922193*hamil[5]*fskin[30]-0.75*hamil[7]*fskin[27]-0.75*hamil[1]*fskin[26]+0.4330127018922193*hamil[7]*fskin[22]+0.4330127018922193*hamil[1]*fskin[19]-0.75*hamil[4]*fskin[16]+0.4330127018922193*hamil[4]*fskin[8]; 
  Ghat[12] = -(0.75*hamil[4]*fskin[31])+0.4330127018922193*hamil[4]*fskin[30]-0.75*hamil[1]*fskin[27]-0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[1]*fskin[22]+0.4330127018922193*hamil[7]*fskin[19]-0.75*hamil[5]*fskin[16]+0.4330127018922193*hamil[5]*fskin[8]; 
  Ghat[13] = -(0.75*hamil[1]*fskin[28])+0.4330127018922193*hamil[1]*fskin[24]-0.75*hamil[4]*fskin[20]-0.75*hamil[5]*fskin[17]+0.4330127018922193*hamil[4]*fskin[13]+0.4330127018922193*hamil[5]*fskin[10]-0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[14] = -(0.75*hamil[1]*fskin[29])+0.4330127018922193*hamil[1]*fskin[25]-0.75*hamil[4]*fskin[21]-0.75*hamil[5]*fskin[18]+0.4330127018922193*hamil[4]*fskin[14]+0.4330127018922193*hamil[5]*fskin[11]-0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[15] = -(0.75*hamil[1]*fskin[31])+0.4330127018922193*hamil[1]*fskin[30]-0.75*hamil[4]*fskin[27]-0.75*hamil[5]*fskin[26]+0.4330127018922193*hamil[4]*fskin[22]+0.4330127018922193*hamil[5]*fskin[19]-0.75*hamil[7]*fskin[16]+0.4330127018922193*hamil[7]*fskin[8]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[6] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[8] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[9] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[10] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[11] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[12] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[13] += 0.7071067811865475*Ghat[8]*dv10*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dv10*dx10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10*dx10; 
  out[16] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[18] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[19] += 0.7071067811865475*Ghat[11]*dv10*dx10; 
  out[20] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[21] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[22] += 0.7071067811865475*Ghat[12]*dv10*dx10; 
  out[23] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[24] += 0.7071067811865475*Ghat[13]*dv10*dx10; 
  out[25] += 0.7071067811865475*Ghat[14]*dv10*dx10; 
  out[26] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[27] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[28] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[29] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[30] += 0.7071067811865475*Ghat[15]*dv10*dx10; 
  out[31] += -(1.224744871391589*Ghat[15]*dv10*dx10); 

  } 
  return 0.0;

} 
