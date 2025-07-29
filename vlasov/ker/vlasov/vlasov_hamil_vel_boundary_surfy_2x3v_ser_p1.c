#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx11 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wv = w[3]; 

  double Ghat[16]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[22]+0.4330127018922193*hamil[7]*fskin[14]+0.75*hamil[6]*fskin[13]+0.75*hamil[4]*fskin[8]+0.4330127018922193*fskin[5]*hamil[6]+0.4330127018922193*fskin[3]*hamil[4]+0.75*fskin[2]*hamil[2]+0.4330127018922193*fskin[0]*hamil[2]; 
  Ghat[1] = 0.75*hamil[7]*fskin[27]+0.4330127018922193*hamil[7]*fskin[21]+0.75*hamil[6]*fskin[20]+0.75*hamil[4]*fskin[16]+0.4330127018922193*hamil[6]*fskin[12]+0.4330127018922193*hamil[4]*fskin[7]+0.75*hamil[2]*fskin[6]+0.4330127018922193*fskin[1]*hamil[2]; 
  Ghat[2] = 0.75*hamil[6]*fskin[22]+0.4330127018922193*hamil[6]*fskin[14]+0.75*hamil[7]*fskin[13]+0.75*hamil[2]*fskin[8]+0.4330127018922193*fskin[5]*hamil[7]+0.75*fskin[2]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[2]*fskin[3]; 
  Ghat[3] = 0.75*hamil[7]*fskin[30]+0.4330127018922193*hamil[7]*fskin[25]+0.75*hamil[6]*fskin[24]+0.75*hamil[4]*fskin[19]+0.4330127018922193*hamil[6]*fskin[15]+0.4330127018922193*hamil[4]*fskin[11]+0.75*hamil[2]*fskin[10]+0.4330127018922193*hamil[2]*fskin[4]; 
  Ghat[4] = 0.75*hamil[4]*fskin[22]+0.4330127018922193*hamil[4]*fskin[14]+0.75*hamil[2]*fskin[13]+0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[3]*hamil[7]+0.75*fskin[2]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[2]*fskin[5]; 
  Ghat[5] = 0.75*hamil[6]*fskin[27]+0.4330127018922193*hamil[6]*fskin[21]+0.75*hamil[7]*fskin[20]+0.75*hamil[2]*fskin[16]+0.4330127018922193*hamil[7]*fskin[12]+0.4330127018922193*hamil[2]*fskin[7]+0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[1]*hamil[4]; 
  Ghat[6] = 0.75*hamil[7]*fskin[31]+0.4330127018922193*hamil[7]*fskin[29]+0.75*hamil[6]*fskin[28]+0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[6]*fskin[23]+0.4330127018922193*hamil[4]*fskin[18]+0.75*hamil[2]*fskin[17]+0.4330127018922193*hamil[2]*fskin[9]; 
  Ghat[7] = 0.75*hamil[6]*fskin[30]+0.4330127018922193*hamil[6]*fskin[25]+0.75*hamil[7]*fskin[24]+0.75*hamil[2]*fskin[19]+0.4330127018922193*hamil[7]*fskin[15]+0.4330127018922193*hamil[2]*fskin[11]+0.75*hamil[4]*fskin[10]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[8] = 0.75*hamil[4]*fskin[27]+0.4330127018922193*hamil[4]*fskin[21]+0.75*hamil[2]*fskin[20]+0.75*hamil[7]*fskin[16]+0.4330127018922193*hamil[2]*fskin[12]+0.4330127018922193*fskin[7]*hamil[7]+0.75*fskin[6]*hamil[6]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[9] = 0.75*hamil[2]*fskin[22]+0.4330127018922193*hamil[2]*fskin[14]+0.75*hamil[4]*fskin[13]+0.75*hamil[6]*fskin[8]+0.75*fskin[2]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[3]*hamil[6]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[10] = 0.75*hamil[4]*fskin[30]+0.4330127018922193*hamil[4]*fskin[25]+0.75*hamil[2]*fskin[24]+0.75*hamil[7]*fskin[19]+0.4330127018922193*hamil[2]*fskin[15]+0.4330127018922193*hamil[7]*fskin[11]+0.75*hamil[6]*fskin[10]+0.4330127018922193*fskin[4]*hamil[6]; 
  Ghat[11] = 0.75*hamil[6]*fskin[31]+0.4330127018922193*hamil[6]*fskin[29]+0.75*hamil[7]*fskin[28]+0.75*hamil[2]*fskin[26]+0.4330127018922193*hamil[7]*fskin[23]+0.4330127018922193*hamil[2]*fskin[18]+0.75*hamil[4]*fskin[17]+0.4330127018922193*hamil[4]*fskin[9]; 
  Ghat[12] = 0.75*hamil[2]*fskin[27]+0.4330127018922193*hamil[2]*fskin[21]+0.75*hamil[4]*fskin[20]+0.75*hamil[6]*fskin[16]+0.4330127018922193*hamil[4]*fskin[12]+0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[1]*hamil[7]+0.4330127018922193*hamil[6]*fskin[7]; 
  Ghat[13] = 0.75*hamil[4]*fskin[31]+0.4330127018922193*hamil[4]*fskin[29]+0.75*hamil[2]*fskin[28]+0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[2]*fskin[23]+0.4330127018922193*hamil[7]*fskin[18]+0.75*hamil[6]*fskin[17]+0.4330127018922193*hamil[6]*fskin[9]; 
  Ghat[14] = 0.75*hamil[2]*fskin[30]+0.4330127018922193*hamil[2]*fskin[25]+0.75*hamil[4]*fskin[24]+0.75*hamil[6]*fskin[19]+0.4330127018922193*hamil[4]*fskin[15]+0.4330127018922193*hamil[6]*fskin[11]+0.75*hamil[7]*fskin[10]+0.4330127018922193*fskin[4]*hamil[7]; 
  Ghat[15] = 0.75*hamil[2]*fskin[31]+0.4330127018922193*hamil[2]*fskin[29]+0.75*hamil[4]*fskin[28]+0.75*hamil[6]*fskin[26]+0.4330127018922193*hamil[4]*fskin[23]+0.4330127018922193*hamil[6]*fskin[18]+0.75*hamil[7]*fskin[17]+0.4330127018922193*hamil[7]*fskin[9]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[22])+0.4330127018922193*hamil[7]*fedge[14]-0.75*hamil[6]*fedge[13]-0.75*hamil[4]*fedge[8]+0.4330127018922193*fedge[5]*hamil[6]+0.4330127018922193*fedge[3]*hamil[4]-0.75*fedge[2]*hamil[2]+0.4330127018922193*fedge[0]*hamil[2]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[27])+0.4330127018922193*hamil[7]*fedge[21]-0.75*hamil[6]*fedge[20]-0.75*hamil[4]*fedge[16]+0.4330127018922193*hamil[6]*fedge[12]+0.4330127018922193*hamil[4]*fedge[7]-0.75*hamil[2]*fedge[6]+0.4330127018922193*fedge[1]*hamil[2]; 
  Ghat[2] = -(0.75*hamil[6]*fedge[22])+0.4330127018922193*hamil[6]*fedge[14]-0.75*hamil[7]*fedge[13]-0.75*hamil[2]*fedge[8]+0.4330127018922193*fedge[5]*hamil[7]-0.75*fedge[2]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[2]*fedge[3]; 
  Ghat[3] = -(0.75*hamil[7]*fedge[30])+0.4330127018922193*hamil[7]*fedge[25]-0.75*hamil[6]*fedge[24]-0.75*hamil[4]*fedge[19]+0.4330127018922193*hamil[6]*fedge[15]+0.4330127018922193*hamil[4]*fedge[11]-0.75*hamil[2]*fedge[10]+0.4330127018922193*hamil[2]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[4]*fedge[22])+0.4330127018922193*hamil[4]*fedge[14]-0.75*hamil[2]*fedge[13]-0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[3]*hamil[7]-0.75*fedge[2]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[2]*fedge[5]; 
  Ghat[5] = -(0.75*hamil[6]*fedge[27])+0.4330127018922193*hamil[6]*fedge[21]-0.75*hamil[7]*fedge[20]-0.75*hamil[2]*fedge[16]+0.4330127018922193*hamil[7]*fedge[12]+0.4330127018922193*hamil[2]*fedge[7]-0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[1]*hamil[4]; 
  Ghat[6] = -(0.75*hamil[7]*fedge[31])+0.4330127018922193*hamil[7]*fedge[29]-0.75*hamil[6]*fedge[28]-0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[6]*fedge[23]+0.4330127018922193*hamil[4]*fedge[18]-0.75*hamil[2]*fedge[17]+0.4330127018922193*hamil[2]*fedge[9]; 
  Ghat[7] = -(0.75*hamil[6]*fedge[30])+0.4330127018922193*hamil[6]*fedge[25]-0.75*hamil[7]*fedge[24]-0.75*hamil[2]*fedge[19]+0.4330127018922193*hamil[7]*fedge[15]+0.4330127018922193*hamil[2]*fedge[11]-0.75*hamil[4]*fedge[10]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[4]*fedge[27])+0.4330127018922193*hamil[4]*fedge[21]-0.75*hamil[2]*fedge[20]-0.75*hamil[7]*fedge[16]+0.4330127018922193*hamil[2]*fedge[12]+0.4330127018922193*fedge[7]*hamil[7]-0.75*fedge[6]*hamil[6]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[9] = -(0.75*hamil[2]*fedge[22])+0.4330127018922193*hamil[2]*fedge[14]-0.75*hamil[4]*fedge[13]-0.75*hamil[6]*fedge[8]-0.75*fedge[2]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[3]*hamil[6]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[10] = -(0.75*hamil[4]*fedge[30])+0.4330127018922193*hamil[4]*fedge[25]-0.75*hamil[2]*fedge[24]-0.75*hamil[7]*fedge[19]+0.4330127018922193*hamil[2]*fedge[15]+0.4330127018922193*hamil[7]*fedge[11]-0.75*hamil[6]*fedge[10]+0.4330127018922193*fedge[4]*hamil[6]; 
  Ghat[11] = -(0.75*hamil[6]*fedge[31])+0.4330127018922193*hamil[6]*fedge[29]-0.75*hamil[7]*fedge[28]-0.75*hamil[2]*fedge[26]+0.4330127018922193*hamil[7]*fedge[23]+0.4330127018922193*hamil[2]*fedge[18]-0.75*hamil[4]*fedge[17]+0.4330127018922193*hamil[4]*fedge[9]; 
  Ghat[12] = -(0.75*hamil[2]*fedge[27])+0.4330127018922193*hamil[2]*fedge[21]-0.75*hamil[4]*fedge[20]-0.75*hamil[6]*fedge[16]+0.4330127018922193*hamil[4]*fedge[12]-0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[1]*hamil[7]+0.4330127018922193*hamil[6]*fedge[7]; 
  Ghat[13] = -(0.75*hamil[4]*fedge[31])+0.4330127018922193*hamil[4]*fedge[29]-0.75*hamil[2]*fedge[28]-0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[2]*fedge[23]+0.4330127018922193*hamil[7]*fedge[18]-0.75*hamil[6]*fedge[17]+0.4330127018922193*hamil[6]*fedge[9]; 
  Ghat[14] = -(0.75*hamil[2]*fedge[30])+0.4330127018922193*hamil[2]*fedge[25]-0.75*hamil[4]*fedge[24]-0.75*hamil[6]*fedge[19]+0.4330127018922193*hamil[4]*fedge[15]+0.4330127018922193*hamil[6]*fedge[11]-0.75*hamil[7]*fedge[10]+0.4330127018922193*fedge[4]*hamil[7]; 
  Ghat[15] = -(0.75*hamil[2]*fedge[31])+0.4330127018922193*hamil[2]*fedge[29]-0.75*hamil[4]*fedge[28]-0.75*hamil[6]*fedge[26]+0.4330127018922193*hamil[4]*fedge[23]+0.4330127018922193*hamil[6]*fedge[18]-0.75*hamil[7]*fedge[17]+0.4330127018922193*hamil[7]*fedge[9]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11*dx11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11*dx11); 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv11*dx11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11*dx11); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv11*dx11); 
  out[6] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[7] += -(0.7071067811865475*Ghat[5]*dv11*dx11); 
  out[8] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv11*dx11); 
  out[10] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv11*dx11); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv11*dx11); 
  out[13] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv11*dx11); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv11*dx11); 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[18] += -(0.7071067811865475*Ghat[11]*dv11*dx11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[21] += -(0.7071067811865475*Ghat[12]*dv11*dx11); 
  out[22] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[23] += -(0.7071067811865475*Ghat[13]*dv11*dx11); 
  out[24] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[25] += -(0.7071067811865475*Ghat[14]*dv11*dx11); 
  out[26] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[27] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[28] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[29] += -(0.7071067811865475*Ghat[15]*dv11*dx11); 
  out[30] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11*dx11); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[22]+0.4330127018922193*hamil[7]*fedge[14]+0.75*hamil[6]*fedge[13]+0.75*hamil[4]*fedge[8]+0.4330127018922193*fedge[5]*hamil[6]+0.4330127018922193*fedge[3]*hamil[4]+0.75*fedge[2]*hamil[2]+0.4330127018922193*fedge[0]*hamil[2]; 
  Ghat[1] = 0.75*hamil[7]*fedge[27]+0.4330127018922193*hamil[7]*fedge[21]+0.75*hamil[6]*fedge[20]+0.75*hamil[4]*fedge[16]+0.4330127018922193*hamil[6]*fedge[12]+0.4330127018922193*hamil[4]*fedge[7]+0.75*hamil[2]*fedge[6]+0.4330127018922193*fedge[1]*hamil[2]; 
  Ghat[2] = 0.75*hamil[6]*fedge[22]+0.4330127018922193*hamil[6]*fedge[14]+0.75*hamil[7]*fedge[13]+0.75*hamil[2]*fedge[8]+0.4330127018922193*fedge[5]*hamil[7]+0.75*fedge[2]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[2]*fedge[3]; 
  Ghat[3] = 0.75*hamil[7]*fedge[30]+0.4330127018922193*hamil[7]*fedge[25]+0.75*hamil[6]*fedge[24]+0.75*hamil[4]*fedge[19]+0.4330127018922193*hamil[6]*fedge[15]+0.4330127018922193*hamil[4]*fedge[11]+0.75*hamil[2]*fedge[10]+0.4330127018922193*hamil[2]*fedge[4]; 
  Ghat[4] = 0.75*hamil[4]*fedge[22]+0.4330127018922193*hamil[4]*fedge[14]+0.75*hamil[2]*fedge[13]+0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[3]*hamil[7]+0.75*fedge[2]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[2]*fedge[5]; 
  Ghat[5] = 0.75*hamil[6]*fedge[27]+0.4330127018922193*hamil[6]*fedge[21]+0.75*hamil[7]*fedge[20]+0.75*hamil[2]*fedge[16]+0.4330127018922193*hamil[7]*fedge[12]+0.4330127018922193*hamil[2]*fedge[7]+0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[1]*hamil[4]; 
  Ghat[6] = 0.75*hamil[7]*fedge[31]+0.4330127018922193*hamil[7]*fedge[29]+0.75*hamil[6]*fedge[28]+0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[6]*fedge[23]+0.4330127018922193*hamil[4]*fedge[18]+0.75*hamil[2]*fedge[17]+0.4330127018922193*hamil[2]*fedge[9]; 
  Ghat[7] = 0.75*hamil[6]*fedge[30]+0.4330127018922193*hamil[6]*fedge[25]+0.75*hamil[7]*fedge[24]+0.75*hamil[2]*fedge[19]+0.4330127018922193*hamil[7]*fedge[15]+0.4330127018922193*hamil[2]*fedge[11]+0.75*hamil[4]*fedge[10]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[8] = 0.75*hamil[4]*fedge[27]+0.4330127018922193*hamil[4]*fedge[21]+0.75*hamil[2]*fedge[20]+0.75*hamil[7]*fedge[16]+0.4330127018922193*hamil[2]*fedge[12]+0.4330127018922193*fedge[7]*hamil[7]+0.75*fedge[6]*hamil[6]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[9] = 0.75*hamil[2]*fedge[22]+0.4330127018922193*hamil[2]*fedge[14]+0.75*hamil[4]*fedge[13]+0.75*hamil[6]*fedge[8]+0.75*fedge[2]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[3]*hamil[6]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[10] = 0.75*hamil[4]*fedge[30]+0.4330127018922193*hamil[4]*fedge[25]+0.75*hamil[2]*fedge[24]+0.75*hamil[7]*fedge[19]+0.4330127018922193*hamil[2]*fedge[15]+0.4330127018922193*hamil[7]*fedge[11]+0.75*hamil[6]*fedge[10]+0.4330127018922193*fedge[4]*hamil[6]; 
  Ghat[11] = 0.75*hamil[6]*fedge[31]+0.4330127018922193*hamil[6]*fedge[29]+0.75*hamil[7]*fedge[28]+0.75*hamil[2]*fedge[26]+0.4330127018922193*hamil[7]*fedge[23]+0.4330127018922193*hamil[2]*fedge[18]+0.75*hamil[4]*fedge[17]+0.4330127018922193*hamil[4]*fedge[9]; 
  Ghat[12] = 0.75*hamil[2]*fedge[27]+0.4330127018922193*hamil[2]*fedge[21]+0.75*hamil[4]*fedge[20]+0.75*hamil[6]*fedge[16]+0.4330127018922193*hamil[4]*fedge[12]+0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[1]*hamil[7]+0.4330127018922193*hamil[6]*fedge[7]; 
  Ghat[13] = 0.75*hamil[4]*fedge[31]+0.4330127018922193*hamil[4]*fedge[29]+0.75*hamil[2]*fedge[28]+0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[2]*fedge[23]+0.4330127018922193*hamil[7]*fedge[18]+0.75*hamil[6]*fedge[17]+0.4330127018922193*hamil[6]*fedge[9]; 
  Ghat[14] = 0.75*hamil[2]*fedge[30]+0.4330127018922193*hamil[2]*fedge[25]+0.75*hamil[4]*fedge[24]+0.75*hamil[6]*fedge[19]+0.4330127018922193*hamil[4]*fedge[15]+0.4330127018922193*hamil[6]*fedge[11]+0.75*hamil[7]*fedge[10]+0.4330127018922193*fedge[4]*hamil[7]; 
  Ghat[15] = 0.75*hamil[2]*fedge[31]+0.4330127018922193*hamil[2]*fedge[29]+0.75*hamil[4]*fedge[28]+0.75*hamil[6]*fedge[26]+0.4330127018922193*hamil[4]*fedge[23]+0.4330127018922193*hamil[6]*fedge[18]+0.75*hamil[7]*fedge[17]+0.4330127018922193*hamil[7]*fedge[9]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[22])+0.4330127018922193*hamil[7]*fskin[14]-0.75*hamil[6]*fskin[13]-0.75*hamil[4]*fskin[8]+0.4330127018922193*fskin[5]*hamil[6]+0.4330127018922193*fskin[3]*hamil[4]-0.75*fskin[2]*hamil[2]+0.4330127018922193*fskin[0]*hamil[2]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[27])+0.4330127018922193*hamil[7]*fskin[21]-0.75*hamil[6]*fskin[20]-0.75*hamil[4]*fskin[16]+0.4330127018922193*hamil[6]*fskin[12]+0.4330127018922193*hamil[4]*fskin[7]-0.75*hamil[2]*fskin[6]+0.4330127018922193*fskin[1]*hamil[2]; 
  Ghat[2] = -(0.75*hamil[6]*fskin[22])+0.4330127018922193*hamil[6]*fskin[14]-0.75*hamil[7]*fskin[13]-0.75*hamil[2]*fskin[8]+0.4330127018922193*fskin[5]*hamil[7]-0.75*fskin[2]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[2]*fskin[3]; 
  Ghat[3] = -(0.75*hamil[7]*fskin[30])+0.4330127018922193*hamil[7]*fskin[25]-0.75*hamil[6]*fskin[24]-0.75*hamil[4]*fskin[19]+0.4330127018922193*hamil[6]*fskin[15]+0.4330127018922193*hamil[4]*fskin[11]-0.75*hamil[2]*fskin[10]+0.4330127018922193*hamil[2]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[4]*fskin[22])+0.4330127018922193*hamil[4]*fskin[14]-0.75*hamil[2]*fskin[13]-0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[3]*hamil[7]-0.75*fskin[2]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[2]*fskin[5]; 
  Ghat[5] = -(0.75*hamil[6]*fskin[27])+0.4330127018922193*hamil[6]*fskin[21]-0.75*hamil[7]*fskin[20]-0.75*hamil[2]*fskin[16]+0.4330127018922193*hamil[7]*fskin[12]+0.4330127018922193*hamil[2]*fskin[7]-0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[1]*hamil[4]; 
  Ghat[6] = -(0.75*hamil[7]*fskin[31])+0.4330127018922193*hamil[7]*fskin[29]-0.75*hamil[6]*fskin[28]-0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[6]*fskin[23]+0.4330127018922193*hamil[4]*fskin[18]-0.75*hamil[2]*fskin[17]+0.4330127018922193*hamil[2]*fskin[9]; 
  Ghat[7] = -(0.75*hamil[6]*fskin[30])+0.4330127018922193*hamil[6]*fskin[25]-0.75*hamil[7]*fskin[24]-0.75*hamil[2]*fskin[19]+0.4330127018922193*hamil[7]*fskin[15]+0.4330127018922193*hamil[2]*fskin[11]-0.75*hamil[4]*fskin[10]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[4]*fskin[27])+0.4330127018922193*hamil[4]*fskin[21]-0.75*hamil[2]*fskin[20]-0.75*hamil[7]*fskin[16]+0.4330127018922193*hamil[2]*fskin[12]+0.4330127018922193*fskin[7]*hamil[7]-0.75*fskin[6]*hamil[6]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[9] = -(0.75*hamil[2]*fskin[22])+0.4330127018922193*hamil[2]*fskin[14]-0.75*hamil[4]*fskin[13]-0.75*hamil[6]*fskin[8]-0.75*fskin[2]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[3]*hamil[6]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[10] = -(0.75*hamil[4]*fskin[30])+0.4330127018922193*hamil[4]*fskin[25]-0.75*hamil[2]*fskin[24]-0.75*hamil[7]*fskin[19]+0.4330127018922193*hamil[2]*fskin[15]+0.4330127018922193*hamil[7]*fskin[11]-0.75*hamil[6]*fskin[10]+0.4330127018922193*fskin[4]*hamil[6]; 
  Ghat[11] = -(0.75*hamil[6]*fskin[31])+0.4330127018922193*hamil[6]*fskin[29]-0.75*hamil[7]*fskin[28]-0.75*hamil[2]*fskin[26]+0.4330127018922193*hamil[7]*fskin[23]+0.4330127018922193*hamil[2]*fskin[18]-0.75*hamil[4]*fskin[17]+0.4330127018922193*hamil[4]*fskin[9]; 
  Ghat[12] = -(0.75*hamil[2]*fskin[27])+0.4330127018922193*hamil[2]*fskin[21]-0.75*hamil[4]*fskin[20]-0.75*hamil[6]*fskin[16]+0.4330127018922193*hamil[4]*fskin[12]-0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[1]*hamil[7]+0.4330127018922193*hamil[6]*fskin[7]; 
  Ghat[13] = -(0.75*hamil[4]*fskin[31])+0.4330127018922193*hamil[4]*fskin[29]-0.75*hamil[2]*fskin[28]-0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[2]*fskin[23]+0.4330127018922193*hamil[7]*fskin[18]-0.75*hamil[6]*fskin[17]+0.4330127018922193*hamil[6]*fskin[9]; 
  Ghat[14] = -(0.75*hamil[2]*fskin[30])+0.4330127018922193*hamil[2]*fskin[25]-0.75*hamil[4]*fskin[24]-0.75*hamil[6]*fskin[19]+0.4330127018922193*hamil[4]*fskin[15]+0.4330127018922193*hamil[6]*fskin[11]-0.75*hamil[7]*fskin[10]+0.4330127018922193*fskin[4]*hamil[7]; 
  Ghat[15] = -(0.75*hamil[2]*fskin[31])+0.4330127018922193*hamil[2]*fskin[29]-0.75*hamil[4]*fskin[28]-0.75*hamil[6]*fskin[26]+0.4330127018922193*hamil[4]*fskin[23]+0.4330127018922193*hamil[6]*fskin[18]-0.75*hamil[7]*fskin[17]+0.4330127018922193*hamil[7]*fskin[9]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv11*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11*dx11; 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += 0.7071067811865475*Ghat[2]*dv11*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11*dx11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11*dx11; 
  out[6] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[7] += 0.7071067811865475*Ghat[5]*dv11*dx11; 
  out[8] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[9] += 0.7071067811865475*Ghat[6]*dv11*dx11; 
  out[10] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[11] += 0.7071067811865475*Ghat[7]*dv11*dx11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11*dx11; 
  out[13] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[14] += 0.7071067811865475*Ghat[9]*dv11*dx11; 
  out[15] += 0.7071067811865475*Ghat[10]*dv11*dx11; 
  out[16] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[17] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[18] += 0.7071067811865475*Ghat[11]*dv11*dx11; 
  out[19] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[20] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[21] += 0.7071067811865475*Ghat[12]*dv11*dx11; 
  out[22] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[23] += 0.7071067811865475*Ghat[13]*dv11*dx11; 
  out[24] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[25] += 0.7071067811865475*Ghat[14]*dv11*dx11; 
  out[26] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[27] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[28] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[29] += 0.7071067811865475*Ghat[15]*dv11*dx11; 
  out[30] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11*dx11); 

  } 
  return 0.0;

} 
