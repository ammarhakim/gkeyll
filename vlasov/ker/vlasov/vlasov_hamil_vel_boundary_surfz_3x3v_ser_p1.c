#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx12 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[5]; 
  const double wv = w[5]; 

  double Ghat[32]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[31]+0.4330127018922193*hamil[7]*fskin[16]+0.75*hamil[6]*fskin[15]+0.75*hamil[5]*fskin[12]+0.4330127018922193*fskin[5]*hamil[6]+0.4330127018922193*fskin[4]*hamil[5]+0.75*fskin[3]*hamil[3]+0.4330127018922193*fskin[0]*hamil[3]; 
  Ghat[1] = 0.75*hamil[7]*fskin[45]+0.4330127018922193*hamil[7]*fskin[29]+0.75*hamil[6]*fskin[27]+0.75*hamil[5]*fskin[24]+0.4330127018922193*hamil[6]*fskin[13]+0.4330127018922193*hamil[5]*fskin[10]+0.75*hamil[3]*fskin[8]+0.4330127018922193*fskin[1]*hamil[3]; 
  Ghat[2] = 0.75*hamil[7]*fskin[46]+0.4330127018922193*hamil[7]*fskin[30]+0.75*hamil[6]*fskin[28]+0.75*hamil[5]*fskin[25]+0.4330127018922193*hamil[6]*fskin[14]+0.4330127018922193*hamil[5]*fskin[11]+0.75*hamil[3]*fskin[9]+0.4330127018922193*fskin[2]*hamil[3]; 
  Ghat[3] = 0.75*hamil[6]*fskin[31]+0.4330127018922193*hamil[6]*fskin[16]+0.75*hamil[7]*fskin[15]+0.75*hamil[3]*fskin[12]+0.4330127018922193*fskin[5]*hamil[7]+0.75*fskin[3]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[3]*fskin[4]; 
  Ghat[4] = 0.75*hamil[5]*fskin[31]+0.4330127018922193*hamil[5]*fskin[16]+0.75*hamil[3]*fskin[15]+0.75*hamil[7]*fskin[12]+0.4330127018922193*fskin[4]*hamil[7]+0.75*fskin[3]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[3]*fskin[5]; 
  Ghat[5] = 0.75*hamil[7]*fskin[56]+0.4330127018922193*hamil[7]*fskin[41]+0.75*hamil[6]*fskin[40]+0.75*hamil[5]*fskin[37]+0.4330127018922193*hamil[6]*fskin[21]+0.4330127018922193*hamil[5]*fskin[20]+0.75*hamil[3]*fskin[19]+0.4330127018922193*hamil[3]*fskin[6]; 
  Ghat[6] = 0.75*hamil[7]*fskin[57]+0.4330127018922193*hamil[7]*fskin[44]+0.75*hamil[6]*fskin[43]+0.75*hamil[5]*fskin[42]+0.4330127018922193*hamil[6]*fskin[26]+0.4330127018922193*hamil[5]*fskin[23]+0.75*hamil[3]*fskin[22]+0.4330127018922193*hamil[3]*fskin[7]; 
  Ghat[7] = 0.75*hamil[6]*fskin[45]+0.4330127018922193*hamil[6]*fskin[29]+0.75*hamil[7]*fskin[27]+0.75*hamil[3]*fskin[24]+0.4330127018922193*hamil[7]*fskin[13]+0.4330127018922193*hamil[3]*fskin[10]+0.75*hamil[5]*fskin[8]+0.4330127018922193*fskin[1]*hamil[5]; 
  Ghat[8] = 0.75*hamil[6]*fskin[46]+0.4330127018922193*hamil[6]*fskin[30]+0.75*hamil[7]*fskin[28]+0.75*hamil[3]*fskin[25]+0.4330127018922193*hamil[7]*fskin[14]+0.4330127018922193*hamil[3]*fskin[11]+0.75*hamil[5]*fskin[9]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[9] = 0.75*hamil[5]*fskin[45]+0.4330127018922193*hamil[5]*fskin[29]+0.75*hamil[3]*fskin[27]+0.75*hamil[7]*fskin[24]+0.4330127018922193*hamil[3]*fskin[13]+0.4330127018922193*hamil[7]*fskin[10]+0.75*hamil[6]*fskin[8]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[10] = 0.75*hamil[5]*fskin[46]+0.4330127018922193*hamil[5]*fskin[30]+0.75*hamil[3]*fskin[28]+0.75*hamil[7]*fskin[25]+0.4330127018922193*hamil[3]*fskin[14]+0.4330127018922193*hamil[7]*fskin[11]+0.75*hamil[6]*fskin[9]+0.4330127018922193*fskin[2]*hamil[6]; 
  Ghat[11] = 0.75*hamil[3]*fskin[31]+0.4330127018922193*hamil[3]*fskin[16]+0.75*hamil[5]*fskin[15]+0.75*hamil[6]*fskin[12]+0.75*fskin[3]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[6]+0.4330127018922193*fskin[5]*hamil[5]; 
  Ghat[12] = 0.75*hamil[7]*fskin[61]+0.4330127018922193*hamil[7]*fskin[54]+0.75*hamil[6]*fskin[52]+0.75*hamil[5]*fskin[49]+0.4330127018922193*hamil[6]*fskin[38]+0.4330127018922193*hamil[5]*fskin[35]+0.75*hamil[3]*fskin[33]+0.4330127018922193*hamil[3]*fskin[17]; 
  Ghat[13] = 0.75*hamil[7]*fskin[62]+0.4330127018922193*hamil[7]*fskin[55]+0.75*hamil[6]*fskin[53]+0.75*hamil[5]*fskin[50]+0.4330127018922193*hamil[6]*fskin[39]+0.4330127018922193*hamil[5]*fskin[36]+0.75*hamil[3]*fskin[34]+0.4330127018922193*hamil[3]*fskin[18]; 
  Ghat[14] = 0.75*hamil[6]*fskin[56]+0.4330127018922193*hamil[6]*fskin[41]+0.75*hamil[7]*fskin[40]+0.75*hamil[3]*fskin[37]+0.4330127018922193*hamil[7]*fskin[21]+0.4330127018922193*hamil[3]*fskin[20]+0.75*hamil[5]*fskin[19]+0.4330127018922193*hamil[5]*fskin[6]; 
  Ghat[15] = 0.75*hamil[5]*fskin[56]+0.4330127018922193*hamil[5]*fskin[41]+0.75*hamil[3]*fskin[40]+0.75*hamil[7]*fskin[37]+0.4330127018922193*hamil[3]*fskin[21]+0.4330127018922193*hamil[7]*fskin[20]+0.75*hamil[6]*fskin[19]+0.4330127018922193*fskin[6]*hamil[6]; 
  Ghat[16] = 0.75*hamil[6]*fskin[57]+0.4330127018922193*hamil[6]*fskin[44]+0.75*hamil[7]*fskin[43]+0.75*hamil[3]*fskin[42]+0.4330127018922193*hamil[7]*fskin[26]+0.4330127018922193*hamil[3]*fskin[23]+0.75*hamil[5]*fskin[22]+0.4330127018922193*hamil[5]*fskin[7]; 
  Ghat[17] = 0.75*hamil[5]*fskin[57]+0.4330127018922193*hamil[5]*fskin[44]+0.75*hamil[3]*fskin[43]+0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[3]*fskin[26]+0.4330127018922193*hamil[7]*fskin[23]+0.75*hamil[6]*fskin[22]+0.4330127018922193*hamil[6]*fskin[7]; 
  Ghat[18] = 0.75*hamil[3]*fskin[45]+0.4330127018922193*hamil[3]*fskin[29]+0.75*hamil[5]*fskin[27]+0.75*hamil[6]*fskin[24]+0.4330127018922193*hamil[5]*fskin[13]+0.4330127018922193*hamil[6]*fskin[10]+0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[1]*hamil[7]; 
  Ghat[19] = 0.75*hamil[3]*fskin[46]+0.4330127018922193*hamil[3]*fskin[30]+0.75*hamil[5]*fskin[28]+0.75*hamil[6]*fskin[25]+0.4330127018922193*hamil[5]*fskin[14]+0.4330127018922193*hamil[6]*fskin[11]+0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[20] = 0.75*hamil[7]*fskin[63]+0.4330127018922193*hamil[7]*fskin[60]+0.75*hamil[6]*fskin[59]+0.75*hamil[5]*fskin[58]+0.4330127018922193*hamil[6]*fskin[51]+0.4330127018922193*hamil[5]*fskin[48]+0.75*hamil[3]*fskin[47]+0.4330127018922193*hamil[3]*fskin[32]; 
  Ghat[21] = 0.75*hamil[6]*fskin[61]+0.4330127018922193*hamil[6]*fskin[54]+0.75*hamil[7]*fskin[52]+0.75*hamil[3]*fskin[49]+0.4330127018922193*hamil[7]*fskin[38]+0.4330127018922193*hamil[3]*fskin[35]+0.75*hamil[5]*fskin[33]+0.4330127018922193*hamil[5]*fskin[17]; 
  Ghat[22] = 0.75*hamil[6]*fskin[62]+0.4330127018922193*hamil[6]*fskin[55]+0.75*hamil[7]*fskin[53]+0.75*hamil[3]*fskin[50]+0.4330127018922193*hamil[7]*fskin[39]+0.4330127018922193*hamil[3]*fskin[36]+0.75*hamil[5]*fskin[34]+0.4330127018922193*hamil[5]*fskin[18]; 
  Ghat[23] = 0.75*hamil[5]*fskin[61]+0.4330127018922193*hamil[5]*fskin[54]+0.75*hamil[3]*fskin[52]+0.75*hamil[7]*fskin[49]+0.4330127018922193*hamil[3]*fskin[38]+0.4330127018922193*hamil[7]*fskin[35]+0.75*hamil[6]*fskin[33]+0.4330127018922193*hamil[6]*fskin[17]; 
  Ghat[24] = 0.75*hamil[5]*fskin[62]+0.4330127018922193*hamil[5]*fskin[55]+0.75*hamil[3]*fskin[53]+0.75*hamil[7]*fskin[50]+0.4330127018922193*hamil[3]*fskin[39]+0.4330127018922193*hamil[7]*fskin[36]+0.75*hamil[6]*fskin[34]+0.4330127018922193*hamil[6]*fskin[18]; 
  Ghat[25] = 0.75*hamil[3]*fskin[56]+0.4330127018922193*hamil[3]*fskin[41]+0.75*hamil[5]*fskin[40]+0.75*hamil[6]*fskin[37]+0.4330127018922193*hamil[5]*fskin[21]+0.4330127018922193*hamil[6]*fskin[20]+0.75*hamil[7]*fskin[19]+0.4330127018922193*fskin[6]*hamil[7]; 
  Ghat[26] = 0.75*hamil[3]*fskin[57]+0.4330127018922193*hamil[3]*fskin[44]+0.75*hamil[5]*fskin[43]+0.75*hamil[6]*fskin[42]+0.4330127018922193*hamil[5]*fskin[26]+0.4330127018922193*hamil[6]*fskin[23]+0.75*hamil[7]*fskin[22]+0.4330127018922193*fskin[7]*hamil[7]; 
  Ghat[27] = 0.75*hamil[6]*fskin[63]+0.4330127018922193*hamil[6]*fskin[60]+0.75*hamil[7]*fskin[59]+0.75*hamil[3]*fskin[58]+0.4330127018922193*hamil[7]*fskin[51]+0.4330127018922193*hamil[3]*fskin[48]+0.75*hamil[5]*fskin[47]+0.4330127018922193*hamil[5]*fskin[32]; 
  Ghat[28] = 0.75*hamil[5]*fskin[63]+0.4330127018922193*hamil[5]*fskin[60]+0.75*hamil[3]*fskin[59]+0.75*hamil[7]*fskin[58]+0.4330127018922193*hamil[3]*fskin[51]+0.4330127018922193*hamil[7]*fskin[48]+0.75*hamil[6]*fskin[47]+0.4330127018922193*hamil[6]*fskin[32]; 
  Ghat[29] = 0.75*hamil[3]*fskin[61]+0.4330127018922193*hamil[3]*fskin[54]+0.75*hamil[5]*fskin[52]+0.75*hamil[6]*fskin[49]+0.4330127018922193*hamil[5]*fskin[38]+0.4330127018922193*hamil[6]*fskin[35]+0.75*hamil[7]*fskin[33]+0.4330127018922193*hamil[7]*fskin[17]; 
  Ghat[30] = 0.75*hamil[3]*fskin[62]+0.4330127018922193*hamil[3]*fskin[55]+0.75*hamil[5]*fskin[53]+0.75*hamil[6]*fskin[50]+0.4330127018922193*hamil[5]*fskin[39]+0.4330127018922193*hamil[6]*fskin[36]+0.75*hamil[7]*fskin[34]+0.4330127018922193*hamil[7]*fskin[18]; 
  Ghat[31] = 0.75*hamil[3]*fskin[63]+0.4330127018922193*hamil[3]*fskin[60]+0.75*hamil[5]*fskin[59]+0.75*hamil[6]*fskin[58]+0.4330127018922193*hamil[5]*fskin[51]+0.4330127018922193*hamil[6]*fskin[48]+0.75*hamil[7]*fskin[47]+0.4330127018922193*hamil[7]*fskin[32]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[31])+0.4330127018922193*hamil[7]*fedge[16]-0.75*hamil[6]*fedge[15]-0.75*hamil[5]*fedge[12]+0.4330127018922193*fedge[5]*hamil[6]+0.4330127018922193*fedge[4]*hamil[5]-0.75*fedge[3]*hamil[3]+0.4330127018922193*fedge[0]*hamil[3]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[45])+0.4330127018922193*hamil[7]*fedge[29]-0.75*hamil[6]*fedge[27]-0.75*hamil[5]*fedge[24]+0.4330127018922193*hamil[6]*fedge[13]+0.4330127018922193*hamil[5]*fedge[10]-0.75*hamil[3]*fedge[8]+0.4330127018922193*fedge[1]*hamil[3]; 
  Ghat[2] = -(0.75*hamil[7]*fedge[46])+0.4330127018922193*hamil[7]*fedge[30]-0.75*hamil[6]*fedge[28]-0.75*hamil[5]*fedge[25]+0.4330127018922193*hamil[6]*fedge[14]+0.4330127018922193*hamil[5]*fedge[11]-0.75*hamil[3]*fedge[9]+0.4330127018922193*fedge[2]*hamil[3]; 
  Ghat[3] = -(0.75*hamil[6]*fedge[31])+0.4330127018922193*hamil[6]*fedge[16]-0.75*hamil[7]*fedge[15]-0.75*hamil[3]*fedge[12]+0.4330127018922193*fedge[5]*hamil[7]-0.75*fedge[3]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[3]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[5]*fedge[31])+0.4330127018922193*hamil[5]*fedge[16]-0.75*hamil[3]*fedge[15]-0.75*hamil[7]*fedge[12]+0.4330127018922193*fedge[4]*hamil[7]-0.75*fedge[3]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[3]*fedge[5]; 
  Ghat[5] = -(0.75*hamil[7]*fedge[56])+0.4330127018922193*hamil[7]*fedge[41]-0.75*hamil[6]*fedge[40]-0.75*hamil[5]*fedge[37]+0.4330127018922193*hamil[6]*fedge[21]+0.4330127018922193*hamil[5]*fedge[20]-0.75*hamil[3]*fedge[19]+0.4330127018922193*hamil[3]*fedge[6]; 
  Ghat[6] = -(0.75*hamil[7]*fedge[57])+0.4330127018922193*hamil[7]*fedge[44]-0.75*hamil[6]*fedge[43]-0.75*hamil[5]*fedge[42]+0.4330127018922193*hamil[6]*fedge[26]+0.4330127018922193*hamil[5]*fedge[23]-0.75*hamil[3]*fedge[22]+0.4330127018922193*hamil[3]*fedge[7]; 
  Ghat[7] = -(0.75*hamil[6]*fedge[45])+0.4330127018922193*hamil[6]*fedge[29]-0.75*hamil[7]*fedge[27]-0.75*hamil[3]*fedge[24]+0.4330127018922193*hamil[7]*fedge[13]+0.4330127018922193*hamil[3]*fedge[10]-0.75*hamil[5]*fedge[8]+0.4330127018922193*fedge[1]*hamil[5]; 
  Ghat[8] = -(0.75*hamil[6]*fedge[46])+0.4330127018922193*hamil[6]*fedge[30]-0.75*hamil[7]*fedge[28]-0.75*hamil[3]*fedge[25]+0.4330127018922193*hamil[7]*fedge[14]+0.4330127018922193*hamil[3]*fedge[11]-0.75*hamil[5]*fedge[9]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[9] = -(0.75*hamil[5]*fedge[45])+0.4330127018922193*hamil[5]*fedge[29]-0.75*hamil[3]*fedge[27]-0.75*hamil[7]*fedge[24]+0.4330127018922193*hamil[3]*fedge[13]+0.4330127018922193*hamil[7]*fedge[10]-0.75*hamil[6]*fedge[8]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[10] = -(0.75*hamil[5]*fedge[46])+0.4330127018922193*hamil[5]*fedge[30]-0.75*hamil[3]*fedge[28]-0.75*hamil[7]*fedge[25]+0.4330127018922193*hamil[3]*fedge[14]+0.4330127018922193*hamil[7]*fedge[11]-0.75*hamil[6]*fedge[9]+0.4330127018922193*fedge[2]*hamil[6]; 
  Ghat[11] = -(0.75*hamil[3]*fedge[31])+0.4330127018922193*hamil[3]*fedge[16]-0.75*hamil[5]*fedge[15]-0.75*hamil[6]*fedge[12]-0.75*fedge[3]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[6]+0.4330127018922193*fedge[5]*hamil[5]; 
  Ghat[12] = -(0.75*hamil[7]*fedge[61])+0.4330127018922193*hamil[7]*fedge[54]-0.75*hamil[6]*fedge[52]-0.75*hamil[5]*fedge[49]+0.4330127018922193*hamil[6]*fedge[38]+0.4330127018922193*hamil[5]*fedge[35]-0.75*hamil[3]*fedge[33]+0.4330127018922193*hamil[3]*fedge[17]; 
  Ghat[13] = -(0.75*hamil[7]*fedge[62])+0.4330127018922193*hamil[7]*fedge[55]-0.75*hamil[6]*fedge[53]-0.75*hamil[5]*fedge[50]+0.4330127018922193*hamil[6]*fedge[39]+0.4330127018922193*hamil[5]*fedge[36]-0.75*hamil[3]*fedge[34]+0.4330127018922193*hamil[3]*fedge[18]; 
  Ghat[14] = -(0.75*hamil[6]*fedge[56])+0.4330127018922193*hamil[6]*fedge[41]-0.75*hamil[7]*fedge[40]-0.75*hamil[3]*fedge[37]+0.4330127018922193*hamil[7]*fedge[21]+0.4330127018922193*hamil[3]*fedge[20]-0.75*hamil[5]*fedge[19]+0.4330127018922193*hamil[5]*fedge[6]; 
  Ghat[15] = -(0.75*hamil[5]*fedge[56])+0.4330127018922193*hamil[5]*fedge[41]-0.75*hamil[3]*fedge[40]-0.75*hamil[7]*fedge[37]+0.4330127018922193*hamil[3]*fedge[21]+0.4330127018922193*hamil[7]*fedge[20]-0.75*hamil[6]*fedge[19]+0.4330127018922193*fedge[6]*hamil[6]; 
  Ghat[16] = -(0.75*hamil[6]*fedge[57])+0.4330127018922193*hamil[6]*fedge[44]-0.75*hamil[7]*fedge[43]-0.75*hamil[3]*fedge[42]+0.4330127018922193*hamil[7]*fedge[26]+0.4330127018922193*hamil[3]*fedge[23]-0.75*hamil[5]*fedge[22]+0.4330127018922193*hamil[5]*fedge[7]; 
  Ghat[17] = -(0.75*hamil[5]*fedge[57])+0.4330127018922193*hamil[5]*fedge[44]-0.75*hamil[3]*fedge[43]-0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[3]*fedge[26]+0.4330127018922193*hamil[7]*fedge[23]-0.75*hamil[6]*fedge[22]+0.4330127018922193*hamil[6]*fedge[7]; 
  Ghat[18] = -(0.75*hamil[3]*fedge[45])+0.4330127018922193*hamil[3]*fedge[29]-0.75*hamil[5]*fedge[27]-0.75*hamil[6]*fedge[24]+0.4330127018922193*hamil[5]*fedge[13]+0.4330127018922193*hamil[6]*fedge[10]-0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[1]*hamil[7]; 
  Ghat[19] = -(0.75*hamil[3]*fedge[46])+0.4330127018922193*hamil[3]*fedge[30]-0.75*hamil[5]*fedge[28]-0.75*hamil[6]*fedge[25]+0.4330127018922193*hamil[5]*fedge[14]+0.4330127018922193*hamil[6]*fedge[11]-0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[20] = -(0.75*hamil[7]*fedge[63])+0.4330127018922193*hamil[7]*fedge[60]-0.75*hamil[6]*fedge[59]-0.75*hamil[5]*fedge[58]+0.4330127018922193*hamil[6]*fedge[51]+0.4330127018922193*hamil[5]*fedge[48]-0.75*hamil[3]*fedge[47]+0.4330127018922193*hamil[3]*fedge[32]; 
  Ghat[21] = -(0.75*hamil[6]*fedge[61])+0.4330127018922193*hamil[6]*fedge[54]-0.75*hamil[7]*fedge[52]-0.75*hamil[3]*fedge[49]+0.4330127018922193*hamil[7]*fedge[38]+0.4330127018922193*hamil[3]*fedge[35]-0.75*hamil[5]*fedge[33]+0.4330127018922193*hamil[5]*fedge[17]; 
  Ghat[22] = -(0.75*hamil[6]*fedge[62])+0.4330127018922193*hamil[6]*fedge[55]-0.75*hamil[7]*fedge[53]-0.75*hamil[3]*fedge[50]+0.4330127018922193*hamil[7]*fedge[39]+0.4330127018922193*hamil[3]*fedge[36]-0.75*hamil[5]*fedge[34]+0.4330127018922193*hamil[5]*fedge[18]; 
  Ghat[23] = -(0.75*hamil[5]*fedge[61])+0.4330127018922193*hamil[5]*fedge[54]-0.75*hamil[3]*fedge[52]-0.75*hamil[7]*fedge[49]+0.4330127018922193*hamil[3]*fedge[38]+0.4330127018922193*hamil[7]*fedge[35]-0.75*hamil[6]*fedge[33]+0.4330127018922193*hamil[6]*fedge[17]; 
  Ghat[24] = -(0.75*hamil[5]*fedge[62])+0.4330127018922193*hamil[5]*fedge[55]-0.75*hamil[3]*fedge[53]-0.75*hamil[7]*fedge[50]+0.4330127018922193*hamil[3]*fedge[39]+0.4330127018922193*hamil[7]*fedge[36]-0.75*hamil[6]*fedge[34]+0.4330127018922193*hamil[6]*fedge[18]; 
  Ghat[25] = -(0.75*hamil[3]*fedge[56])+0.4330127018922193*hamil[3]*fedge[41]-0.75*hamil[5]*fedge[40]-0.75*hamil[6]*fedge[37]+0.4330127018922193*hamil[5]*fedge[21]+0.4330127018922193*hamil[6]*fedge[20]-0.75*hamil[7]*fedge[19]+0.4330127018922193*fedge[6]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[3]*fedge[57])+0.4330127018922193*hamil[3]*fedge[44]-0.75*hamil[5]*fedge[43]-0.75*hamil[6]*fedge[42]+0.4330127018922193*hamil[5]*fedge[26]+0.4330127018922193*hamil[6]*fedge[23]-0.75*hamil[7]*fedge[22]+0.4330127018922193*fedge[7]*hamil[7]; 
  Ghat[27] = -(0.75*hamil[6]*fedge[63])+0.4330127018922193*hamil[6]*fedge[60]-0.75*hamil[7]*fedge[59]-0.75*hamil[3]*fedge[58]+0.4330127018922193*hamil[7]*fedge[51]+0.4330127018922193*hamil[3]*fedge[48]-0.75*hamil[5]*fedge[47]+0.4330127018922193*hamil[5]*fedge[32]; 
  Ghat[28] = -(0.75*hamil[5]*fedge[63])+0.4330127018922193*hamil[5]*fedge[60]-0.75*hamil[3]*fedge[59]-0.75*hamil[7]*fedge[58]+0.4330127018922193*hamil[3]*fedge[51]+0.4330127018922193*hamil[7]*fedge[48]-0.75*hamil[6]*fedge[47]+0.4330127018922193*hamil[6]*fedge[32]; 
  Ghat[29] = -(0.75*hamil[3]*fedge[61])+0.4330127018922193*hamil[3]*fedge[54]-0.75*hamil[5]*fedge[52]-0.75*hamil[6]*fedge[49]+0.4330127018922193*hamil[5]*fedge[38]+0.4330127018922193*hamil[6]*fedge[35]-0.75*hamil[7]*fedge[33]+0.4330127018922193*hamil[7]*fedge[17]; 
  Ghat[30] = -(0.75*hamil[3]*fedge[62])+0.4330127018922193*hamil[3]*fedge[55]-0.75*hamil[5]*fedge[53]-0.75*hamil[6]*fedge[50]+0.4330127018922193*hamil[5]*fedge[39]+0.4330127018922193*hamil[6]*fedge[36]-0.75*hamil[7]*fedge[34]+0.4330127018922193*hamil[7]*fedge[18]; 
  Ghat[31] = -(0.75*hamil[3]*fedge[63])+0.4330127018922193*hamil[3]*fedge[60]-0.75*hamil[5]*fedge[59]-0.75*hamil[6]*fedge[58]+0.4330127018922193*hamil[5]*fedge[51]+0.4330127018922193*hamil[6]*fedge[48]-0.75*hamil[7]*fedge[47]+0.4330127018922193*hamil[7]*fedge[32]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv12*dx12); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv12*dx12); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv12*dx12); 
  out[3] += -(1.224744871391589*Ghat[0]*dv12*dx12); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv12*dx12); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv12*dx12); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv12*dx12); 
  out[7] += -(0.7071067811865475*Ghat[6]*dv12*dx12); 
  out[8] += -(1.224744871391589*Ghat[1]*dv12*dx12); 
  out[9] += -(1.224744871391589*Ghat[2]*dv12*dx12); 
  out[10] += -(0.7071067811865475*Ghat[7]*dv12*dx12); 
  out[11] += -(0.7071067811865475*Ghat[8]*dv12*dx12); 
  out[12] += -(1.224744871391589*Ghat[3]*dv12*dx12); 
  out[13] += -(0.7071067811865475*Ghat[9]*dv12*dx12); 
  out[14] += -(0.7071067811865475*Ghat[10]*dv12*dx12); 
  out[15] += -(1.224744871391589*Ghat[4]*dv12*dx12); 
  out[16] += -(0.7071067811865475*Ghat[11]*dv12*dx12); 
  out[17] += -(0.7071067811865475*Ghat[12]*dv12*dx12); 
  out[18] += -(0.7071067811865475*Ghat[13]*dv12*dx12); 
  out[19] += -(1.224744871391589*Ghat[5]*dv12*dx12); 
  out[20] += -(0.7071067811865475*Ghat[14]*dv12*dx12); 
  out[21] += -(0.7071067811865475*Ghat[15]*dv12*dx12); 
  out[22] += -(1.224744871391589*Ghat[6]*dv12*dx12); 
  out[23] += -(0.7071067811865475*Ghat[16]*dv12*dx12); 
  out[24] += -(1.224744871391589*Ghat[7]*dv12*dx12); 
  out[25] += -(1.224744871391589*Ghat[8]*dv12*dx12); 
  out[26] += -(0.7071067811865475*Ghat[17]*dv12*dx12); 
  out[27] += -(1.224744871391589*Ghat[9]*dv12*dx12); 
  out[28] += -(1.224744871391589*Ghat[10]*dv12*dx12); 
  out[29] += -(0.7071067811865475*Ghat[18]*dv12*dx12); 
  out[30] += -(0.7071067811865475*Ghat[19]*dv12*dx12); 
  out[31] += -(1.224744871391589*Ghat[11]*dv12*dx12); 
  out[32] += -(0.7071067811865475*Ghat[20]*dv12*dx12); 
  out[33] += -(1.224744871391589*Ghat[12]*dv12*dx12); 
  out[34] += -(1.224744871391589*Ghat[13]*dv12*dx12); 
  out[35] += -(0.7071067811865475*Ghat[21]*dv12*dx12); 
  out[36] += -(0.7071067811865475*Ghat[22]*dv12*dx12); 
  out[37] += -(1.224744871391589*Ghat[14]*dv12*dx12); 
  out[38] += -(0.7071067811865475*Ghat[23]*dv12*dx12); 
  out[39] += -(0.7071067811865475*Ghat[24]*dv12*dx12); 
  out[40] += -(1.224744871391589*Ghat[15]*dv12*dx12); 
  out[41] += -(0.7071067811865475*Ghat[25]*dv12*dx12); 
  out[42] += -(1.224744871391589*Ghat[16]*dv12*dx12); 
  out[43] += -(1.224744871391589*Ghat[17]*dv12*dx12); 
  out[44] += -(0.7071067811865475*Ghat[26]*dv12*dx12); 
  out[45] += -(1.224744871391589*Ghat[18]*dv12*dx12); 
  out[46] += -(1.224744871391589*Ghat[19]*dv12*dx12); 
  out[47] += -(1.224744871391589*Ghat[20]*dv12*dx12); 
  out[48] += -(0.7071067811865475*Ghat[27]*dv12*dx12); 
  out[49] += -(1.224744871391589*Ghat[21]*dv12*dx12); 
  out[50] += -(1.224744871391589*Ghat[22]*dv12*dx12); 
  out[51] += -(0.7071067811865475*Ghat[28]*dv12*dx12); 
  out[52] += -(1.224744871391589*Ghat[23]*dv12*dx12); 
  out[53] += -(1.224744871391589*Ghat[24]*dv12*dx12); 
  out[54] += -(0.7071067811865475*Ghat[29]*dv12*dx12); 
  out[55] += -(0.7071067811865475*Ghat[30]*dv12*dx12); 
  out[56] += -(1.224744871391589*Ghat[25]*dv12*dx12); 
  out[57] += -(1.224744871391589*Ghat[26]*dv12*dx12); 
  out[58] += -(1.224744871391589*Ghat[27]*dv12*dx12); 
  out[59] += -(1.224744871391589*Ghat[28]*dv12*dx12); 
  out[60] += -(0.7071067811865475*Ghat[31]*dv12*dx12); 
  out[61] += -(1.224744871391589*Ghat[29]*dv12*dx12); 
  out[62] += -(1.224744871391589*Ghat[30]*dv12*dx12); 
  out[63] += -(1.224744871391589*Ghat[31]*dv12*dx12); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[31]+0.4330127018922193*hamil[7]*fedge[16]+0.75*hamil[6]*fedge[15]+0.75*hamil[5]*fedge[12]+0.4330127018922193*fedge[5]*hamil[6]+0.4330127018922193*fedge[4]*hamil[5]+0.75*fedge[3]*hamil[3]+0.4330127018922193*fedge[0]*hamil[3]; 
  Ghat[1] = 0.75*hamil[7]*fedge[45]+0.4330127018922193*hamil[7]*fedge[29]+0.75*hamil[6]*fedge[27]+0.75*hamil[5]*fedge[24]+0.4330127018922193*hamil[6]*fedge[13]+0.4330127018922193*hamil[5]*fedge[10]+0.75*hamil[3]*fedge[8]+0.4330127018922193*fedge[1]*hamil[3]; 
  Ghat[2] = 0.75*hamil[7]*fedge[46]+0.4330127018922193*hamil[7]*fedge[30]+0.75*hamil[6]*fedge[28]+0.75*hamil[5]*fedge[25]+0.4330127018922193*hamil[6]*fedge[14]+0.4330127018922193*hamil[5]*fedge[11]+0.75*hamil[3]*fedge[9]+0.4330127018922193*fedge[2]*hamil[3]; 
  Ghat[3] = 0.75*hamil[6]*fedge[31]+0.4330127018922193*hamil[6]*fedge[16]+0.75*hamil[7]*fedge[15]+0.75*hamil[3]*fedge[12]+0.4330127018922193*fedge[5]*hamil[7]+0.75*fedge[3]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[3]*fedge[4]; 
  Ghat[4] = 0.75*hamil[5]*fedge[31]+0.4330127018922193*hamil[5]*fedge[16]+0.75*hamil[3]*fedge[15]+0.75*hamil[7]*fedge[12]+0.4330127018922193*fedge[4]*hamil[7]+0.75*fedge[3]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[3]*fedge[5]; 
  Ghat[5] = 0.75*hamil[7]*fedge[56]+0.4330127018922193*hamil[7]*fedge[41]+0.75*hamil[6]*fedge[40]+0.75*hamil[5]*fedge[37]+0.4330127018922193*hamil[6]*fedge[21]+0.4330127018922193*hamil[5]*fedge[20]+0.75*hamil[3]*fedge[19]+0.4330127018922193*hamil[3]*fedge[6]; 
  Ghat[6] = 0.75*hamil[7]*fedge[57]+0.4330127018922193*hamil[7]*fedge[44]+0.75*hamil[6]*fedge[43]+0.75*hamil[5]*fedge[42]+0.4330127018922193*hamil[6]*fedge[26]+0.4330127018922193*hamil[5]*fedge[23]+0.75*hamil[3]*fedge[22]+0.4330127018922193*hamil[3]*fedge[7]; 
  Ghat[7] = 0.75*hamil[6]*fedge[45]+0.4330127018922193*hamil[6]*fedge[29]+0.75*hamil[7]*fedge[27]+0.75*hamil[3]*fedge[24]+0.4330127018922193*hamil[7]*fedge[13]+0.4330127018922193*hamil[3]*fedge[10]+0.75*hamil[5]*fedge[8]+0.4330127018922193*fedge[1]*hamil[5]; 
  Ghat[8] = 0.75*hamil[6]*fedge[46]+0.4330127018922193*hamil[6]*fedge[30]+0.75*hamil[7]*fedge[28]+0.75*hamil[3]*fedge[25]+0.4330127018922193*hamil[7]*fedge[14]+0.4330127018922193*hamil[3]*fedge[11]+0.75*hamil[5]*fedge[9]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[9] = 0.75*hamil[5]*fedge[45]+0.4330127018922193*hamil[5]*fedge[29]+0.75*hamil[3]*fedge[27]+0.75*hamil[7]*fedge[24]+0.4330127018922193*hamil[3]*fedge[13]+0.4330127018922193*hamil[7]*fedge[10]+0.75*hamil[6]*fedge[8]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[10] = 0.75*hamil[5]*fedge[46]+0.4330127018922193*hamil[5]*fedge[30]+0.75*hamil[3]*fedge[28]+0.75*hamil[7]*fedge[25]+0.4330127018922193*hamil[3]*fedge[14]+0.4330127018922193*hamil[7]*fedge[11]+0.75*hamil[6]*fedge[9]+0.4330127018922193*fedge[2]*hamil[6]; 
  Ghat[11] = 0.75*hamil[3]*fedge[31]+0.4330127018922193*hamil[3]*fedge[16]+0.75*hamil[5]*fedge[15]+0.75*hamil[6]*fedge[12]+0.75*fedge[3]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[6]+0.4330127018922193*fedge[5]*hamil[5]; 
  Ghat[12] = 0.75*hamil[7]*fedge[61]+0.4330127018922193*hamil[7]*fedge[54]+0.75*hamil[6]*fedge[52]+0.75*hamil[5]*fedge[49]+0.4330127018922193*hamil[6]*fedge[38]+0.4330127018922193*hamil[5]*fedge[35]+0.75*hamil[3]*fedge[33]+0.4330127018922193*hamil[3]*fedge[17]; 
  Ghat[13] = 0.75*hamil[7]*fedge[62]+0.4330127018922193*hamil[7]*fedge[55]+0.75*hamil[6]*fedge[53]+0.75*hamil[5]*fedge[50]+0.4330127018922193*hamil[6]*fedge[39]+0.4330127018922193*hamil[5]*fedge[36]+0.75*hamil[3]*fedge[34]+0.4330127018922193*hamil[3]*fedge[18]; 
  Ghat[14] = 0.75*hamil[6]*fedge[56]+0.4330127018922193*hamil[6]*fedge[41]+0.75*hamil[7]*fedge[40]+0.75*hamil[3]*fedge[37]+0.4330127018922193*hamil[7]*fedge[21]+0.4330127018922193*hamil[3]*fedge[20]+0.75*hamil[5]*fedge[19]+0.4330127018922193*hamil[5]*fedge[6]; 
  Ghat[15] = 0.75*hamil[5]*fedge[56]+0.4330127018922193*hamil[5]*fedge[41]+0.75*hamil[3]*fedge[40]+0.75*hamil[7]*fedge[37]+0.4330127018922193*hamil[3]*fedge[21]+0.4330127018922193*hamil[7]*fedge[20]+0.75*hamil[6]*fedge[19]+0.4330127018922193*fedge[6]*hamil[6]; 
  Ghat[16] = 0.75*hamil[6]*fedge[57]+0.4330127018922193*hamil[6]*fedge[44]+0.75*hamil[7]*fedge[43]+0.75*hamil[3]*fedge[42]+0.4330127018922193*hamil[7]*fedge[26]+0.4330127018922193*hamil[3]*fedge[23]+0.75*hamil[5]*fedge[22]+0.4330127018922193*hamil[5]*fedge[7]; 
  Ghat[17] = 0.75*hamil[5]*fedge[57]+0.4330127018922193*hamil[5]*fedge[44]+0.75*hamil[3]*fedge[43]+0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[3]*fedge[26]+0.4330127018922193*hamil[7]*fedge[23]+0.75*hamil[6]*fedge[22]+0.4330127018922193*hamil[6]*fedge[7]; 
  Ghat[18] = 0.75*hamil[3]*fedge[45]+0.4330127018922193*hamil[3]*fedge[29]+0.75*hamil[5]*fedge[27]+0.75*hamil[6]*fedge[24]+0.4330127018922193*hamil[5]*fedge[13]+0.4330127018922193*hamil[6]*fedge[10]+0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[1]*hamil[7]; 
  Ghat[19] = 0.75*hamil[3]*fedge[46]+0.4330127018922193*hamil[3]*fedge[30]+0.75*hamil[5]*fedge[28]+0.75*hamil[6]*fedge[25]+0.4330127018922193*hamil[5]*fedge[14]+0.4330127018922193*hamil[6]*fedge[11]+0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[20] = 0.75*hamil[7]*fedge[63]+0.4330127018922193*hamil[7]*fedge[60]+0.75*hamil[6]*fedge[59]+0.75*hamil[5]*fedge[58]+0.4330127018922193*hamil[6]*fedge[51]+0.4330127018922193*hamil[5]*fedge[48]+0.75*hamil[3]*fedge[47]+0.4330127018922193*hamil[3]*fedge[32]; 
  Ghat[21] = 0.75*hamil[6]*fedge[61]+0.4330127018922193*hamil[6]*fedge[54]+0.75*hamil[7]*fedge[52]+0.75*hamil[3]*fedge[49]+0.4330127018922193*hamil[7]*fedge[38]+0.4330127018922193*hamil[3]*fedge[35]+0.75*hamil[5]*fedge[33]+0.4330127018922193*hamil[5]*fedge[17]; 
  Ghat[22] = 0.75*hamil[6]*fedge[62]+0.4330127018922193*hamil[6]*fedge[55]+0.75*hamil[7]*fedge[53]+0.75*hamil[3]*fedge[50]+0.4330127018922193*hamil[7]*fedge[39]+0.4330127018922193*hamil[3]*fedge[36]+0.75*hamil[5]*fedge[34]+0.4330127018922193*hamil[5]*fedge[18]; 
  Ghat[23] = 0.75*hamil[5]*fedge[61]+0.4330127018922193*hamil[5]*fedge[54]+0.75*hamil[3]*fedge[52]+0.75*hamil[7]*fedge[49]+0.4330127018922193*hamil[3]*fedge[38]+0.4330127018922193*hamil[7]*fedge[35]+0.75*hamil[6]*fedge[33]+0.4330127018922193*hamil[6]*fedge[17]; 
  Ghat[24] = 0.75*hamil[5]*fedge[62]+0.4330127018922193*hamil[5]*fedge[55]+0.75*hamil[3]*fedge[53]+0.75*hamil[7]*fedge[50]+0.4330127018922193*hamil[3]*fedge[39]+0.4330127018922193*hamil[7]*fedge[36]+0.75*hamil[6]*fedge[34]+0.4330127018922193*hamil[6]*fedge[18]; 
  Ghat[25] = 0.75*hamil[3]*fedge[56]+0.4330127018922193*hamil[3]*fedge[41]+0.75*hamil[5]*fedge[40]+0.75*hamil[6]*fedge[37]+0.4330127018922193*hamil[5]*fedge[21]+0.4330127018922193*hamil[6]*fedge[20]+0.75*hamil[7]*fedge[19]+0.4330127018922193*fedge[6]*hamil[7]; 
  Ghat[26] = 0.75*hamil[3]*fedge[57]+0.4330127018922193*hamil[3]*fedge[44]+0.75*hamil[5]*fedge[43]+0.75*hamil[6]*fedge[42]+0.4330127018922193*hamil[5]*fedge[26]+0.4330127018922193*hamil[6]*fedge[23]+0.75*hamil[7]*fedge[22]+0.4330127018922193*fedge[7]*hamil[7]; 
  Ghat[27] = 0.75*hamil[6]*fedge[63]+0.4330127018922193*hamil[6]*fedge[60]+0.75*hamil[7]*fedge[59]+0.75*hamil[3]*fedge[58]+0.4330127018922193*hamil[7]*fedge[51]+0.4330127018922193*hamil[3]*fedge[48]+0.75*hamil[5]*fedge[47]+0.4330127018922193*hamil[5]*fedge[32]; 
  Ghat[28] = 0.75*hamil[5]*fedge[63]+0.4330127018922193*hamil[5]*fedge[60]+0.75*hamil[3]*fedge[59]+0.75*hamil[7]*fedge[58]+0.4330127018922193*hamil[3]*fedge[51]+0.4330127018922193*hamil[7]*fedge[48]+0.75*hamil[6]*fedge[47]+0.4330127018922193*hamil[6]*fedge[32]; 
  Ghat[29] = 0.75*hamil[3]*fedge[61]+0.4330127018922193*hamil[3]*fedge[54]+0.75*hamil[5]*fedge[52]+0.75*hamil[6]*fedge[49]+0.4330127018922193*hamil[5]*fedge[38]+0.4330127018922193*hamil[6]*fedge[35]+0.75*hamil[7]*fedge[33]+0.4330127018922193*hamil[7]*fedge[17]; 
  Ghat[30] = 0.75*hamil[3]*fedge[62]+0.4330127018922193*hamil[3]*fedge[55]+0.75*hamil[5]*fedge[53]+0.75*hamil[6]*fedge[50]+0.4330127018922193*hamil[5]*fedge[39]+0.4330127018922193*hamil[6]*fedge[36]+0.75*hamil[7]*fedge[34]+0.4330127018922193*hamil[7]*fedge[18]; 
  Ghat[31] = 0.75*hamil[3]*fedge[63]+0.4330127018922193*hamil[3]*fedge[60]+0.75*hamil[5]*fedge[59]+0.75*hamil[6]*fedge[58]+0.4330127018922193*hamil[5]*fedge[51]+0.4330127018922193*hamil[6]*fedge[48]+0.75*hamil[7]*fedge[47]+0.4330127018922193*hamil[7]*fedge[32]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[31])+0.4330127018922193*hamil[7]*fskin[16]-0.75*hamil[6]*fskin[15]-0.75*hamil[5]*fskin[12]+0.4330127018922193*fskin[5]*hamil[6]+0.4330127018922193*fskin[4]*hamil[5]-0.75*fskin[3]*hamil[3]+0.4330127018922193*fskin[0]*hamil[3]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[45])+0.4330127018922193*hamil[7]*fskin[29]-0.75*hamil[6]*fskin[27]-0.75*hamil[5]*fskin[24]+0.4330127018922193*hamil[6]*fskin[13]+0.4330127018922193*hamil[5]*fskin[10]-0.75*hamil[3]*fskin[8]+0.4330127018922193*fskin[1]*hamil[3]; 
  Ghat[2] = -(0.75*hamil[7]*fskin[46])+0.4330127018922193*hamil[7]*fskin[30]-0.75*hamil[6]*fskin[28]-0.75*hamil[5]*fskin[25]+0.4330127018922193*hamil[6]*fskin[14]+0.4330127018922193*hamil[5]*fskin[11]-0.75*hamil[3]*fskin[9]+0.4330127018922193*fskin[2]*hamil[3]; 
  Ghat[3] = -(0.75*hamil[6]*fskin[31])+0.4330127018922193*hamil[6]*fskin[16]-0.75*hamil[7]*fskin[15]-0.75*hamil[3]*fskin[12]+0.4330127018922193*fskin[5]*hamil[7]-0.75*fskin[3]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[3]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[5]*fskin[31])+0.4330127018922193*hamil[5]*fskin[16]-0.75*hamil[3]*fskin[15]-0.75*hamil[7]*fskin[12]+0.4330127018922193*fskin[4]*hamil[7]-0.75*fskin[3]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[3]*fskin[5]; 
  Ghat[5] = -(0.75*hamil[7]*fskin[56])+0.4330127018922193*hamil[7]*fskin[41]-0.75*hamil[6]*fskin[40]-0.75*hamil[5]*fskin[37]+0.4330127018922193*hamil[6]*fskin[21]+0.4330127018922193*hamil[5]*fskin[20]-0.75*hamil[3]*fskin[19]+0.4330127018922193*hamil[3]*fskin[6]; 
  Ghat[6] = -(0.75*hamil[7]*fskin[57])+0.4330127018922193*hamil[7]*fskin[44]-0.75*hamil[6]*fskin[43]-0.75*hamil[5]*fskin[42]+0.4330127018922193*hamil[6]*fskin[26]+0.4330127018922193*hamil[5]*fskin[23]-0.75*hamil[3]*fskin[22]+0.4330127018922193*hamil[3]*fskin[7]; 
  Ghat[7] = -(0.75*hamil[6]*fskin[45])+0.4330127018922193*hamil[6]*fskin[29]-0.75*hamil[7]*fskin[27]-0.75*hamil[3]*fskin[24]+0.4330127018922193*hamil[7]*fskin[13]+0.4330127018922193*hamil[3]*fskin[10]-0.75*hamil[5]*fskin[8]+0.4330127018922193*fskin[1]*hamil[5]; 
  Ghat[8] = -(0.75*hamil[6]*fskin[46])+0.4330127018922193*hamil[6]*fskin[30]-0.75*hamil[7]*fskin[28]-0.75*hamil[3]*fskin[25]+0.4330127018922193*hamil[7]*fskin[14]+0.4330127018922193*hamil[3]*fskin[11]-0.75*hamil[5]*fskin[9]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[9] = -(0.75*hamil[5]*fskin[45])+0.4330127018922193*hamil[5]*fskin[29]-0.75*hamil[3]*fskin[27]-0.75*hamil[7]*fskin[24]+0.4330127018922193*hamil[3]*fskin[13]+0.4330127018922193*hamil[7]*fskin[10]-0.75*hamil[6]*fskin[8]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[10] = -(0.75*hamil[5]*fskin[46])+0.4330127018922193*hamil[5]*fskin[30]-0.75*hamil[3]*fskin[28]-0.75*hamil[7]*fskin[25]+0.4330127018922193*hamil[3]*fskin[14]+0.4330127018922193*hamil[7]*fskin[11]-0.75*hamil[6]*fskin[9]+0.4330127018922193*fskin[2]*hamil[6]; 
  Ghat[11] = -(0.75*hamil[3]*fskin[31])+0.4330127018922193*hamil[3]*fskin[16]-0.75*hamil[5]*fskin[15]-0.75*hamil[6]*fskin[12]-0.75*fskin[3]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[6]+0.4330127018922193*fskin[5]*hamil[5]; 
  Ghat[12] = -(0.75*hamil[7]*fskin[61])+0.4330127018922193*hamil[7]*fskin[54]-0.75*hamil[6]*fskin[52]-0.75*hamil[5]*fskin[49]+0.4330127018922193*hamil[6]*fskin[38]+0.4330127018922193*hamil[5]*fskin[35]-0.75*hamil[3]*fskin[33]+0.4330127018922193*hamil[3]*fskin[17]; 
  Ghat[13] = -(0.75*hamil[7]*fskin[62])+0.4330127018922193*hamil[7]*fskin[55]-0.75*hamil[6]*fskin[53]-0.75*hamil[5]*fskin[50]+0.4330127018922193*hamil[6]*fskin[39]+0.4330127018922193*hamil[5]*fskin[36]-0.75*hamil[3]*fskin[34]+0.4330127018922193*hamil[3]*fskin[18]; 
  Ghat[14] = -(0.75*hamil[6]*fskin[56])+0.4330127018922193*hamil[6]*fskin[41]-0.75*hamil[7]*fskin[40]-0.75*hamil[3]*fskin[37]+0.4330127018922193*hamil[7]*fskin[21]+0.4330127018922193*hamil[3]*fskin[20]-0.75*hamil[5]*fskin[19]+0.4330127018922193*hamil[5]*fskin[6]; 
  Ghat[15] = -(0.75*hamil[5]*fskin[56])+0.4330127018922193*hamil[5]*fskin[41]-0.75*hamil[3]*fskin[40]-0.75*hamil[7]*fskin[37]+0.4330127018922193*hamil[3]*fskin[21]+0.4330127018922193*hamil[7]*fskin[20]-0.75*hamil[6]*fskin[19]+0.4330127018922193*fskin[6]*hamil[6]; 
  Ghat[16] = -(0.75*hamil[6]*fskin[57])+0.4330127018922193*hamil[6]*fskin[44]-0.75*hamil[7]*fskin[43]-0.75*hamil[3]*fskin[42]+0.4330127018922193*hamil[7]*fskin[26]+0.4330127018922193*hamil[3]*fskin[23]-0.75*hamil[5]*fskin[22]+0.4330127018922193*hamil[5]*fskin[7]; 
  Ghat[17] = -(0.75*hamil[5]*fskin[57])+0.4330127018922193*hamil[5]*fskin[44]-0.75*hamil[3]*fskin[43]-0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[3]*fskin[26]+0.4330127018922193*hamil[7]*fskin[23]-0.75*hamil[6]*fskin[22]+0.4330127018922193*hamil[6]*fskin[7]; 
  Ghat[18] = -(0.75*hamil[3]*fskin[45])+0.4330127018922193*hamil[3]*fskin[29]-0.75*hamil[5]*fskin[27]-0.75*hamil[6]*fskin[24]+0.4330127018922193*hamil[5]*fskin[13]+0.4330127018922193*hamil[6]*fskin[10]-0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[1]*hamil[7]; 
  Ghat[19] = -(0.75*hamil[3]*fskin[46])+0.4330127018922193*hamil[3]*fskin[30]-0.75*hamil[5]*fskin[28]-0.75*hamil[6]*fskin[25]+0.4330127018922193*hamil[5]*fskin[14]+0.4330127018922193*hamil[6]*fskin[11]-0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[20] = -(0.75*hamil[7]*fskin[63])+0.4330127018922193*hamil[7]*fskin[60]-0.75*hamil[6]*fskin[59]-0.75*hamil[5]*fskin[58]+0.4330127018922193*hamil[6]*fskin[51]+0.4330127018922193*hamil[5]*fskin[48]-0.75*hamil[3]*fskin[47]+0.4330127018922193*hamil[3]*fskin[32]; 
  Ghat[21] = -(0.75*hamil[6]*fskin[61])+0.4330127018922193*hamil[6]*fskin[54]-0.75*hamil[7]*fskin[52]-0.75*hamil[3]*fskin[49]+0.4330127018922193*hamil[7]*fskin[38]+0.4330127018922193*hamil[3]*fskin[35]-0.75*hamil[5]*fskin[33]+0.4330127018922193*hamil[5]*fskin[17]; 
  Ghat[22] = -(0.75*hamil[6]*fskin[62])+0.4330127018922193*hamil[6]*fskin[55]-0.75*hamil[7]*fskin[53]-0.75*hamil[3]*fskin[50]+0.4330127018922193*hamil[7]*fskin[39]+0.4330127018922193*hamil[3]*fskin[36]-0.75*hamil[5]*fskin[34]+0.4330127018922193*hamil[5]*fskin[18]; 
  Ghat[23] = -(0.75*hamil[5]*fskin[61])+0.4330127018922193*hamil[5]*fskin[54]-0.75*hamil[3]*fskin[52]-0.75*hamil[7]*fskin[49]+0.4330127018922193*hamil[3]*fskin[38]+0.4330127018922193*hamil[7]*fskin[35]-0.75*hamil[6]*fskin[33]+0.4330127018922193*hamil[6]*fskin[17]; 
  Ghat[24] = -(0.75*hamil[5]*fskin[62])+0.4330127018922193*hamil[5]*fskin[55]-0.75*hamil[3]*fskin[53]-0.75*hamil[7]*fskin[50]+0.4330127018922193*hamil[3]*fskin[39]+0.4330127018922193*hamil[7]*fskin[36]-0.75*hamil[6]*fskin[34]+0.4330127018922193*hamil[6]*fskin[18]; 
  Ghat[25] = -(0.75*hamil[3]*fskin[56])+0.4330127018922193*hamil[3]*fskin[41]-0.75*hamil[5]*fskin[40]-0.75*hamil[6]*fskin[37]+0.4330127018922193*hamil[5]*fskin[21]+0.4330127018922193*hamil[6]*fskin[20]-0.75*hamil[7]*fskin[19]+0.4330127018922193*fskin[6]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[3]*fskin[57])+0.4330127018922193*hamil[3]*fskin[44]-0.75*hamil[5]*fskin[43]-0.75*hamil[6]*fskin[42]+0.4330127018922193*hamil[5]*fskin[26]+0.4330127018922193*hamil[6]*fskin[23]-0.75*hamil[7]*fskin[22]+0.4330127018922193*fskin[7]*hamil[7]; 
  Ghat[27] = -(0.75*hamil[6]*fskin[63])+0.4330127018922193*hamil[6]*fskin[60]-0.75*hamil[7]*fskin[59]-0.75*hamil[3]*fskin[58]+0.4330127018922193*hamil[7]*fskin[51]+0.4330127018922193*hamil[3]*fskin[48]-0.75*hamil[5]*fskin[47]+0.4330127018922193*hamil[5]*fskin[32]; 
  Ghat[28] = -(0.75*hamil[5]*fskin[63])+0.4330127018922193*hamil[5]*fskin[60]-0.75*hamil[3]*fskin[59]-0.75*hamil[7]*fskin[58]+0.4330127018922193*hamil[3]*fskin[51]+0.4330127018922193*hamil[7]*fskin[48]-0.75*hamil[6]*fskin[47]+0.4330127018922193*hamil[6]*fskin[32]; 
  Ghat[29] = -(0.75*hamil[3]*fskin[61])+0.4330127018922193*hamil[3]*fskin[54]-0.75*hamil[5]*fskin[52]-0.75*hamil[6]*fskin[49]+0.4330127018922193*hamil[5]*fskin[38]+0.4330127018922193*hamil[6]*fskin[35]-0.75*hamil[7]*fskin[33]+0.4330127018922193*hamil[7]*fskin[17]; 
  Ghat[30] = -(0.75*hamil[3]*fskin[62])+0.4330127018922193*hamil[3]*fskin[55]-0.75*hamil[5]*fskin[53]-0.75*hamil[6]*fskin[50]+0.4330127018922193*hamil[5]*fskin[39]+0.4330127018922193*hamil[6]*fskin[36]-0.75*hamil[7]*fskin[34]+0.4330127018922193*hamil[7]*fskin[18]; 
  Ghat[31] = -(0.75*hamil[3]*fskin[63])+0.4330127018922193*hamil[3]*fskin[60]-0.75*hamil[5]*fskin[59]-0.75*hamil[6]*fskin[58]+0.4330127018922193*hamil[5]*fskin[51]+0.4330127018922193*hamil[6]*fskin[48]-0.75*hamil[7]*fskin[47]+0.4330127018922193*hamil[7]*fskin[32]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv12*dx12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12*dx12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12*dx12; 
  out[3] += -(1.224744871391589*Ghat[0]*dv12*dx12); 
  out[4] += 0.7071067811865475*Ghat[3]*dv12*dx12; 
  out[5] += 0.7071067811865475*Ghat[4]*dv12*dx12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12*dx12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12*dx12; 
  out[8] += -(1.224744871391589*Ghat[1]*dv12*dx12); 
  out[9] += -(1.224744871391589*Ghat[2]*dv12*dx12); 
  out[10] += 0.7071067811865475*Ghat[7]*dv12*dx12; 
  out[11] += 0.7071067811865475*Ghat[8]*dv12*dx12; 
  out[12] += -(1.224744871391589*Ghat[3]*dv12*dx12); 
  out[13] += 0.7071067811865475*Ghat[9]*dv12*dx12; 
  out[14] += 0.7071067811865475*Ghat[10]*dv12*dx12; 
  out[15] += -(1.224744871391589*Ghat[4]*dv12*dx12); 
  out[16] += 0.7071067811865475*Ghat[11]*dv12*dx12; 
  out[17] += 0.7071067811865475*Ghat[12]*dv12*dx12; 
  out[18] += 0.7071067811865475*Ghat[13]*dv12*dx12; 
  out[19] += -(1.224744871391589*Ghat[5]*dv12*dx12); 
  out[20] += 0.7071067811865475*Ghat[14]*dv12*dx12; 
  out[21] += 0.7071067811865475*Ghat[15]*dv12*dx12; 
  out[22] += -(1.224744871391589*Ghat[6]*dv12*dx12); 
  out[23] += 0.7071067811865475*Ghat[16]*dv12*dx12; 
  out[24] += -(1.224744871391589*Ghat[7]*dv12*dx12); 
  out[25] += -(1.224744871391589*Ghat[8]*dv12*dx12); 
  out[26] += 0.7071067811865475*Ghat[17]*dv12*dx12; 
  out[27] += -(1.224744871391589*Ghat[9]*dv12*dx12); 
  out[28] += -(1.224744871391589*Ghat[10]*dv12*dx12); 
  out[29] += 0.7071067811865475*Ghat[18]*dv12*dx12; 
  out[30] += 0.7071067811865475*Ghat[19]*dv12*dx12; 
  out[31] += -(1.224744871391589*Ghat[11]*dv12*dx12); 
  out[32] += 0.7071067811865475*Ghat[20]*dv12*dx12; 
  out[33] += -(1.224744871391589*Ghat[12]*dv12*dx12); 
  out[34] += -(1.224744871391589*Ghat[13]*dv12*dx12); 
  out[35] += 0.7071067811865475*Ghat[21]*dv12*dx12; 
  out[36] += 0.7071067811865475*Ghat[22]*dv12*dx12; 
  out[37] += -(1.224744871391589*Ghat[14]*dv12*dx12); 
  out[38] += 0.7071067811865475*Ghat[23]*dv12*dx12; 
  out[39] += 0.7071067811865475*Ghat[24]*dv12*dx12; 
  out[40] += -(1.224744871391589*Ghat[15]*dv12*dx12); 
  out[41] += 0.7071067811865475*Ghat[25]*dv12*dx12; 
  out[42] += -(1.224744871391589*Ghat[16]*dv12*dx12); 
  out[43] += -(1.224744871391589*Ghat[17]*dv12*dx12); 
  out[44] += 0.7071067811865475*Ghat[26]*dv12*dx12; 
  out[45] += -(1.224744871391589*Ghat[18]*dv12*dx12); 
  out[46] += -(1.224744871391589*Ghat[19]*dv12*dx12); 
  out[47] += -(1.224744871391589*Ghat[20]*dv12*dx12); 
  out[48] += 0.7071067811865475*Ghat[27]*dv12*dx12; 
  out[49] += -(1.224744871391589*Ghat[21]*dv12*dx12); 
  out[50] += -(1.224744871391589*Ghat[22]*dv12*dx12); 
  out[51] += 0.7071067811865475*Ghat[28]*dv12*dx12; 
  out[52] += -(1.224744871391589*Ghat[23]*dv12*dx12); 
  out[53] += -(1.224744871391589*Ghat[24]*dv12*dx12); 
  out[54] += 0.7071067811865475*Ghat[29]*dv12*dx12; 
  out[55] += 0.7071067811865475*Ghat[30]*dv12*dx12; 
  out[56] += -(1.224744871391589*Ghat[25]*dv12*dx12); 
  out[57] += -(1.224744871391589*Ghat[26]*dv12*dx12); 
  out[58] += -(1.224744871391589*Ghat[27]*dv12*dx12); 
  out[59] += -(1.224744871391589*Ghat[28]*dv12*dx12); 
  out[60] += 0.7071067811865475*Ghat[31]*dv12*dx12; 
  out[61] += -(1.224744871391589*Ghat[29]*dv12*dx12); 
  out[62] += -(1.224744871391589*Ghat[30]*dv12*dx12); 
  out[63] += -(1.224744871391589*Ghat[31]*dv12*dx12); 

  } 
  return 0.0;

} 
