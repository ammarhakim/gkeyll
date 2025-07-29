#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx11 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[4]; 
  const double wv = w[4]; 

  double Ghat[32]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[36]+0.4330127018922193*hamil[7]*fskin[20]+0.75*hamil[6]*fskin[18]+0.75*hamil[4]*fskin[11]+0.4330127018922193*fskin[6]*hamil[6]+0.4330127018922193*fskin[4]*hamil[4]+0.75*fskin[2]*hamil[2]+0.4330127018922193*fskin[0]*hamil[2]; 
  Ghat[1] = 0.75*hamil[7]*fskin[48]+0.4330127018922193*hamil[7]*fskin[35]+0.75*hamil[6]*fskin[32]+0.75*hamil[4]*fskin[23]+0.4330127018922193*hamil[6]*fskin[17]+0.4330127018922193*hamil[4]*fskin[10]+0.75*hamil[2]*fskin[7]+0.4330127018922193*fskin[1]*hamil[2]; 
  Ghat[2] = 0.75*hamil[7]*fskin[50]+0.4330127018922193*hamil[7]*fskin[37]+0.75*hamil[6]*fskin[34]+0.75*hamil[4]*fskin[25]+0.4330127018922193*hamil[6]*fskin[19]+0.4330127018922193*hamil[4]*fskin[12]+0.75*hamil[2]*fskin[9]+0.4330127018922193*hamil[2]*fskin[3]; 
  Ghat[3] = 0.75*hamil[6]*fskin[36]+0.4330127018922193*hamil[6]*fskin[20]+0.75*hamil[7]*fskin[18]+0.75*hamil[2]*fskin[11]+0.4330127018922193*fskin[6]*hamil[7]+0.75*fskin[2]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[2]*fskin[4]; 
  Ghat[4] = 0.75*hamil[7]*fskin[55]+0.4330127018922193*hamil[7]*fskin[41]+0.75*hamil[6]*fskin[39]+0.75*hamil[4]*fskin[30]+0.4330127018922193*hamil[6]*fskin[21]+0.4330127018922193*hamil[4]*fskin[16]+0.75*hamil[2]*fskin[14]+0.4330127018922193*hamil[2]*fskin[5]; 
  Ghat[5] = 0.75*hamil[4]*fskin[36]+0.4330127018922193*hamil[4]*fskin[20]+0.75*hamil[2]*fskin[18]+0.75*hamil[7]*fskin[11]+0.4330127018922193*fskin[4]*hamil[7]+0.75*fskin[2]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[2]*fskin[6]; 
  Ghat[6] = 0.75*hamil[7]*fskin[58]+0.4330127018922193*hamil[7]*fskin[49]+0.75*hamil[6]*fskin[47]+0.75*hamil[4]*fskin[42]+0.4330127018922193*hamil[6]*fskin[33]+0.4330127018922193*hamil[4]*fskin[24]+0.75*hamil[2]*fskin[22]+0.4330127018922193*hamil[2]*fskin[8]; 
  Ghat[7] = 0.75*hamil[6]*fskin[48]+0.4330127018922193*hamil[6]*fskin[35]+0.75*hamil[7]*fskin[32]+0.75*hamil[2]*fskin[23]+0.4330127018922193*hamil[7]*fskin[17]+0.4330127018922193*hamil[2]*fskin[10]+0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[1]*hamil[4]; 
  Ghat[8] = 0.75*hamil[6]*fskin[50]+0.4330127018922193*hamil[6]*fskin[37]+0.75*hamil[7]*fskin[34]+0.75*hamil[2]*fskin[25]+0.4330127018922193*hamil[7]*fskin[19]+0.4330127018922193*hamil[2]*fskin[12]+0.75*hamil[4]*fskin[9]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[9] = 0.75*hamil[7]*fskin[60]+0.4330127018922193*hamil[7]*fskin[54]+0.75*hamil[6]*fskin[51]+0.75*hamil[4]*fskin[44]+0.4330127018922193*hamil[6]*fskin[38]+0.4330127018922193*hamil[4]*fskin[29]+0.75*hamil[2]*fskin[26]+0.4330127018922193*hamil[2]*fskin[13]; 
  Ghat[10] = 0.75*hamil[7]*fskin[62]+0.4330127018922193*hamil[7]*fskin[56]+0.75*hamil[6]*fskin[53]+0.75*hamil[4]*fskin[46]+0.4330127018922193*hamil[6]*fskin[40]+0.4330127018922193*hamil[4]*fskin[31]+0.75*hamil[2]*fskin[28]+0.4330127018922193*hamil[2]*fskin[15]; 
  Ghat[11] = 0.75*hamil[6]*fskin[55]+0.4330127018922193*hamil[6]*fskin[41]+0.75*hamil[7]*fskin[39]+0.75*hamil[2]*fskin[30]+0.4330127018922193*hamil[7]*fskin[21]+0.4330127018922193*hamil[2]*fskin[16]+0.75*hamil[4]*fskin[14]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[12] = 0.75*hamil[4]*fskin[48]+0.4330127018922193*hamil[4]*fskin[35]+0.75*hamil[2]*fskin[32]+0.75*hamil[7]*fskin[23]+0.4330127018922193*hamil[2]*fskin[17]+0.4330127018922193*hamil[7]*fskin[10]+0.75*hamil[6]*fskin[7]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[13] = 0.75*hamil[4]*fskin[50]+0.4330127018922193*hamil[4]*fskin[37]+0.75*hamil[2]*fskin[34]+0.75*hamil[7]*fskin[25]+0.4330127018922193*hamil[2]*fskin[19]+0.4330127018922193*hamil[7]*fskin[12]+0.75*hamil[6]*fskin[9]+0.4330127018922193*fskin[3]*hamil[6]; 
  Ghat[14] = 0.75*hamil[2]*fskin[36]+0.4330127018922193*hamil[2]*fskin[20]+0.75*hamil[4]*fskin[18]+0.75*hamil[6]*fskin[11]+0.75*fskin[2]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[6]+0.4330127018922193*hamil[4]*fskin[6]; 
  Ghat[15] = 0.75*hamil[4]*fskin[55]+0.4330127018922193*hamil[4]*fskin[41]+0.75*hamil[2]*fskin[39]+0.75*hamil[7]*fskin[30]+0.4330127018922193*hamil[2]*fskin[21]+0.4330127018922193*hamil[7]*fskin[16]+0.75*hamil[6]*fskin[14]+0.4330127018922193*fskin[5]*hamil[6]; 
  Ghat[16] = 0.75*hamil[6]*fskin[58]+0.4330127018922193*hamil[6]*fskin[49]+0.75*hamil[7]*fskin[47]+0.75*hamil[2]*fskin[42]+0.4330127018922193*hamil[7]*fskin[33]+0.4330127018922193*hamil[2]*fskin[24]+0.75*hamil[4]*fskin[22]+0.4330127018922193*hamil[4]*fskin[8]; 
  Ghat[17] = 0.75*hamil[7]*fskin[63]+0.4330127018922193*hamil[7]*fskin[61]+0.75*hamil[6]*fskin[59]+0.75*hamil[4]*fskin[57]+0.4330127018922193*hamil[6]*fskin[52]+0.4330127018922193*hamil[4]*fskin[45]+0.75*hamil[2]*fskin[43]+0.4330127018922193*hamil[2]*fskin[27]; 
  Ghat[18] = 0.75*hamil[6]*fskin[60]+0.4330127018922193*hamil[6]*fskin[54]+0.75*hamil[7]*fskin[51]+0.75*hamil[2]*fskin[44]+0.4330127018922193*hamil[7]*fskin[38]+0.4330127018922193*hamil[2]*fskin[29]+0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[4]*fskin[13]; 
  Ghat[19] = 0.75*hamil[6]*fskin[62]+0.4330127018922193*hamil[6]*fskin[56]+0.75*hamil[7]*fskin[53]+0.75*hamil[2]*fskin[46]+0.4330127018922193*hamil[7]*fskin[40]+0.4330127018922193*hamil[2]*fskin[31]+0.75*hamil[4]*fskin[28]+0.4330127018922193*hamil[4]*fskin[15]; 
  Ghat[20] = 0.75*hamil[4]*fskin[58]+0.4330127018922193*hamil[4]*fskin[49]+0.75*hamil[2]*fskin[47]+0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[2]*fskin[33]+0.4330127018922193*hamil[7]*fskin[24]+0.75*hamil[6]*fskin[22]+0.4330127018922193*hamil[6]*fskin[8]; 
  Ghat[21] = 0.75*hamil[2]*fskin[48]+0.4330127018922193*hamil[2]*fskin[35]+0.75*hamil[4]*fskin[32]+0.75*hamil[6]*fskin[23]+0.4330127018922193*hamil[4]*fskin[17]+0.4330127018922193*hamil[6]*fskin[10]+0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[1]*hamil[7]; 
  Ghat[22] = 0.75*hamil[2]*fskin[50]+0.4330127018922193*hamil[2]*fskin[37]+0.75*hamil[4]*fskin[34]+0.75*hamil[6]*fskin[25]+0.4330127018922193*hamil[4]*fskin[19]+0.4330127018922193*hamil[6]*fskin[12]+0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[23] = 0.75*hamil[4]*fskin[60]+0.4330127018922193*hamil[4]*fskin[54]+0.75*hamil[2]*fskin[51]+0.75*hamil[7]*fskin[44]+0.4330127018922193*hamil[2]*fskin[38]+0.4330127018922193*hamil[7]*fskin[29]+0.75*hamil[6]*fskin[26]+0.4330127018922193*hamil[6]*fskin[13]; 
  Ghat[24] = 0.75*hamil[4]*fskin[62]+0.4330127018922193*hamil[4]*fskin[56]+0.75*hamil[2]*fskin[53]+0.75*hamil[7]*fskin[46]+0.4330127018922193*hamil[2]*fskin[40]+0.4330127018922193*hamil[7]*fskin[31]+0.75*hamil[6]*fskin[28]+0.4330127018922193*hamil[6]*fskin[15]; 
  Ghat[25] = 0.75*hamil[2]*fskin[55]+0.4330127018922193*hamil[2]*fskin[41]+0.75*hamil[4]*fskin[39]+0.75*hamil[6]*fskin[30]+0.4330127018922193*hamil[4]*fskin[21]+0.4330127018922193*hamil[6]*fskin[16]+0.75*hamil[7]*fskin[14]+0.4330127018922193*fskin[5]*hamil[7]; 
  Ghat[26] = 0.75*hamil[6]*fskin[63]+0.4330127018922193*hamil[6]*fskin[61]+0.75*hamil[7]*fskin[59]+0.75*hamil[2]*fskin[57]+0.4330127018922193*hamil[7]*fskin[52]+0.4330127018922193*hamil[2]*fskin[45]+0.75*hamil[4]*fskin[43]+0.4330127018922193*hamil[4]*fskin[27]; 
  Ghat[27] = 0.75*hamil[2]*fskin[58]+0.4330127018922193*hamil[2]*fskin[49]+0.75*hamil[4]*fskin[47]+0.75*hamil[6]*fskin[42]+0.4330127018922193*hamil[4]*fskin[33]+0.4330127018922193*hamil[6]*fskin[24]+0.75*hamil[7]*fskin[22]+0.4330127018922193*hamil[7]*fskin[8]; 
  Ghat[28] = 0.75*hamil[4]*fskin[63]+0.4330127018922193*hamil[4]*fskin[61]+0.75*hamil[2]*fskin[59]+0.75*hamil[7]*fskin[57]+0.4330127018922193*hamil[2]*fskin[52]+0.4330127018922193*hamil[7]*fskin[45]+0.75*hamil[6]*fskin[43]+0.4330127018922193*hamil[6]*fskin[27]; 
  Ghat[29] = 0.75*hamil[2]*fskin[60]+0.4330127018922193*hamil[2]*fskin[54]+0.75*hamil[4]*fskin[51]+0.75*hamil[6]*fskin[44]+0.4330127018922193*hamil[4]*fskin[38]+0.4330127018922193*hamil[6]*fskin[29]+0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[7]*fskin[13]; 
  Ghat[30] = 0.75*hamil[2]*fskin[62]+0.4330127018922193*hamil[2]*fskin[56]+0.75*hamil[4]*fskin[53]+0.75*hamil[6]*fskin[46]+0.4330127018922193*hamil[4]*fskin[40]+0.4330127018922193*hamil[6]*fskin[31]+0.75*hamil[7]*fskin[28]+0.4330127018922193*hamil[7]*fskin[15]; 
  Ghat[31] = 0.75*hamil[2]*fskin[63]+0.4330127018922193*hamil[2]*fskin[61]+0.75*hamil[4]*fskin[59]+0.75*hamil[6]*fskin[57]+0.4330127018922193*hamil[4]*fskin[52]+0.4330127018922193*hamil[6]*fskin[45]+0.75*hamil[7]*fskin[43]+0.4330127018922193*hamil[7]*fskin[27]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[36])+0.4330127018922193*hamil[7]*fedge[20]-0.75*hamil[6]*fedge[18]-0.75*hamil[4]*fedge[11]+0.4330127018922193*fedge[6]*hamil[6]+0.4330127018922193*fedge[4]*hamil[4]-0.75*fedge[2]*hamil[2]+0.4330127018922193*fedge[0]*hamil[2]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[48])+0.4330127018922193*hamil[7]*fedge[35]-0.75*hamil[6]*fedge[32]-0.75*hamil[4]*fedge[23]+0.4330127018922193*hamil[6]*fedge[17]+0.4330127018922193*hamil[4]*fedge[10]-0.75*hamil[2]*fedge[7]+0.4330127018922193*fedge[1]*hamil[2]; 
  Ghat[2] = -(0.75*hamil[7]*fedge[50])+0.4330127018922193*hamil[7]*fedge[37]-0.75*hamil[6]*fedge[34]-0.75*hamil[4]*fedge[25]+0.4330127018922193*hamil[6]*fedge[19]+0.4330127018922193*hamil[4]*fedge[12]-0.75*hamil[2]*fedge[9]+0.4330127018922193*hamil[2]*fedge[3]; 
  Ghat[3] = -(0.75*hamil[6]*fedge[36])+0.4330127018922193*hamil[6]*fedge[20]-0.75*hamil[7]*fedge[18]-0.75*hamil[2]*fedge[11]+0.4330127018922193*fedge[6]*hamil[7]-0.75*fedge[2]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[2]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[7]*fedge[55])+0.4330127018922193*hamil[7]*fedge[41]-0.75*hamil[6]*fedge[39]-0.75*hamil[4]*fedge[30]+0.4330127018922193*hamil[6]*fedge[21]+0.4330127018922193*hamil[4]*fedge[16]-0.75*hamil[2]*fedge[14]+0.4330127018922193*hamil[2]*fedge[5]; 
  Ghat[5] = -(0.75*hamil[4]*fedge[36])+0.4330127018922193*hamil[4]*fedge[20]-0.75*hamil[2]*fedge[18]-0.75*hamil[7]*fedge[11]+0.4330127018922193*fedge[4]*hamil[7]-0.75*fedge[2]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[2]*fedge[6]; 
  Ghat[6] = -(0.75*hamil[7]*fedge[58])+0.4330127018922193*hamil[7]*fedge[49]-0.75*hamil[6]*fedge[47]-0.75*hamil[4]*fedge[42]+0.4330127018922193*hamil[6]*fedge[33]+0.4330127018922193*hamil[4]*fedge[24]-0.75*hamil[2]*fedge[22]+0.4330127018922193*hamil[2]*fedge[8]; 
  Ghat[7] = -(0.75*hamil[6]*fedge[48])+0.4330127018922193*hamil[6]*fedge[35]-0.75*hamil[7]*fedge[32]-0.75*hamil[2]*fedge[23]+0.4330127018922193*hamil[7]*fedge[17]+0.4330127018922193*hamil[2]*fedge[10]-0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[1]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[6]*fedge[50])+0.4330127018922193*hamil[6]*fedge[37]-0.75*hamil[7]*fedge[34]-0.75*hamil[2]*fedge[25]+0.4330127018922193*hamil[7]*fedge[19]+0.4330127018922193*hamil[2]*fedge[12]-0.75*hamil[4]*fedge[9]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[9] = -(0.75*hamil[7]*fedge[60])+0.4330127018922193*hamil[7]*fedge[54]-0.75*hamil[6]*fedge[51]-0.75*hamil[4]*fedge[44]+0.4330127018922193*hamil[6]*fedge[38]+0.4330127018922193*hamil[4]*fedge[29]-0.75*hamil[2]*fedge[26]+0.4330127018922193*hamil[2]*fedge[13]; 
  Ghat[10] = -(0.75*hamil[7]*fedge[62])+0.4330127018922193*hamil[7]*fedge[56]-0.75*hamil[6]*fedge[53]-0.75*hamil[4]*fedge[46]+0.4330127018922193*hamil[6]*fedge[40]+0.4330127018922193*hamil[4]*fedge[31]-0.75*hamil[2]*fedge[28]+0.4330127018922193*hamil[2]*fedge[15]; 
  Ghat[11] = -(0.75*hamil[6]*fedge[55])+0.4330127018922193*hamil[6]*fedge[41]-0.75*hamil[7]*fedge[39]-0.75*hamil[2]*fedge[30]+0.4330127018922193*hamil[7]*fedge[21]+0.4330127018922193*hamil[2]*fedge[16]-0.75*hamil[4]*fedge[14]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[12] = -(0.75*hamil[4]*fedge[48])+0.4330127018922193*hamil[4]*fedge[35]-0.75*hamil[2]*fedge[32]-0.75*hamil[7]*fedge[23]+0.4330127018922193*hamil[2]*fedge[17]+0.4330127018922193*hamil[7]*fedge[10]-0.75*hamil[6]*fedge[7]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[13] = -(0.75*hamil[4]*fedge[50])+0.4330127018922193*hamil[4]*fedge[37]-0.75*hamil[2]*fedge[34]-0.75*hamil[7]*fedge[25]+0.4330127018922193*hamil[2]*fedge[19]+0.4330127018922193*hamil[7]*fedge[12]-0.75*hamil[6]*fedge[9]+0.4330127018922193*fedge[3]*hamil[6]; 
  Ghat[14] = -(0.75*hamil[2]*fedge[36])+0.4330127018922193*hamil[2]*fedge[20]-0.75*hamil[4]*fedge[18]-0.75*hamil[6]*fedge[11]-0.75*fedge[2]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[6]+0.4330127018922193*hamil[4]*fedge[6]; 
  Ghat[15] = -(0.75*hamil[4]*fedge[55])+0.4330127018922193*hamil[4]*fedge[41]-0.75*hamil[2]*fedge[39]-0.75*hamil[7]*fedge[30]+0.4330127018922193*hamil[2]*fedge[21]+0.4330127018922193*hamil[7]*fedge[16]-0.75*hamil[6]*fedge[14]+0.4330127018922193*fedge[5]*hamil[6]; 
  Ghat[16] = -(0.75*hamil[6]*fedge[58])+0.4330127018922193*hamil[6]*fedge[49]-0.75*hamil[7]*fedge[47]-0.75*hamil[2]*fedge[42]+0.4330127018922193*hamil[7]*fedge[33]+0.4330127018922193*hamil[2]*fedge[24]-0.75*hamil[4]*fedge[22]+0.4330127018922193*hamil[4]*fedge[8]; 
  Ghat[17] = -(0.75*hamil[7]*fedge[63])+0.4330127018922193*hamil[7]*fedge[61]-0.75*hamil[6]*fedge[59]-0.75*hamil[4]*fedge[57]+0.4330127018922193*hamil[6]*fedge[52]+0.4330127018922193*hamil[4]*fedge[45]-0.75*hamil[2]*fedge[43]+0.4330127018922193*hamil[2]*fedge[27]; 
  Ghat[18] = -(0.75*hamil[6]*fedge[60])+0.4330127018922193*hamil[6]*fedge[54]-0.75*hamil[7]*fedge[51]-0.75*hamil[2]*fedge[44]+0.4330127018922193*hamil[7]*fedge[38]+0.4330127018922193*hamil[2]*fedge[29]-0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[4]*fedge[13]; 
  Ghat[19] = -(0.75*hamil[6]*fedge[62])+0.4330127018922193*hamil[6]*fedge[56]-0.75*hamil[7]*fedge[53]-0.75*hamil[2]*fedge[46]+0.4330127018922193*hamil[7]*fedge[40]+0.4330127018922193*hamil[2]*fedge[31]-0.75*hamil[4]*fedge[28]+0.4330127018922193*hamil[4]*fedge[15]; 
  Ghat[20] = -(0.75*hamil[4]*fedge[58])+0.4330127018922193*hamil[4]*fedge[49]-0.75*hamil[2]*fedge[47]-0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[2]*fedge[33]+0.4330127018922193*hamil[7]*fedge[24]-0.75*hamil[6]*fedge[22]+0.4330127018922193*hamil[6]*fedge[8]; 
  Ghat[21] = -(0.75*hamil[2]*fedge[48])+0.4330127018922193*hamil[2]*fedge[35]-0.75*hamil[4]*fedge[32]-0.75*hamil[6]*fedge[23]+0.4330127018922193*hamil[4]*fedge[17]+0.4330127018922193*hamil[6]*fedge[10]-0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[1]*hamil[7]; 
  Ghat[22] = -(0.75*hamil[2]*fedge[50])+0.4330127018922193*hamil[2]*fedge[37]-0.75*hamil[4]*fedge[34]-0.75*hamil[6]*fedge[25]+0.4330127018922193*hamil[4]*fedge[19]+0.4330127018922193*hamil[6]*fedge[12]-0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[23] = -(0.75*hamil[4]*fedge[60])+0.4330127018922193*hamil[4]*fedge[54]-0.75*hamil[2]*fedge[51]-0.75*hamil[7]*fedge[44]+0.4330127018922193*hamil[2]*fedge[38]+0.4330127018922193*hamil[7]*fedge[29]-0.75*hamil[6]*fedge[26]+0.4330127018922193*hamil[6]*fedge[13]; 
  Ghat[24] = -(0.75*hamil[4]*fedge[62])+0.4330127018922193*hamil[4]*fedge[56]-0.75*hamil[2]*fedge[53]-0.75*hamil[7]*fedge[46]+0.4330127018922193*hamil[2]*fedge[40]+0.4330127018922193*hamil[7]*fedge[31]-0.75*hamil[6]*fedge[28]+0.4330127018922193*hamil[6]*fedge[15]; 
  Ghat[25] = -(0.75*hamil[2]*fedge[55])+0.4330127018922193*hamil[2]*fedge[41]-0.75*hamil[4]*fedge[39]-0.75*hamil[6]*fedge[30]+0.4330127018922193*hamil[4]*fedge[21]+0.4330127018922193*hamil[6]*fedge[16]-0.75*hamil[7]*fedge[14]+0.4330127018922193*fedge[5]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[6]*fedge[63])+0.4330127018922193*hamil[6]*fedge[61]-0.75*hamil[7]*fedge[59]-0.75*hamil[2]*fedge[57]+0.4330127018922193*hamil[7]*fedge[52]+0.4330127018922193*hamil[2]*fedge[45]-0.75*hamil[4]*fedge[43]+0.4330127018922193*hamil[4]*fedge[27]; 
  Ghat[27] = -(0.75*hamil[2]*fedge[58])+0.4330127018922193*hamil[2]*fedge[49]-0.75*hamil[4]*fedge[47]-0.75*hamil[6]*fedge[42]+0.4330127018922193*hamil[4]*fedge[33]+0.4330127018922193*hamil[6]*fedge[24]-0.75*hamil[7]*fedge[22]+0.4330127018922193*hamil[7]*fedge[8]; 
  Ghat[28] = -(0.75*hamil[4]*fedge[63])+0.4330127018922193*hamil[4]*fedge[61]-0.75*hamil[2]*fedge[59]-0.75*hamil[7]*fedge[57]+0.4330127018922193*hamil[2]*fedge[52]+0.4330127018922193*hamil[7]*fedge[45]-0.75*hamil[6]*fedge[43]+0.4330127018922193*hamil[6]*fedge[27]; 
  Ghat[29] = -(0.75*hamil[2]*fedge[60])+0.4330127018922193*hamil[2]*fedge[54]-0.75*hamil[4]*fedge[51]-0.75*hamil[6]*fedge[44]+0.4330127018922193*hamil[4]*fedge[38]+0.4330127018922193*hamil[6]*fedge[29]-0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[7]*fedge[13]; 
  Ghat[30] = -(0.75*hamil[2]*fedge[62])+0.4330127018922193*hamil[2]*fedge[56]-0.75*hamil[4]*fedge[53]-0.75*hamil[6]*fedge[46]+0.4330127018922193*hamil[4]*fedge[40]+0.4330127018922193*hamil[6]*fedge[31]-0.75*hamil[7]*fedge[28]+0.4330127018922193*hamil[7]*fedge[15]; 
  Ghat[31] = -(0.75*hamil[2]*fedge[63])+0.4330127018922193*hamil[2]*fedge[61]-0.75*hamil[4]*fedge[59]-0.75*hamil[6]*fedge[57]+0.4330127018922193*hamil[4]*fedge[52]+0.4330127018922193*hamil[6]*fedge[45]-0.75*hamil[7]*fedge[43]+0.4330127018922193*hamil[7]*fedge[27]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11*dx11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11*dx11); 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv11*dx11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11*dx11); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv11*dx11); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv11*dx11); 
  out[7] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[8] += -(0.7071067811865475*Ghat[6]*dv11*dx11); 
  out[9] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[10] += -(0.7071067811865475*Ghat[7]*dv11*dx11); 
  out[11] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv11*dx11); 
  out[13] += -(0.7071067811865475*Ghat[9]*dv11*dx11); 
  out[14] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv11*dx11); 
  out[16] += -(0.7071067811865475*Ghat[11]*dv11*dx11); 
  out[17] += -(0.7071067811865475*Ghat[12]*dv11*dx11); 
  out[18] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[19] += -(0.7071067811865475*Ghat[13]*dv11*dx11); 
  out[20] += -(0.7071067811865475*Ghat[14]*dv11*dx11); 
  out[21] += -(0.7071067811865475*Ghat[15]*dv11*dx11); 
  out[22] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[23] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[24] += -(0.7071067811865475*Ghat[16]*dv11*dx11); 
  out[25] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[26] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[27] += -(0.7071067811865475*Ghat[17]*dv11*dx11); 
  out[28] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[29] += -(0.7071067811865475*Ghat[18]*dv11*dx11); 
  out[30] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[31] += -(0.7071067811865475*Ghat[19]*dv11*dx11); 
  out[32] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[33] += -(0.7071067811865475*Ghat[20]*dv11*dx11); 
  out[34] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[35] += -(0.7071067811865475*Ghat[21]*dv11*dx11); 
  out[36] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[37] += -(0.7071067811865475*Ghat[22]*dv11*dx11); 
  out[38] += -(0.7071067811865475*Ghat[23]*dv11*dx11); 
  out[39] += -(1.224744871391589*Ghat[15]*dv11*dx11); 
  out[40] += -(0.7071067811865475*Ghat[24]*dv11*dx11); 
  out[41] += -(0.7071067811865475*Ghat[25]*dv11*dx11); 
  out[42] += -(1.224744871391589*Ghat[16]*dv11*dx11); 
  out[43] += -(1.224744871391589*Ghat[17]*dv11*dx11); 
  out[44] += -(1.224744871391589*Ghat[18]*dv11*dx11); 
  out[45] += -(0.7071067811865475*Ghat[26]*dv11*dx11); 
  out[46] += -(1.224744871391589*Ghat[19]*dv11*dx11); 
  out[47] += -(1.224744871391589*Ghat[20]*dv11*dx11); 
  out[48] += -(1.224744871391589*Ghat[21]*dv11*dx11); 
  out[49] += -(0.7071067811865475*Ghat[27]*dv11*dx11); 
  out[50] += -(1.224744871391589*Ghat[22]*dv11*dx11); 
  out[51] += -(1.224744871391589*Ghat[23]*dv11*dx11); 
  out[52] += -(0.7071067811865475*Ghat[28]*dv11*dx11); 
  out[53] += -(1.224744871391589*Ghat[24]*dv11*dx11); 
  out[54] += -(0.7071067811865475*Ghat[29]*dv11*dx11); 
  out[55] += -(1.224744871391589*Ghat[25]*dv11*dx11); 
  out[56] += -(0.7071067811865475*Ghat[30]*dv11*dx11); 
  out[57] += -(1.224744871391589*Ghat[26]*dv11*dx11); 
  out[58] += -(1.224744871391589*Ghat[27]*dv11*dx11); 
  out[59] += -(1.224744871391589*Ghat[28]*dv11*dx11); 
  out[60] += -(1.224744871391589*Ghat[29]*dv11*dx11); 
  out[61] += -(0.7071067811865475*Ghat[31]*dv11*dx11); 
  out[62] += -(1.224744871391589*Ghat[30]*dv11*dx11); 
  out[63] += -(1.224744871391589*Ghat[31]*dv11*dx11); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[36]+0.4330127018922193*hamil[7]*fedge[20]+0.75*hamil[6]*fedge[18]+0.75*hamil[4]*fedge[11]+0.4330127018922193*fedge[6]*hamil[6]+0.4330127018922193*fedge[4]*hamil[4]+0.75*fedge[2]*hamil[2]+0.4330127018922193*fedge[0]*hamil[2]; 
  Ghat[1] = 0.75*hamil[7]*fedge[48]+0.4330127018922193*hamil[7]*fedge[35]+0.75*hamil[6]*fedge[32]+0.75*hamil[4]*fedge[23]+0.4330127018922193*hamil[6]*fedge[17]+0.4330127018922193*hamil[4]*fedge[10]+0.75*hamil[2]*fedge[7]+0.4330127018922193*fedge[1]*hamil[2]; 
  Ghat[2] = 0.75*hamil[7]*fedge[50]+0.4330127018922193*hamil[7]*fedge[37]+0.75*hamil[6]*fedge[34]+0.75*hamil[4]*fedge[25]+0.4330127018922193*hamil[6]*fedge[19]+0.4330127018922193*hamil[4]*fedge[12]+0.75*hamil[2]*fedge[9]+0.4330127018922193*hamil[2]*fedge[3]; 
  Ghat[3] = 0.75*hamil[6]*fedge[36]+0.4330127018922193*hamil[6]*fedge[20]+0.75*hamil[7]*fedge[18]+0.75*hamil[2]*fedge[11]+0.4330127018922193*fedge[6]*hamil[7]+0.75*fedge[2]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[2]*fedge[4]; 
  Ghat[4] = 0.75*hamil[7]*fedge[55]+0.4330127018922193*hamil[7]*fedge[41]+0.75*hamil[6]*fedge[39]+0.75*hamil[4]*fedge[30]+0.4330127018922193*hamil[6]*fedge[21]+0.4330127018922193*hamil[4]*fedge[16]+0.75*hamil[2]*fedge[14]+0.4330127018922193*hamil[2]*fedge[5]; 
  Ghat[5] = 0.75*hamil[4]*fedge[36]+0.4330127018922193*hamil[4]*fedge[20]+0.75*hamil[2]*fedge[18]+0.75*hamil[7]*fedge[11]+0.4330127018922193*fedge[4]*hamil[7]+0.75*fedge[2]*hamil[6]+0.4330127018922193*fedge[0]*hamil[6]+0.4330127018922193*hamil[2]*fedge[6]; 
  Ghat[6] = 0.75*hamil[7]*fedge[58]+0.4330127018922193*hamil[7]*fedge[49]+0.75*hamil[6]*fedge[47]+0.75*hamil[4]*fedge[42]+0.4330127018922193*hamil[6]*fedge[33]+0.4330127018922193*hamil[4]*fedge[24]+0.75*hamil[2]*fedge[22]+0.4330127018922193*hamil[2]*fedge[8]; 
  Ghat[7] = 0.75*hamil[6]*fedge[48]+0.4330127018922193*hamil[6]*fedge[35]+0.75*hamil[7]*fedge[32]+0.75*hamil[2]*fedge[23]+0.4330127018922193*hamil[7]*fedge[17]+0.4330127018922193*hamil[2]*fedge[10]+0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[1]*hamil[4]; 
  Ghat[8] = 0.75*hamil[6]*fedge[50]+0.4330127018922193*hamil[6]*fedge[37]+0.75*hamil[7]*fedge[34]+0.75*hamil[2]*fedge[25]+0.4330127018922193*hamil[7]*fedge[19]+0.4330127018922193*hamil[2]*fedge[12]+0.75*hamil[4]*fedge[9]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[9] = 0.75*hamil[7]*fedge[60]+0.4330127018922193*hamil[7]*fedge[54]+0.75*hamil[6]*fedge[51]+0.75*hamil[4]*fedge[44]+0.4330127018922193*hamil[6]*fedge[38]+0.4330127018922193*hamil[4]*fedge[29]+0.75*hamil[2]*fedge[26]+0.4330127018922193*hamil[2]*fedge[13]; 
  Ghat[10] = 0.75*hamil[7]*fedge[62]+0.4330127018922193*hamil[7]*fedge[56]+0.75*hamil[6]*fedge[53]+0.75*hamil[4]*fedge[46]+0.4330127018922193*hamil[6]*fedge[40]+0.4330127018922193*hamil[4]*fedge[31]+0.75*hamil[2]*fedge[28]+0.4330127018922193*hamil[2]*fedge[15]; 
  Ghat[11] = 0.75*hamil[6]*fedge[55]+0.4330127018922193*hamil[6]*fedge[41]+0.75*hamil[7]*fedge[39]+0.75*hamil[2]*fedge[30]+0.4330127018922193*hamil[7]*fedge[21]+0.4330127018922193*hamil[2]*fedge[16]+0.75*hamil[4]*fedge[14]+0.4330127018922193*hamil[4]*fedge[5]; 
  Ghat[12] = 0.75*hamil[4]*fedge[48]+0.4330127018922193*hamil[4]*fedge[35]+0.75*hamil[2]*fedge[32]+0.75*hamil[7]*fedge[23]+0.4330127018922193*hamil[2]*fedge[17]+0.4330127018922193*hamil[7]*fedge[10]+0.75*hamil[6]*fedge[7]+0.4330127018922193*fedge[1]*hamil[6]; 
  Ghat[13] = 0.75*hamil[4]*fedge[50]+0.4330127018922193*hamil[4]*fedge[37]+0.75*hamil[2]*fedge[34]+0.75*hamil[7]*fedge[25]+0.4330127018922193*hamil[2]*fedge[19]+0.4330127018922193*hamil[7]*fedge[12]+0.75*hamil[6]*fedge[9]+0.4330127018922193*fedge[3]*hamil[6]; 
  Ghat[14] = 0.75*hamil[2]*fedge[36]+0.4330127018922193*hamil[2]*fedge[20]+0.75*hamil[4]*fedge[18]+0.75*hamil[6]*fedge[11]+0.75*fedge[2]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*fedge[4]*hamil[6]+0.4330127018922193*hamil[4]*fedge[6]; 
  Ghat[15] = 0.75*hamil[4]*fedge[55]+0.4330127018922193*hamil[4]*fedge[41]+0.75*hamil[2]*fedge[39]+0.75*hamil[7]*fedge[30]+0.4330127018922193*hamil[2]*fedge[21]+0.4330127018922193*hamil[7]*fedge[16]+0.75*hamil[6]*fedge[14]+0.4330127018922193*fedge[5]*hamil[6]; 
  Ghat[16] = 0.75*hamil[6]*fedge[58]+0.4330127018922193*hamil[6]*fedge[49]+0.75*hamil[7]*fedge[47]+0.75*hamil[2]*fedge[42]+0.4330127018922193*hamil[7]*fedge[33]+0.4330127018922193*hamil[2]*fedge[24]+0.75*hamil[4]*fedge[22]+0.4330127018922193*hamil[4]*fedge[8]; 
  Ghat[17] = 0.75*hamil[7]*fedge[63]+0.4330127018922193*hamil[7]*fedge[61]+0.75*hamil[6]*fedge[59]+0.75*hamil[4]*fedge[57]+0.4330127018922193*hamil[6]*fedge[52]+0.4330127018922193*hamil[4]*fedge[45]+0.75*hamil[2]*fedge[43]+0.4330127018922193*hamil[2]*fedge[27]; 
  Ghat[18] = 0.75*hamil[6]*fedge[60]+0.4330127018922193*hamil[6]*fedge[54]+0.75*hamil[7]*fedge[51]+0.75*hamil[2]*fedge[44]+0.4330127018922193*hamil[7]*fedge[38]+0.4330127018922193*hamil[2]*fedge[29]+0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[4]*fedge[13]; 
  Ghat[19] = 0.75*hamil[6]*fedge[62]+0.4330127018922193*hamil[6]*fedge[56]+0.75*hamil[7]*fedge[53]+0.75*hamil[2]*fedge[46]+0.4330127018922193*hamil[7]*fedge[40]+0.4330127018922193*hamil[2]*fedge[31]+0.75*hamil[4]*fedge[28]+0.4330127018922193*hamil[4]*fedge[15]; 
  Ghat[20] = 0.75*hamil[4]*fedge[58]+0.4330127018922193*hamil[4]*fedge[49]+0.75*hamil[2]*fedge[47]+0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[2]*fedge[33]+0.4330127018922193*hamil[7]*fedge[24]+0.75*hamil[6]*fedge[22]+0.4330127018922193*hamil[6]*fedge[8]; 
  Ghat[21] = 0.75*hamil[2]*fedge[48]+0.4330127018922193*hamil[2]*fedge[35]+0.75*hamil[4]*fedge[32]+0.75*hamil[6]*fedge[23]+0.4330127018922193*hamil[4]*fedge[17]+0.4330127018922193*hamil[6]*fedge[10]+0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[1]*hamil[7]; 
  Ghat[22] = 0.75*hamil[2]*fedge[50]+0.4330127018922193*hamil[2]*fedge[37]+0.75*hamil[4]*fedge[34]+0.75*hamil[6]*fedge[25]+0.4330127018922193*hamil[4]*fedge[19]+0.4330127018922193*hamil[6]*fedge[12]+0.75*hamil[7]*fedge[9]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[23] = 0.75*hamil[4]*fedge[60]+0.4330127018922193*hamil[4]*fedge[54]+0.75*hamil[2]*fedge[51]+0.75*hamil[7]*fedge[44]+0.4330127018922193*hamil[2]*fedge[38]+0.4330127018922193*hamil[7]*fedge[29]+0.75*hamil[6]*fedge[26]+0.4330127018922193*hamil[6]*fedge[13]; 
  Ghat[24] = 0.75*hamil[4]*fedge[62]+0.4330127018922193*hamil[4]*fedge[56]+0.75*hamil[2]*fedge[53]+0.75*hamil[7]*fedge[46]+0.4330127018922193*hamil[2]*fedge[40]+0.4330127018922193*hamil[7]*fedge[31]+0.75*hamil[6]*fedge[28]+0.4330127018922193*hamil[6]*fedge[15]; 
  Ghat[25] = 0.75*hamil[2]*fedge[55]+0.4330127018922193*hamil[2]*fedge[41]+0.75*hamil[4]*fedge[39]+0.75*hamil[6]*fedge[30]+0.4330127018922193*hamil[4]*fedge[21]+0.4330127018922193*hamil[6]*fedge[16]+0.75*hamil[7]*fedge[14]+0.4330127018922193*fedge[5]*hamil[7]; 
  Ghat[26] = 0.75*hamil[6]*fedge[63]+0.4330127018922193*hamil[6]*fedge[61]+0.75*hamil[7]*fedge[59]+0.75*hamil[2]*fedge[57]+0.4330127018922193*hamil[7]*fedge[52]+0.4330127018922193*hamil[2]*fedge[45]+0.75*hamil[4]*fedge[43]+0.4330127018922193*hamil[4]*fedge[27]; 
  Ghat[27] = 0.75*hamil[2]*fedge[58]+0.4330127018922193*hamil[2]*fedge[49]+0.75*hamil[4]*fedge[47]+0.75*hamil[6]*fedge[42]+0.4330127018922193*hamil[4]*fedge[33]+0.4330127018922193*hamil[6]*fedge[24]+0.75*hamil[7]*fedge[22]+0.4330127018922193*hamil[7]*fedge[8]; 
  Ghat[28] = 0.75*hamil[4]*fedge[63]+0.4330127018922193*hamil[4]*fedge[61]+0.75*hamil[2]*fedge[59]+0.75*hamil[7]*fedge[57]+0.4330127018922193*hamil[2]*fedge[52]+0.4330127018922193*hamil[7]*fedge[45]+0.75*hamil[6]*fedge[43]+0.4330127018922193*hamil[6]*fedge[27]; 
  Ghat[29] = 0.75*hamil[2]*fedge[60]+0.4330127018922193*hamil[2]*fedge[54]+0.75*hamil[4]*fedge[51]+0.75*hamil[6]*fedge[44]+0.4330127018922193*hamil[4]*fedge[38]+0.4330127018922193*hamil[6]*fedge[29]+0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[7]*fedge[13]; 
  Ghat[30] = 0.75*hamil[2]*fedge[62]+0.4330127018922193*hamil[2]*fedge[56]+0.75*hamil[4]*fedge[53]+0.75*hamil[6]*fedge[46]+0.4330127018922193*hamil[4]*fedge[40]+0.4330127018922193*hamil[6]*fedge[31]+0.75*hamil[7]*fedge[28]+0.4330127018922193*hamil[7]*fedge[15]; 
  Ghat[31] = 0.75*hamil[2]*fedge[63]+0.4330127018922193*hamil[2]*fedge[61]+0.75*hamil[4]*fedge[59]+0.75*hamil[6]*fedge[57]+0.4330127018922193*hamil[4]*fedge[52]+0.4330127018922193*hamil[6]*fedge[45]+0.75*hamil[7]*fedge[43]+0.4330127018922193*hamil[7]*fedge[27]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[36])+0.4330127018922193*hamil[7]*fskin[20]-0.75*hamil[6]*fskin[18]-0.75*hamil[4]*fskin[11]+0.4330127018922193*fskin[6]*hamil[6]+0.4330127018922193*fskin[4]*hamil[4]-0.75*fskin[2]*hamil[2]+0.4330127018922193*fskin[0]*hamil[2]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[48])+0.4330127018922193*hamil[7]*fskin[35]-0.75*hamil[6]*fskin[32]-0.75*hamil[4]*fskin[23]+0.4330127018922193*hamil[6]*fskin[17]+0.4330127018922193*hamil[4]*fskin[10]-0.75*hamil[2]*fskin[7]+0.4330127018922193*fskin[1]*hamil[2]; 
  Ghat[2] = -(0.75*hamil[7]*fskin[50])+0.4330127018922193*hamil[7]*fskin[37]-0.75*hamil[6]*fskin[34]-0.75*hamil[4]*fskin[25]+0.4330127018922193*hamil[6]*fskin[19]+0.4330127018922193*hamil[4]*fskin[12]-0.75*hamil[2]*fskin[9]+0.4330127018922193*hamil[2]*fskin[3]; 
  Ghat[3] = -(0.75*hamil[6]*fskin[36])+0.4330127018922193*hamil[6]*fskin[20]-0.75*hamil[7]*fskin[18]-0.75*hamil[2]*fskin[11]+0.4330127018922193*fskin[6]*hamil[7]-0.75*fskin[2]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[2]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[7]*fskin[55])+0.4330127018922193*hamil[7]*fskin[41]-0.75*hamil[6]*fskin[39]-0.75*hamil[4]*fskin[30]+0.4330127018922193*hamil[6]*fskin[21]+0.4330127018922193*hamil[4]*fskin[16]-0.75*hamil[2]*fskin[14]+0.4330127018922193*hamil[2]*fskin[5]; 
  Ghat[5] = -(0.75*hamil[4]*fskin[36])+0.4330127018922193*hamil[4]*fskin[20]-0.75*hamil[2]*fskin[18]-0.75*hamil[7]*fskin[11]+0.4330127018922193*fskin[4]*hamil[7]-0.75*fskin[2]*hamil[6]+0.4330127018922193*fskin[0]*hamil[6]+0.4330127018922193*hamil[2]*fskin[6]; 
  Ghat[6] = -(0.75*hamil[7]*fskin[58])+0.4330127018922193*hamil[7]*fskin[49]-0.75*hamil[6]*fskin[47]-0.75*hamil[4]*fskin[42]+0.4330127018922193*hamil[6]*fskin[33]+0.4330127018922193*hamil[4]*fskin[24]-0.75*hamil[2]*fskin[22]+0.4330127018922193*hamil[2]*fskin[8]; 
  Ghat[7] = -(0.75*hamil[6]*fskin[48])+0.4330127018922193*hamil[6]*fskin[35]-0.75*hamil[7]*fskin[32]-0.75*hamil[2]*fskin[23]+0.4330127018922193*hamil[7]*fskin[17]+0.4330127018922193*hamil[2]*fskin[10]-0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[1]*hamil[4]; 
  Ghat[8] = -(0.75*hamil[6]*fskin[50])+0.4330127018922193*hamil[6]*fskin[37]-0.75*hamil[7]*fskin[34]-0.75*hamil[2]*fskin[25]+0.4330127018922193*hamil[7]*fskin[19]+0.4330127018922193*hamil[2]*fskin[12]-0.75*hamil[4]*fskin[9]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[9] = -(0.75*hamil[7]*fskin[60])+0.4330127018922193*hamil[7]*fskin[54]-0.75*hamil[6]*fskin[51]-0.75*hamil[4]*fskin[44]+0.4330127018922193*hamil[6]*fskin[38]+0.4330127018922193*hamil[4]*fskin[29]-0.75*hamil[2]*fskin[26]+0.4330127018922193*hamil[2]*fskin[13]; 
  Ghat[10] = -(0.75*hamil[7]*fskin[62])+0.4330127018922193*hamil[7]*fskin[56]-0.75*hamil[6]*fskin[53]-0.75*hamil[4]*fskin[46]+0.4330127018922193*hamil[6]*fskin[40]+0.4330127018922193*hamil[4]*fskin[31]-0.75*hamil[2]*fskin[28]+0.4330127018922193*hamil[2]*fskin[15]; 
  Ghat[11] = -(0.75*hamil[6]*fskin[55])+0.4330127018922193*hamil[6]*fskin[41]-0.75*hamil[7]*fskin[39]-0.75*hamil[2]*fskin[30]+0.4330127018922193*hamil[7]*fskin[21]+0.4330127018922193*hamil[2]*fskin[16]-0.75*hamil[4]*fskin[14]+0.4330127018922193*hamil[4]*fskin[5]; 
  Ghat[12] = -(0.75*hamil[4]*fskin[48])+0.4330127018922193*hamil[4]*fskin[35]-0.75*hamil[2]*fskin[32]-0.75*hamil[7]*fskin[23]+0.4330127018922193*hamil[2]*fskin[17]+0.4330127018922193*hamil[7]*fskin[10]-0.75*hamil[6]*fskin[7]+0.4330127018922193*fskin[1]*hamil[6]; 
  Ghat[13] = -(0.75*hamil[4]*fskin[50])+0.4330127018922193*hamil[4]*fskin[37]-0.75*hamil[2]*fskin[34]-0.75*hamil[7]*fskin[25]+0.4330127018922193*hamil[2]*fskin[19]+0.4330127018922193*hamil[7]*fskin[12]-0.75*hamil[6]*fskin[9]+0.4330127018922193*fskin[3]*hamil[6]; 
  Ghat[14] = -(0.75*hamil[2]*fskin[36])+0.4330127018922193*hamil[2]*fskin[20]-0.75*hamil[4]*fskin[18]-0.75*hamil[6]*fskin[11]-0.75*fskin[2]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*fskin[4]*hamil[6]+0.4330127018922193*hamil[4]*fskin[6]; 
  Ghat[15] = -(0.75*hamil[4]*fskin[55])+0.4330127018922193*hamil[4]*fskin[41]-0.75*hamil[2]*fskin[39]-0.75*hamil[7]*fskin[30]+0.4330127018922193*hamil[2]*fskin[21]+0.4330127018922193*hamil[7]*fskin[16]-0.75*hamil[6]*fskin[14]+0.4330127018922193*fskin[5]*hamil[6]; 
  Ghat[16] = -(0.75*hamil[6]*fskin[58])+0.4330127018922193*hamil[6]*fskin[49]-0.75*hamil[7]*fskin[47]-0.75*hamil[2]*fskin[42]+0.4330127018922193*hamil[7]*fskin[33]+0.4330127018922193*hamil[2]*fskin[24]-0.75*hamil[4]*fskin[22]+0.4330127018922193*hamil[4]*fskin[8]; 
  Ghat[17] = -(0.75*hamil[7]*fskin[63])+0.4330127018922193*hamil[7]*fskin[61]-0.75*hamil[6]*fskin[59]-0.75*hamil[4]*fskin[57]+0.4330127018922193*hamil[6]*fskin[52]+0.4330127018922193*hamil[4]*fskin[45]-0.75*hamil[2]*fskin[43]+0.4330127018922193*hamil[2]*fskin[27]; 
  Ghat[18] = -(0.75*hamil[6]*fskin[60])+0.4330127018922193*hamil[6]*fskin[54]-0.75*hamil[7]*fskin[51]-0.75*hamil[2]*fskin[44]+0.4330127018922193*hamil[7]*fskin[38]+0.4330127018922193*hamil[2]*fskin[29]-0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[4]*fskin[13]; 
  Ghat[19] = -(0.75*hamil[6]*fskin[62])+0.4330127018922193*hamil[6]*fskin[56]-0.75*hamil[7]*fskin[53]-0.75*hamil[2]*fskin[46]+0.4330127018922193*hamil[7]*fskin[40]+0.4330127018922193*hamil[2]*fskin[31]-0.75*hamil[4]*fskin[28]+0.4330127018922193*hamil[4]*fskin[15]; 
  Ghat[20] = -(0.75*hamil[4]*fskin[58])+0.4330127018922193*hamil[4]*fskin[49]-0.75*hamil[2]*fskin[47]-0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[2]*fskin[33]+0.4330127018922193*hamil[7]*fskin[24]-0.75*hamil[6]*fskin[22]+0.4330127018922193*hamil[6]*fskin[8]; 
  Ghat[21] = -(0.75*hamil[2]*fskin[48])+0.4330127018922193*hamil[2]*fskin[35]-0.75*hamil[4]*fskin[32]-0.75*hamil[6]*fskin[23]+0.4330127018922193*hamil[4]*fskin[17]+0.4330127018922193*hamil[6]*fskin[10]-0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[1]*hamil[7]; 
  Ghat[22] = -(0.75*hamil[2]*fskin[50])+0.4330127018922193*hamil[2]*fskin[37]-0.75*hamil[4]*fskin[34]-0.75*hamil[6]*fskin[25]+0.4330127018922193*hamil[4]*fskin[19]+0.4330127018922193*hamil[6]*fskin[12]-0.75*hamil[7]*fskin[9]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[23] = -(0.75*hamil[4]*fskin[60])+0.4330127018922193*hamil[4]*fskin[54]-0.75*hamil[2]*fskin[51]-0.75*hamil[7]*fskin[44]+0.4330127018922193*hamil[2]*fskin[38]+0.4330127018922193*hamil[7]*fskin[29]-0.75*hamil[6]*fskin[26]+0.4330127018922193*hamil[6]*fskin[13]; 
  Ghat[24] = -(0.75*hamil[4]*fskin[62])+0.4330127018922193*hamil[4]*fskin[56]-0.75*hamil[2]*fskin[53]-0.75*hamil[7]*fskin[46]+0.4330127018922193*hamil[2]*fskin[40]+0.4330127018922193*hamil[7]*fskin[31]-0.75*hamil[6]*fskin[28]+0.4330127018922193*hamil[6]*fskin[15]; 
  Ghat[25] = -(0.75*hamil[2]*fskin[55])+0.4330127018922193*hamil[2]*fskin[41]-0.75*hamil[4]*fskin[39]-0.75*hamil[6]*fskin[30]+0.4330127018922193*hamil[4]*fskin[21]+0.4330127018922193*hamil[6]*fskin[16]-0.75*hamil[7]*fskin[14]+0.4330127018922193*fskin[5]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[6]*fskin[63])+0.4330127018922193*hamil[6]*fskin[61]-0.75*hamil[7]*fskin[59]-0.75*hamil[2]*fskin[57]+0.4330127018922193*hamil[7]*fskin[52]+0.4330127018922193*hamil[2]*fskin[45]-0.75*hamil[4]*fskin[43]+0.4330127018922193*hamil[4]*fskin[27]; 
  Ghat[27] = -(0.75*hamil[2]*fskin[58])+0.4330127018922193*hamil[2]*fskin[49]-0.75*hamil[4]*fskin[47]-0.75*hamil[6]*fskin[42]+0.4330127018922193*hamil[4]*fskin[33]+0.4330127018922193*hamil[6]*fskin[24]-0.75*hamil[7]*fskin[22]+0.4330127018922193*hamil[7]*fskin[8]; 
  Ghat[28] = -(0.75*hamil[4]*fskin[63])+0.4330127018922193*hamil[4]*fskin[61]-0.75*hamil[2]*fskin[59]-0.75*hamil[7]*fskin[57]+0.4330127018922193*hamil[2]*fskin[52]+0.4330127018922193*hamil[7]*fskin[45]-0.75*hamil[6]*fskin[43]+0.4330127018922193*hamil[6]*fskin[27]; 
  Ghat[29] = -(0.75*hamil[2]*fskin[60])+0.4330127018922193*hamil[2]*fskin[54]-0.75*hamil[4]*fskin[51]-0.75*hamil[6]*fskin[44]+0.4330127018922193*hamil[4]*fskin[38]+0.4330127018922193*hamil[6]*fskin[29]-0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[7]*fskin[13]; 
  Ghat[30] = -(0.75*hamil[2]*fskin[62])+0.4330127018922193*hamil[2]*fskin[56]-0.75*hamil[4]*fskin[53]-0.75*hamil[6]*fskin[46]+0.4330127018922193*hamil[4]*fskin[40]+0.4330127018922193*hamil[6]*fskin[31]-0.75*hamil[7]*fskin[28]+0.4330127018922193*hamil[7]*fskin[15]; 
  Ghat[31] = -(0.75*hamil[2]*fskin[63])+0.4330127018922193*hamil[2]*fskin[61]-0.75*hamil[4]*fskin[59]-0.75*hamil[6]*fskin[57]+0.4330127018922193*hamil[4]*fskin[52]+0.4330127018922193*hamil[6]*fskin[45]-0.75*hamil[7]*fskin[43]+0.4330127018922193*hamil[7]*fskin[27]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv11*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11*dx11; 
  out[2] += -(1.224744871391589*Ghat[0]*dv11*dx11); 
  out[3] += 0.7071067811865475*Ghat[2]*dv11*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11*dx11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11*dx11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11*dx11; 
  out[7] += -(1.224744871391589*Ghat[1]*dv11*dx11); 
  out[8] += 0.7071067811865475*Ghat[6]*dv11*dx11; 
  out[9] += -(1.224744871391589*Ghat[2]*dv11*dx11); 
  out[10] += 0.7071067811865475*Ghat[7]*dv11*dx11; 
  out[11] += -(1.224744871391589*Ghat[3]*dv11*dx11); 
  out[12] += 0.7071067811865475*Ghat[8]*dv11*dx11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11*dx11; 
  out[14] += -(1.224744871391589*Ghat[4]*dv11*dx11); 
  out[15] += 0.7071067811865475*Ghat[10]*dv11*dx11; 
  out[16] += 0.7071067811865475*Ghat[11]*dv11*dx11; 
  out[17] += 0.7071067811865475*Ghat[12]*dv11*dx11; 
  out[18] += -(1.224744871391589*Ghat[5]*dv11*dx11); 
  out[19] += 0.7071067811865475*Ghat[13]*dv11*dx11; 
  out[20] += 0.7071067811865475*Ghat[14]*dv11*dx11; 
  out[21] += 0.7071067811865475*Ghat[15]*dv11*dx11; 
  out[22] += -(1.224744871391589*Ghat[6]*dv11*dx11); 
  out[23] += -(1.224744871391589*Ghat[7]*dv11*dx11); 
  out[24] += 0.7071067811865475*Ghat[16]*dv11*dx11; 
  out[25] += -(1.224744871391589*Ghat[8]*dv11*dx11); 
  out[26] += -(1.224744871391589*Ghat[9]*dv11*dx11); 
  out[27] += 0.7071067811865475*Ghat[17]*dv11*dx11; 
  out[28] += -(1.224744871391589*Ghat[10]*dv11*dx11); 
  out[29] += 0.7071067811865475*Ghat[18]*dv11*dx11; 
  out[30] += -(1.224744871391589*Ghat[11]*dv11*dx11); 
  out[31] += 0.7071067811865475*Ghat[19]*dv11*dx11; 
  out[32] += -(1.224744871391589*Ghat[12]*dv11*dx11); 
  out[33] += 0.7071067811865475*Ghat[20]*dv11*dx11; 
  out[34] += -(1.224744871391589*Ghat[13]*dv11*dx11); 
  out[35] += 0.7071067811865475*Ghat[21]*dv11*dx11; 
  out[36] += -(1.224744871391589*Ghat[14]*dv11*dx11); 
  out[37] += 0.7071067811865475*Ghat[22]*dv11*dx11; 
  out[38] += 0.7071067811865475*Ghat[23]*dv11*dx11; 
  out[39] += -(1.224744871391589*Ghat[15]*dv11*dx11); 
  out[40] += 0.7071067811865475*Ghat[24]*dv11*dx11; 
  out[41] += 0.7071067811865475*Ghat[25]*dv11*dx11; 
  out[42] += -(1.224744871391589*Ghat[16]*dv11*dx11); 
  out[43] += -(1.224744871391589*Ghat[17]*dv11*dx11); 
  out[44] += -(1.224744871391589*Ghat[18]*dv11*dx11); 
  out[45] += 0.7071067811865475*Ghat[26]*dv11*dx11; 
  out[46] += -(1.224744871391589*Ghat[19]*dv11*dx11); 
  out[47] += -(1.224744871391589*Ghat[20]*dv11*dx11); 
  out[48] += -(1.224744871391589*Ghat[21]*dv11*dx11); 
  out[49] += 0.7071067811865475*Ghat[27]*dv11*dx11; 
  out[50] += -(1.224744871391589*Ghat[22]*dv11*dx11); 
  out[51] += -(1.224744871391589*Ghat[23]*dv11*dx11); 
  out[52] += 0.7071067811865475*Ghat[28]*dv11*dx11; 
  out[53] += -(1.224744871391589*Ghat[24]*dv11*dx11); 
  out[54] += 0.7071067811865475*Ghat[29]*dv11*dx11; 
  out[55] += -(1.224744871391589*Ghat[25]*dv11*dx11); 
  out[56] += 0.7071067811865475*Ghat[30]*dv11*dx11; 
  out[57] += -(1.224744871391589*Ghat[26]*dv11*dx11); 
  out[58] += -(1.224744871391589*Ghat[27]*dv11*dx11); 
  out[59] += -(1.224744871391589*Ghat[28]*dv11*dx11); 
  out[60] += -(1.224744871391589*Ghat[29]*dv11*dx11); 
  out[61] += 0.7071067811865475*Ghat[31]*dv11*dx11; 
  out[62] += -(1.224744871391589*Ghat[30]*dv11*dx11); 
  out[63] += -(1.224744871391589*Ghat[31]*dv11*dx11); 

  } 
  return 0.0;

} 
