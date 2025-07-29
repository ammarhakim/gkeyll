#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[3]; 
  double wv = w[3]; 

  double Ghat[32]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[38]+0.4330127018922193*hamil[7]*fskin[21]+0.75*hamil[5]*fskin[17]+0.75*hamil[4]*fskin[13]+0.4330127018922193*hamil[5]*fskin[6]+0.4330127018922193*hamil[4]*fskin[5]+0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fskin[51]+0.4330127018922193*hamil[7]*fskin[39]+0.75*hamil[5]*fskin[32]+0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[5]*fskin[18]+0.4330127018922193*hamil[4]*fskin[14]+0.75*hamil[1]*fskin[7]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = 0.75*hamil[7]*fskin[52]+0.4330127018922193*hamil[7]*fskin[40]+0.75*hamil[5]*fskin[33]+0.75*hamil[4]*fskin[27]+0.4330127018922193*hamil[5]*fskin[19]+0.4330127018922193*hamil[4]*fskin[15]+0.75*hamil[1]*fskin[8]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = 0.75*hamil[7]*fskin[54]+0.4330127018922193*hamil[7]*fskin[41]+0.75*hamil[5]*fskin[35]+0.75*hamil[4]*fskin[29]+0.4330127018922193*hamil[5]*fskin[20]+0.4330127018922193*hamil[4]*fskin[16]+0.75*hamil[1]*fskin[10]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = 0.75*hamil[5]*fskin[38]+0.4330127018922193*hamil[5]*fskin[21]+0.75*hamil[7]*fskin[17]+0.75*hamil[1]*fskin[13]+0.4330127018922193*fskin[6]*hamil[7]+0.4330127018922193*hamil[1]*fskin[5]+0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]; 
  Ghat[5] = 0.75*hamil[4]*fskin[38]+0.4330127018922193*hamil[4]*fskin[21]+0.75*hamil[1]*fskin[17]+0.75*hamil[7]*fskin[13]+0.4330127018922193*fskin[5]*hamil[7]+0.4330127018922193*hamil[1]*fskin[6]+0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]; 
  Ghat[6] = 0.75*hamil[7]*fskin[59]+0.4330127018922193*hamil[7]*fskin[53]+0.75*hamil[5]*fskin[47]+0.75*hamil[4]*fskin[43]+0.4330127018922193*hamil[5]*fskin[34]+0.4330127018922193*hamil[4]*fskin[28]+0.75*hamil[1]*fskin[22]+0.4330127018922193*hamil[1]*fskin[9]; 
  Ghat[7] = 0.75*hamil[7]*fskin[60]+0.4330127018922193*hamil[7]*fskin[55]+0.75*hamil[5]*fskin[48]+0.75*hamil[4]*fskin[44]+0.4330127018922193*hamil[5]*fskin[36]+0.4330127018922193*hamil[4]*fskin[30]+0.75*hamil[1]*fskin[23]+0.4330127018922193*hamil[1]*fskin[11]; 
  Ghat[8] = 0.75*hamil[7]*fskin[61]+0.4330127018922193*hamil[7]*fskin[56]+0.75*hamil[5]*fskin[49]+0.75*hamil[4]*fskin[45]+0.4330127018922193*hamil[5]*fskin[37]+0.4330127018922193*hamil[4]*fskin[31]+0.75*hamil[1]*fskin[24]+0.4330127018922193*hamil[1]*fskin[12]; 
  Ghat[9] = 0.75*hamil[5]*fskin[51]+0.4330127018922193*hamil[5]*fskin[39]+0.75*hamil[7]*fskin[32]+0.75*hamil[1]*fskin[26]+0.4330127018922193*hamil[7]*fskin[18]+0.4330127018922193*hamil[1]*fskin[14]+0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[10] = 0.75*hamil[5]*fskin[52]+0.4330127018922193*hamil[5]*fskin[40]+0.75*hamil[7]*fskin[33]+0.75*hamil[1]*fskin[27]+0.4330127018922193*hamil[7]*fskin[19]+0.4330127018922193*hamil[1]*fskin[15]+0.75*hamil[4]*fskin[8]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[11] = 0.75*hamil[5]*fskin[54]+0.4330127018922193*hamil[5]*fskin[41]+0.75*hamil[7]*fskin[35]+0.75*hamil[1]*fskin[29]+0.4330127018922193*hamil[7]*fskin[20]+0.4330127018922193*hamil[1]*fskin[16]+0.75*hamil[4]*fskin[10]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[12] = 0.75*hamil[4]*fskin[51]+0.4330127018922193*hamil[4]*fskin[39]+0.75*hamil[1]*fskin[32]+0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[1]*fskin[18]+0.4330127018922193*hamil[7]*fskin[14]+0.75*hamil[5]*fskin[7]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[13] = 0.75*hamil[4]*fskin[52]+0.4330127018922193*hamil[4]*fskin[40]+0.75*hamil[1]*fskin[33]+0.75*hamil[7]*fskin[27]+0.4330127018922193*hamil[1]*fskin[19]+0.4330127018922193*hamil[7]*fskin[15]+0.75*hamil[5]*fskin[8]+0.4330127018922193*fskin[3]*hamil[5]; 
  Ghat[14] = 0.75*hamil[4]*fskin[54]+0.4330127018922193*hamil[4]*fskin[41]+0.75*hamil[1]*fskin[35]+0.75*hamil[7]*fskin[29]+0.4330127018922193*hamil[1]*fskin[20]+0.4330127018922193*hamil[7]*fskin[16]+0.75*hamil[5]*fskin[10]+0.4330127018922193*fskin[4]*hamil[5]; 
  Ghat[15] = 0.75*hamil[1]*fskin[38]+0.4330127018922193*hamil[1]*fskin[21]+0.75*hamil[4]*fskin[17]+0.75*hamil[5]*fskin[13]+0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*hamil[4]*fskin[6]+0.4330127018922193*fskin[5]*hamil[5]; 
  Ghat[16] = 0.75*hamil[7]*fskin[63]+0.4330127018922193*hamil[7]*fskin[62]+0.75*hamil[5]*fskin[58]+0.75*hamil[4]*fskin[57]+0.4330127018922193*hamil[5]*fskin[50]+0.4330127018922193*hamil[4]*fskin[46]+0.75*hamil[1]*fskin[42]+0.4330127018922193*hamil[1]*fskin[25]; 
  Ghat[17] = 0.75*hamil[5]*fskin[59]+0.4330127018922193*hamil[5]*fskin[53]+0.75*hamil[7]*fskin[47]+0.75*hamil[1]*fskin[43]+0.4330127018922193*hamil[7]*fskin[34]+0.4330127018922193*hamil[1]*fskin[28]+0.75*hamil[4]*fskin[22]+0.4330127018922193*hamil[4]*fskin[9]; 
  Ghat[18] = 0.75*hamil[5]*fskin[60]+0.4330127018922193*hamil[5]*fskin[55]+0.75*hamil[7]*fskin[48]+0.75*hamil[1]*fskin[44]+0.4330127018922193*hamil[7]*fskin[36]+0.4330127018922193*hamil[1]*fskin[30]+0.75*hamil[4]*fskin[23]+0.4330127018922193*hamil[4]*fskin[11]; 
  Ghat[19] = 0.75*hamil[5]*fskin[61]+0.4330127018922193*hamil[5]*fskin[56]+0.75*hamil[7]*fskin[49]+0.75*hamil[1]*fskin[45]+0.4330127018922193*hamil[7]*fskin[37]+0.4330127018922193*hamil[1]*fskin[31]+0.75*hamil[4]*fskin[24]+0.4330127018922193*hamil[4]*fskin[12]; 
  Ghat[20] = 0.75*hamil[4]*fskin[59]+0.4330127018922193*hamil[4]*fskin[53]+0.75*hamil[1]*fskin[47]+0.75*hamil[7]*fskin[43]+0.4330127018922193*hamil[1]*fskin[34]+0.4330127018922193*hamil[7]*fskin[28]+0.75*hamil[5]*fskin[22]+0.4330127018922193*hamil[5]*fskin[9]; 
  Ghat[21] = 0.75*hamil[4]*fskin[60]+0.4330127018922193*hamil[4]*fskin[55]+0.75*hamil[1]*fskin[48]+0.75*hamil[7]*fskin[44]+0.4330127018922193*hamil[1]*fskin[36]+0.4330127018922193*hamil[7]*fskin[30]+0.75*hamil[5]*fskin[23]+0.4330127018922193*hamil[5]*fskin[11]; 
  Ghat[22] = 0.75*hamil[4]*fskin[61]+0.4330127018922193*hamil[4]*fskin[56]+0.75*hamil[1]*fskin[49]+0.75*hamil[7]*fskin[45]+0.4330127018922193*hamil[1]*fskin[37]+0.4330127018922193*hamil[7]*fskin[31]+0.75*hamil[5]*fskin[24]+0.4330127018922193*hamil[5]*fskin[12]; 
  Ghat[23] = 0.75*hamil[1]*fskin[51]+0.4330127018922193*hamil[1]*fskin[39]+0.75*hamil[4]*fskin[32]+0.75*hamil[5]*fskin[26]+0.4330127018922193*hamil[4]*fskin[18]+0.4330127018922193*hamil[5]*fskin[14]+0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[24] = 0.75*hamil[1]*fskin[52]+0.4330127018922193*hamil[1]*fskin[40]+0.75*hamil[4]*fskin[33]+0.75*hamil[5]*fskin[27]+0.4330127018922193*hamil[4]*fskin[19]+0.4330127018922193*hamil[5]*fskin[15]+0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[25] = 0.75*hamil[1]*fskin[54]+0.4330127018922193*hamil[1]*fskin[41]+0.75*hamil[4]*fskin[35]+0.75*hamil[5]*fskin[29]+0.4330127018922193*hamil[4]*fskin[20]+0.4330127018922193*hamil[5]*fskin[16]+0.75*hamil[7]*fskin[10]+0.4330127018922193*fskin[4]*hamil[7]; 
  Ghat[26] = 0.75*hamil[5]*fskin[63]+0.4330127018922193*hamil[5]*fskin[62]+0.75*hamil[7]*fskin[58]+0.75*hamil[1]*fskin[57]+0.4330127018922193*hamil[7]*fskin[50]+0.4330127018922193*hamil[1]*fskin[46]+0.75*hamil[4]*fskin[42]+0.4330127018922193*hamil[4]*fskin[25]; 
  Ghat[27] = 0.75*hamil[4]*fskin[63]+0.4330127018922193*hamil[4]*fskin[62]+0.75*hamil[1]*fskin[58]+0.75*hamil[7]*fskin[57]+0.4330127018922193*hamil[1]*fskin[50]+0.4330127018922193*hamil[7]*fskin[46]+0.75*hamil[5]*fskin[42]+0.4330127018922193*hamil[5]*fskin[25]; 
  Ghat[28] = 0.75*hamil[1]*fskin[59]+0.4330127018922193*hamil[1]*fskin[53]+0.75*hamil[4]*fskin[47]+0.75*hamil[5]*fskin[43]+0.4330127018922193*hamil[4]*fskin[34]+0.4330127018922193*hamil[5]*fskin[28]+0.75*hamil[7]*fskin[22]+0.4330127018922193*hamil[7]*fskin[9]; 
  Ghat[29] = 0.75*hamil[1]*fskin[60]+0.4330127018922193*hamil[1]*fskin[55]+0.75*hamil[4]*fskin[48]+0.75*hamil[5]*fskin[44]+0.4330127018922193*hamil[4]*fskin[36]+0.4330127018922193*hamil[5]*fskin[30]+0.75*hamil[7]*fskin[23]+0.4330127018922193*hamil[7]*fskin[11]; 
  Ghat[30] = 0.75*hamil[1]*fskin[61]+0.4330127018922193*hamil[1]*fskin[56]+0.75*hamil[4]*fskin[49]+0.75*hamil[5]*fskin[45]+0.4330127018922193*hamil[4]*fskin[37]+0.4330127018922193*hamil[5]*fskin[31]+0.75*hamil[7]*fskin[24]+0.4330127018922193*hamil[7]*fskin[12]; 
  Ghat[31] = 0.75*hamil[1]*fskin[63]+0.4330127018922193*hamil[1]*fskin[62]+0.75*hamil[4]*fskin[58]+0.75*hamil[5]*fskin[57]+0.4330127018922193*hamil[4]*fskin[50]+0.4330127018922193*hamil[5]*fskin[46]+0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[7]*fskin[25]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[38])+0.4330127018922193*hamil[7]*fedge[21]-0.75*hamil[5]*fedge[17]-0.75*hamil[4]*fedge[13]+0.4330127018922193*hamil[5]*fedge[6]+0.4330127018922193*hamil[4]*fedge[5]-0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[51])+0.4330127018922193*hamil[7]*fedge[39]-0.75*hamil[5]*fedge[32]-0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[5]*fedge[18]+0.4330127018922193*hamil[4]*fedge[14]-0.75*hamil[1]*fedge[7]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = -(0.75*hamil[7]*fedge[52])+0.4330127018922193*hamil[7]*fedge[40]-0.75*hamil[5]*fedge[33]-0.75*hamil[4]*fedge[27]+0.4330127018922193*hamil[5]*fedge[19]+0.4330127018922193*hamil[4]*fedge[15]-0.75*hamil[1]*fedge[8]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = -(0.75*hamil[7]*fedge[54])+0.4330127018922193*hamil[7]*fedge[41]-0.75*hamil[5]*fedge[35]-0.75*hamil[4]*fedge[29]+0.4330127018922193*hamil[5]*fedge[20]+0.4330127018922193*hamil[4]*fedge[16]-0.75*hamil[1]*fedge[10]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[5]*fedge[38])+0.4330127018922193*hamil[5]*fedge[21]-0.75*hamil[7]*fedge[17]-0.75*hamil[1]*fedge[13]+0.4330127018922193*fedge[6]*hamil[7]+0.4330127018922193*hamil[1]*fedge[5]-0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]; 
  Ghat[5] = -(0.75*hamil[4]*fedge[38])+0.4330127018922193*hamil[4]*fedge[21]-0.75*hamil[1]*fedge[17]-0.75*hamil[7]*fedge[13]+0.4330127018922193*fedge[5]*hamil[7]+0.4330127018922193*hamil[1]*fedge[6]-0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]; 
  Ghat[6] = -(0.75*hamil[7]*fedge[59])+0.4330127018922193*hamil[7]*fedge[53]-0.75*hamil[5]*fedge[47]-0.75*hamil[4]*fedge[43]+0.4330127018922193*hamil[5]*fedge[34]+0.4330127018922193*hamil[4]*fedge[28]-0.75*hamil[1]*fedge[22]+0.4330127018922193*hamil[1]*fedge[9]; 
  Ghat[7] = -(0.75*hamil[7]*fedge[60])+0.4330127018922193*hamil[7]*fedge[55]-0.75*hamil[5]*fedge[48]-0.75*hamil[4]*fedge[44]+0.4330127018922193*hamil[5]*fedge[36]+0.4330127018922193*hamil[4]*fedge[30]-0.75*hamil[1]*fedge[23]+0.4330127018922193*hamil[1]*fedge[11]; 
  Ghat[8] = -(0.75*hamil[7]*fedge[61])+0.4330127018922193*hamil[7]*fedge[56]-0.75*hamil[5]*fedge[49]-0.75*hamil[4]*fedge[45]+0.4330127018922193*hamil[5]*fedge[37]+0.4330127018922193*hamil[4]*fedge[31]-0.75*hamil[1]*fedge[24]+0.4330127018922193*hamil[1]*fedge[12]; 
  Ghat[9] = -(0.75*hamil[5]*fedge[51])+0.4330127018922193*hamil[5]*fedge[39]-0.75*hamil[7]*fedge[32]-0.75*hamil[1]*fedge[26]+0.4330127018922193*hamil[7]*fedge[18]+0.4330127018922193*hamil[1]*fedge[14]-0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[10] = -(0.75*hamil[5]*fedge[52])+0.4330127018922193*hamil[5]*fedge[40]-0.75*hamil[7]*fedge[33]-0.75*hamil[1]*fedge[27]+0.4330127018922193*hamil[7]*fedge[19]+0.4330127018922193*hamil[1]*fedge[15]-0.75*hamil[4]*fedge[8]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[11] = -(0.75*hamil[5]*fedge[54])+0.4330127018922193*hamil[5]*fedge[41]-0.75*hamil[7]*fedge[35]-0.75*hamil[1]*fedge[29]+0.4330127018922193*hamil[7]*fedge[20]+0.4330127018922193*hamil[1]*fedge[16]-0.75*hamil[4]*fedge[10]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[12] = -(0.75*hamil[4]*fedge[51])+0.4330127018922193*hamil[4]*fedge[39]-0.75*hamil[1]*fedge[32]-0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[1]*fedge[18]+0.4330127018922193*hamil[7]*fedge[14]-0.75*hamil[5]*fedge[7]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[13] = -(0.75*hamil[4]*fedge[52])+0.4330127018922193*hamil[4]*fedge[40]-0.75*hamil[1]*fedge[33]-0.75*hamil[7]*fedge[27]+0.4330127018922193*hamil[1]*fedge[19]+0.4330127018922193*hamil[7]*fedge[15]-0.75*hamil[5]*fedge[8]+0.4330127018922193*fedge[3]*hamil[5]; 
  Ghat[14] = -(0.75*hamil[4]*fedge[54])+0.4330127018922193*hamil[4]*fedge[41]-0.75*hamil[1]*fedge[35]-0.75*hamil[7]*fedge[29]+0.4330127018922193*hamil[1]*fedge[20]+0.4330127018922193*hamil[7]*fedge[16]-0.75*hamil[5]*fedge[10]+0.4330127018922193*fedge[4]*hamil[5]; 
  Ghat[15] = -(0.75*hamil[1]*fedge[38])+0.4330127018922193*hamil[1]*fedge[21]-0.75*hamil[4]*fedge[17]-0.75*hamil[5]*fedge[13]-0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*hamil[4]*fedge[6]+0.4330127018922193*fedge[5]*hamil[5]; 
  Ghat[16] = -(0.75*hamil[7]*fedge[63])+0.4330127018922193*hamil[7]*fedge[62]-0.75*hamil[5]*fedge[58]-0.75*hamil[4]*fedge[57]+0.4330127018922193*hamil[5]*fedge[50]+0.4330127018922193*hamil[4]*fedge[46]-0.75*hamil[1]*fedge[42]+0.4330127018922193*hamil[1]*fedge[25]; 
  Ghat[17] = -(0.75*hamil[5]*fedge[59])+0.4330127018922193*hamil[5]*fedge[53]-0.75*hamil[7]*fedge[47]-0.75*hamil[1]*fedge[43]+0.4330127018922193*hamil[7]*fedge[34]+0.4330127018922193*hamil[1]*fedge[28]-0.75*hamil[4]*fedge[22]+0.4330127018922193*hamil[4]*fedge[9]; 
  Ghat[18] = -(0.75*hamil[5]*fedge[60])+0.4330127018922193*hamil[5]*fedge[55]-0.75*hamil[7]*fedge[48]-0.75*hamil[1]*fedge[44]+0.4330127018922193*hamil[7]*fedge[36]+0.4330127018922193*hamil[1]*fedge[30]-0.75*hamil[4]*fedge[23]+0.4330127018922193*hamil[4]*fedge[11]; 
  Ghat[19] = -(0.75*hamil[5]*fedge[61])+0.4330127018922193*hamil[5]*fedge[56]-0.75*hamil[7]*fedge[49]-0.75*hamil[1]*fedge[45]+0.4330127018922193*hamil[7]*fedge[37]+0.4330127018922193*hamil[1]*fedge[31]-0.75*hamil[4]*fedge[24]+0.4330127018922193*hamil[4]*fedge[12]; 
  Ghat[20] = -(0.75*hamil[4]*fedge[59])+0.4330127018922193*hamil[4]*fedge[53]-0.75*hamil[1]*fedge[47]-0.75*hamil[7]*fedge[43]+0.4330127018922193*hamil[1]*fedge[34]+0.4330127018922193*hamil[7]*fedge[28]-0.75*hamil[5]*fedge[22]+0.4330127018922193*hamil[5]*fedge[9]; 
  Ghat[21] = -(0.75*hamil[4]*fedge[60])+0.4330127018922193*hamil[4]*fedge[55]-0.75*hamil[1]*fedge[48]-0.75*hamil[7]*fedge[44]+0.4330127018922193*hamil[1]*fedge[36]+0.4330127018922193*hamil[7]*fedge[30]-0.75*hamil[5]*fedge[23]+0.4330127018922193*hamil[5]*fedge[11]; 
  Ghat[22] = -(0.75*hamil[4]*fedge[61])+0.4330127018922193*hamil[4]*fedge[56]-0.75*hamil[1]*fedge[49]-0.75*hamil[7]*fedge[45]+0.4330127018922193*hamil[1]*fedge[37]+0.4330127018922193*hamil[7]*fedge[31]-0.75*hamil[5]*fedge[24]+0.4330127018922193*hamil[5]*fedge[12]; 
  Ghat[23] = -(0.75*hamil[1]*fedge[51])+0.4330127018922193*hamil[1]*fedge[39]-0.75*hamil[4]*fedge[32]-0.75*hamil[5]*fedge[26]+0.4330127018922193*hamil[4]*fedge[18]+0.4330127018922193*hamil[5]*fedge[14]-0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[24] = -(0.75*hamil[1]*fedge[52])+0.4330127018922193*hamil[1]*fedge[40]-0.75*hamil[4]*fedge[33]-0.75*hamil[5]*fedge[27]+0.4330127018922193*hamil[4]*fedge[19]+0.4330127018922193*hamil[5]*fedge[15]-0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[25] = -(0.75*hamil[1]*fedge[54])+0.4330127018922193*hamil[1]*fedge[41]-0.75*hamil[4]*fedge[35]-0.75*hamil[5]*fedge[29]+0.4330127018922193*hamil[4]*fedge[20]+0.4330127018922193*hamil[5]*fedge[16]-0.75*hamil[7]*fedge[10]+0.4330127018922193*fedge[4]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[5]*fedge[63])+0.4330127018922193*hamil[5]*fedge[62]-0.75*hamil[7]*fedge[58]-0.75*hamil[1]*fedge[57]+0.4330127018922193*hamil[7]*fedge[50]+0.4330127018922193*hamil[1]*fedge[46]-0.75*hamil[4]*fedge[42]+0.4330127018922193*hamil[4]*fedge[25]; 
  Ghat[27] = -(0.75*hamil[4]*fedge[63])+0.4330127018922193*hamil[4]*fedge[62]-0.75*hamil[1]*fedge[58]-0.75*hamil[7]*fedge[57]+0.4330127018922193*hamil[1]*fedge[50]+0.4330127018922193*hamil[7]*fedge[46]-0.75*hamil[5]*fedge[42]+0.4330127018922193*hamil[5]*fedge[25]; 
  Ghat[28] = -(0.75*hamil[1]*fedge[59])+0.4330127018922193*hamil[1]*fedge[53]-0.75*hamil[4]*fedge[47]-0.75*hamil[5]*fedge[43]+0.4330127018922193*hamil[4]*fedge[34]+0.4330127018922193*hamil[5]*fedge[28]-0.75*hamil[7]*fedge[22]+0.4330127018922193*hamil[7]*fedge[9]; 
  Ghat[29] = -(0.75*hamil[1]*fedge[60])+0.4330127018922193*hamil[1]*fedge[55]-0.75*hamil[4]*fedge[48]-0.75*hamil[5]*fedge[44]+0.4330127018922193*hamil[4]*fedge[36]+0.4330127018922193*hamil[5]*fedge[30]-0.75*hamil[7]*fedge[23]+0.4330127018922193*hamil[7]*fedge[11]; 
  Ghat[30] = -(0.75*hamil[1]*fedge[61])+0.4330127018922193*hamil[1]*fedge[56]-0.75*hamil[4]*fedge[49]-0.75*hamil[5]*fedge[45]+0.4330127018922193*hamil[4]*fedge[37]+0.4330127018922193*hamil[5]*fedge[31]-0.75*hamil[7]*fedge[24]+0.4330127018922193*hamil[7]*fedge[12]; 
  Ghat[31] = -(0.75*hamil[1]*fedge[63])+0.4330127018922193*hamil[1]*fedge[62]-0.75*hamil[4]*fedge[58]-0.75*hamil[5]*fedge[57]+0.4330127018922193*hamil[4]*fedge[50]+0.4330127018922193*hamil[5]*fedge[46]-0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[7]*fedge[25]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[8] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[10] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv10*dx10); 
  out[13] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[9]*dv10*dx10); 
  out[15] += -(0.7071067811865475*Ghat[10]*dv10*dx10); 
  out[16] += -(0.7071067811865475*Ghat[11]*dv10*dx10); 
  out[17] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[18] += -(0.7071067811865475*Ghat[12]*dv10*dx10); 
  out[19] += -(0.7071067811865475*Ghat[13]*dv10*dx10); 
  out[20] += -(0.7071067811865475*Ghat[14]*dv10*dx10); 
  out[21] += -(0.7071067811865475*Ghat[15]*dv10*dx10); 
  out[22] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[23] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[24] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[25] += -(0.7071067811865475*Ghat[16]*dv10*dx10); 
  out[26] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[27] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[28] += -(0.7071067811865475*Ghat[17]*dv10*dx10); 
  out[29] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[30] += -(0.7071067811865475*Ghat[18]*dv10*dx10); 
  out[31] += -(0.7071067811865475*Ghat[19]*dv10*dx10); 
  out[32] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[33] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[34] += -(0.7071067811865475*Ghat[20]*dv10*dx10); 
  out[35] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[36] += -(0.7071067811865475*Ghat[21]*dv10*dx10); 
  out[37] += -(0.7071067811865475*Ghat[22]*dv10*dx10); 
  out[38] += -(1.224744871391589*Ghat[15]*dv10*dx10); 
  out[39] += -(0.7071067811865475*Ghat[23]*dv10*dx10); 
  out[40] += -(0.7071067811865475*Ghat[24]*dv10*dx10); 
  out[41] += -(0.7071067811865475*Ghat[25]*dv10*dx10); 
  out[42] += -(1.224744871391589*Ghat[16]*dv10*dx10); 
  out[43] += -(1.224744871391589*Ghat[17]*dv10*dx10); 
  out[44] += -(1.224744871391589*Ghat[18]*dv10*dx10); 
  out[45] += -(1.224744871391589*Ghat[19]*dv10*dx10); 
  out[46] += -(0.7071067811865475*Ghat[26]*dv10*dx10); 
  out[47] += -(1.224744871391589*Ghat[20]*dv10*dx10); 
  out[48] += -(1.224744871391589*Ghat[21]*dv10*dx10); 
  out[49] += -(1.224744871391589*Ghat[22]*dv10*dx10); 
  out[50] += -(0.7071067811865475*Ghat[27]*dv10*dx10); 
  out[51] += -(1.224744871391589*Ghat[23]*dv10*dx10); 
  out[52] += -(1.224744871391589*Ghat[24]*dv10*dx10); 
  out[53] += -(0.7071067811865475*Ghat[28]*dv10*dx10); 
  out[54] += -(1.224744871391589*Ghat[25]*dv10*dx10); 
  out[55] += -(0.7071067811865475*Ghat[29]*dv10*dx10); 
  out[56] += -(0.7071067811865475*Ghat[30]*dv10*dx10); 
  out[57] += -(1.224744871391589*Ghat[26]*dv10*dx10); 
  out[58] += -(1.224744871391589*Ghat[27]*dv10*dx10); 
  out[59] += -(1.224744871391589*Ghat[28]*dv10*dx10); 
  out[60] += -(1.224744871391589*Ghat[29]*dv10*dx10); 
  out[61] += -(1.224744871391589*Ghat[30]*dv10*dx10); 
  out[62] += -(0.7071067811865475*Ghat[31]*dv10*dx10); 
  out[63] += -(1.224744871391589*Ghat[31]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[38]+0.4330127018922193*hamil[7]*fedge[21]+0.75*hamil[5]*fedge[17]+0.75*hamil[4]*fedge[13]+0.4330127018922193*hamil[5]*fedge[6]+0.4330127018922193*hamil[4]*fedge[5]+0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fedge[51]+0.4330127018922193*hamil[7]*fedge[39]+0.75*hamil[5]*fedge[32]+0.75*hamil[4]*fedge[26]+0.4330127018922193*hamil[5]*fedge[18]+0.4330127018922193*hamil[4]*fedge[14]+0.75*hamil[1]*fedge[7]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = 0.75*hamil[7]*fedge[52]+0.4330127018922193*hamil[7]*fedge[40]+0.75*hamil[5]*fedge[33]+0.75*hamil[4]*fedge[27]+0.4330127018922193*hamil[5]*fedge[19]+0.4330127018922193*hamil[4]*fedge[15]+0.75*hamil[1]*fedge[8]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = 0.75*hamil[7]*fedge[54]+0.4330127018922193*hamil[7]*fedge[41]+0.75*hamil[5]*fedge[35]+0.75*hamil[4]*fedge[29]+0.4330127018922193*hamil[5]*fedge[20]+0.4330127018922193*hamil[4]*fedge[16]+0.75*hamil[1]*fedge[10]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = 0.75*hamil[5]*fedge[38]+0.4330127018922193*hamil[5]*fedge[21]+0.75*hamil[7]*fedge[17]+0.75*hamil[1]*fedge[13]+0.4330127018922193*fedge[6]*hamil[7]+0.4330127018922193*hamil[1]*fedge[5]+0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]; 
  Ghat[5] = 0.75*hamil[4]*fedge[38]+0.4330127018922193*hamil[4]*fedge[21]+0.75*hamil[1]*fedge[17]+0.75*hamil[7]*fedge[13]+0.4330127018922193*fedge[5]*hamil[7]+0.4330127018922193*hamil[1]*fedge[6]+0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]; 
  Ghat[6] = 0.75*hamil[7]*fedge[59]+0.4330127018922193*hamil[7]*fedge[53]+0.75*hamil[5]*fedge[47]+0.75*hamil[4]*fedge[43]+0.4330127018922193*hamil[5]*fedge[34]+0.4330127018922193*hamil[4]*fedge[28]+0.75*hamil[1]*fedge[22]+0.4330127018922193*hamil[1]*fedge[9]; 
  Ghat[7] = 0.75*hamil[7]*fedge[60]+0.4330127018922193*hamil[7]*fedge[55]+0.75*hamil[5]*fedge[48]+0.75*hamil[4]*fedge[44]+0.4330127018922193*hamil[5]*fedge[36]+0.4330127018922193*hamil[4]*fedge[30]+0.75*hamil[1]*fedge[23]+0.4330127018922193*hamil[1]*fedge[11]; 
  Ghat[8] = 0.75*hamil[7]*fedge[61]+0.4330127018922193*hamil[7]*fedge[56]+0.75*hamil[5]*fedge[49]+0.75*hamil[4]*fedge[45]+0.4330127018922193*hamil[5]*fedge[37]+0.4330127018922193*hamil[4]*fedge[31]+0.75*hamil[1]*fedge[24]+0.4330127018922193*hamil[1]*fedge[12]; 
  Ghat[9] = 0.75*hamil[5]*fedge[51]+0.4330127018922193*hamil[5]*fedge[39]+0.75*hamil[7]*fedge[32]+0.75*hamil[1]*fedge[26]+0.4330127018922193*hamil[7]*fedge[18]+0.4330127018922193*hamil[1]*fedge[14]+0.75*hamil[4]*fedge[7]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[10] = 0.75*hamil[5]*fedge[52]+0.4330127018922193*hamil[5]*fedge[40]+0.75*hamil[7]*fedge[33]+0.75*hamil[1]*fedge[27]+0.4330127018922193*hamil[7]*fedge[19]+0.4330127018922193*hamil[1]*fedge[15]+0.75*hamil[4]*fedge[8]+0.4330127018922193*fedge[3]*hamil[4]; 
  Ghat[11] = 0.75*hamil[5]*fedge[54]+0.4330127018922193*hamil[5]*fedge[41]+0.75*hamil[7]*fedge[35]+0.75*hamil[1]*fedge[29]+0.4330127018922193*hamil[7]*fedge[20]+0.4330127018922193*hamil[1]*fedge[16]+0.75*hamil[4]*fedge[10]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[12] = 0.75*hamil[4]*fedge[51]+0.4330127018922193*hamil[4]*fedge[39]+0.75*hamil[1]*fedge[32]+0.75*hamil[7]*fedge[26]+0.4330127018922193*hamil[1]*fedge[18]+0.4330127018922193*hamil[7]*fedge[14]+0.75*hamil[5]*fedge[7]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[13] = 0.75*hamil[4]*fedge[52]+0.4330127018922193*hamil[4]*fedge[40]+0.75*hamil[1]*fedge[33]+0.75*hamil[7]*fedge[27]+0.4330127018922193*hamil[1]*fedge[19]+0.4330127018922193*hamil[7]*fedge[15]+0.75*hamil[5]*fedge[8]+0.4330127018922193*fedge[3]*hamil[5]; 
  Ghat[14] = 0.75*hamil[4]*fedge[54]+0.4330127018922193*hamil[4]*fedge[41]+0.75*hamil[1]*fedge[35]+0.75*hamil[7]*fedge[29]+0.4330127018922193*hamil[1]*fedge[20]+0.4330127018922193*hamil[7]*fedge[16]+0.75*hamil[5]*fedge[10]+0.4330127018922193*fedge[4]*hamil[5]; 
  Ghat[15] = 0.75*hamil[1]*fedge[38]+0.4330127018922193*hamil[1]*fedge[21]+0.75*hamil[4]*fedge[17]+0.75*hamil[5]*fedge[13]+0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.4330127018922193*hamil[4]*fedge[6]+0.4330127018922193*fedge[5]*hamil[5]; 
  Ghat[16] = 0.75*hamil[7]*fedge[63]+0.4330127018922193*hamil[7]*fedge[62]+0.75*hamil[5]*fedge[58]+0.75*hamil[4]*fedge[57]+0.4330127018922193*hamil[5]*fedge[50]+0.4330127018922193*hamil[4]*fedge[46]+0.75*hamil[1]*fedge[42]+0.4330127018922193*hamil[1]*fedge[25]; 
  Ghat[17] = 0.75*hamil[5]*fedge[59]+0.4330127018922193*hamil[5]*fedge[53]+0.75*hamil[7]*fedge[47]+0.75*hamil[1]*fedge[43]+0.4330127018922193*hamil[7]*fedge[34]+0.4330127018922193*hamil[1]*fedge[28]+0.75*hamil[4]*fedge[22]+0.4330127018922193*hamil[4]*fedge[9]; 
  Ghat[18] = 0.75*hamil[5]*fedge[60]+0.4330127018922193*hamil[5]*fedge[55]+0.75*hamil[7]*fedge[48]+0.75*hamil[1]*fedge[44]+0.4330127018922193*hamil[7]*fedge[36]+0.4330127018922193*hamil[1]*fedge[30]+0.75*hamil[4]*fedge[23]+0.4330127018922193*hamil[4]*fedge[11]; 
  Ghat[19] = 0.75*hamil[5]*fedge[61]+0.4330127018922193*hamil[5]*fedge[56]+0.75*hamil[7]*fedge[49]+0.75*hamil[1]*fedge[45]+0.4330127018922193*hamil[7]*fedge[37]+0.4330127018922193*hamil[1]*fedge[31]+0.75*hamil[4]*fedge[24]+0.4330127018922193*hamil[4]*fedge[12]; 
  Ghat[20] = 0.75*hamil[4]*fedge[59]+0.4330127018922193*hamil[4]*fedge[53]+0.75*hamil[1]*fedge[47]+0.75*hamil[7]*fedge[43]+0.4330127018922193*hamil[1]*fedge[34]+0.4330127018922193*hamil[7]*fedge[28]+0.75*hamil[5]*fedge[22]+0.4330127018922193*hamil[5]*fedge[9]; 
  Ghat[21] = 0.75*hamil[4]*fedge[60]+0.4330127018922193*hamil[4]*fedge[55]+0.75*hamil[1]*fedge[48]+0.75*hamil[7]*fedge[44]+0.4330127018922193*hamil[1]*fedge[36]+0.4330127018922193*hamil[7]*fedge[30]+0.75*hamil[5]*fedge[23]+0.4330127018922193*hamil[5]*fedge[11]; 
  Ghat[22] = 0.75*hamil[4]*fedge[61]+0.4330127018922193*hamil[4]*fedge[56]+0.75*hamil[1]*fedge[49]+0.75*hamil[7]*fedge[45]+0.4330127018922193*hamil[1]*fedge[37]+0.4330127018922193*hamil[7]*fedge[31]+0.75*hamil[5]*fedge[24]+0.4330127018922193*hamil[5]*fedge[12]; 
  Ghat[23] = 0.75*hamil[1]*fedge[51]+0.4330127018922193*hamil[1]*fedge[39]+0.75*hamil[4]*fedge[32]+0.75*hamil[5]*fedge[26]+0.4330127018922193*hamil[4]*fedge[18]+0.4330127018922193*hamil[5]*fedge[14]+0.75*fedge[7]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]; 
  Ghat[24] = 0.75*hamil[1]*fedge[52]+0.4330127018922193*hamil[1]*fedge[40]+0.75*hamil[4]*fedge[33]+0.75*hamil[5]*fedge[27]+0.4330127018922193*hamil[4]*fedge[19]+0.4330127018922193*hamil[5]*fedge[15]+0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[3]*hamil[7]; 
  Ghat[25] = 0.75*hamil[1]*fedge[54]+0.4330127018922193*hamil[1]*fedge[41]+0.75*hamil[4]*fedge[35]+0.75*hamil[5]*fedge[29]+0.4330127018922193*hamil[4]*fedge[20]+0.4330127018922193*hamil[5]*fedge[16]+0.75*hamil[7]*fedge[10]+0.4330127018922193*fedge[4]*hamil[7]; 
  Ghat[26] = 0.75*hamil[5]*fedge[63]+0.4330127018922193*hamil[5]*fedge[62]+0.75*hamil[7]*fedge[58]+0.75*hamil[1]*fedge[57]+0.4330127018922193*hamil[7]*fedge[50]+0.4330127018922193*hamil[1]*fedge[46]+0.75*hamil[4]*fedge[42]+0.4330127018922193*hamil[4]*fedge[25]; 
  Ghat[27] = 0.75*hamil[4]*fedge[63]+0.4330127018922193*hamil[4]*fedge[62]+0.75*hamil[1]*fedge[58]+0.75*hamil[7]*fedge[57]+0.4330127018922193*hamil[1]*fedge[50]+0.4330127018922193*hamil[7]*fedge[46]+0.75*hamil[5]*fedge[42]+0.4330127018922193*hamil[5]*fedge[25]; 
  Ghat[28] = 0.75*hamil[1]*fedge[59]+0.4330127018922193*hamil[1]*fedge[53]+0.75*hamil[4]*fedge[47]+0.75*hamil[5]*fedge[43]+0.4330127018922193*hamil[4]*fedge[34]+0.4330127018922193*hamil[5]*fedge[28]+0.75*hamil[7]*fedge[22]+0.4330127018922193*hamil[7]*fedge[9]; 
  Ghat[29] = 0.75*hamil[1]*fedge[60]+0.4330127018922193*hamil[1]*fedge[55]+0.75*hamil[4]*fedge[48]+0.75*hamil[5]*fedge[44]+0.4330127018922193*hamil[4]*fedge[36]+0.4330127018922193*hamil[5]*fedge[30]+0.75*hamil[7]*fedge[23]+0.4330127018922193*hamil[7]*fedge[11]; 
  Ghat[30] = 0.75*hamil[1]*fedge[61]+0.4330127018922193*hamil[1]*fedge[56]+0.75*hamil[4]*fedge[49]+0.75*hamil[5]*fedge[45]+0.4330127018922193*hamil[4]*fedge[37]+0.4330127018922193*hamil[5]*fedge[31]+0.75*hamil[7]*fedge[24]+0.4330127018922193*hamil[7]*fedge[12]; 
  Ghat[31] = 0.75*hamil[1]*fedge[63]+0.4330127018922193*hamil[1]*fedge[62]+0.75*hamil[4]*fedge[58]+0.75*hamil[5]*fedge[57]+0.4330127018922193*hamil[4]*fedge[50]+0.4330127018922193*hamil[5]*fedge[46]+0.75*hamil[7]*fedge[42]+0.4330127018922193*hamil[7]*fedge[25]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[38])+0.4330127018922193*hamil[7]*fskin[21]-0.75*hamil[5]*fskin[17]-0.75*hamil[4]*fskin[13]+0.4330127018922193*hamil[5]*fskin[6]+0.4330127018922193*hamil[4]*fskin[5]-0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[51])+0.4330127018922193*hamil[7]*fskin[39]-0.75*hamil[5]*fskin[32]-0.75*hamil[4]*fskin[26]+0.4330127018922193*hamil[5]*fskin[18]+0.4330127018922193*hamil[4]*fskin[14]-0.75*hamil[1]*fskin[7]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = -(0.75*hamil[7]*fskin[52])+0.4330127018922193*hamil[7]*fskin[40]-0.75*hamil[5]*fskin[33]-0.75*hamil[4]*fskin[27]+0.4330127018922193*hamil[5]*fskin[19]+0.4330127018922193*hamil[4]*fskin[15]-0.75*hamil[1]*fskin[8]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = -(0.75*hamil[7]*fskin[54])+0.4330127018922193*hamil[7]*fskin[41]-0.75*hamil[5]*fskin[35]-0.75*hamil[4]*fskin[29]+0.4330127018922193*hamil[5]*fskin[20]+0.4330127018922193*hamil[4]*fskin[16]-0.75*hamil[1]*fskin[10]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[5]*fskin[38])+0.4330127018922193*hamil[5]*fskin[21]-0.75*hamil[7]*fskin[17]-0.75*hamil[1]*fskin[13]+0.4330127018922193*fskin[6]*hamil[7]+0.4330127018922193*hamil[1]*fskin[5]-0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]; 
  Ghat[5] = -(0.75*hamil[4]*fskin[38])+0.4330127018922193*hamil[4]*fskin[21]-0.75*hamil[1]*fskin[17]-0.75*hamil[7]*fskin[13]+0.4330127018922193*fskin[5]*hamil[7]+0.4330127018922193*hamil[1]*fskin[6]-0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]; 
  Ghat[6] = -(0.75*hamil[7]*fskin[59])+0.4330127018922193*hamil[7]*fskin[53]-0.75*hamil[5]*fskin[47]-0.75*hamil[4]*fskin[43]+0.4330127018922193*hamil[5]*fskin[34]+0.4330127018922193*hamil[4]*fskin[28]-0.75*hamil[1]*fskin[22]+0.4330127018922193*hamil[1]*fskin[9]; 
  Ghat[7] = -(0.75*hamil[7]*fskin[60])+0.4330127018922193*hamil[7]*fskin[55]-0.75*hamil[5]*fskin[48]-0.75*hamil[4]*fskin[44]+0.4330127018922193*hamil[5]*fskin[36]+0.4330127018922193*hamil[4]*fskin[30]-0.75*hamil[1]*fskin[23]+0.4330127018922193*hamil[1]*fskin[11]; 
  Ghat[8] = -(0.75*hamil[7]*fskin[61])+0.4330127018922193*hamil[7]*fskin[56]-0.75*hamil[5]*fskin[49]-0.75*hamil[4]*fskin[45]+0.4330127018922193*hamil[5]*fskin[37]+0.4330127018922193*hamil[4]*fskin[31]-0.75*hamil[1]*fskin[24]+0.4330127018922193*hamil[1]*fskin[12]; 
  Ghat[9] = -(0.75*hamil[5]*fskin[51])+0.4330127018922193*hamil[5]*fskin[39]-0.75*hamil[7]*fskin[32]-0.75*hamil[1]*fskin[26]+0.4330127018922193*hamil[7]*fskin[18]+0.4330127018922193*hamil[1]*fskin[14]-0.75*hamil[4]*fskin[7]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[10] = -(0.75*hamil[5]*fskin[52])+0.4330127018922193*hamil[5]*fskin[40]-0.75*hamil[7]*fskin[33]-0.75*hamil[1]*fskin[27]+0.4330127018922193*hamil[7]*fskin[19]+0.4330127018922193*hamil[1]*fskin[15]-0.75*hamil[4]*fskin[8]+0.4330127018922193*fskin[3]*hamil[4]; 
  Ghat[11] = -(0.75*hamil[5]*fskin[54])+0.4330127018922193*hamil[5]*fskin[41]-0.75*hamil[7]*fskin[35]-0.75*hamil[1]*fskin[29]+0.4330127018922193*hamil[7]*fskin[20]+0.4330127018922193*hamil[1]*fskin[16]-0.75*hamil[4]*fskin[10]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[12] = -(0.75*hamil[4]*fskin[51])+0.4330127018922193*hamil[4]*fskin[39]-0.75*hamil[1]*fskin[32]-0.75*hamil[7]*fskin[26]+0.4330127018922193*hamil[1]*fskin[18]+0.4330127018922193*hamil[7]*fskin[14]-0.75*hamil[5]*fskin[7]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[13] = -(0.75*hamil[4]*fskin[52])+0.4330127018922193*hamil[4]*fskin[40]-0.75*hamil[1]*fskin[33]-0.75*hamil[7]*fskin[27]+0.4330127018922193*hamil[1]*fskin[19]+0.4330127018922193*hamil[7]*fskin[15]-0.75*hamil[5]*fskin[8]+0.4330127018922193*fskin[3]*hamil[5]; 
  Ghat[14] = -(0.75*hamil[4]*fskin[54])+0.4330127018922193*hamil[4]*fskin[41]-0.75*hamil[1]*fskin[35]-0.75*hamil[7]*fskin[29]+0.4330127018922193*hamil[1]*fskin[20]+0.4330127018922193*hamil[7]*fskin[16]-0.75*hamil[5]*fskin[10]+0.4330127018922193*fskin[4]*hamil[5]; 
  Ghat[15] = -(0.75*hamil[1]*fskin[38])+0.4330127018922193*hamil[1]*fskin[21]-0.75*hamil[4]*fskin[17]-0.75*hamil[5]*fskin[13]-0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.4330127018922193*hamil[4]*fskin[6]+0.4330127018922193*fskin[5]*hamil[5]; 
  Ghat[16] = -(0.75*hamil[7]*fskin[63])+0.4330127018922193*hamil[7]*fskin[62]-0.75*hamil[5]*fskin[58]-0.75*hamil[4]*fskin[57]+0.4330127018922193*hamil[5]*fskin[50]+0.4330127018922193*hamil[4]*fskin[46]-0.75*hamil[1]*fskin[42]+0.4330127018922193*hamil[1]*fskin[25]; 
  Ghat[17] = -(0.75*hamil[5]*fskin[59])+0.4330127018922193*hamil[5]*fskin[53]-0.75*hamil[7]*fskin[47]-0.75*hamil[1]*fskin[43]+0.4330127018922193*hamil[7]*fskin[34]+0.4330127018922193*hamil[1]*fskin[28]-0.75*hamil[4]*fskin[22]+0.4330127018922193*hamil[4]*fskin[9]; 
  Ghat[18] = -(0.75*hamil[5]*fskin[60])+0.4330127018922193*hamil[5]*fskin[55]-0.75*hamil[7]*fskin[48]-0.75*hamil[1]*fskin[44]+0.4330127018922193*hamil[7]*fskin[36]+0.4330127018922193*hamil[1]*fskin[30]-0.75*hamil[4]*fskin[23]+0.4330127018922193*hamil[4]*fskin[11]; 
  Ghat[19] = -(0.75*hamil[5]*fskin[61])+0.4330127018922193*hamil[5]*fskin[56]-0.75*hamil[7]*fskin[49]-0.75*hamil[1]*fskin[45]+0.4330127018922193*hamil[7]*fskin[37]+0.4330127018922193*hamil[1]*fskin[31]-0.75*hamil[4]*fskin[24]+0.4330127018922193*hamil[4]*fskin[12]; 
  Ghat[20] = -(0.75*hamil[4]*fskin[59])+0.4330127018922193*hamil[4]*fskin[53]-0.75*hamil[1]*fskin[47]-0.75*hamil[7]*fskin[43]+0.4330127018922193*hamil[1]*fskin[34]+0.4330127018922193*hamil[7]*fskin[28]-0.75*hamil[5]*fskin[22]+0.4330127018922193*hamil[5]*fskin[9]; 
  Ghat[21] = -(0.75*hamil[4]*fskin[60])+0.4330127018922193*hamil[4]*fskin[55]-0.75*hamil[1]*fskin[48]-0.75*hamil[7]*fskin[44]+0.4330127018922193*hamil[1]*fskin[36]+0.4330127018922193*hamil[7]*fskin[30]-0.75*hamil[5]*fskin[23]+0.4330127018922193*hamil[5]*fskin[11]; 
  Ghat[22] = -(0.75*hamil[4]*fskin[61])+0.4330127018922193*hamil[4]*fskin[56]-0.75*hamil[1]*fskin[49]-0.75*hamil[7]*fskin[45]+0.4330127018922193*hamil[1]*fskin[37]+0.4330127018922193*hamil[7]*fskin[31]-0.75*hamil[5]*fskin[24]+0.4330127018922193*hamil[5]*fskin[12]; 
  Ghat[23] = -(0.75*hamil[1]*fskin[51])+0.4330127018922193*hamil[1]*fskin[39]-0.75*hamil[4]*fskin[32]-0.75*hamil[5]*fskin[26]+0.4330127018922193*hamil[4]*fskin[18]+0.4330127018922193*hamil[5]*fskin[14]-0.75*fskin[7]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]; 
  Ghat[24] = -(0.75*hamil[1]*fskin[52])+0.4330127018922193*hamil[1]*fskin[40]-0.75*hamil[4]*fskin[33]-0.75*hamil[5]*fskin[27]+0.4330127018922193*hamil[4]*fskin[19]+0.4330127018922193*hamil[5]*fskin[15]-0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[3]*hamil[7]; 
  Ghat[25] = -(0.75*hamil[1]*fskin[54])+0.4330127018922193*hamil[1]*fskin[41]-0.75*hamil[4]*fskin[35]-0.75*hamil[5]*fskin[29]+0.4330127018922193*hamil[4]*fskin[20]+0.4330127018922193*hamil[5]*fskin[16]-0.75*hamil[7]*fskin[10]+0.4330127018922193*fskin[4]*hamil[7]; 
  Ghat[26] = -(0.75*hamil[5]*fskin[63])+0.4330127018922193*hamil[5]*fskin[62]-0.75*hamil[7]*fskin[58]-0.75*hamil[1]*fskin[57]+0.4330127018922193*hamil[7]*fskin[50]+0.4330127018922193*hamil[1]*fskin[46]-0.75*hamil[4]*fskin[42]+0.4330127018922193*hamil[4]*fskin[25]; 
  Ghat[27] = -(0.75*hamil[4]*fskin[63])+0.4330127018922193*hamil[4]*fskin[62]-0.75*hamil[1]*fskin[58]-0.75*hamil[7]*fskin[57]+0.4330127018922193*hamil[1]*fskin[50]+0.4330127018922193*hamil[7]*fskin[46]-0.75*hamil[5]*fskin[42]+0.4330127018922193*hamil[5]*fskin[25]; 
  Ghat[28] = -(0.75*hamil[1]*fskin[59])+0.4330127018922193*hamil[1]*fskin[53]-0.75*hamil[4]*fskin[47]-0.75*hamil[5]*fskin[43]+0.4330127018922193*hamil[4]*fskin[34]+0.4330127018922193*hamil[5]*fskin[28]-0.75*hamil[7]*fskin[22]+0.4330127018922193*hamil[7]*fskin[9]; 
  Ghat[29] = -(0.75*hamil[1]*fskin[60])+0.4330127018922193*hamil[1]*fskin[55]-0.75*hamil[4]*fskin[48]-0.75*hamil[5]*fskin[44]+0.4330127018922193*hamil[4]*fskin[36]+0.4330127018922193*hamil[5]*fskin[30]-0.75*hamil[7]*fskin[23]+0.4330127018922193*hamil[7]*fskin[11]; 
  Ghat[30] = -(0.75*hamil[1]*fskin[61])+0.4330127018922193*hamil[1]*fskin[56]-0.75*hamil[4]*fskin[49]-0.75*hamil[5]*fskin[45]+0.4330127018922193*hamil[4]*fskin[37]+0.4330127018922193*hamil[5]*fskin[31]-0.75*hamil[7]*fskin[24]+0.4330127018922193*hamil[7]*fskin[12]; 
  Ghat[31] = -(0.75*hamil[1]*fskin[63])+0.4330127018922193*hamil[1]*fskin[62]-0.75*hamil[4]*fskin[58]-0.75*hamil[5]*fskin[57]+0.4330127018922193*hamil[4]*fskin[50]+0.4330127018922193*hamil[5]*fskin[46]-0.75*hamil[7]*fskin[42]+0.4330127018922193*hamil[7]*fskin[25]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[7] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[8] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[9] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[10] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[11] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10*dx10; 
  out[13] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[14] += 0.7071067811865475*Ghat[9]*dv10*dx10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10*dx10; 
  out[16] += 0.7071067811865475*Ghat[11]*dv10*dx10; 
  out[17] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[18] += 0.7071067811865475*Ghat[12]*dv10*dx10; 
  out[19] += 0.7071067811865475*Ghat[13]*dv10*dx10; 
  out[20] += 0.7071067811865475*Ghat[14]*dv10*dx10; 
  out[21] += 0.7071067811865475*Ghat[15]*dv10*dx10; 
  out[22] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[23] += -(1.224744871391589*Ghat[7]*dv10*dx10); 
  out[24] += -(1.224744871391589*Ghat[8]*dv10*dx10); 
  out[25] += 0.7071067811865475*Ghat[16]*dv10*dx10; 
  out[26] += -(1.224744871391589*Ghat[9]*dv10*dx10); 
  out[27] += -(1.224744871391589*Ghat[10]*dv10*dx10); 
  out[28] += 0.7071067811865475*Ghat[17]*dv10*dx10; 
  out[29] += -(1.224744871391589*Ghat[11]*dv10*dx10); 
  out[30] += 0.7071067811865475*Ghat[18]*dv10*dx10; 
  out[31] += 0.7071067811865475*Ghat[19]*dv10*dx10; 
  out[32] += -(1.224744871391589*Ghat[12]*dv10*dx10); 
  out[33] += -(1.224744871391589*Ghat[13]*dv10*dx10); 
  out[34] += 0.7071067811865475*Ghat[20]*dv10*dx10; 
  out[35] += -(1.224744871391589*Ghat[14]*dv10*dx10); 
  out[36] += 0.7071067811865475*Ghat[21]*dv10*dx10; 
  out[37] += 0.7071067811865475*Ghat[22]*dv10*dx10; 
  out[38] += -(1.224744871391589*Ghat[15]*dv10*dx10); 
  out[39] += 0.7071067811865475*Ghat[23]*dv10*dx10; 
  out[40] += 0.7071067811865475*Ghat[24]*dv10*dx10; 
  out[41] += 0.7071067811865475*Ghat[25]*dv10*dx10; 
  out[42] += -(1.224744871391589*Ghat[16]*dv10*dx10); 
  out[43] += -(1.224744871391589*Ghat[17]*dv10*dx10); 
  out[44] += -(1.224744871391589*Ghat[18]*dv10*dx10); 
  out[45] += -(1.224744871391589*Ghat[19]*dv10*dx10); 
  out[46] += 0.7071067811865475*Ghat[26]*dv10*dx10; 
  out[47] += -(1.224744871391589*Ghat[20]*dv10*dx10); 
  out[48] += -(1.224744871391589*Ghat[21]*dv10*dx10); 
  out[49] += -(1.224744871391589*Ghat[22]*dv10*dx10); 
  out[50] += 0.7071067811865475*Ghat[27]*dv10*dx10; 
  out[51] += -(1.224744871391589*Ghat[23]*dv10*dx10); 
  out[52] += -(1.224744871391589*Ghat[24]*dv10*dx10); 
  out[53] += 0.7071067811865475*Ghat[28]*dv10*dx10; 
  out[54] += -(1.224744871391589*Ghat[25]*dv10*dx10); 
  out[55] += 0.7071067811865475*Ghat[29]*dv10*dx10; 
  out[56] += 0.7071067811865475*Ghat[30]*dv10*dx10; 
  out[57] += -(1.224744871391589*Ghat[26]*dv10*dx10); 
  out[58] += -(1.224744871391589*Ghat[27]*dv10*dx10); 
  out[59] += -(1.224744871391589*Ghat[28]*dv10*dx10); 
  out[60] += -(1.224744871391589*Ghat[29]*dv10*dx10); 
  out[61] += -(1.224744871391589*Ghat[30]*dv10*dx10); 
  out[62] += 0.7071067811865475*Ghat[31]*dv10*dx10; 
  out[63] += -(1.224744871391589*Ghat[31]*dv10*dx10); 

  } 
  return 0.0;

} 
