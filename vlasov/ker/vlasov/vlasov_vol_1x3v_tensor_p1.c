#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_1x3v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:      q/m*EM fields.
  // cot_vec:   Only used in gen geo.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &field[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &field[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *E2 = &field[4]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &field[6]; 
  const double *B1 = &field[8]; 
  const double *B2 = &field[10]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[16] = {0.0}; 
  double alpha_vdim[48] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[3] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[4] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_vdim[6] = 0.8164965809277261*B2[1]*dv10*dv2; 
  alpha_vdim[8] = -0.8164965809277261*B1[1]*dv10*dv3; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[16] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[17] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[18] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[20] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_vdim[21] = -0.8164965809277261*B2[1]*dv1*dv11; 
  alpha_vdim[24] = 0.8164965809277261*B0[1]*dv11*dv3; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[16]); 

  alpha_vdim[32] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[33] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[34] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[35] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_vdim[37] = 0.8164965809277261*B1[1]*dv1*dv12; 
  alpha_vdim[38] = -0.8164965809277261*B0[1]*dv12*dv2; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[32]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[8]*alpha_vdim[24]+f[5]*alpha_vdim[21]+f[4]*alpha_vdim[20]+f[2]*alpha_vdim[18]+f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[4] += 0.4330127018922193*(f[6]*alpha_vdim[38]+f[5]*alpha_vdim[37]+f[3]*alpha_vdim[35]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5] += 0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(f[4]*alpha_vdim[24]+f[2]*alpha_vdim[21]+f[8]*alpha_vdim[20]+f[5]*alpha_vdim[18]+f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[12]*alpha_vdim[24]+f[1]*alpha_vdim[21]+f[9]*alpha_vdim[20]+f[0]*alpha_vdim[18]+f[5]*alpha_vdim[17]+f[2]*alpha_vdim[16]+alpha_vdim[8]*f[13]+alpha_vdim[4]*f[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.4330127018922193*(f[3]*alpha_vdim[38]+f[2]*alpha_vdim[37]+f[6]*alpha_vdim[35]+f[5]*alpha_vdim[34]+f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[11]*alpha_vdim[38]+f[1]*alpha_vdim[37]+f[7]*alpha_vdim[35]+f[0]*alpha_vdim[34]+f[5]*alpha_vdim[33]+f[2]*alpha_vdim[32]+alpha_vdim[6]*f[13]+alpha_vdim[3]*f[10]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.4330127018922193*(f[1]*alpha_vdim[38]+f[11]*alpha_vdim[37]+f[0]*alpha_vdim[35]+f[7]*alpha_vdim[34]+f[6]*alpha_vdim[33]+f[3]*alpha_vdim[32]+f[1]*alpha_vdim[24]+f[12]*alpha_vdim[21]+f[0]*alpha_vdim[20]+f[9]*alpha_vdim[18]+f[8]*alpha_vdim[17]+f[4]*alpha_vdim[16]); 
  out[11] += 0.4330127018922193*(f[9]*alpha_vdim[24]+f[0]*alpha_vdim[21]+f[12]*alpha_vdim[20]+f[1]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[5]*alpha_vdim[16]+alpha_vdim[4]*f[13]+alpha_vdim[8]*f[10]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[12] += 0.4330127018922193*(f[7]*alpha_vdim[38]+f[0]*alpha_vdim[37]+f[11]*alpha_vdim[35]+f[1]*alpha_vdim[34]+f[2]*alpha_vdim[33]+f[5]*alpha_vdim[32]+alpha_vdim[3]*f[13]+alpha_vdim[6]*f[10]+alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[13] += 0.4330127018922193*(f[0]*alpha_vdim[38]+f[7]*alpha_vdim[37]+f[1]*alpha_vdim[35]+f[11]*alpha_vdim[34]+f[3]*alpha_vdim[33]+f[6]*alpha_vdim[32]+f[0]*alpha_vdim[24]+f[9]*alpha_vdim[21]+f[1]*alpha_vdim[20]+f[12]*alpha_vdim[18]+f[4]*alpha_vdim[17]+f[8]*alpha_vdim[16]+alpha_cdim[2]*f[14]+alpha_cdim[0]*f[10]); 
  out[14] += 0.4330127018922193*(f[5]*alpha_vdim[38]+f[6]*alpha_vdim[37]+f[2]*alpha_vdim[35]+f[3]*alpha_vdim[34]+f[11]*alpha_vdim[33]+f[7]*alpha_vdim[32]+f[5]*alpha_vdim[24]+f[8]*alpha_vdim[21]+f[2]*alpha_vdim[20]+f[4]*alpha_vdim[18]+f[12]*alpha_vdim[17]+f[9]*alpha_vdim[16]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[10]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[15] += 0.4330127018922193*(f[2]*alpha_vdim[38]+f[3]*alpha_vdim[37]+f[5]*alpha_vdim[35]+f[6]*alpha_vdim[34]+f[7]*alpha_vdim[33]+f[11]*alpha_vdim[32]+f[2]*alpha_vdim[24]+f[4]*alpha_vdim[21]+f[5]*alpha_vdim[20]+f[8]*alpha_vdim[18]+f[9]*alpha_vdim[17]+f[12]*alpha_vdim[16]+alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+(alpha_cdim[2]+alpha_vdim[1])*f[10]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]); 

  return cflFreq_mid; 
} 
