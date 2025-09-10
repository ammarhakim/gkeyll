#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[3]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.732050807568877*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[64] = {0.0}; 
  alpha_cdim[0] = 2.0*p0_over_gamma[0]*dx10; 
  alpha_cdim[3] = 2.0*p0_over_gamma[1]*dx10; 
  alpha_cdim[4] = 2.0*p0_over_gamma[2]*dx10; 
  alpha_cdim[10] = 2.0*p0_over_gamma[3]*dx10; 
  alpha_cdim[24] = 2.0*p0_over_gamma[5]*dx10; 
  cflFreq_mid += fabs(0.125*alpha_cdim[0]-0.1397542485937369*alpha_cdim[24]); 

  alpha_cdim[32] = 2.0*p1_over_gamma[0]*dx11; 
  alpha_cdim[35] = 2.0*p1_over_gamma[1]*dx11; 
  alpha_cdim[36] = 2.0*p1_over_gamma[2]*dx11; 
  alpha_cdim[42] = 2.0*p1_over_gamma[3]*dx11; 
  alpha_cdim[48] = 2.0*p1_over_gamma[4]*dx11; 
  cflFreq_mid += fabs(0.125*alpha_cdim[32]-0.1397542485937369*alpha_cdim[48]); 

  out[1] += 0.4330127018922193*(alpha_cdim[24]*f[24]+alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[16]*alpha_cdim[48]+f[10]*alpha_cdim[42]+f[4]*alpha_cdim[36]+f[3]*alpha_cdim[35]+f[0]*alpha_cdim[32]); 
  out[5] += 0.4330127018922193*(f[17]*alpha_cdim[48]+f[13]*alpha_cdim[42]+f[8]*alpha_cdim[36]+f[6]*alpha_cdim[35]+f[1]*alpha_cdim[32]+alpha_cdim[24]*f[26]+alpha_cdim[10]*f[14]+alpha_cdim[4]*f[9]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.4330127018922193*alpha_cdim[24]*f[27]+0.3872983346207416*(alpha_cdim[10]*f[19]+alpha_cdim[3]*f[16])+0.4330127018922193*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[7] += 0.3872983346207416*(f[3]*alpha_cdim[48]+f[19]*alpha_cdim[42])+0.4330127018922193*(f[4]*alpha_cdim[42]+f[10]*alpha_cdim[36])+0.3872983346207416*f[16]*alpha_cdim[35]+0.4330127018922193*(f[0]*alpha_cdim[35]+f[3]*alpha_cdim[32]); 
  out[8] += 0.3872983346207416*(alpha_cdim[10]*f[27]+alpha_cdim[4]*f[24]+f[4]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[9] += 0.4330127018922193*f[19]*alpha_cdim[48]+(0.3872983346207416*f[27]+0.4330127018922193*f[3])*alpha_cdim[42]+0.3872983346207416*f[24]*alpha_cdim[36]+0.4330127018922193*(f[0]*alpha_cdim[36]+f[10]*alpha_cdim[35]+f[4]*alpha_cdim[32]); 
  out[11] += 0.3872983346207416*(f[6]*alpha_cdim[48]+f[21]*alpha_cdim[42])+0.4330127018922193*(f[8]*alpha_cdim[42]+f[13]*alpha_cdim[36])+0.3872983346207416*f[17]*alpha_cdim[35]+0.4330127018922193*(f[1]*alpha_cdim[35]+f[6]*alpha_cdim[32]+alpha_cdim[24]*f[30])+0.3872983346207416*(alpha_cdim[10]*f[22]+alpha_cdim[3]*f[18])+0.4330127018922193*(alpha_cdim[4]*f[14]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[2]*alpha_cdim[3]); 
  out[12] += 0.4330127018922193*f[21]*alpha_cdim[48]+(0.3872983346207416*f[29]+0.4330127018922193*f[6])*alpha_cdim[42]+0.3872983346207416*f[25]*alpha_cdim[36]+0.4330127018922193*(f[1]*alpha_cdim[36]+f[13]*alpha_cdim[35]+f[8]*alpha_cdim[32])+0.3872983346207416*(alpha_cdim[10]*f[30]+alpha_cdim[4]*f[26]+f[9]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[3]*f[14]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[2]*alpha_cdim[4]); 
  out[13] += 0.3872983346207416*(alpha_cdim[4]*f[27]+alpha_cdim[10]*f[24]+f[10]*alpha_cdim[24]+alpha_cdim[3]*f[19]+alpha_cdim[10]*f[16])+0.4330127018922193*(alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[14] += 0.3872983346207416*f[10]*alpha_cdim[48]+(0.3872983346207416*(f[24]+f[16])+0.4330127018922193*f[0])*alpha_cdim[42]+(0.3872983346207416*f[27]+0.4330127018922193*f[3])*alpha_cdim[36]+0.3872983346207416*f[19]*alpha_cdim[35]+0.4330127018922193*(f[4]*alpha_cdim[35]+f[10]*alpha_cdim[32]); 
  out[15] += 0.3872983346207416*f[13]*alpha_cdim[48]+(0.3872983346207416*(f[25]+f[17])+0.4330127018922193*f[1])*alpha_cdim[42]+(0.3872983346207416*f[29]+0.4330127018922193*f[6])*alpha_cdim[36]+0.3872983346207416*f[21]*alpha_cdim[35]+0.4330127018922193*(f[8]*alpha_cdim[35]+f[13]*alpha_cdim[32])+0.3872983346207416*(alpha_cdim[4]*f[30]+alpha_cdim[10]*f[26]+f[14]*alpha_cdim[24]+alpha_cdim[3]*f[22]+alpha_cdim[10]*f[18])+0.4330127018922193*(alpha_cdim[0]*f[14]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+alpha_cdim[4]*f[7]); 
  out[17] += 0.4330127018922193*(alpha_cdim[4]*f[19]+alpha_cdim[0]*f[16])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[3]*f[3]); 
  out[18] += (0.276641667586244*f[16]+0.4330127018922193*f[0])*alpha_cdim[48]+0.3872983346207416*f[10]*alpha_cdim[42]+0.4330127018922193*f[19]*alpha_cdim[36]+0.3872983346207416*f[3]*alpha_cdim[35]+0.4330127018922193*f[16]*alpha_cdim[32]; 
  out[20] += (0.276641667586244*f[17]+0.4330127018922193*f[1])*alpha_cdim[48]+0.3872983346207416*f[13]*alpha_cdim[42]+0.4330127018922193*f[21]*alpha_cdim[36]+0.3872983346207416*f[6]*alpha_cdim[35]+0.4330127018922193*(f[17]*alpha_cdim[32]+alpha_cdim[4]*f[22]+alpha_cdim[0]*f[18])+0.3872983346207416*(alpha_cdim[10]*f[14]+alpha_cdim[3]*f[7]); 
  out[21] += 0.3464101615137755*alpha_cdim[10]*f[27]+0.3872983346207416*f[19]*alpha_cdim[24]+0.4330127018922193*(alpha_cdim[0]*f[19]+alpha_cdim[4]*f[16])+0.3872983346207416*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]); 
  out[22] += (0.276641667586244*f[19]+0.4330127018922193*f[4])*alpha_cdim[48]+(0.3464101615137755*f[27]+0.3872983346207416*f[3])*alpha_cdim[42]+0.4330127018922193*f[16]*alpha_cdim[36]+0.3872983346207416*f[10]*alpha_cdim[35]+0.4330127018922193*f[19]*alpha_cdim[32]; 
  out[23] += (0.276641667586244*f[21]+0.4330127018922193*f[8])*alpha_cdim[48]+(0.3464101615137755*f[29]+0.3872983346207416*f[6])*alpha_cdim[42]+0.4330127018922193*f[17]*alpha_cdim[36]+0.3872983346207416*f[13]*alpha_cdim[35]+0.4330127018922193*f[21]*alpha_cdim[32]+0.3464101615137755*alpha_cdim[10]*f[30]+0.3872983346207416*f[22]*alpha_cdim[24]+0.4330127018922193*(alpha_cdim[0]*f[22]+alpha_cdim[4]*f[18])+0.3872983346207416*(alpha_cdim[3]*f[14]+f[7]*alpha_cdim[10]); 
  out[25] += 0.4330127018922193*alpha_cdim[3]*f[27]+0.276641667586244*alpha_cdim[24]*f[24]+0.4330127018922193*(alpha_cdim[0]*f[24]+f[0]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]); 
  out[26] += 0.3872983346207416*(f[10]*alpha_cdim[42]+f[4]*alpha_cdim[36])+0.4330127018922193*(f[27]*alpha_cdim[35]+f[24]*alpha_cdim[32]); 
  out[28] += 0.3872983346207416*(f[13]*alpha_cdim[42]+f[8]*alpha_cdim[36])+0.4330127018922193*(f[29]*alpha_cdim[35]+f[25]*alpha_cdim[32]+alpha_cdim[3]*f[30])+0.276641667586244*alpha_cdim[24]*f[26]+0.4330127018922193*(alpha_cdim[0]*f[26]+f[2]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[10]*f[14]+alpha_cdim[4]*f[9]); 
  out[29] += 0.276641667586244*alpha_cdim[24]*f[27]+0.4330127018922193*(alpha_cdim[0]*f[27]+alpha_cdim[3]*f[24]+f[3]*alpha_cdim[24])+0.3464101615137755*alpha_cdim[10]*f[19]+0.3872983346207416*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]); 
  out[30] += 0.3872983346207416*f[27]*alpha_cdim[48]+0.3464101615137755*f[19]*alpha_cdim[42]+0.3872983346207416*(f[4]*alpha_cdim[42]+f[10]*alpha_cdim[36])+0.4330127018922193*(f[24]*alpha_cdim[35]+f[27]*alpha_cdim[32]); 
  out[31] += 0.3872983346207416*f[29]*alpha_cdim[48]+0.3464101615137755*f[21]*alpha_cdim[42]+0.3872983346207416*(f[8]*alpha_cdim[42]+f[13]*alpha_cdim[36])+0.4330127018922193*(f[25]*alpha_cdim[35]+f[29]*alpha_cdim[32])+0.276641667586244*alpha_cdim[24]*f[30]+0.4330127018922193*(alpha_cdim[0]*f[30]+alpha_cdim[3]*f[26]+f[7]*alpha_cdim[24])+0.3464101615137755*alpha_cdim[10]*f[22]+0.3872983346207416*(alpha_cdim[4]*f[14]+f[9]*alpha_cdim[10]); 

  return 3.0*cflFreq_mid; 
} 
