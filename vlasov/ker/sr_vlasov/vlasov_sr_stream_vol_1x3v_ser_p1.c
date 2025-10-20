#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.732050807568877*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.732050807568877*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.732050807568877*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.732050807568877*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.732050807568877*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.732050807568877*gamma[19]*dv10; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[40] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[4] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[7] = 1.414213562373095*p0_over_gamma[4]*dx10; 
  alpha_cdim[9] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  alpha_cdim[10] = 1.414213562373095*p0_over_gamma[6]*dx10; 
  alpha_cdim[14] = 1.414213562373095*p0_over_gamma[10]*dx10; 
  alpha_cdim[24] = 1.414213562373095*p0_over_gamma[8]*dx10; 
  alpha_cdim[27] = 1.414213562373095*p0_over_gamma[14]*dx10; 
  alpha_cdim[32] = 1.414213562373095*p0_over_gamma[9]*dx10; 
  alpha_cdim[35] = 1.414213562373095*p0_over_gamma[16]*dx10; 
  cflFreq_mid += fabs(0.125*alpha_cdim[0]-0.1397542485937369*(alpha_cdim[32]+alpha_cdim[24])); 

  out[1] += 0.4330127018922193*(alpha_cdim[35]*f[35]+alpha_cdim[32]*f[32]+alpha_cdim[27]*f[27]+alpha_cdim[24]*f[24]+alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha_cdim[35]*f[38]+alpha_cdim[32]*f[34]+alpha_cdim[27]*f[30]+alpha_cdim[24]*f[26])+0.3872983346207416*(alpha_cdim[14]*f[22]+alpha_cdim[9]*f[19]+alpha_cdim[7]*f[18]+alpha_cdim[2]*f[16])+0.4330127018922193*(alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]+alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[6] += 0.4330127018922193*(alpha_cdim[32]*f[35]+f[32]*alpha_cdim[35])+0.3872983346207416*(alpha_cdim[14]*f[30]+alpha_cdim[10]*f[27]+f[10]*alpha_cdim[27]+alpha_cdim[7]*f[26]+alpha_cdim[3]*f[24]+f[3]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[8] += 0.3872983346207416*(alpha_cdim[14]*f[38]+alpha_cdim[10]*f[35]+f[10]*alpha_cdim[35]+alpha_cdim[9]*f[34]+alpha_cdim[4]*f[32]+f[4]*alpha_cdim[32])+0.4330127018922193*(alpha_cdim[24]*f[27]+f[24]*alpha_cdim[27]+alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14]+alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[11] += 0.4330127018922193*(alpha_cdim[32]*f[38]+f[34]*alpha_cdim[35])+0.3872983346207416*(alpha_cdim[10]*f[30]+alpha_cdim[14]*f[27]+f[14]*alpha_cdim[27]+alpha_cdim[3]*f[26]+alpha_cdim[7]*f[24]+f[7]*alpha_cdim[24]+alpha_cdim[9]*f[22]+alpha_cdim[14]*f[19]+alpha_cdim[2]*f[18]+alpha_cdim[7]*f[16])+0.4330127018922193*(alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14]+alpha_cdim[9]*f[10]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[0]*alpha_cdim[7]+alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]); 
  out[12] += 0.3872983346207416*(alpha_cdim[10]*f[38]+alpha_cdim[14]*f[35]+f[14]*alpha_cdim[35]+alpha_cdim[4]*f[34]+alpha_cdim[9]*f[32]+f[9]*alpha_cdim[32])+0.4330127018922193*(alpha_cdim[24]*f[30]+f[26]*alpha_cdim[27])+0.3872983346207416*(alpha_cdim[7]*f[22]+alpha_cdim[2]*f[19]+alpha_cdim[14]*f[18]+alpha_cdim[9]*f[16])+0.4330127018922193*(alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14]+alpha_cdim[7]*f[10]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[0]*alpha_cdim[9]+alpha_cdim[2]*f[4]+f[2]*alpha_cdim[4]); 
  out[13] += 0.3872983346207416*alpha_cdim[9]*f[38]+(0.3464101615137755*alpha_cdim[27]+0.3872983346207416*alpha_cdim[4])*f[35]+0.3464101615137755*f[27]*alpha_cdim[35]+0.3872983346207416*(f[4]*alpha_cdim[35]+alpha_cdim[14]*f[34]+alpha_cdim[10]*f[32]+f[10]*alpha_cdim[32]+alpha_cdim[7]*f[30]+alpha_cdim[3]*f[27]+f[3]*alpha_cdim[27]+alpha_cdim[14]*f[26]+alpha_cdim[10]*f[24]+f[10]*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[2]*f[14]+f[2]*alpha_cdim[14]+alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[15] += 0.3464101615137755*alpha_cdim[27]*f[38]+0.3872983346207416*(alpha_cdim[4]*f[38]+alpha_cdim[9]*f[35])+0.3464101615137755*f[30]*alpha_cdim[35]+0.3872983346207416*(f[9]*alpha_cdim[35]+alpha_cdim[10]*f[34]+alpha_cdim[14]*f[32]+f[14]*alpha_cdim[32]+alpha_cdim[3]*f[30]+alpha_cdim[7]*f[27]+f[7]*alpha_cdim[27]+alpha_cdim[10]*f[26]+alpha_cdim[14]*f[24]+f[14]*alpha_cdim[24]+alpha_cdim[2]*f[22]+alpha_cdim[7]*f[19]+alpha_cdim[9]*f[18]+alpha_cdim[14]*f[16])+0.4330127018922193*(alpha_cdim[0]*f[14]+f[0]*alpha_cdim[14]+alpha_cdim[2]*f[10]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+f[3]*alpha_cdim[9]+alpha_cdim[4]*f[7]+f[4]*alpha_cdim[7]); 
  out[17] += 0.4330127018922193*(alpha_cdim[10]*f[22]+alpha_cdim[4]*f[19]+alpha_cdim[3]*f[18]+alpha_cdim[0]*f[16])+0.3872983346207416*(alpha_cdim[14]*f[14]+alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[2]*f[2]); 
  out[20] += 0.3464101615137755*alpha_cdim[14]*f[30]+0.3872983346207416*f[22]*alpha_cdim[27]+0.3464101615137755*alpha_cdim[7]*f[26]+0.3872983346207416*f[18]*alpha_cdim[24]+0.4330127018922193*(alpha_cdim[4]*f[22]+alpha_cdim[10]*f[19]+alpha_cdim[0]*f[18]+alpha_cdim[3]*f[16])+0.3872983346207416*(alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]); 
  out[21] += 0.3464101615137755*alpha_cdim[14]*f[38]+0.3872983346207416*f[22]*alpha_cdim[35]+0.3464101615137755*alpha_cdim[9]*f[34]+0.3872983346207416*f[19]*alpha_cdim[32]+0.4330127018922193*(alpha_cdim[3]*f[22]+alpha_cdim[0]*f[19]+alpha_cdim[10]*f[18]+alpha_cdim[4]*f[16])+0.3872983346207416*(alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14]+alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]); 
  out[23] += 0.3464101615137755*alpha_cdim[9]*f[38]+0.3872983346207416*f[19]*alpha_cdim[35]+0.3464101615137755*alpha_cdim[14]*f[34]+0.3872983346207416*f[22]*alpha_cdim[32]+0.3464101615137755*alpha_cdim[7]*f[30]+0.3872983346207416*f[18]*alpha_cdim[27]+0.3464101615137755*alpha_cdim[14]*f[26]+0.3872983346207416*f[22]*alpha_cdim[24]+0.4330127018922193*(alpha_cdim[0]*f[22]+alpha_cdim[3]*f[19]+alpha_cdim[4]*f[18]+alpha_cdim[10]*f[16])+0.3872983346207416*(alpha_cdim[2]*f[14]+f[2]*alpha_cdim[14]+alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]); 
  out[25] += 0.3872983346207416*alpha_cdim[35]*f[35]+0.4330127018922193*alpha_cdim[9]*f[30]+0.276641667586244*alpha_cdim[27]*f[27]+0.4330127018922193*(alpha_cdim[4]*f[27]+f[4]*alpha_cdim[27]+alpha_cdim[2]*f[26])+0.276641667586244*alpha_cdim[24]*f[24]+0.4330127018922193*(alpha_cdim[0]*f[24]+f[0]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[7]*f[7]+alpha_cdim[3]*f[3]); 
  out[28] += 0.3872983346207416*alpha_cdim[35]*f[38]+0.276641667586244*alpha_cdim[27]*f[30]+0.4330127018922193*(alpha_cdim[4]*f[30]+alpha_cdim[9]*f[27]+f[9]*alpha_cdim[27])+0.276641667586244*alpha_cdim[24]*f[26]+0.4330127018922193*(alpha_cdim[0]*f[26]+alpha_cdim[2]*f[24]+f[2]*alpha_cdim[24])+0.3464101615137755*(alpha_cdim[14]*f[22]+alpha_cdim[7]*f[18])+0.3872983346207416*(alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]); 
  out[29] += 0.3464101615137755*(alpha_cdim[14]*f[38]+alpha_cdim[10]*f[35]+f[10]*alpha_cdim[35])+0.3872983346207416*(alpha_cdim[27]*f[32]+f[27]*alpha_cdim[32])+0.4330127018922193*alpha_cdim[2]*f[30]+(0.276641667586244*alpha_cdim[24]+0.4330127018922193*alpha_cdim[0])*f[27]+0.276641667586244*f[24]*alpha_cdim[27]+0.4330127018922193*(f[0]*alpha_cdim[27]+alpha_cdim[9]*f[26]+alpha_cdim[4]*f[24]+f[4]*alpha_cdim[24])+0.3872983346207416*(alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14]+alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]); 
  out[31] += 0.3464101615137755*(alpha_cdim[10]*f[38]+alpha_cdim[14]*f[35]+f[14]*alpha_cdim[35])+0.3872983346207416*alpha_cdim[27]*f[34]+f[30]*(0.3872983346207416*alpha_cdim[32]+0.276641667586244*alpha_cdim[24])+0.4330127018922193*(alpha_cdim[0]*f[30]+alpha_cdim[2]*f[27])+0.276641667586244*f[26]*alpha_cdim[27]+0.4330127018922193*(f[2]*alpha_cdim[27]+alpha_cdim[4]*f[26]+alpha_cdim[9]*f[24]+f[9]*alpha_cdim[24])+0.3464101615137755*(alpha_cdim[7]*f[22]+alpha_cdim[14]*f[18])+0.3872983346207416*(alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14]+alpha_cdim[7]*f[10]+f[7]*alpha_cdim[10]); 
  out[33] += 0.4330127018922193*alpha_cdim[7]*f[38]+0.276641667586244*alpha_cdim[35]*f[35]+0.4330127018922193*(alpha_cdim[3]*f[35]+f[3]*alpha_cdim[35]+alpha_cdim[2]*f[34])+0.276641667586244*alpha_cdim[32]*f[32]+0.4330127018922193*(alpha_cdim[0]*f[32]+f[0]*alpha_cdim[32])+0.3872983346207416*(alpha_cdim[27]*f[27]+alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[4]*f[4]); 
  out[36] += 0.276641667586244*alpha_cdim[35]*f[38]+0.4330127018922193*(alpha_cdim[3]*f[38]+alpha_cdim[7]*f[35]+f[7]*alpha_cdim[35])+0.276641667586244*alpha_cdim[32]*f[34]+0.4330127018922193*(alpha_cdim[0]*f[34]+alpha_cdim[2]*f[32]+f[2]*alpha_cdim[32])+0.3872983346207416*alpha_cdim[27]*f[30]+0.3464101615137755*(alpha_cdim[14]*f[22]+alpha_cdim[9]*f[19])+0.3872983346207416*(alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]); 
  out[37] += 0.4330127018922193*alpha_cdim[2]*f[38]+(0.276641667586244*alpha_cdim[32]+0.3872983346207416*alpha_cdim[24]+0.4330127018922193*alpha_cdim[0])*f[35]+(0.276641667586244*f[32]+0.3872983346207416*f[24])*alpha_cdim[35]+0.4330127018922193*(f[0]*alpha_cdim[35]+alpha_cdim[7]*f[34]+alpha_cdim[3]*f[32]+f[3]*alpha_cdim[32])+0.3464101615137755*(alpha_cdim[14]*f[30]+alpha_cdim[10]*f[27]+f[10]*alpha_cdim[27])+0.3872983346207416*(alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]); 
  out[39] += (0.276641667586244*alpha_cdim[32]+0.3872983346207416*alpha_cdim[24])*f[38]+0.4330127018922193*(alpha_cdim[0]*f[38]+alpha_cdim[2]*f[35])+(0.276641667586244*f[34]+0.3872983346207416*f[26])*alpha_cdim[35]+0.4330127018922193*(f[2]*alpha_cdim[35]+alpha_cdim[3]*f[34]+alpha_cdim[7]*f[32]+f[7]*alpha_cdim[32])+0.3464101615137755*(alpha_cdim[10]*f[30]+alpha_cdim[14]*f[27]+f[14]*alpha_cdim[27]+alpha_cdim[9]*f[22]+alpha_cdim[14]*f[19])+0.3872983346207416*(alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14]+alpha_cdim[9]*f[10]+f[9]*alpha_cdim[10]); 

  return 3.0*cflFreq_mid; 
} 
