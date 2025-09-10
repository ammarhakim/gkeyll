#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
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
  double alpha_cdim[96] = {0.0}; 
  alpha_cdim[0] = 2.0*p0_over_gamma[0]*dx10; 
  alpha_cdim[3] = 2.0*p0_over_gamma[1]*dx10; 
  alpha_cdim[4] = 2.0*p0_over_gamma[2]*dx10; 
  alpha_cdim[10] = 2.0*p0_over_gamma[3]*dx10; 
  alpha_cdim[14] = 2.0*p0_over_gamma[5]*dx10; 
  cflFreq_mid += fabs(0.125*alpha_cdim[0]-0.1397542485937369*alpha_cdim[14]); 

  alpha_cdim[48] = 2.0*p1_over_gamma[0]*dx11; 
  alpha_cdim[51] = 2.0*p1_over_gamma[1]*dx11; 
  alpha_cdim[52] = 2.0*p1_over_gamma[2]*dx11; 
  alpha_cdim[58] = 2.0*p1_over_gamma[3]*dx11; 
  alpha_cdim[61] = 2.0*p1_over_gamma[4]*dx11; 
  cflFreq_mid += fabs(0.125*alpha_cdim[48]-0.1397542485937369*alpha_cdim[61]); 

  out[1] += 0.4330127018922193*(alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[13]*alpha_cdim[61]+f[10]*alpha_cdim[58]+f[4]*alpha_cdim[52]+f[3]*alpha_cdim[51]+f[0]*alpha_cdim[48]); 
  out[5] += 0.4330127018922193*(f[23]*alpha_cdim[61]+f[17]*alpha_cdim[58]+f[8]*alpha_cdim[52]+f[6]*alpha_cdim[51]+f[1]*alpha_cdim[48]+alpha_cdim[14]*f[29]+alpha_cdim[10]*f[18]+alpha_cdim[4]*f[9]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.4330127018922193*alpha_cdim[14]*f[30]+0.3872983346207416*(alpha_cdim[10]*f[27]+alpha_cdim[3]*f[13])+0.4330127018922193*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[7] += 0.3872983346207416*(f[3]*alpha_cdim[61]+f[27]*alpha_cdim[58])+0.4330127018922193*(f[4]*alpha_cdim[58]+f[10]*alpha_cdim[52])+0.3872983346207416*f[13]*alpha_cdim[51]+0.4330127018922193*(f[0]*alpha_cdim[51]+f[3]*alpha_cdim[48]); 
  out[8] += 0.3872983346207416*(alpha_cdim[10]*f[30]+alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14])+0.4330127018922193*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[9] += 0.4330127018922193*f[27]*alpha_cdim[61]+(0.3872983346207416*f[30]+0.4330127018922193*f[3])*alpha_cdim[58]+0.3872983346207416*f[14]*alpha_cdim[52]+0.4330127018922193*(f[0]*alpha_cdim[52]+f[10]*alpha_cdim[51]+f[4]*alpha_cdim[48]); 
  out[11] += 0.9682458365518543*(alpha_cdim[14]*f[28]+alpha_cdim[10]*f[17]+alpha_cdim[4]*f[8]+alpha_cdim[3]*f[6]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(f[24]*alpha_cdim[61]+f[18]*alpha_cdim[58]+f[9]*alpha_cdim[52]+f[7]*alpha_cdim[51]+f[2]*alpha_cdim[48]); 
  out[15] += 0.3872983346207416*(f[6]*alpha_cdim[61]+f[39]*alpha_cdim[58])+0.4330127018922193*(f[8]*alpha_cdim[58]+f[17]*alpha_cdim[52])+0.3872983346207416*f[23]*alpha_cdim[51]+0.4330127018922193*(f[1]*alpha_cdim[51]+f[6]*alpha_cdim[48]+alpha_cdim[14]*f[43])+0.3872983346207416*(alpha_cdim[10]*f[40]+alpha_cdim[3]*f[24])+0.4330127018922193*(alpha_cdim[4]*f[18]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[2]*alpha_cdim[3]); 
  out[16] += 0.4330127018922193*f[39]*alpha_cdim[61]+(0.3872983346207416*f[42]+0.4330127018922193*f[6])*alpha_cdim[58]+0.3872983346207416*f[28]*alpha_cdim[52]+0.4330127018922193*(f[1]*alpha_cdim[52]+f[17]*alpha_cdim[51]+f[8]*alpha_cdim[48])+0.3872983346207416*(alpha_cdim[10]*f[43]+alpha_cdim[4]*f[29])+0.4330127018922193*alpha_cdim[3]*f[18]+0.3872983346207416*f[9]*alpha_cdim[14]+0.4330127018922193*(f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[2]*alpha_cdim[4]); 
  out[17] += 0.3872983346207416*(alpha_cdim[4]*f[30]+alpha_cdim[3]*f[27]+alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[10]*f[13])+0.4330127018922193*(alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[18] += 0.3872983346207416*f[10]*alpha_cdim[61]+(0.3872983346207416*(f[14]+f[13])+0.4330127018922193*f[0])*alpha_cdim[58]+(0.3872983346207416*f[30]+0.4330127018922193*f[3])*alpha_cdim[52]+0.3872983346207416*f[27]*alpha_cdim[51]+0.4330127018922193*(f[4]*alpha_cdim[51]+f[10]*alpha_cdim[48]); 
  out[19] += 0.4330127018922193*(f[37]*alpha_cdim[58]+f[25]*alpha_cdim[52]+f[21]*alpha_cdim[51]+f[11]*alpha_cdim[48])+0.9682458365518543*(alpha_cdim[14]*f[41]+alpha_cdim[10]*f[31]+alpha_cdim[4]*f[16]+alpha_cdim[3]*f[15]+alpha_cdim[0]*f[5]); 
  out[20] += 0.9682458365518543*(f[34]*alpha_cdim[61]+f[31]*alpha_cdim[58]+f[16]*alpha_cdim[52]+f[15]*alpha_cdim[51]+f[5]*alpha_cdim[48])+0.4330127018922193*(alpha_cdim[10]*f[38]+alpha_cdim[4]*f[26]+alpha_cdim[3]*f[22]+alpha_cdim[0]*f[12]); 
  out[21] += 0.9682458365518543*alpha_cdim[14]*f[42]+0.8660254037844386*(alpha_cdim[10]*f[39]+alpha_cdim[3]*f[23])+0.9682458365518543*(alpha_cdim[4]*f[17]+f[8]*alpha_cdim[10]+alpha_cdim[0]*f[6]+f[1]*alpha_cdim[3]); 
  out[22] += 0.8660254037844386*(f[7]*alpha_cdim[61]+f[40]*alpha_cdim[58])+0.9682458365518543*(f[9]*alpha_cdim[58]+f[18]*alpha_cdim[52])+0.8660254037844386*f[24]*alpha_cdim[51]+0.9682458365518543*(f[2]*alpha_cdim[51]+f[7]*alpha_cdim[48]); 
  out[23] += 0.4330127018922193*(alpha_cdim[4]*f[27]+alpha_cdim[0]*f[13])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[3]*f[3]); 
  out[24] += (0.276641667586244*f[13]+0.4330127018922193*f[0])*alpha_cdim[61]+0.3872983346207416*f[10]*alpha_cdim[58]+0.4330127018922193*f[27]*alpha_cdim[52]+0.3872983346207416*f[3]*alpha_cdim[51]+0.4330127018922193*f[13]*alpha_cdim[48]; 
  out[25] += 0.8660254037844386*(alpha_cdim[10]*f[42]+alpha_cdim[4]*f[28])+0.9682458365518543*alpha_cdim[3]*f[17]+0.8660254037844386*f[8]*alpha_cdim[14]+0.9682458365518543*(f[6]*alpha_cdim[10]+alpha_cdim[0]*f[8]+f[1]*alpha_cdim[4]); 
  out[26] += 0.9682458365518543*f[40]*alpha_cdim[61]+(0.8660254037844386*f[43]+0.9682458365518543*f[7])*alpha_cdim[58]+0.8660254037844386*f[29]*alpha_cdim[52]+0.9682458365518543*(f[2]*alpha_cdim[52]+f[18]*alpha_cdim[51]+f[9]*alpha_cdim[48]); 
  out[28] += 0.4330127018922193*alpha_cdim[3]*f[30]+0.276641667586244*alpha_cdim[14]*f[14]+0.4330127018922193*(alpha_cdim[0]*f[14]+f[0]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]); 
  out[29] += 0.3872983346207416*(f[10]*alpha_cdim[58]+f[4]*alpha_cdim[52])+0.4330127018922193*(f[30]*alpha_cdim[51]+f[14]*alpha_cdim[48]); 
  out[31] += 0.3872983346207416*f[17]*alpha_cdim[61]+(0.3872983346207416*(f[28]+f[23])+0.4330127018922193*f[1])*alpha_cdim[58]+(0.3872983346207416*f[42]+0.4330127018922193*f[6])*alpha_cdim[52]+0.3872983346207416*f[39]*alpha_cdim[51]+0.4330127018922193*(f[8]*alpha_cdim[51]+f[17]*alpha_cdim[48])+0.3872983346207416*(alpha_cdim[4]*f[43]+alpha_cdim[3]*f[40]+alpha_cdim[10]*(f[29]+f[24])+alpha_cdim[14]*f[18])+0.4330127018922193*(alpha_cdim[0]*f[18]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+alpha_cdim[4]*f[7]); 
  out[32] += 0.3872983346207416*f[21]*alpha_cdim[61]+0.4330127018922193*(f[25]*alpha_cdim[58]+f[37]*alpha_cdim[52]+f[11]*alpha_cdim[51]+f[21]*alpha_cdim[48])+0.9682458365518543*alpha_cdim[14]*f[47]+0.8660254037844386*(alpha_cdim[10]*f[46]+alpha_cdim[3]*f[34])+0.9682458365518543*(alpha_cdim[4]*f[31]+alpha_cdim[10]*f[16]+alpha_cdim[0]*f[15]+alpha_cdim[3]*f[5]); 
  out[33] += 0.8660254037844386*(f[15]*alpha_cdim[61]+f[46]*alpha_cdim[58])+0.9682458365518543*(f[16]*alpha_cdim[58]+f[31]*alpha_cdim[52])+0.8660254037844386*f[34]*alpha_cdim[51]+0.9682458365518543*(f[5]*alpha_cdim[51]+f[15]*alpha_cdim[48])+0.4330127018922193*(alpha_cdim[4]*f[38]+alpha_cdim[10]*f[26]+alpha_cdim[0]*f[22]+alpha_cdim[3]*f[12]); 
  out[34] += (0.276641667586244*f[23]+0.4330127018922193*f[1])*alpha_cdim[61]+0.3872983346207416*f[17]*alpha_cdim[58]+0.4330127018922193*f[39]*alpha_cdim[52]+0.3872983346207416*f[6]*alpha_cdim[51]+0.4330127018922193*(f[23]*alpha_cdim[48]+alpha_cdim[4]*f[40]+alpha_cdim[0]*f[24])+0.3872983346207416*(alpha_cdim[10]*f[18]+alpha_cdim[3]*f[7]); 
  out[35] += 0.4330127018922193*(f[21]*alpha_cdim[58]+f[11]*alpha_cdim[52]+f[37]*alpha_cdim[51]+f[25]*alpha_cdim[48])+0.8660254037844386*(alpha_cdim[10]*f[47]+alpha_cdim[4]*f[41])+0.9682458365518543*alpha_cdim[3]*f[31]+0.8660254037844386*alpha_cdim[14]*f[16]+0.9682458365518543*(alpha_cdim[0]*f[16]+alpha_cdim[10]*f[15]+alpha_cdim[4]*f[5]); 
  out[36] += 0.9682458365518543*f[46]*alpha_cdim[61]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alpha_cdim[58]+0.8660254037844386*f[41]*alpha_cdim[52]+0.9682458365518543*(f[5]*alpha_cdim[52]+f[31]*alpha_cdim[51]+f[16]*alpha_cdim[48])+0.4330127018922193*alpha_cdim[3]*f[38]+0.3872983346207416*alpha_cdim[14]*f[26]+0.4330127018922193*(alpha_cdim[0]*f[26]+alpha_cdim[10]*f[22]+alpha_cdim[4]*f[12]); 
  out[37] += 0.8660254037844386*(alpha_cdim[4]*f[42]+alpha_cdim[3]*f[39]+alpha_cdim[10]*(f[28]+f[23])+alpha_cdim[14]*f[17])+0.9682458365518543*(alpha_cdim[0]*f[17]+f[1]*alpha_cdim[10]+alpha_cdim[3]*f[8]+alpha_cdim[4]*f[6]); 
  out[38] += 0.8660254037844386*f[18]*alpha_cdim[61]+(0.8660254037844386*(f[29]+f[24])+0.9682458365518543*f[2])*alpha_cdim[58]+(0.8660254037844386*f[43]+0.9682458365518543*f[7])*alpha_cdim[52]+0.8660254037844386*f[40]*alpha_cdim[51]+0.9682458365518543*(f[9]*alpha_cdim[51]+f[18]*alpha_cdim[48]); 
  out[39] += 0.3464101615137755*alpha_cdim[10]*f[30]+0.3872983346207416*alpha_cdim[14]*f[27]+0.4330127018922193*(alpha_cdim[0]*f[27]+alpha_cdim[4]*f[13])+0.3872983346207416*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]); 
  out[40] += (0.276641667586244*f[27]+0.4330127018922193*f[4])*alpha_cdim[61]+(0.3464101615137755*f[30]+0.3872983346207416*f[3])*alpha_cdim[58]+0.4330127018922193*f[13]*alpha_cdim[52]+0.3872983346207416*f[10]*alpha_cdim[51]+0.4330127018922193*f[27]*alpha_cdim[48]; 
  out[41] += 0.3872983346207416*(f[17]*alpha_cdim[58]+f[8]*alpha_cdim[52])+0.4330127018922193*(f[42]*alpha_cdim[51]+f[28]*alpha_cdim[48]+alpha_cdim[3]*f[43])+(0.276641667586244*alpha_cdim[14]+0.4330127018922193*alpha_cdim[0])*f[29]+0.3872983346207416*alpha_cdim[10]*f[18]+0.4330127018922193*f[2]*alpha_cdim[14]+0.3872983346207416*alpha_cdim[4]*f[9]; 
  out[42] += (0.276641667586244*alpha_cdim[14]+0.4330127018922193*alpha_cdim[0])*f[30]+0.3464101615137755*alpha_cdim[10]*f[27]+0.4330127018922193*(alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14])+0.3872983346207416*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]); 
  out[43] += 0.3872983346207416*f[30]*alpha_cdim[61]+0.3464101615137755*f[27]*alpha_cdim[58]+0.3872983346207416*(f[4]*alpha_cdim[58]+f[10]*alpha_cdim[52])+0.4330127018922193*(f[14]*alpha_cdim[51]+f[30]*alpha_cdim[48]); 
  out[44] += 0.3872983346207416*f[37]*alpha_cdim[61]+0.4330127018922193*(f[11]*alpha_cdim[58]+f[21]*alpha_cdim[52]+f[25]*alpha_cdim[51]+f[37]*alpha_cdim[48])+0.8660254037844386*(alpha_cdim[4]*f[47]+alpha_cdim[3]*f[46]+alpha_cdim[10]*(f[41]+f[34])+alpha_cdim[14]*f[31])+0.9682458365518543*(alpha_cdim[0]*f[31]+alpha_cdim[3]*f[16]+alpha_cdim[4]*f[15]+f[5]*alpha_cdim[10]); 
  out[45] += 0.8660254037844386*f[31]*alpha_cdim[61]+(0.8660254037844386*(f[41]+f[34])+0.9682458365518543*f[5])*alpha_cdim[58]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alpha_cdim[52]+0.8660254037844386*f[46]*alpha_cdim[51]+0.9682458365518543*(f[16]*alpha_cdim[51]+f[31]*alpha_cdim[48])+0.3872983346207416*alpha_cdim[14]*f[38]+0.4330127018922193*(alpha_cdim[0]*f[38]+alpha_cdim[3]*f[26]+alpha_cdim[4]*f[22]+alpha_cdim[10]*f[12]); 
  out[46] += (0.276641667586244*f[39]+0.4330127018922193*f[8])*alpha_cdim[61]+(0.3464101615137755*f[42]+0.3872983346207416*f[6])*alpha_cdim[58]+0.4330127018922193*f[23]*alpha_cdim[52]+0.3872983346207416*f[17]*alpha_cdim[51]+0.4330127018922193*f[39]*alpha_cdim[48]+0.3464101615137755*alpha_cdim[10]*f[43]+0.3872983346207416*alpha_cdim[14]*f[40]+0.4330127018922193*(alpha_cdim[0]*f[40]+alpha_cdim[4]*f[24])+0.3872983346207416*(alpha_cdim[3]*f[18]+f[7]*alpha_cdim[10]); 
  out[47] += 0.3872983346207416*f[42]*alpha_cdim[61]+0.3464101615137755*f[39]*alpha_cdim[58]+0.3872983346207416*(f[8]*alpha_cdim[58]+f[17]*alpha_cdim[52])+0.4330127018922193*(f[28]*alpha_cdim[51]+f[42]*alpha_cdim[48])+(0.276641667586244*alpha_cdim[14]+0.4330127018922193*alpha_cdim[0])*f[43]+0.3464101615137755*alpha_cdim[10]*f[40]+0.4330127018922193*alpha_cdim[3]*f[29]+0.3872983346207416*alpha_cdim[4]*f[18]+0.4330127018922193*f[7]*alpha_cdim[14]+0.3872983346207416*f[9]*alpha_cdim[10]; 

  return 5.0*cflFreq_mid; 
} 
