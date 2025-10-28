#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_ext_phiA_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // f:         Input distribution function.
  // out:       Incremented output.

  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *phi = &pots[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[64] = {0.0}; 
  double alpha_vdim[64] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.3094010767585034*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_cdim[32] = 8.0*w1dx1; 
  alpha_cdim[36] = 2.3094010767585034*dv1dx1; 
  cflFreq_mid += 3.0*(fabs(w1dx1)+0.5*dv1dx1); 

  const double *Ax = &pots[4]; 
  const double *Ay = &pots[8]; 

  alpha_vdim[0] = dv10*(dx10*(3.4641016151377544*Ay[1]*wv2-3.4641016151377544*phi[1])-3.4641016151377544*Ax[2]*dx11*wv2); 
  alpha_vdim[1] = -(3.4641016151377544*Ax[3]*dv10*dx11*wv2); 
  alpha_vdim[2] = dv10*dx10*(3.4641016151377544*Ay[3]*wv2-3.4641016151377544*phi[3]); 
  alpha_vdim[4] = dv10*dv2*(Ay[1]*dx10-1.0*Ax[2]*dx11); 
  alpha_vdim[8] = -(1.0*Ax[3]*dv10*dv2*dx11); 
  alpha_vdim[9] = Ay[3]*dv10*dv2*dx10; 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[32] = dv11*(3.4641016151377544*Ax[2]*dx11*wv1-3.4641016151377544*(Ay[1]*dx10*wv1+phi[2]*dx11)); 
  alpha_vdim[33] = dv11*dx11*(3.4641016151377544*Ax[3]*wv1-3.4641016151377544*phi[3]); 
  alpha_vdim[34] = -(3.4641016151377544*Ay[3]*dv11*dx10*wv1); 
  alpha_vdim[35] = dv1*dv11*(Ax[2]*dx11-1.0*Ay[1]*dx10); 
  alpha_vdim[38] = Ax[3]*dv1*dv11*dx11; 
  alpha_vdim[39] = -(1.0*Ay[3]*dv1*dv11*dx10); 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[32]); 

  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[4]*alpha_cdim[36]+f[0]*alpha_cdim[32]); 
  out[3] += 0.4330127018922193*(alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[7]*alpha_vdim[39]+f[6]*alpha_vdim[38]+f[3]*alpha_vdim[35]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5] += 0.4330127018922193*(f[8]*alpha_cdim[36]+f[1]*alpha_cdim[32]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.38729833462074165*alpha_cdim[3]*f[16]+0.4330127018922193*(alpha_vdim[9]*f[12]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[2]*f[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.4330127018922193*(f[10]*alpha_cdim[36]+f[3]*alpha_cdim[32]+alpha_vdim[8]*f[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.4330127018922193*(f[11]*alpha_vdim[39]+f[3]*alpha_vdim[38]+f[6]*alpha_vdim[35]+f[5]*alpha_vdim[34]+f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[3]*f[10]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[3]*alpha_vdim[39]+f[11]*alpha_vdim[38])+0.38729833462074165*f[24]*alpha_cdim[36]+0.4330127018922193*(f[0]*alpha_cdim[36]+f[7]*alpha_vdim[35]+f[0]*alpha_vdim[34]+f[5]*alpha_vdim[33]+f[2]*alpha_vdim[32]+f[4]*alpha_cdim[32]); 
  out[10] += (0.38729833462074165*f[18]+0.4330127018922193*f[2])*alpha_vdim[39]+(0.38729833462074165*f[17]+0.4330127018922193*f[1])*alpha_vdim[38]+0.38729833462074165*f[16]*alpha_vdim[35]+0.4330127018922193*(f[0]*alpha_vdim[35]+f[7]*alpha_vdim[34]+f[6]*alpha_vdim[33]+f[3]*alpha_vdim[32])+0.38729833462074165*(alpha_vdim[9]*f[26]+alpha_vdim[8]*f[25]+alpha_vdim[4]*f[24])+0.4330127018922193*(alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[11] += 0.4330127018922193*(f[13]*alpha_cdim[36]+f[6]*alpha_cdim[32])+0.38729833462074165*alpha_cdim[3]*f[18]+0.4330127018922193*(alpha_vdim[4]*f[12]+alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[5]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[12] += 0.4330127018922193*(f[6]*alpha_vdim[39]+f[7]*alpha_vdim[38])+0.38729833462074165*f[25]*alpha_cdim[36]+0.4330127018922193*(f[1]*alpha_cdim[36]+f[11]*alpha_vdim[35]+f[1]*alpha_vdim[34]+f[2]*alpha_vdim[33]+f[5]*alpha_vdim[32]+f[8]*alpha_cdim[32]+alpha_cdim[3]*f[14]+alpha_cdim[0]*f[9]); 
  out[13] += (0.38729833462074165*f[20]+0.4330127018922193*f[5])*alpha_vdim[39]+(0.38729833462074165*f[16]+0.4330127018922193*f[0])*alpha_vdim[38]+0.38729833462074165*f[17]*alpha_vdim[35]+0.4330127018922193*(f[1]*alpha_vdim[35]+f[11]*alpha_vdim[34]+f[3]*alpha_vdim[33]+f[6]*alpha_vdim[32])+0.38729833462074165*(alpha_vdim[9]*f[28]+alpha_vdim[4]*f[25]+alpha_vdim[8]*f[24]+alpha_cdim[3]*f[19])+0.4330127018922193*(alpha_vdim[2]*f[12]+alpha_cdim[0]*f[10]+f[5]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[3]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[14] += (0.38729833462074165*f[16]+0.4330127018922193*f[0])*alpha_vdim[39]+(0.38729833462074165*f[20]+0.4330127018922193*f[5])*alpha_vdim[38]+(0.38729833462074165*f[27]+0.4330127018922193*f[3])*alpha_cdim[36]+0.38729833462074165*f[18]*alpha_vdim[35]+0.4330127018922193*(f[2]*alpha_vdim[35]+f[3]*alpha_vdim[34]+f[11]*alpha_vdim[33]+f[7]*alpha_vdim[32]+f[10]*alpha_cdim[32])+0.38729833462074165*(alpha_vdim[8]*f[28]+alpha_vdim[4]*f[26]+alpha_vdim[9]*f[24])+0.4330127018922193*(alpha_vdim[1]*f[12]+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+f[5]*alpha_vdim[8]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[15] += (0.38729833462074165*f[17]+0.4330127018922193*f[1])*alpha_vdim[39]+(0.38729833462074165*f[18]+0.4330127018922193*f[2])*alpha_vdim[38]+(0.38729833462074165*f[29]+0.4330127018922193*f[6])*alpha_cdim[36]+0.38729833462074165*f[20]*alpha_vdim[35]+0.4330127018922193*(f[5]*alpha_vdim[35]+f[6]*alpha_vdim[34]+f[7]*alpha_vdim[33]+f[11]*alpha_vdim[32]+f[13]*alpha_cdim[32])+0.38729833462074165*(alpha_vdim[4]*f[28]+alpha_vdim[8]*f[26]+alpha_vdim[9]*f[25]+alpha_cdim[3]*f[22])+0.4330127018922193*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[12]+(alpha_cdim[3]+alpha_vdim[1])*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8]+alpha_vdim[4]*f[5]); 
  out[16] += 0.9682458365518543*(alpha_vdim[9]*f[14]+alpha_vdim[8]*f[13]+alpha_vdim[4]*f[10]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[17] += 0.4330127018922193*alpha_cdim[0]*f[16]+0.9682458365518543*(alpha_vdim[9]*f[15]+alpha_vdim[4]*f[13]+alpha_vdim[2]*f[11]+alpha_vdim[8]*f[10]+alpha_vdim[0]*f[6])+(0.38729833462074165*alpha_cdim[3]+0.9682458365518543*alpha_vdim[1])*f[3]; 
  out[18] += 0.4330127018922193*(f[19]*alpha_cdim[36]+f[16]*alpha_cdim[32])+0.9682458365518543*(alpha_vdim[8]*f[15]+alpha_vdim[4]*f[14]+alpha_vdim[1]*f[11]+alpha_vdim[9]*f[10]+alpha_vdim[0]*f[7]+alpha_vdim[2]*f[3]); 
  out[19] += 0.38729833462074165*(f[7]*alpha_vdim[39]+f[6]*alpha_vdim[38]+f[3]*alpha_vdim[35])+0.4330127018922193*(f[18]*alpha_vdim[34]+f[17]*alpha_vdim[33]+f[16]*alpha_vdim[32])+0.8660254037844386*(alpha_vdim[9]*f[30]+alpha_vdim[8]*f[29]+alpha_vdim[4]*f[27])+0.9682458365518543*(alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[10]+f[7]*alpha_vdim[9]+f[6]*alpha_vdim[8]+f[3]*alpha_vdim[4]); 
  out[20] += 0.4330127018922193*(f[21]*alpha_cdim[36]+f[17]*alpha_cdim[32]+alpha_cdim[0]*f[18])+0.9682458365518543*(alpha_vdim[4]*f[15]+alpha_vdim[8]*f[14]+alpha_vdim[9]*f[13]+alpha_vdim[0]*f[11])+0.38729833462074165*alpha_cdim[3]*f[7]+0.9682458365518543*(alpha_vdim[1]*f[7]+alpha_vdim[2]*f[6]); 
  out[21] += 0.38729833462074165*(f[11]*alpha_vdim[39]+f[3]*alpha_vdim[38]+f[6]*alpha_vdim[35])+0.4330127018922193*(f[20]*alpha_vdim[34]+f[16]*alpha_vdim[33]+f[17]*alpha_vdim[32])+0.8660254037844386*(alpha_vdim[9]*f[31]+alpha_vdim[4]*f[29]+alpha_vdim[8]*f[27])+0.4330127018922193*alpha_cdim[0]*f[19]+0.9682458365518543*(alpha_vdim[2]*f[15]+alpha_vdim[0]*f[13]+alpha_vdim[9]*f[11])+0.38729833462074165*alpha_cdim[3]*f[10]+0.9682458365518543*(alpha_vdim[1]*f[10]+f[3]*alpha_vdim[8]+alpha_vdim[4]*f[6]); 
  out[22] += 0.38729833462074165*(f[3]*alpha_vdim[39]+f[11]*alpha_vdim[38])+0.4330127018922193*f[16]*alpha_cdim[36]+0.38729833462074165*f[7]*alpha_vdim[35]+0.4330127018922193*(f[16]*alpha_vdim[34]+f[20]*alpha_vdim[33]+f[18]*alpha_vdim[32]+f[19]*alpha_cdim[32])+0.8660254037844386*(alpha_vdim[8]*f[31]+alpha_vdim[4]*f[30]+alpha_vdim[9]*f[27])+0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[14]+alpha_vdim[8]*f[11]+alpha_vdim[2]*f[10]+f[3]*alpha_vdim[9]+alpha_vdim[4]*f[7]); 
  out[23] += 0.38729833462074165*(f[6]*alpha_vdim[39]+f[7]*alpha_vdim[38])+0.4330127018922193*f[17]*alpha_cdim[36]+0.38729833462074165*f[11]*alpha_vdim[35]+0.4330127018922193*(f[17]*alpha_vdim[34]+f[18]*alpha_vdim[33]+f[20]*alpha_vdim[32]+f[21]*alpha_cdim[32])+0.8660254037844386*(alpha_vdim[4]*f[31]+alpha_vdim[8]*f[30]+alpha_vdim[9]*f[29])+0.4330127018922193*alpha_cdim[0]*f[22]+0.9682458365518543*alpha_vdim[0]*f[15]+0.38729833462074165*alpha_cdim[3]*f[14]+0.9682458365518543*(alpha_vdim[1]*f[14]+alpha_vdim[2]*f[13]+alpha_vdim[4]*f[11]+f[6]*alpha_vdim[9]+f[7]*alpha_vdim[8]); 
  out[24] += 0.9682458365518543*(f[14]*alpha_vdim[39]+f[13]*alpha_vdim[38]+f[10]*alpha_vdim[35]+f[9]*alpha_vdim[34]+f[8]*alpha_vdim[33]+f[4]*alpha_vdim[32]); 
  out[25] += 0.9682458365518543*(f[15]*alpha_vdim[39]+f[10]*alpha_vdim[38]+f[13]*alpha_vdim[35]+f[12]*alpha_vdim[34]+f[4]*alpha_vdim[33]+f[8]*alpha_vdim[32])+0.4330127018922193*(alpha_cdim[3]*f[27]+alpha_cdim[0]*f[24]); 
  out[26] += 0.9682458365518543*(f[10]*alpha_vdim[39]+f[15]*alpha_vdim[38])+0.38729833462074165*f[4]*alpha_cdim[36]+0.9682458365518543*(f[14]*alpha_vdim[35]+f[4]*alpha_vdim[34]+f[12]*alpha_vdim[33]+f[9]*alpha_vdim[32])+0.4330127018922193*f[24]*alpha_cdim[32]; 
  out[27] += (0.8660254037844386*f[22]+0.9682458365518543*f[9])*alpha_vdim[39]+(0.8660254037844386*f[21]+0.9682458365518543*f[8])*alpha_vdim[38]+0.8660254037844386*f[19]*alpha_vdim[35]+0.9682458365518543*(f[4]*alpha_vdim[35]+f[14]*alpha_vdim[34]+f[13]*alpha_vdim[33]+f[10]*alpha_vdim[32])+0.4330127018922193*(alpha_vdim[2]*f[26]+alpha_vdim[1]*f[25]+alpha_vdim[0]*f[24])+0.38729833462074165*(alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[28] += 0.9682458365518543*(f[13]*alpha_vdim[39]+f[14]*alpha_vdim[38])+0.38729833462074165*f[8]*alpha_cdim[36]+0.9682458365518543*(f[15]*alpha_vdim[35]+f[8]*alpha_vdim[34]+f[9]*alpha_vdim[33]+f[12]*alpha_vdim[32])+0.4330127018922193*(f[25]*alpha_cdim[32]+alpha_cdim[3]*f[30]+alpha_cdim[0]*f[26]); 
  out[29] += (0.8660254037844386*f[23]+0.9682458365518543*f[12])*alpha_vdim[39]+(0.8660254037844386*f[19]+0.9682458365518543*f[4])*alpha_vdim[38]+0.8660254037844386*f[21]*alpha_vdim[35]+0.9682458365518543*(f[8]*alpha_vdim[35]+f[15]*alpha_vdim[34]+f[10]*alpha_vdim[33]+f[13]*alpha_vdim[32])+0.4330127018922193*(alpha_vdim[2]*f[28]+alpha_cdim[0]*f[27]+alpha_vdim[0]*f[25]+(alpha_cdim[3]+alpha_vdim[1])*f[24])+0.38729833462074165*(alpha_vdim[9]*f[12]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]); 
  out[30] += (0.8660254037844386*f[19]+0.9682458365518543*f[4])*alpha_vdim[39]+(0.8660254037844386*f[23]+0.9682458365518543*f[12])*alpha_vdim[38]+0.38729833462074165*f[10]*alpha_cdim[36]+0.8660254037844386*f[22]*alpha_vdim[35]+0.9682458365518543*(f[9]*alpha_vdim[35]+f[10]*alpha_vdim[34]+f[15]*alpha_vdim[33]+f[14]*alpha_vdim[32])+0.4330127018922193*(f[27]*alpha_cdim[32]+alpha_vdim[1]*f[28]+alpha_vdim[0]*f[26]+alpha_vdim[2]*f[24])+0.38729833462074165*(alpha_vdim[8]*f[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]); 
  out[31] += (0.8660254037844386*f[21]+0.9682458365518543*f[8])*alpha_vdim[39]+(0.8660254037844386*f[22]+0.9682458365518543*f[9])*alpha_vdim[38]+0.38729833462074165*f[13]*alpha_cdim[36]+0.8660254037844386*f[23]*alpha_vdim[35]+0.9682458365518543*(f[12]*alpha_vdim[35]+f[13]*alpha_vdim[34]+f[14]*alpha_vdim[33]+f[15]*alpha_vdim[32])+0.4330127018922193*(f[29]*alpha_cdim[32]+alpha_cdim[0]*f[30]+alpha_vdim[0]*f[28]+(alpha_cdim[3]+alpha_vdim[1])*f[26]+alpha_vdim[2]*f[25])+0.38729833462074165*(alpha_vdim[4]*f[12]+alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]); 

  return cflFreq_mid; 
} 

