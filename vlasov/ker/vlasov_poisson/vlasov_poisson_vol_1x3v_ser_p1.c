#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // f:         Input distribution function.
  // out:       Incremented output.

  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 

  const double *phi = &pots[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[40] = {0.0}; 
  double alpha_vdim[120] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.3094010767585034*dv0dx0; 

  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = -(4.898979485566357*phi[1]*dv10*dx10); 

  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 



  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*alpha_vdim[0]*f[0]; 
  out[5] += 0.38729833462074165*alpha_cdim[2]*f[16]+0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[6] += 0.4330127018922193*(alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*alpha_vdim[0]*f[3]; 
  out[8] += 0.4330127018922193*(alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*alpha_vdim[0]*f[4]; 
  out[11] += 0.38729833462074165*alpha_cdim[2]*f[18]+0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+alpha_cdim[2]*f[3]); 
  out[12] += 0.38729833462074165*alpha_cdim[2]*f[19]+0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+alpha_cdim[2]*f[4]); 
  out[13] += 0.4330127018922193*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[10]); 
  out[14] += 0.4330127018922193*alpha_vdim[0]*f[10]; 
  out[15] += 0.38729833462074165*alpha_cdim[2]*f[22]+0.4330127018922193*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+alpha_cdim[2]*f[10]); 
  out[16] += 0.9682458365518543*alpha_vdim[0]*f[2]; 
  out[17] += 0.4330127018922193*alpha_cdim[0]*f[16]+0.9682458365518543*alpha_vdim[0]*f[5]+0.38729833462074165*alpha_cdim[2]*f[2]; 
  out[18] += 0.9682458365518543*alpha_vdim[0]*f[7]; 
  out[19] += 0.9682458365518543*alpha_vdim[0]*f[9]; 
  out[20] += 0.4330127018922193*alpha_cdim[0]*f[18]+0.9682458365518543*alpha_vdim[0]*f[11]+0.38729833462074165*alpha_cdim[2]*f[7]; 
  out[21] += 0.4330127018922193*alpha_cdim[0]*f[19]+0.9682458365518543*alpha_vdim[0]*f[12]+0.38729833462074165*alpha_cdim[2]*f[9]; 
  out[22] += 0.9682458365518543*alpha_vdim[0]*f[14]; 
  out[23] += 0.4330127018922193*alpha_cdim[0]*f[22]+0.9682458365518543*alpha_vdim[0]*f[15]+0.38729833462074165*alpha_cdim[2]*f[14]; 
  out[25] += 0.4330127018922193*(alpha_cdim[2]*f[26]+alpha_cdim[0]*f[24]); 
  out[26] += 0.4330127018922193*alpha_vdim[0]*f[24]; 
  out[28] += 0.4330127018922193*(alpha_cdim[0]*f[26]+alpha_vdim[0]*f[25]+alpha_cdim[2]*f[24]); 
  out[29] += 0.4330127018922193*(alpha_cdim[2]*f[30]+alpha_cdim[0]*f[27]); 
  out[30] += 0.4330127018922193*alpha_vdim[0]*f[27]; 
  out[31] += 0.4330127018922193*(alpha_cdim[0]*f[30]+alpha_vdim[0]*f[29]+alpha_cdim[2]*f[27]); 
  out[33] += 0.4330127018922193*(alpha_cdim[2]*f[34]+alpha_cdim[0]*f[32]); 
  out[34] += 0.4330127018922193*alpha_vdim[0]*f[32]; 
  out[36] += 0.4330127018922193*(alpha_cdim[0]*f[34]+alpha_vdim[0]*f[33]+alpha_cdim[2]*f[32]); 
  out[37] += 0.4330127018922193*(alpha_cdim[2]*f[38]+alpha_cdim[0]*f[35]); 
  out[38] += 0.4330127018922193*alpha_vdim[0]*f[35]; 
  out[39] += 0.4330127018922193*(alpha_cdim[0]*f[38]+alpha_vdim[0]*f[37]+alpha_cdim[2]*f[35]); 

  return cflFreq_mid; 
} 

