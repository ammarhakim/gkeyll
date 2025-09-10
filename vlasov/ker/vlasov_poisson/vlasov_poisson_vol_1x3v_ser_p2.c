#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out) 
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

  double alpha_cdim[48] = {0.0}; 
  double alpha_vdim[144] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.3094010767585034*dv0dx0; 

  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = -(4.898979485566357*phi[1]*dv10*dx10); 
  alpha_vdim[1] = -(10.954451150103324*phi[2]*dv10*dx10); 

  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 



  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.38729833462074165*(alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11])+0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[15] += 0.38729833462074165*(alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21])+0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[16] += 0.38729833462074165*(alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25])+0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[17] += 0.4330127018922193*(alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += 0.4330127018922193*(alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]); 
  out[19] += 0.8660254037844386*alpha_cdim[2]*f[20]+0.4330127018922193*alpha_vdim[0]*f[11]+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.38729833462074165*alpha_vdim[1]); 
  out[20] += 0.8660254037844386*alpha_vdim[1]*f[19]+0.4330127018922193*alpha_cdim[0]*f[12]+0.9682458365518543*alpha_vdim[0]*f[5]+(0.38729833462074165*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += 0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]); 
  out[23] += 0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += 0.4330127018922193*(alpha_vdim[1]*f[23]+alpha_vdim[0]*f[13]); 
  out[25] += 0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.9682458365518543*(alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]); 
  out[28] += 0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += 0.4330127018922193*(alpha_vdim[1]*f[28]+alpha_vdim[0]*f[14]); 
  out[31] += 0.38729833462074165*(alpha_cdim[2]*f[38]+alpha_vdim[1]*f[37])+0.4330127018922193*(alpha_cdim[0]*f[18]+alpha_vdim[0]*f[17]+(alpha_cdim[2]+alpha_vdim[1])*f[10]); 
  out[32] += 0.8660254037844386*alpha_cdim[2]*f[33]+0.4330127018922193*alpha_vdim[0]*f[21]+0.9682458365518543*alpha_cdim[0]*f[15]+(0.9682458365518543*alpha_cdim[2]+0.38729833462074165*alpha_vdim[1])*f[6]; 
  out[33] += 0.8660254037844386*alpha_vdim[1]*f[32]+0.4330127018922193*alpha_cdim[0]*f[22]+0.9682458365518543*alpha_vdim[0]*f[15]+(0.38729833462074165*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[7]; 
  out[34] += 0.4330127018922193*(alpha_cdim[0]*f[24]+alpha_vdim[0]*f[23]+(alpha_cdim[2]+alpha_vdim[1])*f[13]); 
  out[35] += 0.8660254037844386*alpha_cdim[2]*f[36]+0.4330127018922193*alpha_vdim[0]*f[25]+0.9682458365518543*alpha_cdim[0]*f[16]+(0.9682458365518543*alpha_cdim[2]+0.38729833462074165*alpha_vdim[1])*f[8]; 
  out[36] += 0.8660254037844386*alpha_vdim[1]*f[35]+0.4330127018922193*alpha_cdim[0]*f[26]+0.9682458365518543*alpha_vdim[0]*f[16]+(0.38729833462074165*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[9]; 
  out[37] += 0.9682458365518543*(alpha_cdim[2]*f[31]+alpha_cdim[0]*f[17]); 
  out[38] += 0.9682458365518543*(alpha_vdim[1]*f[31]+alpha_vdim[0]*f[18]); 
  out[39] += 0.4330127018922193*(alpha_cdim[2]*f[40]+alpha_cdim[0]*f[27]); 
  out[40] += 0.4330127018922193*(alpha_vdim[1]*f[39]+alpha_vdim[0]*f[27]); 
  out[41] += 0.4330127018922193*(alpha_cdim[0]*f[29]+alpha_vdim[0]*f[28]+(alpha_cdim[2]+alpha_vdim[1])*f[14]); 
  out[42] += 0.4330127018922193*(alpha_cdim[2]*f[43]+alpha_cdim[0]*f[30]); 
  out[43] += 0.4330127018922193*(alpha_vdim[1]*f[42]+alpha_vdim[0]*f[30]); 
  out[44] += 0.8660254037844386*alpha_cdim[2]*f[45]+0.4330127018922193*alpha_vdim[0]*f[37]+0.9682458365518543*alpha_cdim[0]*f[31]+(0.9682458365518543*alpha_cdim[2]+0.38729833462074165*alpha_vdim[1])*f[17]; 
  out[45] += 0.8660254037844386*alpha_vdim[1]*f[44]+0.4330127018922193*alpha_cdim[0]*f[38]+0.9682458365518543*alpha_vdim[0]*f[31]+(0.38729833462074165*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[18]; 
  out[46] += 0.4330127018922193*(alpha_cdim[0]*f[40]+alpha_vdim[0]*f[39]+(alpha_cdim[2]+alpha_vdim[1])*f[27]); 
  out[47] += 0.4330127018922193*(alpha_cdim[0]*f[43]+alpha_vdim[0]*f[42]+(alpha_cdim[2]+alpha_vdim[1])*f[30]); 

  return cflFreq_mid; 
} 

