#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // pots:      potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:     external E and B fields (scaled by q/m).
  // f:         Input distribution function.
  // out:       Incremented output.

  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *phi = &pots[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[40] = {0.0}; 
  double alpha_vdim[120] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.3094010767585034*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[2]; 
  const double *Ez = &EBext[4]; 
  const double *Bx = &EBext[6]; 
  const double *By = &EBext[8]; 
  const double *Bz = &EBext[10]; 

  alpha_vdim[0] = dv10*(-(2.8284271247461907*By[0]*wv3)+2.8284271247461907*Bz[0]*wv2-4.898979485566357*phi[1]*dx10+2.8284271247461907*Ex[0]); 
  alpha_vdim[1] = dv10*(2.8284271247461907*(Bz[1]*wv2+Ex[1])-2.8284271247461907*By[1]*wv3); 
  alpha_vdim[3] = 0.8164965809277261*Bz[0]*dv10*dv2; 
  alpha_vdim[4] = -(0.8164965809277261*By[0]*dv10*dv3); 
  alpha_vdim[6] = 0.8164965809277261*Bz[1]*dv10*dv2; 
  alpha_vdim[8] = -(0.8164965809277261*By[1]*dv10*dv3); 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[40] = dv11*(2.8284271247461907*Bx[0]*wv3-2.8284271247461907*Bz[0]*wv1+2.8284271247461907*Ey[0]); 
  alpha_vdim[41] = dv11*(2.8284271247461907*Bx[1]*wv3-2.8284271247461907*Bz[1]*wv1+2.8284271247461907*Ey[1]); 
  alpha_vdim[42] = -(0.8164965809277261*Bz[0]*dv1*dv11); 
  alpha_vdim[44] = 0.8164965809277261*Bx[0]*dv11*dv3; 
  alpha_vdim[45] = -(0.8164965809277261*Bz[1]*dv1*dv11); 
  alpha_vdim[48] = 0.8164965809277261*Bx[1]*dv11*dv3; 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[40]); 

  alpha_vdim[80] = dv12*(2.8284271247461907*(By[0]*wv1+Ez[0])-2.8284271247461907*Bx[0]*wv2); 
  alpha_vdim[81] = dv12*(2.8284271247461907*(By[1]*wv1+Ez[1])-2.8284271247461907*Bx[1]*wv2); 
  alpha_vdim[82] = 0.8164965809277261*By[0]*dv1*dv12; 
  alpha_vdim[83] = -(0.8164965809277261*Bx[0]*dv12*dv2); 
  alpha_vdim[85] = 0.8164965809277261*By[1]*dv1*dv12; 
  alpha_vdim[86] = -(0.8164965809277261*Bx[1]*dv12*dv2); 
  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[80]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[8]*alpha_vdim[48]+f[5]*alpha_vdim[45]+f[4]*alpha_vdim[44]+f[2]*alpha_vdim[42]+f[1]*alpha_vdim[41]+f[0]*alpha_vdim[40]); 
  out[4] += 0.4330127018922193*(f[6]*alpha_vdim[86]+f[5]*alpha_vdim[85]+f[3]*alpha_vdim[83]+f[2]*alpha_vdim[82]+f[1]*alpha_vdim[81]+f[0]*alpha_vdim[80]); 
  out[5] += 0.38729833462074165*alpha_cdim[2]*f[16]+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(f[4]*alpha_vdim[48]+f[2]*alpha_vdim[45]+f[8]*alpha_vdim[44]+f[5]*alpha_vdim[42]+f[0]*alpha_vdim[41]+f[1]*alpha_vdim[40]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*f[12]*alpha_vdim[48]+0.38729833462074165*f[17]*alpha_vdim[45]+0.4330127018922193*(f[1]*alpha_vdim[45]+f[9]*alpha_vdim[44])+0.38729833462074165*f[16]*alpha_vdim[42]+0.4330127018922193*(f[0]*alpha_vdim[42]+f[5]*alpha_vdim[41]+f[2]*alpha_vdim[40])+0.38729833462074165*(alpha_vdim[6]*f[25]+alpha_vdim[3]*f[24])+0.4330127018922193*(alpha_vdim[8]*f[13]+alpha_vdim[4]*f[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.4330127018922193*(f[3]*alpha_vdim[86]+f[2]*alpha_vdim[85]+f[6]*alpha_vdim[83]+f[5]*alpha_vdim[82]+f[0]*alpha_vdim[81]+f[1]*alpha_vdim[80]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*f[11]*alpha_vdim[86]+0.38729833462074165*f[17]*alpha_vdim[85]+0.4330127018922193*(f[1]*alpha_vdim[85]+f[7]*alpha_vdim[83])+0.38729833462074165*f[16]*alpha_vdim[82]+0.4330127018922193*(f[0]*alpha_vdim[82]+f[5]*alpha_vdim[81]+f[2]*alpha_vdim[80])+0.38729833462074165*(alpha_vdim[8]*f[33]+alpha_vdim[4]*f[32])+0.4330127018922193*(alpha_vdim[6]*f[13]+alpha_vdim[3]*f[10]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.38729833462074165*f[25]*alpha_vdim[86]+0.4330127018922193*(f[1]*alpha_vdim[86]+f[11]*alpha_vdim[85])+0.38729833462074165*f[24]*alpha_vdim[83]+0.4330127018922193*(f[0]*alpha_vdim[83]+f[7]*alpha_vdim[82]+f[6]*alpha_vdim[81]+f[3]*alpha_vdim[80])+0.38729833462074165*f[33]*alpha_vdim[48]+0.4330127018922193*(f[1]*alpha_vdim[48]+f[12]*alpha_vdim[45])+0.38729833462074165*f[32]*alpha_vdim[44]+0.4330127018922193*(f[0]*alpha_vdim[44]+f[9]*alpha_vdim[42]+f[8]*alpha_vdim[41]+f[4]*alpha_vdim[40]); 
  out[11] += 0.4330127018922193*f[9]*alpha_vdim[48]+0.38729833462074165*f[16]*alpha_vdim[45]+0.4330127018922193*(f[0]*alpha_vdim[45]+f[12]*alpha_vdim[44])+0.38729833462074165*f[17]*alpha_vdim[42]+0.4330127018922193*(f[1]*alpha_vdim[42]+f[2]*alpha_vdim[41]+f[5]*alpha_vdim[40])+0.38729833462074165*(alpha_vdim[3]*f[25]+alpha_vdim[6]*f[24]+alpha_cdim[2]*f[18])+0.4330127018922193*(alpha_vdim[4]*f[13]+alpha_vdim[8]*f[10]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[12] += 0.4330127018922193*f[7]*alpha_vdim[86]+0.38729833462074165*f[16]*alpha_vdim[85]+0.4330127018922193*(f[0]*alpha_vdim[85]+f[11]*alpha_vdim[83])+0.38729833462074165*f[17]*alpha_vdim[82]+0.4330127018922193*(f[1]*alpha_vdim[82]+f[2]*alpha_vdim[81]+f[5]*alpha_vdim[80])+0.38729833462074165*(alpha_vdim[4]*f[33]+alpha_vdim[8]*f[32]+alpha_cdim[2]*f[19])+0.4330127018922193*(alpha_vdim[3]*f[13]+alpha_vdim[6]*f[10]+alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[13] += 0.38729833462074165*f[24]*alpha_vdim[86]+0.4330127018922193*(f[0]*alpha_vdim[86]+f[7]*alpha_vdim[85])+0.38729833462074165*f[25]*alpha_vdim[83]+0.4330127018922193*(f[1]*alpha_vdim[83]+f[11]*alpha_vdim[82]+f[3]*alpha_vdim[81]+f[6]*alpha_vdim[80])+0.38729833462074165*f[32]*alpha_vdim[48]+0.4330127018922193*(f[0]*alpha_vdim[48]+f[9]*alpha_vdim[45])+0.38729833462074165*f[33]*alpha_vdim[44]+0.4330127018922193*(f[1]*alpha_vdim[44]+f[12]*alpha_vdim[42]+f[4]*alpha_vdim[41]+f[8]*alpha_vdim[40]+alpha_cdim[2]*f[14]+alpha_cdim[0]*f[10]); 
  out[14] += (0.38729833462074165*f[28]+0.4330127018922193*f[5])*alpha_vdim[86]+(0.38729833462074165*f[20]+0.4330127018922193*f[6])*alpha_vdim[85]+(0.38729833462074165*f[26]+0.4330127018922193*f[2])*alpha_vdim[83]+0.38729833462074165*f[18]*alpha_vdim[82]+0.4330127018922193*(f[3]*alpha_vdim[82]+f[11]*alpha_vdim[81]+f[7]*alpha_vdim[80])+(0.38729833462074165*f[36]+0.4330127018922193*f[5])*alpha_vdim[48]+(0.38729833462074165*f[21]+0.4330127018922193*f[8])*alpha_vdim[45]+(0.38729833462074165*f[34]+0.4330127018922193*f[2])*alpha_vdim[44]+0.38729833462074165*f[19]*alpha_vdim[42]+0.4330127018922193*(f[4]*alpha_vdim[42]+f[12]*alpha_vdim[41]+f[9]*alpha_vdim[40])+0.38729833462074165*(alpha_vdim[8]*f[37]+alpha_vdim[4]*f[35]+alpha_vdim[6]*f[29]+alpha_vdim[3]*f[27])+0.4330127018922193*(alpha_vdim[1]*f[13]+alpha_vdim[0]*f[10]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[15] += (0.38729833462074165*f[26]+0.4330127018922193*f[2])*alpha_vdim[86]+(0.38729833462074165*f[18]+0.4330127018922193*f[3])*alpha_vdim[85]+(0.38729833462074165*f[28]+0.4330127018922193*f[5])*alpha_vdim[83]+0.38729833462074165*f[20]*alpha_vdim[82]+0.4330127018922193*(f[6]*alpha_vdim[82]+f[7]*alpha_vdim[81]+f[11]*alpha_vdim[80])+(0.38729833462074165*f[34]+0.4330127018922193*f[2])*alpha_vdim[48]+(0.38729833462074165*f[19]+0.4330127018922193*f[4])*alpha_vdim[45]+(0.38729833462074165*f[36]+0.4330127018922193*f[5])*alpha_vdim[44]+0.38729833462074165*f[21]*alpha_vdim[42]+0.4330127018922193*(f[8]*alpha_vdim[42]+f[9]*alpha_vdim[41]+f[12]*alpha_vdim[40])+0.38729833462074165*(alpha_vdim[4]*f[37]+alpha_vdim[8]*f[35]+alpha_vdim[3]*f[29]+alpha_vdim[6]*f[27]+alpha_cdim[2]*f[22])+0.4330127018922193*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+(alpha_cdim[2]+alpha_vdim[1])*f[10]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]); 
  out[16] += 0.9682458365518543*(alpha_vdim[8]*f[12]+alpha_vdim[6]*f[11]+alpha_vdim[4]*f[9]+alpha_vdim[3]*f[7]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[17] += 0.4330127018922193*alpha_cdim[0]*f[16]+0.9682458365518543*(alpha_vdim[4]*f[12]+alpha_vdim[3]*f[11]+alpha_vdim[8]*f[9]+alpha_vdim[6]*f[7]+alpha_vdim[0]*f[5])+(0.38729833462074165*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[18] += 0.4330127018922193*f[21]*alpha_vdim[48]+0.38729833462074165*f[5]*alpha_vdim[45]+0.4330127018922193*f[19]*alpha_vdim[44]+0.38729833462074165*f[2]*alpha_vdim[42]+0.4330127018922193*(f[17]*alpha_vdim[41]+f[16]*alpha_vdim[40])+0.8660254037844386*(alpha_vdim[6]*f[28]+alpha_vdim[3]*f[26])+0.9682458365518543*(alpha_vdim[8]*f[15]+alpha_vdim[4]*f[14]+alpha_vdim[1]*f[11]+alpha_vdim[0]*f[7]+f[5]*alpha_vdim[6]+f[2]*alpha_vdim[3]); 
  out[19] += 0.4330127018922193*f[20]*alpha_vdim[86]+0.38729833462074165*f[5]*alpha_vdim[85]+0.4330127018922193*f[18]*alpha_vdim[83]+0.38729833462074165*f[2]*alpha_vdim[82]+0.4330127018922193*(f[17]*alpha_vdim[81]+f[16]*alpha_vdim[80])+0.8660254037844386*(alpha_vdim[8]*f[36]+alpha_vdim[4]*f[34])+0.9682458365518543*(alpha_vdim[6]*f[15]+alpha_vdim[3]*f[14]+alpha_vdim[1]*f[12]+alpha_vdim[0]*f[9]+f[5]*alpha_vdim[8]+f[2]*alpha_vdim[4]); 
  out[20] += 0.4330127018922193*f[19]*alpha_vdim[48]+0.38729833462074165*f[2]*alpha_vdim[45]+0.4330127018922193*f[21]*alpha_vdim[44]+0.38729833462074165*f[5]*alpha_vdim[42]+0.4330127018922193*(f[16]*alpha_vdim[41]+f[17]*alpha_vdim[40])+0.8660254037844386*(alpha_vdim[3]*f[28]+alpha_vdim[6]*f[26])+0.4330127018922193*alpha_cdim[0]*f[18]+0.9682458365518543*(alpha_vdim[4]*f[15]+alpha_vdim[8]*f[14]+alpha_vdim[0]*f[11])+0.38729833462074165*alpha_cdim[2]*f[7]+0.9682458365518543*(alpha_vdim[1]*f[7]+f[2]*alpha_vdim[6]+alpha_vdim[3]*f[5]); 
  out[21] += 0.4330127018922193*f[18]*alpha_vdim[86]+0.38729833462074165*f[2]*alpha_vdim[85]+0.4330127018922193*f[20]*alpha_vdim[83]+0.38729833462074165*f[5]*alpha_vdim[82]+0.4330127018922193*(f[16]*alpha_vdim[81]+f[17]*alpha_vdim[80])+0.8660254037844386*(alpha_vdim[4]*f[36]+alpha_vdim[8]*f[34])+0.4330127018922193*alpha_cdim[0]*f[19]+0.9682458365518543*(alpha_vdim[3]*f[15]+alpha_vdim[6]*f[14]+alpha_vdim[0]*f[12])+0.38729833462074165*alpha_cdim[2]*f[9]+0.9682458365518543*(alpha_vdim[1]*f[9]+f[2]*alpha_vdim[8]+alpha_vdim[4]*f[5]); 
  out[22] += 0.4330127018922193*f[17]*alpha_vdim[86]+0.38729833462074165*f[11]*alpha_vdim[85]+0.4330127018922193*f[16]*alpha_vdim[83]+0.38729833462074165*f[7]*alpha_vdim[82]+0.4330127018922193*(f[20]*alpha_vdim[81]+f[18]*alpha_vdim[80]+f[17]*alpha_vdim[48])+0.38729833462074165*f[12]*alpha_vdim[45]+0.4330127018922193*f[16]*alpha_vdim[44]+0.38729833462074165*f[9]*alpha_vdim[42]+0.4330127018922193*(f[21]*alpha_vdim[41]+f[19]*alpha_vdim[40])+0.8660254037844386*(alpha_vdim[8]*f[39]+alpha_vdim[4]*f[38]+alpha_vdim[6]*f[31]+alpha_vdim[3]*f[30])+0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[14]+alpha_vdim[6]*f[12]+alpha_vdim[8]*f[11]+alpha_vdim[3]*f[9]+alpha_vdim[4]*f[7]); 
  out[23] += 0.4330127018922193*f[16]*alpha_vdim[86]+0.38729833462074165*f[7]*alpha_vdim[85]+0.4330127018922193*f[17]*alpha_vdim[83]+0.38729833462074165*f[11]*alpha_vdim[82]+0.4330127018922193*(f[18]*alpha_vdim[81]+f[20]*alpha_vdim[80]+f[16]*alpha_vdim[48])+0.38729833462074165*f[9]*alpha_vdim[45]+0.4330127018922193*f[17]*alpha_vdim[44]+0.38729833462074165*f[12]*alpha_vdim[42]+0.4330127018922193*(f[19]*alpha_vdim[41]+f[21]*alpha_vdim[40])+0.8660254037844386*(alpha_vdim[4]*f[39]+alpha_vdim[8]*f[38]+alpha_vdim[3]*f[31]+alpha_vdim[6]*f[30])+0.4330127018922193*alpha_cdim[0]*f[22]+0.9682458365518543*alpha_vdim[0]*f[15]+0.38729833462074165*alpha_cdim[2]*f[14]+0.9682458365518543*(alpha_vdim[1]*f[14]+alpha_vdim[3]*f[12]+alpha_vdim[4]*f[11]+alpha_vdim[6]*f[9]+f[7]*alpha_vdim[8]); 
  out[24] += 0.9682458365518543*(f[13]*alpha_vdim[48]+f[11]*alpha_vdim[45]+f[10]*alpha_vdim[44]+f[7]*alpha_vdim[42]+f[6]*alpha_vdim[41]+f[3]*alpha_vdim[40]); 
  out[25] += 0.9682458365518543*(f[10]*alpha_vdim[48]+f[7]*alpha_vdim[45]+f[13]*alpha_vdim[44]+f[11]*alpha_vdim[42]+f[3]*alpha_vdim[41]+f[6]*alpha_vdim[40])+0.4330127018922193*(alpha_cdim[2]*f[26]+alpha_cdim[0]*f[24]); 
  out[26] += 0.9682458365518543*f[15]*alpha_vdim[48]+0.8660254037844386*f[20]*alpha_vdim[45]+0.9682458365518543*(f[6]*alpha_vdim[45]+f[14]*alpha_vdim[44])+0.8660254037844386*f[18]*alpha_vdim[42]+0.9682458365518543*(f[3]*alpha_vdim[42]+f[11]*alpha_vdim[41]+f[7]*alpha_vdim[40])+0.4330127018922193*(alpha_vdim[8]*f[29]+alpha_vdim[4]*f[27]+alpha_vdim[1]*f[25]+alpha_vdim[0]*f[24])+0.38729833462074165*(alpha_vdim[6]*f[6]+alpha_vdim[3]*f[3]); 
  out[27] += 0.38729833462074165*f[6]*alpha_vdim[86]+0.4330127018922193*f[28]*alpha_vdim[85]+0.38729833462074165*f[3]*alpha_vdim[83]+0.4330127018922193*(f[26]*alpha_vdim[82]+f[25]*alpha_vdim[81]+f[24]*alpha_vdim[80])+0.8660254037844386*f[37]*alpha_vdim[48]+0.9682458365518543*(f[6]*alpha_vdim[48]+f[15]*alpha_vdim[45])+0.8660254037844386*f[35]*alpha_vdim[44]+0.9682458365518543*(f[3]*alpha_vdim[44]+f[14]*alpha_vdim[42]+f[13]*alpha_vdim[41]+f[10]*alpha_vdim[40]); 
  out[28] += 0.9682458365518543*f[14]*alpha_vdim[48]+0.8660254037844386*f[18]*alpha_vdim[45]+0.9682458365518543*(f[3]*alpha_vdim[45]+f[15]*alpha_vdim[44])+0.8660254037844386*f[20]*alpha_vdim[42]+0.9682458365518543*(f[6]*alpha_vdim[42]+f[7]*alpha_vdim[41]+f[11]*alpha_vdim[40])+0.4330127018922193*(alpha_vdim[4]*f[29]+alpha_vdim[8]*f[27]+alpha_cdim[0]*f[26]+alpha_vdim[0]*f[25]+(alpha_cdim[2]+alpha_vdim[1])*f[24])+0.38729833462074165*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]); 
  out[29] += 0.38729833462074165*f[3]*alpha_vdim[86]+0.4330127018922193*f[26]*alpha_vdim[85]+0.38729833462074165*f[6]*alpha_vdim[83]+0.4330127018922193*(f[28]*alpha_vdim[82]+f[24]*alpha_vdim[81]+f[25]*alpha_vdim[80])+0.8660254037844386*f[35]*alpha_vdim[48]+0.9682458365518543*(f[3]*alpha_vdim[48]+f[14]*alpha_vdim[45])+0.8660254037844386*f[37]*alpha_vdim[44]+0.9682458365518543*(f[6]*alpha_vdim[44]+f[15]*alpha_vdim[42]+f[10]*alpha_vdim[41]+f[13]*alpha_vdim[40])+0.4330127018922193*(alpha_cdim[2]*f[30]+alpha_cdim[0]*f[27]); 
  out[30] += 0.38729833462074165*f[11]*alpha_vdim[86]+0.4330127018922193*f[25]*alpha_vdim[85]+0.38729833462074165*f[7]*alpha_vdim[83]+0.4330127018922193*(f[24]*alpha_vdim[82]+f[28]*alpha_vdim[81]+f[26]*alpha_vdim[80])+(0.8660254037844386*f[39]+0.9682458365518543*f[11])*alpha_vdim[48]+(0.8660254037844386*f[23]+0.9682458365518543*f[13])*alpha_vdim[45]+(0.8660254037844386*f[38]+0.9682458365518543*f[7])*alpha_vdim[44]+0.8660254037844386*f[22]*alpha_vdim[42]+0.9682458365518543*(f[10]*alpha_vdim[42]+f[15]*alpha_vdim[41]+f[14]*alpha_vdim[40])+0.4330127018922193*(alpha_vdim[1]*f[29]+alpha_vdim[0]*f[27]+alpha_vdim[8]*f[25]+alpha_vdim[4]*f[24])+0.38729833462074165*(alpha_vdim[6]*f[13]+alpha_vdim[3]*f[10]); 
  out[31] += 0.38729833462074165*f[7]*alpha_vdim[86]+0.4330127018922193*f[24]*alpha_vdim[85]+0.38729833462074165*f[11]*alpha_vdim[83]+0.4330127018922193*(f[25]*alpha_vdim[82]+f[26]*alpha_vdim[81]+f[28]*alpha_vdim[80])+(0.8660254037844386*f[38]+0.9682458365518543*f[7])*alpha_vdim[48]+(0.8660254037844386*f[22]+0.9682458365518543*f[10])*alpha_vdim[45]+(0.8660254037844386*f[39]+0.9682458365518543*f[11])*alpha_vdim[44]+0.8660254037844386*f[23]*alpha_vdim[42]+0.9682458365518543*(f[13]*alpha_vdim[42]+f[14]*alpha_vdim[41]+f[15]*alpha_vdim[40])+0.4330127018922193*(alpha_cdim[0]*f[30]+alpha_vdim[0]*f[29]+(alpha_cdim[2]+alpha_vdim[1])*f[27]+alpha_vdim[4]*f[25]+alpha_vdim[8]*f[24])+0.38729833462074165*(alpha_vdim[3]*f[13]+alpha_vdim[6]*f[10]); 
  out[32] += 0.9682458365518543*(f[13]*alpha_vdim[86]+f[12]*alpha_vdim[85]+f[10]*alpha_vdim[83]+f[9]*alpha_vdim[82]+f[8]*alpha_vdim[81]+f[4]*alpha_vdim[80]); 
  out[33] += 0.9682458365518543*(f[10]*alpha_vdim[86]+f[9]*alpha_vdim[85]+f[13]*alpha_vdim[83]+f[12]*alpha_vdim[82]+f[4]*alpha_vdim[81]+f[8]*alpha_vdim[80])+0.4330127018922193*(alpha_cdim[2]*f[34]+alpha_cdim[0]*f[32]); 
  out[34] += 0.9682458365518543*f[15]*alpha_vdim[86]+0.8660254037844386*f[21]*alpha_vdim[85]+0.9682458365518543*(f[8]*alpha_vdim[85]+f[14]*alpha_vdim[83])+0.8660254037844386*f[19]*alpha_vdim[82]+0.9682458365518543*(f[4]*alpha_vdim[82]+f[12]*alpha_vdim[81]+f[9]*alpha_vdim[80])+0.4330127018922193*(alpha_vdim[6]*f[37]+alpha_vdim[3]*f[35]+alpha_vdim[1]*f[33]+alpha_vdim[0]*f[32])+0.38729833462074165*(alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[35] += 0.8660254037844386*f[29]*alpha_vdim[86]+0.9682458365518543*(f[8]*alpha_vdim[86]+f[15]*alpha_vdim[85])+0.8660254037844386*f[27]*alpha_vdim[83]+0.9682458365518543*(f[4]*alpha_vdim[83]+f[14]*alpha_vdim[82]+f[13]*alpha_vdim[81]+f[10]*alpha_vdim[80])+0.38729833462074165*f[8]*alpha_vdim[48]+0.4330127018922193*f[36]*alpha_vdim[45]+0.38729833462074165*f[4]*alpha_vdim[44]+0.4330127018922193*(f[34]*alpha_vdim[42]+f[33]*alpha_vdim[41]+f[32]*alpha_vdim[40]); 
  out[36] += 0.9682458365518543*f[14]*alpha_vdim[86]+0.8660254037844386*f[19]*alpha_vdim[85]+0.9682458365518543*(f[4]*alpha_vdim[85]+f[15]*alpha_vdim[83])+0.8660254037844386*f[21]*alpha_vdim[82]+0.9682458365518543*(f[8]*alpha_vdim[82]+f[9]*alpha_vdim[81]+f[12]*alpha_vdim[80])+0.4330127018922193*(alpha_vdim[3]*f[37]+alpha_vdim[6]*f[35]+alpha_cdim[0]*f[34]+alpha_vdim[0]*f[33]+(alpha_cdim[2]+alpha_vdim[1])*f[32])+0.38729833462074165*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]); 
  out[37] += 0.8660254037844386*f[27]*alpha_vdim[86]+0.9682458365518543*(f[4]*alpha_vdim[86]+f[14]*alpha_vdim[85])+0.8660254037844386*f[29]*alpha_vdim[83]+0.9682458365518543*(f[8]*alpha_vdim[83]+f[15]*alpha_vdim[82]+f[10]*alpha_vdim[81]+f[13]*alpha_vdim[80])+0.38729833462074165*f[4]*alpha_vdim[48]+0.4330127018922193*f[34]*alpha_vdim[45]+0.38729833462074165*f[8]*alpha_vdim[44]+0.4330127018922193*(f[36]*alpha_vdim[42]+f[32]*alpha_vdim[41]+f[33]*alpha_vdim[40]+alpha_cdim[2]*f[38]+alpha_cdim[0]*f[35]); 
  out[38] += (0.8660254037844386*f[31]+0.9682458365518543*f[12])*alpha_vdim[86]+(0.8660254037844386*f[23]+0.9682458365518543*f[13])*alpha_vdim[85]+(0.8660254037844386*f[30]+0.9682458365518543*f[9])*alpha_vdim[83]+0.8660254037844386*f[22]*alpha_vdim[82]+0.9682458365518543*(f[10]*alpha_vdim[82]+f[15]*alpha_vdim[81]+f[14]*alpha_vdim[80])+0.38729833462074165*f[12]*alpha_vdim[48]+0.4330127018922193*f[33]*alpha_vdim[45]+0.38729833462074165*f[9]*alpha_vdim[44]+0.4330127018922193*(f[32]*alpha_vdim[42]+f[36]*alpha_vdim[41]+f[34]*alpha_vdim[40]+alpha_vdim[1]*f[37]+alpha_vdim[0]*f[35]+alpha_vdim[6]*f[33]+alpha_vdim[3]*f[32])+0.38729833462074165*(alpha_vdim[8]*f[13]+alpha_vdim[4]*f[10]); 
  out[39] += (0.8660254037844386*f[30]+0.9682458365518543*f[9])*alpha_vdim[86]+(0.8660254037844386*f[22]+0.9682458365518543*f[10])*alpha_vdim[85]+(0.8660254037844386*f[31]+0.9682458365518543*f[12])*alpha_vdim[83]+0.8660254037844386*f[23]*alpha_vdim[82]+0.9682458365518543*(f[13]*alpha_vdim[82]+f[14]*alpha_vdim[81]+f[15]*alpha_vdim[80])+0.38729833462074165*f[9]*alpha_vdim[48]+0.4330127018922193*f[32]*alpha_vdim[45]+0.38729833462074165*f[12]*alpha_vdim[44]+0.4330127018922193*(f[33]*alpha_vdim[42]+f[34]*alpha_vdim[41]+f[36]*alpha_vdim[40]+alpha_cdim[0]*f[38]+alpha_vdim[0]*f[37]+(alpha_cdim[2]+alpha_vdim[1])*f[35]+alpha_vdim[3]*f[33]+alpha_vdim[6]*f[32])+0.38729833462074165*(alpha_vdim[4]*f[13]+alpha_vdim[8]*f[10]); 

  return cflFreq_mid; 
} 

