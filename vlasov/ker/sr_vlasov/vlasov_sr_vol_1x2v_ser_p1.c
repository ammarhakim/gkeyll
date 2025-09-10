#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 
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
  double alpha_cdim[16] = {0.0}; 
  double alpha_vdim[32] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[6] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[12] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_cdim[0]-0.1976423537605236*alpha_cdim[12]); 

  alpha_vdim[0] = (B2[0]*p1_over_gamma[0]+2.0*E0[0])*dv10; 
  alpha_vdim[1] = (2.0*E0[1]+p1_over_gamma[0]*B2[1])*dv10; 
  alpha_vdim[2] = B2[0]*p1_over_gamma[1]*dv10; 
  alpha_vdim[3] = B2[0]*p1_over_gamma[2]*dv10; 
  alpha_vdim[4] = B2[1]*p1_over_gamma[1]*dv10; 
  alpha_vdim[5] = B2[1]*p1_over_gamma[2]*dv10; 
  alpha_vdim[6] = B2[0]*p1_over_gamma[3]*dv10; 
  alpha_vdim[7] = B2[1]*p1_over_gamma[3]*dv10; 
  alpha_vdim[8] = B2[0]*p1_over_gamma[4]*dv10; 
  alpha_vdim[9] = B2[1]*p1_over_gamma[4]*dv10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[8]); 

  alpha_vdim[16] = (2.0*E1[0]-1.0*B2[0]*p0_over_gamma[0])*dv11; 
  alpha_vdim[17] = (2.0*E1[1]-1.0*p0_over_gamma[0]*B2[1])*dv11; 
  alpha_vdim[18] = -1.0*B2[0]*p0_over_gamma[1]*dv11; 
  alpha_vdim[19] = -1.0*B2[0]*p0_over_gamma[2]*dv11; 
  alpha_vdim[20] = -1.0*B2[1]*p0_over_gamma[1]*dv11; 
  alpha_vdim[21] = -1.0*B2[1]*p0_over_gamma[2]*dv11; 
  alpha_vdim[22] = -1.0*B2[0]*p0_over_gamma[3]*dv11; 
  alpha_vdim[23] = -1.0*B2[1]*p0_over_gamma[3]*dv11; 
  alpha_vdim[28] = -1.0*B2[0]*p0_over_gamma[5]*dv11; 
  alpha_vdim[29] = -1.0*B2[1]*p0_over_gamma[5]*dv11; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[16]-0.1976423537605236*alpha_vdim[28]); 

  out[1] += 0.6123724356957944*(alpha_cdim[12]*f[12]+alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[13]*alpha_vdim[29]+f[12]*alpha_vdim[28]+f[7]*alpha_vdim[23]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[4]*alpha_vdim[20]+f[3]*alpha_vdim[19]+f[2]*alpha_vdim[18]+f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[4] += 0.6123724356957944*alpha_cdim[12]*f[14]+0.5477225575051661*alpha_cdim[6]*f[10]+0.6123724356957944*alpha_vdim[8]*f[9]+f[8]*(0.6123724356957944*alpha_vdim[9]+0.5477225575051661*alpha_cdim[2])+0.6123724356957944*(alpha_vdim[6]*f[7]+f[6]*(alpha_vdim[7]+alpha_cdim[3])+f[3]*alpha_cdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*(alpha_vdim[4]+alpha_cdim[0])+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[12]*alpha_vdim[29]+f[13]*alpha_vdim[28]+f[6]*alpha_vdim[23]+f[7]*alpha_vdim[22]+f[3]*alpha_vdim[21]+f[2]*alpha_vdim[20]+f[5]*alpha_vdim[19]+f[4]*alpha_vdim[18]+f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16])+0.5477225575051661*(alpha_cdim[6]*f[14]+alpha_cdim[3]*f[12]+f[3]*alpha_cdim[12])+0.6123724356957944*(alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[6] += 0.6123724356957944*(f[15]*alpha_vdim[29]+f[14]*alpha_vdim[28])+(0.5477225575051661*f[11]+0.6123724356957944*f[5])*alpha_vdim[23]+0.5477225575051661*f[10]*alpha_vdim[22]+0.6123724356957944*(f[3]*alpha_vdim[22]+f[7]*alpha_vdim[21])+0.5477225575051661*f[9]*alpha_vdim[20]+0.6123724356957944*(f[1]*alpha_vdim[20]+f[6]*alpha_vdim[19])+0.5477225575051661*f[8]*alpha_vdim[18]+0.6123724356957944*(f[0]*alpha_vdim[18]+f[4]*alpha_vdim[17]+f[2]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[7]*f[15]+alpha_vdim[6]*f[14]+alpha_vdim[5]*f[13]+alpha_vdim[3]*f[12])+0.6123724356957944*(alpha_vdim[9]*f[11]+alpha_vdim[8]*f[10]+alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[14]*alpha_vdim[29]+f[15]*alpha_vdim[28])+(0.5477225575051661*f[10]+0.6123724356957944*f[3])*alpha_vdim[23]+0.5477225575051661*f[11]*alpha_vdim[22]+0.6123724356957944*(f[5]*alpha_vdim[22]+f[6]*alpha_vdim[21])+0.5477225575051661*f[8]*alpha_vdim[20]+0.6123724356957944*(f[0]*alpha_vdim[20]+f[7]*alpha_vdim[19])+0.5477225575051661*f[9]*alpha_vdim[18]+0.6123724356957944*(f[1]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[4]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[6]*f[15]+(alpha_vdim[7]+alpha_cdim[3])*f[14]+alpha_vdim[3]*f[13]+(alpha_cdim[6]+alpha_vdim[5])*f[12]+f[6]*alpha_cdim[12])+0.6123724356957944*(alpha_vdim[8]*f[11]+alpha_vdim[9]*f[10])+0.5477225575051661*(alpha_cdim[2]*f[10]+alpha_cdim[6]*f[8])+0.6123724356957944*(alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+(alpha_vdim[4]+alpha_cdim[0])*f[6]+f[4]*alpha_vdim[6]+f[0]*alpha_cdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]+f[2]*alpha_cdim[3]); 
  out[8] += 1.224744871391589*(alpha_vdim[7]*f[11]+alpha_vdim[6]*f[10]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 1.224744871391589*alpha_vdim[6]*f[11]+(1.224744871391589*alpha_vdim[7]+0.6123724356957944*alpha_cdim[3])*f[10]+1.224744871391589*(alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9])+(1.224744871391589*alpha_vdim[4]+0.6123724356957944*alpha_cdim[0])*f[8]+1.224744871391589*f[4]*alpha_vdim[8]+1.369306393762915*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+0.5477225575051661*alpha_cdim[6]*f[6]+1.369306393762915*(alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+0.5477225575051661*alpha_cdim[2]*f[2]+1.369306393762915*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[10] += 0.5477225575051661*(f[7]*alpha_vdim[23]+f[6]*alpha_vdim[22])+0.6123724356957944*f[11]*alpha_vdim[21]+0.5477225575051661*f[4]*alpha_vdim[20]+0.6123724356957944*f[10]*alpha_vdim[19]+0.5477225575051661*f[2]*alpha_vdim[18]+0.6123724356957944*(f[9]*alpha_vdim[17]+f[8]*alpha_vdim[16])+1.224744871391589*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[14]+alpha_vdim[7]*f[13]+alpha_vdim[6]*f[12]+alpha_vdim[4]*f[11]+alpha_vdim[2]*f[10]+alpha_vdim[7]*f[9]+f[7]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[11] += 0.5477225575051661*(f[6]*alpha_vdim[23]+f[7]*alpha_vdim[22])+0.6123724356957944*f[10]*alpha_vdim[21]+0.5477225575051661*f[2]*alpha_vdim[20]+0.6123724356957944*f[11]*alpha_vdim[19]+0.5477225575051661*f[4]*alpha_vdim[18]+0.6123724356957944*(f[8]*alpha_vdim[17]+f[9]*alpha_vdim[16])+1.224744871391589*alpha_vdim[3]*f[15]+0.4898979485566357*alpha_cdim[6]*f[14]+1.224744871391589*(alpha_vdim[5]*f[14]+alpha_vdim[6]*f[13]+alpha_vdim[7]*f[12])+0.5477225575051661*f[10]*alpha_cdim[12]+1.224744871391589*alpha_vdim[2]*f[11]+(1.224744871391589*alpha_vdim[4]+0.6123724356957944*alpha_cdim[0])*f[10]+1.224744871391589*(alpha_vdim[6]*f[9]+f[6]*alpha_vdim[9])+(1.224744871391589*alpha_vdim[7]+0.6123724356957944*alpha_cdim[3])*f[8]+1.224744871391589*f[7]*alpha_vdim[8]+1.369306393762915*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+0.5477225575051661*f[2]*alpha_cdim[6]+1.369306393762915*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[12] += 1.224744871391589*(f[5]*alpha_vdim[29]+f[3]*alpha_vdim[28])+(1.224744871391589*f[15]+1.369306393762915*f[4])*alpha_vdim[23]+(1.224744871391589*f[14]+1.369306393762915*f[2])*alpha_vdim[22]+1.224744871391589*f[13]*alpha_vdim[21]+1.369306393762915*(f[1]*alpha_vdim[21]+f[7]*alpha_vdim[20])+1.224744871391589*f[12]*alpha_vdim[19]+1.369306393762915*(f[0]*alpha_vdim[19]+f[6]*alpha_vdim[18]+f[5]*alpha_vdim[17]+f[3]*alpha_vdim[16]); 
  out[13] += 1.224744871391589*(f[3]*alpha_vdim[29]+f[5]*alpha_vdim[28])+(1.224744871391589*f[14]+1.369306393762915*f[2])*alpha_vdim[23]+(1.224744871391589*f[15]+1.369306393762915*f[4])*alpha_vdim[22]+1.224744871391589*f[12]*alpha_vdim[21]+1.369306393762915*(f[0]*alpha_vdim[21]+f[6]*alpha_vdim[20])+1.224744871391589*f[13]*alpha_vdim[19]+1.369306393762915*(f[1]*alpha_vdim[19]+f[7]*alpha_vdim[18]+f[3]*alpha_vdim[17]+f[5]*alpha_vdim[16])+0.6123724356957944*alpha_cdim[2]*f[14]+0.3912303982179757*alpha_cdim[12]*f[12]+0.6123724356957944*(alpha_cdim[0]*f[12]+f[0]*alpha_cdim[12])+0.5477225575051661*(alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]); 
  out[14] += 1.224744871391589*(f[7]*alpha_vdim[29]+f[6]*alpha_vdim[28])+(1.224744871391589*(f[13]+f[9])+1.369306393762915*f[1])*alpha_vdim[23]+(1.224744871391589*(f[12]+f[8])+1.369306393762915*f[0])*alpha_vdim[22]+(1.224744871391589*f[15]+1.369306393762915*f[4])*alpha_vdim[21]+(1.224744871391589*f[11]+1.369306393762915*f[5])*alpha_vdim[20]+(1.224744871391589*f[14]+1.369306393762915*f[2])*alpha_vdim[19]+1.224744871391589*f[10]*alpha_vdim[18]+1.369306393762915*(f[3]*alpha_vdim[18]+f[7]*alpha_vdim[17]+f[6]*alpha_vdim[16])+0.6123724356957944*(alpha_vdim[4]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[12])+0.5477225575051661*(alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[15] += 1.224744871391589*(f[6]*alpha_vdim[29]+f[7]*alpha_vdim[28])+(1.224744871391589*(f[12]+f[8])+1.369306393762915*f[0])*alpha_vdim[23]+(1.224744871391589*(f[13]+f[9])+1.369306393762915*f[1])*alpha_vdim[22]+(1.224744871391589*f[14]+1.369306393762915*f[2])*alpha_vdim[21]+(1.224744871391589*f[10]+1.369306393762915*f[3])*alpha_vdim[20]+(1.224744871391589*f[15]+1.369306393762915*f[4])*alpha_vdim[19]+1.224744871391589*f[11]*alpha_vdim[18]+1.369306393762915*(f[5]*alpha_vdim[18]+f[6]*alpha_vdim[17]+f[7]*alpha_vdim[16])+0.6123724356957944*alpha_vdim[2]*f[15]+0.3912303982179757*alpha_cdim[12]*f[14]+0.6123724356957944*((alpha_vdim[4]+alpha_cdim[0])*f[14]+alpha_vdim[0]*f[13]+(alpha_cdim[2]+alpha_vdim[1])*f[12]+f[2]*alpha_cdim[12])+0.4898979485566357*alpha_cdim[6]*f[10]+0.5477225575051661*(alpha_vdim[6]*f[7]+f[6]*(alpha_vdim[7]+alpha_cdim[3])+f[3]*alpha_cdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  return cflFreq_mid; 
} 
