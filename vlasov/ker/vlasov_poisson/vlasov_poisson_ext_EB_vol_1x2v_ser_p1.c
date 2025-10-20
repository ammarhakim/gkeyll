#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_ext_EB_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const double *f, double* GKYL_RESTRICT out) 
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

  const double *phi = &pots[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[16] = {0.0}; 
  double alpha_vdim[32] = {0.0}; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.6329931618554527*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[2]; 
  const double *Bx = &EBext[6]; 
  const double *By = &EBext[8]; 
  const double *Bz = &EBext[10]; 

  alpha_vdim[0] = dv10*(2.0*Bz[0]*wv2-3.4641016151377544*phi[1]*dx10+2.0*Ex[0]); 
  alpha_vdim[1] = 2.0*dv10*(Bz[1]*wv2+Ex[1]); 
  alpha_vdim[3] = 0.5773502691896258*Bz[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*Bz[1]*dv10*dv2; 
  cflFreq_mid += 5.0*fabs(0.17677669529663684*alpha_vdim[0]); 

  alpha_vdim[16] = dv11*(2.0*Ey[0]-2.0*Bz[0]*wv1); 
  alpha_vdim[17] = dv11*(2.0*Ey[1]-2.0*Bz[1]*wv1); 
  alpha_vdim[18] = -(0.5773502691896258*Bz[0]*dv1*dv11); 
  alpha_vdim[20] = -(0.5773502691896258*Bz[1]*dv1*dv11); 
  cflFreq_mid += 5.0*fabs(0.17677669529663684*alpha_vdim[16]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[4]*alpha_vdim[20]+f[2]*alpha_vdim[18]+f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[4] += 0.5477225575051661*alpha_cdim[2]*f[8]+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[2]*alpha_vdim[20]+f[4]*alpha_vdim[18]+f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += (0.5477225575051661*f[9]+0.6123724356957944*f[1])*alpha_vdim[20]+0.5477225575051661*f[8]*alpha_vdim[18]+0.6123724356957944*(f[0]*alpha_vdim[18]+f[4]*alpha_vdim[17]+f[2]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[5]*f[13]+alpha_vdim[3]*f[12])+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += (0.5477225575051661*f[8]+0.6123724356957944*f[0])*alpha_vdim[20]+0.5477225575051661*f[9]*alpha_vdim[18]+0.6123724356957944*(f[1]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[4]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[3]*f[13]+alpha_vdim[5]*f[12]+alpha_cdim[2]*f[10])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[8] += 1.369306393762915*(alpha_vdim[5]*f[7]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 0.6123724356957944*alpha_cdim[0]*f[8]+1.369306393762915*(alpha_vdim[3]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[10] += 0.5477225575051661*(f[4]*alpha_vdim[20]+f[2]*alpha_vdim[18])+0.6123724356957944*(f[9]*alpha_vdim[17]+f[8]*alpha_vdim[16])+1.224744871391589*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[14])+1.369306393762915*(alpha_vdim[1]*f[7]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[11] += 0.5477225575051661*(f[2]*alpha_vdim[20]+f[4]*alpha_vdim[18])+0.6123724356957944*(f[8]*alpha_vdim[17]+f[9]*alpha_vdim[16])+1.224744871391589*(alpha_vdim[3]*f[15]+alpha_vdim[5]*f[14])+0.6123724356957944*alpha_cdim[0]*f[10]+1.369306393762915*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[12] += 1.369306393762915*(f[7]*alpha_vdim[20]+f[6]*alpha_vdim[18]+f[5]*alpha_vdim[17]+f[3]*alpha_vdim[16]); 
  out[13] += 1.369306393762915*(f[6]*alpha_vdim[20]+f[7]*alpha_vdim[18]+f[3]*alpha_vdim[17]+f[5]*alpha_vdim[16])+0.6123724356957944*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[12]); 
  out[14] += (1.224744871391589*f[11]+1.369306393762915*f[5])*alpha_vdim[20]+1.224744871391589*f[10]*alpha_vdim[18]+1.369306393762915*(f[3]*alpha_vdim[18]+f[7]*alpha_vdim[17]+f[6]*alpha_vdim[16])+0.6123724356957944*(alpha_vdim[1]*f[13]+alpha_vdim[0]*f[12])+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[15] += (1.224744871391589*f[10]+1.369306393762915*f[3])*alpha_vdim[20]+1.224744871391589*f[11]*alpha_vdim[18]+1.369306393762915*(f[5]*alpha_vdim[18]+f[6]*alpha_vdim[17]+f[7]*alpha_vdim[16])+0.6123724356957944*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+(alpha_cdim[2]+alpha_vdim[1])*f[12])+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  return cflFreq_mid; 
} 

