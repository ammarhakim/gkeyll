#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_1x2v_tensor_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
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

  const double *B2 = &field[10]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[8] = {0.0}; 
  double alpha_vdim[16] = {0.0}; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[8] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[9] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[10] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[12] = -0.5773502691896258*B2[1]*dv1*dv11; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_vdim[8]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[4]*alpha_vdim[12]+f[2]*alpha_vdim[10]+f[1]*alpha_vdim[9]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[2]*alpha_vdim[12]+f[4]*alpha_vdim[10]+f[0]*alpha_vdim[9]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[1]*alpha_vdim[12]+f[0]*alpha_vdim[10]+f[4]*alpha_vdim[9]+f[2]*alpha_vdim[8]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[0]*alpha_vdim[12]+f[1]*alpha_vdim[10]+f[2]*alpha_vdim[9]+f[4]*alpha_vdim[8]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 

  return cflFreq_mid; 
} 
