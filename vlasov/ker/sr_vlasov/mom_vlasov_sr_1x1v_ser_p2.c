#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += (gamma[2]*f[5]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0000000000000002*gamma[2]*f[7]+gamma[1]*f[3]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*gamma[1]*f[6]+gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.4142135623730951*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.4142135623730951*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[2] += volFact*(1.4142135623730951*f[4]*wx1+0.40824829046386296*f[6]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_Ni_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += (p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[4] += (p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
  out[5] += (1.0000000000000002*p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (gamma[2]*f[5]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0000000000000002*gamma[2]*f[7]+gamma[1]*f[3]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*gamma[1]*f[6]+gamma[0]*f[4])*volFact; 
  out[3] += volFact*(1.4142135623730951*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[4] += volFact*(1.4142135623730951*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[5] += volFact*(1.4142135623730951*f[4]*wx1+0.40824829046386296*f[6]*dv1); 
  out[6] += volFact*(p0_over_gamma[1]*f[2]*wx1+f[0]*p0_over_gamma[0]*wx1+0.25819888974716115*p0_over_gamma[1]*f[5]*dv1+0.2886751345948129*p0_over_gamma[0]*f[2]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[7] += volFact*(p0_over_gamma[1]*f[3]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2581988897471611*p0_over_gamma[1]*f[7]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[8] += volFact*(1.0000000000000002*p0_over_gamma[1]*f[6]*wx1+p0_over_gamma[0]*f[4]*wx1+0.28867513459481287*p0_over_gamma[0]*f[6]*dv1+0.2886751345948129*p0_over_gamma[1]*f[4]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*gamma[2]*f[5]+1.4142135623730951*gamma[1]*f[2]+1.4142135623730951*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
} 
