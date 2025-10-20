#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.7320508075688774*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.7320508075688772*gamma[6]*dv11; 

  out[0] += (p0_over_gamma[3]*f[6]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[3]*f[7]+p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[4]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p1_over_gamma[3]*f[6]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (p1_over_gamma[3]*f[7]+p1_over_gamma[2]*f[5]+p1_over_gamma[1]*f[4]+p1_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (gamma[3]*f[6]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[3]*f[7]+gamma[2]*f[5]+gamma[1]*f[4]+gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
  out[2] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2); 
  out[3] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2); 
} 
GKYL_CU_DH void vlasov_sr_Ni_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.7320508075688774*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.7320508075688772*gamma[6]*dv11; 

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (p0_over_gamma[3]*f[6]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[3] += (p0_over_gamma[3]*f[7]+p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[4]+p0_over_gamma[0]*f[1])*volFact; 
  out[4] += (p1_over_gamma[3]*f[6]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (p1_over_gamma[3]*f[7]+p1_over_gamma[2]*f[5]+p1_over_gamma[1]*f[4]+p1_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double dv11 = 2.0/dxv[2]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.7320508075688774*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.7320508075688772*gamma[6]*dv11; 

  out[0] += (gamma[3]*f[6]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[3]*f[7]+gamma[2]*f[5]+gamma[1]*f[4]+gamma[0]*f[1])*volFact; 
  out[2] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[3] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
  out[4] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2); 
  out[5] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2); 
  out[6] += volFact*(p0_over_gamma[3]*f[6]*wx1+p0_over_gamma[2]*f[3]*wx1+p0_over_gamma[1]*f[2]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[2]*f[6]*dv1+0.2886751345948129*f[3]*p0_over_gamma[3]*dv1+0.2886751345948129*p0_over_gamma[0]*f[2]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[7] += volFact*(p0_over_gamma[3]*f[7]*wx1+p0_over_gamma[2]*f[5]*wx1+p0_over_gamma[1]*f[4]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[2]*f[7]*dv1+0.2886751345948129*p0_over_gamma[3]*f[5]*dv1+0.2886751345948129*p0_over_gamma[0]*f[4]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[8] += volFact*(p0_over_gamma[3]*f[6]*wx2+p0_over_gamma[2]*f[3]*wx2+p0_over_gamma[1]*f[2]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2886751345948129*p0_over_gamma[1]*f[6]*dv2+0.25819888974716115*f[3]*p0_over_gamma[5]*dv2+0.2886751345948129*f[2]*p0_over_gamma[3]*dv2+0.2886751345948129*p0_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[9] += volFact*(p0_over_gamma[3]*f[7]*wx2+p0_over_gamma[2]*f[5]*wx2+p0_over_gamma[1]*f[4]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2886751345948129*p0_over_gamma[1]*f[7]*dv2+0.25819888974716115*f[5]*p0_over_gamma[5]*dv2+0.2886751345948129*p0_over_gamma[0]*f[5]*dv2+0.2886751345948129*p0_over_gamma[3]*f[4]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[10] += volFact*(p1_over_gamma[3]*f[6]*wx2+p1_over_gamma[2]*f[3]*wx2+p1_over_gamma[1]*f[2]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2886751345948129*p1_over_gamma[1]*f[6]*dv2+0.2886751345948129*f[2]*p1_over_gamma[3]*dv2+0.2886751345948129*p1_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[11] += volFact*(p1_over_gamma[3]*f[7]*wx2+p1_over_gamma[2]*f[5]*wx2+p1_over_gamma[1]*f[4]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2886751345948129*p1_over_gamma[1]*f[7]*dv2+0.2886751345948129*p1_over_gamma[0]*f[5]*dv2+0.2886751345948129*p1_over_gamma[3]*f[4]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
} 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.4142135623730951*gamma[3]*f[6]+1.4142135623730951*gamma[2]*f[3]+1.4142135623730951*gamma[1]*f[2]+1.4142135623730951*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(2.8284271247461907*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[3] += volFact*(2.8284271247461907*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
} 
