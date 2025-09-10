#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.7320508075688772*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.7320508075688772*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.7320508075688774*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.7320508075688774*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.7320508075688774*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.7320508075688774*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.7320508075688772*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.7320508075688772*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.7320508075688774*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.7320508075688774*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.7320508075688774*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.7320508075688774*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.7320508075688772*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.7320508075688772*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.7320508075688772*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.7320508075688772*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.7320508075688774*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.7320508075688772*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.7320508075688774*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.7320508075688774*gamma[18]*dv12; 

  out[0] += (p0_over_gamma[10]*f[14]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[10]*f[15]+p0_over_gamma[6]*f[13]+p0_over_gamma[5]*f[12]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p1_over_gamma[10]*f[14]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (p1_over_gamma[10]*f[15]+p1_over_gamma[6]*f[13]+p1_over_gamma[5]*f[12]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[4] += (p2_over_gamma[10]*f[14]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[5] += (p2_over_gamma[10]*f[15]+p2_over_gamma[6]*f[13]+p2_over_gamma[5]*f[12]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (gamma[10]*f[14]+gamma[6]*f[10]+gamma[5]*f[9]+gamma[4]*f[7]+gamma[3]*f[4]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[10]*f[15]+gamma[6]*f[13]+gamma[5]*f[12]+gamma[4]*f[11]+gamma[3]*f[8]+gamma[2]*f[6]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += volFact*(2.8284271247461907*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[1] += volFact*(2.8284271247461907*f[1]*wx1+0.8164965809277261*f[5]*dv1); 
  out[2] += volFact*(2.8284271247461907*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[3] += volFact*(2.8284271247461907*f[1]*wx2+0.8164965809277261*f[6]*dv2); 
  out[4] += volFact*(2.8284271247461907*f[0]*wx3+0.8164965809277261*f[4]*dv3); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx3+0.8164965809277261*f[8]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_Ni_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.7320508075688772*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.7320508075688772*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.7320508075688774*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.7320508075688774*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.7320508075688774*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.7320508075688774*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.7320508075688772*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.7320508075688772*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.7320508075688774*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.7320508075688774*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.7320508075688774*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.7320508075688774*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.7320508075688772*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.7320508075688772*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.7320508075688772*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.7320508075688772*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.7320508075688774*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.7320508075688772*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.7320508075688774*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.7320508075688774*gamma[18]*dv12; 

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += (p0_over_gamma[10]*f[14]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[3] += (p0_over_gamma[10]*f[15]+p0_over_gamma[6]*f[13]+p0_over_gamma[5]*f[12]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[4] += (p1_over_gamma[10]*f[14]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (p1_over_gamma[10]*f[15]+p1_over_gamma[6]*f[13]+p1_over_gamma[5]*f[12]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (p2_over_gamma[10]*f[14]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[7] += (p2_over_gamma[10]*f[15]+p2_over_gamma[6]*f[13]+p2_over_gamma[5]*f[12]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double dv11 = 2.0/dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.7320508075688772*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.7320508075688772*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.7320508075688772*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.7320508075688772*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.7320508075688774*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.7320508075688774*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.7320508075688774*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.7320508075688774*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.7320508075688772*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.7320508075688772*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.7320508075688772*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.7320508075688772*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.7320508075688774*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.7320508075688774*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.7320508075688774*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.7320508075688774*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.7320508075688772*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.7320508075688772*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.7320508075688772*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.7320508075688772*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.7320508075688774*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.7320508075688772*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.7320508075688774*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.7320508075688774*gamma[18]*dv12; 

  out[0] += (gamma[10]*f[14]+gamma[6]*f[10]+gamma[5]*f[9]+gamma[4]*f[7]+gamma[3]*f[4]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[10]*f[15]+gamma[6]*f[13]+gamma[5]*f[12]+gamma[4]*f[11]+gamma[3]*f[8]+gamma[2]*f[6]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += volFact*(2.8284271247461907*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[3] += volFact*(2.8284271247461907*f[1]*wx1+0.8164965809277261*f[5]*dv1); 
  out[4] += volFact*(2.8284271247461907*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx2+0.8164965809277261*f[6]*dv2); 
  out[6] += volFact*(2.8284271247461907*f[0]*wx3+0.8164965809277261*f[4]*dv3); 
  out[7] += volFact*(2.8284271247461907*f[1]*wx3+0.8164965809277261*f[8]*dv3); 
  out[8] += volFact*(p0_over_gamma[10]*f[14]*wx1+p0_over_gamma[6]*f[10]*wx1+p0_over_gamma[5]*f[9]*wx1+p0_over_gamma[4]*f[7]*wx1+p0_over_gamma[3]*f[4]*wx1+p0_over_gamma[2]*f[3]*wx1+p0_over_gamma[1]*f[2]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[6]*f[14]*dv1+0.2886751345948129*f[10]*p0_over_gamma[10]*dv1+0.2886751345948129*p0_over_gamma[3]*f[9]*dv1+0.2886751345948129*p0_over_gamma[2]*f[7]*dv1+0.2886751345948129*f[4]*p0_over_gamma[5]*dv1+0.2886751345948129*f[3]*p0_over_gamma[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[2]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[9] += volFact*(p0_over_gamma[10]*f[15]*wx1+p0_over_gamma[6]*f[13]*wx1+p0_over_gamma[5]*f[12]*wx1+p0_over_gamma[4]*f[11]*wx1+p0_over_gamma[3]*f[8]*wx1+p0_over_gamma[2]*f[6]*wx1+p0_over_gamma[1]*f[5]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[6]*f[15]*dv1+0.2886751345948129*p0_over_gamma[10]*f[13]*dv1+0.2886751345948129*p0_over_gamma[3]*f[12]*dv1+0.2886751345948129*p0_over_gamma[2]*f[11]*dv1+0.2886751345948129*p0_over_gamma[5]*f[8]*dv1+0.2886751345948129*p0_over_gamma[4]*f[6]*dv1+0.2886751345948129*p0_over_gamma[0]*f[5]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[10] += volFact*(p0_over_gamma[10]*f[14]*wx2+p0_over_gamma[6]*f[10]*wx2+p0_over_gamma[5]*f[9]*wx2+p0_over_gamma[4]*f[7]*wx2+p0_over_gamma[3]*f[4]*wx2+p0_over_gamma[2]*f[3]*wx2+p0_over_gamma[1]*f[2]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2581988897471611*f[10]*p0_over_gamma[14]*dv2+0.2886751345948129*p0_over_gamma[5]*f[14]*dv2+0.2886751345948129*f[9]*p0_over_gamma[10]*dv2+0.2886751345948129*p0_over_gamma[3]*f[10]*dv2+0.25819888974716115*f[3]*p0_over_gamma[8]*dv2+0.2886751345948129*p0_over_gamma[1]*f[7]*dv2+0.2886751345948129*f[4]*p0_over_gamma[6]*dv2+0.2886751345948129*f[2]*p0_over_gamma[4]*dv2+0.2886751345948129*p0_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[11] += volFact*(p0_over_gamma[10]*f[15]*wx2+p0_over_gamma[6]*f[13]*wx2+p0_over_gamma[5]*f[12]*wx2+p0_over_gamma[4]*f[11]*wx2+p0_over_gamma[3]*f[8]*wx2+p0_over_gamma[2]*f[6]*wx2+p0_over_gamma[1]*f[5]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2886751345948129*p0_over_gamma[5]*f[15]*dv2+0.2581988897471611*f[13]*p0_over_gamma[14]*dv2+0.2886751345948129*p0_over_gamma[3]*f[13]*dv2+0.2886751345948129*p0_over_gamma[10]*f[12]*dv2+0.2886751345948129*p0_over_gamma[1]*f[11]*dv2+0.25819888974716115*f[6]*p0_over_gamma[8]*dv2+0.2886751345948129*p0_over_gamma[6]*f[8]*dv2+0.2886751345948129*p0_over_gamma[0]*f[6]*dv2+0.2886751345948129*p0_over_gamma[4]*f[5]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[12] += volFact*(p0_over_gamma[10]*f[14]*wx3+p0_over_gamma[6]*f[10]*wx3+p0_over_gamma[5]*f[9]*wx3+p0_over_gamma[4]*f[7]*wx3+p0_over_gamma[3]*f[4]*wx3+p0_over_gamma[2]*f[3]*wx3+p0_over_gamma[1]*f[2]*wx3+f[0]*p0_over_gamma[0]*wx3+0.2581988897471611*f[10]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p0_over_gamma[10]*dv3+0.2886751345948129*p0_over_gamma[2]*f[10]*dv3+0.25819888974716115*f[4]*p0_over_gamma[9]*dv3+0.2886751345948129*p0_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p0_over_gamma[6]*dv3+0.2886751345948129*f[2]*p0_over_gamma[5]*dv3+0.2886751345948129*p0_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p0_over_gamma[3]*dv3); 
  out[13] += volFact*(p0_over_gamma[10]*f[15]*wx3+p0_over_gamma[6]*f[13]*wx3+p0_over_gamma[5]*f[12]*wx3+p0_over_gamma[4]*f[11]*wx3+p0_over_gamma[3]*f[8]*wx3+p0_over_gamma[2]*f[6]*wx3+p0_over_gamma[1]*f[5]*wx3+p0_over_gamma[0]*f[1]*wx3+0.2581988897471611*f[13]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[4]*f[15]*dv3+0.2886751345948129*p0_over_gamma[2]*f[13]*dv3+0.2886751345948129*p0_over_gamma[1]*f[12]*dv3+0.2886751345948129*p0_over_gamma[10]*f[11]*dv3+0.25819888974716115*f[8]*p0_over_gamma[9]*dv3+0.2886751345948129*p0_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p0_over_gamma[6]*dv3+0.2886751345948129*f[5]*p0_over_gamma[5]*dv3+0.2886751345948129*f[1]*p0_over_gamma[3]*dv3); 
  out[14] += volFact*(p1_over_gamma[10]*f[14]*wx2+p1_over_gamma[6]*f[10]*wx2+p1_over_gamma[5]*f[9]*wx2+p1_over_gamma[4]*f[7]*wx2+p1_over_gamma[3]*f[4]*wx2+p1_over_gamma[2]*f[3]*wx2+p1_over_gamma[1]*f[2]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2886751345948129*p1_over_gamma[5]*f[14]*dv2+0.2886751345948129*f[9]*p1_over_gamma[10]*dv2+0.2886751345948129*p1_over_gamma[3]*f[10]*dv2+0.2886751345948129*p1_over_gamma[1]*f[7]*dv2+0.2886751345948129*f[4]*p1_over_gamma[6]*dv2+0.2886751345948129*f[2]*p1_over_gamma[4]*dv2+0.2886751345948129*p1_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[15] += volFact*(p1_over_gamma[10]*f[15]*wx2+p1_over_gamma[6]*f[13]*wx2+p1_over_gamma[5]*f[12]*wx2+p1_over_gamma[4]*f[11]*wx2+p1_over_gamma[3]*f[8]*wx2+p1_over_gamma[2]*f[6]*wx2+p1_over_gamma[1]*f[5]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2886751345948129*p1_over_gamma[5]*f[15]*dv2+0.2886751345948129*p1_over_gamma[3]*f[13]*dv2+0.2886751345948129*p1_over_gamma[10]*f[12]*dv2+0.2886751345948129*p1_over_gamma[1]*f[11]*dv2+0.2886751345948129*p1_over_gamma[6]*f[8]*dv2+0.2886751345948129*p1_over_gamma[0]*f[6]*dv2+0.2886751345948129*p1_over_gamma[4]*f[5]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
  out[16] += volFact*(p1_over_gamma[10]*f[14]*wx3+p1_over_gamma[6]*f[10]*wx3+p1_over_gamma[5]*f[9]*wx3+p1_over_gamma[4]*f[7]*wx3+p1_over_gamma[3]*f[4]*wx3+p1_over_gamma[2]*f[3]*wx3+p1_over_gamma[1]*f[2]*wx3+f[0]*p1_over_gamma[0]*wx3+0.2581988897471611*f[9]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p1_over_gamma[10]*dv3+0.2886751345948129*p1_over_gamma[2]*f[10]*dv3+0.25819888974716115*f[4]*p1_over_gamma[9]*dv3+0.2886751345948129*p1_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p1_over_gamma[6]*dv3+0.2886751345948129*f[2]*p1_over_gamma[5]*dv3+0.2886751345948129*p1_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p1_over_gamma[3]*dv3); 
  out[17] += volFact*(p1_over_gamma[10]*f[15]*wx3+p1_over_gamma[6]*f[13]*wx3+p1_over_gamma[5]*f[12]*wx3+p1_over_gamma[4]*f[11]*wx3+p1_over_gamma[3]*f[8]*wx3+p1_over_gamma[2]*f[6]*wx3+p1_over_gamma[1]*f[5]*wx3+p1_over_gamma[0]*f[1]*wx3+0.2581988897471611*f[12]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[4]*f[15]*dv3+0.2886751345948129*p1_over_gamma[2]*f[13]*dv3+0.2886751345948129*p1_over_gamma[1]*f[12]*dv3+0.2886751345948129*p1_over_gamma[10]*f[11]*dv3+0.25819888974716115*f[8]*p1_over_gamma[9]*dv3+0.2886751345948129*p1_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p1_over_gamma[6]*dv3+0.2886751345948129*f[5]*p1_over_gamma[5]*dv3+0.2886751345948129*f[1]*p1_over_gamma[3]*dv3); 
  out[18] += volFact*(p2_over_gamma[10]*f[14]*wx3+p2_over_gamma[6]*f[10]*wx3+p2_over_gamma[5]*f[9]*wx3+p2_over_gamma[4]*f[7]*wx3+p2_over_gamma[3]*f[4]*wx3+p2_over_gamma[2]*f[3]*wx3+p2_over_gamma[1]*f[2]*wx3+f[0]*p2_over_gamma[0]*wx3+0.2886751345948129*p2_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p2_over_gamma[10]*dv3+0.2886751345948129*p2_over_gamma[2]*f[10]*dv3+0.2886751345948129*p2_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p2_over_gamma[6]*dv3+0.2886751345948129*f[2]*p2_over_gamma[5]*dv3+0.2886751345948129*p2_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p2_over_gamma[3]*dv3); 
  out[19] += volFact*(p2_over_gamma[10]*f[15]*wx3+p2_over_gamma[6]*f[13]*wx3+p2_over_gamma[5]*f[12]*wx3+p2_over_gamma[4]*f[11]*wx3+p2_over_gamma[3]*f[8]*wx3+p2_over_gamma[2]*f[6]*wx3+p2_over_gamma[1]*f[5]*wx3+p2_over_gamma[0]*f[1]*wx3+0.2886751345948129*p2_over_gamma[4]*f[15]*dv3+0.2886751345948129*p2_over_gamma[2]*f[13]*dv3+0.2886751345948129*p2_over_gamma[1]*f[12]*dv3+0.2886751345948129*p2_over_gamma[10]*f[11]*dv3+0.2886751345948129*p2_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p2_over_gamma[6]*dv3+0.2886751345948129*f[5]*p2_over_gamma[5]*dv3+0.2886751345948129*f[1]*p2_over_gamma[3]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*gamma[10]*f[14]+1.4142135623730951*gamma[6]*f[10]+1.4142135623730951*gamma[5]*f[9]+1.4142135623730951*gamma[4]*f[7]+1.4142135623730951*gamma[3]*f[4]+1.4142135623730951*gamma[2]*f[3]+1.4142135623730951*gamma[1]*f[2]+1.4142135623730951*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(4.0*f[0]*wx1+1.1547005383792517*f[2]*dv1); 
  out[3] += volFact*(4.0*f[0]*wx2+1.1547005383792517*f[3]*dv2); 
  out[4] += volFact*(4.0*f[0]*wx3+1.1547005383792517*f[4]*dv3); 
} 
