#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
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

  out[0] += (p0_over_gamma[10]*f[25]+p0_over_gamma[6]*f[15]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[5]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[10]*f[29]+p0_over_gamma[6]*f[23]+p0_over_gamma[5]*f[21]+p0_over_gamma[4]*f[18]+p0_over_gamma[3]*f[12]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[10]*f[30]+p0_over_gamma[6]*f[24]+p0_over_gamma[5]*f[22]+p0_over_gamma[4]*f[19]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[8]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[10]*f[31]+p0_over_gamma[6]*f[28]+p0_over_gamma[5]*f[27]+p0_over_gamma[4]*f[26]+p0_over_gamma[3]*f[20]+p0_over_gamma[2]*f[17]+p0_over_gamma[1]*f[16]+p0_over_gamma[0]*f[6])*volFact; 
  out[4] += (p1_over_gamma[10]*f[25]+p1_over_gamma[6]*f[15]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[5]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (p1_over_gamma[10]*f[29]+p1_over_gamma[6]*f[23]+p1_over_gamma[5]*f[21]+p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[12]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (p1_over_gamma[10]*f[30]+p1_over_gamma[6]*f[24]+p1_over_gamma[5]*f[22]+p1_over_gamma[4]*f[19]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[10]+p1_over_gamma[1]*f[8]+p1_over_gamma[0]*f[2])*volFact; 
  out[7] += (p1_over_gamma[10]*f[31]+p1_over_gamma[6]*f[28]+p1_over_gamma[5]*f[27]+p1_over_gamma[4]*f[26]+p1_over_gamma[3]*f[20]+p1_over_gamma[2]*f[17]+p1_over_gamma[1]*f[16]+p1_over_gamma[0]*f[6])*volFact; 
  out[8] += (p2_over_gamma[10]*f[25]+p2_over_gamma[6]*f[15]+p2_over_gamma[5]*f[14]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[5]+p2_over_gamma[2]*f[4]+p2_over_gamma[1]*f[3]+f[0]*p2_over_gamma[0])*volFact; 
  out[9] += (p2_over_gamma[10]*f[29]+p2_over_gamma[6]*f[23]+p2_over_gamma[5]*f[21]+p2_over_gamma[4]*f[18]+p2_over_gamma[3]*f[12]+p2_over_gamma[2]*f[9]+p2_over_gamma[1]*f[7]+p2_over_gamma[0]*f[1])*volFact; 
  out[10] += (p2_over_gamma[10]*f[30]+p2_over_gamma[6]*f[24]+p2_over_gamma[5]*f[22]+p2_over_gamma[4]*f[19]+p2_over_gamma[3]*f[13]+p2_over_gamma[2]*f[10]+p2_over_gamma[1]*f[8]+p2_over_gamma[0]*f[2])*volFact; 
  out[11] += (p2_over_gamma[10]*f[31]+p2_over_gamma[6]*f[28]+p2_over_gamma[5]*f[27]+p2_over_gamma[4]*f[26]+p2_over_gamma[3]*f[20]+p2_over_gamma[2]*f[17]+p2_over_gamma[1]*f[16]+p2_over_gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += (gamma[10]*f[25]+gamma[6]*f[15]+gamma[5]*f[14]+gamma[4]*f[11]+gamma[3]*f[5]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[10]*f[29]+gamma[6]*f[23]+gamma[5]*f[21]+gamma[4]*f[18]+gamma[3]*f[12]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[1])*volFact; 
  out[2] += (gamma[10]*f[30]+gamma[6]*f[24]+gamma[5]*f[22]+gamma[4]*f[19]+gamma[3]*f[13]+gamma[2]*f[10]+gamma[1]*f[8]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[10]*f[31]+gamma[6]*f[28]+gamma[5]*f[27]+gamma[4]*f[26]+gamma[3]*f[20]+gamma[2]*f[17]+gamma[1]*f[16]+gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += volFact*(2.8284271247461907*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
  out[1] += volFact*(2.8284271247461907*f[1]*wx1+0.8164965809277261*f[7]*dv1); 
  out[2] += volFact*(2.8284271247461907*f[2]*wx1+0.8164965809277261*f[8]*dv1); 
  out[3] += volFact*(2.8284271247461907*f[6]*wx1+0.8164965809277261*f[16]*dv1); 
  out[4] += volFact*(2.8284271247461907*f[0]*wx2+0.8164965809277261*f[4]*dv2); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx2+0.8164965809277261*f[9]*dv2); 
  out[6] += volFact*(2.8284271247461907*f[2]*wx2+0.8164965809277261*f[10]*dv2); 
  out[7] += volFact*(2.8284271247461907*f[6]*wx2+0.8164965809277261*f[17]*dv2); 
  out[8] += volFact*(2.8284271247461907*f[0]*wx3+0.8164965809277261*f[5]*dv3); 
  out[9] += volFact*(2.8284271247461907*f[1]*wx3+0.8164965809277261*f[12]*dv3); 
  out[10] += volFact*(2.8284271247461907*f[2]*wx3+0.8164965809277261*f[13]*dv3); 
  out[11] += volFact*(2.8284271247461907*f[6]*wx3+0.8164965809277261*f[20]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_Ni_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
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
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[6]*volFact; 
  out[4] += (p0_over_gamma[10]*f[25]+p0_over_gamma[6]*f[15]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[5]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[5] += (p0_over_gamma[10]*f[29]+p0_over_gamma[6]*f[23]+p0_over_gamma[5]*f[21]+p0_over_gamma[4]*f[18]+p0_over_gamma[3]*f[12]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[1])*volFact; 
  out[6] += (p0_over_gamma[10]*f[30]+p0_over_gamma[6]*f[24]+p0_over_gamma[5]*f[22]+p0_over_gamma[4]*f[19]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[8]+p0_over_gamma[0]*f[2])*volFact; 
  out[7] += (p0_over_gamma[10]*f[31]+p0_over_gamma[6]*f[28]+p0_over_gamma[5]*f[27]+p0_over_gamma[4]*f[26]+p0_over_gamma[3]*f[20]+p0_over_gamma[2]*f[17]+p0_over_gamma[1]*f[16]+p0_over_gamma[0]*f[6])*volFact; 
  out[8] += (p1_over_gamma[10]*f[25]+p1_over_gamma[6]*f[15]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[5]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[9] += (p1_over_gamma[10]*f[29]+p1_over_gamma[6]*f[23]+p1_over_gamma[5]*f[21]+p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[12]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[1])*volFact; 
  out[10] += (p1_over_gamma[10]*f[30]+p1_over_gamma[6]*f[24]+p1_over_gamma[5]*f[22]+p1_over_gamma[4]*f[19]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[10]+p1_over_gamma[1]*f[8]+p1_over_gamma[0]*f[2])*volFact; 
  out[11] += (p1_over_gamma[10]*f[31]+p1_over_gamma[6]*f[28]+p1_over_gamma[5]*f[27]+p1_over_gamma[4]*f[26]+p1_over_gamma[3]*f[20]+p1_over_gamma[2]*f[17]+p1_over_gamma[1]*f[16]+p1_over_gamma[0]*f[6])*volFact; 
  out[12] += (p2_over_gamma[10]*f[25]+p2_over_gamma[6]*f[15]+p2_over_gamma[5]*f[14]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[5]+p2_over_gamma[2]*f[4]+p2_over_gamma[1]*f[3]+f[0]*p2_over_gamma[0])*volFact; 
  out[13] += (p2_over_gamma[10]*f[29]+p2_over_gamma[6]*f[23]+p2_over_gamma[5]*f[21]+p2_over_gamma[4]*f[18]+p2_over_gamma[3]*f[12]+p2_over_gamma[2]*f[9]+p2_over_gamma[1]*f[7]+p2_over_gamma[0]*f[1])*volFact; 
  out[14] += (p2_over_gamma[10]*f[30]+p2_over_gamma[6]*f[24]+p2_over_gamma[5]*f[22]+p2_over_gamma[4]*f[19]+p2_over_gamma[3]*f[13]+p2_over_gamma[2]*f[10]+p2_over_gamma[1]*f[8]+p2_over_gamma[0]*f[2])*volFact; 
  out[15] += (p2_over_gamma[10]*f[31]+p2_over_gamma[6]*f[28]+p2_over_gamma[5]*f[27]+p2_over_gamma[4]*f[26]+p2_over_gamma[3]*f[20]+p2_over_gamma[2]*f[17]+p2_over_gamma[1]*f[16]+p2_over_gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double dv12 = 2.0/dxv[4]; 
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

  out[0] += (gamma[10]*f[25]+gamma[6]*f[15]+gamma[5]*f[14]+gamma[4]*f[11]+gamma[3]*f[5]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (gamma[10]*f[29]+gamma[6]*f[23]+gamma[5]*f[21]+gamma[4]*f[18]+gamma[3]*f[12]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[1])*volFact; 
  out[2] += (gamma[10]*f[30]+gamma[6]*f[24]+gamma[5]*f[22]+gamma[4]*f[19]+gamma[3]*f[13]+gamma[2]*f[10]+gamma[1]*f[8]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[10]*f[31]+gamma[6]*f[28]+gamma[5]*f[27]+gamma[4]*f[26]+gamma[3]*f[20]+gamma[2]*f[17]+gamma[1]*f[16]+gamma[0]*f[6])*volFact; 
  out[4] += volFact*(2.8284271247461907*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx1+0.8164965809277261*f[7]*dv1); 
  out[6] += volFact*(2.8284271247461907*f[2]*wx1+0.8164965809277261*f[8]*dv1); 
  out[7] += volFact*(2.8284271247461907*f[6]*wx1+0.8164965809277261*f[16]*dv1); 
  out[8] += volFact*(2.8284271247461907*f[0]*wx2+0.8164965809277261*f[4]*dv2); 
  out[9] += volFact*(2.8284271247461907*f[1]*wx2+0.8164965809277261*f[9]*dv2); 
  out[10] += volFact*(2.8284271247461907*f[2]*wx2+0.8164965809277261*f[10]*dv2); 
  out[11] += volFact*(2.8284271247461907*f[6]*wx2+0.8164965809277261*f[17]*dv2); 
  out[12] += volFact*(2.8284271247461907*f[0]*wx3+0.8164965809277261*f[5]*dv3); 
  out[13] += volFact*(2.8284271247461907*f[1]*wx3+0.8164965809277261*f[12]*dv3); 
  out[14] += volFact*(2.8284271247461907*f[2]*wx3+0.8164965809277261*f[13]*dv3); 
  out[15] += volFact*(2.8284271247461907*f[6]*wx3+0.8164965809277261*f[20]*dv3); 
  out[16] += volFact*(p0_over_gamma[10]*f[25]*wx1+p0_over_gamma[6]*f[15]*wx1+p0_over_gamma[5]*f[14]*wx1+p0_over_gamma[4]*f[11]*wx1+p0_over_gamma[3]*f[5]*wx1+p0_over_gamma[2]*f[4]*wx1+p0_over_gamma[1]*f[3]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[6]*f[25]*dv1+0.2886751345948129*p0_over_gamma[10]*f[15]*dv1+0.2886751345948129*p0_over_gamma[3]*f[14]*dv1+0.2886751345948129*p0_over_gamma[2]*f[11]*dv1+0.2886751345948129*f[5]*p0_over_gamma[5]*dv1+0.2886751345948129*f[4]*p0_over_gamma[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[17] += volFact*(p0_over_gamma[10]*f[29]*wx1+p0_over_gamma[6]*f[23]*wx1+p0_over_gamma[5]*f[21]*wx1+p0_over_gamma[4]*f[18]*wx1+p0_over_gamma[3]*f[12]*wx1+p0_over_gamma[2]*f[9]*wx1+p0_over_gamma[1]*f[7]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[6]*f[29]*dv1+0.2886751345948129*p0_over_gamma[10]*f[23]*dv1+0.2886751345948129*p0_over_gamma[3]*f[21]*dv1+0.2886751345948129*p0_over_gamma[2]*f[18]*dv1+0.2886751345948129*p0_over_gamma[5]*f[12]*dv1+0.2886751345948129*p0_over_gamma[4]*f[9]*dv1+0.2886751345948129*p0_over_gamma[0]*f[7]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[18] += volFact*(p0_over_gamma[10]*f[30]*wx1+p0_over_gamma[6]*f[24]*wx1+p0_over_gamma[5]*f[22]*wx1+p0_over_gamma[4]*f[19]*wx1+p0_over_gamma[3]*f[13]*wx1+p0_over_gamma[2]*f[10]*wx1+p0_over_gamma[1]*f[8]*wx1+p0_over_gamma[0]*f[2]*wx1+0.2886751345948129*p0_over_gamma[6]*f[30]*dv1+0.2886751345948129*p0_over_gamma[10]*f[24]*dv1+0.2886751345948129*p0_over_gamma[3]*f[22]*dv1+0.2886751345948129*p0_over_gamma[2]*f[19]*dv1+0.2886751345948129*p0_over_gamma[5]*f[13]*dv1+0.2886751345948129*p0_over_gamma[4]*f[10]*dv1+0.2886751345948129*p0_over_gamma[0]*f[8]*dv1+0.2886751345948129*p0_over_gamma[1]*f[2]*dv1); 
  out[19] += volFact*(p0_over_gamma[10]*f[31]*wx1+p0_over_gamma[6]*f[28]*wx1+p0_over_gamma[5]*f[27]*wx1+p0_over_gamma[4]*f[26]*wx1+p0_over_gamma[3]*f[20]*wx1+p0_over_gamma[2]*f[17]*wx1+p0_over_gamma[1]*f[16]*wx1+p0_over_gamma[0]*f[6]*wx1+0.2886751345948129*p0_over_gamma[6]*f[31]*dv1+0.2886751345948129*p0_over_gamma[10]*f[28]*dv1+0.2886751345948129*p0_over_gamma[3]*f[27]*dv1+0.2886751345948129*p0_over_gamma[2]*f[26]*dv1+0.2886751345948129*p0_over_gamma[5]*f[20]*dv1+0.2886751345948129*p0_over_gamma[4]*f[17]*dv1+0.2886751345948129*p0_over_gamma[0]*f[16]*dv1+0.2886751345948129*p0_over_gamma[1]*f[6]*dv1); 
  out[20] += volFact*(p0_over_gamma[10]*f[25]*wx2+p0_over_gamma[6]*f[15]*wx2+p0_over_gamma[5]*f[14]*wx2+p0_over_gamma[4]*f[11]*wx2+p0_over_gamma[3]*f[5]*wx2+p0_over_gamma[2]*f[4]*wx2+p0_over_gamma[1]*f[3]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2886751345948129*p0_over_gamma[5]*f[25]*dv2+0.2581988897471611*p0_over_gamma[14]*f[15]*dv2+0.2886751345948129*p0_over_gamma[3]*f[15]*dv2+0.2886751345948129*p0_over_gamma[10]*f[14]*dv2+0.2886751345948129*p0_over_gamma[1]*f[11]*dv2+0.25819888974716115*f[4]*p0_over_gamma[8]*dv2+0.2886751345948129*f[5]*p0_over_gamma[6]*dv2+0.2886751345948129*f[3]*p0_over_gamma[4]*dv2+0.2886751345948129*p0_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[21] += volFact*(p0_over_gamma[10]*f[29]*wx2+p0_over_gamma[6]*f[23]*wx2+p0_over_gamma[5]*f[21]*wx2+p0_over_gamma[4]*f[18]*wx2+p0_over_gamma[3]*f[12]*wx2+p0_over_gamma[2]*f[9]*wx2+p0_over_gamma[1]*f[7]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2886751345948129*p0_over_gamma[5]*f[29]*dv2+0.2581988897471611*p0_over_gamma[14]*f[23]*dv2+0.2886751345948129*p0_over_gamma[3]*f[23]*dv2+0.2886751345948129*p0_over_gamma[10]*f[21]*dv2+0.2886751345948129*p0_over_gamma[1]*f[18]*dv2+0.2886751345948129*p0_over_gamma[6]*f[12]*dv2+0.25819888974716115*p0_over_gamma[8]*f[9]*dv2+0.2886751345948129*p0_over_gamma[0]*f[9]*dv2+0.2886751345948129*p0_over_gamma[4]*f[7]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[22] += volFact*(p0_over_gamma[10]*f[30]*wx2+p0_over_gamma[6]*f[24]*wx2+p0_over_gamma[5]*f[22]*wx2+p0_over_gamma[4]*f[19]*wx2+p0_over_gamma[3]*f[13]*wx2+p0_over_gamma[2]*f[10]*wx2+p0_over_gamma[1]*f[8]*wx2+p0_over_gamma[0]*f[2]*wx2+0.2886751345948129*p0_over_gamma[5]*f[30]*dv2+0.2581988897471611*p0_over_gamma[14]*f[24]*dv2+0.2886751345948129*p0_over_gamma[3]*f[24]*dv2+0.2886751345948129*p0_over_gamma[10]*f[22]*dv2+0.2886751345948129*p0_over_gamma[1]*f[19]*dv2+0.2886751345948129*p0_over_gamma[6]*f[13]*dv2+0.25819888974716115*p0_over_gamma[8]*f[10]*dv2+0.2886751345948129*p0_over_gamma[0]*f[10]*dv2+0.2886751345948129*p0_over_gamma[4]*f[8]*dv2+0.2886751345948129*f[2]*p0_over_gamma[2]*dv2); 
  out[23] += volFact*(p0_over_gamma[10]*f[31]*wx2+p0_over_gamma[6]*f[28]*wx2+p0_over_gamma[5]*f[27]*wx2+p0_over_gamma[4]*f[26]*wx2+p0_over_gamma[3]*f[20]*wx2+p0_over_gamma[2]*f[17]*wx2+p0_over_gamma[1]*f[16]*wx2+p0_over_gamma[0]*f[6]*wx2+0.2886751345948129*p0_over_gamma[5]*f[31]*dv2+0.2581988897471611*p0_over_gamma[14]*f[28]*dv2+0.2886751345948129*p0_over_gamma[3]*f[28]*dv2+0.2886751345948129*p0_over_gamma[10]*f[27]*dv2+0.2886751345948129*p0_over_gamma[1]*f[26]*dv2+0.2886751345948129*p0_over_gamma[6]*f[20]*dv2+0.25819888974716115*p0_over_gamma[8]*f[17]*dv2+0.2886751345948129*p0_over_gamma[0]*f[17]*dv2+0.2886751345948129*p0_over_gamma[4]*f[16]*dv2+0.2886751345948129*p0_over_gamma[2]*f[6]*dv2); 
  out[24] += volFact*(p0_over_gamma[10]*f[25]*wx3+p0_over_gamma[6]*f[15]*wx3+p0_over_gamma[5]*f[14]*wx3+p0_over_gamma[4]*f[11]*wx3+p0_over_gamma[3]*f[5]*wx3+p0_over_gamma[2]*f[4]*wx3+p0_over_gamma[1]*f[3]*wx3+f[0]*p0_over_gamma[0]*wx3+0.2886751345948129*p0_over_gamma[4]*f[25]*dv3+0.2581988897471611*f[15]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[2]*f[15]*dv3+0.2886751345948129*p0_over_gamma[1]*f[14]*dv3+0.2886751345948129*p0_over_gamma[10]*f[11]*dv3+0.25819888974716115*f[5]*p0_over_gamma[9]*dv3+0.2886751345948129*f[4]*p0_over_gamma[6]*dv3+0.2886751345948129*f[3]*p0_over_gamma[5]*dv3+0.2886751345948129*p0_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p0_over_gamma[3]*dv3); 
  out[25] += volFact*(p0_over_gamma[10]*f[29]*wx3+p0_over_gamma[6]*f[23]*wx3+p0_over_gamma[5]*f[21]*wx3+p0_over_gamma[4]*f[18]*wx3+p0_over_gamma[3]*f[12]*wx3+p0_over_gamma[2]*f[9]*wx3+p0_over_gamma[1]*f[7]*wx3+p0_over_gamma[0]*f[1]*wx3+0.2886751345948129*p0_over_gamma[4]*f[29]*dv3+0.2581988897471611*p0_over_gamma[16]*f[23]*dv3+0.2886751345948129*p0_over_gamma[2]*f[23]*dv3+0.2886751345948129*p0_over_gamma[1]*f[21]*dv3+0.2886751345948129*p0_over_gamma[10]*f[18]*dv3+0.25819888974716115*p0_over_gamma[9]*f[12]*dv3+0.2886751345948129*p0_over_gamma[0]*f[12]*dv3+0.2886751345948129*p0_over_gamma[6]*f[9]*dv3+0.2886751345948129*p0_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p0_over_gamma[3]*dv3); 
  out[26] += volFact*(p0_over_gamma[10]*f[30]*wx3+p0_over_gamma[6]*f[24]*wx3+p0_over_gamma[5]*f[22]*wx3+p0_over_gamma[4]*f[19]*wx3+p0_over_gamma[3]*f[13]*wx3+p0_over_gamma[2]*f[10]*wx3+p0_over_gamma[1]*f[8]*wx3+p0_over_gamma[0]*f[2]*wx3+0.2886751345948129*p0_over_gamma[4]*f[30]*dv3+0.2581988897471611*p0_over_gamma[16]*f[24]*dv3+0.2886751345948129*p0_over_gamma[2]*f[24]*dv3+0.2886751345948129*p0_over_gamma[1]*f[22]*dv3+0.2886751345948129*p0_over_gamma[10]*f[19]*dv3+0.25819888974716115*p0_over_gamma[9]*f[13]*dv3+0.2886751345948129*p0_over_gamma[0]*f[13]*dv3+0.2886751345948129*p0_over_gamma[6]*f[10]*dv3+0.2886751345948129*p0_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p0_over_gamma[3]*dv3); 
  out[27] += volFact*(p0_over_gamma[10]*f[31]*wx3+p0_over_gamma[6]*f[28]*wx3+p0_over_gamma[5]*f[27]*wx3+p0_over_gamma[4]*f[26]*wx3+p0_over_gamma[3]*f[20]*wx3+p0_over_gamma[2]*f[17]*wx3+p0_over_gamma[1]*f[16]*wx3+p0_over_gamma[0]*f[6]*wx3+0.2886751345948129*p0_over_gamma[4]*f[31]*dv3+0.2581988897471611*p0_over_gamma[16]*f[28]*dv3+0.2886751345948129*p0_over_gamma[2]*f[28]*dv3+0.2886751345948129*p0_over_gamma[1]*f[27]*dv3+0.2886751345948129*p0_over_gamma[10]*f[26]*dv3+0.25819888974716115*p0_over_gamma[9]*f[20]*dv3+0.2886751345948129*p0_over_gamma[0]*f[20]*dv3+0.2886751345948129*p0_over_gamma[6]*f[17]*dv3+0.2886751345948129*p0_over_gamma[5]*f[16]*dv3+0.2886751345948129*p0_over_gamma[3]*f[6]*dv3); 
  out[28] += volFact*(p1_over_gamma[10]*f[25]*wx2+p1_over_gamma[6]*f[15]*wx2+p1_over_gamma[5]*f[14]*wx2+p1_over_gamma[4]*f[11]*wx2+p1_over_gamma[3]*f[5]*wx2+p1_over_gamma[2]*f[4]*wx2+p1_over_gamma[1]*f[3]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2886751345948129*p1_over_gamma[5]*f[25]*dv2+0.2886751345948129*p1_over_gamma[3]*f[15]*dv2+0.2886751345948129*p1_over_gamma[10]*f[14]*dv2+0.2886751345948129*p1_over_gamma[1]*f[11]*dv2+0.2886751345948129*f[5]*p1_over_gamma[6]*dv2+0.2886751345948129*f[3]*p1_over_gamma[4]*dv2+0.2886751345948129*p1_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[29] += volFact*(p1_over_gamma[10]*f[29]*wx2+p1_over_gamma[6]*f[23]*wx2+p1_over_gamma[5]*f[21]*wx2+p1_over_gamma[4]*f[18]*wx2+p1_over_gamma[3]*f[12]*wx2+p1_over_gamma[2]*f[9]*wx2+p1_over_gamma[1]*f[7]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2886751345948129*p1_over_gamma[5]*f[29]*dv2+0.2886751345948129*p1_over_gamma[3]*f[23]*dv2+0.2886751345948129*p1_over_gamma[10]*f[21]*dv2+0.2886751345948129*p1_over_gamma[1]*f[18]*dv2+0.2886751345948129*p1_over_gamma[6]*f[12]*dv2+0.2886751345948129*p1_over_gamma[0]*f[9]*dv2+0.2886751345948129*p1_over_gamma[4]*f[7]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
  out[30] += volFact*(p1_over_gamma[10]*f[30]*wx2+p1_over_gamma[6]*f[24]*wx2+p1_over_gamma[5]*f[22]*wx2+p1_over_gamma[4]*f[19]*wx2+p1_over_gamma[3]*f[13]*wx2+p1_over_gamma[2]*f[10]*wx2+p1_over_gamma[1]*f[8]*wx2+p1_over_gamma[0]*f[2]*wx2+0.2886751345948129*p1_over_gamma[5]*f[30]*dv2+0.2886751345948129*p1_over_gamma[3]*f[24]*dv2+0.2886751345948129*p1_over_gamma[10]*f[22]*dv2+0.2886751345948129*p1_over_gamma[1]*f[19]*dv2+0.2886751345948129*p1_over_gamma[6]*f[13]*dv2+0.2886751345948129*p1_over_gamma[0]*f[10]*dv2+0.2886751345948129*p1_over_gamma[4]*f[8]*dv2+0.2886751345948129*f[2]*p1_over_gamma[2]*dv2); 
  out[31] += volFact*(p1_over_gamma[10]*f[31]*wx2+p1_over_gamma[6]*f[28]*wx2+p1_over_gamma[5]*f[27]*wx2+p1_over_gamma[4]*f[26]*wx2+p1_over_gamma[3]*f[20]*wx2+p1_over_gamma[2]*f[17]*wx2+p1_over_gamma[1]*f[16]*wx2+p1_over_gamma[0]*f[6]*wx2+0.2886751345948129*p1_over_gamma[5]*f[31]*dv2+0.2886751345948129*p1_over_gamma[3]*f[28]*dv2+0.2886751345948129*p1_over_gamma[10]*f[27]*dv2+0.2886751345948129*p1_over_gamma[1]*f[26]*dv2+0.2886751345948129*p1_over_gamma[6]*f[20]*dv2+0.2886751345948129*p1_over_gamma[0]*f[17]*dv2+0.2886751345948129*p1_over_gamma[4]*f[16]*dv2+0.2886751345948129*p1_over_gamma[2]*f[6]*dv2); 
  out[32] += volFact*(p1_over_gamma[10]*f[25]*wx3+p1_over_gamma[6]*f[15]*wx3+p1_over_gamma[5]*f[14]*wx3+p1_over_gamma[4]*f[11]*wx3+p1_over_gamma[3]*f[5]*wx3+p1_over_gamma[2]*f[4]*wx3+p1_over_gamma[1]*f[3]*wx3+f[0]*p1_over_gamma[0]*wx3+0.2886751345948129*p1_over_gamma[4]*f[25]*dv3+0.2581988897471611*f[14]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[2]*f[15]*dv3+0.2886751345948129*p1_over_gamma[1]*f[14]*dv3+0.2886751345948129*p1_over_gamma[10]*f[11]*dv3+0.25819888974716115*f[5]*p1_over_gamma[9]*dv3+0.2886751345948129*f[4]*p1_over_gamma[6]*dv3+0.2886751345948129*f[3]*p1_over_gamma[5]*dv3+0.2886751345948129*p1_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p1_over_gamma[3]*dv3); 
  out[33] += volFact*(p1_over_gamma[10]*f[29]*wx3+p1_over_gamma[6]*f[23]*wx3+p1_over_gamma[5]*f[21]*wx3+p1_over_gamma[4]*f[18]*wx3+p1_over_gamma[3]*f[12]*wx3+p1_over_gamma[2]*f[9]*wx3+p1_over_gamma[1]*f[7]*wx3+p1_over_gamma[0]*f[1]*wx3+0.2886751345948129*p1_over_gamma[4]*f[29]*dv3+0.2886751345948129*p1_over_gamma[2]*f[23]*dv3+0.2581988897471611*p1_over_gamma[15]*f[21]*dv3+0.2886751345948129*p1_over_gamma[1]*f[21]*dv3+0.2886751345948129*p1_over_gamma[10]*f[18]*dv3+0.25819888974716115*p1_over_gamma[9]*f[12]*dv3+0.2886751345948129*p1_over_gamma[0]*f[12]*dv3+0.2886751345948129*p1_over_gamma[6]*f[9]*dv3+0.2886751345948129*p1_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p1_over_gamma[3]*dv3); 
  out[34] += volFact*(p1_over_gamma[10]*f[30]*wx3+p1_over_gamma[6]*f[24]*wx3+p1_over_gamma[5]*f[22]*wx3+p1_over_gamma[4]*f[19]*wx3+p1_over_gamma[3]*f[13]*wx3+p1_over_gamma[2]*f[10]*wx3+p1_over_gamma[1]*f[8]*wx3+p1_over_gamma[0]*f[2]*wx3+0.2886751345948129*p1_over_gamma[4]*f[30]*dv3+0.2886751345948129*p1_over_gamma[2]*f[24]*dv3+0.2581988897471611*p1_over_gamma[15]*f[22]*dv3+0.2886751345948129*p1_over_gamma[1]*f[22]*dv3+0.2886751345948129*p1_over_gamma[10]*f[19]*dv3+0.25819888974716115*p1_over_gamma[9]*f[13]*dv3+0.2886751345948129*p1_over_gamma[0]*f[13]*dv3+0.2886751345948129*p1_over_gamma[6]*f[10]*dv3+0.2886751345948129*p1_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p1_over_gamma[3]*dv3); 
  out[35] += volFact*(p1_over_gamma[10]*f[31]*wx3+p1_over_gamma[6]*f[28]*wx3+p1_over_gamma[5]*f[27]*wx3+p1_over_gamma[4]*f[26]*wx3+p1_over_gamma[3]*f[20]*wx3+p1_over_gamma[2]*f[17]*wx3+p1_over_gamma[1]*f[16]*wx3+p1_over_gamma[0]*f[6]*wx3+0.2886751345948129*p1_over_gamma[4]*f[31]*dv3+0.2886751345948129*p1_over_gamma[2]*f[28]*dv3+0.2581988897471611*p1_over_gamma[15]*f[27]*dv3+0.2886751345948129*p1_over_gamma[1]*f[27]*dv3+0.2886751345948129*p1_over_gamma[10]*f[26]*dv3+0.25819888974716115*p1_over_gamma[9]*f[20]*dv3+0.2886751345948129*p1_over_gamma[0]*f[20]*dv3+0.2886751345948129*p1_over_gamma[6]*f[17]*dv3+0.2886751345948129*p1_over_gamma[5]*f[16]*dv3+0.2886751345948129*p1_over_gamma[3]*f[6]*dv3); 
  out[36] += volFact*(p2_over_gamma[10]*f[25]*wx3+p2_over_gamma[6]*f[15]*wx3+p2_over_gamma[5]*f[14]*wx3+p2_over_gamma[4]*f[11]*wx3+p2_over_gamma[3]*f[5]*wx3+p2_over_gamma[2]*f[4]*wx3+p2_over_gamma[1]*f[3]*wx3+f[0]*p2_over_gamma[0]*wx3+0.2886751345948129*p2_over_gamma[4]*f[25]*dv3+0.2886751345948129*p2_over_gamma[2]*f[15]*dv3+0.2886751345948129*p2_over_gamma[1]*f[14]*dv3+0.2886751345948129*p2_over_gamma[10]*f[11]*dv3+0.2886751345948129*f[4]*p2_over_gamma[6]*dv3+0.2886751345948129*f[3]*p2_over_gamma[5]*dv3+0.2886751345948129*p2_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p2_over_gamma[3]*dv3); 
  out[37] += volFact*(p2_over_gamma[10]*f[29]*wx3+p2_over_gamma[6]*f[23]*wx3+p2_over_gamma[5]*f[21]*wx3+p2_over_gamma[4]*f[18]*wx3+p2_over_gamma[3]*f[12]*wx3+p2_over_gamma[2]*f[9]*wx3+p2_over_gamma[1]*f[7]*wx3+p2_over_gamma[0]*f[1]*wx3+0.2886751345948129*p2_over_gamma[4]*f[29]*dv3+0.2886751345948129*p2_over_gamma[2]*f[23]*dv3+0.2886751345948129*p2_over_gamma[1]*f[21]*dv3+0.2886751345948129*p2_over_gamma[10]*f[18]*dv3+0.2886751345948129*p2_over_gamma[0]*f[12]*dv3+0.2886751345948129*p2_over_gamma[6]*f[9]*dv3+0.2886751345948129*p2_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p2_over_gamma[3]*dv3); 
  out[38] += volFact*(p2_over_gamma[10]*f[30]*wx3+p2_over_gamma[6]*f[24]*wx3+p2_over_gamma[5]*f[22]*wx3+p2_over_gamma[4]*f[19]*wx3+p2_over_gamma[3]*f[13]*wx3+p2_over_gamma[2]*f[10]*wx3+p2_over_gamma[1]*f[8]*wx3+p2_over_gamma[0]*f[2]*wx3+0.2886751345948129*p2_over_gamma[4]*f[30]*dv3+0.2886751345948129*p2_over_gamma[2]*f[24]*dv3+0.2886751345948129*p2_over_gamma[1]*f[22]*dv3+0.2886751345948129*p2_over_gamma[10]*f[19]*dv3+0.2886751345948129*p2_over_gamma[0]*f[13]*dv3+0.2886751345948129*p2_over_gamma[6]*f[10]*dv3+0.2886751345948129*p2_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p2_over_gamma[3]*dv3); 
  out[39] += volFact*(p2_over_gamma[10]*f[31]*wx3+p2_over_gamma[6]*f[28]*wx3+p2_over_gamma[5]*f[27]*wx3+p2_over_gamma[4]*f[26]*wx3+p2_over_gamma[3]*f[20]*wx3+p2_over_gamma[2]*f[17]*wx3+p2_over_gamma[1]*f[16]*wx3+p2_over_gamma[0]*f[6]*wx3+0.2886751345948129*p2_over_gamma[4]*f[31]*dv3+0.2886751345948129*p2_over_gamma[2]*f[28]*dv3+0.2886751345948129*p2_over_gamma[1]*f[27]*dv3+0.2886751345948129*p2_over_gamma[10]*f[26]*dv3+0.2886751345948129*p2_over_gamma[0]*f[20]*dv3+0.2886751345948129*p2_over_gamma[6]*f[17]*dv3+0.2886751345948129*p2_over_gamma[5]*f[16]*dv3+0.2886751345948129*p2_over_gamma[3]*f[6]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (2.0*gamma[10]*f[25]+2.0*gamma[6]*f[15]+2.0*gamma[5]*f[14]+2.0*gamma[4]*f[11]+2.0*gamma[3]*f[5]+2.0*gamma[2]*f[4]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(5.656854249492382*f[0]*wx1+1.6329931618554527*f[3]*dv1); 
  out[3] += volFact*(5.656854249492382*f[0]*wx2+1.6329931618554527*f[4]*dv2); 
  out[4] += volFact*(5.656854249492382*f[0]*wx3+1.6329931618554527*f[5]*dv3); 
} 
