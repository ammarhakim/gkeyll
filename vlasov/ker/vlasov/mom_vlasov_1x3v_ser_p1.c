#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += (1.7320508075688772*hamil[7]*f[10]+1.7320508075688772*f[4]*hamil[5]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (1.7320508075688772*hamil[7]*f[13]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*hamil[4]*f[6]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (1.7320508075688772*hamil[7]*f[9]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (1.7320508075688772*hamil[7]*f[12]+1.7320508075688772*hamil[6]*f[8]+1.7320508075688772*hamil[4]*f[5]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
  out[4] += (1.7320508075688772*f[7]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[2]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[6]*hamil[6]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[3])*dv12*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += (1.224744871391589*f[13]*hamil[15]+1.224744871391589*f[10]*hamil[14]+1.224744871391589*f[8]*hamil[12]+1.224744871391589*f[6]*hamil[11]+1.224744871391589*f[4]*hamil[9]+1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[1] += (1.224744871391589*f[10]*hamil[15]+1.224744871391589*f[13]*hamil[14]+1.224744871391589*f[4]*hamil[12]+1.224744871391589*f[3]*hamil[11]+1.224744871391589*f[8]*hamil[9]+1.224744871391589*f[6]*hamil[7]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[2] += (1.224744871391589*f[12]*hamil[15]+1.224744871391589*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+1.224744871391589*f[5]*hamil[11]+1.224744871391589*f[4]*hamil[10]+1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[1]*hamil[6]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (1.224744871391589*f[9]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.224744871391589*f[4]*hamil[13]+1.224744871391589*f[2]*hamil[11]+1.224744871391589*f[8]*hamil[10]+1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[0]*hamil[6]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
  out[4] += (1.224744871391589*f[11]*hamil[15]+1.224744871391589*f[7]*hamil[14]+1.224744871391589*f[6]*hamil[13]+1.224744871391589*f[5]*hamil[12]+1.224744871391589*f[3]*hamil[10]+1.224744871391589*f[2]*hamil[9]+1.224744871391589*f[1]*hamil[8]+1.224744871391589*f[0]*hamil[4])*dv12*volFact; 
  out[5] += (1.224744871391589*f[7]*hamil[15]+1.224744871391589*f[11]*hamil[14]+1.224744871391589*f[3]*hamil[13]+1.224744871391589*f[2]*hamil[12]+1.224744871391589*f[6]*hamil[10]+1.224744871391589*f[5]*hamil[9]+1.224744871391589*f[0]*hamil[8]+1.224744871391589*f[1]*hamil[4])*dv12*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (hamil[7]*f[14]+hamil[6]*f[10]+hamil[5]*f[9]+hamil[4]*f[7]+hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[7]*f[15]+hamil[6]*f[13]+hamil[5]*f[12]+hamil[4]*f[11]+hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[14]*hamil[15]+0.7071067811865475*hamil[14]*f[15]+0.7071067811865475*f[10]*hamil[13]+0.7071067811865475*hamil[10]*f[13]+0.7071067811865475*f[9]*hamil[12]+0.7071067811865475*hamil[9]*f[12]+0.7071067811865475*f[7]*hamil[11]+0.7071067811865475*hamil[7]*f[11]+0.7071067811865475*f[4]*hamil[8]+0.7071067811865475*hamil[4]*f[8]+0.7071067811865475*f[3]*hamil[6]+0.7071067811865475*hamil[3]*f[6]+0.7071067811865475*f[2]*hamil[5]+0.7071067811865475*hamil[2]*f[5]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.8284271247461907*f[0]*wx1_sq+1.6329931618554527*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.8284271247461907*f[1]*wx1_sq+1.6329931618554527*f[5]*dv1*wx1+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.8284271247461907*f[0]*wx1*wx2+0.8164965809277261*f[2]*dv1*wx2+0.8164965809277261*f[3]*dv2*wx1+0.2357022603955158*f[7]*dv1*dv2); 
  out[3] += volFact*(2.8284271247461907*f[1]*wx1*wx2+0.8164965809277261*f[5]*dv1*wx2+0.8164965809277261*f[6]*dv2*wx1+0.2357022603955158*f[11]*dv1*dv2); 
  out[4] += volFact*(2.8284271247461907*f[0]*wx1*wx3+0.8164965809277261*f[2]*dv1*wx3+0.8164965809277261*f[4]*dv3*wx1+0.2357022603955158*f[9]*dv1*dv3); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx1*wx3+0.8164965809277261*f[5]*dv1*wx3+0.8164965809277261*f[8]*dv3*wx1+0.2357022603955158*f[12]*dv1*dv3); 
  out[6] += volFact*(2.8284271247461907*f[0]*wx2_sq+1.6329931618554527*f[3]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq); 
  out[7] += volFact*(2.8284271247461907*f[1]*wx2_sq+1.6329931618554527*f[6]*dv2*wx2+0.2357022603955158*f[1]*dv2_sq); 
  out[8] += volFact*(2.8284271247461907*f[0]*wx2*wx3+0.8164965809277261*f[3]*dv2*wx3+0.8164965809277261*f[4]*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3); 
  out[9] += volFact*(2.8284271247461907*f[1]*wx2*wx3+0.8164965809277261*f[6]*dv2*wx3+0.8164965809277261*f[8]*dv3*wx2+0.2357022603955158*f[13]*dv2*dv3); 
  out[10] += volFact*(2.8284271247461907*f[0]*wx3_sq+1.6329931618554527*f[4]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq); 
  out[11] += volFact*(2.8284271247461907*f[1]*wx3_sq+1.6329931618554527*f[8]*dv3*wx3+0.2357022603955158*f[1]*dv3_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += volFact*(2.4494897427831783*f[2]*dv1*wx1_sq+2.8284271247461907*f[0]*wx1_cu+0.7071067811865475*f[0]*dv1_sq*wx1+0.12247448713915887*f[2]*dv1_cu); 
  out[1] += volFact*(2.4494897427831783*f[5]*dv1*wx1_sq+2.8284271247461907*f[1]*wx1_cu+0.7071067811865475*f[1]*dv1_sq*wx1+0.12247448713915887*f[5]*dv1_cu); 
  out[2] += volFact*(2.8284271247461907*f[0]*wx1_sq*wx2+1.6329931618554527*f[2]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[3]*dv2*wx1_sq+0.4714045207910317*f[7]*dv1*dv2*wx1+0.06804138174397717*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.8284271247461907*f[1]*wx1_sq*wx2+1.6329931618554527*f[5]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[6]*dv2*wx1_sq+0.4714045207910317*f[11]*dv1*dv2*wx1+0.06804138174397717*f[6]*dv1_sq*dv2); 
  out[4] += volFact*(2.8284271247461907*f[0]*wx1_sq*wx3+1.6329931618554527*f[2]*dv1*wx1*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[4]*dv3*wx1_sq+0.4714045207910317*f[9]*dv1*dv3*wx1+0.06804138174397717*f[4]*dv1_sq*dv3); 
  out[5] += volFact*(2.8284271247461907*f[1]*wx1_sq*wx3+1.6329931618554527*f[5]*dv1*wx1*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[8]*dv3*wx1_sq+0.4714045207910317*f[12]*dv1*dv3*wx1+0.06804138174397717*f[8]*dv1_sq*dv3); 
  out[6] += volFact*(2.8284271247461907*f[0]*wx1*wx2_sq+0.8164965809277261*f[2]*dv1*wx2_sq+1.6329931618554527*f[3]*dv2*wx1*wx2+0.4714045207910317*f[7]*dv1*dv2*wx2+0.2357022603955158*f[0]*dv2_sq*wx1+0.06804138174397717*f[2]*dv1*dv2_sq); 
  out[7] += volFact*(2.8284271247461907*f[1]*wx1*wx2_sq+0.8164965809277261*f[5]*dv1*wx2_sq+1.6329931618554527*f[6]*dv2*wx1*wx2+0.4714045207910317*f[11]*dv1*dv2*wx2+0.2357022603955158*f[1]*dv2_sq*wx1+0.06804138174397717*f[5]*dv1*dv2_sq); 
  out[8] += volFact*(2.8284271247461907*f[0]*wx1*wx2*wx3+0.8164965809277261*f[2]*dv1*wx2*wx3+0.8164965809277261*f[3]*dv2*wx1*wx3+0.2357022603955158*f[7]*dv1*dv2*wx3+0.8164965809277261*f[4]*dv3*wx1*wx2+0.2357022603955158*f[9]*dv1*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3*wx1+0.06804138174397717*f[14]*dv1*dv2*dv3); 
  out[9] += volFact*(2.8284271247461907*f[1]*wx1*wx2*wx3+0.8164965809277261*f[5]*dv1*wx2*wx3+0.8164965809277261*f[6]*dv2*wx1*wx3+0.2357022603955158*f[11]*dv1*dv2*wx3+0.8164965809277261*f[8]*dv3*wx1*wx2+0.2357022603955158*f[12]*dv1*dv3*wx2+0.2357022603955158*f[13]*dv2*dv3*wx1+0.06804138174397717*f[15]*dv1*dv2*dv3); 
  out[10] += volFact*(2.8284271247461907*f[0]*wx1*wx3_sq+0.8164965809277261*f[2]*dv1*wx3_sq+1.6329931618554527*f[4]*dv3*wx1*wx3+0.4714045207910317*f[9]*dv1*dv3*wx3+0.2357022603955158*f[0]*dv3_sq*wx1+0.06804138174397717*f[2]*dv1*dv3_sq); 
  out[11] += volFact*(2.8284271247461907*f[1]*wx1*wx3_sq+0.8164965809277261*f[5]*dv1*wx3_sq+1.6329931618554527*f[8]*dv3*wx1*wx3+0.4714045207910317*f[12]*dv1*dv3*wx3+0.2357022603955158*f[1]*dv3_sq*wx1+0.06804138174397717*f[5]*dv1*dv3_sq); 
  out[12] += volFact*(2.4494897427831783*f[3]*dv2*wx2_sq+2.8284271247461907*f[0]*wx2_cu+0.7071067811865475*f[0]*dv2_sq*wx2+0.12247448713915887*f[3]*dv2_cu); 
  out[13] += volFact*(2.4494897427831783*f[6]*dv2*wx2_sq+2.8284271247461907*f[1]*wx2_cu+0.7071067811865475*f[1]*dv2_sq*wx2+0.12247448713915887*f[6]*dv2_cu); 
  out[14] += volFact*(2.8284271247461907*f[0]*wx2_sq*wx3+1.6329931618554527*f[3]*dv2*wx2*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.8164965809277261*f[4]*dv3*wx2_sq+0.4714045207910317*f[10]*dv2*dv3*wx2+0.06804138174397717*f[4]*dv2_sq*dv3); 
  out[15] += volFact*(2.8284271247461907*f[1]*wx2_sq*wx3+1.6329931618554527*f[6]*dv2*wx2*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.8164965809277261*f[8]*dv3*wx2_sq+0.4714045207910317*f[13]*dv2*dv3*wx2+0.06804138174397717*f[8]*dv2_sq*dv3); 
  out[16] += volFact*(2.8284271247461907*f[0]*wx2*wx3_sq+0.8164965809277261*f[3]*dv2*wx3_sq+1.6329931618554527*f[4]*dv3*wx2*wx3+0.4714045207910317*f[10]*dv2*dv3*wx3+0.2357022603955158*f[0]*dv3_sq*wx2+0.06804138174397717*f[3]*dv2*dv3_sq); 
  out[17] += volFact*(2.8284271247461907*f[1]*wx2*wx3_sq+0.8164965809277261*f[6]*dv2*wx3_sq+1.6329931618554527*f[8]*dv3*wx2*wx3+0.4714045207910317*f[13]*dv2*dv3*wx3+0.2357022603955158*f[1]*dv3_sq*wx2+0.06804138174397717*f[6]*dv2*dv3_sq); 
  out[18] += volFact*(2.4494897427831783*f[4]*dv3*wx3_sq+2.8284271247461907*f[0]*wx3_cu+0.7071067811865475*f[0]*dv3_sq*wx3+0.12247448713915887*f[4]*dv3_cu); 
  out[19] += volFact*(2.4494897427831783*f[8]*dv3*wx3_sq+2.8284271247461907*f[1]*wx3_cu+0.7071067811865475*f[1]*dv3_sq*wx3+0.12247448713915887*f[8]*dv3_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += (1.7320508075688772*hamil[7]*f[10]+1.7320508075688772*f[4]*hamil[5]+1.7320508075688772*f[3]*hamil[4]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[3] += (1.7320508075688772*hamil[7]*f[13]+1.7320508075688772*hamil[5]*f[8]+1.7320508075688772*hamil[4]*f[6]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[4] += (1.7320508075688772*hamil[7]*f[9]+1.7320508075688772*f[4]*hamil[6]+1.7320508075688772*f[2]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[5] += (1.7320508075688772*hamil[7]*f[12]+1.7320508075688772*hamil[6]*f[8]+1.7320508075688772*hamil[4]*f[5]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
  out[6] += (1.7320508075688772*f[7]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[2]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv12*volFact; 
  out[7] += (1.7320508075688772*hamil[7]*f[11]+1.7320508075688772*f[6]*hamil[6]+1.7320508075688772*f[5]*hamil[5]+1.7320508075688772*f[1]*hamil[3])*dv12*volFact; 
  out[8] += (hamil[7]*f[14]+hamil[6]*f[10]+hamil[5]*f[9]+hamil[4]*f[7]+hamil[3]*f[4]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[9] += (hamil[7]*f[15]+hamil[6]*f[13]+hamil[5]*f[12]+hamil[4]*f[11]+hamil[3]*f[8]+hamil[2]*f[6]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += (1.224744871391589*f[13]*hamil[15]+1.224744871391589*f[10]*hamil[14]+1.224744871391589*f[8]*hamil[12]+1.224744871391589*f[6]*hamil[11]+1.224744871391589*f[4]*hamil[9]+1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[3] += (1.224744871391589*f[10]*hamil[15]+1.224744871391589*f[13]*hamil[14]+1.224744871391589*f[4]*hamil[12]+1.224744871391589*f[3]*hamil[11]+1.224744871391589*f[8]*hamil[9]+1.224744871391589*f[6]*hamil[7]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[4] += (1.224744871391589*f[12]*hamil[15]+1.224744871391589*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+1.224744871391589*f[5]*hamil[11]+1.224744871391589*f[4]*hamil[10]+1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[1]*hamil[6]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[5] += (1.224744871391589*f[9]*hamil[15]+1.224744871391589*f[12]*hamil[14]+1.224744871391589*f[4]*hamil[13]+1.224744871391589*f[2]*hamil[11]+1.224744871391589*f[8]*hamil[10]+1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[0]*hamil[6]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
  out[6] += (1.224744871391589*f[11]*hamil[15]+1.224744871391589*f[7]*hamil[14]+1.224744871391589*f[6]*hamil[13]+1.224744871391589*f[5]*hamil[12]+1.224744871391589*f[3]*hamil[10]+1.224744871391589*f[2]*hamil[9]+1.224744871391589*f[1]*hamil[8]+1.224744871391589*f[0]*hamil[4])*dv12*volFact; 
  out[7] += (1.224744871391589*f[7]*hamil[15]+1.224744871391589*f[11]*hamil[14]+1.224744871391589*f[3]*hamil[13]+1.224744871391589*f[2]*hamil[12]+1.224744871391589*f[6]*hamil[10]+1.224744871391589*f[5]*hamil[9]+1.224744871391589*f[0]*hamil[8]+1.224744871391589*f[1]*hamil[4])*dv12*volFact; 
  out[8] += (0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[9] += (0.7071067811865475*f[14]*hamil[15]+0.7071067811865475*hamil[14]*f[15]+0.7071067811865475*f[10]*hamil[13]+0.7071067811865475*hamil[10]*f[13]+0.7071067811865475*f[9]*hamil[12]+0.7071067811865475*hamil[9]*f[12]+0.7071067811865475*f[7]*hamil[11]+0.7071067811865475*hamil[7]*f[11]+0.7071067811865475*f[4]*hamil[8]+0.7071067811865475*hamil[4]*f[8]+0.7071067811865475*f[3]*hamil[6]+0.7071067811865475*hamil[3]*f[6]+0.7071067811865475*f[2]*hamil[5]+0.7071067811865475*hamil[2]*f[5]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.4494897427831783*hamil[7]*f[10]+2.4494897427831783*f[4]*hamil[5]+2.4494897427831783*f[3]*hamil[4]+2.4494897427831783*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (2.4494897427831783*hamil[7]*f[9]+2.4494897427831783*f[4]*hamil[6]+2.4494897427831783*f[2]*hamil[4]+2.4494897427831783*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (2.4494897427831783*f[7]*hamil[7]+2.4494897427831783*f[3]*hamil[6]+2.4494897427831783*f[2]*hamil[5]+2.4494897427831783*f[0]*hamil[3])*dv12*volFact; 
  out[4] += (1.4142135623730951*hamil[7]*f[14]+1.4142135623730951*hamil[6]*f[10]+1.4142135623730951*hamil[5]*f[9]+1.4142135623730951*hamil[4]*f[7]+1.4142135623730951*hamil[3]*f[4]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[13]*hamil[15]+1.7320508075688772*f[10]*hamil[14]+1.7320508075688772*f[8]*hamil[12]+1.7320508075688772*f[6]*hamil[11]+1.7320508075688772*f[4]*hamil[9]+1.7320508075688772*f[3]*hamil[7]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[2])*dv10*volFact; 
  out[2] += (1.7320508075688772*f[12]*hamil[15]+1.7320508075688772*f[9]*hamil[14]+1.7320508075688772*f[8]*hamil[13]+1.7320508075688772*f[5]*hamil[11]+1.7320508075688772*f[4]*hamil[10]+1.7320508075688772*f[2]*hamil[7]+1.7320508075688772*f[1]*hamil[6]+1.7320508075688772*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (1.7320508075688772*f[11]*hamil[15]+1.7320508075688772*f[7]*hamil[14]+1.7320508075688772*f[6]*hamil[13]+1.7320508075688772*f[5]*hamil[12]+1.7320508075688772*f[3]*hamil[10]+1.7320508075688772*f[2]*hamil[9]+1.7320508075688772*f[1]*hamil[8]+1.7320508075688772*f[0]*hamil[4])*dv12*volFact; 
  out[4] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
