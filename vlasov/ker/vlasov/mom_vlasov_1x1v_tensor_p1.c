#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[4]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.224744871391589*f[2]*dv1*wx1_sq+1.414213562373095*f[0]*wx1_cu+0.3162277660168379*f[4]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1_cu); 
  out[1] += volFact*(1.224744871391589*f[3]*dv1*wx1_sq+1.414213562373095*f[1]*wx1_cu+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += (3.872983346207417*f[2]*hamil[2]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (3.872983346207417*hamil[2]*f[3]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (hamil[2]*f[4]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0*hamil[2]*f[5]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += (3.872983346207417*f[2]*hamil[2]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[3] += (3.872983346207417*hamil[2]*f[3]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
  out[4] += (hamil[2]*f[4]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[5] += (1.0*hamil[2]*f[5]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (5.477225575051662*f[2]*hamil[2]+2.449489742783178*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (1.414213562373095*hamil[2]*f[4]+1.414213562373095*hamil[1]*f[2]+1.414213562373095*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += (2.738612787525831*f[3]*hamil[5]+2.738612787525831*f[2]*hamil[4]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[1] += (2.738612787525831*f[2]*hamil[5]+2.738612787525831*f[3]*hamil[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[4]*hamil[5]+0.7071067811865475*hamil[4]*f[5]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += (2.738612787525831*f[3]*hamil[5]+2.738612787525831*f[2]*hamil[4]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[3] += (2.738612787525831*f[2]*hamil[5]+2.738612787525831*f[3]*hamil[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[4] += (0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[5] += (0.7071067811865475*f[4]*hamil[5]+0.7071067811865475*hamil[4]*f[5]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_1x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (3.872983346207417*f[3]*hamil[5]+3.872983346207417*f[2]*hamil[4]+1.732050807568877*f[1]*hamil[3]+1.732050807568877*f[0]*hamil[2])*dv10*volFact; 
  out[2] += (f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
