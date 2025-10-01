#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double wx1 = w[2], dv1 = dxv[2]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[5]*dv1*wx1+0.105409255338946*f[9]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[2]*wx1_sq+0.8164965809277261*f[6]*dv1*wx1+0.105409255338946*f[10]*dv1_sq+0.1178511301977579*f[2]*dv1_sq); 
  out[3] += volFact*(1.414213562373095*f[4]*wx1_sq+0.8164965809277261*f[8]*dv1*wx1+0.105409255338946*f[11]*dv1_sq+0.1178511301977579*f[4]*dv1_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double wx1 = w[2], dv1 = dxv[2]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.224744871391589*f[3]*dv1*wx1_sq+1.414213562373095*f[0]*wx1_cu+0.3162277660168379*f[7]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1_cu); 
  out[1] += volFact*(1.224744871391589*f[5]*dv1*wx1_sq+1.414213562373095*f[1]*wx1_cu+0.3162277660168379*f[9]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[5]*dv1_cu); 
  out[2] += volFact*(1.224744871391589*f[6]*dv1*wx1_sq+1.414213562373095*f[2]*wx1_cu+0.3162277660168379*f[10]*dv1_sq*wx1+0.3535533905932737*f[2]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1_cu); 
  out[3] += volFact*(1.224744871391589*f[8]*dv1*wx1_sq+1.414213562373095*f[4]*wx1_cu+0.3162277660168379*f[11]*dv1_sq*wx1+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[8]*dv1_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += (3.872983346207417*hamil[2]*f[3]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (3.872983346207417*hamil[2]*f[5]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (3.872983346207417*hamil[2]*f[6]+1.732050807568877*hamil[1]*f[2])*dv10*volFact; 
  out[3] += (3.872983346207417*hamil[2]*f[8]+1.732050807568877*hamil[1]*f[4])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (hamil[2]*f[7]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0*hamil[2]*f[9]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0*hamil[2]*f[10]+hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[2]*f[11]+hamil[1]*f[8]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += (3.872983346207417*hamil[2]*f[3]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[5] += (3.872983346207417*hamil[2]*f[5]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
  out[6] += (3.872983346207417*hamil[2]*f[6]+1.732050807568877*hamil[1]*f[2])*dv10*volFact; 
  out[7] += (3.872983346207417*hamil[2]*f[8]+1.732050807568877*hamil[1]*f[4])*dv10*volFact; 
  out[8] += (hamil[2]*f[7]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[9] += (1.0*hamil[2]*f[9]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[10] += (1.0*hamil[2]*f[10]+hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[11] += (hamil[2]*f[11]+hamil[1]*f[8]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (7.745966692414834*hamil[2]*f[3]+3.464101615137754*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (2.0*hamil[2]*f[7]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += (1.936491673103709*f[8]*hamil[11]+1.936491673103709*f[6]*hamil[10]+1.936491673103709*f[5]*hamil[9]+0.8660254037844386*f[4]*hamil[8]+1.936491673103709*f[3]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*f[1]*hamil[5]+0.8660254037844386*f[0]*hamil[3])*dv10*volFact; 
  out[1] += (1.936491673103709*f[6]*hamil[11]+1.936491673103709*f[8]*hamil[10]+1.936491673103709*f[3]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+1.936491673103709*f[5]*hamil[7]+0.8660254037844386*f[4]*hamil[6]+0.8660254037844386*f[0]*hamil[5]+0.8660254037844386*f[1]*hamil[3])*dv10*volFact; 
  out[2] += (1.936491673103709*f[5]*hamil[11]+1.936491673103709*f[3]*hamil[10]+1.936491673103709*f[8]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+1.936491673103709*f[6]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[4]*hamil[5]+0.8660254037844386*f[2]*hamil[3])*dv10*volFact; 
  out[3] += (1.936491673103709*f[3]*hamil[11]+1.936491673103709*f[5]*hamil[10]+1.936491673103709*f[6]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+1.936491673103709*hamil[7]*f[8]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[2]*hamil[5]+0.8660254037844386*hamil[3]*f[4])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5000000000000001*f[10]*hamil[11]+0.5000000000000001*hamil[10]*f[11]+0.5000000000000001*f[7]*hamil[9]+0.5000000000000001*hamil[7]*f[9]+0.5*f[6]*hamil[8]+0.5*hamil[6]*f[8]+0.5*f[3]*hamil[5]+0.5*hamil[3]*f[5]+0.5*f[2]*hamil[4]+0.5*hamil[2]*f[4]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5000000000000001*f[9]*hamil[11]+0.5000000000000001*hamil[9]*f[11]+0.5000000000000001*f[7]*hamil[10]+0.5000000000000001*hamil[7]*f[10]+0.5*f[5]*hamil[8]+0.5*hamil[5]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[1]*hamil[4]+0.5*hamil[1]*f[4]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[9]*hamil[10]+0.5*hamil[9]*f[10]+0.5*f[3]*hamil[8]+0.5*hamil[3]*f[8]+0.5*f[5]*hamil[6]+0.5*hamil[5]*f[6]+0.5*f[0]*hamil[4]+0.5*hamil[0]*f[4]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += (1.936491673103709*f[8]*hamil[11]+1.936491673103709*f[6]*hamil[10]+1.936491673103709*f[5]*hamil[9]+0.8660254037844386*f[4]*hamil[8]+1.936491673103709*f[3]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*f[1]*hamil[5]+0.8660254037844386*f[0]*hamil[3])*dv10*volFact; 
  out[5] += (1.936491673103709*f[6]*hamil[11]+1.936491673103709*f[8]*hamil[10]+1.936491673103709*f[3]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+1.936491673103709*f[5]*hamil[7]+0.8660254037844386*f[4]*hamil[6]+0.8660254037844386*f[0]*hamil[5]+0.8660254037844386*f[1]*hamil[3])*dv10*volFact; 
  out[6] += (1.936491673103709*f[5]*hamil[11]+1.936491673103709*f[3]*hamil[10]+1.936491673103709*f[8]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+1.936491673103709*f[6]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[4]*hamil[5]+0.8660254037844386*f[2]*hamil[3])*dv10*volFact; 
  out[7] += (1.936491673103709*f[3]*hamil[11]+1.936491673103709*f[5]*hamil[10]+1.936491673103709*f[6]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+1.936491673103709*hamil[7]*f[8]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[2]*hamil[5]+0.8660254037844386*hamil[3]*f[4])*dv10*volFact; 
  out[8] += (0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[9] += (0.5000000000000001*f[10]*hamil[11]+0.5000000000000001*hamil[10]*f[11]+0.5000000000000001*f[7]*hamil[9]+0.5000000000000001*hamil[7]*f[9]+0.5*f[6]*hamil[8]+0.5*hamil[6]*f[8]+0.5*f[3]*hamil[5]+0.5*hamil[3]*f[5]+0.5*f[2]*hamil[4]+0.5*hamil[2]*f[4]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[10] += (0.5000000000000001*f[9]*hamil[11]+0.5000000000000001*hamil[9]*f[11]+0.5000000000000001*f[7]*hamil[10]+0.5000000000000001*hamil[7]*f[10]+0.5*f[5]*hamil[8]+0.5*hamil[5]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[1]*hamil[4]+0.5*hamil[1]*f[4]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[11] += (0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[9]*hamil[10]+0.5*hamil[9]*f[10]+0.5*f[3]*hamil[8]+0.5*hamil[3]*f[8]+0.5*f[5]*hamil[6]+0.5*hamil[5]*f[6]+0.5*f[0]*hamil[4]+0.5*hamil[0]*f[4]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (3.872983346207417*f[8]*hamil[11]+3.872983346207417*f[6]*hamil[10]+3.872983346207417*f[5]*hamil[9]+1.732050807568877*f[4]*hamil[8]+3.872983346207417*f[3]*hamil[7]+1.732050807568877*f[2]*hamil[6]+1.732050807568877*f[1]*hamil[5]+1.732050807568877*f[0]*hamil[3])*dv10*volFact; 
  out[2] += (f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
