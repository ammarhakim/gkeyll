#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx2 = w[2], dv2 = dxv[2]; 
  double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[7]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[10]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[2]*dv1*wx2+0.5773502691896258*f[3]*dv2*wx1+0.1666666666666667*f[6]*dv1*dv2); 
  out[3] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[4]*dv1*wx2+0.5773502691896258*f[5]*dv2*wx1+0.1666666666666667*f[9]*dv1*dv2); 
  out[4] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+0.149071198499986*f[8]*dv2_sq+0.1666666666666667*f[0]*dv2_sq); 
  out[5] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+0.149071198499986*f[12]*dv2_sq+0.1666666666666667*f[1]*dv2_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double wx2 = w[2], dv2 = dxv[2]; 
  double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(1.732050807568877*f[2]*dv1*wx1_sq+2.0*f[0]*wx1_cu+0.4472135954999579*f[7]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[2]*dv1_cu); 
  out[1] += volFact*(1.732050807568877*f[4]*dv1*wx1_sq+2.0*f[1]*wx1_cu+0.447213595499958*f[10]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[4]*dv1_cu); 
  out[2] += volFact*(2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[2]*dv1*wx1*wx2+0.149071198499986*f[7]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.04303314829119351*f[11]*dv1_sq*dv2+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[4]*dv1*wx1*wx2+0.149071198499986*f[10]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[9]*dv1*dv2*wx1+0.04303314829119353*f[14]*dv1_sq*dv2+0.04811252243246882*f[5]*dv1_sq*dv2); 
  out[4] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.154700538379252*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+0.149071198499986*f[8]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.04303314829119351*f[13]*dv1*dv2_sq+0.04811252243246882*f[2]*dv1*dv2_sq); 
  out[5] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.154700538379252*f[5]*dv2*wx1*wx2+0.3333333333333333*f[9]*dv1*dv2*wx2+0.149071198499986*f[12]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.04303314829119353*f[15]*dv1*dv2_sq+0.04811252243246882*f[4]*dv1*dv2_sq); 
  out[6] += volFact*(1.732050807568877*f[3]*dv2*wx2_sq+2.0*f[0]*wx2_cu+0.4472135954999579*f[8]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.08660254037844387*f[3]*dv2_cu); 
  out[7] += volFact*(1.732050807568877*f[5]*dv2*wx2_sq+2.0*f[1]*wx2_cu+0.447213595499958*f[12]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.08660254037844387*f[5]*dv2_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += (3.872983346207417*hamil[8]*f[13]+1.732050807568877*hamil[7]*f[8]+3.872983346207417*f[6]*hamil[6]+3.872983346207417*f[2]*hamil[4]+1.732050807568877*f[3]*hamil[3]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (3.872983346207417*hamil[8]*f[15]+1.732050807568877*hamil[7]*f[12]+3.872983346207417*hamil[6]*f[9]+1.732050807568877*hamil[3]*f[5]+3.872983346207417*f[4]*hamil[4]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (3.872983346207417*hamil[8]*f[11]+3.872983346207417*f[6]*hamil[7]+1.732050807568877*hamil[6]*f[7]+3.872983346207417*f[3]*hamil[5]+1.732050807568877*f[2]*hamil[3]+1.732050807568877*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (3.872983346207417*hamil[8]*f[14]+1.732050807568877*hamil[6]*f[10]+3.872983346207417*hamil[7]*f[9]+3.872983346207417*f[5]*hamil[5]+1.732050807568877*hamil[3]*f[4]+1.732050807568877*f[1]*hamil[2])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[8]*f[16]+hamil[7]*f[13]+hamil[6]*f[11]+hamil[5]*f[8]+hamil[4]*f[7]+hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[8]*f[17]+1.0*hamil[7]*f[15]+1.0*hamil[6]*f[14]+1.0*hamil[5]*f[12]+1.0*hamil[4]*f[10]+hamil[3]*f[9]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (3.872983346207417*hamil[8]*f[13]+1.732050807568877*hamil[7]*f[8]+3.872983346207417*f[6]*hamil[6]+3.872983346207417*f[2]*hamil[4]+1.732050807568877*f[3]*hamil[3]+1.732050807568877*f[0]*hamil[1])*dv10*volFact; 
  out[3] += (3.872983346207417*hamil[8]*f[15]+1.732050807568877*hamil[7]*f[12]+3.872983346207417*hamil[6]*f[9]+1.732050807568877*hamil[3]*f[5]+3.872983346207417*f[4]*hamil[4]+1.732050807568877*f[1]*hamil[1])*dv10*volFact; 
  out[4] += (3.872983346207417*hamil[8]*f[11]+3.872983346207417*f[6]*hamil[7]+1.732050807568877*hamil[6]*f[7]+3.872983346207417*f[3]*hamil[5]+1.732050807568877*f[2]*hamil[3]+1.732050807568877*f[0]*hamil[2])*dv11*volFact; 
  out[5] += (3.872983346207417*hamil[8]*f[14]+1.732050807568877*hamil[6]*f[10]+3.872983346207417*hamil[7]*f[9]+3.872983346207417*f[5]*hamil[5]+1.732050807568877*hamil[3]*f[4]+1.732050807568877*f[1]*hamil[2])*dv11*volFact; 
  out[6] += (hamil[8]*f[16]+hamil[7]*f[13]+hamil[6]*f[11]+hamil[5]*f[8]+hamil[4]*f[7]+hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (hamil[8]*f[17]+1.0*hamil[7]*f[15]+1.0*hamil[6]*f[14]+1.0*hamil[5]*f[12]+1.0*hamil[4]*f[10]+hamil[3]*f[9]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (5.477225575051662*hamil[8]*f[13]+2.449489742783178*hamil[7]*f[8]+5.477225575051662*f[6]*hamil[6]+5.477225575051662*f[2]*hamil[4]+2.449489742783178*f[3]*hamil[3]+2.449489742783178*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (5.477225575051662*hamil[8]*f[11]+5.477225575051662*f[6]*hamil[7]+2.449489742783178*hamil[6]*f[7]+5.477225575051662*f[3]*hamil[5]+2.449489742783178*f[2]*hamil[3]+2.449489742783178*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (1.414213562373095*hamil[8]*f[16]+1.414213562373095*hamil[7]*f[13]+1.414213562373095*hamil[6]*f[11]+1.414213562373095*hamil[5]*f[8]+1.414213562373095*hamil[4]*f[7]+1.414213562373095*hamil[3]*f[6]+1.414213562373095*hamil[2]*f[3]+1.414213562373095*hamil[1]*f[2]+1.414213562373095*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += (2.738612787525831*f[15]*hamil[17]+2.738612787525831*f[13]*hamil[16]+1.224744871391589*f[12]*hamil[15]+2.738612787525831*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+2.738612787525831*f[6]*hamil[11]+2.738612787525831*f[4]*hamil[10]+1.224744871391589*f[5]*hamil[9]+2.738612787525831*f[2]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[1] += (2.738612787525831*f[13]*hamil[17]+2.738612787525831*f[15]*hamil[16]+1.224744871391589*f[8]*hamil[15]+2.738612787525831*f[6]*hamil[14]+1.224744871391589*f[12]*hamil[13]+2.738612787525831*f[9]*hamil[11]+2.738612787525831*f[2]*hamil[10]+1.224744871391589*f[3]*hamil[9]+2.738612787525831*f[4]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[2] += (2.738612787525831*f[14]*hamil[17]+2.738612787525831*f[11]*hamil[16]+2.738612787525831*f[9]*hamil[15]+1.224744871391589*f[10]*hamil[14]+2.738612787525831*f[6]*hamil[13]+2.738612787525831*f[5]*hamil[12]+1.224744871391589*f[7]*hamil[11]+1.224744871391589*f[4]*hamil[9]+2.738612787525831*f[3]*hamil[8]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (2.738612787525831*f[11]*hamil[17]+2.738612787525831*f[14]*hamil[16]+2.738612787525831*f[6]*hamil[15]+1.224744871391589*f[7]*hamil[14]+2.738612787525831*f[9]*hamil[13]+2.738612787525831*f[3]*hamil[12]+1.224744871391589*f[10]*hamil[11]+1.224744871391589*f[2]*hamil[9]+2.738612787525831*f[5]*hamil[8]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (0.7071067811865475*f[17]*hamil[17]+0.7071067811865475*f[16]*hamil[16]+0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[16]*hamil[17]+0.7071067811865475*hamil[16]*f[17]+0.7071067811865475*f[13]*hamil[15]+0.7071067811865475*hamil[13]*f[15]+0.7071067811865475*f[11]*hamil[14]+0.7071067811865475*hamil[11]*f[14]+0.7071067811865475*f[8]*hamil[12]+0.7071067811865475*hamil[8]*f[12]+0.7071067811865475*f[7]*hamil[10]+0.7071067811865475*hamil[7]*f[10]+0.7071067811865475*f[6]*hamil[9]+0.7071067811865475*hamil[6]*f[9]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (2.738612787525831*f[15]*hamil[17]+2.738612787525831*f[13]*hamil[16]+1.224744871391589*f[12]*hamil[15]+2.738612787525831*f[9]*hamil[14]+1.224744871391589*f[8]*hamil[13]+2.738612787525831*f[6]*hamil[11]+2.738612787525831*f[4]*hamil[10]+1.224744871391589*f[5]*hamil[9]+2.738612787525831*f[2]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[3] += (2.738612787525831*f[13]*hamil[17]+2.738612787525831*f[15]*hamil[16]+1.224744871391589*f[8]*hamil[15]+2.738612787525831*f[6]*hamil[14]+1.224744871391589*f[12]*hamil[13]+2.738612787525831*f[9]*hamil[11]+2.738612787525831*f[2]*hamil[10]+1.224744871391589*f[3]*hamil[9]+2.738612787525831*f[4]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[4] += (2.738612787525831*f[14]*hamil[17]+2.738612787525831*f[11]*hamil[16]+2.738612787525831*f[9]*hamil[15]+1.224744871391589*f[10]*hamil[14]+2.738612787525831*f[6]*hamil[13]+2.738612787525831*f[5]*hamil[12]+1.224744871391589*f[7]*hamil[11]+1.224744871391589*f[4]*hamil[9]+2.738612787525831*f[3]*hamil[8]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[5] += (2.738612787525831*f[11]*hamil[17]+2.738612787525831*f[14]*hamil[16]+2.738612787525831*f[6]*hamil[15]+1.224744871391589*f[7]*hamil[14]+2.738612787525831*f[9]*hamil[13]+2.738612787525831*f[3]*hamil[12]+1.224744871391589*f[10]*hamil[11]+1.224744871391589*f[2]*hamil[9]+2.738612787525831*f[5]*hamil[8]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
  out[6] += (0.7071067811865475*f[17]*hamil[17]+0.7071067811865475*f[16]*hamil[16]+0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[7] += (0.7071067811865475*f[16]*hamil[17]+0.7071067811865475*hamil[16]*f[17]+0.7071067811865475*f[13]*hamil[15]+0.7071067811865475*hamil[13]*f[15]+0.7071067811865475*f[11]*hamil[14]+0.7071067811865475*hamil[11]*f[14]+0.7071067811865475*f[8]*hamil[12]+0.7071067811865475*hamil[8]*f[12]+0.7071067811865475*f[7]*hamil[10]+0.7071067811865475*hamil[7]*f[10]+0.7071067811865475*f[6]*hamil[9]+0.7071067811865475*hamil[6]*f[9]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (3.872983346207417*f[15]*hamil[17]+3.872983346207417*f[13]*hamil[16]+1.732050807568877*f[12]*hamil[15]+3.872983346207417*f[9]*hamil[14]+1.732050807568877*f[8]*hamil[13]+3.872983346207417*f[6]*hamil[11]+3.872983346207417*f[4]*hamil[10]+1.732050807568877*f[5]*hamil[9]+3.872983346207417*f[2]*hamil[7]+1.732050807568877*f[3]*hamil[6]+1.732050807568877*f[1]*hamil[4]+1.732050807568877*f[0]*hamil[2])*dv10*volFact; 
  out[2] += (3.872983346207417*f[14]*hamil[17]+3.872983346207417*f[11]*hamil[16]+3.872983346207417*f[9]*hamil[15]+1.732050807568877*f[10]*hamil[14]+3.872983346207417*f[6]*hamil[13]+3.872983346207417*f[5]*hamil[12]+1.732050807568877*f[7]*hamil[11]+1.732050807568877*f[4]*hamil[9]+3.872983346207417*f[3]*hamil[8]+1.732050807568877*f[2]*hamil[6]+1.732050807568877*f[1]*hamil[5]+1.732050807568877*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (f[17]*hamil[17]+f[16]*hamil[16]+f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
