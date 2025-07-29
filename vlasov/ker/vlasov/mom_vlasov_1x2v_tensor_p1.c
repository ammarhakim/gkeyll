#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (1.7320508075688772*hamil[3]*f[5]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += (1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[1] += (1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[2] += (1.224744871391589*f[4]*hamil[7]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[7]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  out[0] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.1547005383792517*f[2]*dv1*wx1+0.16666666666666666*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.1547005383792517*f[4]*dv1*wx1+0.16666666666666666*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[2]*dv1*wx2+0.5773502691896258*f[3]*dv2*wx1+0.16666666666666666*f[6]*dv1*dv2); 
  out[3] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[4]*dv1*wx2+0.5773502691896258*f[5]*dv2*wx1+0.16666666666666666*f[7]*dv1*dv2); 
  out[4] += volFact*(2.0*f[0]*wx2_sq+1.1547005383792517*f[3]*dv2*wx2+0.16666666666666666*f[0]*dv2_sq); 
  out[5] += volFact*(2.0*f[1]*wx2_sq+1.1547005383792517*f[5]*dv2*wx2+0.16666666666666666*f[1]*dv2_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(1.7320508075688772*f[2]*dv1*wx1_sq+2.0*f[0]*wx1_cu+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[2]*dv1_cu); 
  out[1] += volFact*(1.7320508075688772*f[4]*dv1*wx1_sq+2.0*f[1]*wx1_cu+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[4]*dv1_cu); 
  out[2] += volFact*(2.0*f[0]*wx1_sq*wx2+1.1547005383792517*f[2]*dv1*wx1*wx2+0.16666666666666666*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.0*f[1]*wx1_sq*wx2+1.1547005383792517*f[4]*dv1*wx1*wx2+0.16666666666666666*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[7]*dv1*dv2*wx1+0.04811252243246882*f[5]*dv1_sq*dv2); 
  out[4] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.1547005383792517*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+0.16666666666666666*f[0]*dv2_sq*wx1+0.04811252243246882*f[2]*dv1*dv2_sq); 
  out[5] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.1547005383792517*f[5]*dv2*wx1*wx2+0.3333333333333333*f[7]*dv1*dv2*wx2+0.16666666666666666*f[1]*dv2_sq*wx1+0.04811252243246882*f[4]*dv1*dv2_sq); 
  out[6] += volFact*(1.7320508075688772*f[3]*dv2*wx2_sq+2.0*f[0]*wx2_cu+0.5*f[0]*dv2_sq*wx2+0.08660254037844387*f[3]*dv2_cu); 
  out[7] += volFact*(1.7320508075688772*f[5]*dv2*wx2_sq+2.0*f[1]*wx2_cu+0.5*f[1]*dv2_sq*wx2+0.08660254037844387*f[5]*dv2_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[5]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[4] += (1.7320508075688772*f[2]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
  out[6] += (hamil[3]*f[6]+hamil[2]*f[3]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[7] += (hamil[3]*f[7]+hamil[2]*f[5]+hamil[1]*f[4]+hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.224744871391589*f[5]*hamil[7]+1.224744871391589*f[3]*hamil[6]+1.224744871391589*f[1]*hamil[4]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[3] += (1.224744871391589*f[3]*hamil[7]+1.224744871391589*f[5]*hamil[6]+1.224744871391589*f[0]*hamil[4]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[4] += (1.224744871391589*f[4]*hamil[7]+1.224744871391589*f[2]*hamil[6]+1.224744871391589*f[1]*hamil[5]+1.224744871391589*f[0]*hamil[3])*dv11*volFact; 
  out[5] += (1.224744871391589*f[2]*hamil[7]+1.224744871391589*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[5]+1.224744871391589*f[1]*hamil[3])*dv11*volFact; 
  out[6] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[7] += (0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.4494897427831783*f[3]*hamil[3]+2.4494897427831783*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (2.4494897427831783*f[2]*hamil[3]+2.4494897427831783*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (1.4142135623730951*hamil[3]*f[6]+1.4142135623730951*hamil[2]*f[3]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_1x2v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[5]*hamil[7]+1.7320508075688772*f[3]*hamil[6]+1.7320508075688772*f[1]*hamil[4]+1.7320508075688772*f[0]*hamil[2])*dv10*volFact; 
  out[2] += (1.7320508075688772*f[4]*hamil[7]+1.7320508075688772*f[2]*hamil[6]+1.7320508075688772*f[1]*hamil[5]+1.7320508075688772*f[0]*hamil[3])*dv11*volFact; 
  out[3] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
