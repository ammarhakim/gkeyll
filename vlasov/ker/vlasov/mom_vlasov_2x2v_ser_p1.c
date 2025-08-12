#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double wx1 = w[2], dv1 = dxv[2]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx2 = w[3], dv2 = dxv[3]; 
  double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.1547005383792517*f[3]*dv1*wx1+0.16666666666666666*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.1547005383792517*f[6]*dv1*wx1+0.16666666666666666*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.1547005383792517*f[7]*dv1*wx1+0.16666666666666666*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.1547005383792517*f[11]*dv1*wx1+0.16666666666666666*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[3]*dv1*wx2+0.5773502691896258*f[4]*dv2*wx1+0.16666666666666666*f[10]*dv1*dv2); 
  out[5] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[6]*dv1*wx2+0.5773502691896258*f[8]*dv2*wx1+0.16666666666666666*f[13]*dv1*dv2); 
  out[6] += volFact*(2.0*f[2]*wx1*wx2+0.5773502691896258*f[7]*dv1*wx2+0.5773502691896258*f[9]*dv2*wx1+0.16666666666666666*f[14]*dv1*dv2); 
  out[7] += volFact*(2.0*f[5]*wx1*wx2+0.5773502691896258*f[11]*dv1*wx2+0.5773502691896258*f[12]*dv2*wx1+0.16666666666666666*f[15]*dv1*dv2); 
  out[8] += volFact*(2.0*f[0]*wx2_sq+1.1547005383792517*f[4]*dv2*wx2+0.16666666666666666*f[0]*dv2_sq); 
  out[9] += volFact*(2.0*f[1]*wx2_sq+1.1547005383792517*f[8]*dv2*wx2+0.16666666666666666*f[1]*dv2_sq); 
  out[10] += volFact*(2.0*f[2]*wx2_sq+1.1547005383792517*f[9]*dv2*wx2+0.16666666666666666*f[2]*dv2_sq); 
  out[11] += volFact*(2.0*f[5]*wx2_sq+1.1547005383792517*f[12]*dv2*wx2+0.16666666666666666*f[5]*dv2_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
  const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double wx1 = w[2], dv1 = dxv[2]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double wx2 = w[3], dv2 = dxv[3]; 
  double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(1.7320508075688772*f[3]*dv1*wx1_sq+2.0*f[0]*wx1_cu+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1_cu); 
  out[1] += volFact*(1.7320508075688772*f[6]*dv1*wx1_sq+2.0*f[1]*wx1_cu+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1_cu); 
  out[2] += volFact*(1.7320508075688772*f[7]*dv1*wx1_sq+2.0*f[2]*wx1_cu+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1_cu); 
  out[3] += volFact*(1.7320508075688772*f[11]*dv1*wx1_sq+2.0*f[5]*wx1_cu+0.5*f[5]*dv1_sq*wx1+0.08660254037844387*f[11]*dv1_cu); 
  out[4] += volFact*(2.0*f[0]*wx1_sq*wx2+1.1547005383792517*f[3]*dv1*wx1*wx2+0.16666666666666666*f[0]*dv1_sq*wx2+0.5773502691896258*f[4]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.04811252243246882*f[4]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[1]*wx1_sq*wx2+1.1547005383792517*f[6]*dv1*wx1*wx2+0.16666666666666666*f[1]*dv1_sq*wx2+0.5773502691896258*f[8]*dv2*wx1_sq+0.3333333333333333*f[13]*dv1*dv2*wx1+0.04811252243246882*f[8]*dv1_sq*dv2); 
  out[6] += volFact*(2.0*f[2]*wx1_sq*wx2+1.1547005383792517*f[7]*dv1*wx1*wx2+0.16666666666666666*f[2]*dv1_sq*wx2+0.5773502691896258*f[9]*dv2*wx1_sq+0.3333333333333333*f[14]*dv1*dv2*wx1+0.04811252243246882*f[9]*dv1_sq*dv2); 
  out[7] += volFact*(2.0*f[5]*wx1_sq*wx2+1.1547005383792517*f[11]*dv1*wx1*wx2+0.16666666666666666*f[5]*dv1_sq*wx2+0.5773502691896258*f[12]*dv2*wx1_sq+0.3333333333333333*f[15]*dv1*dv2*wx1+0.04811252243246882*f[12]*dv1_sq*dv2); 
  out[8] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[3]*dv1*wx2_sq+1.1547005383792517*f[4]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+0.16666666666666666*f[0]*dv2_sq*wx1+0.04811252243246882*f[3]*dv1*dv2_sq); 
  out[9] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[6]*dv1*wx2_sq+1.1547005383792517*f[8]*dv2*wx1*wx2+0.3333333333333333*f[13]*dv1*dv2*wx2+0.16666666666666666*f[1]*dv2_sq*wx1+0.04811252243246882*f[6]*dv1*dv2_sq); 
  out[10] += volFact*(2.0*f[2]*wx1*wx2_sq+0.5773502691896258*f[7]*dv1*wx2_sq+1.1547005383792517*f[9]*dv2*wx1*wx2+0.3333333333333333*f[14]*dv1*dv2*wx2+0.16666666666666666*f[2]*dv2_sq*wx1+0.04811252243246882*f[7]*dv1*dv2_sq); 
  out[11] += volFact*(2.0*f[5]*wx1*wx2_sq+0.5773502691896258*f[11]*dv1*wx2_sq+1.1547005383792517*f[12]*dv2*wx1*wx2+0.3333333333333333*f[15]*dv1*dv2*wx2+0.16666666666666666*f[5]*dv2_sq*wx1+0.04811252243246882*f[11]*dv1*dv2_sq); 
  out[12] += volFact*(1.7320508075688772*f[4]*dv2*wx2_sq+2.0*f[0]*wx2_cu+0.5*f[0]*dv2_sq*wx2+0.08660254037844387*f[4]*dv2_cu); 
  out[13] += volFact*(1.7320508075688772*f[8]*dv2*wx2_sq+2.0*f[1]*wx2_cu+0.5*f[1]*dv2_sq*wx2+0.08660254037844387*f[8]*dv2_cu); 
  out[14] += volFact*(1.7320508075688772*f[9]*dv2*wx2_sq+2.0*f[2]*wx2_cu+0.5*f[2]*dv2_sq*wx2+0.08660254037844387*f[9]*dv2_cu); 
  out[15] += volFact*(1.7320508075688772*f[12]*dv2*wx2_sq+2.0*f[5]*wx2_cu+0.5*f[5]*dv2_sq*wx2+0.08660254037844387*f[12]*dv2_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (1.7320508075688772*hamil[3]*f[8]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (1.7320508075688772*hamil[3]*f[9]+1.7320508075688772*hamil[1]*f[2])*dv10*volFact; 
  out[3] += (1.7320508075688772*hamil[3]*f[12]+1.7320508075688772*hamil[1]*f[5])*dv10*volFact; 
  out[4] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[6]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
  out[6] += (1.7320508075688772*hamil[3]*f[7]+1.7320508075688772*f[2]*hamil[2])*dv11*volFact; 
  out[7] += (1.7320508075688772*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[5])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (hamil[3]*f[10]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (hamil[3]*f[13]+hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[2] += (hamil[3]*f[14]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[3]*f[15]+hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.7320508075688772*hamil[3]*f[4]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[5] += (1.7320508075688772*hamil[3]*f[8]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[6] += (1.7320508075688772*hamil[3]*f[9]+1.7320508075688772*hamil[1]*f[2])*dv10*volFact; 
  out[7] += (1.7320508075688772*hamil[3]*f[12]+1.7320508075688772*hamil[1]*f[5])*dv10*volFact; 
  out[8] += (1.7320508075688772*f[3]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv11*volFact; 
  out[9] += (1.7320508075688772*hamil[3]*f[6]+1.7320508075688772*f[1]*hamil[2])*dv11*volFact; 
  out[10] += (1.7320508075688772*hamil[3]*f[7]+1.7320508075688772*f[2]*hamil[2])*dv11*volFact; 
  out[11] += (1.7320508075688772*hamil[3]*f[11]+1.7320508075688772*hamil[2]*f[5])*dv11*volFact; 
  out[12] += (hamil[3]*f[10]+hamil[2]*f[4]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[13] += (hamil[3]*f[13]+hamil[2]*f[8]+hamil[1]*f[6]+hamil[0]*f[1])*volFact; 
  out[14] += (hamil[3]*f[14]+hamil[2]*f[9]+hamil[1]*f[7]+hamil[0]*f[2])*volFact; 
  out[15] += (hamil[3]*f[15]+hamil[2]*f[12]+hamil[1]*f[11]+hamil[0]*f[5])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (3.4641016151377544*hamil[3]*f[4]+3.4641016151377544*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (3.4641016151377544*f[3]*hamil[3]+3.4641016151377544*f[0]*hamil[2])*dv11*volFact; 
  out[3] += (2.0*hamil[3]*f[10]+2.0*hamil[2]*f[4]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += (0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[9]*hamil[14]+0.8660254037844386*f[8]*hamil[13]+0.8660254037844386*f[5]*hamil[11]+0.8660254037844386*f[4]*hamil[10]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[0]*hamil[3])*dv10*volFact; 
  out[1] += (0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[4]*hamil[13]+0.8660254037844386*f[2]*hamil[11]+0.8660254037844386*f[8]*hamil[10]+0.8660254037844386*f[5]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[1]*hamil[3])*dv10*volFact; 
  out[2] += (0.8660254037844386*f[8]*hamil[15]+0.8660254037844386*f[4]*hamil[14]+0.8660254037844386*f[12]*hamil[13]+0.8660254037844386*f[1]*hamil[11]+0.8660254037844386*f[9]*hamil[10]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[5]*hamil[6]+0.8660254037844386*f[2]*hamil[3])*dv10*volFact; 
  out[3] += (0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[9]*hamil[13]+0.8660254037844386*hamil[10]*f[12]+0.8660254037844386*f[0]*hamil[11]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*hamil[3]*f[5])*dv10*volFact; 
  out[4] += (0.8660254037844386*f[11]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[5]*hamil[12]+0.8660254037844386*f[3]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[0]*hamil[4])*dv11*volFact; 
  out[5] += (0.8660254037844386*f[7]*hamil[15]+0.8660254037844386*f[11]*hamil[14]+0.8660254037844386*f[3]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[5]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[1]*hamil[4])*dv11*volFact; 
  out[6] += (0.8660254037844386*f[6]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[11]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[7]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[5]*hamil[8]+0.8660254037844386*f[2]*hamil[4])*dv11*volFact; 
  out[7] += (0.8660254037844386*f[3]*hamil[15]+0.8660254037844386*f[6]*hamil[14]+0.8660254037844386*f[7]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*hamil[10]*f[11]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*hamil[4]*f[5])*dv11*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*f[14]*hamil[15]+0.5*hamil[14]*f[15]+0.5*f[10]*hamil[13]+0.5*hamil[10]*f[13]+0.5*f[9]*hamil[12]+0.5*hamil[9]*f[12]+0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[4]*hamil[8]+0.5*hamil[4]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[2]*hamil[5]+0.5*hamil[2]*f[5]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*f[13]*hamil[15]+0.5*hamil[13]*f[15]+0.5*f[10]*hamil[14]+0.5*hamil[10]*f[14]+0.5*f[8]*hamil[12]+0.5*hamil[8]*f[12]+0.5*f[6]*hamil[11]+0.5*hamil[6]*f[11]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[1]*hamil[5]+0.5*hamil[1]*f[5]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[10]*hamil[15]+0.5*hamil[10]*f[15]+0.5*f[13]*hamil[14]+0.5*hamil[13]*f[14]+0.5*f[4]*hamil[12]+0.5*hamil[4]*f[12]+0.5*f[3]*hamil[11]+0.5*hamil[3]*f[11]+0.5*f[8]*hamil[9]+0.5*hamil[8]*f[9]+0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[0]*hamil[5]+0.5*hamil[0]*f[5]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]*dxv[3]/4; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (0.8660254037844386*f[12]*hamil[15]+0.8660254037844386*f[9]*hamil[14]+0.8660254037844386*f[8]*hamil[13]+0.8660254037844386*f[5]*hamil[11]+0.8660254037844386*f[4]*hamil[10]+0.8660254037844386*f[2]*hamil[7]+0.8660254037844386*f[1]*hamil[6]+0.8660254037844386*f[0]*hamil[3])*dv10*volFact; 
  out[5] += (0.8660254037844386*f[9]*hamil[15]+0.8660254037844386*f[12]*hamil[14]+0.8660254037844386*f[4]*hamil[13]+0.8660254037844386*f[2]*hamil[11]+0.8660254037844386*f[8]*hamil[10]+0.8660254037844386*f[5]*hamil[7]+0.8660254037844386*f[0]*hamil[6]+0.8660254037844386*f[1]*hamil[3])*dv10*volFact; 
  out[6] += (0.8660254037844386*f[8]*hamil[15]+0.8660254037844386*f[4]*hamil[14]+0.8660254037844386*f[12]*hamil[13]+0.8660254037844386*f[1]*hamil[11]+0.8660254037844386*f[9]*hamil[10]+0.8660254037844386*f[0]*hamil[7]+0.8660254037844386*f[5]*hamil[6]+0.8660254037844386*f[2]*hamil[3])*dv10*volFact; 
  out[7] += (0.8660254037844386*f[4]*hamil[15]+0.8660254037844386*f[8]*hamil[14]+0.8660254037844386*f[9]*hamil[13]+0.8660254037844386*hamil[10]*f[12]+0.8660254037844386*f[0]*hamil[11]+0.8660254037844386*f[1]*hamil[7]+0.8660254037844386*f[2]*hamil[6]+0.8660254037844386*hamil[3]*f[5])*dv10*volFact; 
  out[8] += (0.8660254037844386*f[11]*hamil[15]+0.8660254037844386*f[7]*hamil[14]+0.8660254037844386*f[6]*hamil[13]+0.8660254037844386*f[5]*hamil[12]+0.8660254037844386*f[3]*hamil[10]+0.8660254037844386*f[2]*hamil[9]+0.8660254037844386*f[1]*hamil[8]+0.8660254037844386*f[0]*hamil[4])*dv11*volFact; 
  out[9] += (0.8660254037844386*f[7]*hamil[15]+0.8660254037844386*f[11]*hamil[14]+0.8660254037844386*f[3]*hamil[13]+0.8660254037844386*f[2]*hamil[12]+0.8660254037844386*f[6]*hamil[10]+0.8660254037844386*f[5]*hamil[9]+0.8660254037844386*f[0]*hamil[8]+0.8660254037844386*f[1]*hamil[4])*dv11*volFact; 
  out[10] += (0.8660254037844386*f[6]*hamil[15]+0.8660254037844386*f[3]*hamil[14]+0.8660254037844386*f[11]*hamil[13]+0.8660254037844386*f[1]*hamil[12]+0.8660254037844386*f[7]*hamil[10]+0.8660254037844386*f[0]*hamil[9]+0.8660254037844386*f[5]*hamil[8]+0.8660254037844386*f[2]*hamil[4])*dv11*volFact; 
  out[11] += (0.8660254037844386*f[3]*hamil[15]+0.8660254037844386*f[6]*hamil[14]+0.8660254037844386*f[7]*hamil[13]+0.8660254037844386*f[0]*hamil[12]+0.8660254037844386*hamil[10]*f[11]+0.8660254037844386*f[1]*hamil[9]+0.8660254037844386*f[2]*hamil[8]+0.8660254037844386*hamil[4]*f[5])*dv11*volFact; 
  out[12] += (0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[13] += (0.5*f[14]*hamil[15]+0.5*hamil[14]*f[15]+0.5*f[10]*hamil[13]+0.5*hamil[10]*f[13]+0.5*f[9]*hamil[12]+0.5*hamil[9]*f[12]+0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[4]*hamil[8]+0.5*hamil[4]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[2]*hamil[5]+0.5*hamil[2]*f[5]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[14] += (0.5*f[13]*hamil[15]+0.5*hamil[13]*f[15]+0.5*f[10]*hamil[14]+0.5*hamil[10]*f[14]+0.5*f[8]*hamil[12]+0.5*hamil[8]*f[12]+0.5*f[6]*hamil[11]+0.5*hamil[6]*f[11]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[1]*hamil[5]+0.5*hamil[1]*f[5]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[15] += (0.5*f[10]*hamil[15]+0.5*hamil[10]*f[15]+0.5*f[13]*hamil[14]+0.5*hamil[13]*f[14]+0.5*f[4]*hamil[12]+0.5*hamil[4]*f[12]+0.5*f[3]*hamil[11]+0.5*hamil[3]*f[11]+0.5*f[8]*hamil[9]+0.5*hamil[8]*f[9]+0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[0]*hamil[5]+0.5*hamil[0]*f[5]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[12]*hamil[15]+1.7320508075688772*f[9]*hamil[14]+1.7320508075688772*f[8]*hamil[13]+1.7320508075688772*f[5]*hamil[11]+1.7320508075688772*f[4]*hamil[10]+1.7320508075688772*f[2]*hamil[7]+1.7320508075688772*f[1]*hamil[6]+1.7320508075688772*f[0]*hamil[3])*dv10*volFact; 
  out[2] += (1.7320508075688772*f[11]*hamil[15]+1.7320508075688772*f[7]*hamil[14]+1.7320508075688772*f[6]*hamil[13]+1.7320508075688772*f[5]*hamil[12]+1.7320508075688772*f[3]*hamil[10]+1.7320508075688772*f[2]*hamil[9]+1.7320508075688772*f[1]*hamil[8]+1.7320508075688772*f[0]*hamil[4])*dv11*volFact; 
  out[3] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
