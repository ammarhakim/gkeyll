#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
GKYL_CU_DH void vlasov_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
  out[4] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2); 
  out[5] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2); 
  out[6] += volFact*(2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2); 
  out[7] += volFact*(2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2); 
} 
GKYL_CU_DH void vlasov_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[24]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[25]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[26]*dv2_sq+0.1666666666666667*f[2]*dv2_sq+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx2_sq+1.154700538379252*f[12]*dv2*wx2+2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[5]*dv2_sq+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
} 
GKYL_CU_DH void vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[3]*dv1*wx2+0.5773502691896258*f[4]*dv2*wx1+0.1666666666666667*f[10]*dv1*dv2); 
  out[5] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[6]*dv1*wx2+0.5773502691896258*f[8]*dv2*wx1+0.1666666666666667*f[13]*dv1*dv2); 
  out[6] += volFact*(2.0*f[2]*wx1*wx2+0.5773502691896258*f[7]*dv1*wx2+0.5773502691896258*f[9]*dv2*wx1+0.1666666666666667*f[14]*dv1*dv2); 
  out[7] += volFact*(2.0*f[5]*wx1*wx2+0.5773502691896258*f[11]*dv1*wx2+0.5773502691896258*f[12]*dv2*wx1+0.1666666666666667*f[15]*dv1*dv2); 
  out[8] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+0.149071198499986*f[24]*dv2_sq+0.1666666666666667*f[0]*dv2_sq); 
  out[9] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+0.149071198499986*f[25]*dv2_sq+0.1666666666666667*f[1]*dv2_sq); 
  out[10] += volFact*(2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+0.149071198499986*f[26]*dv2_sq+0.1666666666666667*f[2]*dv2_sq); 
  out[11] += volFact*(2.0*f[5]*wx2_sq+1.154700538379252*f[12]*dv2*wx2+0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[5]*dv2_sq); 
} 
GKYL_CU_DH void vlasov_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[3]*dv1*wx2_sq+1.154700538379252*f[4]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.149071198499986*f[24]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.4472135954999579*f[16]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04303314829119351*f[27]*dv1*dv2_sq+0.04811252243246882*f[3]*dv1*dv2_sq+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[6]*dv1*wx2_sq+1.154700538379252*f[8]*dv2*wx1*wx2+0.3333333333333333*f[13]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.149071198499986*f[25]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.447213595499958*f[17]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04303314829119353*f[29]*dv1*dv2_sq+0.04811252243246882*f[6]*dv1*dv2_sq+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx2_sq+0.5773502691896258*f[7]*dv1*wx2_sq+1.154700538379252*f[9]*dv2*wx1*wx2+0.3333333333333333*f[14]*dv1*dv2*wx2+2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.149071198499986*f[26]*dv2_sq*wx1+0.1666666666666667*f[2]*dv2_sq*wx1+0.447213595499958*f[18]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.04303314829119353*f[30]*dv1*dv2_sq+0.04811252243246882*f[7]*dv1*dv2_sq+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx2_sq+0.5773502691896258*f[11]*dv1*wx2_sq+1.154700538379252*f[12]*dv2*wx1*wx2+0.3333333333333333*f[15]*dv1*dv2*wx2+2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.149071198499986*f[28]*dv2_sq*wx1+0.1666666666666667*f[5]*dv2_sq*wx1+0.4472135954999579*f[20]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.04303314829119351*f[31]*dv1*dv2_sq+0.04811252243246882*f[11]*dv1*dv2_sq+0.08660254037844387*f[11]*dv1*dv1_sq); 
  out[4] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[4]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[3]*dv1*wx1*wx2+0.4472135954999579*f[24]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.149071198499986*f[16]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[4]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.08660254037844387*f[4]*dv2*dv2_sq+0.04303314829119351*f[19]*dv1_sq*dv2+0.04811252243246882*f[4]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[8]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[6]*dv1*wx1*wx2+0.447213595499958*f[25]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.149071198499986*f[17]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[8]*dv2*wx1_sq+0.3333333333333333*f[13]*dv1*dv2*wx1+0.08660254037844387*f[8]*dv2*dv2_sq+0.04303314829119353*f[21]*dv1_sq*dv2+0.04811252243246882*f[8]*dv1_sq*dv2); 
  out[6] += volFact*(2.0*f[2]*wx2*wx2_sq+1.732050807568877*f[9]*dv2*wx2_sq+2.0*f[2]*wx1_sq*wx2+1.154700538379252*f[7]*dv1*wx1*wx2+0.447213595499958*f[26]*dv2_sq*wx2+0.5*f[2]*dv2_sq*wx2+0.149071198499986*f[18]*dv1_sq*wx2+0.1666666666666667*f[2]*dv1_sq*wx2+0.5773502691896258*f[9]*dv2*wx1_sq+0.3333333333333333*f[14]*dv1*dv2*wx1+0.08660254037844387*f[9]*dv2*dv2_sq+0.04303314829119353*f[22]*dv1_sq*dv2+0.04811252243246882*f[9]*dv1_sq*dv2); 
  out[7] += volFact*(2.0*f[5]*wx2*wx2_sq+1.732050807568877*f[12]*dv2*wx2_sq+2.0*f[5]*wx1_sq*wx2+1.154700538379252*f[11]*dv1*wx1*wx2+0.4472135954999579*f[28]*dv2_sq*wx2+0.5*f[5]*dv2_sq*wx2+0.149071198499986*f[20]*dv1_sq*wx2+0.1666666666666667*f[5]*dv1_sq*wx2+0.5773502691896258*f[12]*dv2*wx1_sq+0.3333333333333333*f[15]*dv1*dv2*wx1+0.08660254037844387*f[12]*dv2*dv2_sq+0.04303314829119351*f[23]*dv1_sq*dv2+0.04811252243246882*f[12]*dv1_sq*dv2); 
} 
GKYL_CU_DH void vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.4472135954999579*f[16]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.447213595499958*f[17]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.447213595499958*f[18]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.4472135954999579*f[20]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.08660254037844387*f[11]*dv1*dv1_sq); 
  out[4] += volFact*(2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[3]*dv1*wx1*wx2+0.149071198499986*f[16]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[4]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.04303314829119351*f[19]*dv1_sq*dv2+0.04811252243246882*f[4]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[6]*dv1*wx1*wx2+0.149071198499986*f[17]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[8]*dv2*wx1_sq+0.3333333333333333*f[13]*dv1*dv2*wx1+0.04303314829119353*f[21]*dv1_sq*dv2+0.04811252243246882*f[8]*dv1_sq*dv2); 
  out[6] += volFact*(2.0*f[2]*wx1_sq*wx2+1.154700538379252*f[7]*dv1*wx1*wx2+0.149071198499986*f[18]*dv1_sq*wx2+0.1666666666666667*f[2]*dv1_sq*wx2+0.5773502691896258*f[9]*dv2*wx1_sq+0.3333333333333333*f[14]*dv1*dv2*wx1+0.04303314829119353*f[22]*dv1_sq*dv2+0.04811252243246882*f[9]*dv1_sq*dv2); 
  out[7] += volFact*(2.0*f[5]*wx1_sq*wx2+1.154700538379252*f[11]*dv1*wx1*wx2+0.149071198499986*f[20]*dv1_sq*wx2+0.1666666666666667*f[5]*dv1_sq*wx2+0.5773502691896258*f[12]*dv2*wx1_sq+0.3333333333333333*f[15]*dv1*dv2*wx1+0.04303314829119351*f[23]*dv1_sq*dv2+0.04811252243246882*f[12]*dv1_sq*dv2); 
  out[8] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[3]*dv1*wx2_sq+1.154700538379252*f[4]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+0.149071198499986*f[24]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.04303314829119351*f[27]*dv1*dv2_sq+0.04811252243246882*f[3]*dv1*dv2_sq); 
  out[9] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[6]*dv1*wx2_sq+1.154700538379252*f[8]*dv2*wx1*wx2+0.3333333333333333*f[13]*dv1*dv2*wx2+0.149071198499986*f[25]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.04303314829119353*f[29]*dv1*dv2_sq+0.04811252243246882*f[6]*dv1*dv2_sq); 
  out[10] += volFact*(2.0*f[2]*wx1*wx2_sq+0.5773502691896258*f[7]*dv1*wx2_sq+1.154700538379252*f[9]*dv2*wx1*wx2+0.3333333333333333*f[14]*dv1*dv2*wx2+0.149071198499986*f[26]*dv2_sq*wx1+0.1666666666666667*f[2]*dv2_sq*wx1+0.04303314829119353*f[30]*dv1*dv2_sq+0.04811252243246882*f[7]*dv1*dv2_sq); 
  out[11] += volFact*(2.0*f[5]*wx1*wx2_sq+0.5773502691896258*f[11]*dv1*wx2_sq+1.154700538379252*f[12]*dv2*wx1*wx2+0.3333333333333333*f[15]*dv1*dv2*wx2+0.149071198499986*f[28]*dv2_sq*wx1+0.1666666666666667*f[5]*dv2_sq*wx1+0.04303314829119351*f[31]*dv1*dv2_sq+0.04811252243246882*f[11]*dv1*dv2_sq); 
  out[12] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[4]*dv2*wx2_sq+0.4472135954999579*f[24]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.08660254037844387*f[4]*dv2*dv2_sq); 
  out[13] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[8]*dv2*wx2_sq+0.447213595499958*f[25]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.08660254037844387*f[8]*dv2*dv2_sq); 
  out[14] += volFact*(2.0*f[2]*wx2*wx2_sq+1.732050807568877*f[9]*dv2*wx2_sq+0.447213595499958*f[26]*dv2_sq*wx2+0.5*f[2]*dv2_sq*wx2+0.08660254037844387*f[9]*dv2*dv2_sq); 
  out[15] += volFact*(2.0*f[5]*wx2*wx2_sq+1.732050807568877*f[12]*dv2*wx2_sq+0.4472135954999579*f[28]*dv2_sq*wx2+0.5*f[5]*dv2_sq*wx2+0.08660254037844387*f[12]*dv2*dv2_sq); 
} 
GKYL_CU_DH void vlasov_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double tempM0[4], tempM1i[8]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 
  tempM0[3] = 2.0*f[5]*volFact; 

  tempM1i[0] = tempM0[0]*wx1+0.5773502691896258*f[3]*dv1*volFact; 
  tempM1i[1] = tempM0[1]*wx1+0.5773502691896258*f[6]*dv1*volFact; 
  tempM1i[2] = tempM0[2]*wx1+0.5773502691896258*f[7]*dv1*volFact; 
  tempM1i[3] = tempM0[3]*wx1+0.5773502691896258*f[11]*dv1*volFact; 
  tempM1i[4] = tempM0[0]*wx2+0.5773502691896258*f[4]*dv2*volFact; 
  tempM1i[5] = tempM0[1]*wx2+0.5773502691896258*f[8]*dv2*volFact; 
  tempM1i[6] = tempM0[2]*wx2+0.5773502691896258*f[9]*dv2*volFact; 
  tempM1i[7] = tempM0[3]*wx2+0.5773502691896258*f[12]*dv2*volFact; 

  out[0] += tempM0[0]; 
  out[1] += tempM0[1]; 
  out[2] += tempM0[2]; 
  out[3] += tempM0[3]; 
  out[4] += tempM1i[0]; 
  out[5] += tempM1i[1]; 
  out[6] += tempM1i[2]; 
  out[7] += tempM1i[3]; 
  out[8] += tempM1i[4]; 
  out[9] += tempM1i[5]; 
  out[10] += tempM1i[6]; 
  out[11] += tempM1i[7]; 
  out[12] += tempM0[0]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[4]*wx2+2.0*tempM1i[0]*wx1+(0.149071198499986*f[24]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[13] += tempM0[1]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[5]*wx2+2.0*tempM1i[1]*wx1+(0.149071198499986*f[25]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[14] += tempM0[2]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[6]*wx2+2.0*tempM1i[2]*wx1+(0.149071198499986*f[26]*dv2_sq+0.1666666666666667*f[2]*dv2_sq+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[15] += tempM0[3]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[7]*wx2+2.0*tempM1i[3]*wx1+(0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[5]*dv2_sq+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
} 
