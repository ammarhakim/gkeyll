#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += 1.4142135623730951*f[8]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.4142135623730951*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.10540925533894596*f[5]*dv1_sq+0.11785113019775789*f[0]*dv1_sq); 
  out[1] += volFact*(1.4142135623730951*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.11785113019775789*f[1]*dv1_sq); 
  out[2] += volFact*(1.4142135623730951*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.11785113019775789*f[4]*dv1_sq); 
  out[3] += volFact*(1.4142135623730951*f[8]*wx1_sq+0.8164965809277258*f[10]*dv1*wx1+0.11785113019775789*f[8]*dv1_sq); 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double wx1 = w[1], dv1 = dxv[1]; 
  double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.224744871391589*f[2]*dv1*wx1_sq+1.4142135623730951*f[0]*wx1_cu+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.026726124191242432*f[9]*dv1_cu+0.06123724356957942*f[2]*dv1_cu); 
  out[1] += volFact*(1.224744871391589*f[3]*dv1*wx1_sq+1.4142135623730951*f[1]*wx1_cu+0.31622776601683794*f[7]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.026726124191242435*f[11]*dv1_cu+0.06123724356957942*f[3]*dv1_cu); 
  out[2] += volFact*(1.224744871391589*f[6]*dv1*wx1_sq+1.4142135623730951*f[4]*wx1_cu+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1_cu); 
  out[3] += volFact*(1.2247448713915887*f[10]*dv1*wx1_sq+1.4142135623730951*f[8]*wx1_cu+0.3535533905932737*f[8]*dv1_sq*wx1+0.06123724356957942*f[10]*dv1_cu); 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M1i_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += (5.916079783099617*hamil[3]*f[5]+2.6457513110645907*f[0]*hamil[3]+3.872983346207417*f[2]*hamil[2]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[1] += (5.916079783099617*hamil[3]*f[7]+2.6457513110645907*f[1]*hamil[3]+3.872983346207417*hamil[2]*f[3]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[2] += (3.872983346207417*hamil[2]*f[6]+2.6457513110645907*hamil[3]*f[4]+1.7320508075688772*hamil[1]*f[4])*dv10*volFact; 
  out[3] += (3.8729833462074166*hamil[2]*f[10]+2.6457513110645907*hamil[3]*f[8]+1.7320508075688772*hamil[1]*f[8])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_M2_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (hamil[3]*f[9]+hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0*hamil[3]*f[11]+1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
  out[3] += (1.0*hamil[1]*f[10]+hamil[0]*f[8])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_five_moments_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += 1.4142135623730951*f[8]*volFact; 
  out[4] += (5.916079783099617*hamil[3]*f[5]+2.6457513110645907*f[0]*hamil[3]+3.872983346207417*f[2]*hamil[2]+1.7320508075688772*f[0]*hamil[1])*dv10*volFact; 
  out[5] += (5.916079783099617*hamil[3]*f[7]+2.6457513110645907*f[1]*hamil[3]+3.872983346207417*hamil[2]*f[3]+1.7320508075688772*f[1]*hamil[1])*dv10*volFact; 
  out[6] += (3.872983346207417*hamil[2]*f[6]+2.6457513110645907*hamil[3]*f[4]+1.7320508075688772*hamil[1]*f[4])*dv10*volFact; 
  out[7] += (3.8729833462074166*hamil[2]*f[10]+2.6457513110645907*hamil[3]*f[8]+1.7320508075688772*hamil[1]*f[8])*dv10*volFact; 
  out[8] += (hamil[3]*f[9]+hamil[2]*f[5]+hamil[1]*f[2]+f[0]*hamil[0])*volFact; 
  out[9] += (1.0*hamil[3]*f[11]+1.0000000000000002*hamil[2]*f[7]+hamil[1]*f[3]+hamil[0]*f[1])*volFact; 
  out[10] += (1.0000000000000002*hamil[1]*f[6]+hamil[0]*f[4])*volFact; 
  out[11] += (1.0*hamil[1]*f[10]+hamil[0]*f[8])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_int_five_moments_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (8.366600265340756*hamil[3]*f[5]+3.7416573867739418*f[0]*hamil[3]+5.477225575051662*f[2]*hamil[2]+2.4494897427831783*f[0]*hamil[1])*dv10*volFact; 
  out[2] += (1.4142135623730951*hamil[3]*f[9]+1.4142135623730951*hamil[2]*f[5]+1.4142135623730951*hamil[1]*f[2]+1.4142135623730951*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M1i_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += (4.183300132670377*f[7]*hamil[11]+1.8708286933869707*f[1]*hamil[11]+1.2247448713915887*f[8]*hamil[10]+4.183300132670378*f[5]*hamil[9]+1.8708286933869707*f[0]*hamil[9]+2.7386127875258306*f[3]*hamil[7]+1.224744871391589*f[4]*hamil[6]+2.7386127875258306*f[2]*hamil[5]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[1] += (4.183300132670378*f[5]*hamil[11]+1.6733200530681511*f[4]*hamil[11]+1.8708286933869707*f[0]*hamil[11]+1.075705748400954*f[4]*hamil[10]+4.183300132670379*f[7]*hamil[9]+1.8708286933869707*f[1]*hamil[9]+1.0757057484009542*hamil[6]*f[8]+2.4494897427831783*f[6]*hamil[7]+2.7386127875258306*f[2]*hamil[7]+1.0954451150103324*f[1]*hamil[6]+2.7386127875258306*f[3]*hamil[5]+1.0954451150103324*hamil[3]*f[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[2] += (1.6431676725154978*f[8]*hamil[11]+3.7416573867739413*f[7]*hamil[11]+1.6733200530681511*f[1]*hamil[11]+0.7302967433402215*f[8]*hamil[10]+1.075705748400954*f[1]*hamil[10]+2.405351177211819*hamil[7]*f[10]+1.8708286933869707*f[4]*hamil[9]+1.075705748400954*hamil[3]*f[8]+2.4494897427831783*f[3]*hamil[7]+0.7824607964359517*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[6]+2.7386127875258306*hamil[5]*f[6]+1.224744871391589*hamil[2]*f[4]+1.0954451150103324*f[1]*hamil[3])*dv10*volFact; 
  out[3] += (1.6431676725154978*f[4]*hamil[11]+0.7302967433402215*f[4]*hamil[10]+1.2247448713915887*f[0]*hamil[10]+2.7386127875258297*hamil[5]*f[10]+1.8708286933869707*f[8]*hamil[9]+0.7302967433402215*hamil[6]*f[8]+1.224744871391589*hamil[2]*f[8]+2.405351177211819*f[6]*hamil[7]+1.0757057484009542*f[1]*hamil[6]+1.075705748400954*hamil[3]*f[4])*dv10*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_M2_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  out[0] += (0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865474*f[9]*hamil[11]+0.7071067811865474*hamil[9]*f[11]+0.6210590034081187*f[6]*hamil[10]+0.6210590034081187*hamil[6]*f[10]+0.6210590034081186*f[4]*hamil[8]+0.6210590034081186*hamil[4]*f[8]+0.7071067811865475*f[5]*hamil[7]+0.7071067811865475*hamil[5]*f[7]+0.632455532033676*f[3]*hamil[6]+0.632455532033676*hamil[3]*f[6]+0.6324555320336759*f[1]*hamil[4]+0.6324555320336759*hamil[1]*f[4]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[2] += (0.6324555320336759*f[11]*hamil[11]+0.421637021355784*f[10]*hamil[10]+0.6210590034081187*f[3]*hamil[10]+0.6210590034081187*hamil[3]*f[10]+0.421637021355784*f[8]*hamil[8]+0.6210590034081186*f[1]*hamil[8]+0.6210590034081186*hamil[1]*f[8]+0.6324555320336759*f[7]*hamil[7]+0.45175395145262565*f[6]*hamil[6]+0.7071067811865475*f[2]*hamil[6]+0.7071067811865475*hamil[2]*f[6]+0.45175395145262565*f[4]*hamil[4]+0.7071067811865475*f[0]*hamil[4]+0.7071067811865475*hamil[0]*f[4]+0.6324555320336759*f[3]*hamil[3]+0.6324555320336759*f[1]*hamil[1])*volFact; 
  out[3] += (0.4216370213557839*f[6]*hamil[10]+0.7071067811865474*f[2]*hamil[10]+0.4216370213557839*hamil[6]*f[10]+0.7071067811865474*hamil[2]*f[10]+0.421637021355784*f[4]*hamil[8]+0.7071067811865475*f[0]*hamil[8]+0.421637021355784*hamil[4]*f[8]+0.7071067811865475*hamil[0]*f[8]+0.6210590034081187*f[3]*hamil[6]+0.6210590034081187*hamil[3]*f[6]+0.6210590034081186*f[1]*hamil[4]+0.6210590034081186*hamil[1]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_five_moments_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[1]/2; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[4]*volFact; 
  out[3] += 1.4142135623730951*f[8]*volFact; 
  out[4] += (4.183300132670377*f[7]*hamil[11]+1.8708286933869707*f[1]*hamil[11]+1.2247448713915887*f[8]*hamil[10]+4.183300132670378*f[5]*hamil[9]+1.8708286933869707*f[0]*hamil[9]+2.7386127875258306*f[3]*hamil[7]+1.224744871391589*f[4]*hamil[6]+2.7386127875258306*f[2]*hamil[5]+1.224744871391589*f[1]*hamil[3]+1.224744871391589*f[0]*hamil[2])*dv10*volFact; 
  out[5] += (4.183300132670378*f[5]*hamil[11]+1.6733200530681511*f[4]*hamil[11]+1.8708286933869707*f[0]*hamil[11]+1.075705748400954*f[4]*hamil[10]+4.183300132670379*f[7]*hamil[9]+1.8708286933869707*f[1]*hamil[9]+1.0757057484009542*hamil[6]*f[8]+2.4494897427831783*f[6]*hamil[7]+2.7386127875258306*f[2]*hamil[7]+1.0954451150103324*f[1]*hamil[6]+2.7386127875258306*f[3]*hamil[5]+1.0954451150103324*hamil[3]*f[4]+1.224744871391589*f[0]*hamil[3]+1.224744871391589*f[1]*hamil[2])*dv10*volFact; 
  out[6] += (1.6431676725154978*f[8]*hamil[11]+3.7416573867739413*f[7]*hamil[11]+1.6733200530681511*f[1]*hamil[11]+0.7302967433402215*f[8]*hamil[10]+1.075705748400954*f[1]*hamil[10]+2.405351177211819*hamil[7]*f[10]+1.8708286933869707*f[4]*hamil[9]+1.075705748400954*hamil[3]*f[8]+2.4494897427831783*f[3]*hamil[7]+0.7824607964359517*f[4]*hamil[6]+1.224744871391589*f[0]*hamil[6]+2.7386127875258306*hamil[5]*f[6]+1.224744871391589*hamil[2]*f[4]+1.0954451150103324*f[1]*hamil[3])*dv10*volFact; 
  out[7] += (1.6431676725154978*f[4]*hamil[11]+0.7302967433402215*f[4]*hamil[10]+1.2247448713915887*f[0]*hamil[10]+2.7386127875258297*hamil[5]*f[10]+1.8708286933869707*f[8]*hamil[9]+0.7302967433402215*hamil[6]*f[8]+1.224744871391589*hamil[2]*f[8]+2.405351177211819*f[6]*hamil[7]+1.0757057484009542*f[1]*hamil[6]+1.075705748400954*hamil[3]*f[4])*dv10*volFact; 
  out[8] += (0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[9] += (0.7071067811865474*f[9]*hamil[11]+0.7071067811865474*hamil[9]*f[11]+0.6210590034081187*f[6]*hamil[10]+0.6210590034081187*hamil[6]*f[10]+0.6210590034081186*f[4]*hamil[8]+0.6210590034081186*hamil[4]*f[8]+0.7071067811865475*f[5]*hamil[7]+0.7071067811865475*hamil[5]*f[7]+0.632455532033676*f[3]*hamil[6]+0.632455532033676*hamil[3]*f[6]+0.6324555320336759*f[1]*hamil[4]+0.6324555320336759*hamil[1]*f[4]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[10] += (0.6324555320336759*f[11]*hamil[11]+0.421637021355784*f[10]*hamil[10]+0.6210590034081187*f[3]*hamil[10]+0.6210590034081187*hamil[3]*f[10]+0.421637021355784*f[8]*hamil[8]+0.6210590034081186*f[1]*hamil[8]+0.6210590034081186*hamil[1]*f[8]+0.6324555320336759*f[7]*hamil[7]+0.45175395145262565*f[6]*hamil[6]+0.7071067811865475*f[2]*hamil[6]+0.7071067811865475*hamil[2]*f[6]+0.45175395145262565*f[4]*hamil[4]+0.7071067811865475*f[0]*hamil[4]+0.7071067811865475*hamil[0]*f[4]+0.6324555320336759*f[3]*hamil[3]+0.6324555320336759*f[1]*hamil[1])*volFact; 
  out[11] += (0.4216370213557839*f[6]*hamil[10]+0.7071067811865474*f[2]*hamil[10]+0.4216370213557839*hamil[6]*f[10]+0.7071067811865474*hamil[2]*f[10]+0.421637021355784*f[4]*hamil[8]+0.7071067811865475*f[0]*hamil[8]+0.421637021355784*hamil[4]*f[8]+0.7071067811865475*hamil[0]*f[8]+0.6210590034081187*f[3]*hamil[6]+0.6210590034081187*hamil[3]*f[6]+0.6210590034081186*f[1]*hamil[4]+0.6210590034081186*hamil[1]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_gen_int_five_moments_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*0.25; 
  double dv10 = 2.0/dxv[1]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (5.916079783099615*f[7]*hamil[11]+2.6457513110645907*f[1]*hamil[11]+1.7320508075688772*f[8]*hamil[10]+5.916079783099617*f[5]*hamil[9]+2.6457513110645907*f[0]*hamil[9]+3.872983346207417*f[3]*hamil[7]+1.7320508075688774*f[4]*hamil[6]+3.872983346207417*f[2]*hamil[5]+1.7320508075688772*f[1]*hamil[3]+1.7320508075688772*f[0]*hamil[2])*dv10*volFact; 
  out[2] += (f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
