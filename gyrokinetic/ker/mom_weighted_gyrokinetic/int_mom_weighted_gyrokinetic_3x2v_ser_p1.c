#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 

  out[0] += (2.0*wgt[7]*f[16]+2.0*wgt[6]*f[8]+2.0*wgt[5]*f[7]+2.0*wgt[4]*f[6]+2.0*f[3]*wgt[3]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 

  out[0] += (1.4142135623730951*vmap[1]*wgt[7]*f[26]+1.4142135623730951*vmap[1]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[5]*f[18]+1.4142135623730951*vmap[1]*wgt[4]*f[17]+1.4142135623730951*vmap[0]*wgt[7]*f[16]+1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.4142135623730951*vmap[1]*wgt[1]*f[9]+1.4142135623730951*vmap[0]*wgt[6]*f[8]+1.4142135623730951*vmap[0]*wgt[5]*f[7]+1.4142135623730951*vmap[0]*wgt[4]*f[6]+1.4142135623730951*wgt[0]*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[3]*wgt[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[7]*f[43]+0.8944271909999159*vmap1R2*wgt[6]*f[39]+0.8944271909999159*vmap1R2*wgt[5]*f[38]+0.8944271909999159*vmap1R2*wgt[4]*f[37]+0.8944271909999161*vmap1R2*wgt[3]*f[35]+0.8944271909999161*vmap1R2*wgt[2]*f[34]+0.8944271909999161*vmap1R2*wgt[1]*f[33]+0.8944271909999159*wgt[0]*vmap1R2*f[32]+2.0*vmap[0]*vmap[1]*wgt[7]*f[26]+2.0*vmap[0]*vmap[1]*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[5]*f[18]+2.0*vmap[0]*vmap[1]*wgt[4]*f[17]+vmap1R2*wgt[7]*f[16]+vmap0R2*wgt[7]*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+2.0*vmap[0]*vmap[1]*wgt[1]*f[9]+vmap1R2*wgt[6]*f[8]+vmap0R2*wgt[6]*f[8]+vmap1R2*wgt[5]*f[7]+vmap0R2*wgt[5]*f[7]+vmap1R2*wgt[4]*f[6]+vmap0R2*wgt[4]*f[6]+2.0*vmap[0]*wgt[0]*vmap[1]*f[4]+vmap1R2*f[3]*wgt[3]+vmap0R2*f[3]*wgt[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[7]*f[43]+0.8944271909999159*vmap1R2*wgt[6]*f[39]+0.8944271909999159*vmap1R2*wgt[5]*f[38]+0.8944271909999159*vmap1R2*wgt[4]*f[37]+0.8944271909999161*vmap1R2*wgt[3]*f[35]+0.8944271909999161*vmap1R2*wgt[2]*f[34]+0.8944271909999161*vmap1R2*wgt[1]*f[33]+0.8944271909999159*wgt[0]*vmap1R2*f[32]+2.0*vmap[0]*vmap[1]*wgt[7]*f[26]+2.0*vmap[0]*vmap[1]*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[5]*f[18]+2.0*vmap[0]*vmap[1]*wgt[4]*f[17]+vmap1R2*wgt[7]*f[16]+vmap0R2*wgt[7]*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+2.0*vmap[0]*vmap[1]*wgt[1]*f[9]+vmap1R2*wgt[6]*f[8]+vmap0R2*wgt[6]*f[8]+0.3535533905932737*bmag[0]*tmp[7]*wgt[7]+0.3535533905932737*tmp[0]*bmag[7]*wgt[7]+0.3535533905932737*bmag[1]*tmp[6]*wgt[7]+0.3535533905932737*tmp[1]*bmag[6]*wgt[7]+0.3535533905932737*bmag[2]*tmp[5]*wgt[7]+0.3535533905932737*tmp[2]*bmag[5]*wgt[7]+0.3535533905932737*bmag[3]*tmp[4]*wgt[7]+0.3535533905932737*tmp[3]*bmag[4]*wgt[7]+0.3535533905932737*wgt[0]*bmag[7]*tmp[7]+0.3535533905932737*bmag[1]*wgt[6]*tmp[7]+0.3535533905932737*wgt[1]*bmag[6]*tmp[7]+0.3535533905932737*bmag[2]*wgt[5]*tmp[7]+0.3535533905932737*wgt[2]*bmag[5]*tmp[7]+0.3535533905932737*bmag[3]*wgt[4]*tmp[7]+0.3535533905932737*wgt[3]*bmag[4]*tmp[7]+vmap1R2*wgt[5]*f[7]+vmap0R2*wgt[5]*f[7]+0.3535533905932737*tmp[1]*wgt[6]*bmag[7]+0.3535533905932737*wgt[1]*tmp[6]*bmag[7]+0.3535533905932737*tmp[2]*wgt[5]*bmag[7]+0.3535533905932737*wgt[2]*tmp[5]*bmag[7]+0.3535533905932737*tmp[3]*wgt[4]*bmag[7]+0.3535533905932737*wgt[3]*tmp[4]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]*wgt[6]+0.3535533905932737*tmp[0]*bmag[6]*wgt[6]+0.3535533905932737*bmag[4]*tmp[5]*wgt[6]+0.3535533905932737*tmp[4]*bmag[5]*wgt[6]+0.3535533905932737*bmag[2]*tmp[3]*wgt[6]+0.3535533905932737*tmp[2]*bmag[3]*wgt[6]+0.3535533905932737*wgt[0]*bmag[6]*tmp[6]+0.3535533905932737*bmag[4]*wgt[5]*tmp[6]+0.3535533905932737*wgt[4]*bmag[5]*tmp[6]+0.3535533905932737*bmag[2]*wgt[3]*tmp[6]+0.3535533905932737*wgt[2]*bmag[3]*tmp[6]+vmap1R2*wgt[4]*f[6]+vmap0R2*wgt[4]*f[6]+0.3535533905932737*tmp[4]*wgt[5]*bmag[6]+0.3535533905932737*wgt[4]*tmp[5]*bmag[6]+0.3535533905932737*tmp[2]*wgt[3]*bmag[6]+0.3535533905932737*wgt[2]*tmp[3]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]*wgt[5]+0.3535533905932737*tmp[0]*bmag[5]*wgt[5]+0.3535533905932737*bmag[1]*tmp[3]*wgt[5]+0.3535533905932737*tmp[1]*bmag[3]*wgt[5]+0.3535533905932737*wgt[0]*bmag[5]*tmp[5]+0.3535533905932737*bmag[1]*wgt[3]*tmp[5]+0.3535533905932737*wgt[1]*bmag[3]*tmp[5]+0.3535533905932737*tmp[1]*wgt[3]*bmag[5]+0.3535533905932737*wgt[1]*tmp[3]*bmag[5]+0.3535533905932737*bmag[0]*tmp[4]*wgt[4]+0.3535533905932737*tmp[0]*bmag[4]*wgt[4]+0.3535533905932737*bmag[1]*tmp[2]*wgt[4]+0.3535533905932737*tmp[1]*bmag[2]*wgt[4]+0.3535533905932737*wgt[0]*bmag[4]*tmp[4]+0.3535533905932737*bmag[1]*wgt[2]*tmp[4]+0.3535533905932737*wgt[1]*bmag[2]*tmp[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[4]+0.3535533905932737*tmp[1]*wgt[2]*bmag[4]+0.3535533905932737*wgt[1]*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[3]*wgt[3]+vmap1R2*f[3]*wgt[3]+vmap0R2*f[3]*wgt[3]+0.3535533905932737*tmp[0]*bmag[3]*wgt[3]+0.3535533905932737*wgt[0]*bmag[3]*tmp[3]+0.3535533905932737*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.3535533905932737*tmp[0]*bmag[2]*wgt[2]+0.3535533905932737*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.3535533905932737*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.3535533905932737*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.3535533905932737*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.3535533905932737*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[7]*f[16]+2.0*wgt[6]*f[8]+2.0*wgt[5]*f[7]+2.0*wgt[4]*f[6]+2.0*f[3]*wgt[3]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[7]*f[26]+1.4142135623730951*vmap[1]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[5]*f[18]+1.4142135623730951*vmap[1]*wgt[4]*f[17]+1.4142135623730951*vmap[0]*wgt[7]*f[16]+1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.4142135623730951*vmap[1]*wgt[1]*f[9]+1.4142135623730951*vmap[0]*wgt[6]*f[8]+1.4142135623730951*vmap[0]*wgt[5]*f[7]+1.4142135623730951*vmap[0]*wgt[4]*f[6]+1.4142135623730951*wgt[0]*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[3]*wgt[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[7]*f[43]+0.8944271909999159*vmap1R2*wgt[6]*f[39]+0.8944271909999159*vmap1R2*wgt[5]*f[38]+0.8944271909999159*vmap1R2*wgt[4]*f[37]+0.8944271909999161*vmap1R2*wgt[3]*f[35]+0.8944271909999161*vmap1R2*wgt[2]*f[34]+0.8944271909999161*vmap1R2*wgt[1]*f[33]+0.8944271909999159*wgt[0]*vmap1R2*f[32]+2.0*vmap[0]*vmap[1]*wgt[7]*f[26]+2.0*vmap[0]*vmap[1]*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[5]*f[18]+2.0*vmap[0]*vmap[1]*wgt[4]*f[17]+vmap1R2*wgt[7]*f[16]+vmap0R2*wgt[7]*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+2.0*vmap[0]*vmap[1]*wgt[1]*f[9]+vmap1R2*wgt[6]*f[8]+vmap0R2*wgt[6]*f[8]+0.3535533905932737*bmag[0]*tmp[7]*wgt[7]+0.3535533905932737*tmp[0]*bmag[7]*wgt[7]+0.3535533905932737*bmag[1]*tmp[6]*wgt[7]+0.3535533905932737*tmp[1]*bmag[6]*wgt[7]+0.3535533905932737*bmag[2]*tmp[5]*wgt[7]+0.3535533905932737*tmp[2]*bmag[5]*wgt[7]+0.3535533905932737*bmag[3]*tmp[4]*wgt[7]+0.3535533905932737*tmp[3]*bmag[4]*wgt[7]+0.3535533905932737*wgt[0]*bmag[7]*tmp[7]+0.3535533905932737*bmag[1]*wgt[6]*tmp[7]+0.3535533905932737*wgt[1]*bmag[6]*tmp[7]+0.3535533905932737*bmag[2]*wgt[5]*tmp[7]+0.3535533905932737*wgt[2]*bmag[5]*tmp[7]+0.3535533905932737*bmag[3]*wgt[4]*tmp[7]+0.3535533905932737*wgt[3]*bmag[4]*tmp[7]+vmap1R2*wgt[5]*f[7]+vmap0R2*wgt[5]*f[7]+0.3535533905932737*tmp[1]*wgt[6]*bmag[7]+0.3535533905932737*wgt[1]*tmp[6]*bmag[7]+0.3535533905932737*tmp[2]*wgt[5]*bmag[7]+0.3535533905932737*wgt[2]*tmp[5]*bmag[7]+0.3535533905932737*tmp[3]*wgt[4]*bmag[7]+0.3535533905932737*wgt[3]*tmp[4]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]*wgt[6]+0.3535533905932737*tmp[0]*bmag[6]*wgt[6]+0.3535533905932737*bmag[4]*tmp[5]*wgt[6]+0.3535533905932737*tmp[4]*bmag[5]*wgt[6]+0.3535533905932737*bmag[2]*tmp[3]*wgt[6]+0.3535533905932737*tmp[2]*bmag[3]*wgt[6]+0.3535533905932737*wgt[0]*bmag[6]*tmp[6]+0.3535533905932737*bmag[4]*wgt[5]*tmp[6]+0.3535533905932737*wgt[4]*bmag[5]*tmp[6]+0.3535533905932737*bmag[2]*wgt[3]*tmp[6]+0.3535533905932737*wgt[2]*bmag[3]*tmp[6]+vmap1R2*wgt[4]*f[6]+vmap0R2*wgt[4]*f[6]+0.3535533905932737*tmp[4]*wgt[5]*bmag[6]+0.3535533905932737*wgt[4]*tmp[5]*bmag[6]+0.3535533905932737*tmp[2]*wgt[3]*bmag[6]+0.3535533905932737*wgt[2]*tmp[3]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]*wgt[5]+0.3535533905932737*tmp[0]*bmag[5]*wgt[5]+0.3535533905932737*bmag[1]*tmp[3]*wgt[5]+0.3535533905932737*tmp[1]*bmag[3]*wgt[5]+0.3535533905932737*wgt[0]*bmag[5]*tmp[5]+0.3535533905932737*bmag[1]*wgt[3]*tmp[5]+0.3535533905932737*wgt[1]*bmag[3]*tmp[5]+0.3535533905932737*tmp[1]*wgt[3]*bmag[5]+0.3535533905932737*wgt[1]*tmp[3]*bmag[5]+0.3535533905932737*bmag[0]*tmp[4]*wgt[4]+0.3535533905932737*tmp[0]*bmag[4]*wgt[4]+0.3535533905932737*bmag[1]*tmp[2]*wgt[4]+0.3535533905932737*tmp[1]*bmag[2]*wgt[4]+0.3535533905932737*wgt[0]*bmag[4]*tmp[4]+0.3535533905932737*bmag[1]*wgt[2]*tmp[4]+0.3535533905932737*wgt[1]*bmag[2]*tmp[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[4]+0.3535533905932737*tmp[1]*wgt[2]*bmag[4]+0.3535533905932737*wgt[1]*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[3]*wgt[3]+vmap1R2*f[3]*wgt[3]+vmap0R2*f[3]*wgt[3]+0.3535533905932737*tmp[0]*bmag[3]*wgt[3]+0.3535533905932737*wgt[0]*bmag[3]*tmp[3]+0.3535533905932737*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.3535533905932737*tmp[0]*bmag[2]*wgt[2]+0.3535533905932737*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.3535533905932737*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.3535533905932737*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.3535533905932737*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.3535533905932737*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_NOWEIGHT_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[14]+1.4142135623730951*vmap[2]*f[3]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[20]+1.4142135623730951*vmap[2]*f[6]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[21]+1.4142135623730951*vmap[2]*f[7]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[22]+1.4142135623730951*vmap[2]*f[8]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[27]+1.4142135623730951*vmap[2]*f[16]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*m_*volFact; 
  out[2] += (2.0*phi[7]*f[16]*q_+2.0*phi[6]*f[8]*q_+2.0*phi[5]*f[7]*q_+2.0*phi[4]*f[6]*q_+2.0*f[3]*phi[3]*q_+2.0*f[2]*phi[2]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact+(1.264911064067352*vmap1R2*f[32]+2.8284271247461907*vmap[0]*vmap[1]*f[4]+1.4142135623730951*f[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R2)*m_*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 

  out[0] += (0.3535533905932737*bmag[0]*tmp[7]*wgt[7]+0.3535533905932737*tmp[0]*bmag[7]*wgt[7]+0.3535533905932737*bmag[1]*tmp[6]*wgt[7]+0.3535533905932737*tmp[1]*bmag[6]*wgt[7]+0.3535533905932737*bmag[2]*tmp[5]*wgt[7]+0.3535533905932737*tmp[2]*bmag[5]*wgt[7]+0.3535533905932737*bmag[3]*tmp[4]*wgt[7]+0.3535533905932737*tmp[3]*bmag[4]*wgt[7]+0.3535533905932737*wgt[0]*bmag[7]*tmp[7]+0.3535533905932737*bmag[1]*wgt[6]*tmp[7]+0.3535533905932737*wgt[1]*bmag[6]*tmp[7]+0.3535533905932737*bmag[2]*wgt[5]*tmp[7]+0.3535533905932737*wgt[2]*bmag[5]*tmp[7]+0.3535533905932737*bmag[3]*wgt[4]*tmp[7]+0.3535533905932737*wgt[3]*bmag[4]*tmp[7]+0.3535533905932737*tmp[1]*wgt[6]*bmag[7]+0.3535533905932737*wgt[1]*tmp[6]*bmag[7]+0.3535533905932737*tmp[2]*wgt[5]*bmag[7]+0.3535533905932737*wgt[2]*tmp[5]*bmag[7]+0.3535533905932737*tmp[3]*wgt[4]*bmag[7]+0.3535533905932737*wgt[3]*tmp[4]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]*wgt[6]+0.3535533905932737*tmp[0]*bmag[6]*wgt[6]+0.3535533905932737*bmag[4]*tmp[5]*wgt[6]+0.3535533905932737*tmp[4]*bmag[5]*wgt[6]+0.3535533905932737*bmag[2]*tmp[3]*wgt[6]+0.3535533905932737*tmp[2]*bmag[3]*wgt[6]+0.3535533905932737*wgt[0]*bmag[6]*tmp[6]+0.3535533905932737*bmag[4]*wgt[5]*tmp[6]+0.3535533905932737*wgt[4]*bmag[5]*tmp[6]+0.3535533905932737*bmag[2]*wgt[3]*tmp[6]+0.3535533905932737*wgt[2]*bmag[3]*tmp[6]+0.3535533905932737*tmp[4]*wgt[5]*bmag[6]+0.3535533905932737*wgt[4]*tmp[5]*bmag[6]+0.3535533905932737*tmp[2]*wgt[3]*bmag[6]+0.3535533905932737*wgt[2]*tmp[3]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]*wgt[5]+0.3535533905932737*tmp[0]*bmag[5]*wgt[5]+0.3535533905932737*bmag[1]*tmp[3]*wgt[5]+0.3535533905932737*tmp[1]*bmag[3]*wgt[5]+0.3535533905932737*wgt[0]*bmag[5]*tmp[5]+0.3535533905932737*bmag[1]*wgt[3]*tmp[5]+0.3535533905932737*wgt[1]*bmag[3]*tmp[5]+0.3535533905932737*tmp[1]*wgt[3]*bmag[5]+0.3535533905932737*wgt[1]*tmp[3]*bmag[5]+0.3535533905932737*bmag[0]*tmp[4]*wgt[4]+0.3535533905932737*tmp[0]*bmag[4]*wgt[4]+0.3535533905932737*bmag[1]*tmp[2]*wgt[4]+0.3535533905932737*tmp[1]*bmag[2]*wgt[4]+0.3535533905932737*wgt[0]*bmag[4]*tmp[4]+0.3535533905932737*bmag[1]*wgt[2]*tmp[4]+0.3535533905932737*wgt[1]*bmag[2]*tmp[4]+0.3535533905932737*tmp[1]*wgt[2]*bmag[4]+0.3535533905932737*wgt[1]*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[3]*wgt[3]+0.3535533905932737*tmp[0]*bmag[3]*wgt[3]+0.3535533905932737*wgt[0]*bmag[3]*tmp[3]+0.3535533905932737*bmag[0]*tmp[2]*wgt[2]+0.3535533905932737*tmp[0]*bmag[2]*wgt[2]+0.3535533905932737*wgt[0]*bmag[2]*tmp[2]+0.3535533905932737*bmag[0]*tmp[1]*wgt[1]+0.3535533905932737*tmp[0]*bmag[1]*wgt[1]+0.3535533905932737*wgt[0]*bmag[1]*tmp[1]+0.3535533905932737*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.19634954084936207*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[12])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[13])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[14])/m_+(2.8284271247461907*vmap[2]*f[3])/m_; 
  tmp[4] = (2.8284271247461907*vmap[3]*f[20])/m_+(2.8284271247461907*vmap[2]*f[6])/m_; 
  tmp[5] = (2.8284271247461907*vmap[3]*f[21])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
  tmp[6] = (2.8284271247461907*vmap[3]*f[22])/m_+(2.8284271247461907*vmap[2]*f[8])/m_; 
  tmp[7] = (2.8284271247461907*vmap[3]*f[27])/m_+(2.8284271247461907*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[7]*f[16]+2.0*wgt[6]*f[8]+2.0*wgt[5]*f[7]+2.0*wgt[4]*f[6]+2.0*f[3]*wgt[3]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[7]*f[26]+1.4142135623730951*vmap[1]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[5]*f[18]+1.4142135623730951*vmap[1]*wgt[4]*f[17]+1.4142135623730951*vmap[0]*wgt[7]*f[16]+1.4142135623730951*vmap[1]*wgt[3]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[10]+1.4142135623730951*vmap[1]*wgt[1]*f[9]+1.4142135623730951*vmap[0]*wgt[6]*f[8]+1.4142135623730951*vmap[0]*wgt[5]*f[7]+1.4142135623730951*vmap[0]*wgt[4]*f[6]+1.4142135623730951*wgt[0]*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[3]*wgt[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[7]*f[43]+0.8944271909999159*vmap1R2*wgt[6]*f[39]+0.8944271909999159*vmap1R2*wgt[5]*f[38]+0.8944271909999159*vmap1R2*wgt[4]*f[37]+0.8944271909999161*vmap1R2*wgt[3]*f[35]+0.8944271909999161*vmap1R2*wgt[2]*f[34]+0.8944271909999161*vmap1R2*wgt[1]*f[33]+0.8944271909999159*wgt[0]*vmap1R2*f[32]+2.0*vmap[0]*vmap[1]*wgt[7]*f[26]+2.0*vmap[0]*vmap[1]*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[5]*f[18]+2.0*vmap[0]*vmap[1]*wgt[4]*f[17]+vmap1R2*wgt[7]*f[16]+vmap0R2*wgt[7]*f[16]+2.0*vmap[0]*vmap[1]*wgt[3]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[10]+2.0*vmap[0]*vmap[1]*wgt[1]*f[9]+vmap1R2*wgt[6]*f[8]+vmap0R2*wgt[6]*f[8]+vmap1R2*wgt[5]*f[7]+vmap0R2*wgt[5]*f[7]+vmap1R2*wgt[4]*f[6]+vmap0R2*wgt[4]*f[6]+2.0*vmap[0]*wgt[0]*vmap[1]*f[4]+vmap1R2*f[3]*wgt[3]+vmap0R2*f[3]*wgt[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (0.3535533905932737*bmag[0]*tmp[7]*wgt[7]+0.3535533905932737*tmp[0]*bmag[7]*wgt[7]+0.3535533905932737*bmag[1]*tmp[6]*wgt[7]+0.3535533905932737*tmp[1]*bmag[6]*wgt[7]+0.3535533905932737*bmag[2]*tmp[5]*wgt[7]+0.3535533905932737*tmp[2]*bmag[5]*wgt[7]+0.3535533905932737*bmag[3]*tmp[4]*wgt[7]+0.3535533905932737*tmp[3]*bmag[4]*wgt[7]+0.3535533905932737*wgt[0]*bmag[7]*tmp[7]+0.3535533905932737*bmag[1]*wgt[6]*tmp[7]+0.3535533905932737*wgt[1]*bmag[6]*tmp[7]+0.3535533905932737*bmag[2]*wgt[5]*tmp[7]+0.3535533905932737*wgt[2]*bmag[5]*tmp[7]+0.3535533905932737*bmag[3]*wgt[4]*tmp[7]+0.3535533905932737*wgt[3]*bmag[4]*tmp[7]+0.3535533905932737*tmp[1]*wgt[6]*bmag[7]+0.3535533905932737*wgt[1]*tmp[6]*bmag[7]+0.3535533905932737*tmp[2]*wgt[5]*bmag[7]+0.3535533905932737*wgt[2]*tmp[5]*bmag[7]+0.3535533905932737*tmp[3]*wgt[4]*bmag[7]+0.3535533905932737*wgt[3]*tmp[4]*bmag[7]+0.3535533905932737*bmag[0]*tmp[6]*wgt[6]+0.3535533905932737*tmp[0]*bmag[6]*wgt[6]+0.3535533905932737*bmag[4]*tmp[5]*wgt[6]+0.3535533905932737*tmp[4]*bmag[5]*wgt[6]+0.3535533905932737*bmag[2]*tmp[3]*wgt[6]+0.3535533905932737*tmp[2]*bmag[3]*wgt[6]+0.3535533905932737*wgt[0]*bmag[6]*tmp[6]+0.3535533905932737*bmag[4]*wgt[5]*tmp[6]+0.3535533905932737*wgt[4]*bmag[5]*tmp[6]+0.3535533905932737*bmag[2]*wgt[3]*tmp[6]+0.3535533905932737*wgt[2]*bmag[3]*tmp[6]+0.3535533905932737*tmp[4]*wgt[5]*bmag[6]+0.3535533905932737*wgt[4]*tmp[5]*bmag[6]+0.3535533905932737*tmp[2]*wgt[3]*bmag[6]+0.3535533905932737*wgt[2]*tmp[3]*bmag[6]+0.3535533905932737*bmag[0]*tmp[5]*wgt[5]+0.3535533905932737*tmp[0]*bmag[5]*wgt[5]+0.3535533905932737*bmag[1]*tmp[3]*wgt[5]+0.3535533905932737*tmp[1]*bmag[3]*wgt[5]+0.3535533905932737*wgt[0]*bmag[5]*tmp[5]+0.3535533905932737*bmag[1]*wgt[3]*tmp[5]+0.3535533905932737*wgt[1]*bmag[3]*tmp[5]+0.3535533905932737*tmp[1]*wgt[3]*bmag[5]+0.3535533905932737*wgt[1]*tmp[3]*bmag[5]+0.3535533905932737*bmag[0]*tmp[4]*wgt[4]+0.3535533905932737*tmp[0]*bmag[4]*wgt[4]+0.3535533905932737*bmag[1]*tmp[2]*wgt[4]+0.3535533905932737*tmp[1]*bmag[2]*wgt[4]+0.3535533905932737*wgt[0]*bmag[4]*tmp[4]+0.3535533905932737*bmag[1]*wgt[2]*tmp[4]+0.3535533905932737*wgt[1]*bmag[2]*tmp[4]+0.3535533905932737*tmp[1]*wgt[2]*bmag[4]+0.3535533905932737*wgt[1]*tmp[2]*bmag[4]+0.3535533905932737*bmag[0]*tmp[3]*wgt[3]+0.3535533905932737*tmp[0]*bmag[3]*wgt[3]+0.3535533905932737*wgt[0]*bmag[3]*tmp[3]+0.3535533905932737*bmag[0]*tmp[2]*wgt[2]+0.3535533905932737*tmp[0]*bmag[2]*wgt[2]+0.3535533905932737*wgt[0]*bmag[2]*tmp[2]+0.3535533905932737*bmag[0]*tmp[1]*wgt[1]+0.3535533905932737*tmp[0]*bmag[1]*wgt[1]+0.3535533905932737*wgt[0]*bmag[1]*tmp[1]+0.3535533905932737*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

