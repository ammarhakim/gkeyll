#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (2.0*wgt[7]*f[20]+2.0*wgt[6]*f[19]+2.0*wgt[5]*f[12]+2.0*wgt[4]*f[11]+2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 

  out[0] += (1.4142135623730951*vmap[1]*wgt[7]*f[33]+1.4142135623730951*vmap[1]*wgt[6]*f[32]+1.4142135623730951*vmap[1]*wgt[5]*f[22]+1.4142135623730951*vmap[1]*wgt[4]*f[21]+1.4142135623730951*vmap[0]*wgt[7]*f[20]+1.4142135623730951*vmap[0]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[3]*f[15]+1.4142135623730951*vmap[0]*wgt[5]*f[12]+1.4142135623730951*vmap[0]*wgt[4]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*wgt[3]*f[34]+2.0000000000000004*vmap[0]*vmap[1]*wgt[7]*f[33]+2.0000000000000004*vmap[0]*vmap[1]*wgt[6]*f[32]+0.8944271909999161*vmap1R2*wgt[2]*f[24]+0.8944271909999161*vmap1R2*wgt[1]*f[23]+2.0000000000000004*vmap[0]*vmap[1]*wgt[5]*f[22]+2.0000000000000004*vmap[0]*vmap[1]*wgt[4]*f[21]+vmap1R2*wgt[7]*f[20]+vmap0R2*wgt[7]*f[20]+vmap1R2*wgt[6]*f[19]+vmap0R2*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[3]*f[15]+0.8944271909999159*wgt[0]*vmap1R2*f[13]+vmap1R2*wgt[5]*f[12]+vmap0R2*wgt[5]*f[12]+vmap1R2*wgt[4]*f[11]+vmap0R2*wgt[4]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*wgt[3]*f[34]+2.0000000000000004*vmap[0]*vmap[1]*wgt[7]*f[33]+2.0000000000000004*vmap[0]*vmap[1]*wgt[6]*f[32]+0.8944271909999161*vmap1R2*wgt[2]*f[24]+0.8944271909999161*vmap1R2*wgt[1]*f[23]+2.0000000000000004*vmap[0]*vmap[1]*wgt[5]*f[22]+2.0000000000000004*vmap[0]*vmap[1]*wgt[4]*f[21]+vmap1R2*wgt[7]*f[20]+vmap0R2*wgt[7]*f[20]+vmap1R2*wgt[6]*f[19]+vmap0R2*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[3]*f[15]+0.8944271909999159*wgt[0]*vmap1R2*f[13]+vmap1R2*wgt[5]*f[12]+vmap0R2*wgt[5]*f[12]+vmap1R2*wgt[4]*f[11]+vmap0R2*wgt[4]*f[11]+0.31943828249996997*bmag[5]*tmp[7]*wgt[7]+0.4472135954999579*bmag[4]*tmp[7]*wgt[7]+0.5*bmag[0]*tmp[7]*wgt[7]+0.31943828249996997*tmp[5]*bmag[7]*wgt[7]+0.4472135954999579*tmp[4]*bmag[7]*wgt[7]+0.5*tmp[0]*bmag[7]*wgt[7]+0.4*bmag[3]*tmp[6]*wgt[7]+0.4*tmp[3]*bmag[6]*wgt[7]+0.5000000000000001*bmag[1]*tmp[5]*wgt[7]+0.5000000000000001*tmp[1]*bmag[5]*wgt[7]+0.44721359549995804*bmag[2]*tmp[3]*wgt[7]+0.44721359549995804*tmp[2]*bmag[3]*wgt[7]+0.31943828249996997*wgt[5]*bmag[7]*tmp[7]+0.4472135954999579*wgt[4]*bmag[7]*tmp[7]+0.5*wgt[0]*bmag[7]*tmp[7]+0.4*bmag[3]*wgt[6]*tmp[7]+0.4*wgt[3]*bmag[6]*tmp[7]+0.5000000000000001*bmag[1]*wgt[5]*tmp[7]+0.5000000000000001*wgt[1]*bmag[5]*tmp[7]+0.44721359549995804*bmag[2]*wgt[3]*tmp[7]+0.44721359549995804*wgt[2]*bmag[3]*tmp[7]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+0.4*tmp[3]*wgt[6]*bmag[7]+0.4*wgt[3]*tmp[6]*bmag[7]+0.5000000000000001*tmp[1]*wgt[5]*bmag[7]+0.5000000000000001*wgt[1]*tmp[5]*bmag[7]+0.44721359549995804*tmp[2]*wgt[3]*bmag[7]+0.44721359549995804*wgt[2]*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]*wgt[6]+0.31943828249996997*bmag[4]*tmp[6]*wgt[6]+0.5*bmag[0]*tmp[6]*wgt[6]+0.4472135954999579*tmp[5]*bmag[6]*wgt[6]+0.31943828249996997*tmp[4]*bmag[6]*wgt[6]+0.5*tmp[0]*bmag[6]*wgt[6]+0.5000000000000001*bmag[2]*tmp[4]*wgt[6]+0.5000000000000001*tmp[2]*bmag[4]*wgt[6]+0.44721359549995804*bmag[1]*tmp[3]*wgt[6]+0.44721359549995804*tmp[1]*bmag[3]*wgt[6]+0.4472135954999579*wgt[5]*bmag[6]*tmp[6]+0.31943828249996997*wgt[4]*bmag[6]*tmp[6]+0.5*wgt[0]*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*wgt[4]*tmp[6]+0.5000000000000001*wgt[2]*bmag[4]*tmp[6]+0.44721359549995804*bmag[1]*wgt[3]*tmp[6]+0.44721359549995804*wgt[1]*bmag[3]*tmp[6]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+0.5000000000000001*tmp[2]*wgt[4]*bmag[6]+0.5000000000000001*wgt[2]*tmp[4]*bmag[6]+0.44721359549995804*tmp[1]*wgt[3]*bmag[6]+0.44721359549995804*wgt[1]*tmp[3]*bmag[6]+0.31943828249996997*bmag[5]*tmp[5]*wgt[5]+0.5*bmag[0]*tmp[5]*wgt[5]+0.5*tmp[0]*bmag[5]*wgt[5]+0.4472135954999579*bmag[3]*tmp[3]*wgt[5]+0.4472135954999579*bmag[2]*tmp[2]*wgt[5]+0.5*wgt[0]*bmag[5]*tmp[5]+0.4472135954999579*bmag[3]*wgt[3]*tmp[5]+0.4472135954999579*bmag[2]*wgt[2]*tmp[5]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+0.4472135954999579*tmp[3]*wgt[3]*bmag[5]+0.4472135954999579*tmp[2]*wgt[2]*bmag[5]+0.31943828249996997*bmag[4]*tmp[4]*wgt[4]+0.5*bmag[0]*tmp[4]*wgt[4]+0.5*tmp[0]*bmag[4]*wgt[4]+0.4472135954999579*bmag[3]*tmp[3]*wgt[4]+0.4472135954999579*bmag[1]*tmp[1]*wgt[4]+0.5*wgt[0]*bmag[4]*tmp[4]+0.4472135954999579*bmag[3]*wgt[3]*tmp[4]+0.4472135954999579*bmag[1]*wgt[1]*tmp[4]+0.4472135954999579*tmp[3]*wgt[3]*bmag[4]+0.4472135954999579*tmp[1]*wgt[1]*bmag[4]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.5*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.5*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[7]*f[20]+2.0*wgt[6]*f[19]+2.0*wgt[5]*f[12]+2.0*wgt[4]*f[11]+2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[7]*f[33]+1.4142135623730951*vmap[1]*wgt[6]*f[32]+1.4142135623730951*vmap[1]*wgt[5]*f[22]+1.4142135623730951*vmap[1]*wgt[4]*f[21]+1.4142135623730951*vmap[0]*wgt[7]*f[20]+1.4142135623730951*vmap[0]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[3]*f[15]+1.4142135623730951*vmap[0]*wgt[5]*f[12]+1.4142135623730951*vmap[0]*wgt[4]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*wgt[3]*f[34]+2.0000000000000004*vmap[0]*vmap[1]*wgt[7]*f[33]+2.0000000000000004*vmap[0]*vmap[1]*wgt[6]*f[32]+0.8944271909999161*vmap1R2*wgt[2]*f[24]+0.8944271909999161*vmap1R2*wgt[1]*f[23]+2.0000000000000004*vmap[0]*vmap[1]*wgt[5]*f[22]+2.0000000000000004*vmap[0]*vmap[1]*wgt[4]*f[21]+vmap1R2*wgt[7]*f[20]+vmap0R2*wgt[7]*f[20]+vmap1R2*wgt[6]*f[19]+vmap0R2*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[3]*f[15]+0.8944271909999159*wgt[0]*vmap1R2*f[13]+vmap1R2*wgt[5]*f[12]+vmap0R2*wgt[5]*f[12]+vmap1R2*wgt[4]*f[11]+vmap0R2*wgt[4]*f[11]+0.31943828249996997*bmag[5]*tmp[7]*wgt[7]+0.4472135954999579*bmag[4]*tmp[7]*wgt[7]+0.5*bmag[0]*tmp[7]*wgt[7]+0.31943828249996997*tmp[5]*bmag[7]*wgt[7]+0.4472135954999579*tmp[4]*bmag[7]*wgt[7]+0.5*tmp[0]*bmag[7]*wgt[7]+0.4*bmag[3]*tmp[6]*wgt[7]+0.4*tmp[3]*bmag[6]*wgt[7]+0.5000000000000001*bmag[1]*tmp[5]*wgt[7]+0.5000000000000001*tmp[1]*bmag[5]*wgt[7]+0.44721359549995804*bmag[2]*tmp[3]*wgt[7]+0.44721359549995804*tmp[2]*bmag[3]*wgt[7]+0.31943828249996997*wgt[5]*bmag[7]*tmp[7]+0.4472135954999579*wgt[4]*bmag[7]*tmp[7]+0.5*wgt[0]*bmag[7]*tmp[7]+0.4*bmag[3]*wgt[6]*tmp[7]+0.4*wgt[3]*bmag[6]*tmp[7]+0.5000000000000001*bmag[1]*wgt[5]*tmp[7]+0.5000000000000001*wgt[1]*bmag[5]*tmp[7]+0.44721359549995804*bmag[2]*wgt[3]*tmp[7]+0.44721359549995804*wgt[2]*bmag[3]*tmp[7]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+0.4*tmp[3]*wgt[6]*bmag[7]+0.4*wgt[3]*tmp[6]*bmag[7]+0.5000000000000001*tmp[1]*wgt[5]*bmag[7]+0.5000000000000001*wgt[1]*tmp[5]*bmag[7]+0.44721359549995804*tmp[2]*wgt[3]*bmag[7]+0.44721359549995804*wgt[2]*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]*wgt[6]+0.31943828249996997*bmag[4]*tmp[6]*wgt[6]+0.5*bmag[0]*tmp[6]*wgt[6]+0.4472135954999579*tmp[5]*bmag[6]*wgt[6]+0.31943828249996997*tmp[4]*bmag[6]*wgt[6]+0.5*tmp[0]*bmag[6]*wgt[6]+0.5000000000000001*bmag[2]*tmp[4]*wgt[6]+0.5000000000000001*tmp[2]*bmag[4]*wgt[6]+0.44721359549995804*bmag[1]*tmp[3]*wgt[6]+0.44721359549995804*tmp[1]*bmag[3]*wgt[6]+0.4472135954999579*wgt[5]*bmag[6]*tmp[6]+0.31943828249996997*wgt[4]*bmag[6]*tmp[6]+0.5*wgt[0]*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*wgt[4]*tmp[6]+0.5000000000000001*wgt[2]*bmag[4]*tmp[6]+0.44721359549995804*bmag[1]*wgt[3]*tmp[6]+0.44721359549995804*wgt[1]*bmag[3]*tmp[6]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+0.5000000000000001*tmp[2]*wgt[4]*bmag[6]+0.5000000000000001*wgt[2]*tmp[4]*bmag[6]+0.44721359549995804*tmp[1]*wgt[3]*bmag[6]+0.44721359549995804*wgt[1]*tmp[3]*bmag[6]+0.31943828249996997*bmag[5]*tmp[5]*wgt[5]+0.5*bmag[0]*tmp[5]*wgt[5]+0.5*tmp[0]*bmag[5]*wgt[5]+0.4472135954999579*bmag[3]*tmp[3]*wgt[5]+0.4472135954999579*bmag[2]*tmp[2]*wgt[5]+0.5*wgt[0]*bmag[5]*tmp[5]+0.4472135954999579*bmag[3]*wgt[3]*tmp[5]+0.4472135954999579*bmag[2]*wgt[2]*tmp[5]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+0.4472135954999579*tmp[3]*wgt[3]*bmag[5]+0.4472135954999579*tmp[2]*wgt[2]*bmag[5]+0.31943828249996997*bmag[4]*tmp[4]*wgt[4]+0.5*bmag[0]*tmp[4]*wgt[4]+0.5*tmp[0]*bmag[4]*wgt[4]+0.4472135954999579*bmag[3]*tmp[3]*wgt[4]+0.4472135954999579*bmag[1]*tmp[1]*wgt[4]+0.5*wgt[0]*bmag[4]*tmp[4]+0.4472135954999579*bmag[3]*wgt[3]*tmp[4]+0.4472135954999579*bmag[1]*wgt[1]*tmp[4]+0.4472135954999579*tmp[3]*wgt[3]*bmag[4]+0.4472135954999579*tmp[1]*wgt[1]*bmag[4]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+f[1]*vmap1R2*wgt[1]+0.5*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.5*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[16]+1.4142135623730951*vmap[2]*f[5]; 
  tmp[4] = 1.4142135623730951*vmap[3]*f[25]+1.4142135623730951*vmap[2]*f[11]; 
  tmp[5] = 1.4142135623730951*vmap[3]*f[26]+1.4142135623730951*vmap[2]*f[12]; 
  tmp[6] = 1.4142135623730951*vmap[3]*f[35]+1.4142135623730951*vmap[2]*f[19]; 
  tmp[7] = 1.4142135623730951*vmap[3]*f[36]+1.4142135623730951*vmap[2]*f[20]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*m_*volFact; 
  out[2] += (2.0*phi[7]*f[20]*q_+2.0*phi[6]*f[19]*q_+2.0*phi[5]*f[12]*q_+2.0*phi[4]*f[11]*q_+2.0*phi[3]*f[5]*q_+2.0*f[2]*phi[2]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact+(0.8944271909999159*vmap1R2*f[13]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*m_*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 

  out[0] += (0.31943828249996997*bmag[5]*tmp[7]*wgt[7]+0.4472135954999579*bmag[4]*tmp[7]*wgt[7]+0.5*bmag[0]*tmp[7]*wgt[7]+0.31943828249996997*tmp[5]*bmag[7]*wgt[7]+0.4472135954999579*tmp[4]*bmag[7]*wgt[7]+0.5*tmp[0]*bmag[7]*wgt[7]+0.4*bmag[3]*tmp[6]*wgt[7]+0.4*tmp[3]*bmag[6]*wgt[7]+0.5000000000000001*bmag[1]*tmp[5]*wgt[7]+0.5000000000000001*tmp[1]*bmag[5]*wgt[7]+0.44721359549995804*bmag[2]*tmp[3]*wgt[7]+0.44721359549995804*tmp[2]*bmag[3]*wgt[7]+0.31943828249996997*wgt[5]*bmag[7]*tmp[7]+0.4472135954999579*wgt[4]*bmag[7]*tmp[7]+0.5*wgt[0]*bmag[7]*tmp[7]+0.4*bmag[3]*wgt[6]*tmp[7]+0.4*wgt[3]*bmag[6]*tmp[7]+0.5000000000000001*bmag[1]*wgt[5]*tmp[7]+0.5000000000000001*wgt[1]*bmag[5]*tmp[7]+0.44721359549995804*bmag[2]*wgt[3]*tmp[7]+0.44721359549995804*wgt[2]*bmag[3]*tmp[7]+0.4*tmp[3]*wgt[6]*bmag[7]+0.4*wgt[3]*tmp[6]*bmag[7]+0.5000000000000001*tmp[1]*wgt[5]*bmag[7]+0.5000000000000001*wgt[1]*tmp[5]*bmag[7]+0.44721359549995804*tmp[2]*wgt[3]*bmag[7]+0.44721359549995804*wgt[2]*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]*wgt[6]+0.31943828249996997*bmag[4]*tmp[6]*wgt[6]+0.5*bmag[0]*tmp[6]*wgt[6]+0.4472135954999579*tmp[5]*bmag[6]*wgt[6]+0.31943828249996997*tmp[4]*bmag[6]*wgt[6]+0.5*tmp[0]*bmag[6]*wgt[6]+0.5000000000000001*bmag[2]*tmp[4]*wgt[6]+0.5000000000000001*tmp[2]*bmag[4]*wgt[6]+0.44721359549995804*bmag[1]*tmp[3]*wgt[6]+0.44721359549995804*tmp[1]*bmag[3]*wgt[6]+0.4472135954999579*wgt[5]*bmag[6]*tmp[6]+0.31943828249996997*wgt[4]*bmag[6]*tmp[6]+0.5*wgt[0]*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*wgt[4]*tmp[6]+0.5000000000000001*wgt[2]*bmag[4]*tmp[6]+0.44721359549995804*bmag[1]*wgt[3]*tmp[6]+0.44721359549995804*wgt[1]*bmag[3]*tmp[6]+0.5000000000000001*tmp[2]*wgt[4]*bmag[6]+0.5000000000000001*wgt[2]*tmp[4]*bmag[6]+0.44721359549995804*tmp[1]*wgt[3]*bmag[6]+0.44721359549995804*wgt[1]*tmp[3]*bmag[6]+0.31943828249996997*bmag[5]*tmp[5]*wgt[5]+0.5*bmag[0]*tmp[5]*wgt[5]+0.5*tmp[0]*bmag[5]*wgt[5]+0.4472135954999579*bmag[3]*tmp[3]*wgt[5]+0.4472135954999579*bmag[2]*tmp[2]*wgt[5]+0.5*wgt[0]*bmag[5]*tmp[5]+0.4472135954999579*bmag[3]*wgt[3]*tmp[5]+0.4472135954999579*bmag[2]*wgt[2]*tmp[5]+0.4472135954999579*tmp[3]*wgt[3]*bmag[5]+0.4472135954999579*tmp[2]*wgt[2]*bmag[5]+0.31943828249996997*bmag[4]*tmp[4]*wgt[4]+0.5*bmag[0]*tmp[4]*wgt[4]+0.5*tmp[0]*bmag[4]*wgt[4]+0.4472135954999579*bmag[3]*tmp[3]*wgt[4]+0.4472135954999579*bmag[1]*tmp[1]*wgt[4]+0.5*wgt[0]*bmag[4]*tmp[4]+0.4472135954999579*bmag[3]*wgt[3]*tmp[4]+0.4472135954999579*bmag[1]*wgt[1]*tmp[4]+0.4472135954999579*tmp[3]*wgt[3]*bmag[4]+0.4472135954999579*tmp[1]*wgt[1]*bmag[4]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[7]*f[20]+2.0*wgt[6]*f[19]+2.0*wgt[5]*f[12]+2.0*wgt[4]*f[11]+2.0*wgt[3]*f[5]+2.0*f[2]*wgt[2]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[7]*f[33]+1.4142135623730951*vmap[1]*wgt[6]*f[32]+1.4142135623730951*vmap[1]*wgt[5]*f[22]+1.4142135623730951*vmap[1]*wgt[4]*f[21]+1.4142135623730951*vmap[0]*wgt[7]*f[20]+1.4142135623730951*vmap[0]*wgt[6]*f[19]+1.4142135623730951*vmap[1]*wgt[3]*f[15]+1.4142135623730951*vmap[0]*wgt[5]*f[12]+1.4142135623730951*vmap[0]*wgt[4]*f[11]+1.4142135623730951*vmap[1]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[6]+1.4142135623730951*vmap[0]*wgt[3]*f[5]+1.4142135623730951*wgt[0]*vmap[1]*f[3]+1.4142135623730951*vmap[0]*f[2]*wgt[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*wgt[3]*f[34]+2.0000000000000004*vmap[0]*vmap[1]*wgt[7]*f[33]+2.0000000000000004*vmap[0]*vmap[1]*wgt[6]*f[32]+0.8944271909999161*vmap1R2*wgt[2]*f[24]+0.8944271909999161*vmap1R2*wgt[1]*f[23]+2.0000000000000004*vmap[0]*vmap[1]*wgt[5]*f[22]+2.0000000000000004*vmap[0]*vmap[1]*wgt[4]*f[21]+vmap1R2*wgt[7]*f[20]+vmap0R2*wgt[7]*f[20]+vmap1R2*wgt[6]*f[19]+vmap0R2*wgt[6]*f[19]+2.0*vmap[0]*vmap[1]*wgt[3]*f[15]+0.8944271909999159*wgt[0]*vmap1R2*f[13]+vmap1R2*wgt[5]*f[12]+vmap0R2*wgt[5]*f[12]+vmap1R2*wgt[4]*f[11]+vmap0R2*wgt[4]*f[11]+2.0*vmap[0]*vmap[1]*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[6]+vmap1R2*wgt[3]*f[5]+vmap0R2*wgt[3]*f[5]+2.0*vmap[0]*wgt[0]*vmap[1]*f[3]+vmap1R2*f[2]*wgt[2]+vmap0R2*f[2]*wgt[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (0.31943828249996997*bmag[5]*tmp[7]*wgt[7]+0.4472135954999579*bmag[4]*tmp[7]*wgt[7]+0.5*bmag[0]*tmp[7]*wgt[7]+0.31943828249996997*tmp[5]*bmag[7]*wgt[7]+0.4472135954999579*tmp[4]*bmag[7]*wgt[7]+0.5*tmp[0]*bmag[7]*wgt[7]+0.4*bmag[3]*tmp[6]*wgt[7]+0.4*tmp[3]*bmag[6]*wgt[7]+0.5000000000000001*bmag[1]*tmp[5]*wgt[7]+0.5000000000000001*tmp[1]*bmag[5]*wgt[7]+0.44721359549995804*bmag[2]*tmp[3]*wgt[7]+0.44721359549995804*tmp[2]*bmag[3]*wgt[7]+0.31943828249996997*wgt[5]*bmag[7]*tmp[7]+0.4472135954999579*wgt[4]*bmag[7]*tmp[7]+0.5*wgt[0]*bmag[7]*tmp[7]+0.4*bmag[3]*wgt[6]*tmp[7]+0.4*wgt[3]*bmag[6]*tmp[7]+0.5000000000000001*bmag[1]*wgt[5]*tmp[7]+0.5000000000000001*wgt[1]*bmag[5]*tmp[7]+0.44721359549995804*bmag[2]*wgt[3]*tmp[7]+0.44721359549995804*wgt[2]*bmag[3]*tmp[7]+0.4*tmp[3]*wgt[6]*bmag[7]+0.4*wgt[3]*tmp[6]*bmag[7]+0.5000000000000001*tmp[1]*wgt[5]*bmag[7]+0.5000000000000001*wgt[1]*tmp[5]*bmag[7]+0.44721359549995804*tmp[2]*wgt[3]*bmag[7]+0.44721359549995804*wgt[2]*tmp[3]*bmag[7]+0.4472135954999579*bmag[5]*tmp[6]*wgt[6]+0.31943828249996997*bmag[4]*tmp[6]*wgt[6]+0.5*bmag[0]*tmp[6]*wgt[6]+0.4472135954999579*tmp[5]*bmag[6]*wgt[6]+0.31943828249996997*tmp[4]*bmag[6]*wgt[6]+0.5*tmp[0]*bmag[6]*wgt[6]+0.5000000000000001*bmag[2]*tmp[4]*wgt[6]+0.5000000000000001*tmp[2]*bmag[4]*wgt[6]+0.44721359549995804*bmag[1]*tmp[3]*wgt[6]+0.44721359549995804*tmp[1]*bmag[3]*wgt[6]+0.4472135954999579*wgt[5]*bmag[6]*tmp[6]+0.31943828249996997*wgt[4]*bmag[6]*tmp[6]+0.5*wgt[0]*bmag[6]*tmp[6]+0.5000000000000001*bmag[2]*wgt[4]*tmp[6]+0.5000000000000001*wgt[2]*bmag[4]*tmp[6]+0.44721359549995804*bmag[1]*wgt[3]*tmp[6]+0.44721359549995804*wgt[1]*bmag[3]*tmp[6]+0.5000000000000001*tmp[2]*wgt[4]*bmag[6]+0.5000000000000001*wgt[2]*tmp[4]*bmag[6]+0.44721359549995804*tmp[1]*wgt[3]*bmag[6]+0.44721359549995804*wgt[1]*tmp[3]*bmag[6]+0.31943828249996997*bmag[5]*tmp[5]*wgt[5]+0.5*bmag[0]*tmp[5]*wgt[5]+0.5*tmp[0]*bmag[5]*wgt[5]+0.4472135954999579*bmag[3]*tmp[3]*wgt[5]+0.4472135954999579*bmag[2]*tmp[2]*wgt[5]+0.5*wgt[0]*bmag[5]*tmp[5]+0.4472135954999579*bmag[3]*wgt[3]*tmp[5]+0.4472135954999579*bmag[2]*wgt[2]*tmp[5]+0.4472135954999579*tmp[3]*wgt[3]*bmag[5]+0.4472135954999579*tmp[2]*wgt[2]*bmag[5]+0.31943828249996997*bmag[4]*tmp[4]*wgt[4]+0.5*bmag[0]*tmp[4]*wgt[4]+0.5*tmp[0]*bmag[4]*wgt[4]+0.4472135954999579*bmag[3]*tmp[3]*wgt[4]+0.4472135954999579*bmag[1]*tmp[1]*wgt[4]+0.5*wgt[0]*bmag[4]*tmp[4]+0.4472135954999579*bmag[3]*wgt[3]*tmp[4]+0.4472135954999579*bmag[1]*wgt[1]*tmp[4]+0.4472135954999579*tmp[3]*wgt[3]*bmag[4]+0.4472135954999579*tmp[1]*wgt[1]*bmag[4]+0.5*bmag[0]*tmp[3]*wgt[3]+0.5*tmp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp[2]*wgt[3]+0.5*tmp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp[3]+0.5*bmag[1]*wgt[2]*tmp[3]+0.5*wgt[1]*bmag[2]*tmp[3]+0.5*tmp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

