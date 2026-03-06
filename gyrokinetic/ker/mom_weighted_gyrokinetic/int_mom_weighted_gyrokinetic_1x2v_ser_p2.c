#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (2.0*wgt[2]*f[7]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (1.4142135623730951*vmap[1]*wgt[2]*f[11]+1.4142135623730951*vmap[0]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[1]*f[12]+2.0000000000000004*vmap[0]*vmap[1]*wgt[2]*f[11]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+vmap1R2*wgt[2]*f[7]+vmap0R2*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999161*vmap1R2*wgt[1]*f[12]+2.0000000000000004*vmap[0]*vmap[1]*wgt[2]*f[11]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+vmap1R2*wgt[2]*f[7]+vmap0R2*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+0.45175395145262565*bmag[2]*tmp[2]*wgt[2]+0.7071067811865475*bmag[0]*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*bmag[2]*wgt[2]+0.6324555320336759*bmag[1]*tmp[1]*wgt[2]+0.7071067811865475*wgt[0]*bmag[2]*tmp[2]+0.6324555320336759*bmag[1]*wgt[1]*tmp[2]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+0.6324555320336759*tmp[1]*wgt[1]*bmag[2]+f[1]*vmap1R2*wgt[1]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[2]*f[7]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[2]*f[11]+1.4142135623730951*vmap[0]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[1]*f[12]+2.0000000000000004*vmap[0]*vmap[1]*wgt[2]*f[11]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+vmap1R2*wgt[2]*f[7]+vmap0R2*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+0.45175395145262565*bmag[2]*tmp[2]*wgt[2]+0.7071067811865475*bmag[0]*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*bmag[2]*wgt[2]+0.6324555320336759*bmag[1]*tmp[1]*wgt[2]+0.7071067811865475*wgt[0]*bmag[2]*tmp[2]+0.6324555320336759*bmag[1]*wgt[1]*tmp[2]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+0.6324555320336759*tmp[1]*wgt[1]*bmag[2]+f[1]*vmap1R2*wgt[1]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+vmap0R2*f[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+f[0]*vmap0R2*wgt[0]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.0*vmap[1]*f[2]+2.0*f[0]*vmap[0])*m_*volFact; 
  out[2] += (2.0*phi[2]*f[7]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact+(0.6324555320336759*vmap1R2*f[8]+1.4142135623730951*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2)*m_*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 

  out[0] += (0.45175395145262565*bmag[2]*tmp[2]*wgt[2]+0.7071067811865475*bmag[0]*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*bmag[2]*wgt[2]+0.6324555320336759*bmag[1]*tmp[1]*wgt[2]+0.7071067811865475*wgt[0]*bmag[2]*tmp[2]+0.6324555320336759*bmag[1]*wgt[1]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1]*bmag[2]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (2.0*wgt[2]*f[7]+2.0*f[1]*wgt[1]+2.0*f[0]*wgt[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*wgt[2]*f[11]+1.4142135623730951*vmap[0]*wgt[2]*f[7]+1.4142135623730951*vmap[1]*wgt[1]*f[4]+1.4142135623730951*wgt[0]*vmap[1]*f[2]+1.4142135623730951*vmap[0]*f[1]*wgt[1]+1.4142135623730951*f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*wgt[1]*f[12]+2.0000000000000004*vmap[0]*vmap[1]*wgt[2]*f[11]+0.8944271909999159*wgt[0]*vmap1R2*f[8]+vmap1R2*wgt[2]*f[7]+vmap0R2*wgt[2]*f[7]+2.0*vmap[0]*vmap[1]*wgt[1]*f[4]+2.0*vmap[0]*wgt[0]*vmap[1]*f[2]+f[1]*vmap1R2*wgt[1]+vmap0R2*f[1]*wgt[1]+f[0]*wgt[0]*vmap1R2+f[0]*vmap0R2*wgt[0])*volFact; 
  out[3] += (0.45175395145262565*bmag[2]*tmp[2]*wgt[2]+0.7071067811865475*bmag[0]*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*bmag[2]*wgt[2]+0.6324555320336759*bmag[1]*tmp[1]*wgt[2]+0.7071067811865475*wgt[0]*bmag[2]*tmp[2]+0.6324555320336759*bmag[1]*wgt[1]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1]*bmag[2]+0.7071067811865475*bmag[0]*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*bmag[1]*wgt[1]+0.7071067811865475*wgt[0]*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0]*wgt[0])*volFact; 
} 

