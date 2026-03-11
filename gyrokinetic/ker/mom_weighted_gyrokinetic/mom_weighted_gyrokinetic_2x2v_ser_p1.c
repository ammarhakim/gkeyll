#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  out[0] += (wgt[3]*f[5]+f[2]*wgt[2]+f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (wgt[2]*f[5]+f[2]*wgt[3]+f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (wgt[1]*f[5]+f[1]*wgt[3]+f[0]*wgt[2]+wgt[0]*f[2])*volFact; 
  out[3] += (wgt[0]*f[5]+f[0]*wgt[3]+f[1]*wgt[2]+wgt[1]*f[2])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0]; 
  tmp[1] = 1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1]; 
  tmp[2] = 1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2]; 
  tmp[3] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5]; 
  out[0] += (0.5*tmp[3]*wgt[3]+0.5*tmp[2]*wgt[2]+0.5*tmp[1]*wgt[1]+0.5*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.5*tmp[2]*wgt[3]+0.5*wgt[2]*tmp[3]+0.5*tmp[0]*wgt[1]+0.5*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.5*tmp[1]*wgt[3]+0.5*wgt[1]*tmp[3]+0.5*tmp[0]*wgt[2]+0.5*wgt[0]*tmp[2])*volFact; 
  out[3] += (0.5*tmp[0]*wgt[3]+0.5*wgt[0]*tmp[3]+0.5*tmp[1]*wgt[2]+0.5*wgt[1]*tmp[2])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmpPar[4]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmpPar[1] = 0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmpPar[2] = 0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2]; 
  tmpPar[3] = 0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5]; 

  out[0] += (0.5*tmpPar[3]*wgt[3]+0.5*tmpPar[2]*wgt[2]+0.5*tmpPar[1]*wgt[1]+0.5*tmpPar[0]*wgt[0])*volFact; 
  out[1] += (0.5*tmpPar[2]*wgt[3]+0.5*wgt[2]*tmpPar[3]+0.5*tmpPar[0]*wgt[1]+0.5*wgt[0]*tmpPar[1])*volFact; 
  out[2] += (0.5*tmpPar[1]*wgt[3]+0.5*wgt[1]*tmpPar[3]+0.5*tmpPar[0]*wgt[2]+0.5*wgt[0]*tmpPar[2])*volFact; 
  out[3] += (0.5*tmpPar[0]*wgt[3]+0.5*wgt[0]*tmpPar[3]+0.5*tmpPar[1]*wgt[2]+0.5*wgt[1]*tmpPar[2])*volFact; 
  double tmpPerp[4]; 
  tmpPerp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmpPerp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmpPerp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmpPerp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  out[0] += (2.0*(0.25*bmag[0]*tmpPerp[3]*wgt[3]+0.25*tmpPerp[0]*bmag[3]*wgt[3]+0.25*bmag[1]*tmpPerp[2]*wgt[3]+0.25*tmpPerp[1]*bmag[2]*wgt[3]+0.25*wgt[0]*bmag[3]*tmpPerp[3]+0.25*bmag[1]*wgt[2]*tmpPerp[3]+0.25*wgt[1]*bmag[2]*tmpPerp[3]+0.25*tmpPerp[1]*wgt[2]*bmag[3]+0.25*wgt[1]*tmpPerp[2]*bmag[3]+0.25*bmag[0]*tmpPerp[2]*wgt[2]+0.25*tmpPerp[0]*bmag[2]*wgt[2]+0.25*wgt[0]*bmag[2]*tmpPerp[2]+0.25*bmag[0]*tmpPerp[1]*wgt[1]+0.25*tmpPerp[0]*bmag[1]*wgt[1]+0.25*wgt[0]*bmag[1]*tmpPerp[1]+0.25*bmag[0]*tmpPerp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.45*bmag[1]*tmpPerp[3]*wgt[3]+0.45*tmpPerp[1]*bmag[3]*wgt[3]+0.25*bmag[0]*tmpPerp[2]*wgt[3]+0.25*tmpPerp[0]*bmag[2]*wgt[3]+0.45*wgt[1]*bmag[3]*tmpPerp[3]+0.25*bmag[0]*wgt[2]*tmpPerp[3]+0.25*wgt[0]*bmag[2]*tmpPerp[3]+0.25*tmpPerp[0]*wgt[2]*bmag[3]+0.25*wgt[0]*tmpPerp[2]*bmag[3]+0.25*bmag[1]*tmpPerp[2]*wgt[2]+0.25*tmpPerp[1]*bmag[2]*wgt[2]+0.25*wgt[1]*bmag[2]*tmpPerp[2]+0.45*bmag[1]*tmpPerp[1]*wgt[1]+0.25*bmag[0]*tmpPerp[0]*wgt[1]+0.25*bmag[0]*wgt[0]*tmpPerp[1]+0.25*tmpPerp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.45*bmag[2]*tmpPerp[3]*wgt[3]+0.45*tmpPerp[2]*bmag[3]*wgt[3]+0.25*bmag[0]*tmpPerp[1]*wgt[3]+0.25*tmpPerp[0]*bmag[1]*wgt[3]+0.45*wgt[2]*bmag[3]*tmpPerp[3]+0.25*bmag[0]*wgt[1]*tmpPerp[3]+0.25*wgt[0]*bmag[1]*tmpPerp[3]+0.25*tmpPerp[0]*wgt[1]*bmag[3]+0.25*wgt[0]*tmpPerp[1]*bmag[3]+0.45*bmag[2]*tmpPerp[2]*wgt[2]+0.25*bmag[1]*tmpPerp[1]*wgt[2]+0.25*bmag[0]*tmpPerp[0]*wgt[2]+0.25*bmag[1]*wgt[1]*tmpPerp[2]+0.25*bmag[0]*wgt[0]*tmpPerp[2]+0.25*tmpPerp[1]*wgt[1]*bmag[2]+0.25*tmpPerp[0]*wgt[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.81*bmag[3]*tmpPerp[3]*wgt[3]+0.45*bmag[2]*tmpPerp[2]*wgt[3]+0.45*bmag[1]*tmpPerp[1]*wgt[3]+0.25*bmag[0]*tmpPerp[0]*wgt[3]+0.45*bmag[2]*wgt[2]*tmpPerp[3]+0.45*bmag[1]*wgt[1]*tmpPerp[3]+0.25*bmag[0]*wgt[0]*tmpPerp[3]+0.45*tmpPerp[2]*wgt[2]*bmag[3]+0.45*tmpPerp[1]*wgt[1]*bmag[3]+0.25*tmpPerp[0]*wgt[0]*bmag[3]+0.25*bmag[0]*tmpPerp[1]*wgt[2]+0.25*tmpPerp[0]*bmag[1]*wgt[2]+0.25*bmag[0]*wgt[1]*tmpPerp[2]+0.25*wgt[0]*bmag[1]*tmpPerp[2]+0.25*tmpPerp[0]*wgt[1]*bmag[2]+0.25*wgt[0]*tmpPerp[1]*bmag[2])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp[1] = 0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp[2] = 0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2]; 
  tmp[3] = 0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5]; 

  out[0] += (0.5*tmp[3]*wgt[3]+0.5*tmp[2]*wgt[2]+0.5*tmp[1]*wgt[1]+0.5*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.5*tmp[2]*wgt[3]+0.5*wgt[2]*tmp[3]+0.5*tmp[0]*wgt[1]+0.5*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.5*tmp[1]*wgt[3]+0.5*wgt[1]*tmp[3]+0.5*tmp[0]*wgt[2]+0.5*wgt[0]*tmp[2])*volFact; 
  out[3] += (0.5*tmp[0]*wgt[3]+0.5*wgt[0]*tmp[3]+0.5*tmp[1]*wgt[2]+0.5*wgt[1]*tmp[2])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp1[4]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1]; 
  tmp1[2] = 1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2]; 
  tmp1[3] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5]; 
  double tmp2par[4]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp2par[2] = 0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2]; 
  tmp2par[3] = 0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5]; 
  double tmp2perp[4]; 
  tmp2perp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2perp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp2perp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 

  out[0] += (wgt[3]*f[5]+f[2]*wgt[2]+f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (wgt[2]*f[5]+f[2]*wgt[3]+f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (wgt[1]*f[5]+f[1]*wgt[3]+f[0]*wgt[2]+wgt[0]*f[2])*volFact; 
  out[3] += (wgt[0]*f[5]+f[0]*wgt[3]+f[1]*wgt[2]+wgt[1]*f[2])*volFact; 
  out[4] += (0.5*tmp1[3]*wgt[3]+0.5*tmp1[2]*wgt[2]+0.5*tmp1[1]*wgt[1]+0.5*tmp1[0]*wgt[0])*volFact; 
  out[5] += (0.5*tmp1[2]*wgt[3]+0.5*wgt[2]*tmp1[3]+0.5*tmp1[0]*wgt[1]+0.5*wgt[0]*tmp1[1])*volFact; 
  out[6] += (0.5*tmp1[1]*wgt[3]+0.5*wgt[1]*tmp1[3]+0.5*tmp1[0]*wgt[2]+0.5*wgt[0]*tmp1[2])*volFact; 
  out[7] += (0.5*tmp1[0]*wgt[3]+0.5*wgt[0]*tmp1[3]+0.5*tmp1[1]*wgt[2]+0.5*wgt[1]*tmp1[2])*volFact; 
  out[8] += ((0.5*bmag[0]*tmp2perp[3]*wgt[3]+0.5*tmp2perp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp2perp[2]*wgt[3]+0.5*tmp2perp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp2perp[3]+0.5*bmag[1]*wgt[2]*tmp2perp[3]+0.5*wgt[1]*bmag[2]*tmp2perp[3]+0.5*tmp2perp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp2perp[2]*bmag[3]+0.5*bmag[0]*tmp2perp[2]*wgt[2]+0.5*tmp2perp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp2perp[2]+0.5*bmag[0]*tmp2perp[1]*wgt[1]+0.5*tmp2perp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp2perp[1]+0.5*bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_+(0.5*tmp2par[3]*wgt[3]+0.5*tmp2par[2]*wgt[2]+0.5*tmp2par[1]*wgt[1]+0.5*tmp2par[0]*wgt[0])*volFact; 
  out[9] += ((0.9*bmag[1]*tmp2perp[3]*wgt[3]+0.9*tmp2perp[1]*bmag[3]*wgt[3]+0.5*bmag[0]*tmp2perp[2]*wgt[3]+0.5*tmp2perp[0]*bmag[2]*wgt[3]+0.9*wgt[1]*bmag[3]*tmp2perp[3]+0.5*bmag[0]*wgt[2]*tmp2perp[3]+0.5*wgt[0]*bmag[2]*tmp2perp[3]+0.5*tmp2perp[0]*wgt[2]*bmag[3]+0.5*wgt[0]*tmp2perp[2]*bmag[3]+0.5*bmag[1]*tmp2perp[2]*wgt[2]+0.5*tmp2perp[1]*bmag[2]*wgt[2]+0.5*wgt[1]*bmag[2]*tmp2perp[2]+0.9*bmag[1]*tmp2perp[1]*wgt[1]+0.5*bmag[0]*tmp2perp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmp2perp[1]+0.5*tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_+(0.5*tmp2par[2]*wgt[3]+0.5*wgt[2]*tmp2par[3]+0.5*tmp2par[0]*wgt[1]+0.5*wgt[0]*tmp2par[1])*volFact; 
  out[10] += ((0.9*bmag[2]*tmp2perp[3]*wgt[3]+0.9*tmp2perp[2]*bmag[3]*wgt[3]+0.5*bmag[0]*tmp2perp[1]*wgt[3]+0.5*tmp2perp[0]*bmag[1]*wgt[3]+0.9*wgt[2]*bmag[3]*tmp2perp[3]+0.5*bmag[0]*wgt[1]*tmp2perp[3]+0.5*wgt[0]*bmag[1]*tmp2perp[3]+0.5*tmp2perp[0]*wgt[1]*bmag[3]+0.5*wgt[0]*tmp2perp[1]*bmag[3]+0.9*bmag[2]*tmp2perp[2]*wgt[2]+0.5*bmag[1]*tmp2perp[1]*wgt[2]+0.5*bmag[0]*tmp2perp[0]*wgt[2]+0.5*bmag[1]*wgt[1]*tmp2perp[2]+0.5*bmag[0]*wgt[0]*tmp2perp[2]+0.5*tmp2perp[1]*wgt[1]*bmag[2]+0.5*tmp2perp[0]*wgt[0]*bmag[2])*volFact)/m_+(0.5*tmp2par[1]*wgt[3]+0.5*wgt[1]*tmp2par[3]+0.5*tmp2par[0]*wgt[2]+0.5*wgt[0]*tmp2par[2])*volFact; 
  out[11] += ((1.62*bmag[3]*tmp2perp[3]*wgt[3]+0.9*bmag[2]*tmp2perp[2]*wgt[3]+0.9*bmag[1]*tmp2perp[1]*wgt[3]+0.5*bmag[0]*tmp2perp[0]*wgt[3]+0.9*bmag[2]*wgt[2]*tmp2perp[3]+0.9*bmag[1]*wgt[1]*tmp2perp[3]+0.5*bmag[0]*wgt[0]*tmp2perp[3]+0.9*tmp2perp[2]*wgt[2]*bmag[3]+0.9*tmp2perp[1]*wgt[1]*bmag[3]+0.5*tmp2perp[0]*wgt[0]*bmag[3]+0.5*bmag[0]*tmp2perp[1]*wgt[2]+0.5*tmp2perp[0]*bmag[1]*wgt[2]+0.5*bmag[0]*wgt[1]*tmp2perp[2]+0.5*wgt[0]*bmag[1]*tmp2perp[2]+0.5*tmp2perp[0]*wgt[1]*bmag[2]+0.5*wgt[0]*tmp2perp[1]*bmag[2])*volFact)/m_+(0.5*tmp2par[0]*wgt[3]+0.5*wgt[0]*tmp2par[3]+0.5*tmp2par[1]*wgt[2]+0.5*wgt[1]*tmp2par[2])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*m_*volFact; 
  out[6] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*m_*volFact; 
  out[7] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5])*m_*volFact; 
  out[8] += (phi[3]*f[5]+f[2]*phi[2]+f[1]*phi[1]+f[0]*phi[0])*q_*volFact+(0.4472135954999579*vmap1R2*f[16]+vmap[0]*vmap[1]*f[3]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact+(0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact; 
  out[9] += (phi[2]*f[5]+f[2]*phi[3]+f[0]*phi[1]+phi[0]*f[1])*q_*volFact+(0.44721359549995804*vmap1R2*f[17]+vmap[0]*vmap[1]*f[6]+0.5*f[1]*vmap1R2+0.5*vmap0R2*f[1])*m_*volFact+(0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact; 
  out[10] += (phi[1]*f[5]+f[1]*phi[3]+f[0]*phi[2]+phi[0]*f[2])*q_*volFact+(0.44721359549995804*vmap1R2*f[18]+vmap[0]*vmap[1]*f[7]+0.5*vmap1R2*f[2]+0.5*vmap0R2*f[2])*m_*volFact+(0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact; 
  out[11] += (phi[0]*f[5]+f[0]*phi[3]+f[1]*phi[2]+phi[1]*f[2])*q_*volFact+(0.4472135954999579*vmap1R2*f[20]+vmap[0]*vmap[1]*f[11]+0.5*vmap1R2*f[5]+0.5*vmap0R2*f[5])*m_*volFact+(0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PERP_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  out[0] += (2.0*(0.25*bmag[0]*tmp[3]*wgt[3]+0.25*tmp[0]*bmag[3]*wgt[3]+0.25*bmag[1]*tmp[2]*wgt[3]+0.25*tmp[1]*bmag[2]*wgt[3]+0.25*wgt[0]*bmag[3]*tmp[3]+0.25*bmag[1]*wgt[2]*tmp[3]+0.25*wgt[1]*bmag[2]*tmp[3]+0.25*tmp[1]*wgt[2]*bmag[3]+0.25*wgt[1]*tmp[2]*bmag[3]+0.25*bmag[0]*tmp[2]*wgt[2]+0.25*tmp[0]*bmag[2]*wgt[2]+0.25*wgt[0]*bmag[2]*tmp[2]+0.25*bmag[0]*tmp[1]*wgt[1]+0.25*tmp[0]*bmag[1]*wgt[1]+0.25*wgt[0]*bmag[1]*tmp[1]+0.25*bmag[0]*tmp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.45*bmag[1]*tmp[3]*wgt[3]+0.45*tmp[1]*bmag[3]*wgt[3]+0.25*bmag[0]*tmp[2]*wgt[3]+0.25*tmp[0]*bmag[2]*wgt[3]+0.45*wgt[1]*bmag[3]*tmp[3]+0.25*bmag[0]*wgt[2]*tmp[3]+0.25*wgt[0]*bmag[2]*tmp[3]+0.25*tmp[0]*wgt[2]*bmag[3]+0.25*wgt[0]*tmp[2]*bmag[3]+0.25*bmag[1]*tmp[2]*wgt[2]+0.25*tmp[1]*bmag[2]*wgt[2]+0.25*wgt[1]*bmag[2]*tmp[2]+0.45*bmag[1]*tmp[1]*wgt[1]+0.25*bmag[0]*tmp[0]*wgt[1]+0.25*bmag[0]*wgt[0]*tmp[1]+0.25*tmp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.45*bmag[2]*tmp[3]*wgt[3]+0.45*tmp[2]*bmag[3]*wgt[3]+0.25*bmag[0]*tmp[1]*wgt[3]+0.25*tmp[0]*bmag[1]*wgt[3]+0.45*wgt[2]*bmag[3]*tmp[3]+0.25*bmag[0]*wgt[1]*tmp[3]+0.25*wgt[0]*bmag[1]*tmp[3]+0.25*tmp[0]*wgt[1]*bmag[3]+0.25*wgt[0]*tmp[1]*bmag[3]+0.45*bmag[2]*tmp[2]*wgt[2]+0.25*bmag[1]*tmp[1]*wgt[2]+0.25*bmag[0]*tmp[0]*wgt[2]+0.25*bmag[1]*wgt[1]*tmp[2]+0.25*bmag[0]*wgt[0]*tmp[2]+0.25*tmp[1]*wgt[1]*bmag[2]+0.25*tmp[0]*wgt[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.81*bmag[3]*tmp[3]*wgt[3]+0.45*bmag[2]*tmp[2]*wgt[3]+0.45*bmag[1]*tmp[1]*wgt[3]+0.25*bmag[0]*tmp[0]*wgt[3]+0.45*bmag[2]*wgt[2]*tmp[3]+0.45*bmag[1]*wgt[1]*tmp[3]+0.25*bmag[0]*wgt[0]*tmp[3]+0.45*tmp[2]*wgt[2]*bmag[3]+0.45*tmp[1]*wgt[1]*bmag[3]+0.25*tmp[0]*wgt[0]*bmag[3]+0.25*bmag[0]*tmp[1]*wgt[2]+0.25*tmp[0]*bmag[1]*wgt[2]+0.25*bmag[0]*wgt[1]*tmp[2]+0.25*wgt[0]*bmag[1]*tmp[2]+0.25*tmp[0]*wgt[1]*bmag[2]+0.25*wgt[0]*tmp[1]*bmag[2])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp1[4]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1]; 
  tmp1[2] = 1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2]; 
  tmp1[3] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5]; 
  double tmp2par[4]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp2par[2] = 0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2]; 
  tmp2par[3] = 0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5]; 
  double tmp2perp[4]; 
  tmp2perp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2perp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp2perp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 

  out[0] += (wgt[3]*f[5]+f[2]*wgt[2]+f[1]*wgt[1]+f[0]*wgt[0])*volFact; 
  out[1] += (wgt[2]*f[5]+f[2]*wgt[3]+f[0]*wgt[1]+wgt[0]*f[1])*volFact; 
  out[2] += (wgt[1]*f[5]+f[1]*wgt[3]+f[0]*wgt[2]+wgt[0]*f[2])*volFact; 
  out[3] += (wgt[0]*f[5]+f[0]*wgt[3]+f[1]*wgt[2]+wgt[1]*f[2])*volFact; 
  out[4] += (0.5*tmp1[3]*wgt[3]+0.5*tmp1[2]*wgt[2]+0.5*tmp1[1]*wgt[1]+0.5*tmp1[0]*wgt[0])*volFact; 
  out[5] += (0.5*tmp1[2]*wgt[3]+0.5*wgt[2]*tmp1[3]+0.5*tmp1[0]*wgt[1]+0.5*wgt[0]*tmp1[1])*volFact; 
  out[6] += (0.5*tmp1[1]*wgt[3]+0.5*wgt[1]*tmp1[3]+0.5*tmp1[0]*wgt[2]+0.5*wgt[0]*tmp1[2])*volFact; 
  out[7] += (0.5*tmp1[0]*wgt[3]+0.5*wgt[0]*tmp1[3]+0.5*tmp1[1]*wgt[2]+0.5*wgt[1]*tmp1[2])*volFact; 
  out[8] += (0.5*tmp2par[3]*wgt[3]+0.5*tmp2par[2]*wgt[2]+0.5*tmp2par[1]*wgt[1]+0.5*tmp2par[0]*wgt[0])*volFact; 
  out[9] += (0.5*tmp2par[2]*wgt[3]+0.5*wgt[2]*tmp2par[3]+0.5*tmp2par[0]*wgt[1]+0.5*wgt[0]*tmp2par[1])*volFact; 
  out[10] += (0.5*tmp2par[1]*wgt[3]+0.5*wgt[1]*tmp2par[3]+0.5*tmp2par[0]*wgt[2]+0.5*wgt[0]*tmp2par[2])*volFact; 
  out[11] += (0.5*tmp2par[0]*wgt[3]+0.5*wgt[0]*tmp2par[3]+0.5*tmp2par[1]*wgt[2]+0.5*wgt[1]*tmp2par[2])*volFact; 
  out[12] += ((0.5*bmag[0]*tmp2perp[3]*wgt[3]+0.5*tmp2perp[0]*bmag[3]*wgt[3]+0.5*bmag[1]*tmp2perp[2]*wgt[3]+0.5*tmp2perp[1]*bmag[2]*wgt[3]+0.5*wgt[0]*bmag[3]*tmp2perp[3]+0.5*bmag[1]*wgt[2]*tmp2perp[3]+0.5*wgt[1]*bmag[2]*tmp2perp[3]+0.5*tmp2perp[1]*wgt[2]*bmag[3]+0.5*wgt[1]*tmp2perp[2]*bmag[3]+0.5*bmag[0]*tmp2perp[2]*wgt[2]+0.5*tmp2perp[0]*bmag[2]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp2perp[2]+0.5*bmag[0]*tmp2perp[1]*wgt[1]+0.5*tmp2perp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp2perp[1]+0.5*bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_; 
  out[13] += ((0.9*bmag[1]*tmp2perp[3]*wgt[3]+0.9*tmp2perp[1]*bmag[3]*wgt[3]+0.5*bmag[0]*tmp2perp[2]*wgt[3]+0.5*tmp2perp[0]*bmag[2]*wgt[3]+0.9*wgt[1]*bmag[3]*tmp2perp[3]+0.5*bmag[0]*wgt[2]*tmp2perp[3]+0.5*wgt[0]*bmag[2]*tmp2perp[3]+0.5*tmp2perp[0]*wgt[2]*bmag[3]+0.5*wgt[0]*tmp2perp[2]*bmag[3]+0.5*bmag[1]*tmp2perp[2]*wgt[2]+0.5*tmp2perp[1]*bmag[2]*wgt[2]+0.5*wgt[1]*bmag[2]*tmp2perp[2]+0.9*bmag[1]*tmp2perp[1]*wgt[1]+0.5*bmag[0]*tmp2perp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmp2perp[1]+0.5*tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[14] += ((0.9*bmag[2]*tmp2perp[3]*wgt[3]+0.9*tmp2perp[2]*bmag[3]*wgt[3]+0.5*bmag[0]*tmp2perp[1]*wgt[3]+0.5*tmp2perp[0]*bmag[1]*wgt[3]+0.9*wgt[2]*bmag[3]*tmp2perp[3]+0.5*bmag[0]*wgt[1]*tmp2perp[3]+0.5*wgt[0]*bmag[1]*tmp2perp[3]+0.5*tmp2perp[0]*wgt[1]*bmag[3]+0.5*wgt[0]*tmp2perp[1]*bmag[3]+0.9*bmag[2]*tmp2perp[2]*wgt[2]+0.5*bmag[1]*tmp2perp[1]*wgt[2]+0.5*bmag[0]*tmp2perp[0]*wgt[2]+0.5*bmag[1]*wgt[1]*tmp2perp[2]+0.5*bmag[0]*wgt[0]*tmp2perp[2]+0.5*tmp2perp[1]*wgt[1]*bmag[2]+0.5*tmp2perp[0]*wgt[0]*bmag[2])*volFact)/m_; 
  out[15] += ((1.62*bmag[3]*tmp2perp[3]*wgt[3]+0.9*bmag[2]*tmp2perp[2]*wgt[3]+0.9*bmag[1]*tmp2perp[1]*wgt[3]+0.5*bmag[0]*tmp2perp[0]*wgt[3]+0.9*bmag[2]*wgt[2]*tmp2perp[3]+0.9*bmag[1]*wgt[1]*tmp2perp[3]+0.5*bmag[0]*wgt[0]*tmp2perp[3]+0.9*tmp2perp[2]*wgt[2]*bmag[3]+0.9*tmp2perp[1]*wgt[1]*bmag[3]+0.5*tmp2perp[0]*wgt[0]*bmag[3]+0.5*bmag[0]*tmp2perp[1]*wgt[2]+0.5*tmp2perp[0]*bmag[1]*wgt[2]+0.5*bmag[0]*wgt[1]*tmp2perp[2]+0.5*wgt[0]*bmag[1]*tmp2perp[2]+0.5*tmp2perp[0]*wgt[1]*bmag[2]+0.5*wgt[0]*tmp2perp[1]*bmag[2])*volFact)/m_; 
} 

