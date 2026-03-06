#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_weighted_gyrokinetic_M0_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  out[0] += (1.4142135623730951*wgt[2]*f[7]+1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.264911064067352*wgt[1]*f[7]+1.264911064067352*f[1]*wgt[2]+1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
  out[2] += (0.9035079029052515*wgt[2]*f[7]+1.4142135623730951*wgt[0]*f[7]+1.4142135623730951*f[0]*wgt[2]+1.264911064067352*f[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M1_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[3]; 
  tmp[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  tmp[2] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[7]; 
  out[0] += (0.7071067811865475*tmp[2]*wgt[2]+0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmp[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp[2]+0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.45175395145262565*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmpPar[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmpPar[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmpPar[1] = 0.8944271909999161*vmap1R2*f[12]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmpPar[2] = 2.0000000000000004*vmap[0]*vmap[1]*f[11]+vmap1R2*f[7]+vmap0R2*f[7]; 

  out[0] += (0.7071067811865475*tmpPar[2]*wgt[2]+0.7071067811865475*tmpPar[1]*wgt[1]+0.7071067811865475*tmpPar[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmpPar[1]*wgt[2]+0.6324555320336759*wgt[1]*tmpPar[2]+0.7071067811865475*tmpPar[0]*wgt[1]+0.7071067811865475*wgt[0]*tmpPar[1])*volFact; 
  out[2] += (0.45175395145262565*tmpPar[2]*wgt[2]+0.7071067811865475*tmpPar[0]*wgt[2]+0.7071067811865475*wgt[0]*tmpPar[2]+0.6324555320336759*tmpPar[1]*wgt[1])*volFact; 
  double tmpPerp[3]; 
  tmpPerp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmpPerp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmpPerp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 
  out[0] += (2.0*(0.31943828249996997*bmag[2]*tmpPerp[2]*wgt[2]+0.5*bmag[0]*tmpPerp[2]*wgt[2]+0.5*tmpPerp[0]*bmag[2]*wgt[2]+0.4472135954999579*bmag[1]*tmpPerp[1]*wgt[2]+0.5*wgt[0]*bmag[2]*tmpPerp[2]+0.4472135954999579*bmag[1]*wgt[1]*tmpPerp[2]+0.4472135954999579*tmpPerp[1]*wgt[1]*bmag[2]+0.5*bmag[0]*tmpPerp[1]*wgt[1]+0.5*tmpPerp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmpPerp[1]+0.5*bmag[0]*tmpPerp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.7857142857142857*bmag[1]*tmpPerp[2]*wgt[2]+0.7857142857142857*tmpPerp[1]*bmag[2]*wgt[2]+0.4472135954999579*bmag[0]*tmpPerp[1]*wgt[2]+0.4472135954999579*tmpPerp[0]*bmag[1]*wgt[2]+0.7857142857142857*wgt[1]*bmag[2]*tmpPerp[2]+0.4472135954999579*bmag[0]*wgt[1]*tmpPerp[2]+0.4472135954999579*wgt[0]*bmag[1]*tmpPerp[2]+0.4472135954999579*tmpPerp[0]*wgt[1]*bmag[2]+0.4472135954999579*wgt[0]*tmpPerp[1]*bmag[2]+0.9*bmag[1]*tmpPerp[1]*wgt[1]+0.5*bmag[0]*tmpPerp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmpPerp[1]+0.5*tmpPerp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(1.0714285714285714*bmag[2]*tmpPerp[2]*wgt[2]+0.31943828249996997*bmag[0]*tmpPerp[2]*wgt[2]+0.31943828249996997*tmpPerp[0]*bmag[2]*wgt[2]+0.7857142857142857*bmag[1]*tmpPerp[1]*wgt[2]+0.5*bmag[0]*tmpPerp[0]*wgt[2]+0.31943828249996997*wgt[0]*bmag[2]*tmpPerp[2]+0.7857142857142857*bmag[1]*wgt[1]*tmpPerp[2]+0.5*bmag[0]*wgt[0]*tmpPerp[2]+0.7857142857142857*tmpPerp[1]*wgt[1]*bmag[2]+0.5*tmpPerp[0]*wgt[0]*bmag[2]+0.4472135954999579*bmag[0]*tmpPerp[1]*wgt[1]+0.4472135954999579*tmpPerp[0]*bmag[1]*wgt[1]+0.4472135954999579*wgt[0]*bmag[1]*tmpPerp[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PAR_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp[1] = 0.8944271909999161*vmap1R2*f[12]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp[2] = 2.0000000000000004*vmap[0]*vmap[1]*f[11]+vmap1R2*f[7]+vmap0R2*f[7]; 

  out[0] += (0.7071067811865475*tmp[2]*wgt[2]+0.7071067811865475*tmp[1]*wgt[1]+0.7071067811865475*tmp[0]*wgt[0])*volFact; 
  out[1] += (0.6324555320336759*tmp[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp[2]+0.7071067811865475*tmp[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp[1])*volFact; 
  out[2] += (0.45175395145262565*tmp[2]*wgt[2]+0.7071067811865475*tmp[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp[2]+0.6324555320336759*tmp[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[3]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  tmp1[2] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[7]; 
  double tmp2par[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[12]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp2par[2] = 2.0000000000000004*vmap[0]*vmap[1]*f[11]+vmap1R2*f[7]+vmap0R2*f[7]; 
  double tmp2perp[3]; 
  tmp2perp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2perp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 

  out[0] += (1.4142135623730951*wgt[2]*f[7]+1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.264911064067352*wgt[1]*f[7]+1.264911064067352*f[1]*wgt[2]+1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
  out[2] += (0.9035079029052515*wgt[2]*f[7]+1.4142135623730951*wgt[0]*f[7]+1.4142135623730951*f[0]*wgt[2]+1.264911064067352*f[1]*wgt[1])*volFact; 
  out[3] += (0.7071067811865475*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[4] += (0.6324555320336759*tmp1[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp1[2]+0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[5] += (0.45175395145262565*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp1[2]+0.6324555320336759*tmp1[1]*wgt[1])*volFact; 
  out[6] += ((0.6388765649999399*bmag[2]*tmp2perp[2]*wgt[2]+bmag[0]*tmp2perp[2]*wgt[2]+tmp2perp[0]*bmag[2]*wgt[2]+0.8944271909999159*bmag[1]*tmp2perp[1]*wgt[2]+wgt[0]*bmag[2]*tmp2perp[2]+0.8944271909999159*bmag[1]*wgt[1]*tmp2perp[2]+0.8944271909999159*tmp2perp[1]*wgt[1]*bmag[2]+bmag[0]*tmp2perp[1]*wgt[1]+tmp2perp[0]*bmag[1]*wgt[1]+wgt[0]*bmag[1]*tmp2perp[1]+bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_+(0.7071067811865475*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[7] += ((1.5714285714285714*bmag[1]*tmp2perp[2]*wgt[2]+1.5714285714285714*tmp2perp[1]*bmag[2]*wgt[2]+0.8944271909999159*bmag[0]*tmp2perp[1]*wgt[2]+0.8944271909999159*tmp2perp[0]*bmag[1]*wgt[2]+1.5714285714285714*wgt[1]*bmag[2]*tmp2perp[2]+0.8944271909999159*bmag[0]*wgt[1]*tmp2perp[2]+0.8944271909999159*wgt[0]*bmag[1]*tmp2perp[2]+0.8944271909999159*tmp2perp[0]*wgt[1]*bmag[2]+0.8944271909999159*wgt[0]*tmp2perp[1]*bmag[2]+1.8*bmag[1]*tmp2perp[1]*wgt[1]+bmag[0]*tmp2perp[0]*wgt[1]+bmag[0]*wgt[0]*tmp2perp[1]+tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_+(0.6324555320336759*tmp2par[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp2par[2]+0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
  out[8] += ((2.142857142857143*bmag[2]*tmp2perp[2]*wgt[2]+0.6388765649999399*bmag[0]*tmp2perp[2]*wgt[2]+0.6388765649999399*tmp2perp[0]*bmag[2]*wgt[2]+1.5714285714285714*bmag[1]*tmp2perp[1]*wgt[2]+bmag[0]*tmp2perp[0]*wgt[2]+0.6388765649999399*wgt[0]*bmag[2]*tmp2perp[2]+1.5714285714285714*bmag[1]*wgt[1]*tmp2perp[2]+bmag[0]*wgt[0]*tmp2perp[2]+1.5714285714285714*tmp2perp[1]*wgt[1]*bmag[2]+tmp2perp[0]*wgt[0]*bmag[2]+0.8944271909999159*bmag[0]*tmp2perp[1]*wgt[1]+0.8944271909999159*tmp2perp[0]*bmag[1]*wgt[1]+0.8944271909999159*wgt[0]*bmag[1]*tmp2perp[1])*volFact)/m_+(0.45175395145262565*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp2par[2]+0.6324555320336759*tmp2par[1]*wgt[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[3]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1])*m_*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[7])*m_*volFact; 
  out[6] += (1.4142135623730951*phi[2]*f[7]+1.4142135623730951*f[1]*phi[1]+1.4142135623730951*f[0]*phi[0])*q_*volFact+(0.4472135954999579*vmap1R2*f[8]+vmap[0]*vmap[1]*f[2]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact+(0.7071067811865475*bmag[2]*tmp[2]+0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact; 
  out[7] += (1.264911064067352*phi[1]*f[7]+1.264911064067352*f[1]*phi[2]+1.4142135623730951*f[0]*phi[1]+1.4142135623730951*phi[0]*f[1])*q_*volFact+(0.44721359549995804*vmap1R2*f[12]+vmap[0]*vmap[1]*f[4]+0.5*f[1]*vmap1R2+0.5*vmap0R2*f[1])*m_*volFact+(0.6324555320336759*bmag[1]*tmp[2]+0.6324555320336759*tmp[1]*bmag[2]+0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact; 
  out[8] += (0.9035079029052515*phi[2]*f[7]+1.4142135623730951*phi[0]*f[7]+1.4142135623730951*f[0]*phi[2]+1.264911064067352*f[1]*phi[1])*q_*volFact+(1.0000000000000002*vmap[0]*vmap[1]*f[11]+0.5*vmap1R2*f[7]+0.5*vmap0R2*f[7])*m_*volFact+(0.45175395145262565*bmag[2]*tmp[2]+0.7071067811865475*bmag[0]*tmp[2]+0.7071067811865475*tmp[0]*bmag[2]+0.6324555320336759*bmag[1]*tmp[1])*volFact; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M2PERP_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp[3]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 
  out[0] += (2.0*(0.31943828249996997*bmag[2]*tmp[2]*wgt[2]+0.5*bmag[0]*tmp[2]*wgt[2]+0.5*tmp[0]*bmag[2]*wgt[2]+0.4472135954999579*bmag[1]*tmp[1]*wgt[2]+0.5*wgt[0]*bmag[2]*tmp[2]+0.4472135954999579*bmag[1]*wgt[1]*tmp[2]+0.4472135954999579*tmp[1]*wgt[1]*bmag[2]+0.5*bmag[0]*tmp[1]*wgt[1]+0.5*tmp[0]*bmag[1]*wgt[1]+0.5*wgt[0]*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0]*wgt[0])*volFact)/m_; 
  out[1] += (2.0*(0.7857142857142857*bmag[1]*tmp[2]*wgt[2]+0.7857142857142857*tmp[1]*bmag[2]*wgt[2]+0.4472135954999579*bmag[0]*tmp[1]*wgt[2]+0.4472135954999579*tmp[0]*bmag[1]*wgt[2]+0.7857142857142857*wgt[1]*bmag[2]*tmp[2]+0.4472135954999579*bmag[0]*wgt[1]*tmp[2]+0.4472135954999579*wgt[0]*bmag[1]*tmp[2]+0.4472135954999579*tmp[0]*wgt[1]*bmag[2]+0.4472135954999579*wgt[0]*tmp[1]*bmag[2]+0.9*bmag[1]*tmp[1]*wgt[1]+0.5*bmag[0]*tmp[0]*wgt[1]+0.5*bmag[0]*wgt[0]*tmp[1]+0.5*tmp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(1.0714285714285714*bmag[2]*tmp[2]*wgt[2]+0.31943828249996997*bmag[0]*tmp[2]*wgt[2]+0.31943828249996997*tmp[0]*bmag[2]*wgt[2]+0.7857142857142857*bmag[1]*tmp[1]*wgt[2]+0.5*bmag[0]*tmp[0]*wgt[2]+0.31943828249996997*wgt[0]*bmag[2]*tmp[2]+0.7857142857142857*bmag[1]*wgt[1]*tmp[2]+0.5*bmag[0]*wgt[0]*tmp[2]+0.7857142857142857*tmp[1]*wgt[1]*bmag[2]+0.5*tmp[0]*wgt[0]*bmag[2]+0.4472135954999579*bmag[0]*tmp[1]*wgt[1]+0.4472135954999579*tmp[0]*bmag[1]*wgt[1]+0.4472135954999579*wgt[0]*bmag[1]*tmp[1])*volFact)/m_; 
} 

GKYL_CU_DH void mom_weighted_gyrokinetic_M0M1M2PARM2PERP_WEIGHTX_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[1]*dxv[2]/m_; 
  double tmp1[3]; 
  tmp1[0] = 1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0]; 
  tmp1[1] = 1.4142135623730951*vmap[1]*f[4]+1.4142135623730951*vmap[0]*f[1]; 
  tmp1[2] = 1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[7]; 
  double tmp2par[3]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  tmp2par[0] = 0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2; 
  tmp2par[1] = 0.8944271909999161*vmap1R2*f[12]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1]; 
  tmp2par[2] = 2.0000000000000004*vmap[0]*vmap[1]*f[11]+vmap1R2*f[7]+vmap0R2*f[7]; 
  double tmp2perp[3]; 
  tmp2perp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp2perp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp2perp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 

  out[0] += (1.4142135623730951*wgt[2]*f[7]+1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (1.264911064067352*wgt[1]*f[7]+1.264911064067352*f[1]*wgt[2]+1.4142135623730951*f[0]*wgt[1]+1.4142135623730951*wgt[0]*f[1])*volFact; 
  out[2] += (0.9035079029052515*wgt[2]*f[7]+1.4142135623730951*wgt[0]*f[7]+1.4142135623730951*f[0]*wgt[2]+1.264911064067352*f[1]*wgt[1])*volFact; 
  out[3] += (0.7071067811865475*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[1]*wgt[1]+0.7071067811865475*tmp1[0]*wgt[0])*volFact; 
  out[4] += (0.6324555320336759*tmp1[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp1[2]+0.7071067811865475*tmp1[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp1[1])*volFact; 
  out[5] += (0.45175395145262565*tmp1[2]*wgt[2]+0.7071067811865475*tmp1[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp1[2]+0.6324555320336759*tmp1[1]*wgt[1])*volFact; 
  out[6] += (0.7071067811865475*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[1]*wgt[1]+0.7071067811865475*tmp2par[0]*wgt[0])*volFact; 
  out[7] += (0.6324555320336759*tmp2par[1]*wgt[2]+0.6324555320336759*wgt[1]*tmp2par[2]+0.7071067811865475*tmp2par[0]*wgt[1]+0.7071067811865475*wgt[0]*tmp2par[1])*volFact; 
  out[8] += (0.45175395145262565*tmp2par[2]*wgt[2]+0.7071067811865475*tmp2par[0]*wgt[2]+0.7071067811865475*wgt[0]*tmp2par[2]+0.6324555320336759*tmp2par[1]*wgt[1])*volFact; 
  out[9] += ((0.6388765649999399*bmag[2]*tmp2perp[2]*wgt[2]+bmag[0]*tmp2perp[2]*wgt[2]+tmp2perp[0]*bmag[2]*wgt[2]+0.8944271909999159*bmag[1]*tmp2perp[1]*wgt[2]+wgt[0]*bmag[2]*tmp2perp[2]+0.8944271909999159*bmag[1]*wgt[1]*tmp2perp[2]+0.8944271909999159*tmp2perp[1]*wgt[1]*bmag[2]+bmag[0]*tmp2perp[1]*wgt[1]+tmp2perp[0]*bmag[1]*wgt[1]+wgt[0]*bmag[1]*tmp2perp[1]+bmag[0]*tmp2perp[0]*wgt[0])*volFact)/m_; 
  out[10] += ((1.5714285714285714*bmag[1]*tmp2perp[2]*wgt[2]+1.5714285714285714*tmp2perp[1]*bmag[2]*wgt[2]+0.8944271909999159*bmag[0]*tmp2perp[1]*wgt[2]+0.8944271909999159*tmp2perp[0]*bmag[1]*wgt[2]+1.5714285714285714*wgt[1]*bmag[2]*tmp2perp[2]+0.8944271909999159*bmag[0]*wgt[1]*tmp2perp[2]+0.8944271909999159*wgt[0]*bmag[1]*tmp2perp[2]+0.8944271909999159*tmp2perp[0]*wgt[1]*bmag[2]+0.8944271909999159*wgt[0]*tmp2perp[1]*bmag[2]+1.8*bmag[1]*tmp2perp[1]*wgt[1]+bmag[0]*tmp2perp[0]*wgt[1]+bmag[0]*wgt[0]*tmp2perp[1]+tmp2perp[0]*wgt[0]*bmag[1])*volFact)/m_; 
  out[11] += ((2.142857142857143*bmag[2]*tmp2perp[2]*wgt[2]+0.6388765649999399*bmag[0]*tmp2perp[2]*wgt[2]+0.6388765649999399*tmp2perp[0]*bmag[2]*wgt[2]+1.5714285714285714*bmag[1]*tmp2perp[1]*wgt[2]+bmag[0]*tmp2perp[0]*wgt[2]+0.6388765649999399*wgt[0]*bmag[2]*tmp2perp[2]+1.5714285714285714*bmag[1]*wgt[1]*tmp2perp[2]+bmag[0]*wgt[0]*tmp2perp[2]+1.5714285714285714*tmp2perp[1]*wgt[1]*bmag[2]+tmp2perp[0]*wgt[0]*bmag[2]+0.8944271909999159*bmag[0]*tmp2perp[1]*wgt[1]+0.8944271909999159*tmp2perp[0]*bmag[1]*wgt[1]+0.8944271909999159*wgt[0]*bmag[1]*tmp2perp[1])*volFact)/m_; 
} 

