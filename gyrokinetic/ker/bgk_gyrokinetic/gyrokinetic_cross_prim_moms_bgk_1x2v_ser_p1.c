#include <gkyl_gyrokinetic_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
 
GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x2v_ser_p1(const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross) 
{ 
  // m_:              mass. 
  // prim_moms_:      primitive moments of the distribution function. 
  // prim_moms_cross: cross primitive moments. 
 
  const double m_s = m_self; 
  const double m_r = m_other; 
  const double *n_s = &prim_moms_self[0]; 
  const double *upar_s = &prim_moms_self[2]; 
  const double *vtsq_s = &prim_moms_self[4]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *upar_r = &prim_moms_other[2]; 
  const double *vtsq_r = &prim_moms_other[4]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *upar_sr = &prim_moms_cross[2]; 
  double *vtsq_sr = &prim_moms_cross[4]; 
 
  double msNsNusr[2] = {0.0}; 
  double mrNrNurs[2] = {0.0}; 
  double m_n_nu[2] = {0.0}; 
  double m_n_nu_inv[2] = {0.0}; 
  double alphaE[2] = {0.0}; 

  double msNsNusr_inv[2] = {0.0}; 
  double coeff[2] = {0.0}; 
  double dUpar[2] = {0.0}; 
  double cUpar[2] = {0.0}; 

  double dv; 
  double T1[2] = {0.0}; 
  double T2[2] = {0.0}; 
  double T3[2] = {0.0}; 
  double cVtsq[2] = {0.0}; 
  bool negative_cross_temp = false; 

  binop_mul_1d_ser_p1(n_s, nu_sr, msNsNusr); 
  binop_mul_1d_ser_p1(n_r, nu_rs, mrNrNurs); 
  msNsNusr[0] = m_s * msNsNusr[0]; 
  mrNrNurs[0] = m_r * mrNrNurs[0]; 
  m_n_nu[0] = msNsNusr[0] + mrNrNurs[0]; 
  msNsNusr[1] = m_s * msNsNusr[1]; 
  mrNrNurs[1] = m_r * mrNrNurs[1]; 
  m_n_nu[1] = msNsNusr[1] + mrNrNurs[1]; 
  ser_1x_p1_inv(m_n_nu, m_n_nu_inv); 
  binop_mul_1d_ser_p1(msNsNusr, mrNrNurs, alphaE); 
  binop_mul_1d_ser_p1(alphaE, m_n_nu_inv, alphaE); 
  alphaE[0] = alphaE[0] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[1] = alphaE[1] * 2.0 * betaGreenep1 / (m_s+m_r); 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
 
  ser_1x_p1_inv(msNsNusr, msNsNusr_inv); 
  binop_mul_1d_ser_p1(alphaE, msNsNusr_inv, coeff); 
  dUpar[0] = upar_r[0] - upar_s[0]; 
  dUpar[1] = upar_r[1] - upar_s[1]; 
  binop_mul_1d_ser_p1(coeff, dUpar, cUpar); 
  upar_sr[0] = upar_s[0] + cUpar[0]*(m_s+m_r)/2.0; 
  upar_sr[1] = upar_s[1] + cUpar[1]*(m_s+m_r)/2.0; 
 
  dv = 3.0; 
  T1[0] = dv * (m_r*vtsq_r[0]-m_s*vtsq_s[0]); 
  T1[1] = dv * (m_r*vtsq_r[1]-m_s*vtsq_s[1]); 
  binop_mul_1d_ser_p1(dUpar, dUpar, T2); 
  binop_mul_1d_ser_p1(coeff, T2, T3); 
  cVtsq[0] = T1[0] + m_r*T2[0] - (m_s+m_r)*(m_s+m_r)/4.0*T3[0] ; 
  cVtsq[1] = T1[1] + m_r*T2[1] - (m_s+m_r)*(m_s+m_r)/4.0*T3[1] ; 
  binop_mul_1d_ser_p1(coeff, cVtsq, cVtsq); 
  vtsq_sr[0] = vtsq_s[0] + cVtsq[0]/dv; 
  vtsq_sr[1] = vtsq_s[1] + cVtsq[1]/dv; 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  if (-0.5*(2.449489742783178*vtsq_sr[1]-1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.5*(2.449489742783178*vtsq_sr[1]+1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    upar_sr[0] = upar_s[0]; 
    vtsq_sr[0] = vtsq_s[0]; 
    upar_sr[1] = upar_s[1]; 
    vtsq_sr[1] = vtsq_s[1]; 
  } 
} 
