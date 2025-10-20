#include <gkyl_gyrokinetic_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
 
GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_1x1v_ser_p2(const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross) 
{ 
  // m_:              mass. 
  // prim_moms_:      primitive moments of the distribution function. 
  // prim_moms_cross: cross primitive moments. 
 
  const double m_s = m_self; 
  const double m_r = m_other; 
  const double *n_s = &prim_moms_self[0]; 
  const double *upar_s = &prim_moms_self[3]; 
  const double *vtsq_s = &prim_moms_self[6]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *upar_r = &prim_moms_other[3]; 
  const double *vtsq_r = &prim_moms_other[6]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *upar_sr = &prim_moms_cross[3]; 
  double *vtsq_sr = &prim_moms_cross[6]; 
 
  double msNsNusr[3] = {0.0}; 
  double mrNrNurs[3] = {0.0}; 
  double m_n_nu[3] = {0.0}; 
  double m_n_nu_inv[3] = {0.0}; 
  double alphaE[3] = {0.0}; 

  double msNsNusr_inv[3] = {0.0}; 
  double coeff[3] = {0.0}; 
  double dUpar[3] = {0.0}; 
  double cUpar[3] = {0.0}; 

  double dv; 
  double T1[3] = {0.0}; 
  double T2[3] = {0.0}; 
  double T3[3] = {0.0}; 
  double cVtsq[3] = {0.0}; 
  bool negative_cross_temp = false; 

  binop_mul_1d_ser_p2(n_s, nu_sr, msNsNusr); 
  binop_mul_1d_ser_p2(n_r, nu_rs, mrNrNurs); 
  msNsNusr[0] = m_s * msNsNusr[0]; 
  mrNrNurs[0] = m_r * mrNrNurs[0]; 
  m_n_nu[0] = msNsNusr[0] + mrNrNurs[0]; 
  msNsNusr[1] = m_s * msNsNusr[1]; 
  mrNrNurs[1] = m_r * mrNrNurs[1]; 
  m_n_nu[1] = msNsNusr[1] + mrNrNurs[1]; 
  msNsNusr[2] = m_s * msNsNusr[2]; 
  mrNrNurs[2] = m_r * mrNrNurs[2]; 
  m_n_nu[2] = msNsNusr[2] + mrNrNurs[2]; 
  ser_1x_p2_inv(m_n_nu, m_n_nu_inv); 
  binop_mul_1d_ser_p2(msNsNusr, mrNrNurs, alphaE); 
  binop_mul_1d_ser_p2(alphaE, m_n_nu_inv, alphaE); 
  alphaE[0] = alphaE[0] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[1] = alphaE[1] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[2] = alphaE[2] * 2.0 * betaGreenep1 / (m_s+m_r); 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
  n_sr[2] = n_s[2]; 
 
  ser_1x_p2_inv(msNsNusr, msNsNusr_inv); 
  binop_mul_1d_ser_p2(alphaE, msNsNusr_inv, coeff); 
  dUpar[0] = upar_r[0] - upar_s[0]; 
  dUpar[1] = upar_r[1] - upar_s[1]; 
  dUpar[2] = upar_r[2] - upar_s[2]; 
  binop_mul_1d_ser_p2(coeff, dUpar, cUpar); 
  upar_sr[0] = upar_s[0] + cUpar[0]*(m_s+m_r)/2.0; 
  upar_sr[1] = upar_s[1] + cUpar[1]*(m_s+m_r)/2.0; 
  upar_sr[2] = upar_s[2] + cUpar[2]*(m_s+m_r)/2.0; 
 
  dv = 1.0; 
  T1[0] = dv * (m_r*vtsq_r[0]-m_s*vtsq_s[0]); 
  T1[1] = dv * (m_r*vtsq_r[1]-m_s*vtsq_s[1]); 
  T1[2] = dv * (m_r*vtsq_r[2]-m_s*vtsq_s[2]); 
  binop_mul_1d_ser_p2(dUpar, dUpar, T2); 
  binop_mul_1d_ser_p2(coeff, T2, T3); 
  cVtsq[0] = T1[0] + m_r*T2[0] - (m_s+m_r)*(m_s+m_r)/4.0*T3[0] ; 
  cVtsq[1] = T1[1] + m_r*T2[1] - (m_s+m_r)*(m_s+m_r)/4.0*T3[1] ; 
  cVtsq[2] = T1[2] + m_r*T2[2] - (m_s+m_r)*(m_s+m_r)/4.0*T3[2] ; 
  binop_mul_1d_ser_p2(coeff, cVtsq, cVtsq); 
  vtsq_sr[0] = vtsq_s[0] + cVtsq[0]/dv; 
  vtsq_sr[1] = vtsq_s[1] + cVtsq[1]/dv; 
  vtsq_sr[2] = vtsq_s[2] + cVtsq[2]/dv; 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  if (0.7071067811865475*(2.23606797749979*vtsq_sr[2]-1.732050807568877*vtsq_sr[1]+vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.7071067811865475*(2.23606797749979*vtsq_sr[2]+1.732050807568877*vtsq_sr[1]+vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    upar_sr[0] = upar_s[0]; 
    vtsq_sr[0] = vtsq_s[0]; 
    upar_sr[1] = upar_s[1]; 
    vtsq_sr[1] = vtsq_s[1]; 
    upar_sr[2] = upar_s[2]; 
    vtsq_sr[2] = vtsq_s[2]; 
  } 
} 
