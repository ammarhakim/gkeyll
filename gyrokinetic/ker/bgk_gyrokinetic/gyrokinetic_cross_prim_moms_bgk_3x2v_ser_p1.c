#include <gkyl_gyrokinetic_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
 
GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_3x2v_ser_p1(const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, const double *nu_sr, const double *nu_rs, double *prim_moms_cross) 
{ 
  // m_:              mass. 
  // prim_moms_:      primitive moments of the distribution function. 
  // prim_moms_cross: cross primitive moments. 
 
  const double m_s = m_self; 
  const double m_r = m_other; 
  const double *n_s = &prim_moms_self[0]; 
  const double *upar_s = &prim_moms_self[8]; 
  const double *vtsq_s = &prim_moms_self[16]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *upar_r = &prim_moms_other[8]; 
  const double *vtsq_r = &prim_moms_other[16]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *upar_sr = &prim_moms_cross[8]; 
  double *vtsq_sr = &prim_moms_cross[16]; 
 
  double msNsNusr[8] = {0.0}; 
  double mrNrNurs[8] = {0.0}; 
  double m_n_nu[8] = {0.0}; 
  double m_n_nu_inv[8] = {0.0}; 
  double alphaE[8] = {0.0}; 

  double msNsNusr_inv[8] = {0.0}; 
  double coeff[8] = {0.0}; 
  double dUpar[8] = {0.0}; 
  double cUpar[8] = {0.0}; 

  double dv; 
  double T1[8] = {0.0}; 
  double T2[8] = {0.0}; 
  double T3[8] = {0.0}; 
  double cVtsq[8] = {0.0}; 
  bool negative_cross_temp = false; 

  binop_mul_3d_ser_p1(n_s, nu_sr, msNsNusr); 
  binop_mul_3d_ser_p1(n_r, nu_rs, mrNrNurs); 
  msNsNusr[0] = m_s * msNsNusr[0]; 
  mrNrNurs[0] = m_r * mrNrNurs[0]; 
  m_n_nu[0] = msNsNusr[0] + mrNrNurs[0]; 
  msNsNusr[1] = m_s * msNsNusr[1]; 
  mrNrNurs[1] = m_r * mrNrNurs[1]; 
  m_n_nu[1] = msNsNusr[1] + mrNrNurs[1]; 
  msNsNusr[2] = m_s * msNsNusr[2]; 
  mrNrNurs[2] = m_r * mrNrNurs[2]; 
  m_n_nu[2] = msNsNusr[2] + mrNrNurs[2]; 
  msNsNusr[3] = m_s * msNsNusr[3]; 
  mrNrNurs[3] = m_r * mrNrNurs[3]; 
  m_n_nu[3] = msNsNusr[3] + mrNrNurs[3]; 
  msNsNusr[4] = m_s * msNsNusr[4]; 
  mrNrNurs[4] = m_r * mrNrNurs[4]; 
  m_n_nu[4] = msNsNusr[4] + mrNrNurs[4]; 
  msNsNusr[5] = m_s * msNsNusr[5]; 
  mrNrNurs[5] = m_r * mrNrNurs[5]; 
  m_n_nu[5] = msNsNusr[5] + mrNrNurs[5]; 
  msNsNusr[6] = m_s * msNsNusr[6]; 
  mrNrNurs[6] = m_r * mrNrNurs[6]; 
  m_n_nu[6] = msNsNusr[6] + mrNrNurs[6]; 
  msNsNusr[7] = m_s * msNsNusr[7]; 
  mrNrNurs[7] = m_r * mrNrNurs[7]; 
  m_n_nu[7] = msNsNusr[7] + mrNrNurs[7]; 
  ser_3x_p1_inv(m_n_nu, m_n_nu_inv); 
  binop_mul_3d_ser_p1(msNsNusr, mrNrNurs, alphaE); 
  binop_mul_3d_ser_p1(alphaE, m_n_nu_inv, alphaE); 
  alphaE[0] = alphaE[0] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[1] = alphaE[1] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[2] = alphaE[2] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[3] = alphaE[3] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[4] = alphaE[4] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[5] = alphaE[5] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[6] = alphaE[6] * 2.0 * betaGreenep1 / (m_s+m_r); 
  alphaE[7] = alphaE[7] * 2.0 * betaGreenep1 / (m_s+m_r); 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
  n_sr[2] = n_s[2]; 
  n_sr[3] = n_s[3]; 
  n_sr[4] = n_s[4]; 
  n_sr[5] = n_s[5]; 
  n_sr[6] = n_s[6]; 
  n_sr[7] = n_s[7]; 
 
  ser_3x_p1_inv(msNsNusr, msNsNusr_inv); 
  binop_mul_3d_ser_p1(alphaE, msNsNusr_inv, coeff); 
  dUpar[0] = upar_r[0] - upar_s[0]; 
  dUpar[1] = upar_r[1] - upar_s[1]; 
  dUpar[2] = upar_r[2] - upar_s[2]; 
  dUpar[3] = upar_r[3] - upar_s[3]; 
  dUpar[4] = upar_r[4] - upar_s[4]; 
  dUpar[5] = upar_r[5] - upar_s[5]; 
  dUpar[6] = upar_r[6] - upar_s[6]; 
  dUpar[7] = upar_r[7] - upar_s[7]; 
  binop_mul_3d_ser_p1(coeff, dUpar, cUpar); 
  upar_sr[0] = upar_s[0] + cUpar[0]*(m_s+m_r)/2.0; 
  upar_sr[1] = upar_s[1] + cUpar[1]*(m_s+m_r)/2.0; 
  upar_sr[2] = upar_s[2] + cUpar[2]*(m_s+m_r)/2.0; 
  upar_sr[3] = upar_s[3] + cUpar[3]*(m_s+m_r)/2.0; 
  upar_sr[4] = upar_s[4] + cUpar[4]*(m_s+m_r)/2.0; 
  upar_sr[5] = upar_s[5] + cUpar[5]*(m_s+m_r)/2.0; 
  upar_sr[6] = upar_s[6] + cUpar[6]*(m_s+m_r)/2.0; 
  upar_sr[7] = upar_s[7] + cUpar[7]*(m_s+m_r)/2.0; 
 
  dv = 3.0; 
  T1[0] = dv * (m_r*vtsq_r[0]-m_s*vtsq_s[0]); 
  T1[1] = dv * (m_r*vtsq_r[1]-m_s*vtsq_s[1]); 
  T1[2] = dv * (m_r*vtsq_r[2]-m_s*vtsq_s[2]); 
  T1[3] = dv * (m_r*vtsq_r[3]-m_s*vtsq_s[3]); 
  T1[4] = dv * (m_r*vtsq_r[4]-m_s*vtsq_s[4]); 
  T1[5] = dv * (m_r*vtsq_r[5]-m_s*vtsq_s[5]); 
  T1[6] = dv * (m_r*vtsq_r[6]-m_s*vtsq_s[6]); 
  T1[7] = dv * (m_r*vtsq_r[7]-m_s*vtsq_s[7]); 
  binop_mul_3d_ser_p1(dUpar, dUpar, T2); 
  binop_mul_3d_ser_p1(coeff, T2, T3); 
  cVtsq[0] = T1[0] + m_r*T2[0] - (m_s+m_r)*(m_s+m_r)/4.0*T3[0] ; 
  cVtsq[1] = T1[1] + m_r*T2[1] - (m_s+m_r)*(m_s+m_r)/4.0*T3[1] ; 
  cVtsq[2] = T1[2] + m_r*T2[2] - (m_s+m_r)*(m_s+m_r)/4.0*T3[2] ; 
  cVtsq[3] = T1[3] + m_r*T2[3] - (m_s+m_r)*(m_s+m_r)/4.0*T3[3] ; 
  cVtsq[4] = T1[4] + m_r*T2[4] - (m_s+m_r)*(m_s+m_r)/4.0*T3[4] ; 
  cVtsq[5] = T1[5] + m_r*T2[5] - (m_s+m_r)*(m_s+m_r)/4.0*T3[5] ; 
  cVtsq[6] = T1[6] + m_r*T2[6] - (m_s+m_r)*(m_s+m_r)/4.0*T3[6] ; 
  cVtsq[7] = T1[7] + m_r*T2[7] - (m_s+m_r)*(m_s+m_r)/4.0*T3[7] ; 
  binop_mul_3d_ser_p1(coeff, cVtsq, cVtsq); 
  vtsq_sr[0] = vtsq_s[0] + cVtsq[0]/dv; 
  vtsq_sr[1] = vtsq_s[1] + cVtsq[1]/dv; 
  vtsq_sr[2] = vtsq_s[2] + cVtsq[2]/dv; 
  vtsq_sr[3] = vtsq_s[3] + cVtsq[3]/dv; 
  vtsq_sr[4] = vtsq_s[4] + cVtsq[4]/dv; 
  vtsq_sr[5] = vtsq_s[5] + cVtsq[5]/dv; 
  vtsq_sr[6] = vtsq_s[6] + cVtsq[6]/dv; 
  vtsq_sr[7] = vtsq_s[7] + cVtsq[7]/dv; 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  if (-0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*(vtsq_sr[6]+vtsq_sr[5]+vtsq_sr[4])+2.449489742783178*(vtsq_sr[3]+vtsq_sr[2]+vtsq_sr[1])-1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*vtsq_sr[6]-4.242640687119286*(vtsq_sr[5]+vtsq_sr[4])-2.449489742783178*(vtsq_sr[3]+vtsq_sr[2])+2.449489742783178*vtsq_sr[1]+1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*vtsq_sr[6]+4.242640687119286*vtsq_sr[5]-4.242640687119286*vtsq_sr[4]-2.449489742783178*vtsq_sr[3]+2.449489742783178*vtsq_sr[2]-2.449489742783178*vtsq_sr[1]+1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*(vtsq_sr[6]+vtsq_sr[5])-4.242640687119286*vtsq_sr[4]+2.449489742783178*vtsq_sr[3]-2.449489742783178*(vtsq_sr[2]+vtsq_sr[1])-1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*(vtsq_sr[6]+vtsq_sr[5])+4.242640687119286*vtsq_sr[4]+2.449489742783178*vtsq_sr[3]-2.449489742783178*(vtsq_sr[2]+vtsq_sr[1])+1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*vtsq_sr[6]-4.242640687119286*vtsq_sr[5]+4.242640687119286*vtsq_sr[4]-2.449489742783178*vtsq_sr[3]+2.449489742783178*vtsq_sr[2]-2.449489742783178*vtsq_sr[1]-1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*vtsq_sr[6]+4.242640687119286*(vtsq_sr[5]+vtsq_sr[4])-2.449489742783178*(vtsq_sr[3]+vtsq_sr[2])+2.449489742783178*vtsq_sr[1]-1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*(vtsq_sr[6]+vtsq_sr[5]+vtsq_sr[4])+2.449489742783178*(vtsq_sr[3]+vtsq_sr[2]+vtsq_sr[1])+1.414213562373095*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    upar_sr[0] = upar_s[0]; 
    vtsq_sr[0] = vtsq_s[0]; 
    upar_sr[1] = upar_s[1]; 
    vtsq_sr[1] = vtsq_s[1]; 
    upar_sr[2] = upar_s[2]; 
    vtsq_sr[2] = vtsq_s[2]; 
    upar_sr[3] = upar_s[3]; 
    vtsq_sr[3] = vtsq_s[3]; 
    upar_sr[4] = upar_s[4]; 
    vtsq_sr[4] = vtsq_s[4]; 
    upar_sr[5] = upar_s[5]; 
    vtsq_sr[5] = vtsq_s[5]; 
    upar_sr[6] = upar_s[6]; 
    vtsq_sr[6] = vtsq_s[6]; 
    upar_sr[7] = upar_s[7]; 
    vtsq_sr[7] = vtsq_s[7]; 
  } 
} 
