#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[1]*dxv[1];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[6] = {0.}; 
  hamil[0] = 1.414213562373095*phi[0]*q_+0.7071067811865475*vmapSq[0]*m_; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double alphax[6] = {0.}; 
  alphax[0] = ((0.7071067811865475*cmag[1]*jacobtot_inv[1]*hamil[2]+0.7071067811865475*cmag[0]*jacobtot_inv[0]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[1] = ((0.7071067811865475*cmag[0]*jacobtot_inv[1]*hamil[2]+0.7071067811865475*jacobtot_inv[0]*cmag[1]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[2] = ((1.58113883008419*cmag[1]*jacobtot_inv[1]*hamil[4]+1.58113883008419*cmag[0]*jacobtot_inv[0]*hamil[4])*rdx2)/(vmap[1]*m_); 
  alphax[3] = ((1.58113883008419*cmag[0]*jacobtot_inv[1]*hamil[4]+1.58113883008419*jacobtot_inv[0]*cmag[1]*hamil[4])*rdx2)/(vmap[1]*m_); 

  double alphavpar[6] = {0.}; 
  alphavpar[0] = (((-0.7071067811865475*cmag[1]*hamil[1]*jacobtot_inv[1])-0.7071067811865475*cmag[0]*jacobtot_inv[0]*hamil[1])*rdx2)/(vmap[1]*m_); 
  alphavpar[1] = (((-0.7071067811865475*cmag[0]*hamil[1]*jacobtot_inv[1])-0.7071067811865475*jacobtot_inv[0]*cmag[1]*hamil[1])*rdx2)/(vmap[1]*m_); 

  out[1] += 0.8660254037844386*(alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.7745966692414834*alphax[3]*fin[5]+0.7745966692414833*alphax[2]*fin[4]+0.8660254037844386*(alphax[1]*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.936491673103709*(alphavpar[1]*fin[3]+alphavpar[0]*fin[2]); 
  out[5] += 0.8660254037844386*alphax[1]*fin[5]+0.8660254037844387*alphax[0]*fin[4]+(0.7745966692414834*alphax[3]+1.936491673103709*alphavpar[0])*fin[3]+(0.7745966692414834*alphax[2]+1.936491673103709*alphavpar[1])*fin[2]; 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[3] += -(1.224744871391589*(apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[4] += -(2.738612787525831*(apardot[1]*fin[3]+apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[5] += -(2.738612787525831*(apardot[0]*fin[3]+apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.25*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.1767766952966369*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.25*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.1767766952966369*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
