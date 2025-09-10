#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdmu2 = 2.0/dxv[2];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[1]*dxv[1];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[12] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+vmapSq[0]*m_+1.414213562373095*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double alphax[12] = {0.}; 
  alphax[0] = ((0.7071067811865475*cmag[1]*jacobtot_inv[1]*hamil[2]+0.7071067811865475*cmag[0]*jacobtot_inv[0]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[1] = ((0.7071067811865475*cmag[0]*jacobtot_inv[1]*hamil[2]+0.7071067811865475*jacobtot_inv[0]*cmag[1]*hamil[2])*rdx2)/(vmap[1]*m_); 
  alphax[2] = ((1.58113883008419*cmag[1]*jacobtot_inv[1]*hamil[8]+1.58113883008419*cmag[0]*jacobtot_inv[0]*hamil[8])*rdx2)/(vmap[1]*m_); 
  alphax[4] = ((1.58113883008419*cmag[0]*jacobtot_inv[1]*hamil[8]+1.58113883008419*jacobtot_inv[0]*cmag[1]*hamil[8])*rdx2)/(vmap[1]*m_); 

  double alphavpar[12] = {0.}; 
  alphavpar[0] = (((-0.7071067811865475*cmag[1]*hamil[1]*jacobtot_inv[1])-0.7071067811865475*cmag[0]*jacobtot_inv[0]*hamil[1])*rdx2)/(vmap[1]*m_); 
  alphavpar[1] = (((-0.7071067811865475*cmag[0]*hamil[1]*jacobtot_inv[1])-0.7071067811865475*jacobtot_inv[0]*cmag[1]*hamil[1])*rdx2)/(vmap[1]*m_); 
  alphavpar[3] = (((-0.7071067811865475*cmag[1]*jacobtot_inv[1]*hamil[5])-0.7071067811865475*cmag[0]*jacobtot_inv[0]*hamil[5])*rdx2)/(vmap[1]*m_); 
  alphavpar[5] = (((-0.7071067811865475*cmag[0]*jacobtot_inv[1]*hamil[5])-0.7071067811865475*jacobtot_inv[0]*cmag[1]*hamil[5])*rdx2)/(vmap[1]*m_); 

  out[1] += 0.6123724356957944*(alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.6123724356957944*(alphavpar[5]*fin[5]+alphavpar[3]*fin[3]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.5477225575051661*(alphax[4]*fin[9]+alphax[2]*fin[8])+0.6123724356957944*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphax[1]*fin[4]+fin[1]*alphax[4]+alphax[0]*fin[2]+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[5] += 0.6123724356957944*(alphax[4]*fin[7]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[6] += 0.6123724356957944*(alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.5477225575051661*(alphax[4]*fin[11]+alphax[2]*fin[10])+0.6123724356957944*(alphax[1]*fin[7]+alphax[0]*fin[6]+(alphax[4]+alphavpar[0])*fin[5]+fin[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*fin[3]+fin[1]*alphavpar[3]); 
  out[8] += 1.369306393762915*(alphavpar[5]*fin[7]+alphavpar[3]*fin[6]+alphavpar[1]*fin[4]+alphavpar[0]*fin[2]); 
  out[9] += 0.6123724356957944*(alphax[1]*fin[9]+alphax[0]*fin[8])+1.369306393762915*(alphavpar[3]*fin[7]+alphavpar[5]*fin[6])+(0.5477225575051661*alphax[4]+1.369306393762915*alphavpar[0])*fin[4]+(0.5477225575051661*alphax[2]+1.369306393762915*alphavpar[1])*fin[2]; 
  out[10] += 1.369306393762915*(alphavpar[1]*fin[7]+alphavpar[0]*fin[6]+fin[4]*alphavpar[5]+fin[2]*alphavpar[3]); 
  out[11] += 0.6123724356957944*(alphax[1]*fin[11]+alphax[0]*fin[10])+(0.5477225575051661*alphax[4]+1.369306393762915*alphavpar[0])*fin[7]+0.5477225575051661*alphax[2]*fin[6]+1.369306393762915*(alphavpar[1]*fin[6]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]); 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[4] += -(1.224744871391589*(apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[6] += -(1.224744871391589*(apardot[1]*fin[5]+apardot[0]*fin[3])*q_*rdvpar2)/m_; 
  out[7] += -(1.224744871391589*(apardot[0]*fin[5]+apardot[1]*fin[3])*q_*rdvpar2)/m_; 
  out[8] += -(2.738612787525831*(apardot[1]*fin[4]+apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[9] += -(2.738612787525831*(apardot[0]*fin[4]+apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  out[10] += -(2.738612787525831*(apardot[1]*fin[7]+apardot[0]*fin[6])*q_*rdvpar2)/m_; 
  out[11] += -(2.738612787525831*(apardot[0]*fin[7]+apardot[1]*fin[6])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
