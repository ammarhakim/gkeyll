#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[1]*dxv[1];

  const double *bioverJB_x = &bioverJB[0];
  const double *bioverJB_y = &bioverJB[2];
  const double *bioverJB_z = &bioverJB[4];

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0];
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[2];
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[4];

  double hamil[6] = {0.}; 
  hamil[0] = 1.414213562373095*phi[0]*q_+0.7071067811865475*vmapSq[0]*m_; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[4]*hamil[4]; 

  double alphax[6] = {0.}; 
  alphax[0] = (rdx2*((3.535533905932737*dualcurlbhatoverB_z[0]*hamil2[1])/vmap2+(0.7071067811865475*dualcurlbhatoverB_z[0]*hamil2[0])/vmap2))/(m_*q_)+(rtg33inv[0]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[1] = (rdx2*((3.535533905932737*dualcurlbhatoverB_z[1]*hamil2[1])/vmap2+(0.7071067811865475*hamil2[0]*dualcurlbhatoverB_z[1])/vmap2))/(m_*q_)+(rtg33inv[1]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[2] = (3.16227766016838*dualcurlbhatoverB_z[0]*hamil[2]*hamil[4]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[0]*vmap[1]*hamil[4]*rdx2)/(m_*vmap2); 
  alphax[3] = (3.16227766016838*dualcurlbhatoverB_z[1]*hamil[2]*hamil[4]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[1]*vmap[1]*hamil[4]*rdx2)/(m_*vmap2); 
  alphax[4] = (3.16227766016838*dualcurlbhatoverB_z[0]*hamil2[1]*rdx2)/(m_*q_*vmap2); 
  alphax[5] = (3.16227766016838*dualcurlbhatoverB_z[1]*hamil2[1]*rdx2)/(m_*q_*vmap2); 


  out[1] += 0.8660254037844386*(alphax[5]*fin[5]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[3] += 0.7745966692414834*(alphax[3]*fin[5]+fin[3]*alphax[5])+0.7745966692414833*(alphax[2]*fin[4]+fin[2]*alphax[4])+0.8660254037844386*(alphax[1]*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.5532833351724881*alphax[5]*fin[5]+0.8660254037844386*(alphax[1]*fin[5]+fin[1]*alphax[5])+0.5532833351724881*alphax[4]*fin[4]+0.8660254037844387*(alphax[0]*fin[4]+fin[0]*alphax[4])+0.7745966692414834*(alphax[3]*fin[3]+alphax[2]*fin[2]); 

  double alphavpar[6] = {0.}; 
  alphavpar[0] = (-(0.7071067811865475*dualcurlbhatoverB_z[0]*hamil[1]*hamil[2]*rdx2)/(m_*q_*vmap2))-(1.0*rtg33inv[0]*hamil[1]*vmap[1]*rdx2)/(m_*vmap2); 
  alphavpar[1] = (-(0.7071067811865475*dualcurlbhatoverB_z[1]*hamil[1]*hamil[2]*rdx2)/(m_*q_*vmap2))-(1.0*hamil[1]*rtg33inv[1]*vmap[1]*rdx2)/(m_*vmap2); 
  alphavpar[2] = -(1.58113883008419*dualcurlbhatoverB_z[0]*hamil[1]*hamil[4]*rdx2)/(m_*q_*vmap2); 
  alphavpar[3] = -(1.58113883008419*dualcurlbhatoverB_z[1]*hamil[1]*hamil[4]*rdx2)/(m_*q_*vmap2); 


  out[2] += 0.8660254037844386*(alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.8660254037844386*(alphavpar[2]*fin[3]+fin[2]*alphavpar[3]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.732050807568877*alphavpar[3]*fin[5]+1.732050807568877*alphavpar[2]*fin[4]+1.936491673103709*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[5] += 1.732050807568877*alphavpar[2]*fin[5]+1.732050807568877*alphavpar[3]*fin[4]+1.936491673103709*(alphavpar[0]*fin[3]+fin[0]*alphavpar[3]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 

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
  out[4] += -(2.73861278752583*(apardot[1]*fin[3]+apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[5] += -(2.73861278752583*(apardot[0]*fin[3]+apardot[1]*fin[2])*q_*rdvpar2)/m_; 
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
