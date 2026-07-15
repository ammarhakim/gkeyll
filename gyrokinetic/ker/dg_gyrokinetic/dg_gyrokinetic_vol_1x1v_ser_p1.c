#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *yfield,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // yfield: Yushmanov field (gradient of Yushmanov potential).
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
  hamil[0] = 0.7071067811865475*vmapSq[0]*m_; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[4]*hamil[4]; 

  double alphax[6] = {0.}; 
  alphax[0] = (rdx2*((3.5355339059327373*dualcurlbhatoverB_z[0]*hamil2[1])/vmap2+(0.7071067811865475*dualcurlbhatoverB_z[0]*hamil2[0])/vmap2))/(m_*q_)+(rtg33inv[0]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[1] = (rdx2*((3.5355339059327373*dualcurlbhatoverB_z[1]*hamil2[1])/vmap2+(0.7071067811865475*hamil2[0]*dualcurlbhatoverB_z[1])/vmap2))/(m_*q_)+(rtg33inv[1]*vmap[1]*hamil[2]*rdx2)/(m_*vmap2); 
  alphax[2] = (3.1622776601683795*dualcurlbhatoverB_z[0]*hamil[2]*hamil[4]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[0]*vmap[1]*hamil[4]*rdx2)/(m_*vmap2); 
  alphax[3] = (3.1622776601683795*dualcurlbhatoverB_z[1]*hamil[2]*hamil[4]*rdx2)/(m_*q_*vmap2)+(2.23606797749979*rtg33inv[1]*vmap[1]*hamil[4]*rdx2)/(m_*vmap2); 
  alphax[4] = (3.1622776601683795*dualcurlbhatoverB_z[0]*hamil2[1]*rdx2)/(m_*q_*vmap2); 
  alphax[5] = (3.16227766016838*dualcurlbhatoverB_z[1]*hamil2[1]*rdx2)/(m_*q_*vmap2); 


  out[1] += 0.8660254037844386*(alphax[5]*fin[5]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[3] += 0.7745966692414834*(alphax[3]*fin[5]+fin[3]*alphax[5])+0.7745966692414833*(alphax[2]*fin[4]+fin[2]*alphax[4])+0.8660254037844386*(alphax[1]*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[5] += 0.5532833351724881*alphax[5]*fin[5]+0.8660254037844386*(alphax[1]*fin[5]+fin[1]*alphax[5])+0.5532833351724881*alphax[4]*fin[4]+0.8660254037844387*(alphax[0]*fin[4]+fin[0]*alphax[4])+0.7745966692414834*(alphax[3]*fin[3]+alphax[2]*fin[2]); 

  double alphavpar[6] = {0.}; 
  alphavpar[0] = (-((0.9128709291752769*dualcurlbhatoverB_z[1]*hamil[4]*yfield[4])/vmap2)-(0.9128709291752769*dualcurlbhatoverB_z[0]*yfield[3]*hamil[4])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[1]*hamil[2]*yfield[2])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[0]*yfield[1]*hamil[2])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896258*rtg33inv[1]*yfield[2])/vmap2)-(0.5773502691896258*rtg33inv[0]*yfield[1])/vmap2))/m_; 
  alphavpar[1] = (-((0.9128709291752769*dualcurlbhatoverB_z[0]*hamil[4]*yfield[4])/vmap2)-(0.9128709291752769*dualcurlbhatoverB_z[1]*yfield[3]*hamil[4])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[0]*hamil[2]*yfield[2])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[1]*yfield[1]*hamil[2])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896258*rtg33inv[0]*yfield[2])/vmap2)-(0.5773502691896258*rtg33inv[1]*yfield[1])/vmap2))/m_; 
  alphavpar[2] = (-((0.816496580927726*dualcurlbhatoverB_z[1]*hamil[4]*yfield[6])/vmap2)-(0.8164965809277261*dualcurlbhatoverB_z[0]*hamil[4]*yfield[5])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[1]*hamil[2]*yfield[4])/vmap2-(0.9128709291752769*dualcurlbhatoverB_z[1]*yfield[2]*hamil[4])/vmap2-(0.9128709291752769*dualcurlbhatoverB_z[0]*yfield[1]*hamil[4])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[0]*hamil[2]*yfield[3])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896258*rtg33inv[1]*yfield[4])/vmap2)-(0.5773502691896258*rtg33inv[0]*yfield[3])/vmap2))/m_; 
  alphavpar[3] = (-((0.816496580927726*dualcurlbhatoverB_z[0]*hamil[4]*yfield[6])/vmap2)-(0.8164965809277261*dualcurlbhatoverB_z[1]*hamil[4]*yfield[5])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[0]*hamil[2]*yfield[4])/vmap2-(0.9128709291752769*dualcurlbhatoverB_z[0]*yfield[2]*hamil[4])/vmap2-(0.9128709291752769*dualcurlbhatoverB_z[1]*yfield[1]*hamil[4])/vmap2-(0.408248290463863*dualcurlbhatoverB_z[1]*hamil[2]*yfield[3])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896258*rtg33inv[0]*yfield[4])/vmap2)-(0.5773502691896258*rtg33inv[1]*yfield[3])/vmap2))/m_; 
  alphavpar[4] = (-((0.40824829046386296*dualcurlbhatoverB_z[1]*hamil[2]*yfield[6])/vmap2)-(0.408248290463863*dualcurlbhatoverB_z[0]*hamil[2]*yfield[5])/vmap2-(0.8164965809277261*dualcurlbhatoverB_z[1]*hamil[4]*yfield[4])/vmap2-(0.8164965809277261*dualcurlbhatoverB_z[0]*yfield[3]*hamil[4])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896257*rtg33inv[1]*yfield[6])/vmap2)-(0.5773502691896258*rtg33inv[0]*yfield[5])/vmap2))/m_; 
  alphavpar[5] = (-((0.408248290463863*dualcurlbhatoverB_z[0]*hamil[2]*yfield[6])/vmap2)-(0.40824829046386296*dualcurlbhatoverB_z[1]*hamil[2]*yfield[5])/vmap2-(0.816496580927726*dualcurlbhatoverB_z[0]*hamil[4]*yfield[4])/vmap2-(0.816496580927726*dualcurlbhatoverB_z[1]*yfield[3]*hamil[4])/vmap2)/(m_*q_)+(vmap[1]*(-((0.5773502691896258*rtg33inv[0]*yfield[6])/vmap2)-(0.5773502691896257*rtg33inv[1]*yfield[5])/vmap2))/m_; 


  out[2] += 0.8660254037844386*(alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.8660254037844387*(alphavpar[4]*fin[5]+fin[4]*alphavpar[5])+0.8660254037844386*(alphavpar[2]*fin[3]+fin[2]*alphavpar[3]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.7320508075688774*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5])+1.7320508075688772*(alphavpar[2]*fin[4]+fin[2]*alphavpar[4])+1.9364916731037085*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[5] += 1.7320508075688772*(alphavpar[2]*fin[5]+fin[2]*alphavpar[5])+1.7320508075688774*(alphavpar[3]*fin[4]+fin[3]*alphavpar[4])+1.9364916731037085*(alphavpar[0]*fin[3]+fin[0]*alphavpar[3]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 

  return 0.; 
} 
