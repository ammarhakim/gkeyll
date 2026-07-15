#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];
  double rdmu2 = 2.0/dxv[3];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[2]*dxv[2];

  const double *bioverJB_x = &bioverJB[0]; 
  const double *bioverJB_y = &bioverJB[4]; 
  const double *bioverJB_z = &bioverJB[8]; 

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0]; 
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[4]; 
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[8]; 

  double hamil[24] = {0.}; 
  hamil[0] = 1.4142135623730951*vmapSq[0]*m_; 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_; 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[3]*hamil[3]; 
  hamil2[1] = hamil[16]*hamil[16]; 

  double alphax[24] = {0.}; 



  double alphaz[24] = {0.}; 
  alphaz[0] = (0.7071067811865475*rtg33inv[0]*hamil[3]*rdz2)/(vmap[1]*m_); 
  alphaz[1] = (0.7071067811865475*rtg33inv[1]*hamil[3]*rdz2)/(vmap[1]*m_); 
  alphaz[2] = (0.7071067811865475*rtg33inv[2]*hamil[3]*rdz2)/(vmap[1]*m_); 
  alphaz[3] = (1.5811388300841895*rtg33inv[0]*hamil[16]*rdz2)/(vmap[1]*m_); 
  alphaz[5] = (0.7071067811865475*hamil[3]*rtg33inv[3]*rdz2)/(vmap[1]*m_); 
  alphaz[6] = (1.5811388300841895*rtg33inv[1]*hamil[16]*rdz2)/(vmap[1]*m_); 
  alphaz[7] = (1.5811388300841895*rtg33inv[2]*hamil[16]*rdz2)/(vmap[1]*m_); 
  alphaz[11] = (1.5811388300841895*rtg33inv[3]*hamil[16]*rdz2)/(vmap[1]*m_); 


  out[2] += 0.4330127018922193*(alphaz[11]*fin[11]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[5]*fin[5]+alphaz[3]*fin[3]+alphaz[2]*fin[2]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[5] += 0.4330127018922193*(alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[3]*fin[6]+fin[3]*alphaz[6]+alphaz[2]*fin[5]+fin[2]*alphaz[5]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[7] += 0.38729833462074165*alphaz[11]*fin[20]+0.3872983346207417*(alphaz[7]*fin[18]+alphaz[6]*fin[17])+0.38729833462074165*alphaz[3]*fin[16]+0.4330127018922193*(alphaz[5]*fin[11]+fin[5]*alphaz[11]+alphaz[2]*fin[7]+fin[2]*alphaz[7]+alphaz[1]*fin[6]+fin[1]*alphaz[6]+alphaz[0]*fin[3]+fin[0]*alphaz[3]); 
  out[9] += 0.4330127018922193*(alphaz[11]*fin[15]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[5]*fin[12]+alphaz[3]*fin[10]+alphaz[2]*fin[9]+alphaz[1]*fin[8]+alphaz[0]*fin[4]); 
  out[11] += 0.38729833462074165*alphaz[7]*fin[20]+0.3872983346207417*(alphaz[11]*fin[18]+alphaz[3]*fin[17])+0.38729833462074165*alphaz[6]*fin[16]+0.4330127018922193*(alphaz[2]*fin[11]+fin[2]*alphaz[11]+alphaz[5]*fin[7]+fin[5]*alphaz[7]+alphaz[0]*fin[6]+fin[0]*alphaz[6]+alphaz[1]*fin[3]+fin[1]*alphaz[3]); 
  out[12] += 0.4330127018922193*(alphaz[7]*fin[15]+alphaz[11]*fin[14]+alphaz[3]*fin[13]+alphaz[2]*fin[12]+alphaz[6]*fin[10]+alphaz[5]*fin[9]+alphaz[0]*fin[8]+alphaz[1]*fin[4]); 
  out[14] += 0.3872983346207417*alphaz[11]*fin[23]+0.38729833462074165*(alphaz[7]*fin[22]+alphaz[6]*fin[21])+0.3872983346207417*alphaz[3]*fin[19]+0.4330127018922193*(alphaz[5]*fin[15]+alphaz[2]*fin[14]+alphaz[1]*fin[13]+alphaz[11]*fin[12]+alphaz[0]*fin[10]+alphaz[7]*fin[9]+alphaz[6]*fin[8]+alphaz[3]*fin[4]); 
  out[15] += 0.3872983346207417*alphaz[7]*fin[23]+0.38729833462074165*(alphaz[11]*fin[22]+alphaz[3]*fin[21])+0.3872983346207417*alphaz[6]*fin[19]+0.4330127018922193*(alphaz[2]*fin[15]+alphaz[5]*fin[14]+alphaz[0]*fin[13]+alphaz[7]*fin[12]+fin[9]*alphaz[11]+alphaz[1]*fin[10]+alphaz[3]*fin[8]+fin[4]*alphaz[6]); 
  out[18] += 0.43301270189221935*alphaz[5]*fin[20]+0.4330127018922193*(alphaz[2]*fin[18]+alphaz[1]*fin[17])+0.43301270189221935*alphaz[0]*fin[16]+0.3872983346207417*(alphaz[11]*fin[11]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[3]*fin[3]); 
  out[20] += 0.4330127018922193*alphaz[2]*fin[20]+0.43301270189221935*(alphaz[5]*fin[18]+alphaz[0]*fin[17])+0.4330127018922193*alphaz[1]*fin[16]+0.38729833462074165*(alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[3]*fin[6]+fin[3]*alphaz[6]); 
  out[22] += 0.43301270189221935*alphaz[5]*fin[23]+0.4330127018922193*(alphaz[2]*fin[22]+alphaz[1]*fin[21])+0.43301270189221935*alphaz[0]*fin[19]+0.38729833462074165*(alphaz[11]*fin[15]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[3]*fin[10]); 
  out[23] += 0.4330127018922193*alphaz[2]*fin[23]+0.43301270189221935*(alphaz[5]*fin[22]+alphaz[0]*fin[21])+0.4330127018922193*alphaz[1]*fin[19]+0.3872983346207417*(alphaz[7]*fin[15]+alphaz[11]*fin[14]+alphaz[3]*fin[13]+alphaz[6]*fin[10]); 

  double alphavpar[24] = {0.}; 
  alphavpar[0] = (-(0.408248290463863*rtg33inv[3]*yfield[30])-0.408248290463863*rtg33inv[2]*yfield[27]-0.408248290463863*rtg33inv[1]*yfield[26]-0.408248290463863*rtg33inv[0]*yfield[25])/(vmap[1]*m_); 
  alphavpar[1] = (-(0.408248290463863*rtg33inv[2]*yfield[30])-0.408248290463863*rtg33inv[3]*yfield[27]-0.408248290463863*rtg33inv[0]*yfield[26]-0.408248290463863*rtg33inv[1]*yfield[25])/(vmap[1]*m_); 
  alphavpar[2] = (-(0.408248290463863*rtg33inv[1]*yfield[30])-0.408248290463863*rtg33inv[0]*yfield[27]-0.408248290463863*rtg33inv[3]*yfield[26]-0.408248290463863*rtg33inv[2]*yfield[25])/(vmap[1]*m_); 
  alphavpar[3] = (-(0.408248290463863*rtg33inv[3]*yfield[36])-0.408248290463863*rtg33inv[2]*yfield[32]-0.408248290463863*rtg33inv[1]*yfield[31]-0.408248290463863*rtg33inv[0]*yfield[28])/(vmap[1]*m_); 
  alphavpar[4] = (-(0.408248290463863*rtg33inv[3]*yfield[37])-0.408248290463863*rtg33inv[2]*yfield[34]-0.408248290463863*rtg33inv[1]*yfield[33]-0.408248290463863*rtg33inv[0]*yfield[29])/(vmap[1]*m_); 
  alphavpar[5] = (-(0.408248290463863*rtg33inv[0]*yfield[30])-0.408248290463863*rtg33inv[1]*yfield[27]-0.408248290463863*rtg33inv[2]*yfield[26]-0.408248290463863*rtg33inv[3]*yfield[25])/(vmap[1]*m_); 
  alphavpar[6] = (-(0.408248290463863*rtg33inv[2]*yfield[36])-0.408248290463863*rtg33inv[3]*yfield[32]-0.408248290463863*rtg33inv[0]*yfield[31]-0.408248290463863*rtg33inv[1]*yfield[28])/(vmap[1]*m_); 
  alphavpar[7] = (-(0.408248290463863*rtg33inv[1]*yfield[36])-0.408248290463863*rtg33inv[0]*yfield[32]-0.408248290463863*rtg33inv[3]*yfield[31]-0.408248290463863*rtg33inv[2]*yfield[28])/(vmap[1]*m_); 
  alphavpar[8] = (-(0.408248290463863*rtg33inv[2]*yfield[37])-0.408248290463863*rtg33inv[3]*yfield[34]-0.408248290463863*rtg33inv[0]*yfield[33]-0.408248290463863*rtg33inv[1]*yfield[29])/(vmap[1]*m_); 
  alphavpar[9] = (-(0.408248290463863*rtg33inv[1]*yfield[37])-0.408248290463863*rtg33inv[0]*yfield[34]-0.408248290463863*rtg33inv[3]*yfield[33]-0.408248290463863*rtg33inv[2]*yfield[29])/(vmap[1]*m_); 
  alphavpar[10] = (-(0.408248290463863*rtg33inv[3]*yfield[40])-0.408248290463863*rtg33inv[2]*yfield[39]-0.408248290463863*rtg33inv[1]*yfield[38]-0.408248290463863*rtg33inv[0]*yfield[35])/(vmap[1]*m_); 
  alphavpar[11] = (-(0.408248290463863*rtg33inv[0]*yfield[36])-0.408248290463863*rtg33inv[1]*yfield[32]-0.408248290463863*rtg33inv[2]*yfield[31]-0.408248290463863*rtg33inv[3]*yfield[28])/(vmap[1]*m_); 
  alphavpar[12] = (-(0.408248290463863*rtg33inv[0]*yfield[37])-0.408248290463863*rtg33inv[1]*yfield[34]-0.408248290463863*rtg33inv[2]*yfield[33]-0.408248290463863*rtg33inv[3]*yfield[29])/(vmap[1]*m_); 
  alphavpar[13] = (-(0.408248290463863*rtg33inv[2]*yfield[40])-0.408248290463863*rtg33inv[3]*yfield[39]-0.408248290463863*rtg33inv[0]*yfield[38]-0.408248290463863*rtg33inv[1]*yfield[35])/(vmap[1]*m_); 
  alphavpar[14] = (-(0.408248290463863*rtg33inv[1]*yfield[40])-0.408248290463863*rtg33inv[0]*yfield[39]-0.408248290463863*rtg33inv[3]*yfield[38]-0.408248290463863*rtg33inv[2]*yfield[35])/(vmap[1]*m_); 
  alphavpar[15] = (-(0.408248290463863*rtg33inv[0]*yfield[40])-0.408248290463863*rtg33inv[1]*yfield[39]-0.408248290463863*rtg33inv[2]*yfield[38]-0.408248290463863*rtg33inv[3]*yfield[35])/(vmap[1]*m_); 
  alphavpar[16] = (-(0.408248290463863*rtg33inv[3]*yfield[45])-0.40824829046386296*rtg33inv[2]*yfield[43]-0.40824829046386296*rtg33inv[1]*yfield[42]-0.408248290463863*rtg33inv[0]*yfield[41])/(vmap[1]*m_); 
  alphavpar[17] = (-(0.40824829046386296*rtg33inv[2]*yfield[45])-0.408248290463863*rtg33inv[3]*yfield[43]-0.408248290463863*rtg33inv[0]*yfield[42]-0.40824829046386296*rtg33inv[1]*yfield[41])/(vmap[1]*m_); 
  alphavpar[18] = (-(0.40824829046386296*rtg33inv[1]*yfield[45])-0.408248290463863*rtg33inv[0]*yfield[43]-0.408248290463863*rtg33inv[3]*yfield[42]-0.40824829046386296*rtg33inv[2]*yfield[41])/(vmap[1]*m_); 
  alphavpar[19] = (-(0.408248290463863*rtg33inv[3]*yfield[48])-0.40824829046386296*rtg33inv[2]*yfield[47]-0.40824829046386296*rtg33inv[1]*yfield[46]-0.408248290463863*rtg33inv[0]*yfield[44])/(vmap[1]*m_); 
  alphavpar[20] = (-(0.408248290463863*rtg33inv[0]*yfield[45])-0.40824829046386296*rtg33inv[1]*yfield[43]-0.40824829046386296*rtg33inv[2]*yfield[42]-0.408248290463863*rtg33inv[3]*yfield[41])/(vmap[1]*m_); 
  alphavpar[21] = (-(0.40824829046386296*rtg33inv[2]*yfield[48])-0.408248290463863*rtg33inv[3]*yfield[47]-0.408248290463863*rtg33inv[0]*yfield[46]-0.40824829046386296*rtg33inv[1]*yfield[44])/(vmap[1]*m_); 
  alphavpar[22] = (-(0.40824829046386296*rtg33inv[1]*yfield[48])-0.408248290463863*rtg33inv[0]*yfield[47]-0.408248290463863*rtg33inv[3]*yfield[46]-0.40824829046386296*rtg33inv[2]*yfield[44])/(vmap[1]*m_); 
  alphavpar[23] = (-(0.408248290463863*rtg33inv[0]*yfield[48])-0.40824829046386296*rtg33inv[1]*yfield[47]-0.40824829046386296*rtg33inv[2]*yfield[46]-0.408248290463863*rtg33inv[3]*yfield[44])/(vmap[1]*m_); 


  out[3] += 0.4330127018922193*(alphavpar[23]*fin[23]+alphavpar[22]*fin[22]+alphavpar[21]*fin[21]+alphavpar[20]*fin[20]+alphavpar[19]*fin[19]+alphavpar[18]*fin[18]+alphavpar[17]*fin[17]+alphavpar[16]*fin[16]+alphavpar[15]*fin[15]+alphavpar[14]*fin[14]+alphavpar[13]*fin[13]+alphavpar[12]*fin[12]+alphavpar[11]*fin[11]+alphavpar[10]*fin[10]+alphavpar[9]*fin[9]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[6] += 0.43301270189221935*(alphavpar[22]*fin[23]+fin[22]*alphavpar[23]+alphavpar[19]*fin[21]+fin[19]*alphavpar[21]+alphavpar[18]*fin[20]+fin[18]*alphavpar[20]+alphavpar[16]*fin[17]+fin[16]*alphavpar[17])+0.4330127018922193*(alphavpar[14]*fin[15]+fin[14]*alphavpar[15]+alphavpar[10]*fin[13]+fin[10]*alphavpar[13]+alphavpar[9]*fin[12]+fin[9]*alphavpar[12]+alphavpar[7]*fin[11]+fin[7]*alphavpar[11]+alphavpar[4]*fin[8]+fin[4]*alphavpar[8]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.43301270189221935*(alphavpar[21]*fin[23]+fin[21]*alphavpar[23]+alphavpar[19]*fin[22]+fin[19]*alphavpar[22]+alphavpar[17]*fin[20]+fin[17]*alphavpar[20]+alphavpar[16]*fin[18]+fin[16]*alphavpar[18])+0.4330127018922193*(alphavpar[13]*fin[15]+fin[13]*alphavpar[15]+alphavpar[10]*fin[14]+fin[10]*alphavpar[14]+alphavpar[8]*fin[12]+fin[8]*alphavpar[12]+alphavpar[6]*fin[11]+fin[6]*alphavpar[11]+alphavpar[4]*fin[9]+fin[4]*alphavpar[9]+alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[10] += 0.43301270189221935*(alphavpar[20]*fin[23]+fin[20]*alphavpar[23]+alphavpar[18]*fin[22]+fin[18]*alphavpar[22]+alphavpar[17]*fin[21]+fin[17]*alphavpar[21]+alphavpar[16]*fin[19]+fin[16]*alphavpar[19])+0.4330127018922193*(alphavpar[11]*fin[15]+fin[11]*alphavpar[15]+alphavpar[7]*fin[14]+fin[7]*alphavpar[14]+alphavpar[6]*fin[13]+fin[6]*alphavpar[13]+alphavpar[5]*fin[12]+fin[5]*alphavpar[12]+alphavpar[3]*fin[10]+fin[3]*alphavpar[10]+alphavpar[2]*fin[9]+fin[2]*alphavpar[9]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphavpar[19]*fin[23]+fin[19]*alphavpar[23]+alphavpar[21]*fin[22]+fin[21]*alphavpar[22]+alphavpar[16]*fin[20]+fin[16]*alphavpar[20]+alphavpar[17]*fin[18]+fin[17]*alphavpar[18]+alphavpar[10]*fin[15]+fin[10]*alphavpar[15]+alphavpar[13]*fin[14]+fin[13]*alphavpar[14]+alphavpar[4]*fin[12]+fin[4]*alphavpar[12]+alphavpar[3]*fin[11]+fin[3]*alphavpar[11]+alphavpar[8]*fin[9]+fin[8]*alphavpar[9]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[13] += 0.4330127018922193*(alphavpar[18]*fin[23]+fin[18]*alphavpar[23]+alphavpar[20]*fin[22]+fin[20]*alphavpar[22]+alphavpar[16]*fin[21]+fin[16]*alphavpar[21]+alphavpar[17]*fin[19]+fin[17]*alphavpar[19]+alphavpar[7]*fin[15]+fin[7]*alphavpar[15]+alphavpar[11]*fin[14]+fin[11]*alphavpar[14]+alphavpar[3]*fin[13]+fin[3]*alphavpar[13]+alphavpar[2]*fin[12]+fin[2]*alphavpar[12]+alphavpar[6]*fin[10]+fin[6]*alphavpar[10]+alphavpar[5]*fin[9]+fin[5]*alphavpar[9]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*(alphavpar[17]*fin[23]+fin[17]*alphavpar[23]+alphavpar[16]*fin[22]+fin[16]*alphavpar[22]+alphavpar[20]*fin[21]+fin[20]*alphavpar[21]+alphavpar[18]*fin[19]+fin[18]*alphavpar[19]+alphavpar[6]*fin[15]+fin[6]*alphavpar[15]+alphavpar[3]*fin[14]+fin[3]*alphavpar[14]+alphavpar[11]*fin[13]+fin[11]*alphavpar[13]+alphavpar[1]*fin[12]+fin[1]*alphavpar[12]+alphavpar[7]*fin[10]+fin[7]*alphavpar[10]+alphavpar[0]*fin[9]+fin[0]*alphavpar[9]+alphavpar[5]*fin[8]+fin[5]*alphavpar[8]+alphavpar[2]*fin[4]+fin[2]*alphavpar[4]); 
  out[15] += 0.43301270189221935*(alphavpar[16]*fin[23]+fin[16]*alphavpar[23]+alphavpar[17]*fin[22]+fin[17]*alphavpar[22]+alphavpar[18]*fin[21]+fin[18]*alphavpar[21]+alphavpar[19]*fin[20]+fin[19]*alphavpar[20])+0.4330127018922193*(alphavpar[3]*fin[15]+fin[3]*alphavpar[15]+alphavpar[6]*fin[14]+fin[6]*alphavpar[14]+alphavpar[7]*fin[13]+fin[7]*alphavpar[13]+alphavpar[0]*fin[12]+fin[0]*alphavpar[12]+alphavpar[10]*fin[11]+fin[10]*alphavpar[11]+alphavpar[1]*fin[9]+fin[1]*alphavpar[9]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]); 
  out[16] += 0.8660254037844387*(alphavpar[15]*fin[23]+fin[15]*alphavpar[23])+0.8660254037844386*(alphavpar[14]*fin[22]+fin[14]*alphavpar[22]+alphavpar[13]*fin[21]+fin[13]*alphavpar[21]+alphavpar[11]*fin[20]+fin[11]*alphavpar[20])+0.8660254037844387*(alphavpar[10]*fin[19]+fin[10]*alphavpar[19]+alphavpar[7]*fin[18]+fin[7]*alphavpar[18]+alphavpar[6]*fin[17]+fin[6]*alphavpar[17])+0.8660254037844386*(alphavpar[3]*fin[16]+fin[3]*alphavpar[16])+0.9682458365518543*(alphavpar[12]*fin[15]+fin[12]*alphavpar[15]+alphavpar[9]*fin[14]+fin[9]*alphavpar[14]+alphavpar[8]*fin[13]+fin[8]*alphavpar[13]+alphavpar[5]*fin[11]+fin[5]*alphavpar[11]+alphavpar[4]*fin[10]+fin[4]*alphavpar[10]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[17] += 0.8660254037844386*(alphavpar[14]*fin[23]+fin[14]*alphavpar[23])+0.8660254037844387*(alphavpar[15]*fin[22]+fin[15]*alphavpar[22]+alphavpar[10]*fin[21]+fin[10]*alphavpar[21]+alphavpar[7]*fin[20]+fin[7]*alphavpar[20])+0.8660254037844386*(alphavpar[13]*fin[19]+fin[13]*alphavpar[19]+alphavpar[11]*fin[18]+fin[11]*alphavpar[18]+alphavpar[3]*fin[17]+fin[3]*alphavpar[17])+0.8660254037844387*(alphavpar[6]*fin[16]+fin[6]*alphavpar[16])+0.9682458365518543*(alphavpar[9]*fin[15]+fin[9]*alphavpar[15]+alphavpar[12]*fin[14]+fin[12]*alphavpar[14]+alphavpar[4]*fin[13]+fin[4]*alphavpar[13]+alphavpar[2]*fin[11]+fin[2]*alphavpar[11]+alphavpar[8]*fin[10]+fin[8]*alphavpar[10]+alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[18] += 0.8660254037844386*(alphavpar[13]*fin[23]+fin[13]*alphavpar[23])+0.8660254037844387*(alphavpar[10]*fin[22]+fin[10]*alphavpar[22]+alphavpar[15]*fin[21]+fin[15]*alphavpar[21]+alphavpar[6]*fin[20]+fin[6]*alphavpar[20])+0.8660254037844386*(alphavpar[14]*fin[19]+fin[14]*alphavpar[19]+alphavpar[3]*fin[18]+fin[3]*alphavpar[18]+alphavpar[11]*fin[17]+fin[11]*alphavpar[17])+0.8660254037844387*(alphavpar[7]*fin[16]+fin[7]*alphavpar[16])+0.9682458365518543*(alphavpar[8]*fin[15]+fin[8]*alphavpar[15]+alphavpar[4]*fin[14]+fin[4]*alphavpar[14]+alphavpar[12]*fin[13]+fin[12]*alphavpar[13]+alphavpar[1]*fin[11]+fin[1]*alphavpar[11]+alphavpar[9]*fin[10]+fin[9]*alphavpar[10]+alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[19] += 0.8660254037844386*(alphavpar[11]*fin[23]+fin[11]*alphavpar[23])+0.8660254037844387*(alphavpar[7]*fin[22]+fin[7]*alphavpar[22]+alphavpar[6]*fin[21]+fin[6]*alphavpar[21]+alphavpar[15]*fin[20]+fin[15]*alphavpar[20])+0.8660254037844386*(alphavpar[3]*fin[19]+fin[3]*alphavpar[19]+alphavpar[14]*fin[18]+fin[14]*alphavpar[18]+alphavpar[13]*fin[17]+fin[13]*alphavpar[17])+0.8660254037844387*(alphavpar[10]*fin[16]+fin[10]*alphavpar[16])+0.9682458365518543*(alphavpar[5]*fin[15]+fin[5]*alphavpar[15]+alphavpar[2]*fin[14]+fin[2]*alphavpar[14]+alphavpar[1]*fin[13]+fin[1]*alphavpar[13]+alphavpar[11]*fin[12]+fin[11]*alphavpar[12]+alphavpar[0]*fin[10]+fin[0]*alphavpar[10]+alphavpar[7]*fin[9]+fin[7]*alphavpar[9]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4]); 
  out[20] += 0.8660254037844387*(alphavpar[10]*fin[23]+fin[10]*alphavpar[23])+0.8660254037844386*(alphavpar[13]*fin[22]+fin[13]*alphavpar[22]+alphavpar[14]*fin[21]+fin[14]*alphavpar[21]+alphavpar[3]*fin[20]+fin[3]*alphavpar[20])+0.8660254037844387*(alphavpar[15]*fin[19]+fin[15]*alphavpar[19]+alphavpar[6]*fin[18]+fin[6]*alphavpar[18]+alphavpar[7]*fin[17]+fin[7]*alphavpar[17])+0.8660254037844386*(alphavpar[11]*fin[16]+fin[11]*alphavpar[16])+0.9682458365518543*(alphavpar[4]*fin[15]+fin[4]*alphavpar[15]+alphavpar[8]*fin[14]+fin[8]*alphavpar[14]+alphavpar[9]*fin[13]+fin[9]*alphavpar[13]+alphavpar[10]*fin[12]+fin[10]*alphavpar[12]+alphavpar[0]*fin[11]+fin[0]*alphavpar[11]+alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]); 
  out[21] += 0.8660254037844387*(alphavpar[7]*fin[23]+fin[7]*alphavpar[23])+0.8660254037844386*(alphavpar[11]*fin[22]+fin[11]*alphavpar[22]+alphavpar[3]*fin[21]+fin[3]*alphavpar[21]+alphavpar[14]*fin[20]+fin[14]*alphavpar[20])+0.8660254037844387*(alphavpar[6]*fin[19]+fin[6]*alphavpar[19]+alphavpar[15]*fin[18]+fin[15]*alphavpar[18]+alphavpar[10]*fin[17]+fin[10]*alphavpar[17])+0.8660254037844386*(alphavpar[13]*fin[16]+fin[13]*alphavpar[16])+0.9682458365518543*(alphavpar[2]*fin[15]+fin[2]*alphavpar[15]+alphavpar[5]*fin[14]+fin[5]*alphavpar[14]+alphavpar[0]*fin[13]+fin[0]*alphavpar[13]+alphavpar[7]*fin[12]+fin[7]*alphavpar[12]+alphavpar[9]*fin[11]+fin[9]*alphavpar[11]+alphavpar[1]*fin[10]+fin[1]*alphavpar[10]+alphavpar[3]*fin[8]+fin[3]*alphavpar[8]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6]); 
  out[22] += 0.8660254037844387*(alphavpar[6]*fin[23]+fin[6]*alphavpar[23])+0.8660254037844386*(alphavpar[3]*fin[22]+fin[3]*alphavpar[22]+alphavpar[11]*fin[21]+fin[11]*alphavpar[21]+alphavpar[13]*fin[20]+fin[13]*alphavpar[20])+0.8660254037844387*(alphavpar[7]*fin[19]+fin[7]*alphavpar[19]+alphavpar[10]*fin[18]+fin[10]*alphavpar[18]+alphavpar[15]*fin[17]+fin[15]*alphavpar[17])+0.8660254037844386*(alphavpar[14]*fin[16]+fin[14]*alphavpar[16])+0.9682458365518543*(alphavpar[1]*fin[15]+fin[1]*alphavpar[15]+alphavpar[0]*fin[14]+fin[0]*alphavpar[14]+alphavpar[5]*fin[13]+fin[5]*alphavpar[13]+alphavpar[6]*fin[12]+fin[6]*alphavpar[12]+alphavpar[8]*fin[11]+fin[8]*alphavpar[11]+alphavpar[2]*fin[10]+fin[2]*alphavpar[10]+alphavpar[3]*fin[9]+fin[3]*alphavpar[9]+alphavpar[4]*fin[7]+fin[4]*alphavpar[7]); 
  out[23] += 0.8660254037844386*(alphavpar[3]*fin[23]+fin[3]*alphavpar[23])+0.8660254037844387*(alphavpar[6]*fin[22]+fin[6]*alphavpar[22]+alphavpar[7]*fin[21]+fin[7]*alphavpar[21]+alphavpar[10]*fin[20]+fin[10]*alphavpar[20])+0.8660254037844386*(alphavpar[11]*fin[19]+fin[11]*alphavpar[19]+alphavpar[13]*fin[18]+fin[13]*alphavpar[18]+alphavpar[14]*fin[17]+fin[14]*alphavpar[17])+0.8660254037844387*(alphavpar[15]*fin[16]+fin[15]*alphavpar[16])+0.9682458365518543*(alphavpar[0]*fin[15]+fin[0]*alphavpar[15]+alphavpar[1]*fin[14]+fin[1]*alphavpar[14]+alphavpar[2]*fin[13]+fin[2]*alphavpar[13]+alphavpar[3]*fin[12]+fin[3]*alphavpar[12]+alphavpar[4]*fin[11]+fin[4]*alphavpar[11]+alphavpar[5]*fin[10]+fin[5]*alphavpar[10]+alphavpar[6]*fin[9]+fin[6]*alphavpar[9]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8]); 

  return 0.; 
} 
