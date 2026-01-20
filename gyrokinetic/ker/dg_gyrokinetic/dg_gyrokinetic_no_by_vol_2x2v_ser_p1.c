#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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

  const double *dual_xx = &dualcurlbhatoverB[0]; 
  const double *dual_xy = &dualcurlbhatoverB[4]; 
  const double *dual_xz = &dualcurlbhatoverB[8]; 
  const double *dual_yx = &dualcurlbhatoverB[12]; 
  const double *dual_yy = &dualcurlbhatoverB[16]; 
  const double *dual_yz = &dualcurlbhatoverB[20]; 
  const double *dual_zx = &dualcurlbhatoverB[24]; 
  const double *dual_zy = &dualcurlbhatoverB[28]; 
  const double *dual_zz = &dualcurlbhatoverB[32]; 

  double hamil[24] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+1.414213562373095*(vmapSq[0]*m_+bmag[0]*vmap[2]); 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_+1.414213562373095*bmag[2]*vmap[2]; 
  hamil[3] = 1.414213562373095*vmapSq[1]*m_; 
  hamil[4] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[8] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[9] = 1.414213562373095*bmag[2]*vmap[3]; 
  hamil[12] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[16] = 1.414213562373095*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[3]*hamil[3]; 
  hamil2[1] = hamil[16]*hamil[16]; 

  double alphax[24] = {0.}; 
  alphax[0] = 200000.0*dual_xz[0]*rdx2; 
  alphax[1] = 200000.0*dual_xz[1]*rdx2; 
  alphax[2] = 200000.0*dual_xz[2]*rdx2; 
  alphax[5] = 200000.0*dual_xz[3]*rdx2; 


  out[1] += 0.4330127018922193*(alphax[5]*fin[5]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[5] += 0.4330127018922193*(alphax[1]*fin[5]+fin[1]*alphax[5]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[6] += 0.4330127018922193*(alphax[5]*fin[11]+alphax[2]*fin[7]+alphax[1]*fin[6]+alphax[0]*fin[3]); 
  out[8] += 0.4330127018922193*(alphax[5]*fin[12]+alphax[2]*fin[9]+alphax[1]*fin[8]+alphax[0]*fin[4]); 
  out[11] += 0.4330127018922193*(alphax[1]*fin[11]+alphax[0]*fin[7]+alphax[5]*fin[6]+alphax[2]*fin[3]); 
  out[12] += 0.4330127018922193*(alphax[1]*fin[12]+alphax[0]*fin[9]+alphax[5]*fin[8]+alphax[2]*fin[4]); 
  out[13] += 0.4330127018922193*(alphax[5]*fin[15]+alphax[2]*fin[14]+alphax[1]*fin[13]+alphax[0]*fin[10]); 
  out[15] += 0.4330127018922193*(alphax[1]*fin[15]+alphax[0]*fin[14]+alphax[5]*fin[13]+alphax[2]*fin[10]); 
  out[17] += 0.4330127018922194*alphax[5]*fin[20]+0.4330127018922193*(alphax[2]*fin[18]+alphax[1]*fin[17])+0.4330127018922194*alphax[0]*fin[16]; 
  out[20] += 0.4330127018922193*alphax[1]*fin[20]+0.4330127018922194*(alphax[0]*fin[18]+alphax[5]*fin[17])+0.4330127018922193*alphax[2]*fin[16]; 
  out[21] += 0.4330127018922194*alphax[5]*fin[23]+0.4330127018922193*(alphax[2]*fin[22]+alphax[1]*fin[21])+0.4330127018922194*alphax[0]*fin[19]; 
  out[23] += 0.4330127018922193*alphax[1]*fin[23]+0.4330127018922194*(alphax[0]*fin[22]+alphax[5]*fin[21])+0.4330127018922193*alphax[2]*fin[19]; 

  double alphaz[24] = {0.}; 
  alphaz[0] = 200000.0*dual_zz[0]*rdz2; 
  alphaz[1] = 200000.0*dual_zz[1]*rdz2; 
  alphaz[2] = 200000.0*dual_zz[2]*rdz2; 
  alphaz[5] = 200000.0*dual_zz[3]*rdz2; 


  out[2] += 0.4330127018922193*(alphaz[5]*fin[5]+alphaz[2]*fin[2]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[5] += 0.4330127018922193*(alphaz[2]*fin[5]+fin[2]*alphaz[5]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[7] += 0.4330127018922193*(alphaz[5]*fin[11]+alphaz[2]*fin[7]+alphaz[1]*fin[6]+alphaz[0]*fin[3]); 
  out[9] += 0.4330127018922193*(alphaz[5]*fin[12]+alphaz[2]*fin[9]+alphaz[1]*fin[8]+alphaz[0]*fin[4]); 
  out[11] += 0.4330127018922193*(alphaz[2]*fin[11]+alphaz[5]*fin[7]+alphaz[0]*fin[6]+alphaz[1]*fin[3]); 
  out[12] += 0.4330127018922193*(alphaz[2]*fin[12]+alphaz[5]*fin[9]+alphaz[0]*fin[8]+alphaz[1]*fin[4]); 
  out[14] += 0.4330127018922193*(alphaz[5]*fin[15]+alphaz[2]*fin[14]+alphaz[1]*fin[13]+alphaz[0]*fin[10]); 
  out[15] += 0.4330127018922193*(alphaz[2]*fin[15]+alphaz[5]*fin[14]+alphaz[0]*fin[13]+alphaz[1]*fin[10]); 
  out[18] += 0.4330127018922194*alphaz[5]*fin[20]+0.4330127018922193*(alphaz[2]*fin[18]+alphaz[1]*fin[17])+0.4330127018922194*alphaz[0]*fin[16]; 
  out[20] += 0.4330127018922193*alphaz[2]*fin[20]+0.4330127018922194*(alphaz[5]*fin[18]+alphaz[0]*fin[17])+0.4330127018922193*alphaz[1]*fin[16]; 
  out[22] += 0.4330127018922194*alphaz[5]*fin[23]+0.4330127018922193*(alphaz[2]*fin[22]+alphaz[1]*fin[21])+0.4330127018922194*alphaz[0]*fin[19]; 
  out[23] += 0.4330127018922193*alphaz[2]*fin[23]+0.4330127018922194*(alphaz[5]*fin[22]+alphaz[0]*fin[21])+0.4330127018922193*alphaz[1]*fin[19]; 

  double alphavpar[24] = {0.}; 



  return 0.; 
} 
