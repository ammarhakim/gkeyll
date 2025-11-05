#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apar_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *b_i, const double *phi, const double *apar, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_i: covariant components of the field aligned unit vector.
  // apar: parallel component of magnetic vector potential.
  // phi: electrostatic potential .
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
  hamil[0] = 2.0*phi[0]*q_+vmapSq[0]*m_+1.4142135623730951*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.4142135623730951*bmag[1]*vmap[2]; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.4142135623730951*bmag[0]*vmap[3]; 
  hamil[5] = 1.4142135623730951*bmag[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[2]*hamil[2]; 
  hamil2[1] = hamil[8]*hamil[8]; 

  double alphax[12] = {0.}; 



  double alphavpar[12] = {0.}; 



  return 0.; 
} 
