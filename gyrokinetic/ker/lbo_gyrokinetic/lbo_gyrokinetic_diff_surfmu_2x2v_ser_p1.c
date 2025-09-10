#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_mapped_surfmu_2x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvell, const double *jacobvelc, const double *jacobvelr, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double vmap_primeSq = pow(vmap_prime[1],2);

  double fl_over_jacv[24], fc_over_jacv[24], fr_over_jacv[24];
  fl_over_jacv[0] = fl[0]/jacobvell[0]; 
  fl_over_jacv[1] = fl[1]/jacobvell[0]; 
  fl_over_jacv[2] = fl[2]/jacobvell[0]; 
  fl_over_jacv[3] = fl[3]/jacobvell[0]; 
  fl_over_jacv[4] = fl[4]/jacobvell[0]; 
  fl_over_jacv[5] = fl[5]/jacobvell[0]; 
  fl_over_jacv[6] = fl[6]/jacobvell[0]; 
  fl_over_jacv[7] = fl[7]/jacobvell[0]; 
  fl_over_jacv[8] = fl[8]/jacobvell[0]; 
  fl_over_jacv[9] = fl[9]/jacobvell[0]; 
  fl_over_jacv[10] = fl[10]/jacobvell[0]; 
  fl_over_jacv[11] = fl[11]/jacobvell[0]; 
  fl_over_jacv[12] = fl[12]/jacobvell[0]; 
  fl_over_jacv[13] = fl[13]/jacobvell[0]; 
  fl_over_jacv[14] = fl[14]/jacobvell[0]; 
  fl_over_jacv[15] = fl[15]/jacobvell[0]; 
  fl_over_jacv[16] = fl[16]/jacobvell[0]; 
  fl_over_jacv[17] = fl[17]/jacobvell[0]; 
  fl_over_jacv[18] = fl[18]/jacobvell[0]; 
  fl_over_jacv[19] = fl[19]/jacobvell[0]; 
  fl_over_jacv[20] = fl[20]/jacobvell[0]; 
  fl_over_jacv[21] = fl[21]/jacobvell[0]; 
  fl_over_jacv[22] = fl[22]/jacobvell[0]; 
  fl_over_jacv[23] = fl[23]/jacobvell[0]; 

  fc_over_jacv[0] = fc[0]/jacobvelc[0]; 
  fc_over_jacv[1] = fc[1]/jacobvelc[0]; 
  fc_over_jacv[2] = fc[2]/jacobvelc[0]; 
  fc_over_jacv[3] = fc[3]/jacobvelc[0]; 
  fc_over_jacv[4] = fc[4]/jacobvelc[0]; 
  fc_over_jacv[5] = fc[5]/jacobvelc[0]; 
  fc_over_jacv[6] = fc[6]/jacobvelc[0]; 
  fc_over_jacv[7] = fc[7]/jacobvelc[0]; 
  fc_over_jacv[8] = fc[8]/jacobvelc[0]; 
  fc_over_jacv[9] = fc[9]/jacobvelc[0]; 
  fc_over_jacv[10] = fc[10]/jacobvelc[0]; 
  fc_over_jacv[11] = fc[11]/jacobvelc[0]; 
  fc_over_jacv[12] = fc[12]/jacobvelc[0]; 
  fc_over_jacv[13] = fc[13]/jacobvelc[0]; 
  fc_over_jacv[14] = fc[14]/jacobvelc[0]; 
  fc_over_jacv[15] = fc[15]/jacobvelc[0]; 
  fc_over_jacv[16] = fc[16]/jacobvelc[0]; 
  fc_over_jacv[17] = fc[17]/jacobvelc[0]; 
  fc_over_jacv[18] = fc[18]/jacobvelc[0]; 
  fc_over_jacv[19] = fc[19]/jacobvelc[0]; 
  fc_over_jacv[20] = fc[20]/jacobvelc[0]; 
  fc_over_jacv[21] = fc[21]/jacobvelc[0]; 
  fc_over_jacv[22] = fc[22]/jacobvelc[0]; 
  fc_over_jacv[23] = fc[23]/jacobvelc[0]; 

  fr_over_jacv[0] = fr[0]/jacobvelr[0]; 
  fr_over_jacv[1] = fr[1]/jacobvelr[0]; 
  fr_over_jacv[2] = fr[2]/jacobvelr[0]; 
  fr_over_jacv[3] = fr[3]/jacobvelr[0]; 
  fr_over_jacv[4] = fr[4]/jacobvelr[0]; 
  fr_over_jacv[5] = fr[5]/jacobvelr[0]; 
  fr_over_jacv[6] = fr[6]/jacobvelr[0]; 
  fr_over_jacv[7] = fr[7]/jacobvelr[0]; 
  fr_over_jacv[8] = fr[8]/jacobvelr[0]; 
  fr_over_jacv[9] = fr[9]/jacobvelr[0]; 
  fr_over_jacv[10] = fr[10]/jacobvelr[0]; 
  fr_over_jacv[11] = fr[11]/jacobvelr[0]; 
  fr_over_jacv[12] = fr[12]/jacobvelr[0]; 
  fr_over_jacv[13] = fr[13]/jacobvelr[0]; 
  fr_over_jacv[14] = fr[14]/jacobvelr[0]; 
  fr_over_jacv[15] = fr[15]/jacobvelr[0]; 
  fr_over_jacv[16] = fr[16]/jacobvelr[0]; 
  fr_over_jacv[17] = fr[17]/jacobvelr[0]; 
  fr_over_jacv[18] = fr[18]/jacobvelr[0]; 
  fr_over_jacv[19] = fr[19]/jacobvelr[0]; 
  fr_over_jacv[20] = fr[20]/jacobvelr[0]; 
  fr_over_jacv[21] = fr[21]/jacobvelr[0]; 
  fr_over_jacv[22] = fr[22]/jacobvelr[0]; 
  fr_over_jacv[23] = fr[23]/jacobvelr[0]; 

  double dvl = 2.449489742783178*vmapl[3];
  double dvc = 2.449489742783178*vmapc[3];
  double dvr = 2.449489742783178*vmapr[3];

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdv2 = 2.0/dxv[3]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[4] = {0.}; 
  confFac[0] = bmag_inv[3]*nuVtSqSum[3]*m_+bmag_inv[2]*nuVtSqSum[2]*m_+bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = bmag_inv[2]*nuVtSqSum[3]*m_+nuVtSqSum[2]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*bmag_inv[1]*m_; 
  confFac[2] = bmag_inv[1]*nuVtSqSum[3]*m_+nuVtSqSum[1]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[2]*m_+nuVtSqSum[0]*bmag_inv[2]*m_; 
  confFac[3] = bmag_inv[0]*nuVtSqSum[3]*m_+nuVtSqSum[0]*bmag_inv[3]*m_+bmag_inv[1]*nuVtSqSum[2]*m_+nuVtSqSum[1]*bmag_inv[2]*m_; 

  double dfVfac_l = vmap_prime[0]*(0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]); 
  double dfVfac_r = vmap_prime[0]*(1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]); 

  double fVfac_l = (jacobvelc[0]*(0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]))/vmap_primeSq; 
  double fVfac_r = (jacobvelc[0]*(1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]))/vmap_primeSq; 

  double phaseFacl[24] = {0.0}; 
  const double dvlR2 = pow(dvl,2);
  const double dvlR3 = pow(dvl,3);
  const double dvlR4 = pow(dvl,4);
  const double dvcR2 = pow(dvc,2);
  const double dvcR3 = pow(dvc,3);
  const double dvcR4 = pow(dvc,4);

  phaseFacl[0] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[4]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[4]*dfVfac_l*rdv2-9.0*fl_over_jacv[0]*dfVfac_l*rdv2+9.0*fc_over_jacv[0]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[4]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[4]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[4]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[4]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[1] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[8]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[8]*dfVfac_l*rdv2-9.0*fl_over_jacv[1]*dfVfac_l*rdv2+9.0*fc_over_jacv[1]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[8]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[8]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[8]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[8]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[2] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[9]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[9]*dfVfac_l*rdv2-9.0*fl_over_jacv[2]*dfVfac_l*rdv2+9.0*fc_over_jacv[2]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[9]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[9]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[9]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[9]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[3] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[10]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[10]*dfVfac_l*rdv2-9.0*fl_over_jacv[3]*dfVfac_l*rdv2+9.0*fc_over_jacv[3]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[10]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[10]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[10]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[10]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[4] = (dvcR3*((-3.0*fl_over_jacv[4]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[4]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[4]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[4]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[4]*dfVfac_l*rdv2+15.0*fc_over_jacv[4]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[0]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[0]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[4]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[4]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[4]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[4]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[5] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[12]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[12]*dfVfac_l*rdv2-9.0*fl_over_jacv[5]*dfVfac_l*rdv2+9.0*fc_over_jacv[5]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[12]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[12]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[12]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[12]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[6] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[13]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[13]*dfVfac_l*rdv2-9.0*fl_over_jacv[6]*dfVfac_l*rdv2+9.0*fc_over_jacv[6]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[13]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[13]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[13]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[13]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[7] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[14]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[14]*dfVfac_l*rdv2-9.0*fl_over_jacv[7]*dfVfac_l*rdv2+9.0*fc_over_jacv[7]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[14]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[14]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[14]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[14]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[8] = (dvcR3*((-3.0*fl_over_jacv[8]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[8]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[8]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[8]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[8]*dfVfac_l*rdv2+15.0*fc_over_jacv[8]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[1]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[1]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[8]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[8]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[8]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[8]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[9] = (dvcR3*((-3.0*fl_over_jacv[9]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[9]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[9]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[9]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[9]*dfVfac_l*rdv2+15.0*fc_over_jacv[9]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[2]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[2]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[9]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[9]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[9]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[9]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[10] = (dvcR3*((-3.0*fl_over_jacv[10]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[3]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[10]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[3]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[10]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[3]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[10]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[3]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[10]*dfVfac_l*rdv2+15.0*fc_over_jacv[10]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[3]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[3]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[10]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[10]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[10]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[10]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[11] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[15]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[15]*dfVfac_l*rdv2-9.0*fl_over_jacv[11]*dfVfac_l*rdv2+9.0*fc_over_jacv[11]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[15]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[15]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[15]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[15]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[12] = (dvcR3*((-3.0*fl_over_jacv[12]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[5]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[12]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[5]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[12]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[5]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[12]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[5]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[12]*dfVfac_l*rdv2+15.0*fc_over_jacv[12]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[5]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[5]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[12]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[12]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[12]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[12]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[13] = (dvcR3*((-3.0*fl_over_jacv[13]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[6]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[13]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[6]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[13]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[6]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[13]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[6]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[13]*dfVfac_l*rdv2+15.0*fc_over_jacv[13]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[6]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[6]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[13]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[13]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[13]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[13]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[14] = (dvcR3*((-3.0*fl_over_jacv[14]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[7]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[14]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[7]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[14]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[7]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[14]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[7]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[14]*dfVfac_l*rdv2+15.0*fc_over_jacv[14]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[7]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[7]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[14]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[14]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[14]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[14]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[15] = (dvcR3*((-3.0*fl_over_jacv[15]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[11]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[15]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[11]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[15]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[11]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[15]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[11]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[15]*dfVfac_l*rdv2+15.0*fc_over_jacv[15]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[11]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[11]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[15]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[15]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[15]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[15]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[16] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[19]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[19]*dfVfac_l*rdv2-45.0*fl_over_jacv[16]*dfVfac_l*rdv2+45.0*fc_over_jacv[16]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[19]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[19]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[19]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[19]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[17] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[21]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[21]*dfVfac_l*rdv2-45.0*fl_over_jacv[17]*dfVfac_l*rdv2+45.0*fc_over_jacv[17]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[21]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[21]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[21]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[21]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[18] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[22]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[22]*dfVfac_l*rdv2-45.0*fl_over_jacv[18]*dfVfac_l*rdv2+45.0*fc_over_jacv[18]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[22]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[22]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[22]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[22]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[19] = (dvcR3*((-15.0*fl_over_jacv[19]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[16]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[19]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[16]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[19]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[16]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[19]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[16]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[19]*dfVfac_l*rdv2+75.0*fc_over_jacv[19]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[16]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[16]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[19]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[19]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[19]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[19]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[20] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[23]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[23]*dfVfac_l*rdv2-45.0*fl_over_jacv[20]*dfVfac_l*rdv2+45.0*fc_over_jacv[20]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[23]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[23]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[23]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[23]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[21] = (dvcR3*((-15.0*fl_over_jacv[21]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[17]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[21]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[17]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[21]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[17]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[21]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[17]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[21]*dfVfac_l*rdv2+75.0*fc_over_jacv[21]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[17]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[17]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[21]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[21]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[21]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[21]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[22] = (dvcR3*((-15.0*fl_over_jacv[22]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[18]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[22]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[18]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[22]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[18]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[22]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[18]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[22]*dfVfac_l*rdv2+75.0*fc_over_jacv[22]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[18]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[18]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[22]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[22]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[22]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[22]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[23] = (dvcR3*((-15.0*fl_over_jacv[23]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[20]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[23]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[20]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[23]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[20]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[23]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[20]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[23]*dfVfac_l*rdv2+75.0*fc_over_jacv[23]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[20]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[20]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[23]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[23]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[23]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[23]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 

  double phaseFacr[24] = {0.0}; 
  const double dvrR2 = pow(dvr,2);
  const double dvrR3 = pow(dvr,3);
  const double dvrR4 = pow(dvr,4);

  phaseFacr[0] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[4]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[4]*dfVfac_r*rdv2+9.0*fr_over_jacv[0]*dfVfac_r*rdv2-9.0*fc_over_jacv[0]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[4]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[4]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[4]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[4]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[1] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[8]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[8]*dfVfac_r*rdv2+9.0*fr_over_jacv[1]*dfVfac_r*rdv2-9.0*fc_over_jacv[1]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[8]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[8]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[8]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[8]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[2] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[9]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[9]*dfVfac_r*rdv2+9.0*fr_over_jacv[2]*dfVfac_r*rdv2-9.0*fc_over_jacv[2]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[9]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[9]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[9]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[9]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[3] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[10]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[10]*dfVfac_r*rdv2+9.0*fr_over_jacv[3]*dfVfac_r*rdv2-9.0*fc_over_jacv[3]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[10]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[10]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[10]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[10]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[4] = (dvcR2*dvr*(5.0*fr_over_jacv[4]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[4]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[4]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[4]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[4]*dfVfac_r*rdv2)-15.0*fc_over_jacv[4]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[0]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[0]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[4]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[4]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[4]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[4]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[5] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[12]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[12]*dfVfac_r*rdv2+9.0*fr_over_jacv[5]*dfVfac_r*rdv2-9.0*fc_over_jacv[5]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[12]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[12]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[12]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[12]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[6] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[13]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[13]*dfVfac_r*rdv2+9.0*fr_over_jacv[6]*dfVfac_r*rdv2-9.0*fc_over_jacv[6]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[13]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[13]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[13]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[13]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[7] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[14]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[14]*dfVfac_r*rdv2+9.0*fr_over_jacv[7]*dfVfac_r*rdv2-9.0*fc_over_jacv[7]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[14]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[14]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[14]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[14]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[8] = (dvcR2*dvr*(5.0*fr_over_jacv[8]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[8]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[8]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[8]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[8]*dfVfac_r*rdv2)-15.0*fc_over_jacv[8]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[1]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[1]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[8]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[8]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[8]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[8]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[9] = (dvcR2*dvr*(5.0*fr_over_jacv[9]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[9]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[9]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[9]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[9]*dfVfac_r*rdv2)-15.0*fc_over_jacv[9]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[2]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[2]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[9]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[9]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[9]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[9]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[10] = (dvcR2*dvr*(5.0*fr_over_jacv[10]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[3]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[10]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[3]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[10]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[3]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[10]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[3]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[10]*dfVfac_r*rdv2)-15.0*fc_over_jacv[10]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[3]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[3]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[10]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[10]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[10]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[10]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[11] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[15]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[15]*dfVfac_r*rdv2+9.0*fr_over_jacv[11]*dfVfac_r*rdv2-9.0*fc_over_jacv[11]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[15]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[15]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[15]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[15]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[12] = (dvcR2*dvr*(5.0*fr_over_jacv[12]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[5]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[12]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[5]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[12]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[5]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[12]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[5]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[12]*dfVfac_r*rdv2)-15.0*fc_over_jacv[12]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[5]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[5]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[12]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[12]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[12]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[12]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[13] = (dvcR2*dvr*(5.0*fr_over_jacv[13]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[6]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[13]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[6]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[13]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[6]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[13]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[6]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[13]*dfVfac_r*rdv2)-15.0*fc_over_jacv[13]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[6]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[6]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[13]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[13]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[13]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[13]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[14] = (dvcR2*dvr*(5.0*fr_over_jacv[14]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[7]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[14]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[7]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[14]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[7]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[14]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[7]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[14]*dfVfac_r*rdv2)-15.0*fc_over_jacv[14]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[7]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[7]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[14]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[14]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[14]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[14]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[15] = (dvcR2*dvr*(5.0*fr_over_jacv[15]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[11]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[15]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[11]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[15]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[11]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[15]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[11]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[15]*dfVfac_r*rdv2)-15.0*fc_over_jacv[15]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[11]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[11]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[15]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[15]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[15]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[15]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[16] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[19]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[19]*dfVfac_r*rdv2+45.0*fr_over_jacv[16]*dfVfac_r*rdv2-45.0*fc_over_jacv[16]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[19]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[19]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[19]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[19]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[17] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[21]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[21]*dfVfac_r*rdv2+45.0*fr_over_jacv[17]*dfVfac_r*rdv2-45.0*fc_over_jacv[17]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[21]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[21]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[21]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[21]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[18] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[22]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[22]*dfVfac_r*rdv2+45.0*fr_over_jacv[18]*dfVfac_r*rdv2-45.0*fc_over_jacv[18]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[22]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[22]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[22]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[22]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[19] = (dvcR2*dvr*(25.0*fr_over_jacv[19]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[16]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[19]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[16]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[19]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[16]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[19]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[16]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[19]*dfVfac_r*rdv2)-75.0*fc_over_jacv[19]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[16]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[16]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[19]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[19]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[19]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[19]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[20] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[23]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[23]*dfVfac_r*rdv2+45.0*fr_over_jacv[20]*dfVfac_r*rdv2-45.0*fc_over_jacv[20]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[23]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[23]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[23]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[23]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[21] = (dvcR2*dvr*(25.0*fr_over_jacv[21]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[17]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[21]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[17]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[21]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[17]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[21]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[17]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[21]*dfVfac_r*rdv2)-75.0*fc_over_jacv[21]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[17]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[17]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[21]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[21]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[21]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[21]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[22] = (dvcR2*dvr*(25.0*fr_over_jacv[22]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[18]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[22]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[18]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[22]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[18]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[22]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[18]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[22]*dfVfac_r*rdv2)-75.0*fc_over_jacv[22]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[18]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[18]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[22]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[22]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[22]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[22]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[23] = (dvcR2*dvr*(25.0*fr_over_jacv[23]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[20]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[23]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[20]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[23]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[20]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[23]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[20]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[23]*dfVfac_r*rdv2)-75.0*fc_over_jacv[23]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[20]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[20]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[23]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[23]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[23]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[23]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 

  double incrl[24] = {0.0}; 
  incrl[0] = 0.5*confFac[3]*phaseFacl[5]+0.5*confFac[2]*phaseFacl[2]+0.5*confFac[1]*phaseFacl[1]+0.5*confFac[0]*phaseFacl[0]; 
  incrl[1] = 0.5*confFac[2]*phaseFacl[5]+0.5*phaseFacl[2]*confFac[3]+0.5*confFac[0]*phaseFacl[1]+0.5*phaseFacl[0]*confFac[1]; 
  incrl[2] = 0.5*confFac[1]*phaseFacl[5]+0.5*phaseFacl[1]*confFac[3]+0.5*confFac[0]*phaseFacl[2]+0.5*phaseFacl[0]*confFac[2]; 
  incrl[3] = 0.5*confFac[3]*phaseFacl[11]+0.5*confFac[2]*phaseFacl[7]+0.5*confFac[1]*phaseFacl[6]+0.5*confFac[0]*phaseFacl[3]; 
  incrl[4] = 0.5*confFac[3]*phaseFacl[12]+0.5*confFac[2]*phaseFacl[9]+0.5*confFac[1]*phaseFacl[8]+0.5*confFac[0]*phaseFacl[4]; 
  incrl[5] = 0.5*confFac[0]*phaseFacl[5]+0.5*phaseFacl[0]*confFac[3]+0.5*confFac[1]*phaseFacl[2]+0.5*phaseFacl[1]*confFac[2]; 
  incrl[6] = 0.5*confFac[2]*phaseFacl[11]+0.5*confFac[3]*phaseFacl[7]+0.5*confFac[0]*phaseFacl[6]+0.5*confFac[1]*phaseFacl[3]; 
  incrl[7] = 0.5*confFac[1]*phaseFacl[11]+0.5*confFac[0]*phaseFacl[7]+0.5*confFac[3]*phaseFacl[6]+0.5*confFac[2]*phaseFacl[3]; 
  incrl[8] = 0.5*confFac[2]*phaseFacl[12]+0.5*confFac[3]*phaseFacl[9]+0.5*confFac[0]*phaseFacl[8]+0.5*confFac[1]*phaseFacl[4]; 
  incrl[9] = 0.5*confFac[1]*phaseFacl[12]+0.5*confFac[0]*phaseFacl[9]+0.5*confFac[3]*phaseFacl[8]+0.5*confFac[2]*phaseFacl[4]; 
  incrl[10] = 0.5*confFac[3]*phaseFacl[15]+0.5*confFac[2]*phaseFacl[14]+0.5*confFac[1]*phaseFacl[13]+0.5*confFac[0]*phaseFacl[10]; 
  incrl[11] = 0.5*confFac[0]*phaseFacl[11]+0.5*confFac[1]*phaseFacl[7]+0.5*confFac[2]*phaseFacl[6]+0.5*confFac[3]*phaseFacl[3]; 
  incrl[12] = 0.5*confFac[0]*phaseFacl[12]+0.5*confFac[1]*phaseFacl[9]+0.5*confFac[2]*phaseFacl[8]+0.5*confFac[3]*phaseFacl[4]; 
  incrl[13] = 0.5*confFac[2]*phaseFacl[15]+0.5*confFac[3]*phaseFacl[14]+0.5*confFac[0]*phaseFacl[13]+0.5*confFac[1]*phaseFacl[10]; 
  incrl[14] = 0.5*confFac[1]*phaseFacl[15]+0.5*confFac[0]*phaseFacl[14]+0.5*confFac[3]*phaseFacl[13]+0.5*confFac[2]*phaseFacl[10]; 
  incrl[15] = 0.5*confFac[0]*phaseFacl[15]+0.5*confFac[1]*phaseFacl[14]+0.5*confFac[2]*phaseFacl[13]+0.5*confFac[3]*phaseFacl[10]; 
  incrl[16] = 0.5*confFac[3]*phaseFacl[20]+0.5000000000000001*confFac[2]*phaseFacl[18]+0.5000000000000001*confFac[1]*phaseFacl[17]+0.5*confFac[0]*phaseFacl[16]; 
  incrl[17] = 0.5000000000000001*confFac[2]*phaseFacl[20]+0.5*confFac[3]*phaseFacl[18]+0.5*confFac[0]*phaseFacl[17]+0.5000000000000001*confFac[1]*phaseFacl[16]; 
  incrl[18] = 0.5000000000000001*confFac[1]*phaseFacl[20]+0.5*confFac[0]*phaseFacl[18]+0.5*confFac[3]*phaseFacl[17]+0.5000000000000001*confFac[2]*phaseFacl[16]; 
  incrl[19] = 0.5*confFac[3]*phaseFacl[23]+0.5000000000000001*confFac[2]*phaseFacl[22]+0.5000000000000001*confFac[1]*phaseFacl[21]+0.5*confFac[0]*phaseFacl[19]; 
  incrl[20] = 0.5*confFac[0]*phaseFacl[20]+0.5000000000000001*confFac[1]*phaseFacl[18]+0.5000000000000001*confFac[2]*phaseFacl[17]+0.5*confFac[3]*phaseFacl[16]; 
  incrl[21] = 0.5000000000000001*confFac[2]*phaseFacl[23]+0.5*confFac[3]*phaseFacl[22]+0.5*confFac[0]*phaseFacl[21]+0.5000000000000001*confFac[1]*phaseFacl[19]; 
  incrl[22] = 0.5000000000000001*confFac[1]*phaseFacl[23]+0.5*confFac[0]*phaseFacl[22]+0.5*confFac[3]*phaseFacl[21]+0.5000000000000001*confFac[2]*phaseFacl[19]; 
  incrl[23] = 0.5*confFac[0]*phaseFacl[23]+0.5000000000000001*confFac[1]*phaseFacl[22]+0.5000000000000001*confFac[2]*phaseFacl[21]+0.5*confFac[3]*phaseFacl[19]; 

  double incrr[24] = {0.0}; 
  incrr[0] = 0.5*confFac[3]*phaseFacr[5]+0.5*confFac[2]*phaseFacr[2]+0.5*confFac[1]*phaseFacr[1]+0.5*confFac[0]*phaseFacr[0]; 
  incrr[1] = 0.5*confFac[2]*phaseFacr[5]+0.5*phaseFacr[2]*confFac[3]+0.5*confFac[0]*phaseFacr[1]+0.5*phaseFacr[0]*confFac[1]; 
  incrr[2] = 0.5*confFac[1]*phaseFacr[5]+0.5*phaseFacr[1]*confFac[3]+0.5*confFac[0]*phaseFacr[2]+0.5*phaseFacr[0]*confFac[2]; 
  incrr[3] = 0.5*confFac[3]*phaseFacr[11]+0.5*confFac[2]*phaseFacr[7]+0.5*confFac[1]*phaseFacr[6]+0.5*confFac[0]*phaseFacr[3]; 
  incrr[4] = 0.5*confFac[3]*phaseFacr[12]+0.5*confFac[2]*phaseFacr[9]+0.5*confFac[1]*phaseFacr[8]+0.5*confFac[0]*phaseFacr[4]; 
  incrr[5] = 0.5*confFac[0]*phaseFacr[5]+0.5*phaseFacr[0]*confFac[3]+0.5*confFac[1]*phaseFacr[2]+0.5*phaseFacr[1]*confFac[2]; 
  incrr[6] = 0.5*confFac[2]*phaseFacr[11]+0.5*confFac[3]*phaseFacr[7]+0.5*confFac[0]*phaseFacr[6]+0.5*confFac[1]*phaseFacr[3]; 
  incrr[7] = 0.5*confFac[1]*phaseFacr[11]+0.5*confFac[0]*phaseFacr[7]+0.5*confFac[3]*phaseFacr[6]+0.5*confFac[2]*phaseFacr[3]; 
  incrr[8] = 0.5*confFac[2]*phaseFacr[12]+0.5*confFac[3]*phaseFacr[9]+0.5*confFac[0]*phaseFacr[8]+0.5*confFac[1]*phaseFacr[4]; 
  incrr[9] = 0.5*confFac[1]*phaseFacr[12]+0.5*confFac[0]*phaseFacr[9]+0.5*confFac[3]*phaseFacr[8]+0.5*confFac[2]*phaseFacr[4]; 
  incrr[10] = 0.5*confFac[3]*phaseFacr[15]+0.5*confFac[2]*phaseFacr[14]+0.5*confFac[1]*phaseFacr[13]+0.5*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.5*confFac[0]*phaseFacr[11]+0.5*confFac[1]*phaseFacr[7]+0.5*confFac[2]*phaseFacr[6]+0.5*confFac[3]*phaseFacr[3]; 
  incrr[12] = 0.5*confFac[0]*phaseFacr[12]+0.5*confFac[1]*phaseFacr[9]+0.5*confFac[2]*phaseFacr[8]+0.5*confFac[3]*phaseFacr[4]; 
  incrr[13] = 0.5*confFac[2]*phaseFacr[15]+0.5*confFac[3]*phaseFacr[14]+0.5*confFac[0]*phaseFacr[13]+0.5*confFac[1]*phaseFacr[10]; 
  incrr[14] = 0.5*confFac[1]*phaseFacr[15]+0.5*confFac[0]*phaseFacr[14]+0.5*confFac[3]*phaseFacr[13]+0.5*confFac[2]*phaseFacr[10]; 
  incrr[15] = 0.5*confFac[0]*phaseFacr[15]+0.5*confFac[1]*phaseFacr[14]+0.5*confFac[2]*phaseFacr[13]+0.5*confFac[3]*phaseFacr[10]; 
  incrr[16] = 0.5*confFac[3]*phaseFacr[20]+0.5000000000000001*confFac[2]*phaseFacr[18]+0.5000000000000001*confFac[1]*phaseFacr[17]+0.5*confFac[0]*phaseFacr[16]; 
  incrr[17] = 0.5000000000000001*confFac[2]*phaseFacr[20]+0.5*confFac[3]*phaseFacr[18]+0.5*confFac[0]*phaseFacr[17]+0.5000000000000001*confFac[1]*phaseFacr[16]; 
  incrr[18] = 0.5000000000000001*confFac[1]*phaseFacr[20]+0.5*confFac[0]*phaseFacr[18]+0.5*confFac[3]*phaseFacr[17]+0.5000000000000001*confFac[2]*phaseFacr[16]; 
  incrr[19] = 0.5*confFac[3]*phaseFacr[23]+0.5000000000000001*confFac[2]*phaseFacr[22]+0.5000000000000001*confFac[1]*phaseFacr[21]+0.5*confFac[0]*phaseFacr[19]; 
  incrr[20] = 0.5*confFac[0]*phaseFacr[20]+0.5000000000000001*confFac[1]*phaseFacr[18]+0.5000000000000001*confFac[2]*phaseFacr[17]+0.5*confFac[3]*phaseFacr[16]; 
  incrr[21] = 0.5000000000000001*confFac[2]*phaseFacr[23]+0.5*confFac[3]*phaseFacr[22]+0.5*confFac[0]*phaseFacr[21]+0.5000000000000001*confFac[1]*phaseFacr[19]; 
  incrr[22] = 0.5000000000000001*confFac[1]*phaseFacr[23]+0.5*confFac[0]*phaseFacr[22]+0.5*confFac[3]*phaseFacr[21]+0.5000000000000001*confFac[2]*phaseFacr[19]; 
  incrr[23] = 0.5*confFac[0]*phaseFacr[23]+0.5000000000000001*confFac[1]*phaseFacr[22]+0.5000000000000001*confFac[2]*phaseFacr[21]+0.5*confFac[3]*phaseFacr[19]; 

  out[0] += incrr[0]-1.0*incrl[0]; 
  out[1] += incrr[1]-1.0*incrl[1]; 
  out[2] += incrr[2]-1.0*incrl[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += incrr[4]-1.0*incrl[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += incrr[8]-1.0*incrl[8]; 
  out[9] += incrr[9]-1.0*incrl[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 
  out[12] += incrr[12]-1.0*incrl[12]; 
  out[13] += incrr[13]-1.0*incrl[13]; 
  out[14] += incrr[14]-1.0*incrl[14]; 
  out[15] += incrr[15]-1.0*incrl[15]; 
  out[16] += incrr[16]-1.0*incrl[16]; 
  out[17] += incrr[17]-1.0*incrl[17]; 
  out[18] += incrr[18]-1.0*incrl[18]; 
  out[19] += incrr[19]-1.0*incrl[19]; 
  out[20] += incrr[20]-1.0*incrl[20]; 
  out[21] += incrr[21]-1.0*incrl[21]; 
  out[22] += incrr[22]-1.0*incrl[22]; 
  out[23] += incrr[23]-1.0*incrl[23]; 

  return 0.;

} 
GKYL_CU_DH double lbo_gyrokinetic_diff_notmapped_surfmu_2x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvell, const double *jacobvelc, const double *jacobvelr, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdv2 = 2.0/dxv[3]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[4] = {0.}; 
  confFac[0] = bmag_inv[3]*nuVtSqSum[3]*m_+bmag_inv[2]*nuVtSqSum[2]*m_+bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = bmag_inv[2]*nuVtSqSum[3]*m_+nuVtSqSum[2]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*bmag_inv[1]*m_; 
  confFac[2] = bmag_inv[1]*nuVtSqSum[3]*m_+nuVtSqSum[1]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[2]*m_+nuVtSqSum[0]*bmag_inv[2]*m_; 
  confFac[3] = bmag_inv[0]*nuVtSqSum[3]*m_+nuVtSqSum[0]*bmag_inv[3]*m_+bmag_inv[1]*nuVtSqSum[2]*m_+nuVtSqSum[1]*bmag_inv[2]*m_; 

  double dfVfac_l = 0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]; 
  double dfVfac_r = 1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]; 

  double fVfac_l = 0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]; 
  double fVfac_r = 1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]; 

  double phaseFacl[24] = {0.0}; 

  phaseFacl[0] = (-0.5412658773652741*fl[4]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[4]*dfVfac_l*rdv2Sq-0.5625*fl[0]*dfVfac_l*rdv2Sq+0.5625*fc[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[1] = (-0.5412658773652741*fl[8]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[8]*dfVfac_l*rdv2Sq-0.5625*fl[1]*dfVfac_l*rdv2Sq+0.5625*fc[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[2] = (-0.5412658773652741*fl[9]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[9]*dfVfac_l*rdv2Sq-0.5625*fl[2]*dfVfac_l*rdv2Sq+0.5625*fc[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[3] = (-0.5412658773652741*fl[10]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[10]*dfVfac_l*rdv2Sq-0.5625*fl[3]*dfVfac_l*rdv2Sq+0.5625*fc[3]*dfVfac_l*rdv2Sq; 
  phaseFacl[4] = (-0.5*fl[4]*fVfac_l*rdv2Sq)+0.5*fc[4]*fVfac_l*rdv2Sq-0.4330127018922193*fl[0]*fVfac_l*rdv2Sq-0.4330127018922193*fc[0]*fVfac_l*rdv2Sq+0.9375*fl[4]*dfVfac_l*rdv2Sq+0.9375*fc[4]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[0]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[5] = (-0.5412658773652741*fl[12]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[12]*dfVfac_l*rdv2Sq-0.5625*fl[5]*dfVfac_l*rdv2Sq+0.5625*fc[5]*dfVfac_l*rdv2Sq; 
  phaseFacl[6] = (-0.5412658773652741*fl[13]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[13]*dfVfac_l*rdv2Sq-0.5625*fl[6]*dfVfac_l*rdv2Sq+0.5625*fc[6]*dfVfac_l*rdv2Sq; 
  phaseFacl[7] = (-0.5412658773652741*fl[14]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[14]*dfVfac_l*rdv2Sq-0.5625*fl[7]*dfVfac_l*rdv2Sq+0.5625*fc[7]*dfVfac_l*rdv2Sq; 
  phaseFacl[8] = (-0.5*fl[8]*fVfac_l*rdv2Sq)+0.5*fc[8]*fVfac_l*rdv2Sq-0.4330127018922193*fl[1]*fVfac_l*rdv2Sq-0.4330127018922193*fc[1]*fVfac_l*rdv2Sq+0.9375*fl[8]*dfVfac_l*rdv2Sq+0.9375*fc[8]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[1]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[9] = (-0.5*fl[9]*fVfac_l*rdv2Sq)+0.5*fc[9]*fVfac_l*rdv2Sq-0.4330127018922193*fl[2]*fVfac_l*rdv2Sq-0.4330127018922193*fc[2]*fVfac_l*rdv2Sq+0.9375*fl[9]*dfVfac_l*rdv2Sq+0.9375*fc[9]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[2]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[10] = (-0.5*fl[10]*fVfac_l*rdv2Sq)+0.5*fc[10]*fVfac_l*rdv2Sq-0.4330127018922193*fl[3]*fVfac_l*rdv2Sq-0.4330127018922193*fc[3]*fVfac_l*rdv2Sq+0.9375*fl[10]*dfVfac_l*rdv2Sq+0.9375*fc[10]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[3]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[3]*dfVfac_l*rdv2Sq; 
  phaseFacl[11] = (-0.5412658773652741*fl[15]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[15]*dfVfac_l*rdv2Sq-0.5625*fl[11]*dfVfac_l*rdv2Sq+0.5625*fc[11]*dfVfac_l*rdv2Sq; 
  phaseFacl[12] = (-0.5*fl[12]*fVfac_l*rdv2Sq)+0.5*fc[12]*fVfac_l*rdv2Sq-0.4330127018922193*fl[5]*fVfac_l*rdv2Sq-0.4330127018922193*fc[5]*fVfac_l*rdv2Sq+0.9375*fl[12]*dfVfac_l*rdv2Sq+0.9375*fc[12]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[5]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[5]*dfVfac_l*rdv2Sq; 
  phaseFacl[13] = (-0.5*fl[13]*fVfac_l*rdv2Sq)+0.5*fc[13]*fVfac_l*rdv2Sq-0.4330127018922193*fl[6]*fVfac_l*rdv2Sq-0.4330127018922193*fc[6]*fVfac_l*rdv2Sq+0.9375*fl[13]*dfVfac_l*rdv2Sq+0.9375*fc[13]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[6]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[6]*dfVfac_l*rdv2Sq; 
  phaseFacl[14] = (-0.5*fl[14]*fVfac_l*rdv2Sq)+0.5*fc[14]*fVfac_l*rdv2Sq-0.4330127018922193*fl[7]*fVfac_l*rdv2Sq-0.4330127018922193*fc[7]*fVfac_l*rdv2Sq+0.9375*fl[14]*dfVfac_l*rdv2Sq+0.9375*fc[14]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[7]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[7]*dfVfac_l*rdv2Sq; 
  phaseFacl[15] = (-0.5*fl[15]*fVfac_l*rdv2Sq)+0.5*fc[15]*fVfac_l*rdv2Sq-0.4330127018922193*fl[11]*fVfac_l*rdv2Sq-0.4330127018922193*fc[11]*fVfac_l*rdv2Sq+0.9375*fl[15]*dfVfac_l*rdv2Sq+0.9375*fc[15]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[11]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[11]*dfVfac_l*rdv2Sq; 
  phaseFacl[16] = (-0.5412658773652742*fl[19]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[19]*dfVfac_l*rdv2Sq-0.5625*fl[16]*dfVfac_l*rdv2Sq+0.5625*fc[16]*dfVfac_l*rdv2Sq; 
  phaseFacl[17] = (-0.5412658773652742*fl[21]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[21]*dfVfac_l*rdv2Sq-0.5625*fl[17]*dfVfac_l*rdv2Sq+0.5625*fc[17]*dfVfac_l*rdv2Sq; 
  phaseFacl[18] = (-0.5412658773652742*fl[22]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[22]*dfVfac_l*rdv2Sq-0.5625*fl[18]*dfVfac_l*rdv2Sq+0.5625*fc[18]*dfVfac_l*rdv2Sq; 
  phaseFacl[19] = (-0.5*fl[19]*fVfac_l*rdv2Sq)+0.5*fc[19]*fVfac_l*rdv2Sq-0.4330127018922194*fl[16]*fVfac_l*rdv2Sq-0.4330127018922194*fc[16]*fVfac_l*rdv2Sq+0.9375*fl[19]*dfVfac_l*rdv2Sq+0.9375*fc[19]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[16]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[16]*dfVfac_l*rdv2Sq; 
  phaseFacl[20] = (-0.5412658773652742*fl[23]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[23]*dfVfac_l*rdv2Sq-0.5625*fl[20]*dfVfac_l*rdv2Sq+0.5625*fc[20]*dfVfac_l*rdv2Sq; 
  phaseFacl[21] = (-0.5*fl[21]*fVfac_l*rdv2Sq)+0.5*fc[21]*fVfac_l*rdv2Sq-0.4330127018922194*fl[17]*fVfac_l*rdv2Sq-0.4330127018922194*fc[17]*fVfac_l*rdv2Sq+0.9375*fl[21]*dfVfac_l*rdv2Sq+0.9375*fc[21]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[17]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[17]*dfVfac_l*rdv2Sq; 
  phaseFacl[22] = (-0.5*fl[22]*fVfac_l*rdv2Sq)+0.5*fc[22]*fVfac_l*rdv2Sq-0.4330127018922194*fl[18]*fVfac_l*rdv2Sq-0.4330127018922194*fc[18]*fVfac_l*rdv2Sq+0.9375*fl[22]*dfVfac_l*rdv2Sq+0.9375*fc[22]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[18]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[18]*dfVfac_l*rdv2Sq; 
  phaseFacl[23] = (-0.5*fl[23]*fVfac_l*rdv2Sq)+0.5*fc[23]*fVfac_l*rdv2Sq-0.4330127018922194*fl[20]*fVfac_l*rdv2Sq-0.4330127018922194*fc[20]*fVfac_l*rdv2Sq+0.9375*fl[23]*dfVfac_l*rdv2Sq+0.9375*fc[23]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[20]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[20]*dfVfac_l*rdv2Sq; 

  double phaseFacr[24] = {0.0}; 

  phaseFacr[0] = (-0.5412658773652741*fr[4]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[4]*dfVfac_r*rdv2Sq+0.5625*fr[0]*dfVfac_r*rdv2Sq-0.5625*fc[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[1] = (-0.5412658773652741*fr[8]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[8]*dfVfac_r*rdv2Sq+0.5625*fr[1]*dfVfac_r*rdv2Sq-0.5625*fc[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[2] = (-0.5412658773652741*fr[9]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[9]*dfVfac_r*rdv2Sq+0.5625*fr[2]*dfVfac_r*rdv2Sq-0.5625*fc[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[3] = (-0.5412658773652741*fr[10]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[10]*dfVfac_r*rdv2Sq+0.5625*fr[3]*dfVfac_r*rdv2Sq-0.5625*fc[3]*dfVfac_r*rdv2Sq; 
  phaseFacr[4] = 0.5*fr[4]*fVfac_r*rdv2Sq-0.5*fc[4]*fVfac_r*rdv2Sq-0.4330127018922193*fr[0]*fVfac_r*rdv2Sq-0.4330127018922193*fc[0]*fVfac_r*rdv2Sq-0.9375*fr[4]*dfVfac_r*rdv2Sq-0.9375*fc[4]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[0]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[5] = (-0.5412658773652741*fr[12]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[12]*dfVfac_r*rdv2Sq+0.5625*fr[5]*dfVfac_r*rdv2Sq-0.5625*fc[5]*dfVfac_r*rdv2Sq; 
  phaseFacr[6] = (-0.5412658773652741*fr[13]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[13]*dfVfac_r*rdv2Sq+0.5625*fr[6]*dfVfac_r*rdv2Sq-0.5625*fc[6]*dfVfac_r*rdv2Sq; 
  phaseFacr[7] = (-0.5412658773652741*fr[14]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[14]*dfVfac_r*rdv2Sq+0.5625*fr[7]*dfVfac_r*rdv2Sq-0.5625*fc[7]*dfVfac_r*rdv2Sq; 
  phaseFacr[8] = 0.5*fr[8]*fVfac_r*rdv2Sq-0.5*fc[8]*fVfac_r*rdv2Sq-0.4330127018922193*fr[1]*fVfac_r*rdv2Sq-0.4330127018922193*fc[1]*fVfac_r*rdv2Sq-0.9375*fr[8]*dfVfac_r*rdv2Sq-0.9375*fc[8]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[1]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[9] = 0.5*fr[9]*fVfac_r*rdv2Sq-0.5*fc[9]*fVfac_r*rdv2Sq-0.4330127018922193*fr[2]*fVfac_r*rdv2Sq-0.4330127018922193*fc[2]*fVfac_r*rdv2Sq-0.9375*fr[9]*dfVfac_r*rdv2Sq-0.9375*fc[9]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[2]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[10] = 0.5*fr[10]*fVfac_r*rdv2Sq-0.5*fc[10]*fVfac_r*rdv2Sq-0.4330127018922193*fr[3]*fVfac_r*rdv2Sq-0.4330127018922193*fc[3]*fVfac_r*rdv2Sq-0.9375*fr[10]*dfVfac_r*rdv2Sq-0.9375*fc[10]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[3]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[3]*dfVfac_r*rdv2Sq; 
  phaseFacr[11] = (-0.5412658773652741*fr[15]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[15]*dfVfac_r*rdv2Sq+0.5625*fr[11]*dfVfac_r*rdv2Sq-0.5625*fc[11]*dfVfac_r*rdv2Sq; 
  phaseFacr[12] = 0.5*fr[12]*fVfac_r*rdv2Sq-0.5*fc[12]*fVfac_r*rdv2Sq-0.4330127018922193*fr[5]*fVfac_r*rdv2Sq-0.4330127018922193*fc[5]*fVfac_r*rdv2Sq-0.9375*fr[12]*dfVfac_r*rdv2Sq-0.9375*fc[12]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[5]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[5]*dfVfac_r*rdv2Sq; 
  phaseFacr[13] = 0.5*fr[13]*fVfac_r*rdv2Sq-0.5*fc[13]*fVfac_r*rdv2Sq-0.4330127018922193*fr[6]*fVfac_r*rdv2Sq-0.4330127018922193*fc[6]*fVfac_r*rdv2Sq-0.9375*fr[13]*dfVfac_r*rdv2Sq-0.9375*fc[13]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[6]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[6]*dfVfac_r*rdv2Sq; 
  phaseFacr[14] = 0.5*fr[14]*fVfac_r*rdv2Sq-0.5*fc[14]*fVfac_r*rdv2Sq-0.4330127018922193*fr[7]*fVfac_r*rdv2Sq-0.4330127018922193*fc[7]*fVfac_r*rdv2Sq-0.9375*fr[14]*dfVfac_r*rdv2Sq-0.9375*fc[14]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[7]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[7]*dfVfac_r*rdv2Sq; 
  phaseFacr[15] = 0.5*fr[15]*fVfac_r*rdv2Sq-0.5*fc[15]*fVfac_r*rdv2Sq-0.4330127018922193*fr[11]*fVfac_r*rdv2Sq-0.4330127018922193*fc[11]*fVfac_r*rdv2Sq-0.9375*fr[15]*dfVfac_r*rdv2Sq-0.9375*fc[15]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[11]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[11]*dfVfac_r*rdv2Sq; 
  phaseFacr[16] = (-0.5412658773652742*fr[19]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[19]*dfVfac_r*rdv2Sq+0.5625*fr[16]*dfVfac_r*rdv2Sq-0.5625*fc[16]*dfVfac_r*rdv2Sq; 
  phaseFacr[17] = (-0.5412658773652742*fr[21]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[21]*dfVfac_r*rdv2Sq+0.5625*fr[17]*dfVfac_r*rdv2Sq-0.5625*fc[17]*dfVfac_r*rdv2Sq; 
  phaseFacr[18] = (-0.5412658773652742*fr[22]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[22]*dfVfac_r*rdv2Sq+0.5625*fr[18]*dfVfac_r*rdv2Sq-0.5625*fc[18]*dfVfac_r*rdv2Sq; 
  phaseFacr[19] = 0.5*fr[19]*fVfac_r*rdv2Sq-0.5*fc[19]*fVfac_r*rdv2Sq-0.4330127018922194*fr[16]*fVfac_r*rdv2Sq-0.4330127018922194*fc[16]*fVfac_r*rdv2Sq-0.9375*fr[19]*dfVfac_r*rdv2Sq-0.9375*fc[19]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[16]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[16]*dfVfac_r*rdv2Sq; 
  phaseFacr[20] = (-0.5412658773652742*fr[23]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[23]*dfVfac_r*rdv2Sq+0.5625*fr[20]*dfVfac_r*rdv2Sq-0.5625*fc[20]*dfVfac_r*rdv2Sq; 
  phaseFacr[21] = 0.5*fr[21]*fVfac_r*rdv2Sq-0.5*fc[21]*fVfac_r*rdv2Sq-0.4330127018922194*fr[17]*fVfac_r*rdv2Sq-0.4330127018922194*fc[17]*fVfac_r*rdv2Sq-0.9375*fr[21]*dfVfac_r*rdv2Sq-0.9375*fc[21]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[17]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[17]*dfVfac_r*rdv2Sq; 
  phaseFacr[22] = 0.5*fr[22]*fVfac_r*rdv2Sq-0.5*fc[22]*fVfac_r*rdv2Sq-0.4330127018922194*fr[18]*fVfac_r*rdv2Sq-0.4330127018922194*fc[18]*fVfac_r*rdv2Sq-0.9375*fr[22]*dfVfac_r*rdv2Sq-0.9375*fc[22]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[18]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[18]*dfVfac_r*rdv2Sq; 
  phaseFacr[23] = 0.5*fr[23]*fVfac_r*rdv2Sq-0.5*fc[23]*fVfac_r*rdv2Sq-0.4330127018922194*fr[20]*fVfac_r*rdv2Sq-0.4330127018922194*fc[20]*fVfac_r*rdv2Sq-0.9375*fr[23]*dfVfac_r*rdv2Sq-0.9375*fc[23]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[20]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[20]*dfVfac_r*rdv2Sq; 

  double incrl[24] = {0.0}; 
  incrl[0] = 0.5*confFac[3]*phaseFacl[5]+0.5*confFac[2]*phaseFacl[2]+0.5*confFac[1]*phaseFacl[1]+0.5*confFac[0]*phaseFacl[0]; 
  incrl[1] = 0.5*confFac[2]*phaseFacl[5]+0.5*phaseFacl[2]*confFac[3]+0.5*confFac[0]*phaseFacl[1]+0.5*phaseFacl[0]*confFac[1]; 
  incrl[2] = 0.5*confFac[1]*phaseFacl[5]+0.5*phaseFacl[1]*confFac[3]+0.5*confFac[0]*phaseFacl[2]+0.5*phaseFacl[0]*confFac[2]; 
  incrl[3] = 0.5*confFac[3]*phaseFacl[11]+0.5*confFac[2]*phaseFacl[7]+0.5*confFac[1]*phaseFacl[6]+0.5*confFac[0]*phaseFacl[3]; 
  incrl[4] = 0.5*confFac[3]*phaseFacl[12]+0.5*confFac[2]*phaseFacl[9]+0.5*confFac[1]*phaseFacl[8]+0.5*confFac[0]*phaseFacl[4]; 
  incrl[5] = 0.5*confFac[0]*phaseFacl[5]+0.5*phaseFacl[0]*confFac[3]+0.5*confFac[1]*phaseFacl[2]+0.5*phaseFacl[1]*confFac[2]; 
  incrl[6] = 0.5*confFac[2]*phaseFacl[11]+0.5*confFac[3]*phaseFacl[7]+0.5*confFac[0]*phaseFacl[6]+0.5*confFac[1]*phaseFacl[3]; 
  incrl[7] = 0.5*confFac[1]*phaseFacl[11]+0.5*confFac[0]*phaseFacl[7]+0.5*confFac[3]*phaseFacl[6]+0.5*confFac[2]*phaseFacl[3]; 
  incrl[8] = 0.5*confFac[2]*phaseFacl[12]+0.5*confFac[3]*phaseFacl[9]+0.5*confFac[0]*phaseFacl[8]+0.5*confFac[1]*phaseFacl[4]; 
  incrl[9] = 0.5*confFac[1]*phaseFacl[12]+0.5*confFac[0]*phaseFacl[9]+0.5*confFac[3]*phaseFacl[8]+0.5*confFac[2]*phaseFacl[4]; 
  incrl[10] = 0.5*confFac[3]*phaseFacl[15]+0.5*confFac[2]*phaseFacl[14]+0.5*confFac[1]*phaseFacl[13]+0.5*confFac[0]*phaseFacl[10]; 
  incrl[11] = 0.5*confFac[0]*phaseFacl[11]+0.5*confFac[1]*phaseFacl[7]+0.5*confFac[2]*phaseFacl[6]+0.5*confFac[3]*phaseFacl[3]; 
  incrl[12] = 0.5*confFac[0]*phaseFacl[12]+0.5*confFac[1]*phaseFacl[9]+0.5*confFac[2]*phaseFacl[8]+0.5*confFac[3]*phaseFacl[4]; 
  incrl[13] = 0.5*confFac[2]*phaseFacl[15]+0.5*confFac[3]*phaseFacl[14]+0.5*confFac[0]*phaseFacl[13]+0.5*confFac[1]*phaseFacl[10]; 
  incrl[14] = 0.5*confFac[1]*phaseFacl[15]+0.5*confFac[0]*phaseFacl[14]+0.5*confFac[3]*phaseFacl[13]+0.5*confFac[2]*phaseFacl[10]; 
  incrl[15] = 0.5*confFac[0]*phaseFacl[15]+0.5*confFac[1]*phaseFacl[14]+0.5*confFac[2]*phaseFacl[13]+0.5*confFac[3]*phaseFacl[10]; 
  incrl[16] = 0.5*confFac[3]*phaseFacl[20]+0.5000000000000001*confFac[2]*phaseFacl[18]+0.5000000000000001*confFac[1]*phaseFacl[17]+0.5*confFac[0]*phaseFacl[16]; 
  incrl[17] = 0.5000000000000001*confFac[2]*phaseFacl[20]+0.5*confFac[3]*phaseFacl[18]+0.5*confFac[0]*phaseFacl[17]+0.5000000000000001*confFac[1]*phaseFacl[16]; 
  incrl[18] = 0.5000000000000001*confFac[1]*phaseFacl[20]+0.5*confFac[0]*phaseFacl[18]+0.5*confFac[3]*phaseFacl[17]+0.5000000000000001*confFac[2]*phaseFacl[16]; 
  incrl[19] = 0.5*confFac[3]*phaseFacl[23]+0.5000000000000001*confFac[2]*phaseFacl[22]+0.5000000000000001*confFac[1]*phaseFacl[21]+0.5*confFac[0]*phaseFacl[19]; 
  incrl[20] = 0.5*confFac[0]*phaseFacl[20]+0.5000000000000001*confFac[1]*phaseFacl[18]+0.5000000000000001*confFac[2]*phaseFacl[17]+0.5*confFac[3]*phaseFacl[16]; 
  incrl[21] = 0.5000000000000001*confFac[2]*phaseFacl[23]+0.5*confFac[3]*phaseFacl[22]+0.5*confFac[0]*phaseFacl[21]+0.5000000000000001*confFac[1]*phaseFacl[19]; 
  incrl[22] = 0.5000000000000001*confFac[1]*phaseFacl[23]+0.5*confFac[0]*phaseFacl[22]+0.5*confFac[3]*phaseFacl[21]+0.5000000000000001*confFac[2]*phaseFacl[19]; 
  incrl[23] = 0.5*confFac[0]*phaseFacl[23]+0.5000000000000001*confFac[1]*phaseFacl[22]+0.5000000000000001*confFac[2]*phaseFacl[21]+0.5*confFac[3]*phaseFacl[19]; 

  double incrr[24] = {0.0}; 
  incrr[0] = 0.5*confFac[3]*phaseFacr[5]+0.5*confFac[2]*phaseFacr[2]+0.5*confFac[1]*phaseFacr[1]+0.5*confFac[0]*phaseFacr[0]; 
  incrr[1] = 0.5*confFac[2]*phaseFacr[5]+0.5*phaseFacr[2]*confFac[3]+0.5*confFac[0]*phaseFacr[1]+0.5*phaseFacr[0]*confFac[1]; 
  incrr[2] = 0.5*confFac[1]*phaseFacr[5]+0.5*phaseFacr[1]*confFac[3]+0.5*confFac[0]*phaseFacr[2]+0.5*phaseFacr[0]*confFac[2]; 
  incrr[3] = 0.5*confFac[3]*phaseFacr[11]+0.5*confFac[2]*phaseFacr[7]+0.5*confFac[1]*phaseFacr[6]+0.5*confFac[0]*phaseFacr[3]; 
  incrr[4] = 0.5*confFac[3]*phaseFacr[12]+0.5*confFac[2]*phaseFacr[9]+0.5*confFac[1]*phaseFacr[8]+0.5*confFac[0]*phaseFacr[4]; 
  incrr[5] = 0.5*confFac[0]*phaseFacr[5]+0.5*phaseFacr[0]*confFac[3]+0.5*confFac[1]*phaseFacr[2]+0.5*phaseFacr[1]*confFac[2]; 
  incrr[6] = 0.5*confFac[2]*phaseFacr[11]+0.5*confFac[3]*phaseFacr[7]+0.5*confFac[0]*phaseFacr[6]+0.5*confFac[1]*phaseFacr[3]; 
  incrr[7] = 0.5*confFac[1]*phaseFacr[11]+0.5*confFac[0]*phaseFacr[7]+0.5*confFac[3]*phaseFacr[6]+0.5*confFac[2]*phaseFacr[3]; 
  incrr[8] = 0.5*confFac[2]*phaseFacr[12]+0.5*confFac[3]*phaseFacr[9]+0.5*confFac[0]*phaseFacr[8]+0.5*confFac[1]*phaseFacr[4]; 
  incrr[9] = 0.5*confFac[1]*phaseFacr[12]+0.5*confFac[0]*phaseFacr[9]+0.5*confFac[3]*phaseFacr[8]+0.5*confFac[2]*phaseFacr[4]; 
  incrr[10] = 0.5*confFac[3]*phaseFacr[15]+0.5*confFac[2]*phaseFacr[14]+0.5*confFac[1]*phaseFacr[13]+0.5*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.5*confFac[0]*phaseFacr[11]+0.5*confFac[1]*phaseFacr[7]+0.5*confFac[2]*phaseFacr[6]+0.5*confFac[3]*phaseFacr[3]; 
  incrr[12] = 0.5*confFac[0]*phaseFacr[12]+0.5*confFac[1]*phaseFacr[9]+0.5*confFac[2]*phaseFacr[8]+0.5*confFac[3]*phaseFacr[4]; 
  incrr[13] = 0.5*confFac[2]*phaseFacr[15]+0.5*confFac[3]*phaseFacr[14]+0.5*confFac[0]*phaseFacr[13]+0.5*confFac[1]*phaseFacr[10]; 
  incrr[14] = 0.5*confFac[1]*phaseFacr[15]+0.5*confFac[0]*phaseFacr[14]+0.5*confFac[3]*phaseFacr[13]+0.5*confFac[2]*phaseFacr[10]; 
  incrr[15] = 0.5*confFac[0]*phaseFacr[15]+0.5*confFac[1]*phaseFacr[14]+0.5*confFac[2]*phaseFacr[13]+0.5*confFac[3]*phaseFacr[10]; 
  incrr[16] = 0.5*confFac[3]*phaseFacr[20]+0.5000000000000001*confFac[2]*phaseFacr[18]+0.5000000000000001*confFac[1]*phaseFacr[17]+0.5*confFac[0]*phaseFacr[16]; 
  incrr[17] = 0.5000000000000001*confFac[2]*phaseFacr[20]+0.5*confFac[3]*phaseFacr[18]+0.5*confFac[0]*phaseFacr[17]+0.5000000000000001*confFac[1]*phaseFacr[16]; 
  incrr[18] = 0.5000000000000001*confFac[1]*phaseFacr[20]+0.5*confFac[0]*phaseFacr[18]+0.5*confFac[3]*phaseFacr[17]+0.5000000000000001*confFac[2]*phaseFacr[16]; 
  incrr[19] = 0.5*confFac[3]*phaseFacr[23]+0.5000000000000001*confFac[2]*phaseFacr[22]+0.5000000000000001*confFac[1]*phaseFacr[21]+0.5*confFac[0]*phaseFacr[19]; 
  incrr[20] = 0.5*confFac[0]*phaseFacr[20]+0.5000000000000001*confFac[1]*phaseFacr[18]+0.5000000000000001*confFac[2]*phaseFacr[17]+0.5*confFac[3]*phaseFacr[16]; 
  incrr[21] = 0.5000000000000001*confFac[2]*phaseFacr[23]+0.5*confFac[3]*phaseFacr[22]+0.5*confFac[0]*phaseFacr[21]+0.5000000000000001*confFac[1]*phaseFacr[19]; 
  incrr[22] = 0.5000000000000001*confFac[1]*phaseFacr[23]+0.5*confFac[0]*phaseFacr[22]+0.5*confFac[3]*phaseFacr[21]+0.5000000000000001*confFac[2]*phaseFacr[19]; 
  incrr[23] = 0.5*confFac[0]*phaseFacr[23]+0.5000000000000001*confFac[1]*phaseFacr[22]+0.5000000000000001*confFac[2]*phaseFacr[21]+0.5*confFac[3]*phaseFacr[19]; 

  out[0] += incrr[0]-1.0*incrl[0]; 
  out[1] += incrr[1]-1.0*incrl[1]; 
  out[2] += incrr[2]-1.0*incrl[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += incrr[4]-1.0*incrl[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += incrr[8]-1.0*incrl[8]; 
  out[9] += incrr[9]-1.0*incrl[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 
  out[12] += incrr[12]-1.0*incrl[12]; 
  out[13] += incrr[13]-1.0*incrl[13]; 
  out[14] += incrr[14]-1.0*incrl[14]; 
  out[15] += incrr[15]-1.0*incrl[15]; 
  out[16] += incrr[16]-1.0*incrl[16]; 
  out[17] += incrr[17]-1.0*incrl[17]; 
  out[18] += incrr[18]-1.0*incrl[18]; 
  out[19] += incrr[19]-1.0*incrl[19]; 
  out[20] += incrr[20]-1.0*incrl[20]; 
  out[21] += incrr[21]-1.0*incrl[21]; 
  out[22] += incrr[22]-1.0*incrl[22]; 
  out[23] += incrr[23]-1.0*incrl[23]; 

  return 0.;

} 