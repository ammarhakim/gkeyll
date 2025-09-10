#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_mapped_surfmu_1x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvell, const double *jacobvelc, const double *jacobvelr, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double fl_over_jacv[12], fc_over_jacv[12], fr_over_jacv[12];
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

  double dvl = 2.449489742783178*vmapl[3];
  double dvc = 2.449489742783178*vmapc[3];
  double dvr = 2.449489742783178*vmapr[3];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdv2 = 2.0/dxv[2]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[2] = {0.}; 
  confFac[0] = 1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = 1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*m_; 

  double dfVfac_l = vmap_prime[0]*(0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]); 
  double dfVfac_r = vmap_prime[0]*(1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]); 

  double fVfac_l = (jacobvelc[0]*(0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]))/vmap_primeSq; 
  double fVfac_r = (jacobvelc[0]*(1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]))/vmap_primeSq; 

  double phaseFacl[12] = {0.0}; 
  const double dvlR2 = pow(dvl,2);
  const double dvlR3 = pow(dvl,3);
  const double dvlR4 = pow(dvl,4);
  const double dvcR2 = pow(dvc,2);
  const double dvcR3 = pow(dvc,3);
  const double dvcR4 = pow(dvc,4);

  phaseFacl[0] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[3]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[3]*dfVfac_l*rdv2-9.0*fl_over_jacv[0]*dfVfac_l*rdv2+9.0*fc_over_jacv[0]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[3]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[3]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[3]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[3]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[1] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[5]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[5]*dfVfac_l*rdv2-9.0*fl_over_jacv[1]*dfVfac_l*rdv2+9.0*fc_over_jacv[1]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[5]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[5]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[5]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[5]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[2] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[6]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[6]*dfVfac_l*rdv2-9.0*fl_over_jacv[2]*dfVfac_l*rdv2+9.0*fc_over_jacv[2]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[6]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[6]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[6]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[6]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[3] = (dvcR3*((-3.0*fl_over_jacv[3]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[3]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[3]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[3]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[3]*dfVfac_l*rdv2+15.0*fc_over_jacv[3]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[0]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[0]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[3]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[3]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[3]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[3]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[4] = (dvcR2*dvlR2*((-8.660254037844386*fl_over_jacv[7]*dfVfac_l*rdv2)-8.660254037844386*fc_over_jacv[7]*dfVfac_l*rdv2-9.0*fl_over_jacv[4]*dfVfac_l*rdv2+9.0*fc_over_jacv[4]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fc_over_jacv[7]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fc_over_jacv[7]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(1.732050807568877*fl_over_jacv[7]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(1.732050807568877*fl_over_jacv[7]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[5] = (dvcR3*((-3.0*fl_over_jacv[5]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[5]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[5]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[5]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[5]*dfVfac_l*rdv2+15.0*fc_over_jacv[5]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[1]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[1]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[5]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[5]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[5]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[5]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[6] = (dvcR3*((-3.0*fl_over_jacv[6]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[6]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[6]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[6]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[6]*dfVfac_l*rdv2+15.0*fc_over_jacv[6]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[2]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[2]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[6]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[6]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[6]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[6]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[7] = (dvcR3*((-3.0*fl_over_jacv[7]*fVfac_l*rdv2Sq)-1.732050807568877*fl_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvl*((-5.0*fl_over_jacv[7]*fVfac_l*rdv2Sq)-5.196152422706631*fl_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvc*dvlR2*(5.0*fc_over_jacv[7]*fVfac_l*rdv2Sq-5.196152422706631*fc_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvlR3*(3.0*fc_over_jacv[7]*fVfac_l*rdv2Sq-1.732050807568877*fc_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dvlR3+6.0*dvc*dvlR2+6.0*dvcR2*dvl+2.0*dvcR3)+(dvcR2*dvlR2*(15.0*fl_over_jacv[7]*dfVfac_l*rdv2+15.0*fc_over_jacv[7]*dfVfac_l*rdv2+15.58845726811989*fl_over_jacv[4]*dfVfac_l*rdv2-15.58845726811989*fc_over_jacv[4]*dfVfac_l*rdv2))/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fc_over_jacv[7]*dfVfac_l*dvlR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fc_over_jacv[7]*dfVfac_l*dvc*dvlR3*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)+(3.0*fl_over_jacv[7]*dfVfac_l*dvcR3*dvl*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl)-(3.0*fl_over_jacv[7]*dfVfac_l*dvcR4*rdv2)/(dvc*dvlR4+3.0*dvcR2*dvlR3+3.0*dvcR3*dvlR2+dvcR4*dvl); 
  phaseFacl[8] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[10]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[10]*dfVfac_l*rdv2-45.0*fl_over_jacv[8]*dfVfac_l*rdv2+45.0*fc_over_jacv[8]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[10]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[10]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[10]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[10]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[9] = (dvcR2*dvlR2*((-43.30127018922195*fl_over_jacv[11]*dfVfac_l*rdv2)-43.30127018922195*fc_over_jacv[11]*dfVfac_l*rdv2-45.0*fl_over_jacv[9]*dfVfac_l*rdv2+45.0*fc_over_jacv[9]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fc_over_jacv[11]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fc_over_jacv[11]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(8.660254037844387*fl_over_jacv[11]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(8.660254037844387*fl_over_jacv[11]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[10] = (dvcR3*((-15.0*fl_over_jacv[10]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[10]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[10]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[10]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[10]*dfVfac_l*rdv2+75.0*fc_over_jacv[10]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[8]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[8]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[10]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[10]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[10]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[10]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 
  phaseFacl[11] = (dvcR3*((-15.0*fl_over_jacv[11]*fVfac_l*rdv2Sq)-8.660254037844387*fl_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvl*((-25.0*fl_over_jacv[11]*fVfac_l*rdv2Sq)-25.98076211353316*fl_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvc*dvlR2*(25.0*fc_over_jacv[11]*fVfac_l*rdv2Sq-25.98076211353316*fc_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvlR3*(15.0*fc_over_jacv[11]*fVfac_l*rdv2Sq-8.660254037844387*fc_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dvlR3+30.0*dvc*dvlR2+30.0*dvcR2*dvl+10.0*dvcR3)+(dvcR2*dvlR2*(75.0*fl_over_jacv[11]*dfVfac_l*rdv2+75.0*fc_over_jacv[11]*dfVfac_l*rdv2+77.94228634059948*fl_over_jacv[9]*dfVfac_l*rdv2-77.94228634059948*fc_over_jacv[9]*dfVfac_l*rdv2))/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fc_over_jacv[11]*dfVfac_l*dvlR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fc_over_jacv[11]*dfVfac_l*dvc*dvlR3*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)+(15.0*fl_over_jacv[11]*dfVfac_l*dvcR3*dvl*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl)-(15.0*fl_over_jacv[11]*dfVfac_l*dvcR4*rdv2)/(5.0*dvc*dvlR4+15.0*dvcR2*dvlR3+15.0*dvcR3*dvlR2+5.0*dvcR4*dvl); 

  double phaseFacr[12] = {0.0}; 
  const double dvrR2 = pow(dvr,2);
  const double dvrR3 = pow(dvr,3);
  const double dvrR4 = pow(dvr,4);

  phaseFacr[0] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[3]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[3]*dfVfac_r*rdv2+9.0*fr_over_jacv[0]*dfVfac_r*rdv2-9.0*fc_over_jacv[0]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[3]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[3]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[3]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[3]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[1] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[5]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[5]*dfVfac_r*rdv2+9.0*fr_over_jacv[1]*dfVfac_r*rdv2-9.0*fc_over_jacv[1]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[5]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[5]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[5]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[5]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[2] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[6]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[6]*dfVfac_r*rdv2+9.0*fr_over_jacv[2]*dfVfac_r*rdv2-9.0*fc_over_jacv[2]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[6]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[6]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[6]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[6]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[3] = (dvcR2*dvr*(5.0*fr_over_jacv[3]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[3]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[3]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[3]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[3]*dfVfac_r*rdv2)-15.0*fc_over_jacv[3]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[0]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[0]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[3]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[3]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[3]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[3]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[4] = (dvcR2*dvrR2*((-8.660254037844386*fr_over_jacv[7]*dfVfac_r*rdv2)-8.660254037844386*fc_over_jacv[7]*dfVfac_r*rdv2+9.0*fr_over_jacv[4]*dfVfac_r*rdv2-9.0*fc_over_jacv[4]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fc_over_jacv[7]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fc_over_jacv[7]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(1.732050807568877*fr_over_jacv[7]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(1.732050807568877*fr_over_jacv[7]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[5] = (dvcR2*dvr*(5.0*fr_over_jacv[5]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[5]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[5]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[5]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[5]*dfVfac_r*rdv2)-15.0*fc_over_jacv[5]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[1]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[1]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[5]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[5]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[5]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[5]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[6] = (dvcR2*dvr*(5.0*fr_over_jacv[6]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[6]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[6]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[6]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[6]*dfVfac_r*rdv2)-15.0*fc_over_jacv[6]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[2]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[2]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[6]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[6]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[6]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[6]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[7] = (dvcR2*dvr*(5.0*fr_over_jacv[7]*fVfac_r*rdv2Sq-5.196152422706631*fr_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR3*(3.0*fr_over_jacv[7]*fVfac_r*rdv2Sq-1.732050807568877*fr_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvrR3*((-3.0*fc_over_jacv[7]*fVfac_r*rdv2Sq)-1.732050807568877*fc_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvc*dvrR2*((-5.0*fc_over_jacv[7]*fVfac_r*rdv2Sq)-5.196152422706631*fc_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dvrR3+6.0*dvc*dvrR2+6.0*dvcR2*dvr+2.0*dvcR3)+(dvcR2*dvrR2*((-15.0*fr_over_jacv[7]*dfVfac_r*rdv2)-15.0*fc_over_jacv[7]*dfVfac_r*rdv2+15.58845726811989*fr_over_jacv[4]*dfVfac_r*rdv2-15.58845726811989*fc_over_jacv[4]*dfVfac_r*rdv2))/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fc_over_jacv[7]*dfVfac_r*dvrR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fc_over_jacv[7]*dfVfac_r*dvc*dvrR3*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)-(3.0*fr_over_jacv[7]*dfVfac_r*dvcR3*dvr*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr)+(3.0*fr_over_jacv[7]*dfVfac_r*dvcR4*rdv2)/(dvc*dvrR4+3.0*dvcR2*dvrR3+3.0*dvcR3*dvrR2+dvcR4*dvr); 
  phaseFacr[8] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[10]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[10]*dfVfac_r*rdv2+45.0*fr_over_jacv[8]*dfVfac_r*rdv2-45.0*fc_over_jacv[8]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[10]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[10]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[10]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[10]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[9] = (dvcR2*dvrR2*((-43.30127018922195*fr_over_jacv[11]*dfVfac_r*rdv2)-43.30127018922195*fc_over_jacv[11]*dfVfac_r*rdv2+45.0*fr_over_jacv[9]*dfVfac_r*rdv2-45.0*fc_over_jacv[9]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fc_over_jacv[11]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fc_over_jacv[11]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(8.660254037844387*fr_over_jacv[11]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(8.660254037844387*fr_over_jacv[11]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[10] = (dvcR2*dvr*(25.0*fr_over_jacv[10]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[10]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[10]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[10]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[10]*dfVfac_r*rdv2)-75.0*fc_over_jacv[10]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[8]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[8]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[10]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[10]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[10]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[10]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 
  phaseFacr[11] = (dvcR2*dvr*(25.0*fr_over_jacv[11]*fVfac_r*rdv2Sq-25.98076211353316*fr_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR3*(15.0*fr_over_jacv[11]*fVfac_r*rdv2Sq-8.660254037844387*fr_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvrR3*((-15.0*fc_over_jacv[11]*fVfac_r*rdv2Sq)-8.660254037844387*fc_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvc*dvrR2*((-25.0*fc_over_jacv[11]*fVfac_r*rdv2Sq)-25.98076211353316*fc_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dvrR3+30.0*dvc*dvrR2+30.0*dvcR2*dvr+10.0*dvcR3)+(dvcR2*dvrR2*((-75.0*fr_over_jacv[11]*dfVfac_r*rdv2)-75.0*fc_over_jacv[11]*dfVfac_r*rdv2+77.94228634059948*fr_over_jacv[9]*dfVfac_r*rdv2-77.94228634059948*fc_over_jacv[9]*dfVfac_r*rdv2))/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fc_over_jacv[11]*dfVfac_r*dvrR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fc_over_jacv[11]*dfVfac_r*dvc*dvrR3*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)-(15.0*fr_over_jacv[11]*dfVfac_r*dvcR3*dvr*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr)+(15.0*fr_over_jacv[11]*dfVfac_r*dvcR4*rdv2)/(5.0*dvc*dvrR4+15.0*dvcR2*dvrR3+15.0*dvcR3*dvrR2+5.0*dvcR4*dvr); 

  double incrl[12] = {0.0}; 
  incrl[0] = 0.7071067811865475*confFac[1]*phaseFacl[1]+0.7071067811865475*confFac[0]*phaseFacl[0]; 
  incrl[1] = 0.7071067811865475*confFac[0]*phaseFacl[1]+0.7071067811865475*phaseFacl[0]*confFac[1]; 
  incrl[2] = 0.7071067811865475*confFac[1]*phaseFacl[4]+0.7071067811865475*confFac[0]*phaseFacl[2]; 
  incrl[3] = 0.7071067811865475*confFac[1]*phaseFacl[5]+0.7071067811865475*confFac[0]*phaseFacl[3]; 
  incrl[4] = 0.7071067811865475*confFac[0]*phaseFacl[4]+0.7071067811865475*confFac[1]*phaseFacl[2]; 
  incrl[5] = 0.7071067811865475*confFac[0]*phaseFacl[5]+0.7071067811865475*confFac[1]*phaseFacl[3]; 
  incrl[6] = 0.7071067811865475*confFac[1]*phaseFacl[7]+0.7071067811865475*confFac[0]*phaseFacl[6]; 
  incrl[7] = 0.7071067811865475*confFac[0]*phaseFacl[7]+0.7071067811865475*confFac[1]*phaseFacl[6]; 
  incrl[8] = 0.7071067811865475*confFac[1]*phaseFacl[9]+0.7071067811865475*confFac[0]*phaseFacl[8]; 
  incrl[9] = 0.7071067811865475*confFac[0]*phaseFacl[9]+0.7071067811865475*confFac[1]*phaseFacl[8]; 
  incrl[10] = 0.7071067811865475*confFac[1]*phaseFacl[11]+0.7071067811865475*confFac[0]*phaseFacl[10]; 
  incrl[11] = 0.7071067811865475*confFac[0]*phaseFacl[11]+0.7071067811865475*confFac[1]*phaseFacl[10]; 

  double incrr[12] = {0.0}; 
  incrr[0] = 0.7071067811865475*confFac[1]*phaseFacr[1]+0.7071067811865475*confFac[0]*phaseFacr[0]; 
  incrr[1] = 0.7071067811865475*confFac[0]*phaseFacr[1]+0.7071067811865475*phaseFacr[0]*confFac[1]; 
  incrr[2] = 0.7071067811865475*confFac[1]*phaseFacr[4]+0.7071067811865475*confFac[0]*phaseFacr[2]; 
  incrr[3] = 0.7071067811865475*confFac[1]*phaseFacr[5]+0.7071067811865475*confFac[0]*phaseFacr[3]; 
  incrr[4] = 0.7071067811865475*confFac[0]*phaseFacr[4]+0.7071067811865475*confFac[1]*phaseFacr[2]; 
  incrr[5] = 0.7071067811865475*confFac[0]*phaseFacr[5]+0.7071067811865475*confFac[1]*phaseFacr[3]; 
  incrr[6] = 0.7071067811865475*confFac[1]*phaseFacr[7]+0.7071067811865475*confFac[0]*phaseFacr[6]; 
  incrr[7] = 0.7071067811865475*confFac[0]*phaseFacr[7]+0.7071067811865475*confFac[1]*phaseFacr[6]; 
  incrr[8] = 0.7071067811865475*confFac[1]*phaseFacr[9]+0.7071067811865475*confFac[0]*phaseFacr[8]; 
  incrr[9] = 0.7071067811865475*confFac[0]*phaseFacr[9]+0.7071067811865475*confFac[1]*phaseFacr[8]; 
  incrr[10] = 0.7071067811865475*confFac[1]*phaseFacr[11]+0.7071067811865475*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.7071067811865475*confFac[0]*phaseFacr[11]+0.7071067811865475*confFac[1]*phaseFacr[10]; 

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

  return 0.;

} 
GKYL_CU_DH double lbo_gyrokinetic_diff_notmapped_surfmu_1x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvell, const double *jacobvelc, const double *jacobvelr, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdv2 = 2.0/dxv[2]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[2] = {0.}; 
  confFac[0] = 1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = 1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*m_; 

  double dfVfac_l = 0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]; 
  double dfVfac_r = 1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]; 

  double fVfac_l = 0.7071067811865475*vmapc[2]-1.224744871391589*vmapc[3]; 
  double fVfac_r = 1.224744871391589*vmapc[3]+0.7071067811865475*vmapc[2]; 

  double phaseFacl[12] = {0.0}; 

  phaseFacl[0] = (-0.5412658773652741*fl[3]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[3]*dfVfac_l*rdv2Sq-0.5625*fl[0]*dfVfac_l*rdv2Sq+0.5625*fc[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[1] = (-0.5412658773652741*fl[5]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[5]*dfVfac_l*rdv2Sq-0.5625*fl[1]*dfVfac_l*rdv2Sq+0.5625*fc[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[2] = (-0.5412658773652741*fl[6]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[6]*dfVfac_l*rdv2Sq-0.5625*fl[2]*dfVfac_l*rdv2Sq+0.5625*fc[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[3] = (-0.5*fl[3]*fVfac_l*rdv2Sq)+0.5*fc[3]*fVfac_l*rdv2Sq-0.4330127018922193*fl[0]*fVfac_l*rdv2Sq-0.4330127018922193*fc[0]*fVfac_l*rdv2Sq+0.9375*fl[3]*dfVfac_l*rdv2Sq+0.9375*fc[3]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[0]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[4] = (-0.5412658773652741*fl[7]*dfVfac_l*rdv2Sq)-0.5412658773652741*fc[7]*dfVfac_l*rdv2Sq-0.5625*fl[4]*dfVfac_l*rdv2Sq+0.5625*fc[4]*dfVfac_l*rdv2Sq; 
  phaseFacl[5] = (-0.5*fl[5]*fVfac_l*rdv2Sq)+0.5*fc[5]*fVfac_l*rdv2Sq-0.4330127018922193*fl[1]*fVfac_l*rdv2Sq-0.4330127018922193*fc[1]*fVfac_l*rdv2Sq+0.9375*fl[5]*dfVfac_l*rdv2Sq+0.9375*fc[5]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[1]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[6] = (-0.5*fl[6]*fVfac_l*rdv2Sq)+0.5*fc[6]*fVfac_l*rdv2Sq-0.4330127018922193*fl[2]*fVfac_l*rdv2Sq-0.4330127018922193*fc[2]*fVfac_l*rdv2Sq+0.9375*fl[6]*dfVfac_l*rdv2Sq+0.9375*fc[6]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[2]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[7] = (-0.5*fl[7]*fVfac_l*rdv2Sq)+0.5*fc[7]*fVfac_l*rdv2Sq-0.4330127018922193*fl[4]*fVfac_l*rdv2Sq-0.4330127018922193*fc[4]*fVfac_l*rdv2Sq+0.9375*fl[7]*dfVfac_l*rdv2Sq+0.9375*fc[7]*dfVfac_l*rdv2Sq+0.9742785792574932*fl[4]*dfVfac_l*rdv2Sq-0.9742785792574932*fc[4]*dfVfac_l*rdv2Sq; 
  phaseFacl[8] = (-0.5412658773652742*fl[10]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[10]*dfVfac_l*rdv2Sq-0.5625*fl[8]*dfVfac_l*rdv2Sq+0.5625*fc[8]*dfVfac_l*rdv2Sq; 
  phaseFacl[9] = (-0.5412658773652742*fl[11]*dfVfac_l*rdv2Sq)-0.5412658773652742*fc[11]*dfVfac_l*rdv2Sq-0.5625*fl[9]*dfVfac_l*rdv2Sq+0.5625*fc[9]*dfVfac_l*rdv2Sq; 
  phaseFacl[10] = (-0.5*fl[10]*fVfac_l*rdv2Sq)+0.5*fc[10]*fVfac_l*rdv2Sq-0.4330127018922194*fl[8]*fVfac_l*rdv2Sq-0.4330127018922194*fc[8]*fVfac_l*rdv2Sq+0.9375*fl[10]*dfVfac_l*rdv2Sq+0.9375*fc[10]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[8]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[8]*dfVfac_l*rdv2Sq; 
  phaseFacl[11] = (-0.5*fl[11]*fVfac_l*rdv2Sq)+0.5*fc[11]*fVfac_l*rdv2Sq-0.4330127018922194*fl[9]*fVfac_l*rdv2Sq-0.4330127018922194*fc[9]*fVfac_l*rdv2Sq+0.9375*fl[11]*dfVfac_l*rdv2Sq+0.9375*fc[11]*dfVfac_l*rdv2Sq+0.9742785792574935*fl[9]*dfVfac_l*rdv2Sq-0.9742785792574935*fc[9]*dfVfac_l*rdv2Sq; 

  double phaseFacr[12] = {0.0}; 

  phaseFacr[0] = (-0.5412658773652741*fr[3]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[3]*dfVfac_r*rdv2Sq+0.5625*fr[0]*dfVfac_r*rdv2Sq-0.5625*fc[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[1] = (-0.5412658773652741*fr[5]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[5]*dfVfac_r*rdv2Sq+0.5625*fr[1]*dfVfac_r*rdv2Sq-0.5625*fc[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[2] = (-0.5412658773652741*fr[6]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[6]*dfVfac_r*rdv2Sq+0.5625*fr[2]*dfVfac_r*rdv2Sq-0.5625*fc[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[3] = 0.5*fr[3]*fVfac_r*rdv2Sq-0.5*fc[3]*fVfac_r*rdv2Sq-0.4330127018922193*fr[0]*fVfac_r*rdv2Sq-0.4330127018922193*fc[0]*fVfac_r*rdv2Sq-0.9375*fr[3]*dfVfac_r*rdv2Sq-0.9375*fc[3]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[0]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[4] = (-0.5412658773652741*fr[7]*dfVfac_r*rdv2Sq)-0.5412658773652741*fc[7]*dfVfac_r*rdv2Sq+0.5625*fr[4]*dfVfac_r*rdv2Sq-0.5625*fc[4]*dfVfac_r*rdv2Sq; 
  phaseFacr[5] = 0.5*fr[5]*fVfac_r*rdv2Sq-0.5*fc[5]*fVfac_r*rdv2Sq-0.4330127018922193*fr[1]*fVfac_r*rdv2Sq-0.4330127018922193*fc[1]*fVfac_r*rdv2Sq-0.9375*fr[5]*dfVfac_r*rdv2Sq-0.9375*fc[5]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[1]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[6] = 0.5*fr[6]*fVfac_r*rdv2Sq-0.5*fc[6]*fVfac_r*rdv2Sq-0.4330127018922193*fr[2]*fVfac_r*rdv2Sq-0.4330127018922193*fc[2]*fVfac_r*rdv2Sq-0.9375*fr[6]*dfVfac_r*rdv2Sq-0.9375*fc[6]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[2]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[7] = 0.5*fr[7]*fVfac_r*rdv2Sq-0.5*fc[7]*fVfac_r*rdv2Sq-0.4330127018922193*fr[4]*fVfac_r*rdv2Sq-0.4330127018922193*fc[4]*fVfac_r*rdv2Sq-0.9375*fr[7]*dfVfac_r*rdv2Sq-0.9375*fc[7]*dfVfac_r*rdv2Sq+0.9742785792574932*fr[4]*dfVfac_r*rdv2Sq-0.9742785792574932*fc[4]*dfVfac_r*rdv2Sq; 
  phaseFacr[8] = (-0.5412658773652742*fr[10]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[10]*dfVfac_r*rdv2Sq+0.5625*fr[8]*dfVfac_r*rdv2Sq-0.5625*fc[8]*dfVfac_r*rdv2Sq; 
  phaseFacr[9] = (-0.5412658773652742*fr[11]*dfVfac_r*rdv2Sq)-0.5412658773652742*fc[11]*dfVfac_r*rdv2Sq+0.5625*fr[9]*dfVfac_r*rdv2Sq-0.5625*fc[9]*dfVfac_r*rdv2Sq; 
  phaseFacr[10] = 0.5*fr[10]*fVfac_r*rdv2Sq-0.5*fc[10]*fVfac_r*rdv2Sq-0.4330127018922194*fr[8]*fVfac_r*rdv2Sq-0.4330127018922194*fc[8]*fVfac_r*rdv2Sq-0.9375*fr[10]*dfVfac_r*rdv2Sq-0.9375*fc[10]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[8]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[8]*dfVfac_r*rdv2Sq; 
  phaseFacr[11] = 0.5*fr[11]*fVfac_r*rdv2Sq-0.5*fc[11]*fVfac_r*rdv2Sq-0.4330127018922194*fr[9]*fVfac_r*rdv2Sq-0.4330127018922194*fc[9]*fVfac_r*rdv2Sq-0.9375*fr[11]*dfVfac_r*rdv2Sq-0.9375*fc[11]*dfVfac_r*rdv2Sq+0.9742785792574935*fr[9]*dfVfac_r*rdv2Sq-0.9742785792574935*fc[9]*dfVfac_r*rdv2Sq; 

  double incrl[12] = {0.0}; 
  incrl[0] = 0.7071067811865475*confFac[1]*phaseFacl[1]+0.7071067811865475*confFac[0]*phaseFacl[0]; 
  incrl[1] = 0.7071067811865475*confFac[0]*phaseFacl[1]+0.7071067811865475*phaseFacl[0]*confFac[1]; 
  incrl[2] = 0.7071067811865475*confFac[1]*phaseFacl[4]+0.7071067811865475*confFac[0]*phaseFacl[2]; 
  incrl[3] = 0.7071067811865475*confFac[1]*phaseFacl[5]+0.7071067811865475*confFac[0]*phaseFacl[3]; 
  incrl[4] = 0.7071067811865475*confFac[0]*phaseFacl[4]+0.7071067811865475*confFac[1]*phaseFacl[2]; 
  incrl[5] = 0.7071067811865475*confFac[0]*phaseFacl[5]+0.7071067811865475*confFac[1]*phaseFacl[3]; 
  incrl[6] = 0.7071067811865475*confFac[1]*phaseFacl[7]+0.7071067811865475*confFac[0]*phaseFacl[6]; 
  incrl[7] = 0.7071067811865475*confFac[0]*phaseFacl[7]+0.7071067811865475*confFac[1]*phaseFacl[6]; 
  incrl[8] = 0.7071067811865475*confFac[1]*phaseFacl[9]+0.7071067811865475*confFac[0]*phaseFacl[8]; 
  incrl[9] = 0.7071067811865475*confFac[0]*phaseFacl[9]+0.7071067811865475*confFac[1]*phaseFacl[8]; 
  incrl[10] = 0.7071067811865475*confFac[1]*phaseFacl[11]+0.7071067811865475*confFac[0]*phaseFacl[10]; 
  incrl[11] = 0.7071067811865475*confFac[0]*phaseFacl[11]+0.7071067811865475*confFac[1]*phaseFacl[10]; 

  double incrr[12] = {0.0}; 
  incrr[0] = 0.7071067811865475*confFac[1]*phaseFacr[1]+0.7071067811865475*confFac[0]*phaseFacr[0]; 
  incrr[1] = 0.7071067811865475*confFac[0]*phaseFacr[1]+0.7071067811865475*phaseFacr[0]*confFac[1]; 
  incrr[2] = 0.7071067811865475*confFac[1]*phaseFacr[4]+0.7071067811865475*confFac[0]*phaseFacr[2]; 
  incrr[3] = 0.7071067811865475*confFac[1]*phaseFacr[5]+0.7071067811865475*confFac[0]*phaseFacr[3]; 
  incrr[4] = 0.7071067811865475*confFac[0]*phaseFacr[4]+0.7071067811865475*confFac[1]*phaseFacr[2]; 
  incrr[5] = 0.7071067811865475*confFac[0]*phaseFacr[5]+0.7071067811865475*confFac[1]*phaseFacr[3]; 
  incrr[6] = 0.7071067811865475*confFac[1]*phaseFacr[7]+0.7071067811865475*confFac[0]*phaseFacr[6]; 
  incrr[7] = 0.7071067811865475*confFac[0]*phaseFacr[7]+0.7071067811865475*confFac[1]*phaseFacr[6]; 
  incrr[8] = 0.7071067811865475*confFac[1]*phaseFacr[9]+0.7071067811865475*confFac[0]*phaseFacr[8]; 
  incrr[9] = 0.7071067811865475*confFac[0]*phaseFacr[9]+0.7071067811865475*confFac[1]*phaseFacr[8]; 
  incrr[10] = 0.7071067811865475*confFac[1]*phaseFacr[11]+0.7071067811865475*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.7071067811865475*confFac[0]*phaseFacr[11]+0.7071067811865475*confFac[1]*phaseFacr[10]; 

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

  return 0.;

} 