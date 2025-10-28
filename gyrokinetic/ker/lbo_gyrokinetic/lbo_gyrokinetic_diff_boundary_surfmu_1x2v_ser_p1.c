#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_mapped_boundary_surfmu_1x2v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel_edge, const double *jacobvel_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // dxv: Cell spacing. 
  // vmap_edge,vmap_skin: velocity space mapping.
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // jacobvel_edge,jacobvel_skin: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  double vmap_primeSq = pow(vmap_prime[1],2);

  double fedge_over_jacv[12], fskin_over_jacv[12];
  fedge_over_jacv[0] = fedge[0]/jacobvel_edge[0]; 
  fedge_over_jacv[1] = fedge[1]/jacobvel_edge[0]; 
  fedge_over_jacv[2] = fedge[2]/jacobvel_edge[0]; 
  fedge_over_jacv[3] = fedge[3]/jacobvel_edge[0]; 
  fedge_over_jacv[4] = fedge[4]/jacobvel_edge[0]; 
  fedge_over_jacv[5] = fedge[5]/jacobvel_edge[0]; 
  fedge_over_jacv[6] = fedge[6]/jacobvel_edge[0]; 
  fedge_over_jacv[7] = fedge[7]/jacobvel_edge[0]; 
  fedge_over_jacv[8] = fedge[8]/jacobvel_edge[0]; 
  fedge_over_jacv[9] = fedge[9]/jacobvel_edge[0]; 
  fedge_over_jacv[10] = fedge[10]/jacobvel_edge[0]; 
  fedge_over_jacv[11] = fedge[11]/jacobvel_edge[0]; 

  fskin_over_jacv[0] = fskin[0]/jacobvel_skin[0]; 
  fskin_over_jacv[1] = fskin[1]/jacobvel_skin[0]; 
  fskin_over_jacv[2] = fskin[2]/jacobvel_skin[0]; 
  fskin_over_jacv[3] = fskin[3]/jacobvel_skin[0]; 
  fskin_over_jacv[4] = fskin[4]/jacobvel_skin[0]; 
  fskin_over_jacv[5] = fskin[5]/jacobvel_skin[0]; 
  fskin_over_jacv[6] = fskin[6]/jacobvel_skin[0]; 
  fskin_over_jacv[7] = fskin[7]/jacobvel_skin[0]; 
  fskin_over_jacv[8] = fskin[8]/jacobvel_skin[0]; 
  fskin_over_jacv[9] = fskin[9]/jacobvel_skin[0]; 
  fskin_over_jacv[10] = fskin[10]/jacobvel_skin[0]; 
  fskin_over_jacv[11] = fskin[11]/jacobvel_skin[0]; 

  double dv_edge = 2.449489742783178*vmap_edge[3];
  double dv_skin = 2.449489742783178*vmap_skin[3];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdv2 = 2.0/dxv[2]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[2]; 
  confFac[0] = 1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = 1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*m_; 

  double dfVfac_l = vmap_prime[0]*(0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3]); 
  double dfVfac_r = vmap_prime[0]*(1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2]); 

  double fVfac_l = (jacobvel_skin[0]*(0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3]))/vmap_primeSq; 
  double fVfac_r = (jacobvel_skin[0]*(1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2]))/vmap_primeSq; 

  if (edge == -1) { 

  double phaseFacl[12] = {0.0}; 

  phaseFacl[3] = 1.5*fskin_over_jacv[3]*fVfac_l*rdv2Sq-0.8660254037844386*fskin_over_jacv[0]*fVfac_l*rdv2Sq; 
  phaseFacl[5] = 1.5*fskin_over_jacv[5]*fVfac_l*rdv2Sq-0.8660254037844386*fskin_over_jacv[1]*fVfac_l*rdv2Sq; 
  phaseFacl[6] = 1.5*fskin_over_jacv[6]*fVfac_l*rdv2Sq-0.8660254037844386*fskin_over_jacv[2]*fVfac_l*rdv2Sq; 
  phaseFacl[7] = 1.5*fskin_over_jacv[7]*fVfac_l*rdv2Sq-0.8660254037844386*fskin_over_jacv[4]*fVfac_l*rdv2Sq; 
  phaseFacl[10] = 1.5*fskin_over_jacv[10]*fVfac_l*rdv2Sq-0.8660254037844387*fskin_over_jacv[8]*fVfac_l*rdv2Sq; 
  phaseFacl[11] = 1.5*fskin_over_jacv[11]*fVfac_l*rdv2Sq-0.8660254037844387*fskin_over_jacv[9]*fVfac_l*rdv2Sq; 

  double phaseFacr[12] = {0.0}; 
  const double dv_edgeR2 = pow(dv_edge,2);
  const double dv_edgeR3 = pow(dv_edge,3);
  const double dv_edgeR4 = pow(dv_edge,4);
  const double dv_skinR2 = pow(dv_skin,2);
  const double dv_skinR3 = pow(dv_skin,3);
  const double dv_skinR4 = pow(dv_skin,4);

  phaseFacr[0] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[3]*dfVfac_r*rdv2)-8.660254037844386*fedge_over_jacv[3]*dfVfac_r*rdv2-9.0*fskin_over_jacv[0]*dfVfac_r*rdv2+9.0*fedge_over_jacv[0]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[3]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[3]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[3]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[3]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[1] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[5]*dfVfac_r*rdv2)-8.660254037844386*fedge_over_jacv[5]*dfVfac_r*rdv2-9.0*fskin_over_jacv[1]*dfVfac_r*rdv2+9.0*fedge_over_jacv[1]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[5]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[5]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[5]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[5]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[2] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[6]*dfVfac_r*rdv2)-8.660254037844386*fedge_over_jacv[6]*dfVfac_r*rdv2-9.0*fskin_over_jacv[2]*dfVfac_r*rdv2+9.0*fedge_over_jacv[2]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[6]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[6]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[6]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[6]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[3] = (dv_edgeR3*((-3.0*fskin_over_jacv[3]*fVfac_r*rdv2Sq)-1.732050807568877*fskin_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-5.0*fskin_over_jacv[3]*fVfac_r*rdv2Sq)-5.196152422706631*fskin_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*(5.0*fedge_over_jacv[3]*fVfac_r*rdv2Sq-5.196152422706631*fedge_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*(3.0*fedge_over_jacv[3]*fVfac_r*rdv2Sq-1.732050807568877*fedge_over_jacv[0]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-15.0*fskin_over_jacv[3]*dfVfac_r*rdv2)-15.0*fedge_over_jacv[3]*dfVfac_r*rdv2-15.58845726811989*fskin_over_jacv[0]*dfVfac_r*rdv2+15.58845726811989*fedge_over_jacv[0]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[3]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[3]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[3]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[3]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[4] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[7]*dfVfac_r*rdv2)-8.660254037844386*fedge_over_jacv[7]*dfVfac_r*rdv2-9.0*fskin_over_jacv[4]*dfVfac_r*rdv2+9.0*fedge_over_jacv[4]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[7]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[7]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[7]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[7]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[5] = (dv_edgeR3*((-3.0*fskin_over_jacv[5]*fVfac_r*rdv2Sq)-1.732050807568877*fskin_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-5.0*fskin_over_jacv[5]*fVfac_r*rdv2Sq)-5.196152422706631*fskin_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*(5.0*fedge_over_jacv[5]*fVfac_r*rdv2Sq-5.196152422706631*fedge_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*(3.0*fedge_over_jacv[5]*fVfac_r*rdv2Sq-1.732050807568877*fedge_over_jacv[1]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-15.0*fskin_over_jacv[5]*dfVfac_r*rdv2)-15.0*fedge_over_jacv[5]*dfVfac_r*rdv2-15.58845726811989*fskin_over_jacv[1]*dfVfac_r*rdv2+15.58845726811989*fedge_over_jacv[1]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[5]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[5]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[5]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[5]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[6] = (dv_edgeR3*((-3.0*fskin_over_jacv[6]*fVfac_r*rdv2Sq)-1.732050807568877*fskin_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-5.0*fskin_over_jacv[6]*fVfac_r*rdv2Sq)-5.196152422706631*fskin_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*(5.0*fedge_over_jacv[6]*fVfac_r*rdv2Sq-5.196152422706631*fedge_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*(3.0*fedge_over_jacv[6]*fVfac_r*rdv2Sq-1.732050807568877*fedge_over_jacv[2]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-15.0*fskin_over_jacv[6]*dfVfac_r*rdv2)-15.0*fedge_over_jacv[6]*dfVfac_r*rdv2-15.58845726811989*fskin_over_jacv[2]*dfVfac_r*rdv2+15.58845726811989*fedge_over_jacv[2]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[6]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[6]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[6]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[6]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[7] = (dv_edgeR3*((-3.0*fskin_over_jacv[7]*fVfac_r*rdv2Sq)-1.732050807568877*fskin_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-5.0*fskin_over_jacv[7]*fVfac_r*rdv2Sq)-5.196152422706631*fskin_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*(5.0*fedge_over_jacv[7]*fVfac_r*rdv2Sq-5.196152422706631*fedge_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*(3.0*fedge_over_jacv[7]*fVfac_r*rdv2Sq-1.732050807568877*fedge_over_jacv[4]*fVfac_r*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-15.0*fskin_over_jacv[7]*dfVfac_r*rdv2)-15.0*fedge_over_jacv[7]*dfVfac_r*rdv2-15.58845726811989*fskin_over_jacv[4]*dfVfac_r*rdv2+15.58845726811989*fedge_over_jacv[4]*dfVfac_r*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[7]*dfVfac_r*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[7]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[7]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[7]*dfVfac_r*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacr[8] = (dv_edgeR2*dv_skinR2*((-43.30127018922195*fskin_over_jacv[10]*dfVfac_r*rdv2)-43.30127018922195*fedge_over_jacv[10]*dfVfac_r*rdv2-45.0*fskin_over_jacv[8]*dfVfac_r*rdv2+45.0*fedge_over_jacv[8]*dfVfac_r*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fedge_over_jacv[10]*dfVfac_r*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fedge_over_jacv[10]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fskin_over_jacv[10]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fskin_over_jacv[10]*dfVfac_r*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacr[9] = (dv_edgeR2*dv_skinR2*((-43.30127018922195*fskin_over_jacv[11]*dfVfac_r*rdv2)-43.30127018922195*fedge_over_jacv[11]*dfVfac_r*rdv2-45.0*fskin_over_jacv[9]*dfVfac_r*rdv2+45.0*fedge_over_jacv[9]*dfVfac_r*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fedge_over_jacv[11]*dfVfac_r*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fedge_over_jacv[11]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fskin_over_jacv[11]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fskin_over_jacv[11]*dfVfac_r*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacr[10] = (dv_edgeR3*((-15.0*fskin_over_jacv[10]*fVfac_r*rdv2Sq)-8.660254037844387*fskin_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-25.0*fskin_over_jacv[10]*fVfac_r*rdv2Sq)-25.98076211353316*fskin_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edge*dv_skinR2*(25.0*fedge_over_jacv[10]*fVfac_r*rdv2Sq-25.98076211353316*fedge_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_skinR3*(15.0*fedge_over_jacv[10]*fVfac_r*rdv2Sq-8.660254037844387*fedge_over_jacv[8]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-75.0*fskin_over_jacv[10]*dfVfac_r*rdv2)-75.0*fedge_over_jacv[10]*dfVfac_r*rdv2-77.94228634059948*fskin_over_jacv[8]*dfVfac_r*rdv2+77.94228634059948*fedge_over_jacv[8]*dfVfac_r*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fedge_over_jacv[10]*dfVfac_r*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fedge_over_jacv[10]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fskin_over_jacv[10]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fskin_over_jacv[10]*dfVfac_r*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacr[11] = (dv_edgeR3*((-15.0*fskin_over_jacv[11]*fVfac_r*rdv2Sq)-8.660254037844387*fskin_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skin*((-25.0*fskin_over_jacv[11]*fVfac_r*rdv2Sq)-25.98076211353316*fskin_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edge*dv_skinR2*(25.0*fedge_over_jacv[11]*fVfac_r*rdv2Sq-25.98076211353316*fedge_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_skinR3*(15.0*fedge_over_jacv[11]*fVfac_r*rdv2Sq-8.660254037844387*fedge_over_jacv[9]*fVfac_r*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*((-75.0*fskin_over_jacv[11]*dfVfac_r*rdv2)-75.0*fedge_over_jacv[11]*dfVfac_r*rdv2-77.94228634059948*fskin_over_jacv[9]*dfVfac_r*rdv2+77.94228634059948*fedge_over_jacv[9]*dfVfac_r*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fedge_over_jacv[11]*dfVfac_r*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fedge_over_jacv[11]*dfVfac_r*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fskin_over_jacv[11]*dfVfac_r*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fskin_over_jacv[11]*dfVfac_r*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 

  double incrl[12] = {0.0}; 
  incrl[3] = 0.7071067811865475*confFac[1]*phaseFacl[5]+0.7071067811865475*confFac[0]*phaseFacl[3]; 
  incrl[5] = 0.7071067811865475*confFac[0]*phaseFacl[5]+0.7071067811865475*confFac[1]*phaseFacl[3]; 
  incrl[6] = 0.7071067811865475*confFac[1]*phaseFacl[7]+0.7071067811865475*confFac[0]*phaseFacl[6]; 
  incrl[7] = 0.7071067811865475*confFac[0]*phaseFacl[7]+0.7071067811865475*confFac[1]*phaseFacl[6]; 
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

  out[0] += incrr[0]; 
  out[1] += incrr[1]; 
  out[2] += incrr[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += incrr[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += incrr[8]; 
  out[9] += incrr[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 


  } else { 

  double phaseFacl[12] = {0.0}; 
  const double dv_edgeR2 = pow(dv_edge,2);
  const double dv_edgeR3 = pow(dv_edge,3);
  const double dv_edgeR4 = pow(dv_edge,4);
  const double dv_skinR2 = pow(dv_skin,2);
  const double dv_skinR3 = pow(dv_skin,3);
  const double dv_skinR4 = pow(dv_skin,4);

  phaseFacl[0] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[3]*dfVfac_l*rdv2)-8.660254037844386*fedge_over_jacv[3]*dfVfac_l*rdv2+9.0*fskin_over_jacv[0]*dfVfac_l*rdv2-9.0*fedge_over_jacv[0]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[3]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[3]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[3]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[3]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[1] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[5]*dfVfac_l*rdv2)-8.660254037844386*fedge_over_jacv[5]*dfVfac_l*rdv2+9.0*fskin_over_jacv[1]*dfVfac_l*rdv2-9.0*fedge_over_jacv[1]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[5]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[5]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[5]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[5]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[2] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[6]*dfVfac_l*rdv2)-8.660254037844386*fedge_over_jacv[6]*dfVfac_l*rdv2+9.0*fskin_over_jacv[2]*dfVfac_l*rdv2-9.0*fedge_over_jacv[2]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[6]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[6]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[6]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[6]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[3] = (dv_edgeR2*dv_skin*(5.0*fskin_over_jacv[3]*fVfac_l*rdv2Sq-5.196152422706631*fskin_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR3*(3.0*fskin_over_jacv[3]*fVfac_l*rdv2Sq-1.732050807568877*fskin_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*((-3.0*fedge_over_jacv[3]*fVfac_l*rdv2Sq)-1.732050807568877*fedge_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-5.0*fedge_over_jacv[3]*fVfac_l*rdv2Sq)-5.196152422706631*fedge_over_jacv[0]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(15.0*fskin_over_jacv[3]*dfVfac_l*rdv2+15.0*fedge_over_jacv[3]*dfVfac_l*rdv2-15.58845726811989*fskin_over_jacv[0]*dfVfac_l*rdv2+15.58845726811989*fedge_over_jacv[0]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[3]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[3]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[3]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[3]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[4] = (dv_edgeR2*dv_skinR2*((-8.660254037844386*fskin_over_jacv[7]*dfVfac_l*rdv2)-8.660254037844386*fedge_over_jacv[7]*dfVfac_l*rdv2+9.0*fskin_over_jacv[4]*dfVfac_l*rdv2-9.0*fedge_over_jacv[4]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fedge_over_jacv[7]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fedge_over_jacv[7]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(1.732050807568877*fskin_over_jacv[7]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(1.732050807568877*fskin_over_jacv[7]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[5] = (dv_edgeR2*dv_skin*(5.0*fskin_over_jacv[5]*fVfac_l*rdv2Sq-5.196152422706631*fskin_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR3*(3.0*fskin_over_jacv[5]*fVfac_l*rdv2Sq-1.732050807568877*fskin_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*((-3.0*fedge_over_jacv[5]*fVfac_l*rdv2Sq)-1.732050807568877*fedge_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-5.0*fedge_over_jacv[5]*fVfac_l*rdv2Sq)-5.196152422706631*fedge_over_jacv[1]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(15.0*fskin_over_jacv[5]*dfVfac_l*rdv2+15.0*fedge_over_jacv[5]*dfVfac_l*rdv2-15.58845726811989*fskin_over_jacv[1]*dfVfac_l*rdv2+15.58845726811989*fedge_over_jacv[1]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[5]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[5]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[5]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[5]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[6] = (dv_edgeR2*dv_skin*(5.0*fskin_over_jacv[6]*fVfac_l*rdv2Sq-5.196152422706631*fskin_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR3*(3.0*fskin_over_jacv[6]*fVfac_l*rdv2Sq-1.732050807568877*fskin_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*((-3.0*fedge_over_jacv[6]*fVfac_l*rdv2Sq)-1.732050807568877*fedge_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-5.0*fedge_over_jacv[6]*fVfac_l*rdv2Sq)-5.196152422706631*fedge_over_jacv[2]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(15.0*fskin_over_jacv[6]*dfVfac_l*rdv2+15.0*fedge_over_jacv[6]*dfVfac_l*rdv2-15.58845726811989*fskin_over_jacv[2]*dfVfac_l*rdv2+15.58845726811989*fedge_over_jacv[2]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[6]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[6]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[6]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[6]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[7] = (dv_edgeR2*dv_skin*(5.0*fskin_over_jacv[7]*fVfac_l*rdv2Sq-5.196152422706631*fskin_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR3*(3.0*fskin_over_jacv[7]*fVfac_l*rdv2Sq-1.732050807568877*fskin_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_skinR3*((-3.0*fedge_over_jacv[7]*fVfac_l*rdv2Sq)-1.732050807568877*fedge_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-5.0*fedge_over_jacv[7]*fVfac_l*rdv2Sq)-5.196152422706631*fedge_over_jacv[4]*fVfac_l*rdv2Sq))/(2.0*dv_skinR3+6.0*dv_edge*dv_skinR2+6.0*dv_edgeR2*dv_skin+2.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(15.0*fskin_over_jacv[7]*dfVfac_l*rdv2+15.0*fedge_over_jacv[7]*dfVfac_l*rdv2-15.58845726811989*fskin_over_jacv[4]*dfVfac_l*rdv2+15.58845726811989*fedge_over_jacv[4]*dfVfac_l*rdv2))/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fedge_over_jacv[7]*dfVfac_l*dv_skinR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fedge_over_jacv[7]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)+(3.0*fskin_over_jacv[7]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin)-(3.0*fskin_over_jacv[7]*dfVfac_l*dv_edgeR4*rdv2)/(dv_edge*dv_skinR4+3.0*dv_edgeR2*dv_skinR3+3.0*dv_edgeR3*dv_skinR2+dv_edgeR4*dv_skin); 
  phaseFacl[8] = (dv_edgeR2*dv_skinR2*((-43.30127018922195*fskin_over_jacv[10]*dfVfac_l*rdv2)-43.30127018922195*fedge_over_jacv[10]*dfVfac_l*rdv2+45.0*fskin_over_jacv[8]*dfVfac_l*rdv2-45.0*fedge_over_jacv[8]*dfVfac_l*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fedge_over_jacv[10]*dfVfac_l*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fedge_over_jacv[10]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fskin_over_jacv[10]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fskin_over_jacv[10]*dfVfac_l*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacl[9] = (dv_edgeR2*dv_skinR2*((-43.30127018922195*fskin_over_jacv[11]*dfVfac_l*rdv2)-43.30127018922195*fedge_over_jacv[11]*dfVfac_l*rdv2+45.0*fskin_over_jacv[9]*dfVfac_l*rdv2-45.0*fedge_over_jacv[9]*dfVfac_l*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fedge_over_jacv[11]*dfVfac_l*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fedge_over_jacv[11]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(8.660254037844387*fskin_over_jacv[11]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(8.660254037844387*fskin_over_jacv[11]*dfVfac_l*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacl[10] = (dv_edgeR2*dv_skin*(25.0*fskin_over_jacv[10]*fVfac_l*rdv2Sq-25.98076211353316*fskin_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR3*(15.0*fskin_over_jacv[10]*fVfac_l*rdv2Sq-8.660254037844387*fskin_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_skinR3*((-15.0*fedge_over_jacv[10]*fVfac_l*rdv2Sq)-8.660254037844387*fedge_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-25.0*fedge_over_jacv[10]*fVfac_l*rdv2Sq)-25.98076211353316*fedge_over_jacv[8]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(75.0*fskin_over_jacv[10]*dfVfac_l*rdv2+75.0*fedge_over_jacv[10]*dfVfac_l*rdv2-77.94228634059948*fskin_over_jacv[8]*dfVfac_l*rdv2+77.94228634059948*fedge_over_jacv[8]*dfVfac_l*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fedge_over_jacv[10]*dfVfac_l*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fedge_over_jacv[10]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fskin_over_jacv[10]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fskin_over_jacv[10]*dfVfac_l*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 
  phaseFacl[11] = (dv_edgeR2*dv_skin*(25.0*fskin_over_jacv[11]*fVfac_l*rdv2Sq-25.98076211353316*fskin_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR3*(15.0*fskin_over_jacv[11]*fVfac_l*rdv2Sq-8.660254037844387*fskin_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_skinR3*((-15.0*fedge_over_jacv[11]*fVfac_l*rdv2Sq)-8.660254037844387*fedge_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edge*dv_skinR2*((-25.0*fedge_over_jacv[11]*fVfac_l*rdv2Sq)-25.98076211353316*fedge_over_jacv[9]*fVfac_l*rdv2Sq))/(10.0*dv_skinR3+30.0*dv_edge*dv_skinR2+30.0*dv_edgeR2*dv_skin+10.0*dv_edgeR3)+(dv_edgeR2*dv_skinR2*(75.0*fskin_over_jacv[11]*dfVfac_l*rdv2+75.0*fedge_over_jacv[11]*dfVfac_l*rdv2-77.94228634059948*fskin_over_jacv[9]*dfVfac_l*rdv2+77.94228634059948*fedge_over_jacv[9]*dfVfac_l*rdv2))/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fedge_over_jacv[11]*dfVfac_l*dv_skinR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fedge_over_jacv[11]*dfVfac_l*dv_edge*dv_skinR3*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)+(15.0*fskin_over_jacv[11]*dfVfac_l*dv_edgeR3*dv_skin*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin)-(15.0*fskin_over_jacv[11]*dfVfac_l*dv_edgeR4*rdv2)/(5.0*dv_edge*dv_skinR4+15.0*dv_edgeR2*dv_skinR3+15.0*dv_edgeR3*dv_skinR2+5.0*dv_edgeR4*dv_skin); 

  double phaseFacr[12] = {0.0}; 

  phaseFacr[3] = (-1.5*fskin_over_jacv[3]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin_over_jacv[0]*fVfac_r*rdv2Sq; 
  phaseFacr[5] = (-1.5*fskin_over_jacv[5]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin_over_jacv[1]*fVfac_r*rdv2Sq; 
  phaseFacr[6] = (-1.5*fskin_over_jacv[6]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin_over_jacv[2]*fVfac_r*rdv2Sq; 
  phaseFacr[7] = (-1.5*fskin_over_jacv[7]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin_over_jacv[4]*fVfac_r*rdv2Sq; 
  phaseFacr[10] = (-1.5*fskin_over_jacv[10]*fVfac_r*rdv2Sq)-0.8660254037844387*fskin_over_jacv[8]*fVfac_r*rdv2Sq; 
  phaseFacr[11] = (-1.5*fskin_over_jacv[11]*fVfac_r*rdv2Sq)-0.8660254037844387*fskin_over_jacv[9]*fVfac_r*rdv2Sq; 

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
  incrr[3] = 0.7071067811865475*confFac[1]*phaseFacr[5]+0.7071067811865475*confFac[0]*phaseFacr[3]; 
  incrr[5] = 0.7071067811865475*confFac[0]*phaseFacr[5]+0.7071067811865475*confFac[1]*phaseFacr[3]; 
  incrr[6] = 0.7071067811865475*confFac[1]*phaseFacr[7]+0.7071067811865475*confFac[0]*phaseFacr[6]; 
  incrr[7] = 0.7071067811865475*confFac[0]*phaseFacr[7]+0.7071067811865475*confFac[1]*phaseFacr[6]; 
  incrr[10] = 0.7071067811865475*confFac[1]*phaseFacr[11]+0.7071067811865475*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.7071067811865475*confFac[0]*phaseFacr[11]+0.7071067811865475*confFac[1]*phaseFacr[10]; 

  out[0] += -1.0*incrl[0]; 
  out[1] += -1.0*incrl[1]; 
  out[2] += -1.0*incrl[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += -1.0*incrl[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += -1.0*incrl[8]; 
  out[9] += -1.0*incrl[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 

  } 

  return 0.;

} 

GKYL_CU_DH double lbo_gyrokinetic_diff_notmapped_boundary_surfmu_1x2v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel_edge, const double *jacobvel_skin, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // dxv: Cell spacing. 
  // vmap_edge,vmap_skin: velocity space mapping.
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // jacobvel_edge,jacobvel_skin: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdv2 = 2.0/dxv[2]; 
  double rdv2Sq = rdv2*rdv2; 

  double confFac[2]; 
  confFac[0] = 1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*m_; 
  confFac[1] = 1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*m_; 

  double dfVfac_l = 0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3]; 
  double dfVfac_r = 1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2]; 

  double fVfac_l = 0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3]; 
  double fVfac_r = 1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2]; 

  if (edge == -1) { 

  double phaseFacl[12] = {0.0}; 

  phaseFacl[3] = 1.5*fskin[3]*fVfac_l*rdv2Sq-0.8660254037844386*fskin[0]*fVfac_l*rdv2Sq; 
  phaseFacl[5] = 1.5*fskin[5]*fVfac_l*rdv2Sq-0.8660254037844386*fskin[1]*fVfac_l*rdv2Sq; 
  phaseFacl[6] = 1.5*fskin[6]*fVfac_l*rdv2Sq-0.8660254037844386*fskin[2]*fVfac_l*rdv2Sq; 
  phaseFacl[7] = 1.5*fskin[7]*fVfac_l*rdv2Sq-0.8660254037844386*fskin[4]*fVfac_l*rdv2Sq; 
  phaseFacl[10] = 1.5*fskin[10]*fVfac_l*rdv2Sq-0.8660254037844387*fskin[8]*fVfac_l*rdv2Sq; 
  phaseFacl[11] = 1.5*fskin[11]*fVfac_l*rdv2Sq-0.8660254037844387*fskin[9]*fVfac_l*rdv2Sq; 

  double phaseFacr[12] = {0.0}; 

  phaseFacr[0] = (-0.5412658773652741*fskin[3]*dfVfac_r*rdv2Sq)-0.5412658773652741*fedge[3]*dfVfac_r*rdv2Sq-0.5625*fskin[0]*dfVfac_r*rdv2Sq+0.5625*fedge[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[1] = (-0.5412658773652741*fskin[5]*dfVfac_r*rdv2Sq)-0.5412658773652741*fedge[5]*dfVfac_r*rdv2Sq-0.5625*fskin[1]*dfVfac_r*rdv2Sq+0.5625*fedge[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[2] = (-0.5412658773652741*fskin[6]*dfVfac_r*rdv2Sq)-0.5412658773652741*fedge[6]*dfVfac_r*rdv2Sq-0.5625*fskin[2]*dfVfac_r*rdv2Sq+0.5625*fedge[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[3] = (-0.5*fskin[3]*fVfac_r*rdv2Sq)+0.5*fedge[3]*fVfac_r*rdv2Sq-0.4330127018922193*fskin[0]*fVfac_r*rdv2Sq-0.4330127018922193*fedge[0]*fVfac_r*rdv2Sq-0.9375*fskin[3]*dfVfac_r*rdv2Sq-0.9375*fedge[3]*dfVfac_r*rdv2Sq-0.9742785792574932*fskin[0]*dfVfac_r*rdv2Sq+0.9742785792574932*fedge[0]*dfVfac_r*rdv2Sq; 
  phaseFacr[4] = (-0.5412658773652741*fskin[7]*dfVfac_r*rdv2Sq)-0.5412658773652741*fedge[7]*dfVfac_r*rdv2Sq-0.5625*fskin[4]*dfVfac_r*rdv2Sq+0.5625*fedge[4]*dfVfac_r*rdv2Sq; 
  phaseFacr[5] = (-0.5*fskin[5]*fVfac_r*rdv2Sq)+0.5*fedge[5]*fVfac_r*rdv2Sq-0.4330127018922193*fskin[1]*fVfac_r*rdv2Sq-0.4330127018922193*fedge[1]*fVfac_r*rdv2Sq-0.9375*fskin[5]*dfVfac_r*rdv2Sq-0.9375*fedge[5]*dfVfac_r*rdv2Sq-0.9742785792574932*fskin[1]*dfVfac_r*rdv2Sq+0.9742785792574932*fedge[1]*dfVfac_r*rdv2Sq; 
  phaseFacr[6] = (-0.5*fskin[6]*fVfac_r*rdv2Sq)+0.5*fedge[6]*fVfac_r*rdv2Sq-0.4330127018922193*fskin[2]*fVfac_r*rdv2Sq-0.4330127018922193*fedge[2]*fVfac_r*rdv2Sq-0.9375*fskin[6]*dfVfac_r*rdv2Sq-0.9375*fedge[6]*dfVfac_r*rdv2Sq-0.9742785792574932*fskin[2]*dfVfac_r*rdv2Sq+0.9742785792574932*fedge[2]*dfVfac_r*rdv2Sq; 
  phaseFacr[7] = (-0.5*fskin[7]*fVfac_r*rdv2Sq)+0.5*fedge[7]*fVfac_r*rdv2Sq-0.4330127018922193*fskin[4]*fVfac_r*rdv2Sq-0.4330127018922193*fedge[4]*fVfac_r*rdv2Sq-0.9375*fskin[7]*dfVfac_r*rdv2Sq-0.9375*fedge[7]*dfVfac_r*rdv2Sq-0.9742785792574932*fskin[4]*dfVfac_r*rdv2Sq+0.9742785792574932*fedge[4]*dfVfac_r*rdv2Sq; 
  phaseFacr[8] = (-0.5412658773652742*fskin[10]*dfVfac_r*rdv2Sq)-0.5412658773652742*fedge[10]*dfVfac_r*rdv2Sq-0.5625*fskin[8]*dfVfac_r*rdv2Sq+0.5625*fedge[8]*dfVfac_r*rdv2Sq; 
  phaseFacr[9] = (-0.5412658773652742*fskin[11]*dfVfac_r*rdv2Sq)-0.5412658773652742*fedge[11]*dfVfac_r*rdv2Sq-0.5625*fskin[9]*dfVfac_r*rdv2Sq+0.5625*fedge[9]*dfVfac_r*rdv2Sq; 
  phaseFacr[10] = (-0.5*fskin[10]*fVfac_r*rdv2Sq)+0.5*fedge[10]*fVfac_r*rdv2Sq-0.4330127018922194*fskin[8]*fVfac_r*rdv2Sq-0.4330127018922194*fedge[8]*fVfac_r*rdv2Sq-0.9375*fskin[10]*dfVfac_r*rdv2Sq-0.9375*fedge[10]*dfVfac_r*rdv2Sq-0.9742785792574935*fskin[8]*dfVfac_r*rdv2Sq+0.9742785792574935*fedge[8]*dfVfac_r*rdv2Sq; 
  phaseFacr[11] = (-0.5*fskin[11]*fVfac_r*rdv2Sq)+0.5*fedge[11]*fVfac_r*rdv2Sq-0.4330127018922194*fskin[9]*fVfac_r*rdv2Sq-0.4330127018922194*fedge[9]*fVfac_r*rdv2Sq-0.9375*fskin[11]*dfVfac_r*rdv2Sq-0.9375*fedge[11]*dfVfac_r*rdv2Sq-0.9742785792574935*fskin[9]*dfVfac_r*rdv2Sq+0.9742785792574935*fedge[9]*dfVfac_r*rdv2Sq; 

  double incrl[12] = {0.0}; 
  incrl[3] = 0.7071067811865475*confFac[1]*phaseFacl[5]+0.7071067811865475*confFac[0]*phaseFacl[3]; 
  incrl[5] = 0.7071067811865475*confFac[0]*phaseFacl[5]+0.7071067811865475*confFac[1]*phaseFacl[3]; 
  incrl[6] = 0.7071067811865475*confFac[1]*phaseFacl[7]+0.7071067811865475*confFac[0]*phaseFacl[6]; 
  incrl[7] = 0.7071067811865475*confFac[0]*phaseFacl[7]+0.7071067811865475*confFac[1]*phaseFacl[6]; 
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

  out[0] += incrr[0]; 
  out[1] += incrr[1]; 
  out[2] += incrr[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += incrr[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += incrr[8]; 
  out[9] += incrr[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 


  } else { 

  double phaseFacl[12] = {0.0}; 

  phaseFacl[0] = (-0.5412658773652741*fskin[3]*dfVfac_l*rdv2Sq)-0.5412658773652741*fedge[3]*dfVfac_l*rdv2Sq+0.5625*fskin[0]*dfVfac_l*rdv2Sq-0.5625*fedge[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[1] = (-0.5412658773652741*fskin[5]*dfVfac_l*rdv2Sq)-0.5412658773652741*fedge[5]*dfVfac_l*rdv2Sq+0.5625*fskin[1]*dfVfac_l*rdv2Sq-0.5625*fedge[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[2] = (-0.5412658773652741*fskin[6]*dfVfac_l*rdv2Sq)-0.5412658773652741*fedge[6]*dfVfac_l*rdv2Sq+0.5625*fskin[2]*dfVfac_l*rdv2Sq-0.5625*fedge[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[3] = 0.5*fskin[3]*fVfac_l*rdv2Sq-0.5*fedge[3]*fVfac_l*rdv2Sq-0.4330127018922193*fskin[0]*fVfac_l*rdv2Sq-0.4330127018922193*fedge[0]*fVfac_l*rdv2Sq+0.9375*fskin[3]*dfVfac_l*rdv2Sq+0.9375*fedge[3]*dfVfac_l*rdv2Sq-0.9742785792574932*fskin[0]*dfVfac_l*rdv2Sq+0.9742785792574932*fedge[0]*dfVfac_l*rdv2Sq; 
  phaseFacl[4] = (-0.5412658773652741*fskin[7]*dfVfac_l*rdv2Sq)-0.5412658773652741*fedge[7]*dfVfac_l*rdv2Sq+0.5625*fskin[4]*dfVfac_l*rdv2Sq-0.5625*fedge[4]*dfVfac_l*rdv2Sq; 
  phaseFacl[5] = 0.5*fskin[5]*fVfac_l*rdv2Sq-0.5*fedge[5]*fVfac_l*rdv2Sq-0.4330127018922193*fskin[1]*fVfac_l*rdv2Sq-0.4330127018922193*fedge[1]*fVfac_l*rdv2Sq+0.9375*fskin[5]*dfVfac_l*rdv2Sq+0.9375*fedge[5]*dfVfac_l*rdv2Sq-0.9742785792574932*fskin[1]*dfVfac_l*rdv2Sq+0.9742785792574932*fedge[1]*dfVfac_l*rdv2Sq; 
  phaseFacl[6] = 0.5*fskin[6]*fVfac_l*rdv2Sq-0.5*fedge[6]*fVfac_l*rdv2Sq-0.4330127018922193*fskin[2]*fVfac_l*rdv2Sq-0.4330127018922193*fedge[2]*fVfac_l*rdv2Sq+0.9375*fskin[6]*dfVfac_l*rdv2Sq+0.9375*fedge[6]*dfVfac_l*rdv2Sq-0.9742785792574932*fskin[2]*dfVfac_l*rdv2Sq+0.9742785792574932*fedge[2]*dfVfac_l*rdv2Sq; 
  phaseFacl[7] = 0.5*fskin[7]*fVfac_l*rdv2Sq-0.5*fedge[7]*fVfac_l*rdv2Sq-0.4330127018922193*fskin[4]*fVfac_l*rdv2Sq-0.4330127018922193*fedge[4]*fVfac_l*rdv2Sq+0.9375*fskin[7]*dfVfac_l*rdv2Sq+0.9375*fedge[7]*dfVfac_l*rdv2Sq-0.9742785792574932*fskin[4]*dfVfac_l*rdv2Sq+0.9742785792574932*fedge[4]*dfVfac_l*rdv2Sq; 
  phaseFacl[8] = (-0.5412658773652742*fskin[10]*dfVfac_l*rdv2Sq)-0.5412658773652742*fedge[10]*dfVfac_l*rdv2Sq+0.5625*fskin[8]*dfVfac_l*rdv2Sq-0.5625*fedge[8]*dfVfac_l*rdv2Sq; 
  phaseFacl[9] = (-0.5412658773652742*fskin[11]*dfVfac_l*rdv2Sq)-0.5412658773652742*fedge[11]*dfVfac_l*rdv2Sq+0.5625*fskin[9]*dfVfac_l*rdv2Sq-0.5625*fedge[9]*dfVfac_l*rdv2Sq; 
  phaseFacl[10] = 0.5*fskin[10]*fVfac_l*rdv2Sq-0.5*fedge[10]*fVfac_l*rdv2Sq-0.4330127018922194*fskin[8]*fVfac_l*rdv2Sq-0.4330127018922194*fedge[8]*fVfac_l*rdv2Sq+0.9375*fskin[10]*dfVfac_l*rdv2Sq+0.9375*fedge[10]*dfVfac_l*rdv2Sq-0.9742785792574935*fskin[8]*dfVfac_l*rdv2Sq+0.9742785792574935*fedge[8]*dfVfac_l*rdv2Sq; 
  phaseFacl[11] = 0.5*fskin[11]*fVfac_l*rdv2Sq-0.5*fedge[11]*fVfac_l*rdv2Sq-0.4330127018922194*fskin[9]*fVfac_l*rdv2Sq-0.4330127018922194*fedge[9]*fVfac_l*rdv2Sq+0.9375*fskin[11]*dfVfac_l*rdv2Sq+0.9375*fedge[11]*dfVfac_l*rdv2Sq-0.9742785792574935*fskin[9]*dfVfac_l*rdv2Sq+0.9742785792574935*fedge[9]*dfVfac_l*rdv2Sq; 

  double phaseFacr[12] = {0.0}; 

  phaseFacr[3] = (-1.5*fskin[3]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin[0]*fVfac_r*rdv2Sq; 
  phaseFacr[5] = (-1.5*fskin[5]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin[1]*fVfac_r*rdv2Sq; 
  phaseFacr[6] = (-1.5*fskin[6]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin[2]*fVfac_r*rdv2Sq; 
  phaseFacr[7] = (-1.5*fskin[7]*fVfac_r*rdv2Sq)-0.8660254037844386*fskin[4]*fVfac_r*rdv2Sq; 
  phaseFacr[10] = (-1.5*fskin[10]*fVfac_r*rdv2Sq)-0.8660254037844387*fskin[8]*fVfac_r*rdv2Sq; 
  phaseFacr[11] = (-1.5*fskin[11]*fVfac_r*rdv2Sq)-0.8660254037844387*fskin[9]*fVfac_r*rdv2Sq; 

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
  incrr[3] = 0.7071067811865475*confFac[1]*phaseFacr[5]+0.7071067811865475*confFac[0]*phaseFacr[3]; 
  incrr[5] = 0.7071067811865475*confFac[0]*phaseFacr[5]+0.7071067811865475*confFac[1]*phaseFacr[3]; 
  incrr[6] = 0.7071067811865475*confFac[1]*phaseFacr[7]+0.7071067811865475*confFac[0]*phaseFacr[6]; 
  incrr[7] = 0.7071067811865475*confFac[0]*phaseFacr[7]+0.7071067811865475*confFac[1]*phaseFacr[6]; 
  incrr[10] = 0.7071067811865475*confFac[1]*phaseFacr[11]+0.7071067811865475*confFac[0]*phaseFacr[10]; 
  incrr[11] = 0.7071067811865475*confFac[0]*phaseFacr[11]+0.7071067811865475*confFac[1]*phaseFacr[10]; 

  out[0] += -1.0*incrl[0]; 
  out[1] += -1.0*incrl[1]; 
  out[2] += -1.0*incrl[2]; 
  out[3] += incrr[3]-1.0*incrl[3]; 
  out[4] += -1.0*incrl[4]; 
  out[5] += incrr[5]-1.0*incrl[5]; 
  out[6] += incrr[6]-1.0*incrl[6]; 
  out[7] += incrr[7]-1.0*incrl[7]; 
  out[8] += -1.0*incrl[8]; 
  out[9] += -1.0*incrl[9]; 
  out[10] += incrr[10]-1.0*incrl[10]; 
  out[11] += incrr[11]-1.0*incrl[11]; 

  } 

  return 0.;

} 
