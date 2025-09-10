#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfz_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *cmag, const double *b_i, const double *phi, double* GKYL_RESTRICT alpha_surf,
    double* GKYL_RESTRICT sgn_alpha_surf) 
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
  // phi: electrostatic potential.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double rdx2 = 2.0/dxv[0];
  double rdy2 = 2.0/dxv[1];
  double rdz2 = 2.0/dxv[2];
  double rdvpar2 = 2.0/dxv[3];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*(phi[0]*q_+vmapSq[0]*m_)+1.414213562373095*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[4] = 2.0*vmapSq[1]*m_; 
  hamil[5] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_+1.414213562373095*vmap[2]*bmag[5]; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[14] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = 1.414213562373095*vmap[3]*bmag[5]; 
  hamil[32] = 2.0*vmapSq[2]*m_; 

  double *alphaL = &alpha_surf[48];
  double *sgn_alpha_surfL = &sgn_alpha_surf[48];
  alphaL[0] = (((-0.7954951288348656*b_x[3]*jacobtot_inv[5]*hamil[16])+0.4592793267718456*b_x[0]*jacobtot_inv[5]*hamil[16]-0.7954951288348656*jacobtot_inv[3]*b_x[5]*hamil[16]+0.4592793267718456*jacobtot_inv[0]*b_x[5]*hamil[16]+0.4592793267718456*b_x[1]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*jacobtot_inv[1]*b_x[3]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[16]-0.7954951288348656*b_x[5]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*b_x[1]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*jacobtot_inv[1]*b_x[5]*hamil[8]-0.7954951288348656*b_x[3]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*b_x[0]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*jacobtot_inv[0]*b_x[3]*hamil[8]-0.2651650429449552*b_x[1]*jacobtot_inv[1]*hamil[8]-0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[8]+0.4592793267718456*b_x[3]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[6]+0.4592793267718456*jacobtot_inv[3]*b_x[5]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[6]-0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[6]+0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[6]+0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[6]+0.4592793267718456*hamil[2]*b_x[5]*jacobtot_inv[5]-0.2651650429449552*b_x[1]*hamil[2]*jacobtot_inv[5]-0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[5]+0.4592793267718456*hamil[2]*b_x[3]*jacobtot_inv[3]-0.2651650429449552*b_x[0]*hamil[2]*jacobtot_inv[3]-0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[3]+0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[2]+0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[2])*rdy2+(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[32]-0.592927061281571*jacobtot_inv[0]*b_y[5]*hamil[32]-0.592927061281571*b_y[1]*jacobtot_inv[3]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*b_y[1]*hamil[32]+0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[7]-0.4592793267718456*b_y[1]*jacobtot_inv[5]*hamil[7]-0.4592793267718456*jacobtot_inv[1]*b_y[5]*hamil[7]+0.7954951288348656*b_y[3]*jacobtot_inv[3]*hamil[7]-0.4592793267718456*b_y[0]*jacobtot_inv[3]*hamil[7]-0.4592793267718456*jacobtot_inv[0]*b_y[3]*hamil[7]+0.2651650429449552*b_y[1]*jacobtot_inv[1]*hamil[7]+0.2651650429449552*b_y[0]*jacobtot_inv[0]*hamil[7]-0.4592793267718456*hamil[1]*b_y[5]*jacobtot_inv[5]+0.2651650429449552*b_y[1]*hamil[1]*jacobtot_inv[5]+0.2651650429449552*hamil[1]*jacobtot_inv[1]*b_y[5]-0.4592793267718456*hamil[1]*b_y[3]*jacobtot_inv[3]+0.2651650429449552*b_y[0]*hamil[1]*jacobtot_inv[3]+0.2651650429449552*jacobtot_inv[0]*hamil[1]*b_y[3]-0.1530931089239486*b_y[1]*hamil[1]*jacobtot_inv[1]-0.1530931089239486*b_y[0]*jacobtot_inv[0]*hamil[1])*rdx2+(vmap[0]*(0.4592793267718456*jacobtot_inv[3]*hamil[4]*b_y[5]-0.2651650429449552*jacobtot_inv[0]*hamil[4]*b_y[5]-0.2651650429449552*b_y[1]*jacobtot_inv[3]*hamil[4]+0.1530931089239486*jacobtot_inv[0]*b_y[1]*hamil[4])*rdx2)/vmap[1])/q_+(0.375*hamil[4]*cmag[5]*jacobtot_inv[5]-0.2165063509461096*cmag[1]*hamil[4]*jacobtot_inv[5]-0.2165063509461096*jacobtot_inv[1]*hamil[4]*cmag[5]+0.375*cmag[3]*jacobtot_inv[3]*hamil[4]-0.2165063509461096*cmag[0]*jacobtot_inv[3]*hamil[4]-0.2165063509461096*jacobtot_inv[0]*cmag[3]*hamil[4]+0.125*cmag[1]*jacobtot_inv[1]*hamil[4]+0.125*cmag[0]*jacobtot_inv[0]*hamil[4])/(vmap[1]*m_); 
  alphaL[1] = (((-1.431891231902758*b_x[5]*jacobtot_inv[5]*hamil[16])+0.826702788189322*b_x[1]*jacobtot_inv[5]*hamil[16]+0.826702788189322*jacobtot_inv[1]*b_x[5]*hamil[16]-0.7954951288348656*b_x[3]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*b_x[0]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*jacobtot_inv[0]*b_x[3]*hamil[16]-0.4772970773009194*b_x[1]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[16]-0.7954951288348656*b_x[3]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*b_x[0]*jacobtot_inv[5]*hamil[8]-0.7954951288348656*jacobtot_inv[3]*b_x[5]*hamil[8]+0.4592793267718456*jacobtot_inv[0]*b_x[5]*hamil[8]+0.4592793267718456*b_x[1]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*jacobtot_inv[1]*b_x[3]*hamil[8]-0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[8]+0.826702788189322*b_x[5]*jacobtot_inv[5]*hamil[6]-0.4772970773009194*b_x[1]*jacobtot_inv[5]*hamil[6]-0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[6]+0.4592793267718456*b_x[3]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[6]+0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[6]+0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[6]+0.4592793267718456*hamil[2]*b_x[3]*jacobtot_inv[5]-0.2651650429449552*b_x[0]*hamil[2]*jacobtot_inv[5]+0.4592793267718456*hamil[2]*jacobtot_inv[3]*b_x[5]-0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[5]-0.2651650429449552*b_x[1]*hamil[2]*jacobtot_inv[3]-0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[3]+0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[2]+0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[2])*rdy2+(1.026979795322186*b_y[5]*jacobtot_inv[5]*hamil[32]-0.592927061281571*b_y[1]*jacobtot_inv[5]*hamil[32]-0.592927061281571*jacobtot_inv[1]*b_y[5]*hamil[32]+0.3423265984407287*b_y[1]*jacobtot_inv[1]*hamil[32]+0.7954951288348656*b_y[3]*jacobtot_inv[5]*hamil[7]-0.4592793267718456*b_y[0]*jacobtot_inv[5]*hamil[7]+0.7954951288348656*jacobtot_inv[3]*b_y[5]*hamil[7]-0.4592793267718456*jacobtot_inv[0]*b_y[5]*hamil[7]-0.4592793267718456*b_y[1]*jacobtot_inv[3]*hamil[7]-0.4592793267718456*jacobtot_inv[1]*b_y[3]*hamil[7]+0.2651650429449552*b_y[0]*jacobtot_inv[1]*hamil[7]+0.2651650429449552*jacobtot_inv[0]*b_y[1]*hamil[7]-0.4592793267718456*hamil[1]*b_y[3]*jacobtot_inv[5]+0.2651650429449552*b_y[0]*hamil[1]*jacobtot_inv[5]-0.4592793267718456*hamil[1]*jacobtot_inv[3]*b_y[5]+0.2651650429449552*jacobtot_inv[0]*hamil[1]*b_y[5]+0.2651650429449552*b_y[1]*hamil[1]*jacobtot_inv[3]+0.2651650429449552*hamil[1]*jacobtot_inv[1]*b_y[3]-0.1530931089239486*b_y[0]*hamil[1]*jacobtot_inv[1]-0.1530931089239486*jacobtot_inv[0]*b_y[1]*hamil[1])*rdx2+(vmap[0]*(0.4592793267718456*hamil[4]*b_y[5]*jacobtot_inv[5]-0.2651650429449552*b_y[1]*hamil[4]*jacobtot_inv[5]-0.2651650429449552*jacobtot_inv[1]*hamil[4]*b_y[5]+0.1530931089239486*b_y[1]*jacobtot_inv[1]*hamil[4])*rdx2)/vmap[1])/q_+(0.375*cmag[3]*hamil[4]*jacobtot_inv[5]-0.2165063509461096*cmag[0]*hamil[4]*jacobtot_inv[5]+0.375*jacobtot_inv[3]*hamil[4]*cmag[5]-0.2165063509461096*jacobtot_inv[0]*hamil[4]*cmag[5]-0.2165063509461096*cmag[1]*jacobtot_inv[3]*hamil[4]-0.2165063509461096*jacobtot_inv[1]*cmag[3]*hamil[4]+0.125*cmag[0]*jacobtot_inv[1]*hamil[4]+0.125*jacobtot_inv[0]*cmag[1]*hamil[4])/(vmap[1]*m_); 
  alphaL[2] = ((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[16]-0.4592793267718456*b_y[1]*jacobtot_inv[5]*hamil[16]-0.4592793267718456*jacobtot_inv[1]*b_y[5]*hamil[16]+0.7954951288348656*b_y[3]*jacobtot_inv[3]*hamil[16]-0.4592793267718456*b_y[0]*jacobtot_inv[3]*hamil[16]-0.4592793267718456*jacobtot_inv[0]*b_y[3]*hamil[16]+0.2651650429449552*b_y[1]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*b_y[0]*jacobtot_inv[0]*hamil[16]-0.4592793267718456*b_y[5]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_y[1]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_y[5]*hamil[6]-0.4592793267718456*b_y[3]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*b_y[0]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_y[3]*hamil[6]-0.1530931089239486*b_y[1]*jacobtot_inv[1]*hamil[6]-0.1530931089239486*b_y[0]*jacobtot_inv[0]*hamil[6])*rdx2)/q_; 
  alphaL[3] = ((vmap[0]*(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[32]-0.592927061281571*jacobtot_inv[0]*b_y[5]*hamil[32]-0.592927061281571*b_y[1]*jacobtot_inv[3]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*b_y[1]*hamil[32])*rdx2)/vmap[1]+(0.4592793267718456*jacobtot_inv[3]*hamil[4]*b_y[5]-0.2651650429449552*jacobtot_inv[0]*hamil[4]*b_y[5]-0.2651650429449552*b_y[1]*jacobtot_inv[3]*hamil[4]+0.1530931089239486*jacobtot_inv[0]*b_y[1]*hamil[4])*rdx2)/q_+(0.8385254915624212*cmag[5]*jacobtot_inv[5]*hamil[32]-0.4841229182759271*cmag[1]*jacobtot_inv[5]*hamil[32]-0.4841229182759271*jacobtot_inv[1]*cmag[5]*hamil[32]+0.8385254915624212*cmag[3]*jacobtot_inv[3]*hamil[32]-0.4841229182759271*cmag[0]*jacobtot_inv[3]*hamil[32]-0.4841229182759271*jacobtot_inv[0]*cmag[3]*hamil[32]+0.2795084971874737*cmag[1]*jacobtot_inv[1]*hamil[32]+0.2795084971874737*cmag[0]*jacobtot_inv[0]*hamil[32])/(vmap[1]*m_); 
  alphaL[4] = ((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[21]-0.4592793267718456*b_y[1]*jacobtot_inv[5]*hamil[21]-0.4592793267718456*jacobtot_inv[1]*b_y[5]*hamil[21]+0.7954951288348656*b_y[3]*jacobtot_inv[3]*hamil[21]-0.4592793267718456*b_y[0]*jacobtot_inv[3]*hamil[21]-0.4592793267718456*jacobtot_inv[0]*b_y[3]*hamil[21]+0.2651650429449552*b_y[1]*jacobtot_inv[1]*hamil[21]+0.2651650429449552*b_y[0]*jacobtot_inv[0]*hamil[21]-0.4592793267718456*b_y[5]*jacobtot_inv[5]*hamil[12]+0.2651650429449552*b_y[1]*jacobtot_inv[5]*hamil[12]+0.2651650429449552*jacobtot_inv[1]*b_y[5]*hamil[12]-0.4592793267718456*b_y[3]*jacobtot_inv[3]*hamil[12]+0.2651650429449552*b_y[0]*jacobtot_inv[3]*hamil[12]+0.2651650429449552*jacobtot_inv[0]*b_y[3]*hamil[12]-0.1530931089239486*b_y[1]*jacobtot_inv[1]*hamil[12]-0.1530931089239486*b_y[0]*jacobtot_inv[0]*hamil[12])*rdx2)/q_; 
  alphaL[5] = ((0.7954951288348656*b_y[3]*jacobtot_inv[5]*hamil[16]-0.4592793267718456*b_y[0]*jacobtot_inv[5]*hamil[16]+0.7954951288348656*jacobtot_inv[3]*b_y[5]*hamil[16]-0.4592793267718456*jacobtot_inv[0]*b_y[5]*hamil[16]-0.4592793267718456*b_y[1]*jacobtot_inv[3]*hamil[16]-0.4592793267718456*jacobtot_inv[1]*b_y[3]*hamil[16]+0.2651650429449552*b_y[0]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_y[1]*hamil[16]-0.4592793267718456*b_y[3]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_y[0]*jacobtot_inv[5]*hamil[6]-0.4592793267718456*jacobtot_inv[3]*b_y[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_y[5]*hamil[6]+0.2651650429449552*b_y[1]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_y[3]*hamil[6]-0.1530931089239486*b_y[0]*jacobtot_inv[1]*hamil[6]-0.1530931089239486*jacobtot_inv[0]*b_y[1]*hamil[6])*rdx2)/q_; 
  alphaL[6] = ((vmap[0]*(1.026979795322186*b_y[5]*jacobtot_inv[5]*hamil[32]-0.592927061281571*b_y[1]*jacobtot_inv[5]*hamil[32]-0.592927061281571*jacobtot_inv[1]*b_y[5]*hamil[32]+0.3423265984407287*b_y[1]*jacobtot_inv[1]*hamil[32])*rdx2)/vmap[1]+(0.4592793267718456*hamil[4]*b_y[5]*jacobtot_inv[5]-0.2651650429449552*b_y[1]*hamil[4]*jacobtot_inv[5]-0.2651650429449552*jacobtot_inv[1]*hamil[4]*b_y[5]+0.1530931089239486*b_y[1]*jacobtot_inv[1]*hamil[4])*rdx2)/q_+(0.8385254915624212*cmag[3]*jacobtot_inv[5]*hamil[32]-0.4841229182759271*cmag[0]*jacobtot_inv[5]*hamil[32]+0.8385254915624212*jacobtot_inv[3]*cmag[5]*hamil[32]-0.4841229182759271*jacobtot_inv[0]*cmag[5]*hamil[32]-0.4841229182759271*cmag[1]*jacobtot_inv[3]*hamil[32]-0.4841229182759271*jacobtot_inv[1]*cmag[3]*hamil[32]+0.2795084971874737*cmag[0]*jacobtot_inv[1]*hamil[32]+0.2795084971874737*jacobtot_inv[0]*cmag[1]*hamil[32])/(vmap[1]*m_); 
  alphaL[8] = ((0.7954951288348656*b_y[3]*jacobtot_inv[5]*hamil[21]-0.4592793267718456*b_y[0]*jacobtot_inv[5]*hamil[21]+0.7954951288348656*jacobtot_inv[3]*b_y[5]*hamil[21]-0.4592793267718456*jacobtot_inv[0]*b_y[5]*hamil[21]-0.4592793267718456*b_y[1]*jacobtot_inv[3]*hamil[21]-0.4592793267718456*jacobtot_inv[1]*b_y[3]*hamil[21]+0.2651650429449552*b_y[0]*jacobtot_inv[1]*hamil[21]+0.2651650429449552*jacobtot_inv[0]*b_y[1]*hamil[21]-0.4592793267718456*b_y[3]*jacobtot_inv[5]*hamil[12]+0.2651650429449552*b_y[0]*jacobtot_inv[5]*hamil[12]-0.4592793267718456*jacobtot_inv[3]*b_y[5]*hamil[12]+0.2651650429449552*jacobtot_inv[0]*b_y[5]*hamil[12]+0.2651650429449552*b_y[1]*jacobtot_inv[3]*hamil[12]+0.2651650429449552*jacobtot_inv[1]*b_y[3]*hamil[12]-0.1530931089239486*b_y[0]*jacobtot_inv[1]*hamil[12]-0.1530931089239486*jacobtot_inv[0]*b_y[1]*hamil[12])*rdx2)/q_; 
  alphaL[16] = ((0.9185586535436913*jacobtot_inv[3]*b_y[5]*hamil[32]-0.5303300858899105*jacobtot_inv[0]*b_y[5]*hamil[32]-0.5303300858899105*b_y[1]*jacobtot_inv[3]*hamil[32]+0.3061862178478971*jacobtot_inv[0]*b_y[1]*hamil[32])*rdx2)/q_; 
  alphaL[17] = ((0.9185586535436916*b_y[5]*jacobtot_inv[5]*hamil[32]-0.5303300858899104*b_y[1]*jacobtot_inv[5]*hamil[32]-0.5303300858899104*jacobtot_inv[1]*b_y[5]*hamil[32]+0.3061862178478971*b_y[1]*jacobtot_inv[1]*hamil[32])*rdx2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]+0.25*alphaL[8]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*(alphaL[8]+alphaL[5])-0.25*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]+0.25*alphaL[8]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]-0.25*alphaL[8]+0.3354101966249678*alphaL[6]+0.25*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[8]+0.25*(alphaL[5]+alphaL[4])-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]-0.25*alphaL[8]-0.3354101966249678*alphaL[6]+0.25*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]+0.25*alphaL[8]+0.3354101966249678*alphaL[6]-0.25*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[8]-0.25*(alphaL[5]+alphaL[4])+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]+0.25*alphaL[8]-0.3354101966249678*alphaL[6]-0.25*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]-0.25*alphaL[8]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*(alphaL[8]+alphaL[5])+0.25*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[17])+0.2236067977499786*alphaL[16]-0.25*alphaL[8]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])-0.25*alphaL[8]-0.3354101966249678*alphaL[6]-0.25*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[17]+alphaL[16]))-0.25*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])-0.25*alphaL[8]+0.3354101966249678*alphaL[6]-0.25*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])+0.25*alphaL[8]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[17]+alphaL[16]))+0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])+0.25*alphaL[8]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])-0.25*alphaL[8]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[17]+alphaL[16]))-0.25*alphaL[8]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])-0.25*alphaL[8]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])+0.25*alphaL[8]-0.3354101966249678*alphaL[6]+0.25*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2]+alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[17]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[17]+alphaL[16])+0.25*alphaL[8]+0.3354101966249678*alphaL[6]+0.25*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
