#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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

  double *alphaR = &alpha_surf[24];
  double *sgn_alpha_surfR = &sgn_alpha_surf[24];
  alphaR[0] = ((0.3423265984407287*jacobtot_inv[1]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*b_x[3]*hamil[32]-0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[16]-0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[16]-0.2651650429449552*b_x[5]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[8]-0.2651650429449552*b_x[1]*jacobtot_inv[1]*hamil[8]-0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[8]-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[7]-0.1530931089239486*hamil[3]*b_x[5]*jacobtot_inv[5]-0.1530931089239486*b_x[3]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[3])*rdz2+(vmap[0]*(0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[4])*rdz2)/vmap[1]+((-0.3423265984407287*jacobtot_inv[3]*b_z[5]*hamil[32])-0.3423265984407287*jacobtot_inv[0]*b_z[1]*hamil[32]+0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[16]+0.2651650429449552*b_z[0]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_z[3]*hamil[16]+0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[7]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[7]+0.2651650429449552*b_z[5]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_z[3]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[6]+0.2651650429449552*b_z[0]*jacobtot_inv[0]*hamil[6]+0.1530931089239486*hamil[1]*b_z[5]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*b_z[3]*jacobtot_inv[3]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[1]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[1])*rdx2+(vmap[0]*((-0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_z[5])-0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[4])*rdx2)/vmap[1])/q_; 
  alphaR[1] = ((0.3423265984407287*jacobtot_inv[0]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[1]*b_x[3]*hamil[32]-0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[16]-0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[16]-0.4772970773009194*b_x[1]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[16]-0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[8]-0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[8]-0.2755675960631073*b_x[5]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[7]-0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[7]-0.1530931089239486*b_x[3]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*hamil[3]*jacobtot_inv[3]*b_x[5]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[3])*rdz2+(vmap[0]*(0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[4])*rdz2)/vmap[1]+((-0.3423265984407287*b_z[5]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[1]*hamil[32]+0.2651650429449552*b_z[0]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[16]+0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_z[3]*hamil[16]+0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[7]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[7]+0.2651650429449552*b_z[3]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[3]*b_z[5]*hamil[6]+0.2651650429449552*b_z[0]*jacobtot_inv[1]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[6]+0.1530931089239486*hamil[1]*b_z[3]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*jacobtot_inv[3]*b_z[5]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[1]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[1])*rdx2+(vmap[0]*((-0.1530931089239486*hamil[4]*b_z[5]*jacobtot_inv[5])-0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[4])*rdx2)/vmap[1])/q_; 
  alphaR[2] = ((0.3423265984407287*b_x[5]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*b_x[3]*jacobtot_inv[3]*hamil[32]-0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[16]-0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[16]-0.2651650429449552*b_x[1]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*jacobtot_inv[1]*b_x[5]*hamil[8]-0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[8]-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[7]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[7]-0.1530931089239486*b_x[1]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[1]*hamil[3]*b_x[5]-0.1530931089239486*b_x[0]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[3])*rdz2+(vmap[0]*(0.1530931089239486*hamil[4]*b_x[5]*jacobtot_inv[5]+0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[4])*rdz2)/vmap[1]+((-0.3423265984407287*jacobtot_inv[0]*b_z[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[3]*hamil[32]+0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[16]+0.4772970773009194*b_z[3]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*b_z[0]*jacobtot_inv[0]*hamil[16]+0.2755675960631073*b_z[5]*jacobtot_inv[5]*hamil[7]+0.2755675960631073*b_z[3]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[7]+0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[6]+0.2651650429449552*b_z[0]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_z[3]*hamil[6]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*jacobtot_inv[1]*b_z[5]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[3]+0.1530931089239486*jacobtot_inv[0]*hamil[1]*b_z[3])*rdx2+(vmap[0]*((-0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_z[5])-0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[4])*rdx2)/vmap[1])/q_; 
  alphaR[3] = ((vmap[0]*(0.3423265984407287*jacobtot_inv[1]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*b_x[3]*hamil[32])*rdz2)/vmap[1]+(0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[4])*rdz2+(vmap[0]*((-0.3423265984407287*jacobtot_inv[3]*b_z[5]*hamil[32])-0.3423265984407287*jacobtot_inv[0]*b_z[1]*hamil[32])*rdx2)/vmap[1]+((-0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_z[5])-0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[4])*rdx2)/q_; 
  alphaR[4] = (((-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[21]-0.1530931089239486*b_x[5]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[14])*rdz2+(0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[21]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[21]+0.1530931089239486*b_z[5]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*b_z[3]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[12])*rdx2)/q_; 
  alphaR[5] = ((0.3423265984407287*b_x[3]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*jacobtot_inv[3]*b_x[5]*hamil[32]-0.4772970773009194*b_x[1]*jacobtot_inv[5]*hamil[16]-0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[16]-0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[8]-0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[8]-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[8]-0.2755675960631073*b_x[1]*jacobtot_inv[5]*hamil[7]-0.2755675960631073*jacobtot_inv[1]*b_x[5]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[7]-0.1530931089239486*b_x[0]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[0]*hamil[3]*b_x[5]-0.1530931089239486*b_x[1]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[3])*rdz2+(vmap[0]*(0.1530931089239486*b_x[3]*hamil[4]*jacobtot_inv[5]+0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_x[5])*rdz2)/vmap[1]+((-0.3423265984407287*b_z[1]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*jacobtot_inv[1]*b_z[5]*hamil[32]+0.4772970773009194*b_z[3]*jacobtot_inv[5]*hamil[16]+0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[16]+0.2651650429449552*b_z[0]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[16]+0.2755675960631073*b_z[3]*jacobtot_inv[5]*hamil[7]+0.2755675960631073*jacobtot_inv[3]*b_z[5]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[7]+0.2651650429449552*b_z[0]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[6]+0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_z[3]*hamil[6]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[5]+0.1530931089239486*jacobtot_inv[0]*hamil[1]*b_z[5]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[3]+0.1530931089239486*hamil[1]*jacobtot_inv[1]*b_z[3])*rdx2+(vmap[0]*((-0.1530931089239486*b_z[1]*hamil[4]*jacobtot_inv[5])-0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_z[5])*rdx2)/vmap[1])/q_; 
  alphaR[6] = ((vmap[0]*(0.3423265984407287*jacobtot_inv[0]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[1]*b_x[3]*hamil[32])*rdz2)/vmap[1]+(0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[4])*rdz2+(vmap[0]*((-0.3423265984407287*b_z[5]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[1]*hamil[32])*rdx2)/vmap[1]+((-0.1530931089239486*hamil[4]*b_z[5]*jacobtot_inv[5])-0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[4])*rdx2)/q_; 
  alphaR[7] = ((vmap[0]*(0.3423265984407287*b_x[5]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*b_x[3]*jacobtot_inv[3]*hamil[32])*rdz2)/vmap[1]+(0.1530931089239486*hamil[4]*b_x[5]*jacobtot_inv[5]+0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[4])*rdz2+(vmap[0]*((-0.3423265984407287*jacobtot_inv[0]*b_z[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[3]*hamil[32])*rdx2)/vmap[1]+((-0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_z[5])-0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[4])*rdx2)/q_; 
  alphaR[8] = (((-0.2755675960631073*b_x[5]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[21]-0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[21]-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[14])*rdz2+(0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[21]+0.1530931089239486*b_z[3]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[3]*b_z[5]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[12])*rdx2)/q_; 
  alphaR[9] = (((-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[21]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[21]-0.1530931089239486*b_x[1]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*b_x[5]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[14])*rdz2+(0.2755675960631073*b_z[5]*jacobtot_inv[5]*hamil[21]+0.2755675960631073*b_z[3]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[12])*rdx2)/q_; 
  alphaR[11] = ((vmap[0]*(0.3423265984407287*b_x[3]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*jacobtot_inv[3]*b_x[5]*hamil[32])*rdz2)/vmap[1]+(0.1530931089239486*b_x[3]*hamil[4]*jacobtot_inv[5]+0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_x[5])*rdz2+(vmap[0]*((-0.3423265984407287*b_z[1]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*jacobtot_inv[1]*b_z[5]*hamil[32])*rdx2)/vmap[1]+((-0.1530931089239486*b_z[1]*hamil[4]*jacobtot_inv[5])-0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_z[5])*rdx2)/q_; 
  alphaR[12] = (((-0.2755675960631073*b_x[1]*jacobtot_inv[5]*hamil[21])-0.2755675960631073*jacobtot_inv[1]*b_x[5]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[14]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[14])*rdz2+(0.2755675960631073*b_z[3]*jacobtot_inv[5]*hamil[21]+0.2755675960631073*jacobtot_inv[3]*b_z[5]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[12]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[12])*rdx2)/q_; 
  alphaR[16] = ((0.3061862178478971*jacobtot_inv[1]*b_x[5]*hamil[32]+0.3061862178478971*jacobtot_inv[0]*b_x[3]*hamil[32])*rdz2+((-0.3061862178478971*jacobtot_inv[3]*b_z[5]*hamil[32])-0.3061862178478971*jacobtot_inv[0]*b_z[1]*hamil[32])*rdx2)/q_; 
  alphaR[17] = ((0.3061862178478971*jacobtot_inv[0]*b_x[5]*hamil[32]+0.3061862178478971*jacobtot_inv[1]*b_x[3]*hamil[32])*rdz2+((-0.3061862178478971*b_z[5]*jacobtot_inv[5]*hamil[32])-0.3061862178478971*b_z[1]*jacobtot_inv[1]*hamil[32])*rdx2)/q_; 
  alphaR[18] = ((0.3061862178478971*b_x[5]*jacobtot_inv[5]*hamil[32]+0.3061862178478971*b_x[3]*jacobtot_inv[3]*hamil[32])*rdz2+((-0.3061862178478971*jacobtot_inv[0]*b_z[5]*hamil[32])-0.3061862178478971*b_z[1]*jacobtot_inv[3]*hamil[32])*rdx2)/q_; 
  alphaR[20] = ((0.3061862178478971*b_x[3]*jacobtot_inv[5]*hamil[32]+0.3061862178478971*jacobtot_inv[3]*b_x[5]*hamil[32])*rdz2+((-0.3061862178478971*b_z[1]*jacobtot_inv[5]*hamil[32])-0.3061862178478971*jacobtot_inv[1]*b_z[5]*hamil[32])*rdx2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.2236067977499786*alphaR[20]-0.2236067977499786*(alphaR[18]+alphaR[17])+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*(alphaR[9]+alphaR[8])+0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.2795084971874732*alphaR[20])+0.2795084971874732*(alphaR[18]+alphaR[17])-0.2795084971874732*alphaR[16]-0.25*alphaR[12]+0.25*(alphaR[9]+alphaR[8]+alphaR[5])-0.25*(alphaR[4]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR[20]-0.2236067977499786*(alphaR[18]+alphaR[17])+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*(alphaR[9]+alphaR[8])-0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR[20]-0.2236067977499786*(alphaR[18]+alphaR[17])+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*(alphaR[9]+alphaR[8])+0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[3]-0.25*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaR[20])+0.2795084971874732*(alphaR[18]+alphaR[17])-0.2795084971874732*alphaR[16]+0.25*alphaR[12]-0.25*(alphaR[9]+alphaR[8])+0.25*(alphaR[5]+alphaR[4])-0.25*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR[20]-0.2236067977499786*(alphaR[18]+alphaR[17])+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*(alphaR[9]+alphaR[8])-0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[3]-0.25*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR[20]-0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.25*alphaR[12]-0.25*alphaR[9]+0.25*alphaR[8]-0.25*(alphaR[5]+alphaR[4])+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR[20]-0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[12]+0.25*alphaR[9]-0.25*(alphaR[8]+alphaR[5])+0.25*(alphaR[4]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaR[20]+alphaR[18]))+0.2236067977499786*(alphaR[17]+alphaR[16])+0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaR[20]+alphaR[18])-0.2795084971874732*(alphaR[17]+alphaR[16])+0.25*(alphaR[12]+alphaR[9])-0.25*(alphaR[8]+alphaR[5]+alphaR[4]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaR[20]+alphaR[18]))+0.2236067977499786*(alphaR[17]+alphaR[16])+0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaR[20]+alphaR[18]))+0.2236067977499786*(alphaR[17]+alphaR[16])-0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaR[20]+alphaR[18])-0.2795084971874732*(alphaR[17]+alphaR[16])-0.25*(alphaR[12]+alphaR[9])+0.25*alphaR[8]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaR[20]+alphaR[18]))+0.2236067977499786*(alphaR[17]+alphaR[16])-0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16])-0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*(alphaR[9]+alphaR[8])-0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16]))-0.25*(alphaR[12]+alphaR[9]+alphaR[8])+0.25*alphaR[5]-0.25*alphaR[4]+0.25*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16])-0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*(alphaR[9]+alphaR[8])+0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16])+0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*(alphaR[9]+alphaR[8])-0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[3]+0.25*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[12]+alphaR[9]+alphaR[8]+alphaR[5]+alphaR[4]+alphaR[2]+alphaR[1]+alphaR[0])-0.2795084971874732*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16]) > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[20]+alphaR[18]+alphaR[17]+alphaR[16])+0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*(alphaR[9]+alphaR[8])+0.3354101966249678*(alphaR[7]+alphaR[6])+0.25*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[3]+0.25*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
