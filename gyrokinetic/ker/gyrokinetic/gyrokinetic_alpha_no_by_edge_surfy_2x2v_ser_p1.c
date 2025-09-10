#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

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

  double *alphaR = &alpha_surf[12];
  double *sgn_alpha_surfR = &sgn_alpha_surf[12];
  alphaR[0] = (0.75*cmag[3]*hamil[3]*jacobtot_inv[3]+0.4330127018922193*cmag[1]*hamil[3]*jacobtot_inv[3]+0.4330127018922193*jacobtot_inv[1]*cmag[3]*hamil[3]+0.75*cmag[2]*jacobtot_inv[2]*hamil[3]+0.4330127018922193*cmag[0]*jacobtot_inv[2]*hamil[3]+0.4330127018922193*jacobtot_inv[0]*cmag[2]*hamil[3]+0.25*cmag[1]*jacobtot_inv[1]*hamil[3]+0.25*cmag[0]*jacobtot_inv[0]*hamil[3])/(vmap[1]*m_); 
  alphaR[1] = (0.75*cmag[2]*hamil[3]*jacobtot_inv[3]+0.4330127018922193*cmag[0]*hamil[3]*jacobtot_inv[3]+0.75*jacobtot_inv[2]*cmag[3]*hamil[3]+0.4330127018922193*jacobtot_inv[0]*cmag[3]*hamil[3]+0.4330127018922193*cmag[1]*jacobtot_inv[2]*hamil[3]+0.4330127018922193*jacobtot_inv[1]*cmag[2]*hamil[3]+0.25*cmag[0]*jacobtot_inv[1]*hamil[3]+0.25*jacobtot_inv[0]*cmag[1]*hamil[3])/(vmap[1]*m_); 
  alphaR[2] = (1.677050983124842*cmag[3]*jacobtot_inv[3]*hamil[16]+0.9682458365518543*cmag[1]*jacobtot_inv[3]*hamil[16]+0.9682458365518543*jacobtot_inv[1]*cmag[3]*hamil[16]+1.677050983124842*cmag[2]*jacobtot_inv[2]*hamil[16]+0.9682458365518543*cmag[0]*jacobtot_inv[2]*hamil[16]+0.9682458365518543*jacobtot_inv[0]*cmag[2]*hamil[16]+0.5590169943749475*cmag[1]*jacobtot_inv[1]*hamil[16]+0.5590169943749475*cmag[0]*jacobtot_inv[0]*hamil[16])/(vmap[1]*m_); 
  alphaR[4] = (1.677050983124842*cmag[2]*jacobtot_inv[3]*hamil[16]+0.9682458365518543*cmag[0]*jacobtot_inv[3]*hamil[16]+1.677050983124842*jacobtot_inv[2]*cmag[3]*hamil[16]+0.9682458365518543*jacobtot_inv[0]*cmag[3]*hamil[16]+0.9682458365518543*cmag[1]*jacobtot_inv[2]*hamil[16]+0.9682458365518543*jacobtot_inv[1]*cmag[2]*hamil[16]+0.5590169943749475*cmag[0]*jacobtot_inv[1]*hamil[16]+0.5590169943749475*jacobtot_inv[0]*cmag[1]*hamil[16])/(vmap[1]*m_); 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4743416490252568*alphaR[4]-0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.3535533905932734*alphaR[0]-0.3535533905932734*alphaR[1] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252568*alphaR[4])+0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaR[1]+alphaR[0])-0.4743416490252568*(alphaR[4]+alphaR[2]) > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*(alphaR[4]+alphaR[2])+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
