#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_tensor_5x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_tensor_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
   double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvx = w[2];
  double rdvx2 = 2.0/dxv[2];
  double wvy = w[3];
  double rdvy2 = 2.0/dxv[3];
  double wvz = w[4];
  double rdvz2 = 2.0/dxv[4];

  double *alphaL = &alpha_surf[16];
  double *sgn_alpha_surfL = &sgn_alpha_surf[16];
  alphaL[0] = 1.224744871391589*hamil[4]*rdvy2-2.121320343559642*hamil[10]*rdvy2; 
  alphaL[1] = 1.224744871391589*hamil[9]*rdvy2-2.121320343559642*hamil[17]*rdvy2; 
  alphaL[2] = 1.224744871391589*hamil[11]*rdvy2-2.121320343559642*hamil[19]*rdvy2; 
  alphaL[4] = 1.224744871391589*hamil[15]*rdvy2-2.121320343559642*hamil[24]*rdvy2; 
  alphaL[5] = 1.224744871391589*hamil[18]*rdvy2-2.121320343559642*hamil[26]*rdvy2; 
  alphaL[8] = 1.224744871391589*hamil[23]*rdvy2-2.121320343559642*hamil[28]*rdvy2; 
  alphaL[9] = 1.224744871391589*hamil[25]*rdvy2-2.121320343559642*hamil[30]*rdvy2; 
  alphaL[12] = 1.224744871391589*hamil[29]*rdvy2-2.121320343559642*hamil[31]*rdvy2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.25*alphaL[12])+0.25*(alphaL[9]+alphaL[8]+alphaL[5])-0.25*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.25*alphaL[12]-0.25*(alphaL[9]+alphaL[8])+0.25*(alphaL[5]+alphaL[4])-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL[12])+0.25*(alphaL[9]+alphaL[8]+alphaL[5])-0.25*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL[12]-0.25*(alphaL[9]+alphaL[8])+0.25*(alphaL[5]+alphaL[4])-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL[12]-0.25*alphaL[9]+0.25*alphaL[8]-0.25*(alphaL[5]+alphaL[4])+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL[12])+0.25*alphaL[9]-0.25*(alphaL[8]+alphaL[5])+0.25*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL[12]-0.25*alphaL[9]+0.25*alphaL[8]-0.25*(alphaL[5]+alphaL[4])+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL[12])+0.25*alphaL[9]-0.25*(alphaL[8]+alphaL[5])+0.25*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[12]+alphaL[9])-0.25*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL[12]+alphaL[9]))+0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[12]+alphaL[9])-0.25*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL[12]+alphaL[9]))+0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL[12]+alphaL[9]+alphaL[8]))+0.25*alphaL[5]-0.25*alphaL[4]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[12]+alphaL[9]+alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL[12]+alphaL[9]+alphaL[8]))+0.25*alphaL[5]-0.25*alphaL[4]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[12]+alphaL[9]+alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
