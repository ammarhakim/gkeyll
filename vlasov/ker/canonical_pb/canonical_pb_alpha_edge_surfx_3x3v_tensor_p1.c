#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_tensor_6x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_tensor_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_3x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvx = w[3];
  double rdvx2 = 2.0/dxv[3];
  double wvy = w[4];
  double rdvy2 = 2.0/dxv[4];
  double wvz = w[5];
  double rdvz2 = 2.0/dxv[5];

  double *alphaR = &alpha_surf[0];
  double *sgn_alpha_surfR = &sgn_alpha_surf[0];
  alphaR[0] = 2.121320343559642*hamil[10]*rdvx2+1.224744871391589*hamil[4]*rdvx2; 
  alphaR[1] = 2.121320343559642*hamil[23]*rdvx2+1.224744871391589*hamil[11]*rdvx2; 
  alphaR[2] = 2.121320343559642*hamil[24]*rdvx2+1.224744871391589*hamil[12]*rdvx2; 
  alphaR[4] = 2.121320343559642*hamil[29]*rdvx2+1.224744871391589*hamil[16]*rdvx2; 
  alphaR[5] = 2.121320343559642*hamil[35]*rdvx2+1.224744871391589*hamil[20]*rdvx2; 
  alphaR[6] = 2.121320343559642*hamil[42]*rdvx2+1.224744871391589*hamil[25]*rdvx2; 
  alphaR[9] = 2.121320343559642*hamil[44]*rdvx2+1.224744871391589*hamil[30]*rdvx2; 
  alphaR[10] = 2.121320343559642*hamil[45]*rdvx2+1.224744871391589*hamil[31]*rdvx2; 
  alphaR[12] = 2.121320343559642*hamil[48]*rdvx2+1.224744871391589*hamil[36]*rdvx2; 
  alphaR[13] = 2.121320343559642*hamil[49]*rdvx2+1.224744871391589*hamil[37]*rdvx2; 
  alphaR[15] = 2.121320343559642*hamil[54]*rdvx2+1.224744871391589*hamil[41]*rdvx2; 
  alphaR[17] = 2.121320343559642*hamil[57]*rdvx2+1.224744871391589*hamil[46]*rdvx2; 
  alphaR[20] = 2.121320343559642*hamil[58]*rdvx2+1.224744871391589*hamil[50]*rdvx2; 
  alphaR[23] = 2.121320343559642*hamil[60]*rdvx2+1.224744871391589*hamil[55]*rdvx2; 
  alphaR[24] = 2.121320343559642*hamil[61]*rdvx2+1.224744871391589*hamil[56]*rdvx2; 
  alphaR[28] = 2.121320343559642*hamil[63]*rdvx2+1.224744871391589*hamil[62]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[20]+alphaR[17])+0.1767766952966367*(alphaR[15]+alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9]+alphaR[6])-0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[20])-0.1767766952966367*(alphaR[17]+alphaR[15]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[6]+alphaR[5])-0.1767766952966367*(alphaR[4]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*(alphaR[24]+alphaR[23])-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[15]+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[10]+alphaR[9])+0.1767766952966367*alphaR[6]-0.1767766952966367*alphaR[5]+0.1767766952966367*alphaR[4]-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[20]+alphaR[17]+alphaR[15])-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[20]+alphaR[17])+0.1767766952966367*(alphaR[15]+alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9]+alphaR[6])-0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*(alphaR[24]+alphaR[23]+alphaR[20])-0.1767766952966367*(alphaR[17]+alphaR[15]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[6]+alphaR[5])-0.1767766952966367*(alphaR[4]+alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*(alphaR[24]+alphaR[23])-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[15]+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[10]+alphaR[9])+0.1767766952966367*alphaR[6]-0.1767766952966367*alphaR[5]+0.1767766952966367*alphaR[4]-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*(alphaR[24]+alphaR[23])+0.1767766952966367*(alphaR[20]+alphaR[17]+alphaR[15])-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])-0.1767766952966367*(alphaR[2]+alphaR[1])+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*alphaR[24]-0.1767766952966367*alphaR[23]+0.1767766952966367*(alphaR[20]+alphaR[17]+alphaR[15])-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[24]+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[15]+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[10])+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[6]+0.1767766952966367*alphaR[5]-0.1767766952966367*alphaR[4]+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[24]+0.1767766952966367*(alphaR[23]+alphaR[20])-0.1767766952966367*(alphaR[17]+alphaR[15]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[10])-0.1767766952966367*(alphaR[9]+alphaR[6]+alphaR[5])+0.1767766952966367*(alphaR[4]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*alphaR[24]-0.1767766952966367*(alphaR[23]+alphaR[20]+alphaR[17])+0.1767766952966367*(alphaR[15]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[6])+0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*alphaR[24]-0.1767766952966367*alphaR[23]+0.1767766952966367*(alphaR[20]+alphaR[17]+alphaR[15])-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*(alphaR[6]+alphaR[5]+alphaR[4])+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[24]+0.1767766952966367*alphaR[23]-0.1767766952966367*alphaR[20]+0.1767766952966367*alphaR[17]-0.1767766952966367*alphaR[15]+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[10])+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[6]+0.1767766952966367*alphaR[5]-0.1767766952966367*alphaR[4]+0.1767766952966367*alphaR[2]-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaR[28]-0.1767766952966367*alphaR[24]+0.1767766952966367*(alphaR[23]+alphaR[20])-0.1767766952966367*(alphaR[17]+alphaR[15]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[10])-0.1767766952966367*(alphaR[9]+alphaR[6]+alphaR[5])+0.1767766952966367*(alphaR[4]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaR[28])+0.1767766952966367*alphaR[24]-0.1767766952966367*(alphaR[23]+alphaR[20]+alphaR[17])+0.1767766952966367*(alphaR[15]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[6])+0.1767766952966367*(alphaR[5]+alphaR[4]+alphaR[2])-0.1767766952966367*alphaR[1]+0.1767766952966367*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]))+0.1767766952966367*(alphaR[23]+alphaR[20]+alphaR[17]+alphaR[15]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24])-0.1767766952966367*(alphaR[23]+alphaR[20])+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[15]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[10])-0.1767766952966367*(alphaR[9]+alphaR[6])+0.1767766952966367*alphaR[5]-0.1767766952966367*(alphaR[4]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24])-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[17]+alphaR[15])+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[10])+0.1767766952966367*alphaR[9]-0.1767766952966367*(alphaR[6]+alphaR[5])+0.1767766952966367*alphaR[4]-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]))+0.1767766952966367*alphaR[23]-0.1767766952966367*(alphaR[20]+alphaR[17])+0.1767766952966367*alphaR[15]-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[6]+0.1767766952966367*(alphaR[5]+alphaR[4])-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]))+0.1767766952966367*(alphaR[23]+alphaR[20]+alphaR[17]+alphaR[15]+alphaR[13])-0.1767766952966367*alphaR[12]+0.1767766952966367*alphaR[10]-0.1767766952966367*(alphaR[9]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24])-0.1767766952966367*(alphaR[23]+alphaR[20])+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[15]+alphaR[13])+0.1767766952966367*(alphaR[12]+alphaR[10])-0.1767766952966367*(alphaR[9]+alphaR[6])+0.1767766952966367*alphaR[5]-0.1767766952966367*(alphaR[4]+alphaR[2])+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24])-0.1767766952966367*alphaR[23]+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[17]+alphaR[15])+0.1767766952966367*alphaR[13]-0.1767766952966367*(alphaR[12]+alphaR[10])+0.1767766952966367*alphaR[9]-0.1767766952966367*(alphaR[6]+alphaR[5])+0.1767766952966367*alphaR[4]-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]))+0.1767766952966367*alphaR[23]-0.1767766952966367*(alphaR[20]+alphaR[17])+0.1767766952966367*alphaR[15]-0.1767766952966367*alphaR[13]+0.1767766952966367*alphaR[12]-0.1767766952966367*alphaR[10]+0.1767766952966367*alphaR[9]-0.1767766952966367*alphaR[6]+0.1767766952966367*(alphaR[5]+alphaR[4])-0.1767766952966367*alphaR[2]+0.1767766952966367*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[20]+alphaR[17])+0.1767766952966367*alphaR[15]-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9])+0.1767766952966367*alphaR[6]-0.1767766952966367*(alphaR[5]+alphaR[4])+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]))+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[17]+alphaR[15])+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[6]+alphaR[5])-0.1767766952966367*alphaR[4]+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]+alphaR[20]))+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[15]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[6])-0.1767766952966367*alphaR[5]+0.1767766952966367*(alphaR[4]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]+alphaR[20]+alphaR[17]+alphaR[15]+alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[27] = 1.0; 
  else  
    sgn_alpha_surfR[27] = -1.0; 
  
  if (sgn_alpha_surfR[27] == sgn_alpha_surfR[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23])-0.1767766952966367*(alphaR[20]+alphaR[17])+0.1767766952966367*alphaR[15]-0.1767766952966367*(alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9])+0.1767766952966367*alphaR[6]-0.1767766952966367*(alphaR[5]+alphaR[4])+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[28] = 1.0; 
  else  
    sgn_alpha_surfR[28] = -1.0; 
  
  if (sgn_alpha_surfR[28] == sgn_alpha_surfR[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]))+0.1767766952966367*alphaR[20]-0.1767766952966367*(alphaR[17]+alphaR[15])+0.1767766952966367*(alphaR[13]+alphaR[12])-0.1767766952966367*(alphaR[10]+alphaR[9])+0.1767766952966367*(alphaR[6]+alphaR[5])-0.1767766952966367*alphaR[4]+0.1767766952966367*(alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[29] = 1.0; 
  else  
    sgn_alpha_surfR[29] = -1.0; 
  
  if (sgn_alpha_surfR[29] == sgn_alpha_surfR[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]+alphaR[20]))+0.1767766952966367*alphaR[17]-0.1767766952966367*(alphaR[15]+alphaR[13]+alphaR[12])+0.1767766952966367*(alphaR[10]+alphaR[9]+alphaR[6])-0.1767766952966367*alphaR[5]+0.1767766952966367*(alphaR[4]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[30] = 1.0; 
  else  
    sgn_alpha_surfR[30] = -1.0; 
  
  if (sgn_alpha_surfR[30] == sgn_alpha_surfR[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaR[28]+alphaR[24]+alphaR[23]+alphaR[20]+alphaR[17]+alphaR[15]+alphaR[13]+alphaR[12]+alphaR[10]+alphaR[9]+alphaR[6]+alphaR[5]+alphaR[4]+alphaR[2]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[31] = 1.0; 
  else  
    sgn_alpha_surfR[31] = -1.0; 
  
  if (sgn_alpha_surfR[31] == sgn_alpha_surfR[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
