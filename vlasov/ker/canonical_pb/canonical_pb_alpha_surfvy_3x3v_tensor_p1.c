#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_tensor_6x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_tensor_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvy_3x3v_tensor_p1(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[128];
  double *sgn_alpha_surfL = &sgn_alpha_surf[128];
  alphaL[0] = 2.121320343559642*hamil[14]*rdy2-1.224744871391589*hamil[2]*rdy2; 
  alphaL[1] = 2.121320343559642*hamil[26]*rdy2-1.224744871391589*hamil[7]*rdy2; 
  alphaL[3] = 2.121320343559642*hamil[28]*rdy2-1.224744871391589*hamil[9]*rdy2; 
  alphaL[4] = 2.121320343559642*hamil[30]*rdy2-1.224744871391589*hamil[11]*rdy2; 
  alphaL[5] = 2.121320343559642*hamil[39]*rdy2-1.224744871391589*hamil[18]*rdy2; 
  alphaL[7] = 2.121320343559642*hamil[43]*rdy2-1.224744871391589*hamil[22]*rdy2; 
  alphaL[9] = 2.121320343559642*hamil[44]*rdy2-1.224744871391589*hamil[23]*rdy2; 
  alphaL[11] = 2.121320343559642*hamil[46]*rdy2-1.224744871391589*hamil[25]*rdy2; 
  alphaL[12] = 2.121320343559642*hamil[51]*rdy2-1.224744871391589*hamil[32]*rdy2; 
  alphaL[14] = 2.121320343559642*hamil[53]*rdy2-1.224744871391589*hamil[34]*rdy2; 
  alphaL[15] = 2.121320343559642*hamil[55]*rdy2-1.224744871391589*hamil[36]*rdy2; 
  alphaL[18] = 2.121320343559642*hamil[57]*rdy2-1.224744871391589*hamil[42]*rdy2; 
  alphaL[21] = 2.121320343559642*hamil[59]*rdy2-1.224744871391589*hamil[47]*rdy2; 
  alphaL[23] = 2.121320343559642*hamil[60]*rdy2-1.224744871391589*hamil[48]*rdy2; 
  alphaL[25] = 2.121320343559642*hamil[62]*rdy2-1.224744871391589*hamil[50]*rdy2; 
  alphaL[29] = 2.121320343559642*hamil[63]*rdy2-1.224744871391589*hamil[58]*rdy2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*(alphaL[25]+alphaL[23]+alphaL[21]+alphaL[18])+0.1767766952966367*(alphaL[15]+alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9]+alphaL[7])-0.1767766952966367*(alphaL[5]+alphaL[4]+alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*(alphaL[25]+alphaL[23]+alphaL[21])-0.1767766952966367*(alphaL[18]+alphaL[15]+alphaL[14]+alphaL[12])+0.1767766952966367*(alphaL[11]+alphaL[9]+alphaL[7]+alphaL[5])-0.1767766952966367*(alphaL[4]+alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*(alphaL[25]+alphaL[23])-0.1767766952966367*alphaL[21]+0.1767766952966367*alphaL[18]-0.1767766952966367*alphaL[15]+0.1767766952966367*(alphaL[14]+alphaL[12])-0.1767766952966367*(alphaL[11]+alphaL[9])+0.1767766952966367*alphaL[7]-0.1767766952966367*alphaL[5]+0.1767766952966367*alphaL[4]-0.1767766952966367*(alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*(alphaL[25]+alphaL[23])+0.1767766952966367*(alphaL[21]+alphaL[18]+alphaL[15])-0.1767766952966367*(alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9])+0.1767766952966367*(alphaL[7]+alphaL[5]+alphaL[4])-0.1767766952966367*(alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*alphaL[25]-0.1767766952966367*alphaL[23]+0.1767766952966367*(alphaL[21]+alphaL[18]+alphaL[15])-0.1767766952966367*alphaL[14]+0.1767766952966367*alphaL[12]-0.1767766952966367*alphaL[11]+0.1767766952966367*alphaL[9]-0.1767766952966367*(alphaL[7]+alphaL[5]+alphaL[4])+0.1767766952966367*alphaL[3]-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*alphaL[25]+0.1767766952966367*alphaL[23]-0.1767766952966367*alphaL[21]+0.1767766952966367*alphaL[18]-0.1767766952966367*alphaL[15]+0.1767766952966367*alphaL[14]-0.1767766952966367*(alphaL[12]+alphaL[11])+0.1767766952966367*alphaL[9]-0.1767766952966367*alphaL[7]+0.1767766952966367*alphaL[5]-0.1767766952966367*alphaL[4]+0.1767766952966367*alphaL[3]-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*alphaL[25]+0.1767766952966367*(alphaL[23]+alphaL[21])-0.1767766952966367*(alphaL[18]+alphaL[15]+alphaL[14])+0.1767766952966367*(alphaL[12]+alphaL[11])-0.1767766952966367*(alphaL[9]+alphaL[7]+alphaL[5])+0.1767766952966367*(alphaL[4]+alphaL[3])-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*alphaL[25]-0.1767766952966367*(alphaL[23]+alphaL[21]+alphaL[18])+0.1767766952966367*(alphaL[15]+alphaL[14])-0.1767766952966367*alphaL[12]+0.1767766952966367*alphaL[11]-0.1767766952966367*(alphaL[9]+alphaL[7])+0.1767766952966367*(alphaL[5]+alphaL[4]+alphaL[3])-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*(alphaL[25]+alphaL[23]+alphaL[21]+alphaL[18])+0.1767766952966367*(alphaL[15]+alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9]+alphaL[7])-0.1767766952966367*(alphaL[5]+alphaL[4]+alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*(alphaL[25]+alphaL[23]+alphaL[21])-0.1767766952966367*(alphaL[18]+alphaL[15]+alphaL[14]+alphaL[12])+0.1767766952966367*(alphaL[11]+alphaL[9]+alphaL[7]+alphaL[5])-0.1767766952966367*(alphaL[4]+alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*(alphaL[25]+alphaL[23])-0.1767766952966367*alphaL[21]+0.1767766952966367*alphaL[18]-0.1767766952966367*alphaL[15]+0.1767766952966367*(alphaL[14]+alphaL[12])-0.1767766952966367*(alphaL[11]+alphaL[9])+0.1767766952966367*alphaL[7]-0.1767766952966367*alphaL[5]+0.1767766952966367*alphaL[4]-0.1767766952966367*(alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*(alphaL[25]+alphaL[23])+0.1767766952966367*(alphaL[21]+alphaL[18]+alphaL[15])-0.1767766952966367*(alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9])+0.1767766952966367*(alphaL[7]+alphaL[5]+alphaL[4])-0.1767766952966367*(alphaL[3]+alphaL[1])+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*alphaL[25]-0.1767766952966367*alphaL[23]+0.1767766952966367*(alphaL[21]+alphaL[18]+alphaL[15])-0.1767766952966367*alphaL[14]+0.1767766952966367*alphaL[12]-0.1767766952966367*alphaL[11]+0.1767766952966367*alphaL[9]-0.1767766952966367*(alphaL[7]+alphaL[5]+alphaL[4])+0.1767766952966367*alphaL[3]-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*alphaL[25]+0.1767766952966367*alphaL[23]-0.1767766952966367*alphaL[21]+0.1767766952966367*alphaL[18]-0.1767766952966367*alphaL[15]+0.1767766952966367*alphaL[14]-0.1767766952966367*(alphaL[12]+alphaL[11])+0.1767766952966367*alphaL[9]-0.1767766952966367*alphaL[7]+0.1767766952966367*alphaL[5]-0.1767766952966367*alphaL[4]+0.1767766952966367*alphaL[3]-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alphaL[29]-0.1767766952966367*alphaL[25]+0.1767766952966367*(alphaL[23]+alphaL[21])-0.1767766952966367*(alphaL[18]+alphaL[15]+alphaL[14])+0.1767766952966367*(alphaL[12]+alphaL[11])-0.1767766952966367*(alphaL[9]+alphaL[7]+alphaL[5])+0.1767766952966367*(alphaL[4]+alphaL[3])-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alphaL[29])+0.1767766952966367*alphaL[25]-0.1767766952966367*(alphaL[23]+alphaL[21]+alphaL[18])+0.1767766952966367*(alphaL[15]+alphaL[14])-0.1767766952966367*alphaL[12]+0.1767766952966367*alphaL[11]-0.1767766952966367*(alphaL[9]+alphaL[7])+0.1767766952966367*(alphaL[5]+alphaL[4]+alphaL[3])-0.1767766952966367*alphaL[1]+0.1767766952966367*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]))+0.1767766952966367*(alphaL[23]+alphaL[21]+alphaL[18]+alphaL[15]+alphaL[14])-0.1767766952966367*alphaL[12]+0.1767766952966367*alphaL[11]-0.1767766952966367*(alphaL[9]+alphaL[7]+alphaL[5]+alphaL[4]+alphaL[3])+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25])-0.1767766952966367*(alphaL[23]+alphaL[21])+0.1767766952966367*alphaL[18]-0.1767766952966367*(alphaL[15]+alphaL[14])+0.1767766952966367*(alphaL[12]+alphaL[11])-0.1767766952966367*(alphaL[9]+alphaL[7])+0.1767766952966367*alphaL[5]-0.1767766952966367*(alphaL[4]+alphaL[3])+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25])-0.1767766952966367*alphaL[23]+0.1767766952966367*alphaL[21]-0.1767766952966367*(alphaL[18]+alphaL[15])+0.1767766952966367*alphaL[14]-0.1767766952966367*(alphaL[12]+alphaL[11])+0.1767766952966367*alphaL[9]-0.1767766952966367*(alphaL[7]+alphaL[5])+0.1767766952966367*alphaL[4]-0.1767766952966367*alphaL[3]+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]))+0.1767766952966367*alphaL[23]-0.1767766952966367*(alphaL[21]+alphaL[18])+0.1767766952966367*alphaL[15]-0.1767766952966367*alphaL[14]+0.1767766952966367*alphaL[12]-0.1767766952966367*alphaL[11]+0.1767766952966367*alphaL[9]-0.1767766952966367*alphaL[7]+0.1767766952966367*(alphaL[5]+alphaL[4])-0.1767766952966367*alphaL[3]+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23])-0.1767766952966367*(alphaL[21]+alphaL[18])+0.1767766952966367*alphaL[15]-0.1767766952966367*(alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9])+0.1767766952966367*alphaL[7]-0.1767766952966367*(alphaL[5]+alphaL[4])+0.1767766952966367*(alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]))+0.1767766952966367*alphaL[21]-0.1767766952966367*(alphaL[18]+alphaL[15])+0.1767766952966367*(alphaL[14]+alphaL[12])-0.1767766952966367*(alphaL[11]+alphaL[9])+0.1767766952966367*(alphaL[7]+alphaL[5])-0.1767766952966367*alphaL[4]+0.1767766952966367*(alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]+alphaL[21]))+0.1767766952966367*alphaL[18]-0.1767766952966367*(alphaL[15]+alphaL[14]+alphaL[12])+0.1767766952966367*(alphaL[11]+alphaL[9]+alphaL[7])-0.1767766952966367*alphaL[5]+0.1767766952966367*(alphaL[4]+alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]+alphaL[21]+alphaL[18]+alphaL[15]+alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9]+alphaL[7]+alphaL[5]+alphaL[4]+alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]))+0.1767766952966367*(alphaL[23]+alphaL[21]+alphaL[18]+alphaL[15]+alphaL[14])-0.1767766952966367*alphaL[12]+0.1767766952966367*alphaL[11]-0.1767766952966367*(alphaL[9]+alphaL[7]+alphaL[5]+alphaL[4]+alphaL[3])+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25])-0.1767766952966367*(alphaL[23]+alphaL[21])+0.1767766952966367*alphaL[18]-0.1767766952966367*(alphaL[15]+alphaL[14])+0.1767766952966367*(alphaL[12]+alphaL[11])-0.1767766952966367*(alphaL[9]+alphaL[7])+0.1767766952966367*alphaL[5]-0.1767766952966367*(alphaL[4]+alphaL[3])+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25])-0.1767766952966367*alphaL[23]+0.1767766952966367*alphaL[21]-0.1767766952966367*(alphaL[18]+alphaL[15])+0.1767766952966367*alphaL[14]-0.1767766952966367*(alphaL[12]+alphaL[11])+0.1767766952966367*alphaL[9]-0.1767766952966367*(alphaL[7]+alphaL[5])+0.1767766952966367*alphaL[4]-0.1767766952966367*alphaL[3]+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]))+0.1767766952966367*alphaL[23]-0.1767766952966367*(alphaL[21]+alphaL[18])+0.1767766952966367*alphaL[15]-0.1767766952966367*alphaL[14]+0.1767766952966367*alphaL[12]-0.1767766952966367*alphaL[11]+0.1767766952966367*alphaL[9]-0.1767766952966367*alphaL[7]+0.1767766952966367*(alphaL[5]+alphaL[4])-0.1767766952966367*alphaL[3]+0.1767766952966367*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23])-0.1767766952966367*(alphaL[21]+alphaL[18])+0.1767766952966367*alphaL[15]-0.1767766952966367*(alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9])+0.1767766952966367*alphaL[7]-0.1767766952966367*(alphaL[5]+alphaL[4])+0.1767766952966367*(alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]))+0.1767766952966367*alphaL[21]-0.1767766952966367*(alphaL[18]+alphaL[15])+0.1767766952966367*(alphaL[14]+alphaL[12])-0.1767766952966367*(alphaL[11]+alphaL[9])+0.1767766952966367*(alphaL[7]+alphaL[5])-0.1767766952966367*alphaL[4]+0.1767766952966367*(alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]+alphaL[21]))+0.1767766952966367*alphaL[18]-0.1767766952966367*(alphaL[15]+alphaL[14]+alphaL[12])+0.1767766952966367*(alphaL[11]+alphaL[9]+alphaL[7])-0.1767766952966367*alphaL[5]+0.1767766952966367*(alphaL[4]+alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alphaL[29]+alphaL[25]+alphaL[23]+alphaL[21]+alphaL[18]+alphaL[15]+alphaL[14]+alphaL[12]+alphaL[11]+alphaL[9]+alphaL[7]+alphaL[5]+alphaL[4]+alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
