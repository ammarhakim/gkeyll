#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[112];
  double *sgn_alpha_surfL = &sgn_alpha_surf[144];
  alphaL[0] = (-2.738612787525831*hamil[50]*rdy2)+2.121320343559642*hamil[10]*rdy2-1.224744871391589*hamil[2]*rdy2; 
  alphaL[1] = (-2.738612787525831*hamil[53]*rdy2)+2.121320343559642*hamil[17]*rdy2-1.224744871391589*hamil[6]*rdy2; 
  alphaL[3] = (-2.738612787525831*hamil[55]*rdy2)+2.121320343559642*hamil[19]*rdy2-1.224744871391589*hamil[8]*rdy2; 
  alphaL[4] = (-2.738612787525831*hamil[57]*rdy2)+2.121320343559642*hamil[24]*rdy2-1.224744871391589*hamil[13]*rdy2; 
  alphaL[6] = (-2.738612787525831*hamil[59]*rdy2)+2.121320343559642*hamil[26]*rdy2-1.224744871391589*hamil[16]*rdy2; 
  alphaL[8] = (-2.738612787525831*hamil[60]*rdy2)+2.121320343559642*hamil[28]*rdy2-1.224744871391589*hamil[20]*rdy2; 
  alphaL[10] = (-2.738612787525831*hamil[62]*rdy2)+2.121320343559642*hamil[30]*rdy2-1.224744871391589*hamil[22]*rdy2; 
  alphaL[13] = (-2.738612787525831*hamil[63]*rdy2)+2.121320343559642*hamil[31]*rdy2-1.224744871391589*hamil[27]*rdy2; 
  alphaL[16] = 2.121320343559642*hamil[39]*rdy2-1.224744871391589*hamil[34]*rdy2; 
  alphaL[17] = 2.121320343559642*hamil[43]*rdy2-1.224744871391589*hamil[37]*rdy2; 
  alphaL[19] = 2.121320343559642*hamil[46]*rdy2-1.224744871391589*hamil[41]*rdy2; 
  alphaL[21] = 2.121320343559642*hamil[47]*rdy2-1.224744871391589*hamil[44]*rdy2; 
  alphaL[24] = 2.121320343559642*hamil[73]*rdy2-1.224744871391589*hamil[66]*rdy2; 
  alphaL[25] = 2.121320343559642*hamil[76]*rdy2-1.224744871391589*hamil[69]*rdy2; 
  alphaL[27] = 2.121320343559642*hamil[78]*rdy2-1.224744871391589*hamil[71]*rdy2; 
  alphaL[29] = 2.121320343559642*hamil[79]*rdy2-1.224744871391589*hamil[75]*rdy2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.3*alphaL[21]-0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.45*alphaL[13]+0.45*alphaL[10]+0.3354101966249678*(alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.375*alphaL[29])+0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.3*alphaL[21]+0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.45*alphaL[13]-0.45*alphaL[10]-0.3354101966249678*alphaL[8]+0.3354101966249678*(alphaL[6]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[25])+0.2236067977499786*alphaL[24]-0.375*alphaL[21]+0.375*alphaL[19]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[25])+0.2236067977499786*alphaL[24]+0.375*alphaL[21]-0.375*alphaL[19]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[29])+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.3*alphaL[21]-0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.45*alphaL[13]-0.45*alphaL[10]+0.3354101966249678*alphaL[8]-0.3354101966249678*(alphaL[6]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[29]-0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[29])+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.3*alphaL[21]+0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.45*alphaL[13]+0.45*alphaL[10]-0.3354101966249678*(alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.3*alphaL[21]-0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.45*alphaL[13]+0.45*alphaL[10]+0.3354101966249678*(alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL[29])+0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.3*alphaL[21]+0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.45*alphaL[13]-0.45*alphaL[10]-0.3354101966249678*alphaL[8]+0.3354101966249678*(alphaL[6]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[25])+0.2236067977499786*alphaL[24]-0.375*alphaL[21]+0.375*alphaL[19]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[25])+0.2236067977499786*alphaL[24]+0.375*alphaL[21]-0.375*alphaL[19]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[29])+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.3*alphaL[21]-0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.45*alphaL[13]-0.45*alphaL[10]+0.3354101966249678*alphaL[8]-0.3354101966249678*(alphaL[6]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[29]-0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[29])+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.3*alphaL[21]+0.3*alphaL[19]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.45*alphaL[13]+0.45*alphaL[10]-0.3354101966249678*(alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[29]+alphaL[27]))+0.2236067977499786*(alphaL[25]+alphaL[24])-0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[29]+alphaL[27]))+0.2236067977499786*(alphaL[25]+alphaL[24])+0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[25]+alphaL[24])+0.375*(alphaL[21]+alphaL[19])-0.2795084971874732*(alphaL[17]+alphaL[16])-0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[25]+alphaL[24]+alphaL[17]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[25]+alphaL[24])-0.375*(alphaL[21]+alphaL[19])-0.2795084971874732*(alphaL[17]+alphaL[16])+0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[29]+alphaL[27]))-0.2795084971874732*(alphaL[25]+alphaL[24])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[29]+alphaL[27]))+0.2236067977499786*(alphaL[25]+alphaL[24])-0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[29]+alphaL[27]))+0.2236067977499786*(alphaL[25]+alphaL[24])+0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[25]+alphaL[24])+0.375*(alphaL[21]+alphaL[19])-0.2795084971874732*(alphaL[17]+alphaL[16])-0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[25]+alphaL[24]+alphaL[17]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[25]+alphaL[24])-0.375*(alphaL[21]+alphaL[19])-0.2795084971874732*(alphaL[17]+alphaL[16])+0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[32] = 1.0; 
  else  
    sgn_alpha_surfL[32] = -1.0; 
  
  if (sgn_alpha_surfL[32] == sgn_alpha_surfL[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])-0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[33] = 1.0; 
  else  
    sgn_alpha_surfL[33] = -1.0; 
  
  if (sgn_alpha_surfL[33] == sgn_alpha_surfL[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[29]+alphaL[27]))-0.2795084971874732*(alphaL[25]+alphaL[24])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[34] = 1.0; 
  else  
    sgn_alpha_surfL[34] = -1.0; 
  
  if (sgn_alpha_surfL[34] == sgn_alpha_surfL[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.3*(alphaL[21]+alphaL[19])+0.2236067977499786*(alphaL[17]+alphaL[16])+0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[35] = 1.0; 
  else  
    sgn_alpha_surfL[35] = -1.0; 
  
  if (sgn_alpha_surfL[35] == sgn_alpha_surfL[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
