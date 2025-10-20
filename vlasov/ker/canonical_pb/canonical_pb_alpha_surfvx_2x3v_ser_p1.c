#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[80];
  double *sgn_alpha_surfL = &sgn_alpha_surf[108];
  alphaL[0] = (-2.738612787525831*hamil[33]*rdx2)+2.121320343559642*hamil[7]*rdx2-1.224744871391589*hamil[1]*rdx2; 
  alphaL[2] = (-2.738612787525831*hamil[37]*rdx2)+2.121320343559642*hamil[16]*rdx2-1.224744871391589*hamil[6]*rdx2; 
  alphaL[3] = (-2.738612787525831*hamil[38]*rdx2)+2.121320343559642*hamil[18]*rdx2-1.224744871391589*hamil[9]*rdx2; 
  alphaL[4] = (-2.738612787525831*hamil[40]*rdx2)+2.121320343559642*hamil[21]*rdx2-1.224744871391589*hamil[12]*rdx2; 
  alphaL[7] = (-2.738612787525831*hamil[43]*rdx2)+2.121320343559642*hamil[26]*rdx2-1.224744871391589*hamil[17]*rdx2; 
  alphaL[9] = (-2.738612787525831*hamil[44]*rdx2)+2.121320343559642*hamil[27]*rdx2-1.224744871391589*hamil[20]*rdx2; 
  alphaL[10] = (-2.738612787525831*hamil[45]*rdx2)+2.121320343559642*hamil[29]*rdx2-1.224744871391589*hamil[23]*rdx2; 
  alphaL[14] = (-2.738612787525831*hamil[47]*rdx2)+2.121320343559642*hamil[31]*rdx2-1.224744871391589*hamil[28]*rdx2; 
  alphaL[16] = 2.121320343559642*hamil[54]*rdx2-1.224744871391589*hamil[49]*rdx2; 
  alphaL[18] = 2.121320343559642*hamil[59]*rdx2-1.224744871391589*hamil[53]*rdx2; 
  alphaL[19] = 2.121320343559642*hamil[61]*rdx2-1.224744871391589*hamil[56]*rdx2; 
  alphaL[22] = 2.121320343559642*hamil[63]*rdx2-1.224744871391589*hamil[60]*rdx2; 
  alphaL[24] = 2.121320343559642*hamil[70]*rdx2-1.224744871391589*hamil[65]*rdx2; 
  alphaL[26] = 2.121320343559642*hamil[75]*rdx2-1.224744871391589*hamil[69]*rdx2; 
  alphaL[27] = 2.121320343559642*hamil[77]*rdx2-1.224744871391589*hamil[72]*rdx2; 
  alphaL[30] = 2.121320343559642*hamil[79]*rdx2-1.224744871391589*hamil[76]*rdx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3*alphaL[30]-0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]+0.3*alphaL[22]-0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.45*alphaL[14]+0.45*alphaL[10]+0.3354101966249678*(alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.375*alphaL[30])+0.375*alphaL[27]+0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[30]-0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]-0.3*alphaL[22]+0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.45*alphaL[14]-0.45*alphaL[10]-0.3354101966249678*alphaL[9]+0.3354101966249678*(alphaL[7]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[26])+0.2236067977499786*alphaL[24]-0.375*alphaL[22]+0.375*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.3354101966249678*alphaL[9]-0.3354101966249678*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[26])+0.2236067977499786*alphaL[24]+0.375*alphaL[22]-0.375*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.3354101966249678*alphaL[9]+0.3354101966249678*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[30])+0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]+0.3*alphaL[22]-0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.45*alphaL[14]-0.45*alphaL[10]+0.3354101966249678*alphaL[9]-0.3354101966249678*(alphaL[7]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[30]-0.375*alphaL[27]+0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[30])+0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]-0.3*alphaL[22]+0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.45*alphaL[14]+0.45*alphaL[10]-0.3354101966249678*(alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[30]+alphaL[27]))+0.2236067977499786*(alphaL[26]+alphaL[24])-0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.45*(alphaL[14]+alphaL[10])-0.3354101966249678*(alphaL[9]+alphaL[7]+alphaL[4]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[27])-0.2795084971874732*(alphaL[26]+alphaL[24])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.3354101966249678*(alphaL[7]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[30]+alphaL[27]))+0.2236067977499786*(alphaL[26]+alphaL[24])+0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.45*(alphaL[14]+alphaL[10])+0.3354101966249678*alphaL[9]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[26]+alphaL[24])+0.375*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])-0.3354101966249678*(alphaL[9]+alphaL[4])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[2]+alphaL[0])-0.2795084971874732*(alphaL[26]+alphaL[24]+alphaL[18]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[26]+alphaL[24])-0.375*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])+0.3354101966249678*(alphaL[9]+alphaL[4])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[30]+alphaL[27])+0.2236067977499786*(alphaL[26]+alphaL[24])-0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.45*(alphaL[14]+alphaL[10])-0.3354101966249678*alphaL[9]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[27]))-0.2795084971874732*(alphaL[26]+alphaL[24])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.3354101966249678*(alphaL[7]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[30]+alphaL[27])+0.2236067977499786*(alphaL[26]+alphaL[24])+0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.45*(alphaL[14]+alphaL[10])+0.3354101966249678*(alphaL[9]+alphaL[7]+alphaL[4]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[30]-0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]+0.3*alphaL[22]-0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.45*alphaL[14]+0.45*alphaL[10]+0.3354101966249678*(alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL[30])+0.375*alphaL[27]+0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[30]-0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]-0.3*alphaL[22]+0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.45*alphaL[14]-0.45*alphaL[10]-0.3354101966249678*alphaL[9]+0.3354101966249678*(alphaL[7]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[26])+0.2236067977499786*alphaL[24]-0.375*alphaL[22]+0.375*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.3354101966249678*alphaL[9]-0.3354101966249678*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[26])+0.2236067977499786*alphaL[24]+0.375*alphaL[22]-0.375*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.3354101966249678*alphaL[9]+0.3354101966249678*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[30])+0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]+0.3*alphaL[22]-0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.45*alphaL[14]-0.45*alphaL[10]+0.3354101966249678*alphaL[9]-0.3354101966249678*(alphaL[7]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[30]-0.375*alphaL[27]+0.2795084971874732*alphaL[26]-0.2795084971874732*alphaL[24]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[30])+0.3*alphaL[27]-0.2236067977499786*alphaL[26]+0.2236067977499786*alphaL[24]-0.3*alphaL[22]+0.3*alphaL[19]-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.45*alphaL[14]+0.45*alphaL[10]-0.3354101966249678*(alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[30]+alphaL[27]))+0.2236067977499786*(alphaL[26]+alphaL[24])-0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.45*(alphaL[14]+alphaL[10])-0.3354101966249678*(alphaL[9]+alphaL[7]+alphaL[4]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[27])-0.2795084971874732*(alphaL[26]+alphaL[24])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.3354101966249678*(alphaL[7]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[30]+alphaL[27]))+0.2236067977499786*(alphaL[26]+alphaL[24])+0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.45*(alphaL[14]+alphaL[10])+0.3354101966249678*alphaL[9]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[26]+alphaL[24])+0.375*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])-0.3354101966249678*(alphaL[9]+alphaL[4])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[2]+alphaL[0])-0.2795084971874732*(alphaL[26]+alphaL[24]+alphaL[18]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[26]+alphaL[24])-0.375*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])+0.3354101966249678*(alphaL[9]+alphaL[4])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[32] = 1.0; 
  else  
    sgn_alpha_surfL[32] = -1.0; 
  
  if (sgn_alpha_surfL[32] == sgn_alpha_surfL[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[30]+alphaL[27])+0.2236067977499786*(alphaL[26]+alphaL[24])-0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])-0.45*(alphaL[14]+alphaL[10])-0.3354101966249678*alphaL[9]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[33] = 1.0; 
  else  
    sgn_alpha_surfL[33] = -1.0; 
  
  if (sgn_alpha_surfL[33] == sgn_alpha_surfL[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[27]))-0.2795084971874732*(alphaL[26]+alphaL[24])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.3354101966249678*(alphaL[7]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[34] = 1.0; 
  else  
    sgn_alpha_surfL[34] = -1.0; 
  
  if (sgn_alpha_surfL[34] == sgn_alpha_surfL[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[30]+alphaL[27])+0.2236067977499786*(alphaL[26]+alphaL[24])+0.3*(alphaL[22]+alphaL[19])+0.2236067977499786*(alphaL[18]+alphaL[16])+0.45*(alphaL[14]+alphaL[10])+0.3354101966249678*(alphaL[9]+alphaL[7]+alphaL[4]+alphaL[3])+0.25*(alphaL[2]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[35] = 1.0; 
  else  
    sgn_alpha_surfL[35] = -1.0; 
  
  if (sgn_alpha_surfL[35] == sgn_alpha_surfL[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
