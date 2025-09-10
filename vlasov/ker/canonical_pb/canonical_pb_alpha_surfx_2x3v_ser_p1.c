#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[0];
  double *sgn_alpha_surfL = &sgn_alpha_surf[0];
  alphaL[0] = 1.224744871391589*hamil[3]*rdvx2-2.121320343559642*hamil[7]*rdvx2; 
  alphaL[1] = 1.224744871391589*hamil[8]*rdvx2-2.121320343559642*hamil[16]*rdvx2; 
  alphaL[2] = 2.738612787525831*hamil[32]*rdvx2-4.743416490252569*hamil[33]*rdvx2; 
  alphaL[3] = 1.224744871391589*hamil[11]*rdvx2-2.121320343559642*hamil[18]*rdvx2; 
  alphaL[4] = 1.224744871391589*hamil[14]*rdvx2-2.121320343559642*hamil[21]*rdvx2; 
  alphaL[5] = 2.738612787525831*hamil[34]*rdvx2-4.743416490252569*hamil[37]*rdvx2; 
  alphaL[6] = 1.224744871391589*hamil[19]*rdvx2-2.121320343559642*hamil[26]*rdvx2; 
  alphaL[7] = 2.738612787525831*hamil[35]*rdvx2-4.743416490252569*hamil[38]*rdvx2; 
  alphaL[8] = 1.224744871391589*hamil[22]*rdvx2-2.121320343559642*hamil[27]*rdvx2; 
  alphaL[9] = 2.738612787525831*hamil[36]*rdvx2-4.743416490252569*hamil[40]*rdvx2; 
  alphaL[10] = 1.224744871391589*hamil[25]*rdvx2-2.121320343559642*hamil[29]*rdvx2; 
  alphaL[11] = 2.738612787525831*hamil[39]*rdvx2-4.743416490252569*hamil[43]*rdvx2; 
  alphaL[12] = 2.738612787525831*hamil[41]*rdvx2-4.743416490252569*hamil[44]*rdvx2; 
  alphaL[13] = 1.224744871391589*hamil[30]*rdvx2-2.121320343559642*hamil[31]*rdvx2; 
  alphaL[14] = 2.738612787525831*hamil[42]*rdvx2-4.743416490252569*hamil[45]*rdvx2; 
  alphaL[15] = 2.738612787525831*hamil[46]*rdvx2-4.743416490252569*hamil[47]*rdvx2; 
  alphaL[24] = 1.224744871391589*hamil[51]*rdvx2-2.121320343559642*hamil[54]*rdvx2; 
  alphaL[25] = 1.224744871391589*hamil[55]*rdvx2-2.121320343559642*hamil[59]*rdvx2; 
  alphaL[27] = 1.224744871391589*hamil[58]*rdvx2-2.121320343559642*hamil[61]*rdvx2; 
  alphaL[29] = 1.224744871391589*hamil[62]*rdvx2-2.121320343559642*hamil[63]*rdvx2; 
  alphaL[32] = 1.224744871391589*hamil[67]*rdvx2-2.121320343559642*hamil[70]*rdvx2; 
  alphaL[33] = 1.224744871391589*hamil[71]*rdvx2-2.121320343559642*hamil[75]*rdvx2; 
  alphaL[35] = 1.224744871391589*hamil[74]*rdvx2-2.121320343559642*hamil[77]*rdvx2; 
  alphaL[37] = 1.224744871391589*hamil[78]*rdvx2-2.121320343559642*hamil[79]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.603738353924943*alphaL[15]-0.603738353924943*alphaL[14]-0.45*(alphaL[13]+alphaL[12]+alphaL[11])+0.45*(alphaL[10]+alphaL[9])+0.3354101966249678*alphaL[8]+0.45*alphaL[7]+0.3354101966249678*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.375*alphaL[37])+0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.45*alphaL[11]+0.45*alphaL[7]+0.3354101966249678*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.603738353924943*alphaL[15]+0.603738353924943*alphaL[14]+0.45*(alphaL[13]+alphaL[12])-0.45*(alphaL[11]+alphaL[10]+alphaL[9])-0.3354101966249678*alphaL[8]+0.45*alphaL[7]+0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[4])-0.3354101966249678*(alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]-0.375*alphaL[29]+0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.45*alphaL[12]+0.45*alphaL[9]+0.3354101966249678*(alphaL[8]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.3354101966249678*alphaL[5]-0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]+0.375*alphaL[29]-0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.45*alphaL[12]-0.45*alphaL[9]-0.3354101966249678*alphaL[8]+0.3354101966249678*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.603738353924943*alphaL[15]+0.603738353924943*alphaL[14]+0.45*alphaL[13]-0.45*alphaL[12]+0.45*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[9]+0.3354101966249678*alphaL[8]-0.45*alphaL[7]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[37]-0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.45*alphaL[11]-0.45*alphaL[7]-0.3354101966249678*alphaL[6]+0.3354101966249678*(alphaL[5]+alphaL[3])-0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.603738353924943*alphaL[15]-0.603738353924943*alphaL[14]-0.45*alphaL[13]+0.45*(alphaL[12]+alphaL[11]+alphaL[10])-0.45*alphaL[9]-0.3354101966249678*alphaL[8]-0.45*alphaL[7]-0.3354101966249678*alphaL[6]+0.3354101966249678*(alphaL[5]+alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.45*alphaL[13]+0.45*alphaL[10]+0.3354101966249678*(alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL[37])+0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.45*alphaL[13]-0.45*alphaL[10]-0.3354101966249678*alphaL[8]+0.3354101966249678*(alphaL[6]+alphaL[4])-0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]-0.375*alphaL[29]+0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]+0.375*alphaL[29]-0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[4]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.45*alphaL[13]-0.45*alphaL[10]+0.3354101966249678*alphaL[8]-0.3354101966249678*(alphaL[6]+alphaL[4])+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[37]-0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[3]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.45*alphaL[13]+0.45*alphaL[10]-0.3354101966249678*(alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.603738353924943*alphaL[15]+0.603738353924943*alphaL[14]-0.45*alphaL[13]+0.45*(alphaL[12]+alphaL[11]+alphaL[10])-0.45*alphaL[9]+0.3354101966249678*alphaL[8]-0.45*alphaL[7]+0.3354101966249678*alphaL[6]-0.3354101966249678*(alphaL[5]+alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL[37])+0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.45*alphaL[11]-0.45*alphaL[7]+0.3354101966249678*alphaL[6]-0.3354101966249678*(alphaL[5]+alphaL[3])+0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaL[37]-0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.603738353924943*alphaL[15]-0.603738353924943*alphaL[14]+0.45*alphaL[13]-0.45*alphaL[12]+0.45*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[9]-0.3354101966249678*alphaL[8]-0.45*alphaL[7]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]-0.375*alphaL[29]+0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]+0.45*alphaL[12]-0.45*alphaL[9]+0.3354101966249678*alphaL[8]-0.3354101966249678*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.3354101966249678*alphaL[5]+0.3354101966249678*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL[33])+0.2236067977499786*alphaL[32]+0.375*alphaL[29]-0.375*alphaL[27]+0.2795084971874732*alphaL[25]-0.2795084971874732*alphaL[24]-0.45*alphaL[12]+0.45*alphaL[9]-0.3354101966249678*(alphaL[8]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]+0.3*alphaL[29]-0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]+0.603738353924943*alphaL[15]-0.603738353924943*alphaL[14]+0.45*(alphaL[13]+alphaL[12])-0.45*(alphaL[11]+alphaL[10]+alphaL[9])+0.3354101966249678*alphaL[8]+0.45*alphaL[7]-0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[4])+0.3354101966249678*(alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL[37]-0.375*alphaL[35]+0.2795084971874732*alphaL[33]-0.2795084971874732*alphaL[32]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.45*alphaL[11]+0.45*alphaL[7]-0.3354101966249678*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaL[37])+0.3*alphaL[35]-0.2236067977499786*alphaL[33]+0.2236067977499786*alphaL[32]-0.3*alphaL[29]+0.3*alphaL[27]-0.2236067977499786*alphaL[25]+0.2236067977499786*alphaL[24]-0.603738353924943*alphaL[15]+0.603738353924943*alphaL[14]-0.45*(alphaL[13]+alphaL[12]+alphaL[11])+0.45*(alphaL[10]+alphaL[9])-0.3354101966249678*alphaL[8]+0.45*alphaL[7]-0.3354101966249678*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.603738353924943*(alphaL[15]+alphaL[14])+0.45*(alphaL[13]+alphaL[12]+alphaL[11]+alphaL[10]+alphaL[9])-0.3354101966249678*alphaL[8]+0.45*alphaL[7]-0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[4]+alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[37]+alphaL[35])-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.45*(alphaL[11]+alphaL[7])-0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.603738353924943*(alphaL[15]+alphaL[14])-0.45*(alphaL[13]+alphaL[12])+0.45*alphaL[11]-0.45*(alphaL[10]+alphaL[9])+0.3354101966249678*alphaL[8]+0.45*alphaL[7]-0.3354101966249678*(alphaL[6]+alphaL[5])+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])+0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])+0.45*(alphaL[12]+alphaL[9])-0.3354101966249678*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[33]+alphaL[32]+alphaL[25]+alphaL[24]))-0.3354101966249678*(alphaL[5]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])-0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])-0.45*(alphaL[12]+alphaL[9])+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[32] = 1.0; 
  else  
    sgn_alpha_surfL[32] = -1.0; 
  
  if (sgn_alpha_surfL[32] == sgn_alpha_surfL[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.603738353924943*(alphaL[15]+alphaL[14])-0.45*alphaL[13]+0.45*alphaL[12]-0.45*(alphaL[11]+alphaL[10])+0.45*alphaL[9]-0.3354101966249678*alphaL[8]-0.45*alphaL[7]+0.3354101966249678*alphaL[6]-0.3354101966249678*(alphaL[5]+alphaL[4])+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[33] = 1.0; 
  else  
    sgn_alpha_surfL[33] = -1.0; 
  
  if (sgn_alpha_surfL[33] == sgn_alpha_surfL[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[37]+alphaL[35]))-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.45*(alphaL[11]+alphaL[7])+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[34] = 1.0; 
  else  
    sgn_alpha_surfL[34] = -1.0; 
  
  if (sgn_alpha_surfL[34] == sgn_alpha_surfL[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.603738353924943*(alphaL[15]+alphaL[14])+0.45*alphaL[13]-0.45*(alphaL[12]+alphaL[11])+0.45*alphaL[10]-0.45*alphaL[9]+0.3354101966249678*alphaL[8]-0.45*alphaL[7]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[5]+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[35] = 1.0; 
  else  
    sgn_alpha_surfL[35] = -1.0; 
  
  if (sgn_alpha_surfL[35] == sgn_alpha_surfL[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[36] = 1.0; 
  else  
    sgn_alpha_surfL[36] = -1.0; 
  
  if (sgn_alpha_surfL[36] == sgn_alpha_surfL[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[37]+alphaL[35])-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[37] = 1.0; 
  else  
    sgn_alpha_surfL[37] = -1.0; 
  
  if (sgn_alpha_surfL[37] == sgn_alpha_surfL[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*alphaL[8]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[38] = 1.0; 
  else  
    sgn_alpha_surfL[38] = -1.0; 
  
  if (sgn_alpha_surfL[38] == sgn_alpha_surfL[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])+0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])-0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[39] = 1.0; 
  else  
    sgn_alpha_surfL[39] = -1.0; 
  
  if (sgn_alpha_surfL[39] == sgn_alpha_surfL[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[33]+alphaL[32]+alphaL[25]+alphaL[24]) > 0.) 
    sgn_alpha_surfL[40] = 1.0; 
  else  
    sgn_alpha_surfL[40] = -1.0; 
  
  if (sgn_alpha_surfL[40] == sgn_alpha_surfL[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])-0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])+0.3354101966249678*(alphaL[8]+alphaL[4])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[41] = 1.0; 
  else  
    sgn_alpha_surfL[41] = -1.0; 
  
  if (sgn_alpha_surfL[41] == sgn_alpha_surfL[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.45*(alphaL[13]+alphaL[10])-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[42] = 1.0; 
  else  
    sgn_alpha_surfL[42] = -1.0; 
  
  if (sgn_alpha_surfL[42] == sgn_alpha_surfL[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[37]+alphaL[35]))-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.3354101966249678*(alphaL[6]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[43] = 1.0; 
  else  
    sgn_alpha_surfL[43] = -1.0; 
  
  if (sgn_alpha_surfL[43] == sgn_alpha_surfL[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.45*(alphaL[13]+alphaL[10])+0.3354101966249678*(alphaL[8]+alphaL[6]+alphaL[4]+alphaL[3])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[44] = 1.0; 
  else  
    sgn_alpha_surfL[44] = -1.0; 
  
  if (sgn_alpha_surfL[44] == sgn_alpha_surfL[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.603738353924943*(alphaL[15]+alphaL[14])+0.45*alphaL[13]-0.45*(alphaL[12]+alphaL[11])+0.45*alphaL[10]-0.45*alphaL[9]-0.3354101966249678*alphaL[8]-0.45*alphaL[7]-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[5]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[45] = 1.0; 
  else  
    sgn_alpha_surfL[45] = -1.0; 
  
  if (sgn_alpha_surfL[45] == sgn_alpha_surfL[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[37]+alphaL[35])-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.45*(alphaL[11]+alphaL[7])-0.3354101966249678*alphaL[6]+0.3354101966249678*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[46] = 1.0; 
  else  
    sgn_alpha_surfL[46] = -1.0; 
  
  if (sgn_alpha_surfL[46] == sgn_alpha_surfL[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaL[37]+alphaL[35]))+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.603738353924943*(alphaL[15]+alphaL[14])-0.45*alphaL[13]+0.45*alphaL[12]-0.45*(alphaL[11]+alphaL[10])+0.45*alphaL[9]+0.3354101966249678*alphaL[8]-0.45*alphaL[7]-0.3354101966249678*alphaL[6]+0.3354101966249678*(alphaL[5]+alphaL[4])-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[47] = 1.0; 
  else  
    sgn_alpha_surfL[47] = -1.0; 
  
  if (sgn_alpha_surfL[47] == sgn_alpha_surfL[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])+0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])-0.45*(alphaL[12]+alphaL[9])-0.3354101966249678*alphaL[8]+0.3354101966249678*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[48] = 1.0; 
  else  
    sgn_alpha_surfL[48] = -1.0; 
  
  if (sgn_alpha_surfL[48] == sgn_alpha_surfL[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[33]+alphaL[32]+alphaL[25]+alphaL[24]))+0.3354101966249678*(alphaL[5]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[49] = 1.0; 
  else  
    sgn_alpha_surfL[49] = -1.0; 
  
  if (sgn_alpha_surfL[49] == sgn_alpha_surfL[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL[33]+alphaL[32])-0.375*(alphaL[29]+alphaL[27])-0.2795084971874732*(alphaL[25]+alphaL[24])+0.45*(alphaL[12]+alphaL[9])+0.3354101966249678*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[50] = 1.0; 
  else  
    sgn_alpha_surfL[50] = -1.0; 
  
  if (sgn_alpha_surfL[50] == sgn_alpha_surfL[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])-0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])-0.603738353924943*(alphaL[15]+alphaL[14])-0.45*(alphaL[13]+alphaL[12])+0.45*alphaL[11]-0.45*(alphaL[10]+alphaL[9])-0.3354101966249678*alphaL[8]+0.45*alphaL[7]+0.3354101966249678*(alphaL[6]+alphaL[5])-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[51] = 1.0; 
  else  
    sgn_alpha_surfL[51] = -1.0; 
  
  if (sgn_alpha_surfL[51] == sgn_alpha_surfL[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[37]+alphaL[35]))-0.2795084971874732*(alphaL[33]+alphaL[32])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.45*(alphaL[11]+alphaL[7])+0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[52] = 1.0; 
  else  
    sgn_alpha_surfL[52] = -1.0; 
  
  if (sgn_alpha_surfL[52] == sgn_alpha_surfL[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaL[37]+alphaL[35])+0.2236067977499786*(alphaL[33]+alphaL[32])+0.3*(alphaL[29]+alphaL[27])+0.2236067977499786*(alphaL[25]+alphaL[24])+0.603738353924943*(alphaL[15]+alphaL[14])+0.45*(alphaL[13]+alphaL[12]+alphaL[11]+alphaL[10]+alphaL[9])+0.3354101966249678*alphaL[8]+0.45*alphaL[7]+0.3354101966249678*(alphaL[6]+alphaL[5]+alphaL[4]+alphaL[3]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[53] = 1.0; 
  else  
    sgn_alpha_surfL[53] = -1.0; 
  
  if (sgn_alpha_surfL[53] == sgn_alpha_surfL[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
