#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
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

  double *alphaR = &alpha_surf[48];
  double *sgn_alpha_surfR = &sgn_alpha_surf[81];
  alphaR[0] = 2.738612787525831*hamil[38]*rdvy2+2.121320343559642*hamil[10]*rdvy2+1.224744871391589*hamil[4]*rdvy2; 
  alphaR[1] = 2.738612787525831*hamil[60]*rdvy2+2.121320343559642*hamil[22]*rdvy2+1.224744871391589*hamil[9]*rdvy2; 
  alphaR[2] = 2.738612787525831*hamil[62]*rdvy2+2.121320343559642*hamil[24]*rdvy2+1.224744871391589*hamil[11]*rdvy2; 
  alphaR[3] = 4.743416490252569*hamil[41]*rdvy2+2.738612787525831*hamil[19]*rdvy2; 
  alphaR[4] = 2.738612787525831*hamil[75]*rdvy2+2.121320343559642*hamil[29]*rdvy2+1.224744871391589*hamil[15]*rdvy2; 
  alphaR[5] = 2.738612787525831*hamil[88]*rdvy2+2.121320343559642*hamil[51]*rdvy2+1.224744871391589*hamil[23]*rdvy2; 
  alphaR[6] = 4.743416490252569*hamil[65]*rdvy2+2.738612787525831*hamil[40]*rdvy2; 
  alphaR[7] = 4.743416490252569*hamil[67]*rdvy2+2.738612787525831*hamil[42]*rdvy2; 
  alphaR[8] = 2.738612787525831*hamil[95]*rdvy2+2.121320343559642*hamil[53]*rdvy2+1.224744871391589*hamil[28]*rdvy2; 
  alphaR[9] = 2.738612787525831*hamil[97]*rdvy2+2.121320343559642*hamil[55]*rdvy2+1.224744871391589*hamil[30]*rdvy2; 
  alphaR[10] = 4.743416490252569*hamil[78]*rdvy2+2.738612787525831*hamil[46]*rdvy2; 
  alphaR[11] = 2.121320343559642*hamil[59]*rdvy2+1.224744871391589*hamil[37]*rdvy2; 
  alphaR[12] = 2.121320343559642*hamil[64]*rdvy2+1.224744871391589*hamil[39]*rdvy2; 
  alphaR[14] = 2.121320343559642*hamil[84]*rdvy2+1.224744871391589*hamil[50]*rdvy2; 
  alphaR[15] = 4.743416490252569*hamil[90]*rdvy2+2.738612787525831*hamil[66]*rdvy2; 
  alphaR[16] = 2.738612787525831*hamil[108]*rdvy2+2.121320343559642*hamil[86]*rdvy2+1.224744871391589*hamil[54]*rdvy2; 
  alphaR[17] = 4.743416490252569*hamil[100]*rdvy2+2.738612787525831*hamil[77]*rdvy2; 
  alphaR[18] = 4.743416490252569*hamil[102]*rdvy2+2.738612787525831*hamil[79]*rdvy2; 
  alphaR[19] = 2.121320343559642*hamil[87]*rdvy2+1.224744871391589*hamil[61]*rdvy2; 
  alphaR[20] = 2.121320343559642*hamil[89]*rdvy2+1.224744871391589*hamil[63]*rdvy2; 
  alphaR[25] = 2.121320343559642*hamil[94]*rdvy2+1.224744871391589*hamil[74]*rdvy2; 
  alphaR[26] = 2.121320343559642*hamil[99]*rdvy2+1.224744871391589*hamil[76]*rdvy2; 
  alphaR[28] = 2.121320343559642*hamil[104]*rdvy2+1.224744871391589*hamil[83]*rdvy2; 
  alphaR[29] = 2.121320343559642*hamil[106]*rdvy2+1.224744871391589*hamil[85]*rdvy2; 
  alphaR[31] = 4.743416490252569*hamil[110]*rdvy2+2.738612787525831*hamil[101]*rdvy2; 
  alphaR[35] = 2.121320343559642*hamil[107]*rdvy2+1.224744871391589*hamil[96]*rdvy2; 
  alphaR[36] = 2.121320343559642*hamil[109]*rdvy2+1.224744871391589*hamil[98]*rdvy2; 
  alphaR[41] = 2.121320343559642*hamil[111]*rdvy2+1.224744871391589*hamil[105]*rdvy2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])+0.81*alphaR[31]-0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])-0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.5031152949374518*alphaR[41])+0.375*(alphaR[29]+alphaR[28])-0.3*(alphaR[20]+alphaR[19])-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])-0.81*alphaR[31]-0.3*(alphaR[29]+alphaR[28])+0.3*(alphaR[26]+alphaR[25])-0.3*(alphaR[20]+alphaR[19])+0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16])-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[10]+alphaR[9]+alphaR[8])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])-0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])-0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[9]+alphaR[8]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[41])+0.375*(alphaR[29]+alphaR[28])-0.3*(alphaR[20]+alphaR[19])-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*alphaR[5]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])-0.3*(alphaR[29]+alphaR[28])+0.3*(alphaR[26]+alphaR[25])-0.3*(alphaR[20]+alphaR[19])+0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[9]+alphaR[8])+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])-0.81*alphaR[31]-0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])+0.603738353924943*(alphaR[18]+alphaR[17])-0.603738353924943*alphaR[16]+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[10]+0.45*(alphaR[9]+alphaR[8])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[41])+0.375*(alphaR[29]+alphaR[28])-0.3*(alphaR[20]+alphaR[19])+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])+0.81*alphaR[31]-0.3*(alphaR[29]+alphaR[28])+0.3*(alphaR[26]+alphaR[25])-0.3*(alphaR[20]+alphaR[19])-0.603738353924943*(alphaR[18]+alphaR[17])+0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])-0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]+0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*(alphaR[10]+alphaR[8]+alphaR[6])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[20])-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]+0.45*alphaR[6]-0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]-0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]+0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*(alphaR[10]+alphaR[8])+0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])-0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]+0.375*alphaR[20]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]-0.3354101966249678*(alphaR[4]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[20])-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]-0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]+0.375*alphaR[20]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])-0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]+0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[20])-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]-0.45*alphaR[6]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]-0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]+0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[8]+alphaR[6])+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]-0.81*alphaR[31]+0.3*alphaR[29]-0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])+0.3*alphaR[19]+0.603738353924943*alphaR[18]-0.603738353924943*alphaR[17]+0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]-0.375*alphaR[29]+0.375*alphaR[28]-0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]+0.81*alphaR[31]+0.3*alphaR[29]-0.3*alphaR[28]+0.3*(alphaR[26]+alphaR[25])-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[18]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[16]+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[7])+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]+0.3*alphaR[29]-0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])+0.3*alphaR[19]+0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]-0.375*alphaR[29]+0.375*alphaR[28]-0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*alphaR[5]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]+0.3*alphaR[29]-0.3*alphaR[28]+0.3*(alphaR[26]+alphaR[25])-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]+0.81*alphaR[31]+0.3*alphaR[29]-0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])+0.3*alphaR[19]-0.603738353924943*alphaR[18]+0.603738353924943*(alphaR[17]+alphaR[16])-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[10]+alphaR[9])+0.45*(alphaR[8]+alphaR[7])-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]-0.375*alphaR[29]+0.375*alphaR[28]-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]-0.81*alphaR[31]+0.3*alphaR[29]-0.3*alphaR[28]+0.3*(alphaR[26]+alphaR[25])-0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[18]-0.603738353924943*(alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])-0.45*alphaR[8]+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])-0.3*(alphaR[29]+alphaR[26])+0.375*(alphaR[25]+alphaR[19])-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*(alphaR[10]+alphaR[9]+alphaR[7])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[27] = 1.0; 
  else  
    sgn_alpha_surfR[27] = -1.0; 
  
  if (sgn_alpha_surfR[27] == sgn_alpha_surfR[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[29]+alphaR[19])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[28] = 1.0; 
  else  
    sgn_alpha_surfR[28] = -1.0; 
  
  if (sgn_alpha_surfR[28] == sgn_alpha_surfR[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]-0.3*alphaR[29]+0.3*alphaR[26]-0.375*alphaR[25]+0.375*alphaR[19]+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*(alphaR[10]+alphaR[9])+0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[29] = 1.0; 
  else  
    sgn_alpha_surfR[29] = -1.0; 
  
  if (sgn_alpha_surfR[29] == sgn_alpha_surfR[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])-0.3*(alphaR[29]+alphaR[26])+0.375*(alphaR[25]+alphaR[19])+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*alphaR[9]-0.3354101966249678*(alphaR[4]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[30] = 1.0; 
  else  
    sgn_alpha_surfR[30] = -1.0; 
  
  if (sgn_alpha_surfR[30] == sgn_alpha_surfR[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[29]+alphaR[19])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[31] = 1.0; 
  else  
    sgn_alpha_surfR[31] = -1.0; 
  
  if (sgn_alpha_surfR[31] == sgn_alpha_surfR[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]-0.3*alphaR[29]+0.3*alphaR[26]-0.375*alphaR[25]+0.375*alphaR[19]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*alphaR[9]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[32] = 1.0; 
  else  
    sgn_alpha_surfR[32] = -1.0; 
  
  if (sgn_alpha_surfR[32] == sgn_alpha_surfR[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])-0.3*(alphaR[29]+alphaR[26])+0.375*(alphaR[25]+alphaR[19])+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]-0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[33] = 1.0; 
  else  
    sgn_alpha_surfR[33] = -1.0; 
  
  if (sgn_alpha_surfR[33] == sgn_alpha_surfR[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[29]+alphaR[19])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[34] = 1.0; 
  else  
    sgn_alpha_surfR[34] = -1.0; 
  
  if (sgn_alpha_surfR[34] == sgn_alpha_surfR[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]-0.3*alphaR[29]+0.3*alphaR[26]-0.375*alphaR[25]+0.375*alphaR[19]-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[7])+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[35] = 1.0; 
  else  
    sgn_alpha_surfR[35] = -1.0; 
  
  if (sgn_alpha_surfR[35] == sgn_alpha_surfR[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[26]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[36] = 1.0; 
  else  
    sgn_alpha_surfR[36] = -1.0; 
  
  if (sgn_alpha_surfR[36] == sgn_alpha_surfR[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[14]+alphaR[12]+alphaR[11]))-0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[37] = 1.0; 
  else  
    sgn_alpha_surfR[37] = -1.0; 
  
  if (sgn_alpha_surfR[37] == sgn_alpha_surfR[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[26]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])-0.45*alphaR[10]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[38] = 1.0; 
  else  
    sgn_alpha_surfR[38] = -1.0; 
  
  if (sgn_alpha_surfR[38] == sgn_alpha_surfR[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[26]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])-0.3354101966249678*alphaR[4]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[39] = 1.0; 
  else  
    sgn_alpha_surfR[39] = -1.0; 
  
  if (sgn_alpha_surfR[39] == sgn_alpha_surfR[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaR[0]-0.2795084971874732*(alphaR[14]+alphaR[12]+alphaR[11]) > 0.) 
    sgn_alpha_surfR[40] = 1.0; 
  else  
    sgn_alpha_surfR[40] = -1.0; 
  
  if (sgn_alpha_surfR[40] == sgn_alpha_surfR[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[26]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])+0.3354101966249678*alphaR[4]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[41] = 1.0; 
  else  
    sgn_alpha_surfR[41] = -1.0; 
  
  if (sgn_alpha_surfR[41] == sgn_alpha_surfR[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[26]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])-0.45*alphaR[10]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[42] = 1.0; 
  else  
    sgn_alpha_surfR[42] = -1.0; 
  
  if (sgn_alpha_surfR[42] == sgn_alpha_surfR[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[14]+alphaR[12]+alphaR[11]))+0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[43] = 1.0; 
  else  
    sgn_alpha_surfR[43] = -1.0; 
  
  if (sgn_alpha_surfR[43] == sgn_alpha_surfR[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[26]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[12]+alphaR[11])+0.45*alphaR[10]+0.3354101966249678*(alphaR[4]+alphaR[3])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[44] = 1.0; 
  else  
    sgn_alpha_surfR[44] = -1.0; 
  
  if (sgn_alpha_surfR[44] == sgn_alpha_surfR[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]+0.3*alphaR[29]-0.3*alphaR[26]+0.375*alphaR[25]-0.375*alphaR[19]+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[7])-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[45] = 1.0; 
  else  
    sgn_alpha_surfR[45] = -1.0; 
  
  if (sgn_alpha_surfR[45] == sgn_alpha_surfR[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[29]+alphaR[19]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[46] = 1.0; 
  else  
    sgn_alpha_surfR[46] = -1.0; 
  
  if (sgn_alpha_surfR[46] == sgn_alpha_surfR[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])+0.3*(alphaR[29]+alphaR[26])-0.375*(alphaR[25]+alphaR[19])-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]-0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[47] = 1.0; 
  else  
    sgn_alpha_surfR[47] = -1.0; 
  
  if (sgn_alpha_surfR[47] == sgn_alpha_surfR[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]+0.3*alphaR[29]-0.3*alphaR[26]+0.375*alphaR[25]-0.375*alphaR[19]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*alphaR[9]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[48] = 1.0; 
  else  
    sgn_alpha_surfR[48] = -1.0; 
  
  if (sgn_alpha_surfR[48] == sgn_alpha_surfR[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[29]+alphaR[19]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[49] = 1.0; 
  else  
    sgn_alpha_surfR[49] = -1.0; 
  
  if (sgn_alpha_surfR[49] == sgn_alpha_surfR[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])+0.3*(alphaR[29]+alphaR[26])-0.375*(alphaR[25]+alphaR[19])+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*alphaR[9]+0.3354101966249678*(alphaR[4]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[50] = 1.0; 
  else  
    sgn_alpha_surfR[50] = -1.0; 
  
  if (sgn_alpha_surfR[50] == sgn_alpha_surfR[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[35]+0.3*alphaR[29]-0.3*alphaR[26]+0.375*alphaR[25]-0.375*alphaR[19]-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]-0.45*(alphaR[10]+alphaR[9])+0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[51] = 1.0; 
  else  
    sgn_alpha_surfR[51] = -1.0; 
  
  if (sgn_alpha_surfR[51] == sgn_alpha_surfR[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[29]+alphaR[19]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]+0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[52] = 1.0; 
  else  
    sgn_alpha_surfR[52] = -1.0; 
  
  if (sgn_alpha_surfR[52] == sgn_alpha_surfR[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[35])+0.3*(alphaR[29]+alphaR[26])-0.375*(alphaR[25]+alphaR[19])+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[12])-0.2795084971874732*alphaR[11]+0.45*(alphaR[10]+alphaR[9]+alphaR[7])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[53] = 1.0; 
  else  
    sgn_alpha_surfR[53] = -1.0; 
  
  if (sgn_alpha_surfR[53] == sgn_alpha_surfR[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]-0.81*alphaR[31]-0.3*alphaR[29]+0.3*alphaR[28]-0.3*(alphaR[26]+alphaR[25])+0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[18]+0.603738353924943*(alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])-0.45*alphaR[8]+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[54] = 1.0; 
  else  
    sgn_alpha_surfR[54] = -1.0; 
  
  if (sgn_alpha_surfR[54] == sgn_alpha_surfR[53]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]+0.375*alphaR[29]-0.375*alphaR[28]+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[55] = 1.0; 
  else  
    sgn_alpha_surfR[55] = -1.0; 
  
  if (sgn_alpha_surfR[55] == sgn_alpha_surfR[54]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]+0.81*alphaR[31]-0.3*alphaR[29]+0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])-0.3*alphaR[19]+0.603738353924943*alphaR[18]-0.603738353924943*(alphaR[17]+alphaR[16])+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[10]+alphaR[9])+0.45*(alphaR[8]+alphaR[7])-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[56] = 1.0; 
  else  
    sgn_alpha_surfR[56] = -1.0; 
  
  if (sgn_alpha_surfR[56] == sgn_alpha_surfR[55]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]-0.3*alphaR[29]+0.3*alphaR[28]-0.3*(alphaR[26]+alphaR[25])+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[57] = 1.0; 
  else  
    sgn_alpha_surfR[57] = -1.0; 
  
  if (sgn_alpha_surfR[57] == sgn_alpha_surfR[56]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]+0.375*alphaR[29]-0.375*alphaR[28]+0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*alphaR[5]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[58] = 1.0; 
  else  
    sgn_alpha_surfR[58] = -1.0; 
  
  if (sgn_alpha_surfR[58] == sgn_alpha_surfR[57]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]-0.3*alphaR[29]+0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])-0.3*alphaR[19]-0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[59] = 1.0; 
  else  
    sgn_alpha_surfR[59] = -1.0; 
  
  if (sgn_alpha_surfR[59] == sgn_alpha_surfR[58]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[41]+alphaR[36]))+0.4024922359499623*alphaR[35]+0.81*alphaR[31]-0.3*alphaR[29]+0.3*alphaR[28]-0.3*(alphaR[26]+alphaR[25])+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[18]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[16]-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[7])+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[60] = 1.0; 
  else  
    sgn_alpha_surfR[60] = -1.0; 
  
  if (sgn_alpha_surfR[60] == sgn_alpha_surfR[59]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[41]+0.375*alphaR[29]-0.375*alphaR[28]+0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[61] = 1.0; 
  else  
    sgn_alpha_surfR[61] = -1.0; 
  
  if (sgn_alpha_surfR[61] == sgn_alpha_surfR[60]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[41])+0.4024922359499623*alphaR[36]-0.4024922359499623*alphaR[35]-0.81*alphaR[31]-0.3*alphaR[29]+0.3*(alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20])-0.3*alphaR[19]-0.603738353924943*alphaR[18]+0.603738353924943*alphaR[17]-0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[62] = 1.0; 
  else  
    sgn_alpha_surfR[62] = -1.0; 
  
  if (sgn_alpha_surfR[62] == sgn_alpha_surfR[61]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]+0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]-0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[8]+alphaR[6])-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[63] = 1.0; 
  else  
    sgn_alpha_surfR[63] = -1.0; 
  
  if (sgn_alpha_surfR[63] == sgn_alpha_surfR[62]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[20]))-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]-0.45*alphaR[6]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[64] = 1.0; 
  else  
    sgn_alpha_surfR[64] = -1.0; 
  
  if (sgn_alpha_surfR[64] == sgn_alpha_surfR[63]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])+0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]-0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[65] = 1.0; 
  else  
    sgn_alpha_surfR[65] = -1.0; 
  
  if (sgn_alpha_surfR[65] == sgn_alpha_surfR[64]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]+0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]-0.375*alphaR[20]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[66] = 1.0; 
  else  
    sgn_alpha_surfR[66] = -1.0; 
  
  if (sgn_alpha_surfR[66] == sgn_alpha_surfR[65]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[20]))-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[67] = 1.0; 
  else  
    sgn_alpha_surfR[67] = -1.0; 
  
  if (sgn_alpha_surfR[67] == sgn_alpha_surfR[66]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])+0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]-0.375*alphaR[20]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]+0.3354101966249678*(alphaR[4]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[68] = 1.0; 
  else  
    sgn_alpha_surfR[68] = -1.0; 
  
  if (sgn_alpha_surfR[68] == sgn_alpha_surfR[67]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[36]+0.3*alphaR[28]+0.375*alphaR[26]-0.3*alphaR[25]-0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*(alphaR[10]+alphaR[8])+0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[69] = 1.0; 
  else  
    sgn_alpha_surfR[69] = -1.0; 
  
  if (sgn_alpha_surfR[69] == sgn_alpha_surfR[68]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[20]))-0.2795084971874732*(alphaR[14]+alphaR[12])+0.2236067977499786*alphaR[11]+0.45*alphaR[6]+0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[70] = 1.0; 
  else  
    sgn_alpha_surfR[70] = -1.0; 
  
  if (sgn_alpha_surfR[70] == sgn_alpha_surfR[69]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[36])+0.3*alphaR[28]-0.375*alphaR[26]+0.3*alphaR[25]-0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*(alphaR[10]+alphaR[8]+alphaR[6])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[71] = 1.0; 
  else  
    sgn_alpha_surfR[71] = -1.0; 
  
  if (sgn_alpha_surfR[71] == sgn_alpha_surfR[70]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])+0.81*alphaR[31]+0.3*(alphaR[29]+alphaR[28])-0.3*(alphaR[26]+alphaR[25])+0.3*(alphaR[20]+alphaR[19])+0.603738353924943*(alphaR[18]+alphaR[17])-0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[72] = 1.0; 
  else  
    sgn_alpha_surfR[72] = -1.0; 
  
  if (sgn_alpha_surfR[72] == sgn_alpha_surfR[71]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[41])-0.375*(alphaR[29]+alphaR[28])+0.3*(alphaR[20]+alphaR[19])-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[73] = 1.0; 
  else  
    sgn_alpha_surfR[73] = -1.0; 
  
  if (sgn_alpha_surfR[73] == sgn_alpha_surfR[72]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])-0.81*alphaR[31]+0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])-0.603738353924943*(alphaR[18]+alphaR[17])+0.603738353924943*alphaR[16]-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*alphaR[10]+0.45*(alphaR[9]+alphaR[8])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[74] = 1.0; 
  else  
    sgn_alpha_surfR[74] = -1.0; 
  
  if (sgn_alpha_surfR[74] == sgn_alpha_surfR[73]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])+0.3*(alphaR[29]+alphaR[28])-0.3*(alphaR[26]+alphaR[25])+0.3*(alphaR[20]+alphaR[19])-0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[9]+alphaR[8])+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[75] = 1.0; 
  else  
    sgn_alpha_surfR[75] = -1.0; 
  
  if (sgn_alpha_surfR[75] == sgn_alpha_surfR[74]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[41])-0.375*(alphaR[29]+alphaR[28])+0.3*(alphaR[20]+alphaR[19])-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*alphaR[5]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[76] = 1.0; 
  else  
    sgn_alpha_surfR[76] = -1.0; 
  
  if (sgn_alpha_surfR[76] == sgn_alpha_surfR[75]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])+0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])+0.603738353924943*alphaR[16]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[9]+alphaR[8]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[77] = 1.0; 
  else  
    sgn_alpha_surfR[77] = -1.0; 
  
  if (sgn_alpha_surfR[77] == sgn_alpha_surfR[76]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[41]-0.4024922359499623*(alphaR[36]+alphaR[35])-0.81*alphaR[31]+0.3*(alphaR[29]+alphaR[28])-0.3*(alphaR[26]+alphaR[25])+0.3*(alphaR[20]+alphaR[19])-0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16])+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])-0.45*(alphaR[10]+alphaR[9]+alphaR[8])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[78] = 1.0; 
  else  
    sgn_alpha_surfR[78] = -1.0; 
  
  if (sgn_alpha_surfR[78] == sgn_alpha_surfR[77]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[41])-0.375*(alphaR[29]+alphaR[28])+0.3*(alphaR[20]+alphaR[19])+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[12]+alphaR[11])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[79] = 1.0; 
  else  
    sgn_alpha_surfR[79] = -1.0; 
  
  if (sgn_alpha_surfR[79] == sgn_alpha_surfR[78]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[41]+alphaR[36]+alphaR[35])+0.81*alphaR[31]+0.3*(alphaR[29]+alphaR[28]+alphaR[26]+alphaR[25]+alphaR[20]+alphaR[19])+0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[80] = 1.0; 
  else  
    sgn_alpha_surfR[80] = -1.0; 
  
  if (sgn_alpha_surfR[80] == sgn_alpha_surfR[79]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
