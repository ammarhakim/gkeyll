#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[96];
  double *sgn_alpha_surfL = &sgn_alpha_surf[162];
  alphaL[0] = (-2.738612787525831*hamil[35]*rdx2)+2.121320343559642*hamil[7]*rdx2-1.224744871391589*hamil[1]*rdx2; 
  alphaL[1] = 4.743416490252569*hamil[33]*rdx2-2.738612787525831*hamil[16]*rdx2; 
  alphaL[2] = (-2.738612787525831*hamil[58]*rdx2)+2.121320343559642*hamil[21]*rdx2-1.224744871391589*hamil[6]*rdx2; 
  alphaL[3] = (-2.738612787525831*hamil[63]*rdx2)+2.121320343559642*hamil[23]*rdx2-1.224744871391589*hamil[9]*rdx2; 
  alphaL[4] = (-2.738612787525831*hamil[72]*rdx2)+2.121320343559642*hamil[26]*rdx2-1.224744871391589*hamil[12]*rdx2; 
  alphaL[5] = 4.743416490252569*hamil[56]*rdx2-2.738612787525831*hamil[31]*rdx2; 
  alphaL[6] = 4.743416490252569*hamil[61]*rdx2-2.738612787525831*hamil[37]*rdx2; 
  alphaL[7] = (-2.738612787525831*hamil[89]*rdx2)+2.121320343559642*hamil[51]*rdx2-1.224744871391589*hamil[22]*rdx2; 
  alphaL[8] = 4.743416490252569*hamil[70]*rdx2-2.738612787525831*hamil[43]*rdx2; 
  alphaL[9] = (-2.738612787525831*hamil[93]*rdx2)+2.121320343559642*hamil[52]*rdx2-1.224744871391589*hamil[25]*rdx2; 
  alphaL[10] = (-2.738612787525831*hamil[98]*rdx2)+2.121320343559642*hamil[54]*rdx2-1.224744871391589*hamil[28]*rdx2; 
  alphaL[12] = 2.121320343559642*hamil[57]*rdx2-1.224744871391589*hamil[32]*rdx2; 
  alphaL[13] = 2.121320343559642*hamil[66]*rdx2-1.224744871391589*hamil[40]*rdx2; 
  alphaL[14] = 2.121320343559642*hamil[81]*rdx2-1.224744871391589*hamil[47]*rdx2; 
  alphaL[15] = 4.743416490252569*hamil[87]*rdx2-2.738612787525831*hamil[59]*rdx2; 
  alphaL[16] = 4.743416490252569*hamil[91]*rdx2-2.738612787525831*hamil[68]*rdx2; 
  alphaL[17] = 4.743416490252569*hamil[96]*rdx2-2.738612787525831*hamil[74]*rdx2; 
  alphaL[18] = (-2.738612787525831*hamil[109]*rdx2)+2.121320343559642*hamil[86]*rdx2-1.224744871391589*hamil[53]*rdx2; 
  alphaL[22] = 2.121320343559642*hamil[88]*rdx2-1.224744871391589*hamil[60]*rdx2; 
  alphaL[24] = 2.121320343559642*hamil[90]*rdx2-1.224744871391589*hamil[65]*rdx2; 
  alphaL[26] = 2.121320343559642*hamil[92]*rdx2-1.224744871391589*hamil[69]*rdx2; 
  alphaL[27] = 2.121320343559642*hamil[101]*rdx2-1.224744871391589*hamil[77]*rdx2; 
  alphaL[29] = 2.121320343559642*hamil[103]*rdx2-1.224744871391589*hamil[80]*rdx2; 
  alphaL[30] = 2.121320343559642*hamil[105]*rdx2-1.224744871391589*hamil[83]*rdx2; 
  alphaL[31] = 4.743416490252569*hamil[107]*rdx2-2.738612787525831*hamil[94]*rdx2; 
  alphaL[38] = 2.121320343559642*hamil[108]*rdx2-1.224744871391589*hamil[95]*rdx2; 
  alphaL[40] = 2.121320343559642*hamil[110]*rdx2-1.224744871391589*hamil[100]*rdx2; 
  alphaL[43] = 2.121320343559642*hamil[111]*rdx2-1.224744871391589*hamil[104]*rdx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])+0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])-0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.5031152949374518*alphaL[43])+0.375*(alphaL[30]+alphaL[29])-0.3*(alphaL[24]+alphaL[22])-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])-0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[29])+0.3*(alphaL[27]+alphaL[26])-0.3*(alphaL[24]+alphaL[22])+0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16])-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9]+alphaL[8])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])-0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]+0.375*alphaL[24]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*(alphaL[9]+alphaL[8]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[29]+alphaL[24])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]+0.45*alphaL[5]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]-0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]+0.375*alphaL[24]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*(alphaL[9]+alphaL[8])+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]-0.81*alphaL[31]+0.3*alphaL[30]-0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])+0.3*alphaL[22]+0.603738353924943*(alphaL[18]+alphaL[17])-0.603738353924943*alphaL[16]+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*(alphaL[9]+alphaL[8])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]-0.375*alphaL[30]+0.375*alphaL[29]-0.3*alphaL[24]+0.3*alphaL[22]+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]+0.81*alphaL[31]+0.3*alphaL[30]-0.3*alphaL[29]+0.3*(alphaL[27]+alphaL[26])-0.3*alphaL[24]+0.3*alphaL[22]-0.603738353924943*(alphaL[18]+alphaL[17])+0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[26]+alphaL[22])-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*(alphaL[10]+alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[22])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.45*alphaL[6]-0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[26]+0.375*alphaL[22]+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*(alphaL[10]+alphaL[8])+0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[26])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])+0.45*alphaL[8]-0.3354101966249678*(alphaL[4]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[12]))-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[26]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])-0.45*alphaL[8]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[26]-0.375*alphaL[22]+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[22]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.45*alphaL[6]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[26]+alphaL[22])-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*alphaL[10]-0.45*(alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]-0.81*alphaL[31]-0.3*alphaL[30]+0.3*alphaL[29]-0.3*(alphaL[27]+alphaL[26])+0.3*alphaL[24]-0.3*alphaL[22]+0.603738353924943*alphaL[18]-0.603738353924943*alphaL[17]+0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]+0.375*alphaL[30]-0.375*alphaL[29]+0.3*alphaL[24]-0.3*alphaL[22]+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]+0.81*alphaL[31]-0.3*alphaL[30]+0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])-0.3*alphaL[22]-0.603738353924943*alphaL[18]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[16]+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[7])+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]+0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]-0.375*alphaL[24]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[29]+alphaL[24]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]-0.45*alphaL[5]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])+0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]-0.375*alphaL[24]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])+0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[29])-0.3*(alphaL[27]+alphaL[26])+0.3*(alphaL[24]+alphaL[22])-0.603738353924943*alphaL[18]+0.603738353924943*(alphaL[17]+alphaL[16])-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9])+0.45*(alphaL[8]+alphaL[7])-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[43])-0.375*(alphaL[30]+alphaL[29])+0.3*(alphaL[24]+alphaL[22])-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])-0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])+0.603738353924943*alphaL[18]-0.603738353924943*(alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9])-0.45*alphaL[8]+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])-0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[43])+0.375*(alphaL[30]+alphaL[29])-0.3*(alphaL[24]+alphaL[22])-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*alphaL[7]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])-0.3*(alphaL[30]+alphaL[29])+0.3*(alphaL[27]+alphaL[26])-0.3*(alphaL[24]+alphaL[22])+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9])+0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])-0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]+0.375*alphaL[24]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*alphaL[9]-0.3354101966249678*(alphaL[4]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[29]+alphaL[24])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]-0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]+0.375*alphaL[24]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*alphaL[9]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[32] = 1.0; 
  else  
    sgn_alpha_surfL[32] = -1.0; 
  
  if (sgn_alpha_surfL[32] == sgn_alpha_surfL[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]+0.3*alphaL[30]-0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])+0.3*alphaL[22]+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[33] = 1.0; 
  else  
    sgn_alpha_surfL[33] = -1.0; 
  
  if (sgn_alpha_surfL[33] == sgn_alpha_surfL[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]-0.375*alphaL[30]+0.375*alphaL[29]-0.3*alphaL[24]+0.3*alphaL[22]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*alphaL[7]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[34] = 1.0; 
  else  
    sgn_alpha_surfL[34] = -1.0; 
  
  if (sgn_alpha_surfL[34] == sgn_alpha_surfL[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]+0.3*alphaL[30]-0.3*alphaL[29]+0.3*(alphaL[27]+alphaL[26])-0.3*alphaL[24]+0.3*alphaL[22]-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[35] = 1.0; 
  else  
    sgn_alpha_surfL[35] = -1.0; 
  
  if (sgn_alpha_surfL[35] == sgn_alpha_surfL[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[26]+alphaL[22])+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*alphaL[10]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[36] = 1.0; 
  else  
    sgn_alpha_surfL[36] = -1.0; 
  
  if (sgn_alpha_surfL[36] == sgn_alpha_surfL[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[22])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[37] = 1.0; 
  else  
    sgn_alpha_surfL[37] = -1.0; 
  
  if (sgn_alpha_surfL[37] == sgn_alpha_surfL[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[26]+0.375*alphaL[22]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*alphaL[10]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[38] = 1.0; 
  else  
    sgn_alpha_surfL[38] = -1.0; 
  
  if (sgn_alpha_surfL[38] == sgn_alpha_surfL[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[26])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])-0.3354101966249678*alphaL[4]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[39] = 1.0; 
  else  
    sgn_alpha_surfL[39] = -1.0; 
  
  if (sgn_alpha_surfL[39] == sgn_alpha_surfL[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL[0]-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[12]) > 0.) 
    sgn_alpha_surfL[40] = 1.0; 
  else  
    sgn_alpha_surfL[40] = -1.0; 
  
  if (sgn_alpha_surfL[40] == sgn_alpha_surfL[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[26]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])+0.3354101966249678*alphaL[4]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[41] = 1.0; 
  else  
    sgn_alpha_surfL[41] = -1.0; 
  
  if (sgn_alpha_surfL[41] == sgn_alpha_surfL[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[26]-0.375*alphaL[22]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*alphaL[10]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[42] = 1.0; 
  else  
    sgn_alpha_surfL[42] = -1.0; 
  
  if (sgn_alpha_surfL[42] == sgn_alpha_surfL[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[22]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[43] = 1.0; 
  else  
    sgn_alpha_surfL[43] = -1.0; 
  
  if (sgn_alpha_surfL[43] == sgn_alpha_surfL[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[26]+alphaL[22])+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*alphaL[10]+0.3354101966249678*(alphaL[4]+alphaL[3])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[44] = 1.0; 
  else  
    sgn_alpha_surfL[44] = -1.0; 
  
  if (sgn_alpha_surfL[44] == sgn_alpha_surfL[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]-0.3*alphaL[30]+0.3*alphaL[29]-0.3*(alphaL[27]+alphaL[26])+0.3*alphaL[24]-0.3*alphaL[22]+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[45] = 1.0; 
  else  
    sgn_alpha_surfL[45] = -1.0; 
  
  if (sgn_alpha_surfL[45] == sgn_alpha_surfL[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]+0.375*alphaL[30]-0.375*alphaL[29]+0.3*alphaL[24]-0.3*alphaL[22]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*alphaL[7]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[46] = 1.0; 
  else  
    sgn_alpha_surfL[46] = -1.0; 
  
  if (sgn_alpha_surfL[46] == sgn_alpha_surfL[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]-0.3*alphaL[30]+0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])-0.3*alphaL[22]-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[47] = 1.0; 
  else  
    sgn_alpha_surfL[47] = -1.0; 
  
  if (sgn_alpha_surfL[47] == sgn_alpha_surfL[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]+0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]-0.375*alphaL[24]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*alphaL[9]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[48] = 1.0; 
  else  
    sgn_alpha_surfL[48] = -1.0; 
  
  if (sgn_alpha_surfL[48] == sgn_alpha_surfL[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[29]+alphaL[24]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[49] = 1.0; 
  else  
    sgn_alpha_surfL[49] = -1.0; 
  
  if (sgn_alpha_surfL[49] == sgn_alpha_surfL[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])+0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]-0.375*alphaL[24]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*alphaL[9]+0.3354101966249678*(alphaL[4]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[50] = 1.0; 
  else  
    sgn_alpha_surfL[50] = -1.0; 
  
  if (sgn_alpha_surfL[50] == sgn_alpha_surfL[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])+0.3*(alphaL[30]+alphaL[29])-0.3*(alphaL[27]+alphaL[26])+0.3*(alphaL[24]+alphaL[22])-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9])+0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[51] = 1.0; 
  else  
    sgn_alpha_surfL[51] = -1.0; 
  
  if (sgn_alpha_surfL[51] == sgn_alpha_surfL[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[43])-0.375*(alphaL[30]+alphaL[29])+0.3*(alphaL[24]+alphaL[22])-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*alphaL[7]+0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[52] = 1.0; 
  else  
    sgn_alpha_surfL[52] = -1.0; 
  
  if (sgn_alpha_surfL[52] == sgn_alpha_surfL[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])+0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[53] = 1.0; 
  else  
    sgn_alpha_surfL[53] = -1.0; 
  
  if (sgn_alpha_surfL[53] == sgn_alpha_surfL[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])-0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])-0.603738353924943*alphaL[18]+0.603738353924943*(alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9])-0.45*alphaL[8]+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[54] = 1.0; 
  else  
    sgn_alpha_surfL[54] = -1.0; 
  
  if (sgn_alpha_surfL[54] == sgn_alpha_surfL[53]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[43])+0.375*(alphaL[30]+alphaL[29])-0.3*(alphaL[24]+alphaL[22])+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[55] = 1.0; 
  else  
    sgn_alpha_surfL[55] = -1.0; 
  
  if (sgn_alpha_surfL[55] == sgn_alpha_surfL[54]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])+0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[29])+0.3*(alphaL[27]+alphaL[26])-0.3*(alphaL[24]+alphaL[22])+0.603738353924943*alphaL[18]-0.603738353924943*(alphaL[17]+alphaL[16])+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9])+0.45*(alphaL[8]+alphaL[7])-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[56] = 1.0; 
  else  
    sgn_alpha_surfL[56] = -1.0; 
  
  if (sgn_alpha_surfL[56] == sgn_alpha_surfL[55]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])-0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]+0.375*alphaL[24]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[57] = 1.0; 
  else  
    sgn_alpha_surfL[57] = -1.0; 
  
  if (sgn_alpha_surfL[57] == sgn_alpha_surfL[56]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[29]+alphaL[24])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]-0.45*alphaL[5]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[58] = 1.0; 
  else  
    sgn_alpha_surfL[58] = -1.0; 
  
  if (sgn_alpha_surfL[58] == sgn_alpha_surfL[57]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]-0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]+0.375*alphaL[24]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[59] = 1.0; 
  else  
    sgn_alpha_surfL[59] = -1.0; 
  
  if (sgn_alpha_surfL[59] == sgn_alpha_surfL[58]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]+0.81*alphaL[31]+0.3*alphaL[30]-0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])+0.3*alphaL[22]+0.603738353924943*alphaL[18]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[16]-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[7])+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[60] = 1.0; 
  else  
    sgn_alpha_surfL[60] = -1.0; 
  
  if (sgn_alpha_surfL[60] == sgn_alpha_surfL[59]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]-0.375*alphaL[30]+0.375*alphaL[29]-0.3*alphaL[24]+0.3*alphaL[22]-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[61] = 1.0; 
  else  
    sgn_alpha_surfL[61] = -1.0; 
  
  if (sgn_alpha_surfL[61] == sgn_alpha_surfL[60]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]-0.81*alphaL[31]+0.3*alphaL[30]-0.3*alphaL[29]+0.3*(alphaL[27]+alphaL[26])-0.3*alphaL[24]+0.3*alphaL[22]-0.603738353924943*alphaL[18]+0.603738353924943*alphaL[17]-0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[62] = 1.0; 
  else  
    sgn_alpha_surfL[62] = -1.0; 
  
  if (sgn_alpha_surfL[62] == sgn_alpha_surfL[61]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[26]+alphaL[22])+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*alphaL[10]-0.45*(alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[63] = 1.0; 
  else  
    sgn_alpha_surfL[63] = -1.0; 
  
  if (sgn_alpha_surfL[63] == sgn_alpha_surfL[62]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[22])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.45*alphaL[6]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[64] = 1.0; 
  else  
    sgn_alpha_surfL[64] = -1.0; 
  
  if (sgn_alpha_surfL[64] == sgn_alpha_surfL[63]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[26]+0.375*alphaL[22]-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[65] = 1.0; 
  else  
    sgn_alpha_surfL[65] = -1.0; 
  
  if (sgn_alpha_surfL[65] == sgn_alpha_surfL[64]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[26])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])-0.45*alphaL[8]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[66] = 1.0; 
  else  
    sgn_alpha_surfL[66] = -1.0; 
  
  if (sgn_alpha_surfL[66] == sgn_alpha_surfL[65]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[12]))+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[67] = 1.0; 
  else  
    sgn_alpha_surfL[67] = -1.0; 
  
  if (sgn_alpha_surfL[67] == sgn_alpha_surfL[66]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[26]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[12])+0.45*alphaL[8]+0.3354101966249678*(alphaL[4]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[68] = 1.0; 
  else  
    sgn_alpha_surfL[68] = -1.0; 
  
  if (sgn_alpha_surfL[68] == sgn_alpha_surfL[67]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[38]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[26]-0.375*alphaL[22]-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]-0.45*(alphaL[10]+alphaL[8])+0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[69] = 1.0; 
  else  
    sgn_alpha_surfL[69] = -1.0; 
  
  if (sgn_alpha_surfL[69] == sgn_alpha_surfL[68]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[22]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.45*alphaL[6]+0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[70] = 1.0; 
  else  
    sgn_alpha_surfL[70] = -1.0; 
  
  if (sgn_alpha_surfL[70] == sgn_alpha_surfL[69]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[38])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[26]+alphaL[22])+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[12]+0.45*(alphaL[10]+alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[71] = 1.0; 
  else  
    sgn_alpha_surfL[71] = -1.0; 
  
  if (sgn_alpha_surfL[71] == sgn_alpha_surfL[70]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[43]+alphaL[40]))+0.4024922359499623*alphaL[38]+0.81*alphaL[31]-0.3*alphaL[30]+0.3*alphaL[29]-0.3*(alphaL[27]+alphaL[26])+0.3*alphaL[24]-0.3*alphaL[22]+0.603738353924943*(alphaL[18]+alphaL[17])-0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[72] = 1.0; 
  else  
    sgn_alpha_surfL[72] = -1.0; 
  
  if (sgn_alpha_surfL[72] == sgn_alpha_surfL[71]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[43]+0.375*alphaL[30]-0.375*alphaL[29]+0.3*alphaL[24]-0.3*alphaL[22]-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[73] = 1.0; 
  else  
    sgn_alpha_surfL[73] = -1.0; 
  
  if (sgn_alpha_surfL[73] == sgn_alpha_surfL[72]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[43])+0.4024922359499623*alphaL[40]-0.4024922359499623*alphaL[38]-0.81*alphaL[31]-0.3*alphaL[30]+0.3*(alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24])-0.3*alphaL[22]-0.603738353924943*(alphaL[18]+alphaL[17])+0.603738353924943*alphaL[16]-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*alphaL[10]+0.45*(alphaL[9]+alphaL[8])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[74] = 1.0; 
  else  
    sgn_alpha_surfL[74] = -1.0; 
  
  if (sgn_alpha_surfL[74] == sgn_alpha_surfL[73]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[40]+0.3*alphaL[29]+0.375*alphaL[27]-0.3*alphaL[26]-0.375*alphaL[24]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.45*(alphaL[9]+alphaL[8])+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[75] = 1.0; 
  else  
    sgn_alpha_surfL[75] = -1.0; 
  
  if (sgn_alpha_surfL[75] == sgn_alpha_surfL[74]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[29]+alphaL[24]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[12]+0.45*alphaL[5]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[76] = 1.0; 
  else  
    sgn_alpha_surfL[76] = -1.0; 
  
  if (sgn_alpha_surfL[76] == sgn_alpha_surfL[75]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[40])+0.3*alphaL[29]-0.375*alphaL[27]+0.3*alphaL[26]-0.375*alphaL[24]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.45*(alphaL[9]+alphaL[8]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[77] = 1.0; 
  else  
    sgn_alpha_surfL[77] = -1.0; 
  
  if (sgn_alpha_surfL[77] == sgn_alpha_surfL[76]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[43]-0.4024922359499623*(alphaL[40]+alphaL[38])-0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[29])-0.3*(alphaL[27]+alphaL[26])+0.3*(alphaL[24]+alphaL[22])-0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16])+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])-0.45*(alphaL[10]+alphaL[9]+alphaL[8])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[78] = 1.0; 
  else  
    sgn_alpha_surfL[78] = -1.0; 
  
  if (sgn_alpha_surfL[78] == sgn_alpha_surfL[77]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[43])-0.375*(alphaL[30]+alphaL[29])+0.3*(alphaL[24]+alphaL[22])+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[12])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[79] = 1.0; 
  else  
    sgn_alpha_surfL[79] = -1.0; 
  
  if (sgn_alpha_surfL[79] == sgn_alpha_surfL[78]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[43]+alphaL[40]+alphaL[38])+0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[29]+alphaL[27]+alphaL[26]+alphaL[24]+alphaL[22])+0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[12])+0.45*(alphaL[10]+alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[80] = 1.0; 
  else  
    sgn_alpha_surfL[80] = -1.0; 
  
  if (sgn_alpha_surfL[80] == sgn_alpha_surfL[79]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
