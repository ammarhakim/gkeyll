#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
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

  double *alphaR = &alpha_surf[0];
  double *sgn_alpha_surfR = &sgn_alpha_surf[0];
  alphaR[0] = 2.738612787525831*hamil[33]*rdvx2+2.121320343559642*hamil[7]*rdvx2+1.224744871391589*hamil[3]*rdvx2; 
  alphaR[1] = 2.738612787525831*hamil[56]*rdvx2+2.121320343559642*hamil[21]*rdvx2+1.224744871391589*hamil[8]*rdvx2; 
  alphaR[2] = 4.743416490252569*hamil[35]*rdvx2+2.738612787525831*hamil[18]*rdvx2; 
  alphaR[3] = 2.738612787525831*hamil[61]*rdvx2+2.121320343559642*hamil[23]*rdvx2+1.224744871391589*hamil[11]*rdvx2; 
  alphaR[4] = 2.738612787525831*hamil[70]*rdvx2+2.121320343559642*hamil[26]*rdvx2+1.224744871391589*hamil[14]*rdvx2; 
  alphaR[5] = 4.743416490252569*hamil[58]*rdvx2+2.738612787525831*hamil[36]*rdvx2; 
  alphaR[6] = 2.738612787525831*hamil[87]*rdvx2+2.121320343559642*hamil[51]*rdvx2+1.224744871391589*hamil[24]*rdvx2; 
  alphaR[7] = 4.743416490252569*hamil[63]*rdvx2+2.738612787525831*hamil[39]*rdvx2; 
  alphaR[8] = 2.738612787525831*hamil[91]*rdvx2+2.121320343559642*hamil[52]*rdvx2+1.224744871391589*hamil[27]*rdvx2; 
  alphaR[9] = 4.743416490252569*hamil[72]*rdvx2+2.738612787525831*hamil[45]*rdvx2; 
  alphaR[10] = 2.738612787525831*hamil[96]*rdvx2+2.121320343559642*hamil[54]*rdvx2+1.224744871391589*hamil[30]*rdvx2; 
  alphaR[11] = 2.121320343559642*hamil[57]*rdvx2+1.224744871391589*hamil[34]*rdvx2; 
  alphaR[13] = 2.121320343559642*hamil[66]*rdvx2+1.224744871391589*hamil[42]*rdvx2; 
  alphaR[14] = 2.121320343559642*hamil[81]*rdvx2+1.224744871391589*hamil[49]*rdvx2; 
  alphaR[15] = 4.743416490252569*hamil[89]*rdvx2+2.738612787525831*hamil[64]*rdvx2; 
  alphaR[16] = 4.743416490252569*hamil[93]*rdvx2+2.738612787525831*hamil[73]*rdvx2; 
  alphaR[17] = 2.738612787525831*hamil[107]*rdvx2+2.121320343559642*hamil[86]*rdvx2+1.224744871391589*hamil[55]*rdvx2; 
  alphaR[18] = 4.743416490252569*hamil[98]*rdvx2+2.738612787525831*hamil[76]*rdvx2; 
  alphaR[21] = 2.121320343559642*hamil[88]*rdvx2+1.224744871391589*hamil[62]*rdvx2; 
  alphaR[23] = 2.121320343559642*hamil[90]*rdvx2+1.224744871391589*hamil[67]*rdvx2; 
  alphaR[25] = 2.121320343559642*hamil[92]*rdvx2+1.224744871391589*hamil[71]*rdvx2; 
  alphaR[27] = 2.121320343559642*hamil[101]*rdvx2+1.224744871391589*hamil[79]*rdvx2; 
  alphaR[28] = 2.121320343559642*hamil[103]*rdvx2+1.224744871391589*hamil[82]*rdvx2; 
  alphaR[30] = 2.121320343559642*hamil[105]*rdvx2+1.224744871391589*hamil[85]*rdvx2; 
  alphaR[31] = 4.743416490252569*hamil[109]*rdvx2+2.738612787525831*hamil[99]*rdvx2; 
  alphaR[37] = 2.121320343559642*hamil[108]*rdvx2+1.224744871391589*hamil[97]*rdvx2; 
  alphaR[39] = 2.121320343559642*hamil[110]*rdvx2+1.224744871391589*hamil[102]*rdvx2; 
  alphaR[42] = 2.121320343559642*hamil[111]*rdvx2+1.224744871391589*hamil[106]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])+0.81*alphaR[31]-0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])-0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.5031152949374518*alphaR[42])+0.375*(alphaR[30]+alphaR[28])-0.3*(alphaR[23]+alphaR[21])-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])-0.81*alphaR[31]-0.3*(alphaR[30]+alphaR[28])+0.3*(alphaR[27]+alphaR[25])-0.3*(alphaR[23]+alphaR[21])+0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16])-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[9]+alphaR[8])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])-0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[23]-0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*(alphaR[9]+alphaR[8]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[23])-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]+0.45*alphaR[5]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]-0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[23]+0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*(alphaR[9]+alphaR[8])+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]-0.81*alphaR[31]+0.3*alphaR[30]-0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])+0.3*alphaR[21]+0.603738353924943*(alphaR[18]+alphaR[17])-0.603738353924943*alphaR[16]+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*(alphaR[9]+alphaR[8])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]-0.375*alphaR[30]+0.375*alphaR[28]-0.3*alphaR[23]+0.3*alphaR[21]+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]+0.81*alphaR[31]+0.3*alphaR[30]-0.3*alphaR[28]+0.3*(alphaR[27]+alphaR[25])-0.3*alphaR[23]+0.3*alphaR[21]-0.603738353924943*(alphaR[18]+alphaR[17])+0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])-0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])-0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[8]+alphaR[6])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[42])+0.375*(alphaR[30]+alphaR[28])-0.3*(alphaR[23]+alphaR[21])-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*alphaR[6]-0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])-0.3*(alphaR[30]+alphaR[28])+0.3*(alphaR[27]+alphaR[25])-0.3*(alphaR[23]+alphaR[21])+0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[8])+0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])-0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[23]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]-0.3354101966249678*(alphaR[4]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[23])-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]-0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[23]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.3*alphaR[30]-0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])+0.3*alphaR[21]+0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]-0.375*alphaR[30]+0.375*alphaR[28]-0.3*alphaR[23]+0.3*alphaR[21]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*alphaR[6]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]+0.3*alphaR[30]-0.3*alphaR[28]+0.3*(alphaR[27]+alphaR[25])-0.3*alphaR[23]+0.3*alphaR[21]-0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[8]+alphaR[6])+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])-0.81*alphaR[31]-0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])+0.603738353924943*alphaR[18]-0.603738353924943*alphaR[17]+0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[42])+0.375*(alphaR[30]+alphaR[28])-0.3*(alphaR[23]+alphaR[21])+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])+0.81*alphaR[31]-0.3*(alphaR[30]+alphaR[28])+0.3*(alphaR[27]+alphaR[25])-0.3*(alphaR[23]+alphaR[21])-0.603738353924943*alphaR[18]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[16]+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[7])+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])-0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[23]+0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[28]+alphaR[23])-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]-0.45*alphaR[5]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]-0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[23]-0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.81*alphaR[31]+0.3*alphaR[30]-0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])+0.3*alphaR[21]-0.603738353924943*alphaR[18]+0.603738353924943*(alphaR[17]+alphaR[16])-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[9])+0.45*(alphaR[8]+alphaR[7])-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]-0.375*alphaR[30]+0.375*alphaR[28]-0.3*alphaR[23]+0.3*alphaR[21]-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]-0.81*alphaR[31]+0.3*alphaR[30]-0.3*alphaR[28]+0.3*(alphaR[27]+alphaR[25])-0.3*alphaR[23]+0.3*alphaR[21]+0.603738353924943*alphaR[18]-0.603738353924943*(alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])-0.45*alphaR[8]+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])-0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])-0.3*(alphaR[30]+alphaR[27])+0.375*(alphaR[25]+alphaR[21])-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*(alphaR[10]+alphaR[9]+alphaR[7])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[27] = 1.0; 
  else  
    sgn_alpha_surfR[27] = -1.0; 
  
  if (sgn_alpha_surfR[27] == sgn_alpha_surfR[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[30]+alphaR[21])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[28] = 1.0; 
  else  
    sgn_alpha_surfR[28] = -1.0; 
  
  if (sgn_alpha_surfR[28] == sgn_alpha_surfR[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]-0.3*alphaR[30]+0.3*alphaR[27]-0.375*alphaR[25]+0.375*alphaR[21]+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*(alphaR[10]+alphaR[9])+0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[29] = 1.0; 
  else  
    sgn_alpha_surfR[29] = -1.0; 
  
  if (sgn_alpha_surfR[29] == sgn_alpha_surfR[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[27]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])+0.45*alphaR[9]-0.3354101966249678*(alphaR[4]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[30] = 1.0; 
  else  
    sgn_alpha_surfR[30] = -1.0; 
  
  if (sgn_alpha_surfR[30] == sgn_alpha_surfR[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[14]+alphaR[13]+alphaR[11]))-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[31] = 1.0; 
  else  
    sgn_alpha_surfR[31] = -1.0; 
  
  if (sgn_alpha_surfR[31] == sgn_alpha_surfR[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[27]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])-0.45*alphaR[9]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[32] = 1.0; 
  else  
    sgn_alpha_surfR[32] = -1.0; 
  
  if (sgn_alpha_surfR[32] == sgn_alpha_surfR[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]+0.3*alphaR[30]-0.3*alphaR[27]+0.375*alphaR[25]-0.375*alphaR[21]+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]-0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[33] = 1.0; 
  else  
    sgn_alpha_surfR[33] = -1.0; 
  
  if (sgn_alpha_surfR[33] == sgn_alpha_surfR[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[30]+alphaR[21]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[34] = 1.0; 
  else  
    sgn_alpha_surfR[34] = -1.0; 
  
  if (sgn_alpha_surfR[34] == sgn_alpha_surfR[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])+0.3*(alphaR[30]+alphaR[27])-0.375*(alphaR[25]+alphaR[21])-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[7])+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[35] = 1.0; 
  else  
    sgn_alpha_surfR[35] = -1.0; 
  
  if (sgn_alpha_surfR[35] == sgn_alpha_surfR[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])-0.3*(alphaR[30]+alphaR[27])+0.375*(alphaR[25]+alphaR[21])+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[36] = 1.0; 
  else  
    sgn_alpha_surfR[36] = -1.0; 
  
  if (sgn_alpha_surfR[36] == sgn_alpha_surfR[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[30]+alphaR[21])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[37] = 1.0; 
  else  
    sgn_alpha_surfR[37] = -1.0; 
  
  if (sgn_alpha_surfR[37] == sgn_alpha_surfR[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]-0.3*alphaR[30]+0.3*alphaR[27]-0.375*alphaR[25]+0.375*alphaR[21]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[38] = 1.0; 
  else  
    sgn_alpha_surfR[38] = -1.0; 
  
  if (sgn_alpha_surfR[38] == sgn_alpha_surfR[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[27]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])-0.3354101966249678*alphaR[4]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[39] = 1.0; 
  else  
    sgn_alpha_surfR[39] = -1.0; 
  
  if (sgn_alpha_surfR[39] == sgn_alpha_surfR[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaR[0]-0.2795084971874732*(alphaR[14]+alphaR[13]+alphaR[11]) > 0.) 
    sgn_alpha_surfR[40] = 1.0; 
  else  
    sgn_alpha_surfR[40] = -1.0; 
  
  if (sgn_alpha_surfR[40] == sgn_alpha_surfR[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[27]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])+0.3354101966249678*alphaR[4]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[41] = 1.0; 
  else  
    sgn_alpha_surfR[41] = -1.0; 
  
  if (sgn_alpha_surfR[41] == sgn_alpha_surfR[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]+0.3*alphaR[30]-0.3*alphaR[27]+0.375*alphaR[25]-0.375*alphaR[21]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[42] = 1.0; 
  else  
    sgn_alpha_surfR[42] = -1.0; 
  
  if (sgn_alpha_surfR[42] == sgn_alpha_surfR[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[30]+alphaR[21]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[3]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[43] = 1.0; 
  else  
    sgn_alpha_surfR[43] = -1.0; 
  
  if (sgn_alpha_surfR[43] == sgn_alpha_surfR[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])+0.3*(alphaR[30]+alphaR[27])-0.375*(alphaR[25]+alphaR[21])+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]+0.3354101966249678*(alphaR[4]+alphaR[3])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[44] = 1.0; 
  else  
    sgn_alpha_surfR[44] = -1.0; 
  
  if (sgn_alpha_surfR[44] == sgn_alpha_surfR[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])-0.3*(alphaR[30]+alphaR[27])+0.375*(alphaR[25]+alphaR[21])+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[7])-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[45] = 1.0; 
  else  
    sgn_alpha_surfR[45] = -1.0; 
  
  if (sgn_alpha_surfR[45] == sgn_alpha_surfR[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[30]+alphaR[21])-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[46] = 1.0; 
  else  
    sgn_alpha_surfR[46] = -1.0; 
  
  if (sgn_alpha_surfR[46] == sgn_alpha_surfR[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]-0.3*alphaR[30]+0.3*alphaR[27]-0.375*alphaR[25]+0.375*alphaR[21]-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]-0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[47] = 1.0; 
  else  
    sgn_alpha_surfR[47] = -1.0; 
  
  if (sgn_alpha_surfR[47] == sgn_alpha_surfR[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[27]+alphaR[25])+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])-0.45*alphaR[9]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[48] = 1.0; 
  else  
    sgn_alpha_surfR[48] = -1.0; 
  
  if (sgn_alpha_surfR[48] == sgn_alpha_surfR[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[14]+alphaR[13]+alphaR[11]))+0.3354101966249678*alphaR[2]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[49] = 1.0; 
  else  
    sgn_alpha_surfR[49] = -1.0; 
  
  if (sgn_alpha_surfR[49] == sgn_alpha_surfR[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[27]+alphaR[25]))+0.2236067977499786*alphaR[14]-0.2795084971874732*(alphaR[13]+alphaR[11])+0.45*alphaR[9]+0.3354101966249678*(alphaR[4]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[50] = 1.0; 
  else  
    sgn_alpha_surfR[50] = -1.0; 
  
  if (sgn_alpha_surfR[50] == sgn_alpha_surfR[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[37]+0.3*alphaR[30]-0.3*alphaR[27]+0.375*alphaR[25]-0.375*alphaR[21]-0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]-0.45*(alphaR[10]+alphaR[9])+0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[51] = 1.0; 
  else  
    sgn_alpha_surfR[51] = -1.0; 
  
  if (sgn_alpha_surfR[51] == sgn_alpha_surfR[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[30]+alphaR[21]))-0.2795084971874732*alphaR[14]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]+0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[52] = 1.0; 
  else  
    sgn_alpha_surfR[52] = -1.0; 
  
  if (sgn_alpha_surfR[52] == sgn_alpha_surfR[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[37])+0.3*(alphaR[30]+alphaR[27])-0.375*(alphaR[25]+alphaR[21])+0.603738353924943*alphaR[18]+0.2236067977499786*(alphaR[14]+alphaR[13])-0.2795084971874732*alphaR[11]+0.45*(alphaR[10]+alphaR[9]+alphaR[7])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[53] = 1.0; 
  else  
    sgn_alpha_surfR[53] = -1.0; 
  
  if (sgn_alpha_surfR[53] == sgn_alpha_surfR[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]-0.81*alphaR[31]-0.3*alphaR[30]+0.3*alphaR[28]-0.3*(alphaR[27]+alphaR[25])+0.3*alphaR[23]-0.3*alphaR[21]-0.603738353924943*alphaR[18]+0.603738353924943*(alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])-0.45*alphaR[8]+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[54] = 1.0; 
  else  
    sgn_alpha_surfR[54] = -1.0; 
  
  if (sgn_alpha_surfR[54] == sgn_alpha_surfR[53]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]+0.375*alphaR[30]-0.375*alphaR[28]+0.3*alphaR[23]-0.3*alphaR[21]+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*alphaR[7]-0.45*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[55] = 1.0; 
  else  
    sgn_alpha_surfR[55] = -1.0; 
  
  if (sgn_alpha_surfR[55] == sgn_alpha_surfR[54]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.81*alphaR[31]-0.3*alphaR[30]+0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])-0.3*alphaR[21]+0.603738353924943*alphaR[18]-0.603738353924943*(alphaR[17]+alphaR[16])+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[9])+0.45*(alphaR[8]+alphaR[7])-0.45*(alphaR[6]+alphaR[5])+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[56] = 1.0; 
  else  
    sgn_alpha_surfR[56] = -1.0; 
  
  if (sgn_alpha_surfR[56] == sgn_alpha_surfR[55]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]+0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[23]+0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[2])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[57] = 1.0; 
  else  
    sgn_alpha_surfR[57] = -1.0; 
  
  if (sgn_alpha_surfR[57] == sgn_alpha_surfR[56]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[23]))-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]-0.45*alphaR[5]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[58] = 1.0; 
  else  
    sgn_alpha_surfR[58] = -1.0; 
  
  if (sgn_alpha_surfR[58] == sgn_alpha_surfR[57]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])+0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[23]-0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[59] = 1.0; 
  else  
    sgn_alpha_surfR[59] = -1.0; 
  
  if (sgn_alpha_surfR[59] == sgn_alpha_surfR[58]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])+0.81*alphaR[31]+0.3*(alphaR[30]+alphaR[28])-0.3*(alphaR[27]+alphaR[25])+0.3*(alphaR[23]+alphaR[21])+0.603738353924943*alphaR[18]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[16]-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[9]-0.45*(alphaR[8]+alphaR[7])+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[60] = 1.0; 
  else  
    sgn_alpha_surfR[60] = -1.0; 
  
  if (sgn_alpha_surfR[60] == sgn_alpha_surfR[59]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[42])-0.375*(alphaR[30]+alphaR[28])+0.3*(alphaR[23]+alphaR[21])-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[61] = 1.0; 
  else  
    sgn_alpha_surfR[61] = -1.0; 
  
  if (sgn_alpha_surfR[61] == sgn_alpha_surfR[60]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])-0.81*alphaR[31]+0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])-0.603738353924943*alphaR[18]+0.603738353924943*alphaR[17]-0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[62] = 1.0; 
  else  
    sgn_alpha_surfR[62] = -1.0; 
  
  if (sgn_alpha_surfR[62] == sgn_alpha_surfR[61]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]-0.3*alphaR[30]+0.3*alphaR[28]-0.3*(alphaR[27]+alphaR[25])+0.3*alphaR[23]-0.3*alphaR[21]+0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[8]+alphaR[6])-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[63] = 1.0; 
  else  
    sgn_alpha_surfR[63] = -1.0; 
  
  if (sgn_alpha_surfR[63] == sgn_alpha_surfR[62]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]+0.375*alphaR[30]-0.375*alphaR[28]+0.3*alphaR[23]-0.3*alphaR[21]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*alphaR[6]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[64] = 1.0; 
  else  
    sgn_alpha_surfR[64] = -1.0; 
  
  if (sgn_alpha_surfR[64] == sgn_alpha_surfR[63]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]-0.3*alphaR[30]+0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])-0.3*alphaR[21]-0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[65] = 1.0; 
  else  
    sgn_alpha_surfR[65] = -1.0; 
  
  if (sgn_alpha_surfR[65] == sgn_alpha_surfR[64]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]+0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[23]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[66] = 1.0; 
  else  
    sgn_alpha_surfR[66] = -1.0; 
  
  if (sgn_alpha_surfR[66] == sgn_alpha_surfR[65]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[23]))-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]+0.3354101966249678*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[67] = 1.0; 
  else  
    sgn_alpha_surfR[67] = -1.0; 
  
  if (sgn_alpha_surfR[67] == sgn_alpha_surfR[66]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])+0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[23]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]+0.3354101966249678*(alphaR[4]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[68] = 1.0; 
  else  
    sgn_alpha_surfR[68] = -1.0; 
  
  if (sgn_alpha_surfR[68] == sgn_alpha_surfR[67]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])+0.3*(alphaR[30]+alphaR[28])-0.3*(alphaR[27]+alphaR[25])+0.3*(alphaR[23]+alphaR[21])-0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[8])+0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[69] = 1.0; 
  else  
    sgn_alpha_surfR[69] = -1.0; 
  
  if (sgn_alpha_surfR[69] == sgn_alpha_surfR[68]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[42])-0.375*(alphaR[30]+alphaR[28])+0.3*(alphaR[23]+alphaR[21])-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*alphaR[6]+0.3354101966249678*(alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[70] = 1.0; 
  else  
    sgn_alpha_surfR[70] = -1.0; 
  
  if (sgn_alpha_surfR[70] == sgn_alpha_surfR[69]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])+0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])+0.603738353924943*alphaR[17]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[8]+alphaR[6])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[71] = 1.0; 
  else  
    sgn_alpha_surfR[71] = -1.0; 
  
  if (sgn_alpha_surfR[71] == sgn_alpha_surfR[70]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaR[42]+alphaR[39]))+0.4024922359499623*alphaR[37]+0.81*alphaR[31]-0.3*alphaR[30]+0.3*alphaR[28]-0.3*(alphaR[27]+alphaR[25])+0.3*alphaR[23]-0.3*alphaR[21]+0.603738353924943*(alphaR[18]+alphaR[17])-0.603738353924943*(alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*alphaR[10]-0.45*(alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[72] = 1.0; 
  else  
    sgn_alpha_surfR[72] = -1.0; 
  
  if (sgn_alpha_surfR[72] == sgn_alpha_surfR[71]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[42]+0.375*alphaR[30]-0.375*alphaR[28]+0.3*alphaR[23]-0.3*alphaR[21]-0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[73] = 1.0; 
  else  
    sgn_alpha_surfR[73] = -1.0; 
  
  if (sgn_alpha_surfR[73] == sgn_alpha_surfR[72]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaR[42])+0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]-0.81*alphaR[31]-0.3*alphaR[30]+0.3*(alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23])-0.3*alphaR[21]-0.603738353924943*(alphaR[18]+alphaR[17])+0.603738353924943*alphaR[16]-0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*alphaR[10]+0.45*(alphaR[9]+alphaR[8])-0.45*(alphaR[7]+alphaR[6])+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[74] = 1.0; 
  else  
    sgn_alpha_surfR[74] = -1.0; 
  
  if (sgn_alpha_surfR[74] == sgn_alpha_surfR[73]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaR[39]+0.3*alphaR[28]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[23]-0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]-0.45*(alphaR[9]+alphaR[8])+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[75] = 1.0; 
  else  
    sgn_alpha_surfR[75] = -1.0; 
  
  if (sgn_alpha_surfR[75] == sgn_alpha_surfR[74]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[28]+alphaR[23]))-0.2795084971874732*(alphaR[14]+alphaR[13])+0.2236067977499786*alphaR[11]+0.45*alphaR[5]+0.3354101966249678*(alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[76] = 1.0; 
  else  
    sgn_alpha_surfR[76] = -1.0; 
  
  if (sgn_alpha_surfR[76] == sgn_alpha_surfR[75]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[39])+0.3*alphaR[28]-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[23]+0.603738353924943*alphaR[16]+0.2236067977499786*alphaR[14]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[11]+0.45*(alphaR[9]+alphaR[8]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[77] = 1.0; 
  else  
    sgn_alpha_surfR[77] = -1.0; 
  
  if (sgn_alpha_surfR[77] == sgn_alpha_surfR[76]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaR[42]-0.4024922359499623*(alphaR[39]+alphaR[37])-0.81*alphaR[31]+0.3*(alphaR[30]+alphaR[28])-0.3*(alphaR[27]+alphaR[25])+0.3*(alphaR[23]+alphaR[21])-0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16])+0.603738353924943*alphaR[15]+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])-0.45*(alphaR[10]+alphaR[9]+alphaR[8])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[78] = 1.0; 
  else  
    sgn_alpha_surfR[78] = -1.0; 
  
  if (sgn_alpha_surfR[78] == sgn_alpha_surfR[77]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaR[42])-0.375*(alphaR[30]+alphaR[28])+0.3*(alphaR[23]+alphaR[21])+0.603738353924943*alphaR[15]-0.2795084971874732*alphaR[14]+0.2236067977499786*(alphaR[13]+alphaR[11])+0.45*(alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[79] = 1.0; 
  else  
    sgn_alpha_surfR[79] = -1.0; 
  
  if (sgn_alpha_surfR[79] == sgn_alpha_surfR[78]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaR[42]+alphaR[39]+alphaR[37])+0.81*alphaR[31]+0.3*(alphaR[30]+alphaR[28]+alphaR[27]+alphaR[25]+alphaR[23]+alphaR[21])+0.603738353924943*(alphaR[18]+alphaR[17]+alphaR[16]+alphaR[15])+0.2236067977499786*(alphaR[14]+alphaR[13]+alphaR[11])+0.45*(alphaR[10]+alphaR[9]+alphaR[8]+alphaR[7]+alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1])+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[80] = 1.0; 
  else  
    sgn_alpha_surfR[80] = -1.0; 
  
  if (sgn_alpha_surfR[80] == sgn_alpha_surfR[79]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
