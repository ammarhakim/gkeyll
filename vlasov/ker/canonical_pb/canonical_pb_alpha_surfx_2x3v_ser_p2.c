#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
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
  alphaL[0] = 2.738612787525831*hamil[33]*rdvx2-2.121320343559642*hamil[7]*rdvx2+1.224744871391589*hamil[3]*rdvx2; 
  alphaL[1] = 2.738612787525831*hamil[56]*rdvx2-2.121320343559642*hamil[21]*rdvx2+1.224744871391589*hamil[8]*rdvx2; 
  alphaL[2] = 2.738612787525831*hamil[18]*rdvx2-4.743416490252569*hamil[35]*rdvx2; 
  alphaL[3] = 2.738612787525831*hamil[61]*rdvx2-2.121320343559642*hamil[23]*rdvx2+1.224744871391589*hamil[11]*rdvx2; 
  alphaL[4] = 2.738612787525831*hamil[70]*rdvx2-2.121320343559642*hamil[26]*rdvx2+1.224744871391589*hamil[14]*rdvx2; 
  alphaL[5] = 2.738612787525831*hamil[36]*rdvx2-4.743416490252569*hamil[58]*rdvx2; 
  alphaL[6] = 2.738612787525831*hamil[87]*rdvx2-2.121320343559642*hamil[51]*rdvx2+1.224744871391589*hamil[24]*rdvx2; 
  alphaL[7] = 2.738612787525831*hamil[39]*rdvx2-4.743416490252569*hamil[63]*rdvx2; 
  alphaL[8] = 2.738612787525831*hamil[91]*rdvx2-2.121320343559642*hamil[52]*rdvx2+1.224744871391589*hamil[27]*rdvx2; 
  alphaL[9] = 2.738612787525831*hamil[45]*rdvx2-4.743416490252569*hamil[72]*rdvx2; 
  alphaL[10] = 2.738612787525831*hamil[96]*rdvx2-2.121320343559642*hamil[54]*rdvx2+1.224744871391589*hamil[30]*rdvx2; 
  alphaL[11] = 1.224744871391589*hamil[34]*rdvx2-2.121320343559642*hamil[57]*rdvx2; 
  alphaL[13] = 1.224744871391589*hamil[42]*rdvx2-2.121320343559642*hamil[66]*rdvx2; 
  alphaL[14] = 1.224744871391589*hamil[49]*rdvx2-2.121320343559642*hamil[81]*rdvx2; 
  alphaL[15] = 2.738612787525831*hamil[64]*rdvx2-4.743416490252569*hamil[89]*rdvx2; 
  alphaL[16] = 2.738612787525831*hamil[73]*rdvx2-4.743416490252569*hamil[93]*rdvx2; 
  alphaL[17] = 2.738612787525831*hamil[107]*rdvx2-2.121320343559642*hamil[86]*rdvx2+1.224744871391589*hamil[55]*rdvx2; 
  alphaL[18] = 2.738612787525831*hamil[76]*rdvx2-4.743416490252569*hamil[98]*rdvx2; 
  alphaL[21] = 1.224744871391589*hamil[62]*rdvx2-2.121320343559642*hamil[88]*rdvx2; 
  alphaL[23] = 1.224744871391589*hamil[67]*rdvx2-2.121320343559642*hamil[90]*rdvx2; 
  alphaL[25] = 1.224744871391589*hamil[71]*rdvx2-2.121320343559642*hamil[92]*rdvx2; 
  alphaL[27] = 1.224744871391589*hamil[79]*rdvx2-2.121320343559642*hamil[101]*rdvx2; 
  alphaL[28] = 1.224744871391589*hamil[82]*rdvx2-2.121320343559642*hamil[103]*rdvx2; 
  alphaL[30] = 1.224744871391589*hamil[85]*rdvx2-2.121320343559642*hamil[105]*rdvx2; 
  alphaL[31] = 2.738612787525831*hamil[99]*rdvx2-4.743416490252569*hamil[109]*rdvx2; 
  alphaL[37] = 1.224744871391589*hamil[97]*rdvx2-2.121320343559642*hamil[108]*rdvx2; 
  alphaL[39] = 1.224744871391589*hamil[102]*rdvx2-2.121320343559642*hamil[110]*rdvx2; 
  alphaL[42] = 1.224744871391589*hamil[106]*rdvx2-2.121320343559642*hamil[111]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])+0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])-0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.5031152949374518*alphaL[42])+0.375*(alphaL[30]+alphaL[28])-0.3*(alphaL[23]+alphaL[21])-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])-0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[28])+0.3*(alphaL[27]+alphaL[25])-0.3*(alphaL[23]+alphaL[21])+0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16])-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[9]+alphaL[8])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])-0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[23]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*(alphaL[9]+alphaL[8]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[28]+alphaL[23])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]+0.45*alphaL[5]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]-0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[23]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*(alphaL[9]+alphaL[8])+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]-0.81*alphaL[31]+0.3*alphaL[30]-0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])+0.3*alphaL[21]+0.603738353924943*(alphaL[18]+alphaL[17])-0.603738353924943*alphaL[16]+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*(alphaL[9]+alphaL[8])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]-0.375*alphaL[30]+0.375*alphaL[28]-0.3*alphaL[23]+0.3*alphaL[21]+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]+0.81*alphaL[31]+0.3*alphaL[30]-0.3*alphaL[28]+0.3*(alphaL[27]+alphaL[25])-0.3*alphaL[23]+0.3*alphaL[21]-0.603738353924943*(alphaL[18]+alphaL[17])+0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])-0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[42])+0.375*(alphaL[30]+alphaL[28])-0.3*(alphaL[23]+alphaL[21])-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*alphaL[6]-0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])-0.3*(alphaL[30]+alphaL[28])+0.3*(alphaL[27]+alphaL[25])-0.3*(alphaL[23]+alphaL[21])+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[8])+0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])-0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[23]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]-0.3354101966249678*(alphaL[4]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[28]+alphaL[23])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]-0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[23]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.3*alphaL[30]-0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])+0.3*alphaL[21]+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]-0.375*alphaL[30]+0.375*alphaL[28]-0.3*alphaL[23]+0.3*alphaL[21]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*alphaL[6]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]+0.3*alphaL[30]-0.3*alphaL[28]+0.3*(alphaL[27]+alphaL[25])-0.3*alphaL[23]+0.3*alphaL[21]-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*(alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])-0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])+0.603738353924943*alphaL[18]-0.603738353924943*alphaL[17]+0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[42])+0.375*(alphaL[30]+alphaL[28])-0.3*(alphaL[23]+alphaL[21])+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])+0.81*alphaL[31]-0.3*(alphaL[30]+alphaL[28])+0.3*(alphaL[27]+alphaL[25])-0.3*(alphaL[23]+alphaL[21])-0.603738353924943*alphaL[18]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[16]+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[7])+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])-0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[23]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[28]+alphaL[23])-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]-0.45*alphaL[5]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]-0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[23]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.81*alphaL[31]+0.3*alphaL[30]-0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])+0.3*alphaL[21]-0.603738353924943*alphaL[18]+0.603738353924943*(alphaL[17]+alphaL[16])-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[9])+0.45*(alphaL[8]+alphaL[7])-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]-0.375*alphaL[30]+0.375*alphaL[28]-0.3*alphaL[23]+0.3*alphaL[21]-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]-0.81*alphaL[31]+0.3*alphaL[30]-0.3*alphaL[28]+0.3*(alphaL[27]+alphaL[25])-0.3*alphaL[23]+0.3*alphaL[21]+0.603738353924943*alphaL[18]-0.603738353924943*(alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[9])-0.45*alphaL[8]+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])-0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[25]+alphaL[21])-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*(alphaL[10]+alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[21])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]+0.45*alphaL[7]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[25]+0.375*alphaL[21]+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*(alphaL[10]+alphaL[9])+0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[25])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])+0.45*alphaL[9]-0.3354101966249678*(alphaL[4]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[11]))-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[25]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])-0.45*alphaL[9]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[32] = 1.0; 
  else  
    sgn_alpha_surfL[32] = -1.0; 
  
  if (sgn_alpha_surfL[32] == sgn_alpha_surfL[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[25]-0.375*alphaL[21]+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[9]-0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[33] = 1.0; 
  else  
    sgn_alpha_surfL[33] = -1.0; 
  
  if (sgn_alpha_surfL[33] == sgn_alpha_surfL[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[21]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]-0.45*alphaL[7]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[34] = 1.0; 
  else  
    sgn_alpha_surfL[34] = -1.0; 
  
  if (sgn_alpha_surfL[34] == sgn_alpha_surfL[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[25]+alphaL[21])-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[35] = 1.0; 
  else  
    sgn_alpha_surfL[35] = -1.0; 
  
  if (sgn_alpha_surfL[35] == sgn_alpha_surfL[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[25]+alphaL[21])+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[36] = 1.0; 
  else  
    sgn_alpha_surfL[36] = -1.0; 
  
  if (sgn_alpha_surfL[36] == sgn_alpha_surfL[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[21])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[37] = 1.0; 
  else  
    sgn_alpha_surfL[37] = -1.0; 
  
  if (sgn_alpha_surfL[37] == sgn_alpha_surfL[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[25]+0.375*alphaL[21]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[38] = 1.0; 
  else  
    sgn_alpha_surfL[38] = -1.0; 
  
  if (sgn_alpha_surfL[38] == sgn_alpha_surfL[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[25])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])-0.3354101966249678*alphaL[4]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[39] = 1.0; 
  else  
    sgn_alpha_surfL[39] = -1.0; 
  
  if (sgn_alpha_surfL[39] == sgn_alpha_surfL[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL[0]-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[11]) > 0.) 
    sgn_alpha_surfL[40] = 1.0; 
  else  
    sgn_alpha_surfL[40] = -1.0; 
  
  if (sgn_alpha_surfL[40] == sgn_alpha_surfL[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[25]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])+0.3354101966249678*alphaL[4]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[41] = 1.0; 
  else  
    sgn_alpha_surfL[41] = -1.0; 
  
  if (sgn_alpha_surfL[41] == sgn_alpha_surfL[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[25]-0.375*alphaL[21]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*alphaL[10]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[42] = 1.0; 
  else  
    sgn_alpha_surfL[42] = -1.0; 
  
  if (sgn_alpha_surfL[42] == sgn_alpha_surfL[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[21]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[3]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[43] = 1.0; 
  else  
    sgn_alpha_surfL[43] = -1.0; 
  
  if (sgn_alpha_surfL[43] == sgn_alpha_surfL[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[25]+alphaL[21])+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*alphaL[10]+0.3354101966249678*(alphaL[4]+alphaL[3])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[44] = 1.0; 
  else  
    sgn_alpha_surfL[44] = -1.0; 
  
  if (sgn_alpha_surfL[44] == sgn_alpha_surfL[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])-0.3*(alphaL[30]+alphaL[27])+0.375*(alphaL[25]+alphaL[21])+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[7])-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[45] = 1.0; 
  else  
    sgn_alpha_surfL[45] = -1.0; 
  
  if (sgn_alpha_surfL[45] == sgn_alpha_surfL[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[30]+alphaL[21])-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]-0.45*alphaL[7]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[46] = 1.0; 
  else  
    sgn_alpha_surfL[46] = -1.0; 
  
  if (sgn_alpha_surfL[46] == sgn_alpha_surfL[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]-0.3*alphaL[30]+0.3*alphaL[27]-0.375*alphaL[25]+0.375*alphaL[21]-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[9]-0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[47] = 1.0; 
  else  
    sgn_alpha_surfL[47] = -1.0; 
  
  if (sgn_alpha_surfL[47] == sgn_alpha_surfL[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaL[27]+alphaL[25])+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])-0.45*alphaL[9]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[48] = 1.0; 
  else  
    sgn_alpha_surfL[48] = -1.0; 
  
  if (sgn_alpha_surfL[48] == sgn_alpha_surfL[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL[14]+alphaL[13]+alphaL[11]))+0.3354101966249678*alphaL[2]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[49] = 1.0; 
  else  
    sgn_alpha_surfL[49] = -1.0; 
  
  if (sgn_alpha_surfL[49] == sgn_alpha_surfL[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[27]+alphaL[25]))+0.2236067977499786*alphaL[14]-0.2795084971874732*(alphaL[13]+alphaL[11])+0.45*alphaL[9]+0.3354101966249678*(alphaL[4]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[50] = 1.0; 
  else  
    sgn_alpha_surfL[50] = -1.0; 
  
  if (sgn_alpha_surfL[50] == sgn_alpha_surfL[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[37]+0.3*alphaL[30]-0.3*alphaL[27]+0.375*alphaL[25]-0.375*alphaL[21]-0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]-0.45*(alphaL[10]+alphaL[9])+0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[51] = 1.0; 
  else  
    sgn_alpha_surfL[51] = -1.0; 
  
  if (sgn_alpha_surfL[51] == sgn_alpha_surfL[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[30]+alphaL[21]))-0.2795084971874732*alphaL[14]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[11]+0.45*alphaL[7]+0.3354101966249678*(alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[52] = 1.0; 
  else  
    sgn_alpha_surfL[52] = -1.0; 
  
  if (sgn_alpha_surfL[52] == sgn_alpha_surfL[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[37])+0.3*(alphaL[30]+alphaL[27])-0.375*(alphaL[25]+alphaL[21])+0.603738353924943*alphaL[18]+0.2236067977499786*(alphaL[14]+alphaL[13])-0.2795084971874732*alphaL[11]+0.45*(alphaL[10]+alphaL[9]+alphaL[7])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[53] = 1.0; 
  else  
    sgn_alpha_surfL[53] = -1.0; 
  
  if (sgn_alpha_surfL[53] == sgn_alpha_surfL[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]-0.81*alphaL[31]-0.3*alphaL[30]+0.3*alphaL[28]-0.3*(alphaL[27]+alphaL[25])+0.3*alphaL[23]-0.3*alphaL[21]-0.603738353924943*alphaL[18]+0.603738353924943*(alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[9])-0.45*alphaL[8]+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[54] = 1.0; 
  else  
    sgn_alpha_surfL[54] = -1.0; 
  
  if (sgn_alpha_surfL[54] == sgn_alpha_surfL[53]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]+0.375*alphaL[30]-0.375*alphaL[28]+0.3*alphaL[23]-0.3*alphaL[21]+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*alphaL[7]-0.45*(alphaL[6]+alphaL[5])-0.3354101966249678*(alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[55] = 1.0; 
  else  
    sgn_alpha_surfL[55] = -1.0; 
  
  if (sgn_alpha_surfL[55] == sgn_alpha_surfL[54]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.81*alphaL[31]-0.3*alphaL[30]+0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])-0.3*alphaL[21]+0.603738353924943*alphaL[18]-0.603738353924943*(alphaL[17]+alphaL[16])+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[9])+0.45*(alphaL[8]+alphaL[7])-0.45*(alphaL[6]+alphaL[5])+0.3354101966249678*alphaL[4]-0.3354101966249678*(alphaL[3]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[56] = 1.0; 
  else  
    sgn_alpha_surfL[56] = -1.0; 
  
  if (sgn_alpha_surfL[56] == sgn_alpha_surfL[55]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]+0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[23]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[5])-0.3354101966249678*(alphaL[4]+alphaL[2])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[57] = 1.0; 
  else  
    sgn_alpha_surfL[57] = -1.0; 
  
  if (sgn_alpha_surfL[57] == sgn_alpha_surfL[56]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[28]+alphaL[23]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]-0.45*alphaL[5]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[58] = 1.0; 
  else  
    sgn_alpha_surfL[58] = -1.0; 
  
  if (sgn_alpha_surfL[58] == sgn_alpha_surfL[57]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])+0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[23]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[59] = 1.0; 
  else  
    sgn_alpha_surfL[59] = -1.0; 
  
  if (sgn_alpha_surfL[59] == sgn_alpha_surfL[58]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])+0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[28])-0.3*(alphaL[27]+alphaL[25])+0.3*(alphaL[23]+alphaL[21])+0.603738353924943*alphaL[18]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[16]-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*alphaL[9]-0.45*(alphaL[8]+alphaL[7])+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[60] = 1.0; 
  else  
    sgn_alpha_surfL[60] = -1.0; 
  
  if (sgn_alpha_surfL[60] == sgn_alpha_surfL[59]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[42])-0.375*(alphaL[30]+alphaL[28])+0.3*(alphaL[23]+alphaL[21])-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[61] = 1.0; 
  else  
    sgn_alpha_surfL[61] = -1.0; 
  
  if (sgn_alpha_surfL[61] == sgn_alpha_surfL[60]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])-0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])-0.603738353924943*alphaL[18]+0.603738353924943*alphaL[17]-0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*alphaL[9]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*(alphaL[4]+alphaL[3])-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[62] = 1.0; 
  else  
    sgn_alpha_surfL[62] = -1.0; 
  
  if (sgn_alpha_surfL[62] == sgn_alpha_surfL[61]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]-0.3*alphaL[30]+0.3*alphaL[28]-0.3*(alphaL[27]+alphaL[25])+0.3*alphaL[23]-0.3*alphaL[21]+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*(alphaL[8]+alphaL[6])-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[63] = 1.0; 
  else  
    sgn_alpha_surfL[63] = -1.0; 
  
  if (sgn_alpha_surfL[63] == sgn_alpha_surfL[62]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]+0.375*alphaL[30]-0.375*alphaL[28]+0.3*alphaL[23]-0.3*alphaL[21]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*alphaL[6]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[64] = 1.0; 
  else  
    sgn_alpha_surfL[64] = -1.0; 
  
  if (sgn_alpha_surfL[64] == sgn_alpha_surfL[63]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]-0.3*alphaL[30]+0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])-0.3*alphaL[21]-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[65] = 1.0; 
  else  
    sgn_alpha_surfL[65] = -1.0; 
  
  if (sgn_alpha_surfL[65] == sgn_alpha_surfL[64]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]+0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[23]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[66] = 1.0; 
  else  
    sgn_alpha_surfL[66] = -1.0; 
  
  if (sgn_alpha_surfL[66] == sgn_alpha_surfL[65]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[28]+alphaL[23]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]+0.3354101966249678*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[67] = 1.0; 
  else  
    sgn_alpha_surfL[67] = -1.0; 
  
  if (sgn_alpha_surfL[67] == sgn_alpha_surfL[66]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])+0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[23]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]+0.3354101966249678*(alphaL[4]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[68] = 1.0; 
  else  
    sgn_alpha_surfL[68] = -1.0; 
  
  if (sgn_alpha_surfL[68] == sgn_alpha_surfL[67]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])+0.3*(alphaL[30]+alphaL[28])-0.3*(alphaL[27]+alphaL[25])+0.3*(alphaL[23]+alphaL[21])-0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[8])+0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[69] = 1.0; 
  else  
    sgn_alpha_surfL[69] = -1.0; 
  
  if (sgn_alpha_surfL[69] == sgn_alpha_surfL[68]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[42])-0.375*(alphaL[30]+alphaL[28])+0.3*(alphaL[23]+alphaL[21])-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*alphaL[6]+0.3354101966249678*(alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[70] = 1.0; 
  else  
    sgn_alpha_surfL[70] = -1.0; 
  
  if (sgn_alpha_surfL[70] == sgn_alpha_surfL[69]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])+0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])+0.603738353924943*alphaL[17]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[8]+alphaL[6])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[71] = 1.0; 
  else  
    sgn_alpha_surfL[71] = -1.0; 
  
  if (sgn_alpha_surfL[71] == sgn_alpha_surfL[70]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*(alphaL[42]+alphaL[39]))+0.4024922359499623*alphaL[37]+0.81*alphaL[31]-0.3*alphaL[30]+0.3*alphaL[28]-0.3*(alphaL[27]+alphaL[25])+0.3*alphaL[23]-0.3*alphaL[21]+0.603738353924943*(alphaL[18]+alphaL[17])-0.603738353924943*(alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*alphaL[10]-0.45*(alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*(alphaL[4]+alphaL[3])+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[72] = 1.0; 
  else  
    sgn_alpha_surfL[72] = -1.0; 
  
  if (sgn_alpha_surfL[72] == sgn_alpha_surfL[71]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[42]+0.375*alphaL[30]-0.375*alphaL[28]+0.3*alphaL[23]-0.3*alphaL[21]-0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[73] = 1.0; 
  else  
    sgn_alpha_surfL[73] = -1.0; 
  
  if (sgn_alpha_surfL[73] == sgn_alpha_surfL[72]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL[42])+0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]-0.81*alphaL[31]-0.3*alphaL[30]+0.3*(alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23])-0.3*alphaL[21]-0.603738353924943*(alphaL[18]+alphaL[17])+0.603738353924943*alphaL[16]-0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*alphaL[10]+0.45*(alphaL[9]+alphaL[8])-0.45*(alphaL[7]+alphaL[6])+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[74] = 1.0; 
  else  
    sgn_alpha_surfL[74] = -1.0; 
  
  if (sgn_alpha_surfL[74] == sgn_alpha_surfL[73]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL[39]+0.3*alphaL[28]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[23]-0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]-0.45*(alphaL[9]+alphaL[8])+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[75] = 1.0; 
  else  
    sgn_alpha_surfL[75] = -1.0; 
  
  if (sgn_alpha_surfL[75] == sgn_alpha_surfL[74]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaL[28]+alphaL[23]))-0.2795084971874732*(alphaL[14]+alphaL[13])+0.2236067977499786*alphaL[11]+0.45*alphaL[5]+0.3354101966249678*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[76] = 1.0; 
  else  
    sgn_alpha_surfL[76] = -1.0; 
  
  if (sgn_alpha_surfL[76] == sgn_alpha_surfL[75]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[39])+0.3*alphaL[28]-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[23]+0.603738353924943*alphaL[16]+0.2236067977499786*alphaL[14]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[11]+0.45*(alphaL[9]+alphaL[8]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[77] = 1.0; 
  else  
    sgn_alpha_surfL[77] = -1.0; 
  
  if (sgn_alpha_surfL[77] == sgn_alpha_surfL[76]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL[42]-0.4024922359499623*(alphaL[39]+alphaL[37])-0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[28])-0.3*(alphaL[27]+alphaL[25])+0.3*(alphaL[23]+alphaL[21])-0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16])+0.603738353924943*alphaL[15]+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])-0.45*(alphaL[10]+alphaL[9]+alphaL[8])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])-0.3354101966249678*alphaL[4]+0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[78] = 1.0; 
  else  
    sgn_alpha_surfL[78] = -1.0; 
  
  if (sgn_alpha_surfL[78] == sgn_alpha_surfL[77]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL[42])-0.375*(alphaL[30]+alphaL[28])+0.3*(alphaL[23]+alphaL[21])+0.603738353924943*alphaL[15]-0.2795084971874732*alphaL[14]+0.2236067977499786*(alphaL[13]+alphaL[11])+0.45*(alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[79] = 1.0; 
  else  
    sgn_alpha_surfL[79] = -1.0; 
  
  if (sgn_alpha_surfL[79] == sgn_alpha_surfL[78]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*(alphaL[42]+alphaL[39]+alphaL[37])+0.81*alphaL[31]+0.3*(alphaL[30]+alphaL[28]+alphaL[27]+alphaL[25]+alphaL[23]+alphaL[21])+0.603738353924943*(alphaL[18]+alphaL[17]+alphaL[16]+alphaL[15])+0.2236067977499786*(alphaL[14]+alphaL[13]+alphaL[11])+0.45*(alphaL[10]+alphaL[9]+alphaL[8]+alphaL[7]+alphaL[6]+alphaL[5])+0.3354101966249678*(alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[80] = 1.0; 
  else  
    sgn_alpha_surfL[80] = -1.0; 
  
  if (sgn_alpha_surfL[80] == sgn_alpha_surfL[79]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
