#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
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

  double *alphaR = &alpha_surf[40];
  double *sgn_alpha_surfR = &sgn_alpha_surf[54];
  alphaR[0] = 2.121320343559642*hamil[10]*rdvy2+1.224744871391589*hamil[4]*rdvy2; 
  alphaR[1] = 2.121320343559642*hamil[17]*rdvy2+1.224744871391589*hamil[9]*rdvy2; 
  alphaR[2] = 2.121320343559642*hamil[19]*rdvy2+1.224744871391589*hamil[11]*rdvy2; 
  alphaR[3] = 4.743416490252569*hamil[50]*rdvy2+2.738612787525831*hamil[48]*rdvy2; 
  alphaR[4] = 2.121320343559642*hamil[24]*rdvy2+1.224744871391589*hamil[15]*rdvy2; 
  alphaR[5] = 2.121320343559642*hamil[26]*rdvy2+1.224744871391589*hamil[18]*rdvy2; 
  alphaR[6] = 4.743416490252569*hamil[53]*rdvy2+2.738612787525831*hamil[49]*rdvy2; 
  alphaR[7] = 4.743416490252569*hamil[55]*rdvy2+2.738612787525831*hamil[51]*rdvy2; 
  alphaR[8] = 2.121320343559642*hamil[28]*rdvy2+1.224744871391589*hamil[23]*rdvy2; 
  alphaR[9] = 2.121320343559642*hamil[30]*rdvy2+1.224744871391589*hamil[25]*rdvy2; 
  alphaR[10] = 4.743416490252569*hamil[57]*rdvy2+2.738612787525831*hamil[52]*rdvy2; 
  alphaR[11] = 4.743416490252569*hamil[59]*rdvy2+2.738612787525831*hamil[54]*rdvy2; 
  alphaR[12] = 2.121320343559642*hamil[31]*rdvy2+1.224744871391589*hamil[29]*rdvy2; 
  alphaR[13] = 4.743416490252569*hamil[60]*rdvy2+2.738612787525831*hamil[56]*rdvy2; 
  alphaR[14] = 4.743416490252569*hamil[62]*rdvy2+2.738612787525831*hamil[58]*rdvy2; 
  alphaR[15] = 4.743416490252569*hamil[63]*rdvy2+2.738612787525831*hamil[61]*rdvy2; 
  alphaR[16] = 2.121320343559642*hamil[39]*rdvy2+1.224744871391589*hamil[35]*rdvy2; 
  alphaR[17] = 2.121320343559642*hamil[43]*rdvy2+1.224744871391589*hamil[38]*rdvy2; 
  alphaR[19] = 2.121320343559642*hamil[46]*rdvy2+1.224744871391589*hamil[42]*rdvy2; 
  alphaR[21] = 2.121320343559642*hamil[47]*rdvy2+1.224744871391589*hamil[45]*rdvy2; 
  alphaR[32] = 2.121320343559642*hamil[73]*rdvy2+1.224744871391589*hamil[68]*rdvy2; 
  alphaR[33] = 2.121320343559642*hamil[76]*rdvy2+1.224744871391589*hamil[72]*rdvy2; 
  alphaR[34] = 2.121320343559642*hamil[78]*rdvy2+1.224744871391589*hamil[74]*rdvy2; 
  alphaR[36] = 2.121320343559642*hamil[79]*rdvy2+1.224744871391589*hamil[77]*rdvy2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.603738353924943*alphaR[15]-0.603738353924943*alphaR[14]-0.45*(alphaR[13]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])+0.3354101966249678*alphaR[8]+0.45*alphaR[7]+0.3354101966249678*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.375*alphaR[36])+0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.45*alphaR[11]+0.45*alphaR[7]+0.3354101966249678*(alphaR[6]+alphaR[5])-0.3354101966249678*(alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.603738353924943*alphaR[15]+0.603738353924943*alphaR[14]+0.45*(alphaR[13]+alphaR[12])-0.45*(alphaR[11]+alphaR[10]+alphaR[9])-0.3354101966249678*alphaR[8]+0.45*alphaR[7]+0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[4])-0.3354101966249678*(alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.45*alphaR[12]+0.45*alphaR[9]+0.3354101966249678*(alphaR[8]+alphaR[5])-0.3354101966249678*(alphaR[4]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaR[36])+0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.3354101966249678*alphaR[5]-0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.45*alphaR[12]-0.45*alphaR[9]-0.3354101966249678*alphaR[8]+0.3354101966249678*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.603738353924943*alphaR[15]+0.603738353924943*alphaR[14]+0.45*alphaR[13]-0.45*alphaR[12]+0.45*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]+0.3354101966249678*alphaR[8]-0.45*alphaR[7]-0.3354101966249678*alphaR[6]+0.3354101966249678*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaR[36])+0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.45*alphaR[11]-0.45*alphaR[7]-0.3354101966249678*alphaR[6]+0.3354101966249678*(alphaR[5]+alphaR[3])-0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*alphaR[36]-0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.603738353924943*alphaR[15]-0.603738353924943*alphaR[14]-0.45*alphaR[13]+0.45*(alphaR[12]+alphaR[11]+alphaR[10])-0.45*alphaR[9]-0.3354101966249678*alphaR[8]-0.45*alphaR[7]-0.3354101966249678*alphaR[6]+0.3354101966249678*(alphaR[5]+alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]-0.375*alphaR[21]+0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.45*alphaR[13]+0.45*alphaR[10]+0.3354101966249678*(alphaR[8]+alphaR[6])-0.3354101966249678*(alphaR[4]+alphaR[3])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]+0.375*alphaR[21]-0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.45*alphaR[13]-0.45*alphaR[10]-0.3354101966249678*alphaR[8]+0.3354101966249678*(alphaR[6]+alphaR[4])-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]-0.375*alphaR[21]+0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.3354101966249678*alphaR[8]-0.3354101966249678*alphaR[4]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]+0.375*alphaR[21]-0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.3354101966249678*alphaR[8]+0.3354101966249678*alphaR[4]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]-0.375*alphaR[21]+0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.45*alphaR[13]-0.45*alphaR[10]+0.3354101966249678*alphaR[8]-0.3354101966249678*(alphaR[6]+alphaR[4])+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.3354101966249678*alphaR[6]+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR[33])+0.2236067977499786*alphaR[32]+0.375*alphaR[21]-0.375*alphaR[19]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.45*alphaR[13]+0.45*alphaR[10]-0.3354101966249678*(alphaR[8]+alphaR[6])+0.3354101966249678*(alphaR[4]+alphaR[3])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.603738353924943*alphaR[15]+0.603738353924943*alphaR[14]-0.45*alphaR[13]+0.45*(alphaR[12]+alphaR[11]+alphaR[10])-0.45*alphaR[9]+0.3354101966249678*alphaR[8]-0.45*alphaR[7]+0.3354101966249678*alphaR[6]-0.3354101966249678*(alphaR[5]+alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaR[36]-0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.45*alphaR[11]-0.45*alphaR[7]+0.3354101966249678*alphaR[6]-0.3354101966249678*(alphaR[5]+alphaR[3])+0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.603738353924943*alphaR[15]-0.603738353924943*alphaR[14]+0.45*alphaR[13]-0.45*alphaR[12]+0.45*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[9]-0.3354101966249678*alphaR[8]-0.45*alphaR[7]+0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.45*alphaR[12]-0.45*alphaR[9]+0.3354101966249678*alphaR[8]-0.3354101966249678*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaR[36]-0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.3354101966249678*alphaR[5]+0.3354101966249678*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.45*alphaR[12]+0.45*alphaR[9]-0.3354101966249678*(alphaR[8]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]+0.3*alphaR[21]-0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.603738353924943*alphaR[15]-0.603738353924943*alphaR[14]+0.45*(alphaR[13]+alphaR[12])-0.45*(alphaR[11]+alphaR[10]+alphaR[9])+0.3354101966249678*alphaR[8]+0.45*alphaR[7]-0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[4])+0.3354101966249678*(alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaR[36]-0.375*alphaR[34]+0.2795084971874732*alphaR[33]-0.2795084971874732*alphaR[32]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.45*alphaR[11]+0.45*alphaR[7]-0.3354101966249678*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*alphaR[36])+0.3*alphaR[34]-0.2236067977499786*alphaR[33]+0.2236067977499786*alphaR[32]-0.3*alphaR[21]+0.3*alphaR[19]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.603738353924943*alphaR[15]+0.603738353924943*alphaR[14]-0.45*(alphaR[13]+alphaR[12]+alphaR[11])+0.45*(alphaR[10]+alphaR[9])-0.3354101966249678*alphaR[8]+0.45*alphaR[7]-0.3354101966249678*(alphaR[6]+alphaR[5])+0.3354101966249678*(alphaR[4]+alphaR[3]+alphaR[2])-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.603738353924943*(alphaR[15]+alphaR[14])+0.45*(alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9])-0.3354101966249678*alphaR[8]+0.45*alphaR[7]-0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[27] = 1.0; 
  else  
    sgn_alpha_surfR[27] = -1.0; 
  
  if (sgn_alpha_surfR[27] == sgn_alpha_surfR[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[36]+alphaR[34])-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.45*(alphaR[11]+alphaR[7])-0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[28] = 1.0; 
  else  
    sgn_alpha_surfR[28] = -1.0; 
  
  if (sgn_alpha_surfR[28] == sgn_alpha_surfR[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.603738353924943*(alphaR[15]+alphaR[14])-0.45*(alphaR[13]+alphaR[12])+0.45*alphaR[11]-0.45*(alphaR[10]+alphaR[9])+0.3354101966249678*alphaR[8]+0.45*alphaR[7]-0.3354101966249678*(alphaR[6]+alphaR[5])+0.3354101966249678*alphaR[4]-0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[29] = 1.0; 
  else  
    sgn_alpha_surfR[29] = -1.0; 
  
  if (sgn_alpha_surfR[29] == sgn_alpha_surfR[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.45*(alphaR[12]+alphaR[9])-0.3354101966249678*(alphaR[8]+alphaR[5]+alphaR[4]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[30] = 1.0; 
  else  
    sgn_alpha_surfR[30] = -1.0; 
  
  if (sgn_alpha_surfR[30] == sgn_alpha_surfR[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[36]+alphaR[34])-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.3354101966249678*(alphaR[5]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[31] = 1.0; 
  else  
    sgn_alpha_surfR[31] = -1.0; 
  
  if (sgn_alpha_surfR[31] == sgn_alpha_surfR[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.45*(alphaR[12]+alphaR[9])+0.3354101966249678*alphaR[8]-0.3354101966249678*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[32] = 1.0; 
  else  
    sgn_alpha_surfR[32] = -1.0; 
  
  if (sgn_alpha_surfR[32] == sgn_alpha_surfR[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.603738353924943*(alphaR[15]+alphaR[14])-0.45*alphaR[13]+0.45*alphaR[12]-0.45*(alphaR[11]+alphaR[10])+0.45*alphaR[9]-0.3354101966249678*alphaR[8]-0.45*alphaR[7]+0.3354101966249678*alphaR[6]-0.3354101966249678*(alphaR[5]+alphaR[4])+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[33] = 1.0; 
  else  
    sgn_alpha_surfR[33] = -1.0; 
  
  if (sgn_alpha_surfR[33] == sgn_alpha_surfR[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*(alphaR[36]+alphaR[34])-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.45*(alphaR[11]+alphaR[7])+0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[34] = 1.0; 
  else  
    sgn_alpha_surfR[34] = -1.0; 
  
  if (sgn_alpha_surfR[34] == sgn_alpha_surfR[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3*(alphaR[36]+alphaR[34]))+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.603738353924943*(alphaR[15]+alphaR[14])+0.45*alphaR[13]-0.45*(alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]+0.3354101966249678*alphaR[8]-0.45*alphaR[7]+0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[5]+0.3354101966249678*(alphaR[4]+alphaR[3])-0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[35] = 1.0; 
  else  
    sgn_alpha_surfR[35] = -1.0; 
  
  if (sgn_alpha_surfR[35] == sgn_alpha_surfR[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])+0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])+0.45*(alphaR[13]+alphaR[10])-0.3354101966249678*(alphaR[8]+alphaR[6]+alphaR[4]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[36] = 1.0; 
  else  
    sgn_alpha_surfR[36] = -1.0; 
  
  if (sgn_alpha_surfR[36] == sgn_alpha_surfR[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[33]+alphaR[32]+alphaR[17]+alphaR[16]))-0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[37] = 1.0; 
  else  
    sgn_alpha_surfR[37] = -1.0; 
  
  if (sgn_alpha_surfR[37] == sgn_alpha_surfR[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])-0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])-0.45*(alphaR[13]+alphaR[10])+0.3354101966249678*alphaR[8]-0.3354101966249678*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[38] = 1.0; 
  else  
    sgn_alpha_surfR[38] = -1.0; 
  
  if (sgn_alpha_surfR[38] == sgn_alpha_surfR[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])+0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])-0.3354101966249678*(alphaR[8]+alphaR[4])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[39] = 1.0; 
  else  
    sgn_alpha_surfR[39] = -1.0; 
  
  if (sgn_alpha_surfR[39] == sgn_alpha_surfR[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0])-0.2795084971874732*(alphaR[33]+alphaR[32]+alphaR[17]+alphaR[16]) > 0.) 
    sgn_alpha_surfR[40] = 1.0; 
  else  
    sgn_alpha_surfR[40] = -1.0; 
  
  if (sgn_alpha_surfR[40] == sgn_alpha_surfR[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])-0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])+0.3354101966249678*(alphaR[8]+alphaR[4])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[41] = 1.0; 
  else  
    sgn_alpha_surfR[41] = -1.0; 
  
  if (sgn_alpha_surfR[41] == sgn_alpha_surfR[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])+0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])-0.45*(alphaR[13]+alphaR[10])-0.3354101966249678*alphaR[8]+0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[42] = 1.0; 
  else  
    sgn_alpha_surfR[42] = -1.0; 
  
  if (sgn_alpha_surfR[42] == sgn_alpha_surfR[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaR[33]+alphaR[32]+alphaR[17]+alphaR[16]))+0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[43] = 1.0; 
  else  
    sgn_alpha_surfR[43] = -1.0; 
  
  if (sgn_alpha_surfR[43] == sgn_alpha_surfR[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaR[33]+alphaR[32])-0.375*(alphaR[21]+alphaR[19])-0.2795084971874732*(alphaR[17]+alphaR[16])+0.45*(alphaR[13]+alphaR[10])+0.3354101966249678*(alphaR[8]+alphaR[6]+alphaR[4]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[44] = 1.0; 
  else  
    sgn_alpha_surfR[44] = -1.0; 
  
  if (sgn_alpha_surfR[44] == sgn_alpha_surfR[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.603738353924943*(alphaR[15]+alphaR[14])+0.45*alphaR[13]-0.45*(alphaR[12]+alphaR[11])+0.45*alphaR[10]-0.45*alphaR[9]-0.3354101966249678*alphaR[8]-0.45*alphaR[7]-0.3354101966249678*alphaR[6]+0.3354101966249678*alphaR[5]-0.3354101966249678*(alphaR[4]+alphaR[3])+0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[45] = 1.0; 
  else  
    sgn_alpha_surfR[45] = -1.0; 
  
  if (sgn_alpha_surfR[45] == sgn_alpha_surfR[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[36]+alphaR[34]))-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.45*(alphaR[11]+alphaR[7])-0.3354101966249678*alphaR[6]+0.3354101966249678*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[46] = 1.0; 
  else  
    sgn_alpha_surfR[46] = -1.0; 
  
  if (sgn_alpha_surfR[46] == sgn_alpha_surfR[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.603738353924943*(alphaR[15]+alphaR[14])-0.45*alphaR[13]+0.45*alphaR[12]-0.45*(alphaR[11]+alphaR[10])+0.45*alphaR[9]+0.3354101966249678*alphaR[8]-0.45*alphaR[7]-0.3354101966249678*alphaR[6]+0.3354101966249678*(alphaR[5]+alphaR[4])-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[47] = 1.0; 
  else  
    sgn_alpha_surfR[47] = -1.0; 
  
  if (sgn_alpha_surfR[47] == sgn_alpha_surfR[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.45*(alphaR[12]+alphaR[9])-0.3354101966249678*alphaR[8]+0.3354101966249678*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[48] = 1.0; 
  else  
    sgn_alpha_surfR[48] = -1.0; 
  
  if (sgn_alpha_surfR[48] == sgn_alpha_surfR[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[36]+alphaR[34]))-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.3354101966249678*(alphaR[5]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[49] = 1.0; 
  else  
    sgn_alpha_surfR[49] = -1.0; 
  
  if (sgn_alpha_surfR[49] == sgn_alpha_surfR[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.45*(alphaR[12]+alphaR[9])+0.3354101966249678*(alphaR[8]+alphaR[5]+alphaR[4]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[50] = 1.0; 
  else  
    sgn_alpha_surfR[50] = -1.0; 
  
  if (sgn_alpha_surfR[50] == sgn_alpha_surfR[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])-0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])-0.603738353924943*(alphaR[15]+alphaR[14])-0.45*(alphaR[13]+alphaR[12])+0.45*alphaR[11]-0.45*(alphaR[10]+alphaR[9])-0.3354101966249678*alphaR[8]+0.45*alphaR[7]+0.3354101966249678*(alphaR[6]+alphaR[5])-0.3354101966249678*alphaR[4]+0.3354101966249678*(alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[51] = 1.0; 
  else  
    sgn_alpha_surfR[51] = -1.0; 
  
  if (sgn_alpha_surfR[51] == sgn_alpha_surfR[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*(alphaR[36]+alphaR[34]))-0.2795084971874732*(alphaR[33]+alphaR[32])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.45*(alphaR[11]+alphaR[7])+0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[52] = 1.0; 
  else  
    sgn_alpha_surfR[52] = -1.0; 
  
  if (sgn_alpha_surfR[52] == sgn_alpha_surfR[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3*(alphaR[36]+alphaR[34])+0.2236067977499786*(alphaR[33]+alphaR[32])+0.3*(alphaR[21]+alphaR[19])+0.2236067977499786*(alphaR[17]+alphaR[16])+0.603738353924943*(alphaR[15]+alphaR[14])+0.45*(alphaR[13]+alphaR[12]+alphaR[11]+alphaR[10]+alphaR[9])+0.3354101966249678*alphaR[8]+0.45*alphaR[7]+0.3354101966249678*(alphaR[6]+alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[53] = 1.0; 
  else  
    sgn_alpha_surfR[53] = -1.0; 
  
  if (sgn_alpha_surfR[53] == sgn_alpha_surfR[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
