#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_tensor_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x3v_tensor_p2(const double *w, const double *dxv, const double *hamil,
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
  double wvx = w[1];
  double rdvx2 = 2.0/dxv[1];
  double wvy = w[2];
  double rdvy2 = 2.0/dxv[2];
  double wvz = w[3];
  double rdvz2 = 2.0/dxv[3];

  double *alphaL = &alpha_surf[27];
  double *sgn_alpha_surfL = &sgn_alpha_surf[27];
  alphaL[0] = (-2.738612787525831*hamil[20]*rdx2)+2.121320343559642*hamil[5]*rdx2-1.224744871391589*hamil[1]*rdx2; 
  alphaL[1] = (-6.123724356957944*hamil[44]*rdx2)+4.743416490252569*hamil[19]*rdx2-2.738612787525831*hamil[11]*rdx2; 
  alphaL[2] = (-2.738612787525831*hamil[33]*rdx2)+2.121320343559642*hamil[15]*rdx2-1.224744871391589*hamil[6]*rdx2; 
  alphaL[3] = (-2.738612787525831*hamil[36]*rdx2)+2.121320343559642*hamil[16]*rdx2-1.224744871391589*hamil[8]*rdx2; 
  alphaL[4] = (-6.123724356957944*hamil[54]*rdx2)+4.743416490252569*hamil[32]*rdx2-2.738612787525831*hamil[21]*rdx2; 
  alphaL[5] = (-6.123724356957944*hamil[57]*rdx2)+4.743416490252569*hamil[35]*rdx2-2.738612787525831*hamil[25]*rdx2; 
  alphaL[6] = (-2.738612787525831*hamil[51]*rdx2)+2.121320343559642*hamil[31]*rdx2-1.224744871391589*hamil[17]*rdx2; 
  alphaL[8] = (-2.738612787525831*hamil[56]*rdx2)+2.121320343559642*hamil[34]*rdx2-1.224744871391589*hamil[23]*rdx2; 
  alphaL[9] = (-2.738612787525831*hamil[61]*rdx2)+2.121320343559642*hamil[41]*rdx2-1.224744871391589*hamil[28]*rdx2; 
  alphaL[10] = (-6.123724356957944*hamil[66]*rdx2)+4.743416490252569*hamil[50]*rdx2-2.738612787525831*hamil[37]*rdx2; 
  alphaL[12] = (-6.123724356957945*hamil[72]*rdx2)+4.743416490252569*hamil[55]*rdx2-2.738612787525831*hamil[45]*rdx2; 
  alphaL[14] = (-2.738612787525831*hamil[68]*rdx2)+2.121320343559642*hamil[52]*rdx2-1.224744871391589*hamil[39]*rdx2; 
  alphaL[15] = (-6.123724356957945*hamil[73]*rdx2)+4.743416490252569*hamil[60]*rdx2-2.738612787525831*hamil[47]*rdx2; 
  alphaL[16] = (-2.738612787525831*hamil[70]*rdx2)+2.121320343559642*hamil[53]*rdx2-1.224744871391589*hamil[42]*rdx2; 
  alphaL[18] = (-6.123724356957945*hamil[76]*rdx2)+4.743416490252569*hamil[67]*rdx2-2.738612787525831*hamil[58]*rdx2; 
  alphaL[19] = (-6.123724356957945*hamil[77]*rdx2)+4.743416490252569*hamil[69]*rdx2-2.738612787525831*hamil[62]*rdx2; 
  alphaL[22] = (-2.738612787525831*hamil[79]*rdx2)+2.121320343559642*hamil[71]*rdx2-1.224744871391589*hamil[64]*rdx2; 
  alphaL[25] = (-6.123724356957944*hamil[80]*rdx2)+4.743416490252569*hamil[78]*rdx2-2.738612787525831*hamil[74]*rdx2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.3794733192202045*alphaL[25])+0.2828427124746191*alphaL[22]+0.5692099788303081*(alphaL[19]+alphaL[18])-0.4242640687119282*(alphaL[16]+alphaL[15])-0.4242640687119286*alphaL[14]-0.4242640687119282*alphaL[12]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.4743416490252561*alphaL[25]-0.3535533905932734*alphaL[22]-0.7115124735378848*alphaL[19]+0.5303300858899102*(alphaL[16]+alphaL[15])-0.4242640687119282*alphaL[12]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202045*alphaL[25])+0.2828427124746191*alphaL[22]+0.5692099788303081*alphaL[19]-0.5692099788303081*alphaL[18]-0.4242640687119282*(alphaL[16]+alphaL[15])+0.4242640687119286*alphaL[14]-0.4242640687119282*alphaL[12]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252561*alphaL[25]-0.3535533905932734*alphaL[22]-0.7115124735378848*alphaL[18]-0.4242640687119282*alphaL[15]+0.5303300858899102*(alphaL[14]+alphaL[12])+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.6363961030678927*alphaL[5]-0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.592927061281571*alphaL[25])+0.4419417382415923*alphaL[22]+0.5303300858899102*(alphaL[15]+alphaL[12])-0.3952847075210471*(alphaL[9]+alphaL[8])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252561*alphaL[25]-0.3535533905932734*alphaL[22]+0.7115124735378848*alphaL[18]-0.4242640687119282*alphaL[15]-0.5303300858899102*alphaL[14]+0.5303300858899102*alphaL[12]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.6363961030678927*alphaL[5]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202045*alphaL[25])+0.2828427124746191*alphaL[22]-0.5692099788303081*alphaL[19]+0.5692099788303081*alphaL[18]+0.4242640687119282*alphaL[16]-0.4242640687119282*alphaL[15]-0.4242640687119286*alphaL[14]-0.4242640687119282*alphaL[12]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252561*alphaL[25]-0.3535533905932734*alphaL[22]+0.7115124735378848*alphaL[19]-0.5303300858899102*alphaL[16]+0.5303300858899102*alphaL[15]-0.4242640687119282*alphaL[12]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202045*alphaL[25])+0.2828427124746191*alphaL[22]-0.5692099788303081*(alphaL[19]+alphaL[18])+0.4242640687119282*alphaL[16]-0.4242640687119282*alphaL[15]+0.4242640687119286*alphaL[14]-0.4242640687119282*alphaL[12]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[22]-0.4242640687119282*alphaL[16]-0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[22])+0.5303300858899102*alphaL[16]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[22]-0.4242640687119282*alphaL[16]+0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[22])+0.5303300858899102*alphaL[14]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4419417382415923*alphaL[22]-0.3952847075210471*(alphaL[9]+alphaL[8])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[22])-0.5303300858899102*alphaL[14]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[22]+0.4242640687119282*alphaL[16]-0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[22])-0.5303300858899102*alphaL[16]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[22]+0.4242640687119282*alphaL[16]+0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]+0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202045*alphaL[25]+0.2828427124746191*alphaL[22]-0.5692099788303081*(alphaL[19]+alphaL[18])-0.4242640687119282*alphaL[16]+0.4242640687119282*alphaL[15]-0.4242640687119286*alphaL[14]+0.4242640687119282*alphaL[12]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252561*alphaL[25])-0.3535533905932734*alphaL[22]+0.7115124735378848*alphaL[19]+0.5303300858899102*alphaL[16]-0.5303300858899102*alphaL[15]+0.4242640687119282*alphaL[12]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202045*alphaL[25]+0.2828427124746191*alphaL[22]-0.5692099788303081*alphaL[19]+0.5692099788303081*alphaL[18]-0.4242640687119282*alphaL[16]+0.4242640687119282*alphaL[15]+0.4242640687119286*alphaL[14]+0.4242640687119282*alphaL[12]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252561*alphaL[25])-0.3535533905932734*alphaL[22]+0.7115124735378848*alphaL[18]+0.4242640687119282*alphaL[15]+0.5303300858899102*alphaL[14]-0.5303300858899102*alphaL[12]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.6363961030678927*alphaL[5]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.592927061281571*alphaL[25]+0.4419417382415923*alphaL[22]-0.5303300858899102*(alphaL[15]+alphaL[12])-0.3952847075210471*(alphaL[9]+alphaL[8])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252561*alphaL[25])-0.3535533905932734*alphaL[22]-0.7115124735378848*alphaL[18]+0.4242640687119282*alphaL[15]-0.5303300858899102*(alphaL[14]+alphaL[12])+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.6363961030678927*alphaL[5]+0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202045*alphaL[25]+0.2828427124746191*alphaL[22]+0.5692099788303081*alphaL[19]-0.5692099788303081*alphaL[18]+0.4242640687119282*(alphaL[16]+alphaL[15])-0.4242640687119286*alphaL[14]+0.4242640687119282*alphaL[12]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252561*alphaL[25])-0.3535533905932734*alphaL[22]-0.7115124735378848*alphaL[19]-0.5303300858899102*(alphaL[16]+alphaL[15])+0.4242640687119282*alphaL[12]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202045*alphaL[25]+0.2828427124746191*alphaL[22]+0.5692099788303081*(alphaL[19]+alphaL[18])+0.4242640687119282*(alphaL[16]+alphaL[15])+0.4242640687119286*alphaL[14]+0.4242640687119282*alphaL[12]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
