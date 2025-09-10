#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_tensor_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_tensor_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_tensor_p2(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[27];
  double *sgn_alpha_surfL = &sgn_alpha_surf[27];
  alphaL[0] = 2.738612787525831*hamil[26]*rdvy2-2.121320343559642*hamil[9]*rdvy2+1.224744871391589*hamil[4]*rdvy2; 
  alphaL[1] = 2.738612787525831*hamil[36]*rdvy2-2.121320343559642*hamil[16]*rdvy2+1.224744871391589*hamil[8]*rdvy2; 
  alphaL[2] = 2.738612787525831*hamil[38]*rdvy2-2.121320343559642*hamil[18]*rdvy2+1.224744871391589*hamil[10]*rdvy2; 
  alphaL[3] = 6.123724356957944*hamil[48]*rdvy2-4.743416490252569*hamil[29]*rdvy2+2.738612787525831*hamil[14]*rdvy2; 
  alphaL[4] = 2.738612787525831*hamil[51]*rdvy2-2.121320343559642*hamil[31]*rdvy2+1.224744871391589*hamil[17]*rdvy2; 
  alphaL[5] = 6.123724356957944*hamil[61]*rdvy2-4.743416490252569*hamil[41]*rdvy2+2.738612787525831*hamil[28]*rdvy2; 
  alphaL[6] = 6.123724356957944*hamil[63]*rdvy2-4.743416490252569*hamil[43]*rdvy2+2.738612787525831*hamil[30]*rdvy2; 
  alphaL[7] = 2.738612787525831*hamil[57]*rdvy2-2.121320343559642*hamil[35]*rdvy2+1.224744871391589*hamil[25]*rdvy2; 
  alphaL[8] = 2.738612787525831*hamil[59]*rdvy2-2.121320343559642*hamil[40]*rdvy2+1.224744871391589*hamil[27]*rdvy2; 
  alphaL[10] = 6.123724356957944*hamil[70]*rdvy2-4.743416490252569*hamil[53]*rdvy2+2.738612787525831*hamil[42]*rdvy2; 
  alphaL[11] = 2.738612787525831*hamil[66]*rdvy2-2.121320343559642*hamil[50]*rdvy2+1.224744871391589*hamil[37]*rdvy2; 
  alphaL[12] = 2.738612787525831*hamil[68]*rdvy2-2.121320343559642*hamil[52]*rdvy2+1.224744871391589*hamil[39]*rdvy2; 
  alphaL[13] = 6.123724356957945*hamil[73]*rdvy2-4.743416490252569*hamil[60]*rdvy2+2.738612787525831*hamil[47]*rdvy2; 
  alphaL[14] = 6.123724356957945*hamil[75]*rdvy2-4.743416490252569*hamil[65]*rdvy2+2.738612787525831*hamil[49]*rdvy2; 
  alphaL[17] = 6.123724356957945*hamil[77]*rdvy2-4.743416490252569*hamil[69]*rdvy2+2.738612787525831*hamil[62]*rdvy2; 
  alphaL[18] = 6.123724356957945*hamil[79]*rdvy2-4.743416490252569*hamil[71]*rdvy2+2.738612787525831*hamil[64]*rdvy2; 
  alphaL[20] = 2.738612787525831*hamil[76]*rdvy2-2.121320343559642*hamil[67]*rdvy2+1.224744871391589*hamil[58]*rdvy2; 
  alphaL[23] = 6.123724356957944*hamil[80]*rdvy2-4.743416490252569*hamil[78]*rdvy2+2.738612787525831*hamil[74]*rdvy2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.3794733192202054*alphaL[23])+0.2828427124746191*alphaL[20]+0.5692099788303081*(alphaL[18]+alphaL[17])-0.4242640687119286*(alphaL[14]+alphaL[13])-0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.2828427124746191*alphaL[20]-0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202054*alphaL[23]+0.2828427124746191*alphaL[20]-0.5692099788303081*(alphaL[18]+alphaL[17])+0.4242640687119286*(alphaL[14]+alphaL[13])-0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaL[23]-0.3535533905932734*alphaL[20]-0.7115124735378848*alphaL[18]+0.5303300858899102*alphaL[14]-0.4242640687119286*alphaL[13]+0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]-0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[20])+0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252568*alphaL[23])-0.3535533905932734*alphaL[20]+0.7115124735378848*alphaL[18]-0.5303300858899102*alphaL[14]+0.4242640687119286*alphaL[13]+0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202054*alphaL[23])+0.2828427124746191*alphaL[20]+0.5692099788303081*alphaL[18]-0.5692099788303081*alphaL[17]-0.4242640687119286*(alphaL[14]+alphaL[13])-0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[20]-0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202054*alphaL[23]+0.2828427124746191*alphaL[20]-0.5692099788303081*alphaL[18]+0.5692099788303081*alphaL[17]+0.4242640687119286*(alphaL[14]+alphaL[13])-0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaL[23]-0.3535533905932734*alphaL[20]-0.7115124735378848*alphaL[17]-0.4242640687119286*alphaL[14]+0.5303300858899102*(alphaL[13]+alphaL[11])+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.6363961030678927*alphaL[6]-0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[20])+0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252568*alphaL[23])-0.3535533905932734*alphaL[20]+0.7115124735378848*alphaL[17]+0.4242640687119286*alphaL[14]-0.5303300858899102*alphaL[13]+0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.6363961030678927*alphaL[6]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.592927061281571*alphaL[23])+0.4419417382415923*alphaL[20]+0.5303300858899102*(alphaL[14]+alphaL[13])-0.3952847075210471*(alphaL[8]+alphaL[7])-0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4419417382415923*alphaL[20]-0.3952847075210471*(alphaL[8]+alphaL[7])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.592927061281571*alphaL[23]+0.4419417382415923*alphaL[20]-0.5303300858899102*(alphaL[14]+alphaL[13])-0.3952847075210471*(alphaL[8]+alphaL[7])+0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaL[23]-0.3535533905932734*alphaL[20]+0.7115124735378848*alphaL[17]-0.4242640687119286*alphaL[14]+0.5303300858899102*alphaL[13]-0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.6363961030678927*alphaL[6]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[20])-0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252568*alphaL[23])-0.3535533905932734*alphaL[20]-0.7115124735378848*alphaL[17]+0.4242640687119286*alphaL[14]-0.5303300858899102*(alphaL[13]+alphaL[11])+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.6363961030678927*alphaL[6]+0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202054*alphaL[23])+0.2828427124746191*alphaL[20]-0.5692099788303081*alphaL[18]+0.5692099788303081*alphaL[17]-0.4242640687119286*(alphaL[14]+alphaL[13])+0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[20]+0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202054*alphaL[23]+0.2828427124746191*alphaL[20]+0.5692099788303081*alphaL[18]-0.5692099788303081*alphaL[17]+0.4242640687119286*(alphaL[14]+alphaL[13])+0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaL[23]-0.3535533905932734*alphaL[20]+0.7115124735378848*alphaL[18]+0.5303300858899102*alphaL[14]-0.4242640687119286*alphaL[13]-0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3535533905932734*alphaL[20])-0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4743416490252568*alphaL[23])-0.3535533905932734*alphaL[20]-0.7115124735378848*alphaL[18]-0.5303300858899102*alphaL[14]+0.4242640687119286*alphaL[13]-0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]+0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3794733192202054*alphaL[23])+0.2828427124746191*alphaL[20]-0.5692099788303081*(alphaL[18]+alphaL[17])-0.4242640687119286*(alphaL[14]+alphaL[13])+0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2828427124746191*alphaL[20]+0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3794733192202054*alphaL[23]+0.2828427124746191*alphaL[20]+0.5692099788303081*(alphaL[18]+alphaL[17])+0.4242640687119286*(alphaL[14]+alphaL[13])+0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[8]+alphaL[7])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
