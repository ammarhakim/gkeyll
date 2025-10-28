#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *hamil,
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

  double *alphaL = &alpha_surf[0];
  double *sgn_alpha_surfL = &sgn_alpha_surf[0];
  alphaL[0] = 2.738612787525831*hamil[19]*rdvx2-2.121320343559642*hamil[5]*rdvx2+1.224744871391589*hamil[2]*rdvx2; 
  alphaL[1] = 2.738612787525831*hamil[12]*rdvx2-4.743416490252569*hamil[20]*rdvx2; 
  alphaL[2] = 2.738612787525831*hamil[32]*rdvx2-2.121320343559642*hamil[15]*rdvx2+1.224744871391589*hamil[7]*rdvx2; 
  alphaL[3] = 2.738612787525831*hamil[35]*rdvx2-2.121320343559642*hamil[16]*rdvx2+1.224744871391589*hamil[9]*rdvx2; 
  alphaL[4] = 2.738612787525831*hamil[22]*rdvx2-4.743416490252569*hamil[33]*rdvx2; 
  alphaL[5] = 2.738612787525831*hamil[26]*rdvx2-4.743416490252569*hamil[36]*rdvx2; 
  alphaL[6] = 2.738612787525831*hamil[44]*rdvx2-2.121320343559642*hamil[31]*rdvx2+1.224744871391589*hamil[18]*rdvx2; 
  alphaL[8] = 1.224744871391589*hamil[24]*rdvx2-2.121320343559642*hamil[34]*rdvx2; 
  alphaL[9] = 1.224744871391589*hamil[29]*rdvx2-2.121320343559642*hamil[41]*rdvx2; 
  alphaL[10] = 2.738612787525831*hamil[38]*rdvx2-4.743416490252569*hamil[45]*rdvx2; 
  alphaL[14] = 1.224744871391589*hamil[40]*rdvx2-2.121320343559642*hamil[46]*rdvx2; 
  alphaL[16] = 1.224744871391589*hamil[43]*rdvx2-2.121320343559642*hamil[47]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.4242640687119282*alphaL[16])-0.4242640687119286*alphaL[14]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.5303300858899102*alphaL[16]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaL[16])+0.4242640687119286*alphaL[14]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[14]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.6363961030678927*alphaL[5]-0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210471*(alphaL[9]+alphaL[8]))-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[14])+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.6363961030678927*alphaL[5]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]-0.4242640687119286*alphaL[14]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[16])-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]+0.4242640687119286*alphaL[14]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaL[16])-0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[16]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaL[16])+0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[14]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaL[0]-0.3952847075210471*(alphaL[9]+alphaL[8]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[14])+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]-0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[16])-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]+0.4242640687119286*alphaL[14]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]+0.4743416490252568*(alphaL[3]+alphaL[2])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaL[16])-0.4242640687119286*alphaL[14]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*alphaL[6]-0.6363961030678927*(alphaL[5]+alphaL[4])-0.4743416490252568*(alphaL[3]+alphaL[2])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[16]-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaL[16])+0.4242640687119286*alphaL[14]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*alphaL[6]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[14]+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]-0.6363961030678927*alphaL[5]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210471*(alphaL[9]+alphaL[8]))+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[14])+0.3162277660168378*alphaL[9]-0.3952847075210471*alphaL[8]+0.6363961030678927*alphaL[5]+0.4743416490252568*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]-0.4242640687119286*alphaL[14]-0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])-0.6363961030678927*(alphaL[6]+alphaL[5])+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[16])-0.3952847075210471*alphaL[9]+0.3162277660168378*alphaL[8]+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaL[16]+0.4242640687119286*alphaL[14]+0.8538149682454614*alphaL[10]+0.3162277660168378*(alphaL[9]+alphaL[8])+0.6363961030678927*(alphaL[6]+alphaL[5]+alphaL[4])+0.4743416490252568*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
