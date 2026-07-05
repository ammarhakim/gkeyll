#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x3v_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  double G1[8] = {0.0}; 
  double Ghat[8] = {0.0}; 
  G1[0] = 0.3535533905932737*Fhat_l_nodal_c[7]+0.3535533905932737*Fhat_l_nodal_c[6]+0.3535533905932737*Fhat_l_nodal_c[5]+0.3535533905932737*Fhat_l_nodal_c[4]+0.3535533905932737*Fhat_l_nodal_c[3]+0.3535533905932737*Fhat_l_nodal_c[2]+0.3535533905932737*Fhat_l_nodal_c[1]+0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[1] = 0.3535533905932737*Fhat_l_nodal_c[7]+0.3535533905932737*Fhat_l_nodal_c[6]+0.3535533905932737*Fhat_l_nodal_c[5]+0.3535533905932737*Fhat_l_nodal_c[4]-0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_l_nodal_c[2]-0.3535533905932737*Fhat_l_nodal_c[1]-0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[2] = 0.3535533905932737*Fhat_l_nodal_c[7]+0.3535533905932737*Fhat_l_nodal_c[6]-0.3535533905932737*Fhat_l_nodal_c[5]-0.3535533905932737*Fhat_l_nodal_c[4]+0.3535533905932737*Fhat_l_nodal_c[3]+0.3535533905932737*Fhat_l_nodal_c[2]-0.3535533905932737*Fhat_l_nodal_c[1]-0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[3] = 0.3535533905932737*Fhat_l_nodal_c[7]-0.3535533905932737*Fhat_l_nodal_c[6]+0.3535533905932737*Fhat_l_nodal_c[5]-0.3535533905932737*Fhat_l_nodal_c[4]+0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_l_nodal_c[2]+0.3535533905932737*Fhat_l_nodal_c[1]-0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[4] = 0.3535533905932737*Fhat_l_nodal_c[7]+0.3535533905932737*Fhat_l_nodal_c[6]-0.3535533905932737*Fhat_l_nodal_c[5]-0.3535533905932737*Fhat_l_nodal_c[4]-0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_l_nodal_c[2]+0.3535533905932737*Fhat_l_nodal_c[1]+0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[5] = 0.3535533905932737*Fhat_l_nodal_c[7]-0.3535533905932737*Fhat_l_nodal_c[6]+0.3535533905932737*Fhat_l_nodal_c[5]-0.3535533905932737*Fhat_l_nodal_c[4]-0.3535533905932737*Fhat_l_nodal_c[3]+0.3535533905932737*Fhat_l_nodal_c[2]-0.3535533905932737*Fhat_l_nodal_c[1]+0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[6] = 0.3535533905932737*Fhat_l_nodal_c[7]-0.3535533905932737*Fhat_l_nodal_c[6]-0.3535533905932737*Fhat_l_nodal_c[5]+0.3535533905932737*Fhat_l_nodal_c[4]+0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_l_nodal_c[2]-0.3535533905932737*Fhat_l_nodal_c[1]+0.3535533905932737*Fhat_l_nodal_c[0]; 
  G1[7] = 0.3535533905932737*Fhat_l_nodal_c[7]-0.3535533905932737*Fhat_l_nodal_c[6]-0.3535533905932737*Fhat_l_nodal_c[5]+0.3535533905932737*Fhat_l_nodal_c[4]-0.3535533905932737*Fhat_l_nodal_c[3]+0.3535533905932737*Fhat_l_nodal_c[2]+0.3535533905932737*Fhat_l_nodal_c[1]-0.3535533905932737*Fhat_l_nodal_c[0]; 
  Ghat[0] = G1[0]; 
  Ghat[1] = G1[1]; 
  Ghat[2] = G1[2]; 
  Ghat[3] = G1[3]; 
  Ghat[4] = G1[4]; 
  Ghat[5] = G1[5]; 
  Ghat[6] = G1[6]; 
  Ghat[7] = G1[7]; 
  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -(1.224744871391589*Ghat[1]*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dx10); 
  out[7] += 0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -(1.224744871391589*Ghat[3]*dx10); 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -(1.224744871391589*Ghat[4]*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dx10); 
  out[14] += 0.7071067811865475*Ghat[7]*dx10; 
  out[15] += -(1.224744871391589*Ghat[7]*dx10); 

  G1[0] = 0.3535533905932737*Fhat_r_nodal_c[7]+0.3535533905932737*Fhat_r_nodal_c[6]+0.3535533905932737*Fhat_r_nodal_c[5]+0.3535533905932737*Fhat_r_nodal_c[4]+0.3535533905932737*Fhat_r_nodal_c[3]+0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*Fhat_r_nodal_c[1]+0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[1] = 0.3535533905932737*Fhat_r_nodal_c[7]+0.3535533905932737*Fhat_r_nodal_c[6]+0.3535533905932737*Fhat_r_nodal_c[5]+0.3535533905932737*Fhat_r_nodal_c[4]-0.3535533905932737*Fhat_r_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]-0.3535533905932737*Fhat_r_nodal_c[1]-0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[2] = 0.3535533905932737*Fhat_r_nodal_c[7]+0.3535533905932737*Fhat_r_nodal_c[6]-0.3535533905932737*Fhat_r_nodal_c[5]-0.3535533905932737*Fhat_r_nodal_c[4]+0.3535533905932737*Fhat_r_nodal_c[3]+0.3535533905932737*Fhat_r_nodal_c[2]-0.3535533905932737*Fhat_r_nodal_c[1]-0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[3] = 0.3535533905932737*Fhat_r_nodal_c[7]-0.3535533905932737*Fhat_r_nodal_c[6]+0.3535533905932737*Fhat_r_nodal_c[5]-0.3535533905932737*Fhat_r_nodal_c[4]+0.3535533905932737*Fhat_r_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*Fhat_r_nodal_c[1]-0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[4] = 0.3535533905932737*Fhat_r_nodal_c[7]+0.3535533905932737*Fhat_r_nodal_c[6]-0.3535533905932737*Fhat_r_nodal_c[5]-0.3535533905932737*Fhat_r_nodal_c[4]-0.3535533905932737*Fhat_r_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*Fhat_r_nodal_c[1]+0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[5] = 0.3535533905932737*Fhat_r_nodal_c[7]-0.3535533905932737*Fhat_r_nodal_c[6]+0.3535533905932737*Fhat_r_nodal_c[5]-0.3535533905932737*Fhat_r_nodal_c[4]-0.3535533905932737*Fhat_r_nodal_c[3]+0.3535533905932737*Fhat_r_nodal_c[2]-0.3535533905932737*Fhat_r_nodal_c[1]+0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[6] = 0.3535533905932737*Fhat_r_nodal_c[7]-0.3535533905932737*Fhat_r_nodal_c[6]-0.3535533905932737*Fhat_r_nodal_c[5]+0.3535533905932737*Fhat_r_nodal_c[4]+0.3535533905932737*Fhat_r_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]-0.3535533905932737*Fhat_r_nodal_c[1]+0.3535533905932737*Fhat_r_nodal_c[0]; 
  G1[7] = 0.3535533905932737*Fhat_r_nodal_c[7]-0.3535533905932737*Fhat_r_nodal_c[6]-0.3535533905932737*Fhat_r_nodal_c[5]+0.3535533905932737*Fhat_r_nodal_c[4]-0.3535533905932737*Fhat_r_nodal_c[3]+0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*Fhat_r_nodal_c[1]-0.3535533905932737*Fhat_r_nodal_c[0]; 
  Ghat[0] = G1[0]; 
  Ghat[1] = G1[1]; 
  Ghat[2] = G1[2]; 
  Ghat[3] = G1[3]; 
  Ghat[4] = G1[4]; 
  Ghat[5] = G1[5]; 
  Ghat[6] = G1[6]; 
  Ghat[7] = G1[7]; 
  out[0] += -(0.7071067811865475*Ghat[0]*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dx10); 
  out[7] += -(0.7071067811865475*Ghat[4]*dx10); 
  out[8] += -(1.224744871391589*Ghat[3]*dx10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dx10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dx10); 
  out[11] += -(1.224744871391589*Ghat[4]*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dx10); 
  out[14] += -(0.7071067811865475*Ghat[7]*dx10); 
  out[15] += -(1.224744871391589*Ghat[7]*dx10); 

  return 0.0;

} 
