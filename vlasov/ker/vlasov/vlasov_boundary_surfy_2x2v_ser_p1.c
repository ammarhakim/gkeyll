#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 

  const double *Fhat_nodal = &flux[8]; 
  double G1[8] = {0.0}; 
  double Ghat[8] = {0.0}; 
  G1[0] = 0.5*Fhat_nodal[3]+0.5*Fhat_nodal[2]+0.5*Fhat_nodal[1]+0.5*Fhat_nodal[0]; 
  G1[1] = 0.5*Fhat_nodal[3]+0.5*Fhat_nodal[2]-0.5*Fhat_nodal[1]-0.5*Fhat_nodal[0]; 
  G1[2] = 0.5*Fhat_nodal[3]-0.5*Fhat_nodal[2]+0.5*Fhat_nodal[1]-0.5*Fhat_nodal[0]; 
  G1[3] = 0.5*Fhat_nodal[3]-0.5*Fhat_nodal[2]-0.5*Fhat_nodal[1]+0.5*Fhat_nodal[0]; 
  G1[4] = 0.5*Fhat_nodal[7]+0.5*Fhat_nodal[6]+0.5*Fhat_nodal[5]+0.5*Fhat_nodal[4]; 
  G1[5] = 0.5*Fhat_nodal[7]+0.5*Fhat_nodal[6]-0.5*Fhat_nodal[5]-0.5*Fhat_nodal[4]; 
  G1[6] = 0.5*Fhat_nodal[7]-0.5*Fhat_nodal[6]+0.5*Fhat_nodal[5]-0.5*Fhat_nodal[4]; 
  G1[7] = 0.5*Fhat_nodal[7]-0.5*Fhat_nodal[6]-0.5*Fhat_nodal[5]+0.5*Fhat_nodal[4]; 
  Ghat[0] = 0.7071067811865475*G1[4]+0.7071067811865475*G1[0]; 
  Ghat[1] = 0.7071067811865475*G1[4]-0.7071067811865475*G1[0]; 
  Ghat[2] = 0.7071067811865475*G1[5]+0.7071067811865475*G1[1]; 
  Ghat[3] = 0.7071067811865475*G1[6]+0.7071067811865475*G1[2]; 
  Ghat[4] = 0.7071067811865475*G1[5]-0.7071067811865475*G1[1]; 
  Ghat[5] = 0.7071067811865475*G1[6]-0.7071067811865475*G1[2]; 
  Ghat[6] = 0.7071067811865475*G1[7]+0.7071067811865475*G1[3]; 
  Ghat[7] = 0.7071067811865475*G1[7]-0.7071067811865475*G1[3]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*Ghat[0]*dx11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dx11); 
  out[2] += -(1.224744871391589*Ghat[0]*dx11); 
  out[3] += -(0.7071067811865475*Ghat[2]*dx11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dx11); 
  out[5] += -(1.224744871391589*Ghat[1]*dx11); 
  out[6] += -(0.7071067811865475*Ghat[4]*dx11); 
  out[7] += -(1.224744871391589*Ghat[2]*dx11); 
  out[8] += -(0.7071067811865475*Ghat[5]*dx11); 
  out[9] += -(1.224744871391589*Ghat[3]*dx11); 
  out[10] += -(0.7071067811865475*Ghat[6]*dx11); 
  out[11] += -(1.224744871391589*Ghat[4]*dx11); 
  out[12] += -(1.224744871391589*Ghat[5]*dx11); 
  out[13] += -(0.7071067811865475*Ghat[7]*dx11); 
  out[14] += -(1.224744871391589*Ghat[6]*dx11); 
  out[15] += -(1.224744871391589*Ghat[7]*dx11); 

  } else { 

  out[0] += 0.7071067811865475*Ghat[0]*dx11; 
  out[1] += 0.7071067811865475*Ghat[1]*dx11; 
  out[2] += -(1.224744871391589*Ghat[0]*dx11); 
  out[3] += 0.7071067811865475*Ghat[2]*dx11; 
  out[4] += 0.7071067811865475*Ghat[3]*dx11; 
  out[5] += -(1.224744871391589*Ghat[1]*dx11); 
  out[6] += 0.7071067811865475*Ghat[4]*dx11; 
  out[7] += -(1.224744871391589*Ghat[2]*dx11); 
  out[8] += 0.7071067811865475*Ghat[5]*dx11; 
  out[9] += -(1.224744871391589*Ghat[3]*dx11); 
  out[10] += 0.7071067811865475*Ghat[6]*dx11; 
  out[11] += -(1.224744871391589*Ghat[4]*dx11); 
  out[12] += -(1.224744871391589*Ghat[5]*dx11); 
  out[13] += 0.7071067811865475*Ghat[7]*dx11; 
  out[14] += -(1.224744871391589*Ghat[6]*dx11); 
  out[15] += -(1.224744871391589*Ghat[7]*dx11); 

  } 
  return 0.0;

} 
