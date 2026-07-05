#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[2]; 

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

  out[0] += -(0.7071067811865475*Ghat[0]*dv11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv11); 
  out[3] += -(1.224744871391589*Ghat[0]*dv11); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv11); 
  out[6] += -(1.224744871391589*Ghat[1]*dv11); 
  out[7] += -(1.224744871391589*Ghat[2]*dv11); 
  out[8] += -(0.7071067811865475*Ghat[5]*dv11); 
  out[9] += -(0.7071067811865475*Ghat[6]*dv11); 
  out[10] += -(1.224744871391589*Ghat[3]*dv11); 
  out[11] += -(1.224744871391589*Ghat[4]*dv11); 
  out[12] += -(0.7071067811865475*Ghat[7]*dv11); 
  out[13] += -(1.224744871391589*Ghat[5]*dv11); 
  out[14] += -(1.224744871391589*Ghat[6]*dv11); 
  out[15] += -(1.224744871391589*Ghat[7]*dv11); 

  } else { 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -(1.224744871391589*Ghat[0]*dv11); 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -(1.224744871391589*Ghat[1]*dv11); 
  out[7] += -(1.224744871391589*Ghat[2]*dv11); 
  out[8] += 0.7071067811865475*Ghat[5]*dv11; 
  out[9] += 0.7071067811865475*Ghat[6]*dv11; 
  out[10] += -(1.224744871391589*Ghat[3]*dv11); 
  out[11] += -(1.224744871391589*Ghat[4]*dv11); 
  out[12] += 0.7071067811865475*Ghat[7]*dv11; 
  out[13] += -(1.224744871391589*Ghat[5]*dv11); 
  out[14] += -(1.224744871391589*Ghat[6]*dv11); 
  out[15] += -(1.224744871391589*Ghat[7]*dv11); 

  } 
  return 0.0;

} 
