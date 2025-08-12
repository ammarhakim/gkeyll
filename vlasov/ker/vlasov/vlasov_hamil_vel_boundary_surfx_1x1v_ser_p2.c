#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double fUpwind[3]; 
  double Ghat[3]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = 1.5811388300841895*fskin[4]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  fUpwind[1] = 1.5811388300841898*fskin[6]+1.224744871391589*fskin[3]+0.7071067811865475*fskin[2]; 
  fUpwind[2] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[5]; 

  } else { 

  fUpwind[0] = 1.5811388300841895*fedge[4]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  fUpwind[1] = 1.5811388300841898*fedge[6]-1.224744871391589*fedge[3]+0.7071067811865475*fedge[2]; 
  fUpwind[2] = 0.7071067811865475*fedge[5]-1.224744871391589*fedge[7]; 

  } 
  Ghat[0] = 2.7386127875258306*fUpwind[1]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 2.4494897427831783*fUpwind[2]*hamil[2]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 2.4494897427831783*fUpwind[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += -(1.5811388300841895*Ghat[0]*dv10*dx10); 
  out[5] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[6] += -(1.5811388300841898*Ghat[1]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = 1.5811388300841895*fedge[4]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  fUpwind[1] = 1.5811388300841898*fedge[6]+1.224744871391589*fedge[3]+0.7071067811865475*fedge[2]; 
  fUpwind[2] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[5]; 

  } else { 

  fUpwind[0] = 1.5811388300841895*fskin[4]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  fUpwind[1] = 1.5811388300841898*fskin[6]-1.224744871391589*fskin[3]+0.7071067811865475*fskin[2]; 
  fUpwind[2] = 0.7071067811865475*fskin[5]-1.224744871391589*fskin[7]; 

  } 
  Ghat[0] = 2.7386127875258306*fUpwind[1]*hamil[2]+1.224744871391589*fUpwind[0]*hamil[1]; 
  Ghat[1] = 2.4494897427831783*fUpwind[2]*hamil[2]+2.7386127875258306*fUpwind[0]*hamil[2]+1.224744871391589*fUpwind[1]*hamil[1]; 
  Ghat[2] = 2.4494897427831783*fUpwind[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += 1.5811388300841895*Ghat[0]*dv10*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[6] += 1.5811388300841898*Ghat[1]*dv10*dx10; 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 

  } 
  return 0.0;

} 
