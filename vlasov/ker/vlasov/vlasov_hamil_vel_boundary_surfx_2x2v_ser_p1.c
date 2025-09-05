#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double fUpwind[8]; 
  double Ghat[8]; 

  if (edge == -1) { 

  if (wv>0) { 

  fUpwind[0] = 1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  fUpwind[1] = 1.224744871391589*fskin[5]+0.7071067811865475*fskin[2]; 
  fUpwind[2] = 1.224744871391589*fskin[6]+0.7071067811865475*fskin[3]; 
  fUpwind[3] = 1.224744871391589*fskin[8]+0.7071067811865475*fskin[4]; 
  fUpwind[4] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[7]; 
  fUpwind[5] = 1.224744871391589*fskin[12]+0.7071067811865475*fskin[9]; 
  fUpwind[6] = 1.224744871391589*fskin[13]+0.7071067811865475*fskin[10]; 
  fUpwind[7] = 1.224744871391589*fskin[15]+0.7071067811865475*fskin[14]; 

  } else { 

  fUpwind[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[1]; 
  fUpwind[1] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[5]; 
  fUpwind[2] = 0.7071067811865475*fedge[3]-1.224744871391589*fedge[6]; 
  fUpwind[3] = 0.7071067811865475*fedge[4]-1.224744871391589*fedge[8]; 
  fUpwind[4] = 0.7071067811865475*fedge[7]-1.224744871391589*fedge[11]; 
  fUpwind[5] = 0.7071067811865475*fedge[9]-1.224744871391589*fedge[12]; 
  fUpwind[6] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[13]; 
  fUpwind[7] = 0.7071067811865475*fedge[14]-1.224744871391589*fedge[15]; 

  } 
  Ghat[0] = 0.8660254037844386*fUpwind[3]*hamil[3]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[3]*fUpwind[5]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 0.8660254037844386*hamil[3]*fUpwind[6]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 0.8660254037844386*fUpwind[0]*hamil[3]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 0.8660254037844386*hamil[3]*fUpwind[7]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 0.8660254037844386*hamil[1]*fUpwind[5]+0.8660254037844386*fUpwind[1]*hamil[3]; 
  Ghat[6] = 0.8660254037844386*hamil[1]*fUpwind[6]+0.8660254037844386*fUpwind[2]*hamil[3]; 
  Ghat[7] = 0.8660254037844386*hamil[1]*fUpwind[7]+0.8660254037844386*hamil[3]*fUpwind[4]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[11] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  fUpwind[0] = 1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  fUpwind[1] = 1.224744871391589*fedge[5]+0.7071067811865475*fedge[2]; 
  fUpwind[2] = 1.224744871391589*fedge[6]+0.7071067811865475*fedge[3]; 
  fUpwind[3] = 1.224744871391589*fedge[8]+0.7071067811865475*fedge[4]; 
  fUpwind[4] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[7]; 
  fUpwind[5] = 1.224744871391589*fedge[12]+0.7071067811865475*fedge[9]; 
  fUpwind[6] = 1.224744871391589*fedge[13]+0.7071067811865475*fedge[10]; 
  fUpwind[7] = 1.224744871391589*fedge[15]+0.7071067811865475*fedge[14]; 

  } else { 

  fUpwind[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[1]; 
  fUpwind[1] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[5]; 
  fUpwind[2] = 0.7071067811865475*fskin[3]-1.224744871391589*fskin[6]; 
  fUpwind[3] = 0.7071067811865475*fskin[4]-1.224744871391589*fskin[8]; 
  fUpwind[4] = 0.7071067811865475*fskin[7]-1.224744871391589*fskin[11]; 
  fUpwind[5] = 0.7071067811865475*fskin[9]-1.224744871391589*fskin[12]; 
  fUpwind[6] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[13]; 
  fUpwind[7] = 0.7071067811865475*fskin[14]-1.224744871391589*fskin[15]; 

  } 
  Ghat[0] = 0.8660254037844386*fUpwind[3]*hamil[3]+0.8660254037844386*fUpwind[0]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[3]*fUpwind[5]+0.8660254037844386*fUpwind[1]*hamil[1]; 
  Ghat[2] = 0.8660254037844386*hamil[3]*fUpwind[6]+0.8660254037844386*hamil[1]*fUpwind[2]; 
  Ghat[3] = 0.8660254037844386*fUpwind[0]*hamil[3]+0.8660254037844386*hamil[1]*fUpwind[3]; 
  Ghat[4] = 0.8660254037844386*hamil[3]*fUpwind[7]+0.8660254037844386*hamil[1]*fUpwind[4]; 
  Ghat[5] = 0.8660254037844386*hamil[1]*fUpwind[5]+0.8660254037844386*fUpwind[1]*hamil[3]; 
  Ghat[6] = 0.8660254037844386*hamil[1]*fUpwind[6]+0.8660254037844386*fUpwind[2]*hamil[3]; 
  Ghat[7] = 0.8660254037844386*hamil[1]*fUpwind[7]+0.8660254037844386*hamil[3]*fUpwind[4]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[11] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[14] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[15] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } 
  return 0.0;

} 
