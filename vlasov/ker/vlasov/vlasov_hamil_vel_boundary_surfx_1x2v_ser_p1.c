#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double Ghat[4]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 1.060660171779821*hamil[3]*fskin[5]+0.6123724356957944*fskin[3]*hamil[3]+1.060660171779821*fskin[1]*hamil[1]+0.6123724356957944*fskin[0]*hamil[1]; 
  Ghat[1] = 1.060660171779821*hamil[3]*fskin[7]+0.6123724356957944*hamil[3]*fskin[6]+1.060660171779821*hamil[1]*fskin[4]+0.6123724356957944*hamil[1]*fskin[2]; 
  Ghat[2] = 1.060660171779821*hamil[1]*fskin[5]+1.060660171779821*fskin[1]*hamil[3]+0.6123724356957944*fskin[0]*hamil[3]+0.6123724356957944*hamil[1]*fskin[3]; 
  Ghat[3] = 1.060660171779821*hamil[1]*fskin[7]+0.6123724356957944*hamil[1]*fskin[6]+1.060660171779821*hamil[3]*fskin[4]+0.6123724356957944*fskin[2]*hamil[3]; 

  } else { 

  Ghat[0] = -(1.060660171779821*hamil[3]*fedge[5])+0.6123724356957944*fedge[3]*hamil[3]-1.060660171779821*fedge[1]*hamil[1]+0.6123724356957944*fedge[0]*hamil[1]; 
  Ghat[1] = -(1.060660171779821*hamil[3]*fedge[7])+0.6123724356957944*hamil[3]*fedge[6]-1.060660171779821*hamil[1]*fedge[4]+0.6123724356957944*hamil[1]*fedge[2]; 
  Ghat[2] = -(1.060660171779821*hamil[1]*fedge[5])-1.060660171779821*fedge[1]*hamil[3]+0.6123724356957944*fedge[0]*hamil[3]+0.6123724356957944*hamil[1]*fedge[3]; 
  Ghat[3] = -(1.060660171779821*hamil[1]*fedge[7])+0.6123724356957944*hamil[1]*fedge[6]-1.060660171779821*hamil[3]*fedge[4]+0.6123724356957944*fedge[2]*hamil[3]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[6] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[3]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 1.060660171779821*hamil[3]*fedge[5]+0.6123724356957944*fedge[3]*hamil[3]+1.060660171779821*fedge[1]*hamil[1]+0.6123724356957944*fedge[0]*hamil[1]; 
  Ghat[1] = 1.060660171779821*hamil[3]*fedge[7]+0.6123724356957944*hamil[3]*fedge[6]+1.060660171779821*hamil[1]*fedge[4]+0.6123724356957944*hamil[1]*fedge[2]; 
  Ghat[2] = 1.060660171779821*hamil[1]*fedge[5]+1.060660171779821*fedge[1]*hamil[3]+0.6123724356957944*fedge[0]*hamil[3]+0.6123724356957944*hamil[1]*fedge[3]; 
  Ghat[3] = 1.060660171779821*hamil[1]*fedge[7]+0.6123724356957944*hamil[1]*fedge[6]+1.060660171779821*hamil[3]*fedge[4]+0.6123724356957944*fedge[2]*hamil[3]; 

  } else { 

  Ghat[0] = -(1.060660171779821*hamil[3]*fskin[5])+0.6123724356957944*fskin[3]*hamil[3]-1.060660171779821*fskin[1]*hamil[1]+0.6123724356957944*fskin[0]*hamil[1]; 
  Ghat[1] = -(1.060660171779821*hamil[3]*fskin[7])+0.6123724356957944*hamil[3]*fskin[6]-1.060660171779821*hamil[1]*fskin[4]+0.6123724356957944*hamil[1]*fskin[2]; 
  Ghat[2] = -(1.060660171779821*hamil[1]*fskin[5])-1.060660171779821*fskin[1]*hamil[3]+0.6123724356957944*fskin[0]*hamil[3]+0.6123724356957944*hamil[1]*fskin[3]; 
  Ghat[3] = -(1.060660171779821*hamil[1]*fskin[7])+0.6123724356957944*hamil[1]*fskin[6]-1.060660171779821*hamil[3]*fskin[4]+0.6123724356957944*fskin[2]*hamil[3]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[6] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[7] += -(1.224744871391589*Ghat[3]*dv10*dx10); 

  } 
  return 0.0;

} 
