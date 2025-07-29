#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_2x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double Ghat[4]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 1.5*fskin[1]*hamil[1]+0.8660254037844386*fskin[0]*hamil[1]; 
  Ghat[1] = 1.5*hamil[1]*fskin[4]+0.8660254037844386*hamil[1]*fskin[2]; 
  Ghat[2] = 1.5*hamil[1]*fskin[5]+0.8660254037844386*hamil[1]*fskin[3]; 
  Ghat[3] = 1.5*hamil[1]*fskin[7]+0.8660254037844386*hamil[1]*fskin[6]; 

  } else { 

  Ghat[0] = 0.8660254037844386*fedge[0]*hamil[1]-1.5*fedge[1]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[1]*fedge[2]-1.5*hamil[1]*fedge[4]; 
  Ghat[2] = 0.8660254037844386*hamil[1]*fedge[3]-1.5*hamil[1]*fedge[5]; 
  Ghat[3] = 0.8660254037844386*hamil[1]*fedge[6]-1.5*hamil[1]*fedge[7]; 

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

  Ghat[0] = 1.5*fedge[1]*hamil[1]+0.8660254037844386*fedge[0]*hamil[1]; 
  Ghat[1] = 1.5*hamil[1]*fedge[4]+0.8660254037844386*hamil[1]*fedge[2]; 
  Ghat[2] = 1.5*hamil[1]*fedge[5]+0.8660254037844386*hamil[1]*fedge[3]; 
  Ghat[3] = 1.5*hamil[1]*fedge[7]+0.8660254037844386*hamil[1]*fedge[6]; 

  } else { 

  Ghat[0] = 0.8660254037844386*fskin[0]*hamil[1]-1.5*fskin[1]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[1]*fskin[2]-1.5*hamil[1]*fskin[4]; 
  Ghat[2] = 0.8660254037844386*hamil[1]*fskin[3]-1.5*hamil[1]*fskin[5]; 
  Ghat[3] = 0.8660254037844386*hamil[1]*fskin[6]-1.5*hamil[1]*fskin[7]; 

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
