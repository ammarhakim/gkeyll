#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 

  double Ghat[2]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 1.5*fskin[1]*hamil[1]+0.8660254037844386*fskin[0]*hamil[1]; 
  Ghat[1] = 1.5*hamil[1]*fskin[3]+0.8660254037844386*hamil[1]*fskin[2]; 

  } else { 

  Ghat[0] = 0.8660254037844386*fedge[0]*hamil[1]-1.5*fedge[1]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[1]*fedge[2]-1.5*hamil[1]*fedge[3]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 1.5*fedge[1]*hamil[1]+0.8660254037844386*fedge[0]*hamil[1]; 
  Ghat[1] = 1.5*hamil[1]*fedge[3]+0.8660254037844386*hamil[1]*fedge[2]; 

  } else { 

  Ghat[0] = 0.8660254037844386*fskin[0]*hamil[1]-1.5*fskin[1]*hamil[1]; 
  Ghat[1] = 0.8660254037844386*hamil[1]*fskin[2]-1.5*hamil[1]*fskin[3]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } 
  return 0.0;

} 
