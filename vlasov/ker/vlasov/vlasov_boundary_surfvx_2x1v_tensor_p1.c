#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_2x1v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  const double dv10 = 2.0/dxv[2]; 

  const double *Ghat = &flux[0]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += -(1.224744871391589*Ghat[3]*dv10); 

  } else { 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -(1.224744871391589*Ghat[0]*dv10); 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -(1.224744871391589*Ghat[1]*dv10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += -(1.224744871391589*Ghat[3]*dv10); 

  } 
  return 0.0;

} 
