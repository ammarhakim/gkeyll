#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  const double dv10 = 2.0/dxv[1]; 

  const double *Ghat = &flux[0]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 

  } else { 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 

  } 
  return 0.0;

} 
