#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  const double dv11 = 2.0/dxv[2]; 

  const double *Ghat = &flux[4]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*flux[0]*dv11); 
  out[1] += -(0.7071067811865475*flux[1]*dv11); 
  out[2] += -(0.7071067811865475*flux[2]*dv11); 
  out[3] += -(1.224744871391589*flux[0]*dv11); 
  out[4] += -(0.7071067811865475*flux[3]*dv11); 
  out[5] += -(1.224744871391589*flux[1]*dv11); 
  out[6] += -(1.224744871391589*flux[2]*dv11); 
  out[7] += -(1.224744871391589*flux[3]*dv11); 

  } else { 

  out[0] += 0.7071067811865475*flux[0]*dv11; 
  out[1] += 0.7071067811865475*flux[1]*dv11; 
  out[2] += 0.7071067811865475*flux[2]*dv11; 
  out[3] += -(1.224744871391589*flux[0]*dv11); 
  out[4] += 0.7071067811865475*flux[3]*dv11; 
  out[5] += -(1.224744871391589*flux[1]*dv11); 
  out[6] += -(1.224744871391589*flux[2]*dv11); 
  out[7] += -(1.224744871391589*flux[3]*dv11); 

  } 
  return 0.0;

} 
