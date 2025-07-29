#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  const double dv11 = 2.0/dxv[3]; 

  const double *Ghat = &flux[8]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*flux[0]*dv11); 
  out[1] += -(0.7071067811865475*flux[1]*dv11); 
  out[2] += -(0.7071067811865475*flux[2]*dv11); 
  out[3] += -(0.7071067811865475*flux[3]*dv11); 
  out[4] += -(1.224744871391589*flux[0]*dv11); 
  out[5] += -(0.7071067811865475*flux[4]*dv11); 
  out[6] += -(0.7071067811865475*flux[5]*dv11); 
  out[7] += -(0.7071067811865475*flux[6]*dv11); 
  out[8] += -(1.224744871391589*flux[1]*dv11); 
  out[9] += -(1.224744871391589*flux[2]*dv11); 
  out[10] += -(1.224744871391589*flux[3]*dv11); 
  out[11] += -(0.7071067811865475*flux[7]*dv11); 
  out[12] += -(1.224744871391589*flux[4]*dv11); 
  out[13] += -(1.224744871391589*flux[5]*dv11); 
  out[14] += -(1.224744871391589*flux[6]*dv11); 
  out[15] += -(1.224744871391589*flux[7]*dv11); 

  } else { 

  out[0] += 0.7071067811865475*flux[0]*dv11; 
  out[1] += 0.7071067811865475*flux[1]*dv11; 
  out[2] += 0.7071067811865475*flux[2]*dv11; 
  out[3] += 0.7071067811865475*flux[3]*dv11; 
  out[4] += -(1.224744871391589*flux[0]*dv11); 
  out[5] += 0.7071067811865475*flux[4]*dv11; 
  out[6] += 0.7071067811865475*flux[5]*dv11; 
  out[7] += 0.7071067811865475*flux[6]*dv11; 
  out[8] += -(1.224744871391589*flux[1]*dv11); 
  out[9] += -(1.224744871391589*flux[2]*dv11); 
  out[10] += -(1.224744871391589*flux[3]*dv11); 
  out[11] += 0.7071067811865475*flux[7]*dv11; 
  out[12] += -(1.224744871391589*flux[4]*dv11); 
  out[13] += -(1.224744871391589*flux[5]*dv11); 
  out[14] += -(1.224744871391589*flux[6]*dv11); 
  out[15] += -(1.224744871391589*flux[7]*dv11); 

  } 
  return 0.0;

} 
