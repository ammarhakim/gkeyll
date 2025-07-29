#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvy_1x2v_tensor_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  const double dv11 = 2.0/dxv[2]; 

  const double *Ghat_l = &flux_l[4]; 
  const double *Ghat_r = &flux_r[4]; 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11); 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11); 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11); 
  out[7] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11); 

  return 0.0;

} 
