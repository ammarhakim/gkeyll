#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  double wv = w[1]; 

  double fUpwind_r[3] = {0.0}; 
  double fUpwind_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 
  double Ghat_l[3] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = 1.58113883008419*fc[4]+1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_r[1] = 1.58113883008419*fc[6]+1.224744871391589*fc[3]+0.7071067811865475*fc[2]; 
  fUpwind_r[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[5]; 

  fUpwind_l[0] = 1.58113883008419*fl[4]+1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fUpwind_l[1] = 1.58113883008419*fl[6]+1.224744871391589*fl[3]+0.7071067811865475*fl[2]; 
  fUpwind_l[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[5]; 

  } else { 

  fUpwind_r[0] = 1.58113883008419*fr[4]-1.224744871391589*fr[1]+0.7071067811865475*fr[0]; 
  fUpwind_r[1] = 1.58113883008419*fr[6]-1.224744871391589*fr[3]+0.7071067811865475*fr[2]; 
  fUpwind_r[2] = 0.7071067811865475*fr[5]-1.224744871391589*fr[7]; 

  fUpwind_l[0] = 1.58113883008419*fc[4]-1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_l[1] = 1.58113883008419*fc[6]-1.224744871391589*fc[3]+0.7071067811865475*fc[2]; 
  fUpwind_l[2] = 0.7071067811865475*fc[5]-1.224744871391589*fc[7]; 

  } 
  Ghat_l[0] = 2.738612787525831*fUpwind_l[1]*hamil[2]+1.224744871391589*fUpwind_l[0]*hamil[1]; 
  Ghat_l[1] = (2.449489742783178*fUpwind_l[2]+2.738612787525831*fUpwind_l[0])*hamil[2]+1.224744871391589*fUpwind_l[1]*hamil[1]; 
  Ghat_l[2] = 2.449489742783178*fUpwind_l[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind_l[2]; 

  Ghat_r[0] = 2.738612787525831*fUpwind_r[1]*hamil[2]+1.224744871391589*fUpwind_r[0]*hamil[1]; 
  Ghat_r[1] = (2.449489742783178*fUpwind_r[2]+2.738612787525831*fUpwind_r[0])*hamil[2]+1.224744871391589*fUpwind_r[1]*hamil[1]; 
  Ghat_r[2] = 2.449489742783178*fUpwind_r[1]*hamil[2]+1.224744871391589*hamil[1]*fUpwind_r[2]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10; 

  return fabs(1.767766952966369*dv10*dx10*(1.732050807568877*hamil[1]));

} 
