#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double fUpwind_r[8] = {0.0}; 
  double fUpwind_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 
  double Ghat_l[8] = {0.0}; 
  if (wv>0) { 

  fUpwind_r[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_r[1] = 1.224744871391589*fc[5]+0.7071067811865475*fc[2]; 
  fUpwind_r[2] = 1.224744871391589*fc[6]+0.7071067811865475*fc[3]; 
  fUpwind_r[3] = 1.224744871391589*fc[8]+0.7071067811865475*fc[4]; 
  fUpwind_r[4] = 1.224744871391589*fc[11]+0.7071067811865475*fc[7]; 
  fUpwind_r[5] = 1.224744871391589*fc[12]+0.7071067811865475*fc[9]; 
  fUpwind_r[6] = 1.224744871391589*fc[13]+0.7071067811865475*fc[10]; 
  fUpwind_r[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[14]; 

  fUpwind_l[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fUpwind_l[1] = 1.224744871391589*fl[5]+0.7071067811865475*fl[2]; 
  fUpwind_l[2] = 1.224744871391589*fl[6]+0.7071067811865475*fl[3]; 
  fUpwind_l[3] = 1.224744871391589*fl[8]+0.7071067811865475*fl[4]; 
  fUpwind_l[4] = 1.224744871391589*fl[11]+0.7071067811865475*fl[7]; 
  fUpwind_l[5] = 1.224744871391589*fl[12]+0.7071067811865475*fl[9]; 
  fUpwind_l[6] = 1.224744871391589*fl[13]+0.7071067811865475*fl[10]; 
  fUpwind_l[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[14]; 

  } else { 

  fUpwind_r[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  fUpwind_r[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[5]; 
  fUpwind_r[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[6]; 
  fUpwind_r[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[8]; 
  fUpwind_r[4] = 0.7071067811865475*fr[7]-1.224744871391589*fr[11]; 
  fUpwind_r[5] = 0.7071067811865475*fr[9]-1.224744871391589*fr[12]; 
  fUpwind_r[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[13]; 
  fUpwind_r[7] = 0.7071067811865475*fr[14]-1.224744871391589*fr[15]; 

  fUpwind_l[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  fUpwind_l[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[5]; 
  fUpwind_l[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[6]; 
  fUpwind_l[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[8]; 
  fUpwind_l[4] = 0.7071067811865475*fc[7]-1.224744871391589*fc[11]; 
  fUpwind_l[5] = 0.7071067811865475*fc[9]-1.224744871391589*fc[12]; 
  fUpwind_l[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[13]; 
  fUpwind_l[7] = 0.7071067811865475*fc[14]-1.224744871391589*fc[15]; 

  } 
  Ghat_l[0] = 0.8660254037844386*(fUpwind_l[3]*hamil[3]+fUpwind_l[0]*hamil[1]); 
  Ghat_l[1] = 0.8660254037844386*(hamil[3]*fUpwind_l[5]+fUpwind_l[1]*hamil[1]); 
  Ghat_l[2] = 0.8660254037844386*(hamil[3]*fUpwind_l[6]+hamil[1]*fUpwind_l[2]); 
  Ghat_l[3] = 0.8660254037844386*(fUpwind_l[0]*hamil[3]+hamil[1]*fUpwind_l[3]); 
  Ghat_l[4] = 0.8660254037844386*(hamil[3]*fUpwind_l[7]+hamil[1]*fUpwind_l[4]); 
  Ghat_l[5] = 0.8660254037844386*(hamil[1]*fUpwind_l[5]+fUpwind_l[1]*hamil[3]); 
  Ghat_l[6] = 0.8660254037844386*(hamil[1]*fUpwind_l[6]+fUpwind_l[2]*hamil[3]); 
  Ghat_l[7] = 0.8660254037844386*(hamil[1]*fUpwind_l[7]+hamil[3]*fUpwind_l[4]); 

  Ghat_r[0] = 0.8660254037844386*(fUpwind_r[3]*hamil[3]+fUpwind_r[0]*hamil[1]); 
  Ghat_r[1] = 0.8660254037844386*(hamil[3]*fUpwind_r[5]+fUpwind_r[1]*hamil[1]); 
  Ghat_r[2] = 0.8660254037844386*(hamil[3]*fUpwind_r[6]+hamil[1]*fUpwind_r[2]); 
  Ghat_r[3] = 0.8660254037844386*(fUpwind_r[0]*hamil[3]+hamil[1]*fUpwind_r[3]); 
  Ghat_r[4] = 0.8660254037844386*(hamil[3]*fUpwind_r[7]+hamil[1]*fUpwind_r[4]); 
  Ghat_r[5] = 0.8660254037844386*(hamil[1]*fUpwind_r[5]+fUpwind_r[1]*hamil[3]); 
  Ghat_r[6] = 0.8660254037844386*(hamil[1]*fUpwind_r[6]+fUpwind_r[2]*hamil[3]); 
  Ghat_r[7] = 0.8660254037844386*(hamil[1]*fUpwind_r[7]+hamil[3]*fUpwind_r[4]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10*dx10; 

  return fabs(0.75*dv10*dx10*(1.732050807568877*hamil[1]));

} 
