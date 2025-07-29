#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_surfx_2x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double wv = w[2]; 

  double Ghat_r[4]; 
  double Ghat_l[4]; 
  if (wv>0) { 

  Ghat_r[0] = (1.5*fc[1]+0.8660254037844386*fc[0])*hamil[1]; 
  Ghat_r[1] = hamil[1]*(1.5*fc[4]+0.8660254037844386*fc[2]); 
  Ghat_r[2] = hamil[1]*(1.5*fc[5]+0.8660254037844386*fc[3]); 
  Ghat_r[3] = hamil[1]*(1.5*fc[7]+0.8660254037844386*fc[6]); 

  Ghat_l[0] = (1.5*fl[1]+0.8660254037844386*fl[0])*hamil[1]; 
  Ghat_l[1] = hamil[1]*(1.5*fl[4]+0.8660254037844386*fl[2]); 
  Ghat_l[2] = hamil[1]*(1.5*fl[5]+0.8660254037844386*fl[3]); 
  Ghat_l[3] = hamil[1]*(1.5*fl[7]+0.8660254037844386*fl[6]); 

  } else { 

  Ghat_r[0] = 0.8660254037844386*fr[0]*hamil[1]-1.5*fr[1]*hamil[1]; 
  Ghat_r[1] = 0.8660254037844386*hamil[1]*fr[2]-1.5*hamil[1]*fr[4]; 
  Ghat_r[2] = 0.8660254037844386*hamil[1]*fr[3]-1.5*hamil[1]*fr[5]; 
  Ghat_r[3] = 0.8660254037844386*hamil[1]*fr[6]-1.5*hamil[1]*fr[7]; 

  Ghat_l[0] = 0.8660254037844386*fc[0]*hamil[1]-1.5*fc[1]*hamil[1]; 
  Ghat_l[1] = 0.8660254037844386*hamil[1]*fc[2]-1.5*hamil[1]*fc[4]; 
  Ghat_l[2] = 0.8660254037844386*hamil[1]*fc[3]-1.5*hamil[1]*fc[5]; 
  Ghat_l[3] = 0.8660254037844386*hamil[1]*fc[6]-1.5*hamil[1]*fc[7]; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10*dx10; 
  out[4] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 
  out[5] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10*dx10); 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10*dx10; 
  out[7] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10*dx10); 

  return fabs(1.0606601717798212*dv10*dx10*1.7320508075688772*hamil[1]);

} 
