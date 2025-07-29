#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 

  double Ghat_r[2]; 
  double Ghat_l[2]; 
  if (wv>0) { 

  Ghat_r[0] = (1.5*fc[1]+0.8660254037844386*fc[0])*hamil[1]; 
  Ghat_r[1] = hamil[1]*(1.5*fc[3]+0.8660254037844386*fc[2]); 

  Ghat_l[0] = (1.5*fl[1]+0.8660254037844386*fl[0])*hamil[1]; 
  Ghat_l[1] = hamil[1]*(1.5*fl[3]+0.8660254037844386*fl[2]); 

  } else { 

  Ghat_r[0] = 0.8660254037844386*fr[0]*hamil[1]-1.5*fr[1]*hamil[1]; 
  Ghat_r[1] = 0.8660254037844386*hamil[1]*fr[2]-1.5*hamil[1]*fr[3]; 

  Ghat_l[0] = 0.8660254037844386*fc[0]*hamil[1]-1.5*fc[1]*hamil[1]; 
  Ghat_l[1] = 0.8660254037844386*hamil[1]*fc[2]-1.5*hamil[1]*fc[3]; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 

  return fabs(1.0606601717798212*dv10*dx10*1.7320508075688772*hamil[1]);

} 
