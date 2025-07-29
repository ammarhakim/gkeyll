#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[1] += (1.5*f[3]*hamil[3]+1.5*f[0]*hamil[1])*dv10*dx10; 
  out[4] += (1.5*hamil[3]*f[6]+1.5*hamil[1]*f[2])*dv10*dx10; 
  out[5] += (1.5*f[0]*hamil[3]+1.5*hamil[1]*f[3])*dv10*dx10; 
  out[7] += (1.5*hamil[1]*f[6]+1.5*f[2]*hamil[3])*dv10*dx10; 
} 
