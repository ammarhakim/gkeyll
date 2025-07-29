#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  out[1] += (1.5*hamil[3]*f[4]+1.5*f[0]*hamil[1])*dv10*dx10; 
  out[2] += 1.5*f[3]*hamil[3]*dv11*dx11+1.5*f[0]*hamil[2]*dv11*dx11; 
  out[5] += (1.5*hamil[3]*f[6]+1.5*f[1]*hamil[2])*dv11*dx11+(1.5*hamil[3]*f[9]+1.5*hamil[1]*f[2])*dv10*dx10; 
  out[6] += (1.5*hamil[3]*f[10]+1.5*hamil[1]*f[3])*dv10*dx10; 
  out[7] += 1.5*f[0]*hamil[3]*dv11*dx11+1.5*hamil[2]*f[3]*dv11*dx11; 
  out[8] += (1.5*hamil[1]*f[4]+1.5*f[0]*hamil[3])*dv10*dx10; 
  out[9] += 1.5*hamil[3]*f[10]*dv11*dx11+1.5*hamil[2]*f[4]*dv11*dx11; 
  out[11] += (1.5*hamil[2]*f[6]+1.5*f[1]*hamil[3])*dv11*dx11+(1.5*hamil[3]*f[14]+1.5*hamil[1]*f[7])*dv10*dx10; 
  out[12] += (1.5*hamil[3]*f[13]+1.5*hamil[2]*f[8])*dv11*dx11+(1.5*hamil[1]*f[9]+1.5*f[2]*hamil[3])*dv10*dx10; 
  out[13] += (1.5*hamil[1]*f[10]+1.5*f[3]*hamil[3])*dv10*dx10; 
  out[14] += 1.5*hamil[2]*f[10]*dv11*dx11+1.5*hamil[3]*f[4]*dv11*dx11; 
  out[15] += (1.5*hamil[2]*f[13]+1.5*hamil[3]*f[8])*dv11*dx11+(1.5*hamil[1]*f[14]+1.5*hamil[3]*f[7])*dv10*dx10; 
} 
