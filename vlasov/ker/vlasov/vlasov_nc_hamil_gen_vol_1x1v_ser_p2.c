#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wx1 = w[1]; 
  const double dv10_sq = dv10*dv10; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  out[1] += (3.354101966249684*hamil[2]*poisson_tensor_conf_0[2]*f[6]+1.5*hamil[1]*poisson_tensor_conf_0[2]*f[4]+3.354101966249685*poisson_tensor_conf_0[1]*hamil[2]*f[3]+3.354101966249685*poisson_tensor_conf_0[0]*f[2]*hamil[2]+1.5*f[1]*hamil[1]*poisson_tensor_conf_0[1]+1.5*f[0]*poisson_tensor_conf_0[0]*hamil[1])*dv10*dx10; 
  out[3] += (3.0*poisson_tensor_conf_0[1]*hamil[2]*f[7]+1.5*hamil[1]*poisson_tensor_conf_0[2]*f[6]+3.0*poisson_tensor_conf_0[0]*hamil[2]*f[5]+3.354101966249685*hamil[2]*poisson_tensor_conf_0[2]*f[4]+1.5*hamil[1]*poisson_tensor_conf_0[1]*f[3]+3.354101966249685*f[1]*poisson_tensor_conf_0[1]*hamil[2]+3.354101966249685*f[0]*poisson_tensor_conf_0[0]*hamil[2]+1.5*poisson_tensor_conf_0[0]*hamil[1]*f[2])*dv10*dx10; 
  out[4] += (6.708203932499369*poisson_tensor_conf_0[1]*hamil[2]*f[6]+3.0*hamil[1]*poisson_tensor_conf_0[1]*f[4]+6.708203932499369*hamil[2]*poisson_tensor_conf_0[2]*f[3]+7.5*poisson_tensor_conf_0[0]*hamil[2]*f[3]+3.0*f[1]*hamil[1]*poisson_tensor_conf_0[2]+7.5*poisson_tensor_conf_0[1]*f[2]*hamil[2]+3.354101966249685*f[0]*hamil[1]*poisson_tensor_conf_0[1]+3.354101966249685*poisson_tensor_conf_0[0]*f[1]*hamil[1])*dv10*dx10; 
  out[6] += (6.0*hamil[2]*poisson_tensor_conf_0[2]*f[7]+6.708203932499369*poisson_tensor_conf_0[0]*hamil[2]*f[7]+3.0*hamil[1]*poisson_tensor_conf_0[1]*f[6]+6.708203932499369*poisson_tensor_conf_0[1]*hamil[2]*f[5]+6.708203932499369*poisson_tensor_conf_0[1]*hamil[2]*f[4]+3.0*hamil[1]*poisson_tensor_conf_0[2]*f[3]+3.354101966249684*poisson_tensor_conf_0[0]*hamil[1]*f[3]+6.708203932499369*f[1]*hamil[2]*poisson_tensor_conf_0[2]+7.500000000000001*f[0]*poisson_tensor_conf_0[1]*hamil[2]+7.500000000000001*poisson_tensor_conf_0[0]*f[1]*hamil[2]+3.354101966249684*hamil[1]*poisson_tensor_conf_0[1]*f[2])*dv10*dx10; 
  out[7] += (1.5*hamil[1]*poisson_tensor_conf_0[1]*f[7]+3.0*hamil[2]*poisson_tensor_conf_0[2]*f[6]+1.5*poisson_tensor_conf_0[0]*hamil[1]*f[5]+3.0*poisson_tensor_conf_0[1]*hamil[2]*f[3]+3.0*poisson_tensor_conf_0[0]*f[2]*hamil[2])*dv10*dx10; 

} 
