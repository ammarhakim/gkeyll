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

printf(
  "[VOL DBG] dxv[0]=%.17g dxv[1]=%.17g | dx10=%.17g dv10=%.17g dv10_sq=%.17g | "
  "w[1]=%.17g | hamil[1]=%.17g hamil[2]=%.17g | "
  "f[0..7]=[%.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g] | "
  "poisson_tensor_conf_0=%p -> {%.17g, %.17g, %.17g} | "
  "factor(dv10*dx10)=%.17g | out[1]=%.17g out[3]=%.17g out[4]=%.17g out[6]=%.17g out[7]=%.17g\n",
  dxv[0], dxv[1],
  dx10, dv10, dv10_sq,
  w[1],
  hamil[1], hamil[2],
  f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7],
  (void*)poisson_tensor_conf_0,
  poisson_tensor_conf_0[0], poisson_tensor_conf_0[1], poisson_tensor_conf_0[2],
  dv10*dx10,
  out[1], out[3], out[4], out[6], out[7]
);

} 
