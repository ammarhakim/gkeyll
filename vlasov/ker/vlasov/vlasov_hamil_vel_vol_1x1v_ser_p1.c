#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[1] += 2.121320343559642*f[0]*hamil[1]*dv10*dx10; 
  out[3] += 2.121320343559642*hamil[1]*f[2]*dv10*dx10; 

} 
