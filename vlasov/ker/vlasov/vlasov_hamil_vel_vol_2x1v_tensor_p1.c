#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  out[1] += 2.1213203435596424*f[0]*hamil[1]*dv10*dx10; 
  out[4] += 2.1213203435596424*hamil[1]*f[2]*dv10*dx10; 
  out[5] += 2.1213203435596424*hamil[1]*f[3]*dv10*dx10; 
  out[7] += 2.1213203435596424*hamil[1]*f[6]*dv10*dx10; 
} 
