#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_alpha_1x1v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv0 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wx1 = w[1]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_0[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_0[1]*dv10; 
  out[2] += 3.872983346207417*pt_conf_0[0]*hamil[2]*dv10; 
  out[3] += 3.872983346207417*pt_conf_0[1]*hamil[2]*dv10; 
  } 
  if (dir == 1) { 
  } 

} 
