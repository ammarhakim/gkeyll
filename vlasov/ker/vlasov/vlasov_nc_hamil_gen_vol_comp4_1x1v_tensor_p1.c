#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp4_1x1v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  out[0] += 1.936491673103709*alpha[2]*dv10; 
  out[1] += 1.936491673103709*alpha[3]*dv10; 
  out[2] += (1.732050807568877*alpha[4]+1.936491673103709*alpha[0])*dv10; 
  out[3] += (1.732050807568877*alpha[5]+1.936491673103709*alpha[1])*dv10; 
  out[4] += 1.732050807568877*alpha[2]*dv10; 
  out[5] += 1.732050807568877*alpha[3]*dv10; 
  } 

} 
