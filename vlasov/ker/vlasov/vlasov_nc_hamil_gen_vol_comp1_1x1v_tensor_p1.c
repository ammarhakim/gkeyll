#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_1x1v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  
  if (dir == 0) { 
  out[0] += 0.8660254037844386*alpha[0]*dx10; 
  out[1] += 0.8660254037844386*alpha[1]*dx10; 
  out[2] += 0.8660254037844386*alpha[2]*dx10; 
  out[3] += 0.8660254037844386*alpha[3]*dx10; 
  out[4] += 0.8660254037844386*alpha[4]*dx10; 
  out[5] += 0.8660254037844386*alpha[5]*dx10; 
  } 

  if (dir == 1) { 
  } 

} 
