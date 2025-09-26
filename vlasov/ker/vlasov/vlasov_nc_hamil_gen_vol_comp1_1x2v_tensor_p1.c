#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_1x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  
  if (dir == 0) { 
  out[0] += 0.6123724356957944*alpha[0]*dx10; 
  out[1] += 0.6123724356957944*alpha[1]*dx10; 
  out[2] += 0.6123724356957944*alpha[2]*dx10; 
  out[3] += 0.6123724356957944*alpha[3]*dx10; 
  out[4] += 0.6123724356957944*alpha[4]*dx10; 
  out[5] += 0.6123724356957944*alpha[5]*dx10; 
  out[6] += 0.6123724356957944*alpha[6]*dx10; 
  out[7] += 0.6123724356957944*alpha[7]*dx10; 
  out[8] += 0.6123724356957944*alpha[8]*dx10; 
  out[9] += 0.6123724356957944*alpha[9]*dx10; 
  out[10] += 0.6123724356957944*alpha[10]*dx10; 
  out[11] += 0.6123724356957944*alpha[11]*dx10; 
  out[12] += 0.6123724356957944*alpha[12]*dx10; 
  out[13] += 0.6123724356957944*alpha[13]*dx10; 
  out[14] += 0.6123724356957944*alpha[14]*dx10; 
  out[15] += 0.6123724356957944*alpha[15]*dx10; 
  out[16] += 0.6123724356957944*alpha[16]*dx10; 
  out[17] += 0.6123724356957944*alpha[17]*dx10; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  } 

} 
