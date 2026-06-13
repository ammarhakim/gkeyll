#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp1_2x3v_ser_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  
  if (dir == 0) { 
  out[0] += 0.3061862178478971*alpha[0]*dx10; 
  out[1] += 0.3061862178478971*alpha[1]*dx10; 
  out[2] += 0.3061862178478971*alpha[2]*dx10; 
  out[3] += 0.3061862178478971*alpha[3]*dx10; 
  out[4] += 0.3061862178478971*alpha[4]*dx10; 
  out[5] += 0.3061862178478971*alpha[5]*dx10; 
  out[6] += 0.3061862178478971*alpha[6]*dx10; 
  out[7] += 0.3061862178478971*alpha[7]*dx10; 
  out[8] += 0.3061862178478971*alpha[8]*dx10; 
  out[9] += 0.3061862178478971*alpha[9]*dx10; 
  out[10] += 0.3061862178478971*alpha[10]*dx10; 
  out[11] += 0.3061862178478971*alpha[11]*dx10; 
  out[12] += 0.3061862178478971*alpha[12]*dx10; 
  out[13] += 0.3061862178478971*alpha[13]*dx10; 
  out[14] += 0.3061862178478971*alpha[14]*dx10; 
  out[15] += 0.3061862178478971*alpha[15]*dx10; 
  out[16] += 0.3061862178478971*alpha[16]*dx10; 
  out[17] += 0.3061862178478971*alpha[17]*dx10; 
  out[18] += 0.3061862178478971*alpha[18]*dx10; 
  out[19] += 0.3061862178478971*alpha[19]*dx10; 
  out[20] += 0.3061862178478971*alpha[20]*dx10; 
  out[21] += 0.3061862178478971*alpha[21]*dx10; 
  out[22] += 0.3061862178478971*alpha[22]*dx10; 
  out[23] += 0.3061862178478971*alpha[23]*dx10; 
  out[24] += 0.3061862178478971*alpha[24]*dx10; 
  out[25] += 0.3061862178478971*alpha[25]*dx10; 
  out[26] += 0.3061862178478971*alpha[26]*dx10; 
  out[27] += 0.3061862178478971*alpha[27]*dx10; 
  out[28] += 0.3061862178478971*alpha[28]*dx10; 
  out[29] += 0.3061862178478971*alpha[29]*dx10; 
  out[30] += 0.3061862178478971*alpha[30]*dx10; 
  out[31] += 0.3061862178478971*alpha[31]*dx10; 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  } 

  if (dir == 3) { 
  } 

  if (dir == 4) { 
  } 

} 
