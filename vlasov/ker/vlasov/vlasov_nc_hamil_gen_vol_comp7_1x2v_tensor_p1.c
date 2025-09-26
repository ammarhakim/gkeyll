#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp7_1x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  out[0] += 1.369306393762915*alpha[2]*dv10; 
  out[1] += 1.369306393762915*alpha[4]*dv10; 
  out[2] += (1.224744871391589*alpha[7]+1.369306393762915*alpha[0])*dv10; 
  out[3] += 1.369306393762915*alpha[6]*dv10; 
  out[4] += (1.224744871391589*alpha[10]+1.369306393762915*alpha[1])*dv10; 
  out[5] += 1.369306393762915*alpha[9]*dv10; 
  out[6] += (1.224744871391589*alpha[11]+1.369306393762915*alpha[3])*dv10; 
  out[7] += 1.224744871391589*alpha[2]*dv10; 
  out[8] += 1.369306393762915*alpha[13]*dv10; 
  out[9] += (1.224744871391589*alpha[14]+1.369306393762915*alpha[5])*dv10; 
  out[10] += 1.224744871391589*alpha[4]*dv10; 
  out[11] += 1.224744871391589*alpha[6]*dv10; 
  out[12] += 1.369306393762915*alpha[15]*dv10; 
  out[13] += (1.224744871391589*alpha[16]+1.369306393762915*alpha[8])*dv10; 
  out[14] += 1.224744871391589*alpha[9]*dv10; 
  out[15] += (1.224744871391589*alpha[17]+1.369306393762915*alpha[12])*dv10; 
  out[16] += 1.224744871391589*alpha[13]*dv10; 
  out[17] += 1.224744871391589*alpha[15]*dv10; 
  } 

  if (dir == 2) { 
  } 

} 
