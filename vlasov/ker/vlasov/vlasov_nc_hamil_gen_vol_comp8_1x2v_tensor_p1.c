#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_comp8_1x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *alpha, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  
  if (dir == 0) { 
  } 

  if (dir == 1) { 
  } 

  if (dir == 2) { 
  out[0] += 1.369306393762915*alpha[3]*dv11; 
  out[1] += 1.369306393762915*alpha[5]*dv11; 
  out[2] += 1.369306393762915*alpha[6]*dv11; 
  out[3] += (1.224744871391589*alpha[8]+1.369306393762915*alpha[0])*dv11; 
  out[4] += 1.369306393762915*alpha[9]*dv11; 
  out[5] += (1.224744871391589*alpha[12]+1.369306393762915*alpha[1])*dv11; 
  out[6] += (1.224744871391589*alpha[13]+1.369306393762915*alpha[2])*dv11; 
  out[7] += 1.369306393762915*alpha[11]*dv11; 
  out[8] += 1.224744871391589*alpha[3]*dv11; 
  out[9] += (1.224744871391589*alpha[15]+1.369306393762915*alpha[4])*dv11; 
  out[10] += 1.369306393762915*alpha[14]*dv11; 
  out[11] += (1.224744871391589*alpha[16]+1.369306393762915*alpha[7])*dv11; 
  out[12] += 1.224744871391589*alpha[5]*dv11; 
  out[13] += 1.224744871391589*alpha[6]*dv11; 
  out[14] += (1.224744871391589*alpha[17]+1.369306393762915*alpha[10])*dv11; 
  out[15] += 1.224744871391589*alpha[9]*dv11; 
  out[16] += 1.224744871391589*alpha[11]*dv11; 
  out[17] += 1.224744871391589*alpha[14]*dv11; 
  } 

} 
