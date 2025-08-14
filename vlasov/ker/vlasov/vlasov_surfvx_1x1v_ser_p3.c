#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Ghat_l = &flux_l[0]; 
  const double *Ghat_r = &flux_r[0]; 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10); 
  out[3] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10); 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[5] += (1.5811388300841895*Ghat_l[0]-1.5811388300841895*Ghat_r[0])*dv10; 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10); 
  out[7] += (1.5811388300841898*Ghat_l[1]-1.5811388300841898*Ghat_r[1])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[9] += -(1.8708286933869707*(Ghat_r[0]+Ghat_l[0])*dv10); 
  out[10] += -(1.2247448713915887*(Ghat_r[3]+Ghat_l[3])*dv10); 
  out[11] += -(1.8708286933869707*(Ghat_r[1]+Ghat_l[1])*dv10); 

  return 0.0;

} 
