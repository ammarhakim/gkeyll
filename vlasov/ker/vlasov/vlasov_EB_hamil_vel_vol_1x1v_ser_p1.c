#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_EB_hamil_vel_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  double alpha_vdim[4] = {0.0}; 

  alpha_vdim[0] = dv10*(1.4142135623730951*Ex[0]); 
  alpha_vdim[1] = dv10*(1.4142135623730951*Ex[1]); 
  out[2] += 0.8660254037844386*alpha_vdim[1]*f[1]+0.8660254037844386*alpha_vdim[0]*f[0]; 
  out[3] += 0.8660254037844386*alpha_vdim[0]*f[1]+0.8660254037844386*f[0]*alpha_vdim[1]; 

} 
