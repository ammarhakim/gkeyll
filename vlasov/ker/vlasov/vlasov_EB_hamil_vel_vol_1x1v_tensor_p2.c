#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_EB_hamil_vel_vol_1x1v_tensor_p2(const double *w, const double *dxv, 
  const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  double alpha_vdim[9] = {0.0}; 

  alpha_vdim[0] = dv10*(1.4142135623730951*Ex[0]); 
  alpha_vdim[1] = dv10*(1.4142135623730951*Ex[1]); 
  alpha_vdim[4] = dv10*(1.4142135623730951*Ex[2]); 
  out[2] += 0.8660254037844386*alpha_vdim[4]*f[4]+0.8660254037844386*alpha_vdim[1]*f[1]+0.8660254037844386*alpha_vdim[0]*f[0]; 
  out[3] += 0.7745966692414833*alpha_vdim[1]*f[4]+0.7745966692414833*f[1]*alpha_vdim[4]+0.8660254037844386*alpha_vdim[0]*f[1]+0.8660254037844386*f[0]*alpha_vdim[1]; 
  out[5] += 1.9364916731037085*alpha_vdim[4]*f[6]+1.9364916731037085*alpha_vdim[1]*f[3]+1.9364916731037085*alpha_vdim[0]*f[2]; 
  out[6] += 0.5532833351724881*alpha_vdim[4]*f[4]+0.8660254037844387*alpha_vdim[0]*f[4]+0.8660254037844387*f[0]*alpha_vdim[4]+0.7745966692414834*alpha_vdim[1]*f[1]; 
  out[7] += 1.7320508075688772*alpha_vdim[1]*f[6]+1.7320508075688774*f[3]*alpha_vdim[4]+1.9364916731037085*alpha_vdim[0]*f[3]+1.9364916731037085*alpha_vdim[1]*f[2]; 
  out[8] += 1.237179148263484*alpha_vdim[4]*f[6]+1.9364916731037085*alpha_vdim[0]*f[6]+1.9364916731037085*f[2]*alpha_vdim[4]+1.7320508075688772*alpha_vdim[1]*f[3]; 

} 
