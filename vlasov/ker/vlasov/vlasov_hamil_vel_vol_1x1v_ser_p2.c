#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[1] += (4.743416490252569*f[2]*hamil[2]+2.1213203435596424*f[0]*hamil[1])*dv10*dx10; 
  out[3] += (4.242640687119286*hamil[2]*f[5]+4.743416490252569*f[0]*hamil[2]+2.1213203435596424*hamil[1]*f[2])*dv10*dx10; 
  out[4] += (10.606601717798211*hamil[2]*f[3]+4.743416490252569*f[1]*hamil[1])*dv10*dx10; 
  out[6] += (9.48683298050514*hamil[2]*f[7]+4.743416490252569*hamil[1]*f[3]+10.606601717798213*f[1]*hamil[2])*dv10*dx10; 
  out[7] += (2.1213203435596424*hamil[1]*f[5]+4.242640687119286*f[2]*hamil[2])*dv10*dx10; 
} 
