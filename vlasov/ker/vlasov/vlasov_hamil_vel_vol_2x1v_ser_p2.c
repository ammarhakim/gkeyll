#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_vol_2x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  out[1] += (4.743416490252569*hamil[2]*f[3]+2.121320343559642*f[0]*hamil[1])*dv10*dx10; 
  out[4] += (4.743416490252569*hamil[2]*f[6]+2.121320343559642*hamil[1]*f[2])*dv10*dx10; 
  out[5] += (4.242640687119286*hamil[2]*f[9]+2.121320343559642*hamil[1]*f[3]+4.743416490252569*f[0]*hamil[2])*dv10*dx10; 
  out[7] += (10.60660171779821*hamil[2]*f[5]+4.743416490252569*f[1]*hamil[1])*dv10*dx10; 
  out[10] += (4.242640687119286*hamil[2]*f[16]+2.121320343559642*hamil[1]*f[6]+4.743416490252569*f[2]*hamil[2])*dv10*dx10; 
  out[11] += (10.60660171779821*hamil[2]*f[10]+4.743416490252569*hamil[1]*f[4])*dv10*dx10; 
  out[12] += (4.743416490252569*hamil[2]*f[14]+2.121320343559642*hamil[1]*f[8])*dv10*dx10; 
  out[13] += (9.48683298050514*hamil[2]*f[15]+4.743416490252569*hamil[1]*f[5]+10.60660171779821*f[1]*hamil[2])*dv10*dx10; 
  out[15] += (2.121320343559642*hamil[1]*f[9]+4.242640687119286*hamil[2]*f[3])*dv10*dx10; 
  out[17] += (9.48683298050514*hamil[2]*f[19]+4.743416490252569*hamil[1]*f[10]+10.60660171779821*hamil[2]*f[4])*dv10*dx10; 
  out[18] += (2.121320343559642*hamil[1]*f[14]+4.743416490252569*hamil[2]*f[8])*dv10*dx10; 
  out[19] += (2.121320343559642*hamil[1]*f[16]+4.242640687119286*hamil[2]*f[6])*dv10*dx10; 

} 
