#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[2] += dv10*dx10*((-4.743416490252569*f[1]*phi[2])-2.121320343559642*f[0]*phi[1]); 
  out[3] += dv10*dx10*(phi[2]*((-4.242640687119286*f[4])-4.743416490252569*f[0])-2.121320343559642*f[1]*phi[1]); 
  out[5] += dv10*dx10*((-10.60660171779821*phi[2]*f[3])-4.743416490252569*phi[1]*f[2]); 
  out[6] += dv10*dx10*((-2.121320343559642*phi[1]*f[4])-4.242640687119286*f[1]*phi[2]); 
  out[7] += dv10*dx10*(phi[2]*((-9.48683298050514*f[6])-10.60660171779821*f[2])-4.743416490252569*phi[1]*f[3]); 

} 
