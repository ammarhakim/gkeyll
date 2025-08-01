#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[2] += (-(4.743416490252569*f[1]*phi[2])-2.1213203435596424*f[0]*phi[1])*dv10*dx10; 
  out[3] += (-(4.242640687119286*phi[2]*f[4])-4.743416490252569*f[0]*phi[2]-2.1213203435596424*f[1]*phi[1])*dv10*dx10; 
  out[5] += (-(10.606601717798211*phi[2]*f[3])-4.743416490252569*phi[1]*f[2])*dv10*dx10; 
  out[6] += (-(2.1213203435596424*phi[1]*f[4])-4.242640687119286*f[1]*phi[2])*dv10*dx10; 
  out[7] += (-(9.48683298050514*phi[2]*f[6])-4.743416490252569*phi[1]*f[3]-10.606601717798213*f[2]*phi[2])*dv10*dx10; 

} 
