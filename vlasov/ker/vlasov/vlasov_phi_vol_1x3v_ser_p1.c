#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[2] += -(2.1213203435596424*f[0]*phi[1]*dv10*dx10); 
  out[5] += -(2.1213203435596424*f[1]*phi[1]*dv10*dx10); 
  out[7] += -(2.1213203435596424*phi[1]*f[3]*dv10*dx10); 
  out[9] += -(2.1213203435596424*phi[1]*f[4]*dv10*dx10); 
  out[11] += -(2.1213203435596424*phi[1]*f[6]*dv10*dx10); 
  out[12] += -(2.1213203435596424*phi[1]*f[8]*dv10*dx10); 
  out[14] += -(2.1213203435596424*phi[1]*f[10]*dv10*dx10); 
  out[15] += -(2.1213203435596424*phi[1]*f[13]*dv10*dx10); 
} 
