#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  out[2] += -(2.1213203435596424*f[0]*phi[1]*dv10*dx10); 
  out[4] += -(2.1213203435596424*f[1]*phi[1]*dv10*dx10); 
  out[6] += -(2.1213203435596424*phi[1]*f[3]*dv10*dx10); 
  out[7] += -(2.1213203435596424*phi[1]*f[5]*dv10*dx10); 
} 
