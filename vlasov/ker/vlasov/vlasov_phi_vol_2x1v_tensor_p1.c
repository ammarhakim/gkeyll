#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x1v_tensor_p1(const double *w, const double *dxv, 
  const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  out[3] += (-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1])*dv10*dx10; 
  out[5] += (-(1.5*phi[3]*f[4])-1.5*f[1]*phi[1])*dv10*dx10; 
  out[6] += (-(1.5*f[0]*phi[3])-1.5*phi[1]*f[2])*dv10*dx10; 
  out[7] += (-(1.5*phi[1]*f[4])-1.5*f[1]*phi[3])*dv10*dx10; 
} 
