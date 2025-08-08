#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  out[3] += dv10*dx10*(-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1]); 
  out[6] += dv10*dx10*(-(1.5*phi[3]*f[5])-1.5*f[1]*phi[1]); 
  out[7] += dv10*dx10*(-(1.5*f[0]*phi[3])-1.5*phi[1]*f[2]); 
  out[10] += dv10*dx10*(-(1.5*phi[3]*f[9])-1.5*phi[1]*f[4]); 
  out[11] += dv10*dx10*(-(1.5*phi[1]*f[5])-1.5*f[1]*phi[3]); 
  out[13] += dv10*dx10*(-(1.5*phi[3]*f[12])-1.5*phi[1]*f[8]); 
  out[14] += dv10*dx10*(-(1.5*phi[1]*f[9])-1.5*phi[3]*f[4]); 
  out[15] += dv10*dx10*(-(1.5*phi[1]*f[12])-1.5*phi[3]*f[8]); 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  out[4] += dv11*dx11*(-(1.5*f[1]*phi[3])-1.5*f[0]*phi[2]); 
  out[8] += dv11*dx11*(-(1.5*f[0]*phi[3])-1.5*f[1]*phi[2]); 
  out[9] += dv11*dx11*(-(1.5*phi[3]*f[5])-1.5*f[2]*phi[2]); 
  out[10] += dv11*dx11*(-(1.5*phi[3]*f[6])-1.5*phi[2]*f[3]); 
  out[12] += dv11*dx11*(-(1.5*phi[2]*f[5])-1.5*f[2]*phi[3]); 
  out[13] += dv11*dx11*(-(1.5*phi[2]*f[6])-1.5*f[3]*phi[3]); 
  out[14] += dv11*dx11*(-(1.5*phi[3]*f[11])-1.5*phi[2]*f[7]); 
  out[15] += dv11*dx11*(-(1.5*phi[2]*f[11])-1.5*phi[3]*f[7]); 

} 
