#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[3] += dv10*dx10*(-(1.5*f[2]*phi[3]*jacob_vx_inv)-1.5*f[0]*phi[1]*jacob_vx_inv); 
  out[6] += dv10*dx10*(-(1.5*phi[3]*f[5]*jacob_vx_inv)-1.5*f[1]*phi[1]*jacob_vx_inv); 
  out[7] += dv10*dx10*(-(1.5*f[0]*phi[3]*jacob_vx_inv)-1.5*phi[1]*f[2]*jacob_vx_inv); 
  out[10] += dv10*dx10*(-(1.5*phi[3]*f[9]*jacob_vx_inv)-1.5*phi[1]*f[4]*jacob_vx_inv); 
  out[11] += dv10*dx10*(-(1.5*phi[1]*f[5]*jacob_vx_inv)-1.5*f[1]*phi[3]*jacob_vx_inv); 
  out[13] += dv10*dx10*(-(1.5*phi[3]*f[12]*jacob_vx_inv)-1.5*phi[1]*f[8]*jacob_vx_inv); 
  out[14] += dv10*dx10*(-(1.5*phi[1]*f[9]*jacob_vx_inv)-1.5*phi[3]*f[4]*jacob_vx_inv); 
  out[15] += dv10*dx10*(-(1.5*phi[1]*f[12]*jacob_vx_inv)-1.5*phi[3]*f[8]*jacob_vx_inv); 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[4] += dv11*dx11*(-(1.5*f[1]*phi[3]*jacob_vy_inv)-1.5*f[0]*phi[2]*jacob_vy_inv); 
  out[8] += dv11*dx11*(-(1.5*f[0]*phi[3]*jacob_vy_inv)-1.5*f[1]*phi[2]*jacob_vy_inv); 
  out[9] += dv11*dx11*(-(1.5*phi[3]*f[5]*jacob_vy_inv)-1.5*f[2]*phi[2]*jacob_vy_inv); 
  out[10] += dv11*dx11*(-(1.5*phi[3]*f[6]*jacob_vy_inv)-1.5*phi[2]*f[3]*jacob_vy_inv); 
  out[12] += dv11*dx11*(-(1.5*phi[2]*f[5]*jacob_vy_inv)-1.5*f[2]*phi[3]*jacob_vy_inv); 
  out[13] += dv11*dx11*(-(1.5*phi[2]*f[6]*jacob_vy_inv)-1.5*f[3]*phi[3]*jacob_vy_inv); 
  out[14] += dv11*dx11*(-(1.5*phi[3]*f[11]*jacob_vy_inv)-1.5*phi[2]*f[7]*jacob_vy_inv); 
  out[15] += dv11*dx11*(-(1.5*phi[2]*f[11]*jacob_vy_inv)-1.5*phi[3]*f[7]*jacob_vy_inv); 

} 
