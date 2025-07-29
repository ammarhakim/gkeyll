#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x2v_tensor_p1(const double *w, const double *dxv, 
  const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  out[3] += (-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1])*dv10*dx10; 
  out[4] += -(1.5*f[1]*phi[3]*dv11*dx11)-1.5*f[0]*phi[2]*dv11*dx11; 
  out[6] += (-(1.5*phi[3]*f[5])-1.5*f[1]*phi[1])*dv10*dx10; 
  out[7] += (-(1.5*f[0]*phi[3])-1.5*phi[1]*f[2])*dv10*dx10; 
  out[8] += -(1.5*f[0]*phi[3]*dv11*dx11)-1.5*f[1]*phi[2]*dv11*dx11; 
  out[9] += -(1.5*phi[3]*f[5]*dv11*dx11)-1.5*f[2]*phi[2]*dv11*dx11; 
  out[10] += (-(1.5*phi[3]*f[6])-1.5*phi[2]*f[3])*dv11*dx11+(-(1.5*phi[3]*f[9])-1.5*phi[1]*f[4])*dv10*dx10; 
  out[11] += (-(1.5*phi[1]*f[5])-1.5*f[1]*phi[3])*dv10*dx10; 
  out[12] += -(1.5*phi[2]*f[5]*dv11*dx11)-1.5*f[2]*phi[3]*dv11*dx11; 
  out[13] += (-(1.5*phi[2]*f[6])-1.5*f[3]*phi[3])*dv11*dx11+(-(1.5*phi[3]*f[12])-1.5*phi[1]*f[8])*dv10*dx10; 
  out[14] += (-(1.5*phi[3]*f[11])-1.5*phi[2]*f[7])*dv11*dx11+(-(1.5*phi[1]*f[9])-1.5*phi[3]*f[4])*dv10*dx10; 
  out[15] += (-(1.5*phi[2]*f[11])-1.5*phi[3]*f[7])*dv11*dx11+(-(1.5*phi[1]*f[12])-1.5*phi[3]*f[8])*dv10*dx10; 
} 
