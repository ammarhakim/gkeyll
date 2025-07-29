#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x3v_ser_p1(const double *w, const double *dxv, 
  const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  out[3] += (-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1])*dv10*dx10; 
  out[4] += -(1.5*f[1]*phi[3]*dv11*dx11)-1.5*f[0]*phi[2]*dv11*dx11; 
  out[7] += (-(1.5*phi[3]*f[6])-1.5*f[1]*phi[1])*dv10*dx10; 
  out[8] += (-(1.5*f[0]*phi[3])-1.5*phi[1]*f[2])*dv10*dx10; 
  out[9] += -(1.5*f[0]*phi[3]*dv11*dx11)-1.5*f[1]*phi[2]*dv11*dx11; 
  out[10] += -(1.5*phi[3]*f[6]*dv11*dx11)-1.5*f[2]*phi[2]*dv11*dx11; 
  out[11] += (-(1.5*phi[3]*f[7])-1.5*phi[2]*f[3])*dv11*dx11+(-(1.5*phi[3]*f[10])-1.5*phi[1]*f[4])*dv10*dx10; 
  out[14] += (-(1.5*phi[3]*f[13])-1.5*phi[1]*f[5])*dv10*dx10; 
  out[15] += -(1.5*phi[3]*f[12]*dv11*dx11)-1.5*phi[2]*f[5]*dv11*dx11; 
  out[16] += (-(1.5*phi[1]*f[6])-1.5*f[1]*phi[3])*dv10*dx10; 
  out[17] += -(1.5*phi[2]*f[6]*dv11*dx11)-1.5*f[2]*phi[3]*dv11*dx11; 
  out[18] += (-(1.5*phi[2]*f[7])-1.5*f[3]*phi[3])*dv11*dx11+(-(1.5*phi[3]*f[17])-1.5*phi[1]*f[9])*dv10*dx10; 
  out[19] += (-(1.5*phi[3]*f[16])-1.5*phi[2]*f[8])*dv11*dx11+(-(1.5*phi[1]*f[10])-1.5*phi[3]*f[4])*dv10*dx10; 
  out[21] += (-(1.5*phi[3]*f[20])-1.5*phi[1]*f[12])*dv10*dx10; 
  out[22] += (-(1.5*phi[1]*f[13])-1.5*phi[3]*f[5])*dv10*dx10; 
  out[23] += -(1.5*phi[2]*f[12]*dv11*dx11)-1.5*phi[3]*f[5]*dv11*dx11; 
  out[24] += -(1.5*phi[3]*f[20]*dv11*dx11)-1.5*phi[2]*f[13]*dv11*dx11; 
  out[25] += (-(1.5*phi[3]*f[21])-1.5*phi[2]*f[14])*dv11*dx11+(-(1.5*phi[3]*f[24])-1.5*phi[1]*f[15])*dv10*dx10; 
  out[26] += (-(1.5*phi[2]*f[16])-1.5*phi[3]*f[8])*dv11*dx11+(-(1.5*phi[1]*f[17])-1.5*phi[3]*f[9])*dv10*dx10; 
  out[27] += (-(1.5*phi[1]*f[20])-1.5*phi[3]*f[12])*dv10*dx10; 
  out[28] += -(1.5*phi[2]*f[20]*dv11*dx11)-1.5*phi[3]*f[13]*dv11*dx11; 
  out[29] += (-(1.5*phi[2]*f[21])-1.5*phi[3]*f[14])*dv11*dx11+(-(1.5*phi[3]*f[28])-1.5*phi[1]*f[23])*dv10*dx10; 
  out[30] += (-(1.5*phi[3]*f[27])-1.5*phi[2]*f[22])*dv11*dx11+(-(1.5*phi[1]*f[24])-1.5*phi[3]*f[15])*dv10*dx10; 
  out[31] += (-(1.5*phi[2]*f[27])-1.5*phi[3]*f[22])*dv11*dx11+(-(1.5*phi[1]*f[28])-1.5*phi[3]*f[23])*dv10*dx10; 
} 
