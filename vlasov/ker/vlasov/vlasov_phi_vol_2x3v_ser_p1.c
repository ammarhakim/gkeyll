#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x3v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[3] += dv10*dx10*(-(1.5*f[2]*phi[3]*jacob_cx_inv*jacob_vx_inv)-1.5*f[0]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[7] += dv10*dx10*(-(1.5*phi[3]*f[6]*jacob_cx_inv*jacob_vx_inv)-1.5*f[1]*phi[1]*jacob_cx_inv*jacob_vx_inv); 
  out[8] += dv10*dx10*(-(1.5*f[0]*phi[3]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[2]*jacob_cx_inv*jacob_vx_inv); 
  out[11] += dv10*dx10*(-(1.5*phi[3]*f[10]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[4]*jacob_cx_inv*jacob_vx_inv); 
  out[14] += dv10*dx10*(-(1.5*phi[3]*f[13]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[5]*jacob_cx_inv*jacob_vx_inv); 
  out[16] += dv10*dx10*(-(1.5*phi[1]*f[6]*jacob_cx_inv*jacob_vx_inv)-1.5*f[1]*phi[3]*jacob_cx_inv*jacob_vx_inv); 
  out[18] += dv10*dx10*(-(1.5*phi[3]*f[17]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[9]*jacob_cx_inv*jacob_vx_inv); 
  out[19] += dv10*dx10*(-(1.5*phi[1]*f[10]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[4]*jacob_cx_inv*jacob_vx_inv); 
  out[21] += dv10*dx10*(-(1.5*phi[3]*f[20]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[12]*jacob_cx_inv*jacob_vx_inv); 
  out[22] += dv10*dx10*(-(1.5*phi[1]*f[13]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[5]*jacob_cx_inv*jacob_vx_inv); 
  out[25] += dv10*dx10*(-(1.5*phi[3]*f[24]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[15]*jacob_cx_inv*jacob_vx_inv); 
  out[26] += dv10*dx10*(-(1.5*phi[1]*f[17]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[9]*jacob_cx_inv*jacob_vx_inv); 
  out[27] += dv10*dx10*(-(1.5*phi[1]*f[20]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[12]*jacob_cx_inv*jacob_vx_inv); 
  out[29] += dv10*dx10*(-(1.5*phi[3]*f[28]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[1]*f[23]*jacob_cx_inv*jacob_vx_inv); 
  out[30] += dv10*dx10*(-(1.5*phi[1]*f[24]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[15]*jacob_cx_inv*jacob_vx_inv); 
  out[31] += dv10*dx10*(-(1.5*phi[1]*f[28]*jacob_cx_inv*jacob_vx_inv)-1.5*phi[3]*f[23]*jacob_cx_inv*jacob_vx_inv); 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  out[4] += dv11*dx11*(-(1.5*f[1]*phi[3]*jacob_cy_inv*jacob_vy_inv)-1.5*f[0]*phi[2]*jacob_cy_inv*jacob_vy_inv); 
  out[9] += dv11*dx11*(-(1.5*f[0]*phi[3]*jacob_cy_inv*jacob_vy_inv)-1.5*f[1]*phi[2]*jacob_cy_inv*jacob_vy_inv); 
  out[10] += dv11*dx11*(-(1.5*phi[3]*f[6]*jacob_cy_inv*jacob_vy_inv)-1.5*f[2]*phi[2]*jacob_cy_inv*jacob_vy_inv); 
  out[11] += dv11*dx11*(-(1.5*phi[3]*f[7]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[3]*jacob_cy_inv*jacob_vy_inv); 
  out[15] += dv11*dx11*(-(1.5*phi[3]*f[12]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[5]*jacob_cy_inv*jacob_vy_inv); 
  out[17] += dv11*dx11*(-(1.5*phi[2]*f[6]*jacob_cy_inv*jacob_vy_inv)-1.5*f[2]*phi[3]*jacob_cy_inv*jacob_vy_inv); 
  out[18] += dv11*dx11*(-(1.5*phi[2]*f[7]*jacob_cy_inv*jacob_vy_inv)-1.5*f[3]*phi[3]*jacob_cy_inv*jacob_vy_inv); 
  out[19] += dv11*dx11*(-(1.5*phi[3]*f[16]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[8]*jacob_cy_inv*jacob_vy_inv); 
  out[23] += dv11*dx11*(-(1.5*phi[2]*f[12]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[3]*f[5]*jacob_cy_inv*jacob_vy_inv); 
  out[24] += dv11*dx11*(-(1.5*phi[3]*f[20]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[13]*jacob_cy_inv*jacob_vy_inv); 
  out[25] += dv11*dx11*(-(1.5*phi[3]*f[21]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[14]*jacob_cy_inv*jacob_vy_inv); 
  out[26] += dv11*dx11*(-(1.5*phi[2]*f[16]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[3]*f[8]*jacob_cy_inv*jacob_vy_inv); 
  out[28] += dv11*dx11*(-(1.5*phi[2]*f[20]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[3]*f[13]*jacob_cy_inv*jacob_vy_inv); 
  out[29] += dv11*dx11*(-(1.5*phi[2]*f[21]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[3]*f[14]*jacob_cy_inv*jacob_vy_inv); 
  out[30] += dv11*dx11*(-(1.5*phi[3]*f[27]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[2]*f[22]*jacob_cy_inv*jacob_vy_inv); 
  out[31] += dv11*dx11*(-(1.5*phi[2]*f[27]*jacob_cy_inv*jacob_vy_inv)-1.5*phi[3]*f[22]*jacob_cy_inv*jacob_vy_inv); 

} 
