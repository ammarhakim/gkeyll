#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_1x2v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += dv10*(1.224744871391589*Ex[2]*f[7]*jacob_vx_inv+1.224744871391589*Ex[1]*f[1]*jacob_vx_inv+1.224744871391589*Ex[0]*f[0]*jacob_vx_inv); 
  out[4] += dv10*(Ex[1]*(1.0954451150103324*f[7]*jacob_vx_inv+1.224744871391589*f[0]*jacob_vx_inv)+1.0954451150103324*f[1]*Ex[2]*jacob_vx_inv+1.224744871391589*Ex[0]*f[1]*jacob_vx_inv); 
  out[6] += dv10*(1.224744871391589*Ex[2]*f[13]*jacob_vx_inv+1.224744871391589*Ex[1]*f[5]*jacob_vx_inv+1.224744871391589*Ex[0]*f[3]*jacob_vx_inv); 
  out[8] += dv10*(2.7386127875258306*Ex[2]*f[11]*jacob_vx_inv+2.7386127875258306*Ex[1]*f[4]*jacob_vx_inv+2.7386127875258306*Ex[0]*f[2]*jacob_vx_inv); 
  out[10] += dv10*(Ex[1]*(1.0954451150103324*f[13]*jacob_vx_inv+1.224744871391589*f[3]*jacob_vx_inv)+1.0954451150103324*Ex[2]*f[5]*jacob_vx_inv+1.224744871391589*Ex[0]*f[5]*jacob_vx_inv); 
  out[11] += dv10*(Ex[2]*(0.7824607964359517*f[7]*jacob_vx_inv+1.224744871391589*f[0]*jacob_vx_inv)+1.224744871391589*Ex[0]*f[7]*jacob_vx_inv+1.0954451150103324*Ex[1]*f[1]*jacob_vx_inv); 
  out[12] += dv10*(Ex[1]*(2.4494897427831783*f[11]*jacob_vx_inv+2.7386127875258306*f[2]*jacob_vx_inv)+2.4494897427831783*Ex[2]*f[4]*jacob_vx_inv+2.7386127875258306*Ex[0]*f[4]*jacob_vx_inv); 
  out[14] += dv10*(2.7386127875258306*Ex[2]*f[17]*jacob_vx_inv+2.7386127875258306*Ex[1]*f[10]*jacob_vx_inv+2.7386127875258306*Ex[0]*f[6]*jacob_vx_inv); 
  out[16] += dv10*(1.224744871391589*Ex[1]*f[15]*jacob_vx_inv+1.224744871391589*Ex[0]*f[9]*jacob_vx_inv); 
  out[17] += dv10*(Ex[2]*(0.7824607964359517*f[13]*jacob_vx_inv+1.224744871391589*f[3]*jacob_vx_inv)+1.224744871391589*Ex[0]*f[13]*jacob_vx_inv+1.0954451150103324*Ex[1]*f[5]*jacob_vx_inv); 
  out[18] += dv10*(Ex[1]*(2.4494897427831783*f[17]*jacob_vx_inv+2.7386127875258306*f[6]*jacob_vx_inv)+2.4494897427831783*Ex[2]*f[10]*jacob_vx_inv+2.7386127875258306*Ex[0]*f[10]*jacob_vx_inv); 
  out[19] += dv10*(1.0954451150103324*Ex[2]*f[15]*jacob_vx_inv+1.224744871391589*Ex[0]*f[15]*jacob_vx_inv+1.224744871391589*Ex[1]*f[9]*jacob_vx_inv); 

  double dv11 = 2.0/dxv[2]; 
  const double *Ey = &qmem[3]; 
  const double *jacob_vy = &jacob_vel[3]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[3] += dv11*(1.224744871391589*Ey[2]*f[7]*jacob_vy_inv+1.224744871391589*Ey[1]*f[1]*jacob_vy_inv+1.224744871391589*Ey[0]*f[0]*jacob_vy_inv); 
  out[5] += dv11*(Ey[1]*(1.0954451150103324*f[7]*jacob_vy_inv+1.224744871391589*f[0]*jacob_vy_inv)+1.0954451150103324*f[1]*Ey[2]*jacob_vy_inv+1.224744871391589*Ey[0]*f[1]*jacob_vy_inv); 
  out[6] += dv11*(1.224744871391589*Ey[2]*f[11]*jacob_vy_inv+1.224744871391589*Ey[1]*f[4]*jacob_vy_inv+1.224744871391589*Ey[0]*f[2]*jacob_vy_inv); 
  out[9] += dv11*(2.7386127875258306*Ey[2]*f[13]*jacob_vy_inv+2.7386127875258306*Ey[1]*f[5]*jacob_vy_inv+2.7386127875258306*Ey[0]*f[3]*jacob_vy_inv); 
  out[10] += dv11*(Ey[1]*(1.0954451150103324*f[11]*jacob_vy_inv+1.224744871391589*f[2]*jacob_vy_inv)+1.0954451150103324*Ey[2]*f[4]*jacob_vy_inv+1.224744871391589*Ey[0]*f[4]*jacob_vy_inv); 
  out[13] += dv11*(Ey[2]*(0.7824607964359517*f[7]*jacob_vy_inv+1.224744871391589*f[0]*jacob_vy_inv)+1.224744871391589*Ey[0]*f[7]*jacob_vy_inv+1.0954451150103324*Ey[1]*f[1]*jacob_vy_inv); 
  out[14] += dv11*(1.224744871391589*Ey[1]*f[12]*jacob_vy_inv+1.224744871391589*Ey[0]*f[8]*jacob_vy_inv); 
  out[15] += dv11*(Ey[1]*(2.4494897427831783*f[13]*jacob_vy_inv+2.7386127875258306*f[3]*jacob_vy_inv)+2.4494897427831783*Ey[2]*f[5]*jacob_vy_inv+2.7386127875258306*Ey[0]*f[5]*jacob_vy_inv); 
  out[16] += dv11*(2.7386127875258306*Ey[2]*f[17]*jacob_vy_inv+2.7386127875258306*Ey[1]*f[10]*jacob_vy_inv+2.7386127875258306*Ey[0]*f[6]*jacob_vy_inv); 
  out[17] += dv11*(Ey[2]*(0.7824607964359517*f[11]*jacob_vy_inv+1.224744871391589*f[2]*jacob_vy_inv)+1.224744871391589*Ey[0]*f[11]*jacob_vy_inv+1.0954451150103324*Ey[1]*f[4]*jacob_vy_inv); 
  out[18] += dv11*(1.0954451150103324*Ey[2]*f[12]*jacob_vy_inv+1.224744871391589*Ey[0]*f[12]*jacob_vy_inv+1.224744871391589*Ey[1]*f[8]*jacob_vy_inv); 
  out[19] += dv11*(Ey[1]*(2.4494897427831783*f[17]*jacob_vy_inv+2.7386127875258306*f[6]*jacob_vy_inv)+2.4494897427831783*Ey[2]*f[10]*jacob_vy_inv+2.7386127875258306*Ey[0]*f[10]*jacob_vy_inv); 

} 
