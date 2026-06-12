#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_2x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[3] += dv10*(0.8660254037844386*Ex[3]*f[6]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[2]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[1]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[0]*jacob_vx_inv); 
  out[7] += dv10*(0.8660254037844386*Ex[2]*f[6]*jacob_vx_inv+0.8660254037844386*f[2]*Ex[3]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[1]*jacob_vx_inv+0.8660254037844386*f[0]*Ex[1]*jacob_vx_inv); 
  out[8] += dv10*(0.8660254037844386*Ex[1]*f[6]*jacob_vx_inv+0.8660254037844386*f[1]*Ex[3]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[2]*jacob_vx_inv+0.8660254037844386*f[0]*Ex[2]*jacob_vx_inv); 
  out[11] += dv10*(0.8660254037844386*Ex[3]*f[17]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[10]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[9]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[4]*jacob_vx_inv); 
  out[14] += dv10*(0.8660254037844386*Ex[3]*f[20]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[13]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[12]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[5]*jacob_vx_inv); 
  out[16] += dv10*(0.8660254037844386*Ex[0]*f[6]*jacob_vx_inv+0.8660254037844386*f[0]*Ex[3]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[2]*jacob_vx_inv+0.8660254037844386*f[1]*Ex[2]*jacob_vx_inv); 
  out[18] += dv10*(0.8660254037844386*Ex[2]*f[17]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[10]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[9]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[4]*jacob_vx_inv); 
  out[19] += dv10*(0.8660254037844386*Ex[1]*f[17]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[10]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[9]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[4]*jacob_vx_inv); 
  out[21] += dv10*(0.8660254037844386*Ex[2]*f[20]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[13]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[12]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[5]*jacob_vx_inv); 
  out[22] += dv10*(0.8660254037844386*Ex[1]*f[20]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[13]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[12]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[5]*jacob_vx_inv); 
  out[25] += dv10*(0.8660254037844386*Ex[3]*f[28]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[24]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[23]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[15]*jacob_vx_inv); 
  out[26] += dv10*(0.8660254037844386*Ex[0]*f[17]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[10]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[9]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[4]*jacob_vx_inv); 
  out[27] += dv10*(0.8660254037844386*Ex[0]*f[20]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[13]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[12]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[5]*jacob_vx_inv); 
  out[29] += dv10*(0.8660254037844386*Ex[2]*f[28]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[24]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[23]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[15]*jacob_vx_inv); 
  out[30] += dv10*(0.8660254037844386*Ex[1]*f[28]*jacob_vx_inv+0.8660254037844386*Ex[0]*f[24]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[23]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[15]*jacob_vx_inv); 
  out[31] += dv10*(0.8660254037844386*Ex[0]*f[28]*jacob_vx_inv+0.8660254037844386*Ex[1]*f[24]*jacob_vx_inv+0.8660254037844386*Ex[2]*f[23]*jacob_vx_inv+0.8660254037844386*Ex[3]*f[15]*jacob_vx_inv); 

  double dv11 = 2.0/dxv[3]; 
  const double *Ey = &qmem[4]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[4] += dv11*(0.8660254037844386*Ey[3]*f[6]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[2]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[1]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[0]*jacob_vy_inv); 
  out[9] += dv11*(0.8660254037844386*Ey[2]*f[6]*jacob_vy_inv+0.8660254037844386*f[2]*Ey[3]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[1]*jacob_vy_inv+0.8660254037844386*f[0]*Ey[1]*jacob_vy_inv); 
  out[10] += dv11*(0.8660254037844386*Ey[1]*f[6]*jacob_vy_inv+0.8660254037844386*f[1]*Ey[3]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[2]*jacob_vy_inv+0.8660254037844386*f[0]*Ey[2]*jacob_vy_inv); 
  out[11] += dv11*(0.8660254037844386*Ey[3]*f[16]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[8]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[7]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[3]*jacob_vy_inv); 
  out[15] += dv11*(0.8660254037844386*Ey[3]*f[20]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[13]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[12]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[5]*jacob_vy_inv); 
  out[17] += dv11*(0.8660254037844386*Ey[0]*f[6]*jacob_vy_inv+0.8660254037844386*f[0]*Ey[3]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[2]*jacob_vy_inv+0.8660254037844386*f[1]*Ey[2]*jacob_vy_inv); 
  out[18] += dv11*(0.8660254037844386*Ey[2]*f[16]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[8]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[7]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[3]*jacob_vy_inv); 
  out[19] += dv11*(0.8660254037844386*Ey[1]*f[16]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[8]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[7]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[3]*jacob_vy_inv); 
  out[23] += dv11*(0.8660254037844386*Ey[2]*f[20]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[13]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[12]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[5]*jacob_vy_inv); 
  out[24] += dv11*(0.8660254037844386*Ey[1]*f[20]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[13]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[12]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[5]*jacob_vy_inv); 
  out[25] += dv11*(0.8660254037844386*Ey[3]*f[27]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[22]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[21]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[14]*jacob_vy_inv); 
  out[26] += dv11*(0.8660254037844386*Ey[0]*f[16]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[8]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[7]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[3]*jacob_vy_inv); 
  out[28] += dv11*(0.8660254037844386*Ey[0]*f[20]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[13]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[12]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[5]*jacob_vy_inv); 
  out[29] += dv11*(0.8660254037844386*Ey[2]*f[27]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[22]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[21]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[14]*jacob_vy_inv); 
  out[30] += dv11*(0.8660254037844386*Ey[1]*f[27]*jacob_vy_inv+0.8660254037844386*Ey[0]*f[22]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[21]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[14]*jacob_vy_inv); 
  out[31] += dv11*(0.8660254037844386*Ey[0]*f[27]*jacob_vy_inv+0.8660254037844386*Ey[1]*f[22]*jacob_vy_inv+0.8660254037844386*Ey[2]*f[21]*jacob_vy_inv+0.8660254037844386*Ey[3]*f[14]*jacob_vy_inv); 

  double dv12 = 2.0/dxv[4]; 
  const double *Ez = &qmem[8]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[5] += dv12*(0.8660254037844386*Ez[3]*f[6]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[2]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[1]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[0]*jacob_vz_inv); 
  out[12] += dv12*(0.8660254037844386*Ez[2]*f[6]*jacob_vz_inv+0.8660254037844386*f[2]*Ez[3]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[1]*jacob_vz_inv+0.8660254037844386*f[0]*Ez[1]*jacob_vz_inv); 
  out[13] += dv12*(0.8660254037844386*Ez[1]*f[6]*jacob_vz_inv+0.8660254037844386*f[1]*Ez[3]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[2]*jacob_vz_inv+0.8660254037844386*f[0]*Ez[2]*jacob_vz_inv); 
  out[14] += dv12*(0.8660254037844386*Ez[3]*f[16]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[8]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[7]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[3]*jacob_vz_inv); 
  out[15] += dv12*(0.8660254037844386*Ez[3]*f[17]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[10]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[9]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[4]*jacob_vz_inv); 
  out[20] += dv12*(0.8660254037844386*Ez[0]*f[6]*jacob_vz_inv+0.8660254037844386*f[0]*Ez[3]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[2]*jacob_vz_inv+0.8660254037844386*f[1]*Ez[2]*jacob_vz_inv); 
  out[21] += dv12*(0.8660254037844386*Ez[2]*f[16]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[8]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[7]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[3]*jacob_vz_inv); 
  out[22] += dv12*(0.8660254037844386*Ez[1]*f[16]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[8]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[7]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[3]*jacob_vz_inv); 
  out[23] += dv12*(0.8660254037844386*Ez[2]*f[17]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[10]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[9]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[4]*jacob_vz_inv); 
  out[24] += dv12*(0.8660254037844386*Ez[1]*f[17]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[10]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[9]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[4]*jacob_vz_inv); 
  out[25] += dv12*(0.8660254037844386*Ez[3]*f[26]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[19]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[18]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[11]*jacob_vz_inv); 
  out[27] += dv12*(0.8660254037844386*Ez[0]*f[16]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[8]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[7]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[3]*jacob_vz_inv); 
  out[28] += dv12*(0.8660254037844386*Ez[0]*f[17]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[10]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[9]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[4]*jacob_vz_inv); 
  out[29] += dv12*(0.8660254037844386*Ez[2]*f[26]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[19]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[18]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[11]*jacob_vz_inv); 
  out[30] += dv12*(0.8660254037844386*Ez[1]*f[26]*jacob_vz_inv+0.8660254037844386*Ez[0]*f[19]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[18]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[11]*jacob_vz_inv); 
  out[31] += dv12*(0.8660254037844386*Ez[0]*f[26]*jacob_vz_inv+0.8660254037844386*Ez[1]*f[19]*jacob_vz_inv+0.8660254037844386*Ez[2]*f[18]*jacob_vz_inv+0.8660254037844386*Ez[3]*f[11]*jacob_vz_inv); 

} 
