#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_EB_hamil_vel_vol_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Ey = &qmem[2]; 
  const double *Bz = &qmem[10]; 
  double alpha_vdim[8] = {0.0}; 

  alpha_vdim[0] = dv10*(1.7320508075688772*Bz[0]*hamil[2]*dv11+2.0*Ex[0]); 
  alpha_vdim[1] = dv10*(1.7320508075688772*Bz[1]*hamil[2]*dv11+2.0*Ex[1]); 
  alpha_vdim[2] = dv10*(1.7320508075688772*Bz[0]*hamil[3]*dv11); 
  alpha_vdim[4] = dv10*(1.7320508075688772*Bz[1]*hamil[3]*dv11); 
  out[2] += 0.6123724356957944*alpha_vdim[4]*f[4]+0.6123724356957944*alpha_vdim[2]*f[2]+0.6123724356957944*alpha_vdim[1]*f[1]+0.6123724356957944*alpha_vdim[0]*f[0]; 
  out[4] += 0.6123724356957944*alpha_vdim[2]*f[4]+0.6123724356957944*f[2]*alpha_vdim[4]+0.6123724356957944*alpha_vdim[0]*f[1]+0.6123724356957944*f[0]*alpha_vdim[1]; 
  out[6] += 0.6123724356957944*alpha_vdim[4]*f[7]+0.6123724356957944*alpha_vdim[2]*f[6]+0.6123724356957944*alpha_vdim[1]*f[5]+0.6123724356957944*alpha_vdim[0]*f[3]; 
  out[7] += 0.6123724356957944*alpha_vdim[2]*f[7]+0.6123724356957944*alpha_vdim[4]*f[6]+0.6123724356957944*alpha_vdim[0]*f[5]+0.6123724356957944*alpha_vdim[1]*f[3]; 

  alpha_vdim[0] = dv11*(2.0*Ey[0]-1.7320508075688772*Bz[0]*hamil[1]*dv10); 
  alpha_vdim[1] = dv11*(2.0*Ey[1]-1.7320508075688772*Bz[1]*hamil[1]*dv10); 
  alpha_vdim[3] = dv11*(-(1.7320508075688772*Bz[0]*hamil[3]*dv10)); 
  alpha_vdim[5] = dv11*(-(1.7320508075688772*Bz[1]*hamil[3]*dv10)); 
  out[3] += 0.6123724356957944*alpha_vdim[5]*f[5]+0.6123724356957944*alpha_vdim[3]*f[3]+0.6123724356957944*alpha_vdim[1]*f[1]+0.6123724356957944*alpha_vdim[0]*f[0]; 
  out[5] += 0.6123724356957944*alpha_vdim[3]*f[5]+0.6123724356957944*f[3]*alpha_vdim[5]+0.6123724356957944*alpha_vdim[0]*f[1]+0.6123724356957944*f[0]*alpha_vdim[1]; 
  out[6] += 0.6123724356957944*alpha_vdim[5]*f[7]+0.6123724356957944*alpha_vdim[3]*f[6]+0.6123724356957944*alpha_vdim[1]*f[4]+0.6123724356957944*alpha_vdim[0]*f[2]; 
  out[7] += 0.6123724356957944*alpha_vdim[3]*f[7]+0.6123724356957944*alpha_vdim[5]*f[6]+0.6123724356957944*alpha_vdim[0]*f[4]+0.6123724356957944*alpha_vdim[1]*f[2]; 

} 
