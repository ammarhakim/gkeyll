#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_alpha_quad_vx_1x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[15]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf[4]; 
  double Bz_quad = 0.0;
  Bz_quad = 0.9681844646844028*Bz[2]-1.0546722811938851*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[0] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(6.458522336955395*hamil[8])+5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]-2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[1] += 2.0/(dxv[2]*jacob_vel_surf_vy[1])*(-(2.549857826886422*hamil[8])+1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]-1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[2] += 2.0/(dxv[2]*jacob_vel_surf_vy[2])*(2.549857826886422*hamil[8]-1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]+1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[3] += 2.0/(dxv[2]*jacob_vel_surf_vy[3])*(6.458522336955395*hamil[8]-5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]+2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

  Bz_quad = -(0.5164305132317774*Bz[2])-0.41639003950091297*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[4] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(6.458522336955395*hamil[8])+5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]-2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[5] += 2.0/(dxv[2]*jacob_vel_surf_vy[1])*(-(2.549857826886422*hamil[8])+1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]-1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[6] += 2.0/(dxv[2]*jacob_vel_surf_vy[2])*(2.549857826886422*hamil[8]-1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]+1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[7] += 2.0/(dxv[2]*jacob_vel_surf_vy[3])*(6.458522336955395*hamil[8]-5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]+2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

  Bz_quad = -(0.5164305132317774*Bz[2])+0.41639003950091297*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[8] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(6.458522336955395*hamil[8])+5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]-2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[9] += 2.0/(dxv[2]*jacob_vel_surf_vy[1])*(-(2.549857826886422*hamil[8])+1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]-1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[10] += 2.0/(dxv[2]*jacob_vel_surf_vy[2])*(2.549857826886422*hamil[8]-1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]+1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[11] += 2.0/(dxv[2]*jacob_vel_surf_vy[3])*(6.458522336955395*hamil[8]-5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]+2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

  Bz_quad = 0.9681844646844028*Bz[2]+1.0546722811938851*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[12] += 2.0/(dxv[2]*jacob_vel_surf_vy[0])*(-(6.458522336955395*hamil[8])+5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]-2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[13] += 2.0/(dxv[2]*jacob_vel_surf_vy[1])*(-(2.549857826886422*hamil[8])+1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]-1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[14] += 2.0/(dxv[2]*jacob_vel_surf_vy[2])*(2.549857826886422*hamil[8]-1.9751113797455495*hamil[7]+1.9364916731037085*hamil[6]+1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;
  alpha_quad[15] += 2.0/(dxv[2]*jacob_vel_surf_vy[3])*(6.458522336955395*hamil[8]-5.002749890427371*hamil[7]+1.9364916731037085*hamil[6]+2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2])*Bz_quad;

} 
