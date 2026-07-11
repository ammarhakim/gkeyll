#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_vel_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[40]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf[0]; 
  double Bz_quad = 0.0;
  Bz_quad = -(0.5999999999999995*Bz[7])-0.5999999999999999*Bz[6]+0.4472135954999579*Bz[5]+0.4472135954999579*Bz[4]+0.9*Bz[3]-0.6708203932499369*Bz[2]-0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[0] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[1] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[2] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = 0.75*Bz[7]-0.5590169943749475*Bz[5]+0.4472135954999579*Bz[4]-0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[3] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[4] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[5] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = -(0.5999999999999995*Bz[7])+0.5999999999999999*Bz[6]+0.4472135954999579*Bz[5]+0.4472135954999579*Bz[4]-0.9*Bz[3]+0.6708203932499369*Bz[2]-0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[6] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[7] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[8] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = 0.75*Bz[6]+0.4472135954999579*Bz[5]-0.5590169943749475*Bz[4]-0.6708203932499369*Bz[2]+0.5*Bz[0];
  alpha_quad[9] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[10] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[11] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = -(0.5590169943749475*Bz[5])-0.5590169943749475*Bz[4]+0.5*Bz[0];
  alpha_quad[12] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[13] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[14] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = -(0.75*Bz[6])+0.4472135954999579*Bz[5]-0.5590169943749475*Bz[4]+0.6708203932499369*Bz[2]+0.5*Bz[0];
  alpha_quad[15] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[16] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[17] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = 0.5999999999999995*Bz[7]-0.5999999999999999*Bz[6]+0.4472135954999579*Bz[5]+0.4472135954999579*Bz[4]-0.9*Bz[3]-0.6708203932499369*Bz[2]+0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[18] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[19] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[20] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = -(0.75*Bz[7])-0.5590169943749475*Bz[5]+0.4472135954999579*Bz[4]+0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[21] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[22] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[23] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

  Bz_quad = 0.5999999999999995*Bz[7]+0.5999999999999999*Bz[6]+0.4472135954999579*Bz[5]+0.4472135954999579*Bz[4]+0.9*Bz[3]+0.6708203932499369*Bz[2]+0.6708203932499369*Bz[1]+0.5*Bz[0];
  alpha_quad[24] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]+4.499999999999999*hamil[6]-2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[25] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;
  alpha_quad[26] -= 2.0/(dxv[2]*jacob_vel_surf_vx[0])*(1.9364916731037085*hamil[7]-4.499999999999999*hamil[6]+2.5980762113533156*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*Bz_quad;

} 
