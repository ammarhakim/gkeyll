#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_alpha_quad_vx_1x1v_ser_p2(const double *dxv, 
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vx = &rad[0]; 

  alpha_quad[0] += 1.5811388300841895*rad_vx[2]-1.224744871391589*rad_vx[1]+0.7071067811865475*rad_vx[0]; 
  alpha_quad[1] += 1.5811388300841895*rad_vx[2]-1.224744871391589*rad_vx[1]+0.7071067811865475*rad_vx[0]; 
  alpha_quad[2] += 1.5811388300841895*rad_vx[2]-1.224744871391589*rad_vx[1]+0.7071067811865475*rad_vx[0]; 
  alpha_quad[3] += 1.5811388300841895*rad_vx[2]-1.224744871391589*rad_vx[1]+0.7071067811865475*rad_vx[0]; 

} 
