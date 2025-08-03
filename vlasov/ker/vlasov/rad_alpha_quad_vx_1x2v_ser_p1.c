#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vx = &rad[0]; 

  alpha_quad[0] += 0.8660254037844386*rad_vx[3]-0.5*rad_vx[2]-0.8660254037844386*rad_vx[1]+0.5*rad_vx[0]; 
  alpha_quad[1] += -(0.8660254037844386*rad_vx[3])+0.5*rad_vx[2]-0.8660254037844386*rad_vx[1]+0.5*rad_vx[0]; 
  alpha_quad[2] += 0.8660254037844386*rad_vx[3]-0.5*rad_vx[2]-0.8660254037844386*rad_vx[1]+0.5*rad_vx[0]; 
  alpha_quad[3] += -(0.8660254037844386*rad_vx[3])+0.5*rad_vx[2]-0.8660254037844386*rad_vx[1]+0.5*rad_vx[0]; 

} 
