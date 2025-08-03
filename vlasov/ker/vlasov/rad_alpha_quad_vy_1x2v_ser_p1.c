#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_alpha_quad_vy_1x2v_ser_p1(const double *dxv, 
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *rad_vy = &rad[4]; 

  alpha_quad[0] += 0.8660254037844386*rad_vy[3]-0.8660254037844386*rad_vy[2]-0.5*rad_vy[1]+0.5*rad_vy[0]; 
  alpha_quad[1] += -(0.8660254037844386*rad_vy[3])-0.8660254037844386*rad_vy[2]+0.5*rad_vy[1]+0.5*rad_vy[0]; 
  alpha_quad[2] += 0.8660254037844386*rad_vy[3]-0.8660254037844386*rad_vy[2]-0.5*rad_vy[1]+0.5*rad_vy[0]; 
  alpha_quad[3] += -(0.8660254037844386*rad_vy[3])-0.8660254037844386*rad_vy[2]+0.5*rad_vy[1]+0.5*rad_vy[0]; 

} 
