#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vy_3x3v_ser_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[8]; 

  double force_quad = 0.0;
  force_quad = -(0.3535533905932737*Ey[7])+0.3535533905932737*Ey[6]+0.3535533905932737*Ey[5]+0.3535533905932737*Ey[4]-0.3535533905932737*Ey[3]-0.3535533905932737*Ey[2]-0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[0] += force_quad;
  alpha_quad[1] += force_quad;
  alpha_quad[2] += force_quad;
  alpha_quad[3] += force_quad;

  force_quad = 0.3535533905932737*Ey[7]-0.3535533905932737*Ey[6]-0.3535533905932737*Ey[5]+0.3535533905932737*Ey[4]+0.3535533905932737*Ey[3]-0.3535533905932737*Ey[2]-0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[4] += force_quad;
  alpha_quad[5] += force_quad;
  alpha_quad[6] += force_quad;
  alpha_quad[7] += force_quad;

  force_quad = 0.3535533905932737*Ey[7]-0.3535533905932737*Ey[6]+0.3535533905932737*Ey[5]-0.3535533905932737*Ey[4]-0.3535533905932737*Ey[3]+0.3535533905932737*Ey[2]-0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[8] += force_quad;
  alpha_quad[9] += force_quad;
  alpha_quad[10] += force_quad;
  alpha_quad[11] += force_quad;

  force_quad = -(0.3535533905932737*Ey[7])+0.3535533905932737*Ey[6]-0.3535533905932737*Ey[5]-0.3535533905932737*Ey[4]+0.3535533905932737*Ey[3]+0.3535533905932737*Ey[2]-0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[12] += force_quad;
  alpha_quad[13] += force_quad;
  alpha_quad[14] += force_quad;
  alpha_quad[15] += force_quad;

  force_quad = 0.3535533905932737*Ey[7]+0.3535533905932737*Ey[6]-0.3535533905932737*Ey[5]-0.3535533905932737*Ey[4]-0.3535533905932737*Ey[3]-0.3535533905932737*Ey[2]+0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[16] += force_quad;
  alpha_quad[17] += force_quad;
  alpha_quad[18] += force_quad;
  alpha_quad[19] += force_quad;

  force_quad = -(0.3535533905932737*Ey[7])-0.3535533905932737*Ey[6]+0.3535533905932737*Ey[5]-0.3535533905932737*Ey[4]+0.3535533905932737*Ey[3]-0.3535533905932737*Ey[2]+0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[20] += force_quad;
  alpha_quad[21] += force_quad;
  alpha_quad[22] += force_quad;
  alpha_quad[23] += force_quad;

  force_quad = -(0.3535533905932737*Ey[7])-0.3535533905932737*Ey[6]-0.3535533905932737*Ey[5]+0.3535533905932737*Ey[4]-0.3535533905932737*Ey[3]+0.3535533905932737*Ey[2]+0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[24] += force_quad;
  alpha_quad[25] += force_quad;
  alpha_quad[26] += force_quad;
  alpha_quad[27] += force_quad;

  force_quad = 0.3535533905932737*Ey[7]+0.3535533905932737*Ey[6]+0.3535533905932737*Ey[5]+0.3535533905932737*Ey[4]+0.3535533905932737*Ey[3]+0.3535533905932737*Ey[2]+0.3535533905932737*Ey[1]+0.3535533905932737*Ey[0];
  alpha_quad[28] += force_quad;
  alpha_quad[29] += force_quad;
  alpha_quad[30] += force_quad;
  alpha_quad[31] += force_quad;

} 
