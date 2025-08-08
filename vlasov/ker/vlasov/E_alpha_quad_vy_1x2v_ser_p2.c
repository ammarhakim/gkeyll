#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[3]; 

  double force_quad = 0.0;
  force_quad = 0.9681844646844028*Ey[2]-1.0546722811938851*Ey[1]+0.7071067811865475*Ey[0];
  alpha_quad[0] += force_quad;
  alpha_quad[1] += force_quad;
  alpha_quad[2] += force_quad;
  alpha_quad[3] += force_quad;

  force_quad = -(0.5164305132317774*Ey[2])-0.41639003950091297*Ey[1]+0.7071067811865475*Ey[0];
  alpha_quad[4] += force_quad;
  alpha_quad[5] += force_quad;
  alpha_quad[6] += force_quad;
  alpha_quad[7] += force_quad;

  force_quad = -(0.5164305132317774*Ey[2])+0.41639003950091297*Ey[1]+0.7071067811865475*Ey[0];
  alpha_quad[8] += force_quad;
  alpha_quad[9] += force_quad;
  alpha_quad[10] += force_quad;
  alpha_quad[11] += force_quad;

  force_quad = 0.9681844646844028*Ey[2]+1.0546722811938851*Ey[1]+0.7071067811865475*Ey[0];
  alpha_quad[12] += force_quad;
  alpha_quad[13] += force_quad;
  alpha_quad[14] += force_quad;
  alpha_quad[15] += force_quad;

} 
