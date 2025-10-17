#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vx_1x3v_tensor_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ex = &qmem[0]; 

  double force_quad = 0.0;
  force_quad = 0.9681844646844028*Ex[2]-1.0546722811938851*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[0] += force_quad;
  alpha_quad[1] += force_quad;
  alpha_quad[2] += force_quad;
  alpha_quad[3] += force_quad;
  alpha_quad[4] += force_quad;
  alpha_quad[5] += force_quad;
  alpha_quad[6] += force_quad;
  alpha_quad[7] += force_quad;
  alpha_quad[8] += force_quad;
  alpha_quad[9] += force_quad;
  alpha_quad[10] += force_quad;
  alpha_quad[11] += force_quad;
  alpha_quad[12] += force_quad;
  alpha_quad[13] += force_quad;
  alpha_quad[14] += force_quad;
  alpha_quad[15] += force_quad;

  force_quad = -(0.5164305132317774*Ex[2])-0.41639003950091297*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[16] += force_quad;
  alpha_quad[17] += force_quad;
  alpha_quad[18] += force_quad;
  alpha_quad[19] += force_quad;
  alpha_quad[20] += force_quad;
  alpha_quad[21] += force_quad;
  alpha_quad[22] += force_quad;
  alpha_quad[23] += force_quad;
  alpha_quad[24] += force_quad;
  alpha_quad[25] += force_quad;
  alpha_quad[26] += force_quad;
  alpha_quad[27] += force_quad;
  alpha_quad[28] += force_quad;
  alpha_quad[29] += force_quad;
  alpha_quad[30] += force_quad;
  alpha_quad[31] += force_quad;

  force_quad = -(0.5164305132317774*Ex[2])+0.41639003950091297*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[32] += force_quad;
  alpha_quad[33] += force_quad;
  alpha_quad[34] += force_quad;
  alpha_quad[35] += force_quad;
  alpha_quad[36] += force_quad;
  alpha_quad[37] += force_quad;
  alpha_quad[38] += force_quad;
  alpha_quad[39] += force_quad;
  alpha_quad[40] += force_quad;
  alpha_quad[41] += force_quad;
  alpha_quad[42] += force_quad;
  alpha_quad[43] += force_quad;
  alpha_quad[44] += force_quad;
  alpha_quad[45] += force_quad;
  alpha_quad[46] += force_quad;
  alpha_quad[47] += force_quad;

  force_quad = 0.9681844646844028*Ex[2]+1.0546722811938851*Ex[1]+0.7071067811865475*Ex[0];
  alpha_quad[48] += force_quad;
  alpha_quad[49] += force_quad;
  alpha_quad[50] += force_quad;
  alpha_quad[51] += force_quad;
  alpha_quad[52] += force_quad;
  alpha_quad[53] += force_quad;
  alpha_quad[54] += force_quad;
  alpha_quad[55] += force_quad;
  alpha_quad[56] += force_quad;
  alpha_quad[57] += force_quad;
  alpha_quad[58] += force_quad;
  alpha_quad[59] += force_quad;
  alpha_quad[60] += force_quad;
  alpha_quad[61] += force_quad;
  alpha_quad[62] += force_quad;
  alpha_quad[63] += force_quad;

} 
