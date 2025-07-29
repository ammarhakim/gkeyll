#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vy_2x2v_ser_p1(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  alpha_quad[0] += -dx11*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[1] += -dx11*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[2] += -dx11*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[3] += -dx11*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[4] += -dx11*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[5] += -dx11*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[6] += -dx11*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 
  alpha_quad[7] += -dx11*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2]); 

} 
