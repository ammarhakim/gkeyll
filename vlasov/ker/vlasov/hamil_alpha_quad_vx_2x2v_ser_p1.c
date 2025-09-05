#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(-(0.75*hamil[15])+0.75*hamil[13]+0.4330127018922193*hamil[12]+0.75*hamil[11]-0.4330127018922193*hamil[8]-0.75*hamil[6]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[1] += -dx10*(0.75*hamil[15]-0.75*hamil[13]-0.4330127018922193*hamil[12]+0.75*hamil[11]+0.4330127018922193*hamil[8]-0.75*hamil[6]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[2] += -dx10*(0.75*hamil[15]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[11]-0.4330127018922193*hamil[8]-0.75*hamil[6]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[3] += -dx10*(-(0.75*hamil[15])-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[11]+0.4330127018922193*hamil[8]-0.75*hamil[6]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[4] += -dx10*(-(0.75*hamil[15])+0.75*hamil[13]+0.4330127018922193*hamil[12]+0.75*hamil[11]-0.4330127018922193*hamil[8]-0.75*hamil[6]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[5] += -dx10*(0.75*hamil[15]-0.75*hamil[13]-0.4330127018922193*hamil[12]+0.75*hamil[11]+0.4330127018922193*hamil[8]-0.75*hamil[6]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[6] += -dx10*(0.75*hamil[15]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[11]-0.4330127018922193*hamil[8]-0.75*hamil[6]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 
  alpha_quad[7] += -dx10*(-(0.75*hamil[15])-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[11]+0.4330127018922193*hamil[8]-0.75*hamil[6]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1]); 

} 
