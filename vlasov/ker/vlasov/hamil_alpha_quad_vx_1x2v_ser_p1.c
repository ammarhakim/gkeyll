#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_alpha_quad_vx_1x2v_ser_p1(const double *dxv, 
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dx10 = 2.0/dxv[0]; 

  alpha_quad[0] += -dx10*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[5]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[1]); 
  alpha_quad[1] += -dx10*(-(1.060660171779821*hamil[7])+0.6123724356957944*hamil[5]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[1]); 
  alpha_quad[2] += -dx10*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[5]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[1]); 
  alpha_quad[3] += -dx10*(-(1.060660171779821*hamil[7])+0.6123724356957944*hamil[5]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[1]); 

} 
