#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_phase_alpha_quad_x_1x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge, 
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  if (hamil_pt_edge == -1) {
    alpha_quad[0] += (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[1] += (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*((-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[2] += (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*((-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[3] += (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*((-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*((-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
  }
  else if (hamil_pt_edge == 1) {
    alpha_quad[0] += (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*((-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*((-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[1] += (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*((-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[2] += (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*((-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
    alpha_quad[3] += (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[5]+0.6123724356957944*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2])*dv10;
  }

} 
