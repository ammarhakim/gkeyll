#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[2]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[4]; 
  const double *poisson_tensor_conf_x3 = &poisson_tensor_conf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[10]; 
  alpha_quad[0] += (1.060660171779821*hamil[7]-1.060660171779821*hamil[6]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7071067811865475*poisson_tensor_conf_1[1])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7071067811865475*poisson_tensor_conf_0[1])*(wx1-0.2886751345948129*dv0))-1.0*(0.7071067811865475*poisson_tensor_conf_x1[0]-0.7071067811865475*poisson_tensor_conf_x1[1])*(1.060660171779821*hamil[7]-1.060660171779821*hamil[5]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[1])*dx10; 
  alpha_quad[1] += (1.060660171779821*hamil[7]-1.060660171779821*hamil[6]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7071067811865475*poisson_tensor_conf_1[1])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7071067811865475*poisson_tensor_conf_0[1])*(wx1+0.2886751345948129*dv0))-1.0*(0.7071067811865475*poisson_tensor_conf_x1[0]-0.7071067811865475*poisson_tensor_conf_x1[1])*((-1.060660171779821*hamil[7])-1.060660171779821*hamil[5]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[1])*dx10; 
  alpha_quad[2] += ((-1.060660171779821*hamil[7])-1.060660171779821*hamil[6]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv0))-1.0*(0.7071067811865475*poisson_tensor_conf_x1[1]+0.7071067811865475*poisson_tensor_conf_x1[0])*(1.060660171779821*hamil[7]-1.060660171779821*hamil[5]-0.6123724356957944*hamil[4]+0.6123724356957944*hamil[1])*dx10; 
  alpha_quad[3] += ((-1.060660171779821*hamil[7])-1.060660171779821*hamil[6]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[2])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv0))-1.0*(0.7071067811865475*poisson_tensor_conf_x1[1]+0.7071067811865475*poisson_tensor_conf_x1[0])*((-1.060660171779821*hamil[7])-1.060660171779821*hamil[5]+0.6123724356957944*hamil[4]+0.6123724356957944*hamil[1])*dx10; 

} 
