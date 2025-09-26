#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  double Ghat[3]; 
  double Ghat_nodal[8]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[4] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[5] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[6] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[7] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.24597051986529*Ghat_nodal[3]+0.4611362613212576*Ghat_nodal[2]+0.4611362613212576*Ghat_nodal[1]+0.24597051986529*Ghat_nodal[0]; 
  Ghat[1] = 0.3668728630454641*Ghat_nodal[3]+0.271546746793545*Ghat_nodal[2]-0.271546746793545*Ghat_nodal[1]-0.3668728630454641*Ghat_nodal[0]; 
  Ghat[2] = 0.3367876570272824*Ghat_nodal[3]-0.3367876570272824*Ghat_nodal[2]-0.3367876570272824*Ghat_nodal[1]+0.3367876570272824*Ghat_nodal[0]; 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += -0.7071067811865475*Ghat[2]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-2.738612787525831*hamil[2])*dv10; 
  f_l_quad = (-0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[4] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[5] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[6] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(2.738612787525831*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[7] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.24597051986529*Ghat_nodal[3]+0.4611362613212576*Ghat_nodal[2]+0.4611362613212576*Ghat_nodal[1]+0.24597051986529*Ghat_nodal[0]; 
  Ghat[1] = 0.3668728630454641*Ghat_nodal[3]+0.271546746793545*Ghat_nodal[2]-0.271546746793545*Ghat_nodal[1]-0.3668728630454641*Ghat_nodal[0]; 
  Ghat[2] = 0.3367876570272824*Ghat_nodal[3]-0.3367876570272824*Ghat_nodal[2]-0.3367876570272824*Ghat_nodal[1]+0.3367876570272824*Ghat_nodal[0]; 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += 0.7071067811865475*Ghat[2]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 

  } 
  return 0.0;

} 
