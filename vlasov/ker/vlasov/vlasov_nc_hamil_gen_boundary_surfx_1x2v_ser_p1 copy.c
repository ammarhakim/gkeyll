#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Ghat[4]; 
  double Ghat_nodal[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3]); 
  f_l_quad = 0.6123724356957944*fskin[7]+0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = (-0.6123724356957944*fedge[7])+0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])+(1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3]); 
  f_l_quad = (-0.6123724356957944*fskin[7])-0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = 0.6123724356957944*fedge[7]-0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3]); 
  f_l_quad = (-0.6123724356957944*fskin[7])-0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = 0.6123724356957944*fedge[7]-0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Ghat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1]); 
  f_l_quad = 0.6123724356957944*fskin[7]+0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = (-0.6123724356957944*fedge[7])+0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Ghat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.5*Ghat_nodal[3]+0.5*Ghat_nodal[2]+0.5*Ghat_nodal[1]+0.5*Ghat_nodal[0]; 
  Ghat[1] = 0.5*Ghat_nodal[3]+0.5*Ghat_nodal[2]-0.5*Ghat_nodal[1]-0.5*Ghat_nodal[0]; 
  Ghat[2] = 0.5*Ghat_nodal[3]-0.5*Ghat_nodal[2]+0.5*Ghat_nodal[1]-0.5*Ghat_nodal[0]; 
  Ghat[3] = 0.5*Ghat_nodal[3]-0.5*Ghat_nodal[2]-0.5*Ghat_nodal[1]+0.5*Ghat_nodal[0]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dv10*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dv10*dx10; 
  out[6] += -0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10*dx10; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3]); 
  f_l_quad = 0.6123724356957944*fedge[7]+0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = (-0.6123724356957944*fskin[7])+0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])+(0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3]); 
  f_l_quad = (-0.6123724356957944*fedge[7])-0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = 0.6123724356957944*fskin[7]-0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3]); 
  f_l_quad = (-0.6123724356957944*fedge[7])-0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = 0.6123724356957944*fskin[7]-0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Ghat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1]); 
  f_l_quad = 0.6123724356957944*fedge[7]+0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = (-0.6123724356957944*fskin[7])+0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Ghat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.5*Ghat_nodal[3]+0.5*Ghat_nodal[2]+0.5*Ghat_nodal[1]+0.5*Ghat_nodal[0]; 
  Ghat[1] = 0.5*Ghat_nodal[3]+0.5*Ghat_nodal[2]-0.5*Ghat_nodal[1]-0.5*Ghat_nodal[0]; 
  Ghat[2] = 0.5*Ghat_nodal[3]-0.5*Ghat_nodal[2]+0.5*Ghat_nodal[1]-0.5*Ghat_nodal[0]; 
  Ghat[3] = 0.5*Ghat_nodal[3]-0.5*Ghat_nodal[2]-0.5*Ghat_nodal[1]+0.5*Ghat_nodal[0]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dv10*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dv10*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10*dx10; 

  } 
  return 0.0;

} 
