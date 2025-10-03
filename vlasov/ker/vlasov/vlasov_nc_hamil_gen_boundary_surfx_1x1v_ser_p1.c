#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  double Fhat_nodal[2]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 1.224744871391589*hamil[1]*(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  f_l_quad = (-0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.224744871391589*hamil[1]*(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += -0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.5*Fhat_nodal[0]-0.5*Fhat_nodal[1])*dx10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dx10; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 1.224744871391589*hamil[1]*(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*dv10; 
  f_l_quad = (-0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.224744871391589*hamil[1]*(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*dv10; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += 0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.5*Fhat_nodal[1]-0.5*Fhat_nodal[0])*dx10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dx10; 

  } 
  return 0.0;

} 
