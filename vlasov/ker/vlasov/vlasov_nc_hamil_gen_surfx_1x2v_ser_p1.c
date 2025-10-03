#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  double Fhat_l_nodal[4] = {0.0}; 
  double Fhat_r_nodal[4] = {0.0}; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.6123724356957944*fl[7]+0.3535533905932737*fl[6]-0.6123724356957944*fl[5]-0.6123724356957944*fl[4]-0.3535533905932737*fl[3]-0.3535533905932737*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  f_cl_quad = (-0.6123724356957944*fc[7])+0.3535533905932737*fc[6]+0.6123724356957944*fc[5]+0.6123724356957944*fc[4]-0.3535533905932737*fc[3]-0.3535533905932737*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_cr_quad = 0.6123724356957944*fc[7]+0.3535533905932737*fc[6]-0.6123724356957944*fc[5]-0.6123724356957944*fc[4]-0.3535533905932737*fc[3]-0.3535533905932737*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_rl_quad = (-0.6123724356957944*fr[7])+0.3535533905932737*fr[6]+0.6123724356957944*fr[5]+0.6123724356957944*fr[4]-0.3535533905932737*fr[3]-0.3535533905932737*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  Fhat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[2]-0.8660254037844386*hamil[3])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6123724356957944*fl[7])-0.3535533905932737*fl[6]+0.6123724356957944*fl[5]-0.6123724356957944*fl[4]+0.3535533905932737*fl[3]-0.3535533905932737*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  f_cl_quad = 0.6123724356957944*fc[7]-0.3535533905932737*fc[6]-0.6123724356957944*fc[5]+0.6123724356957944*fc[4]+0.3535533905932737*fc[3]-0.3535533905932737*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_cr_quad = (-0.6123724356957944*fc[7])-0.3535533905932737*fc[6]+0.6123724356957944*fc[5]-0.6123724356957944*fc[4]+0.3535533905932737*fc[3]-0.3535533905932737*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_rl_quad = 0.6123724356957944*fr[7]-0.3535533905932737*fr[6]-0.6123724356957944*fr[5]+0.6123724356957944*fr[4]+0.3535533905932737*fr[3]-0.3535533905932737*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  Fhat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[1]-0.8660254037844386*hamil[3])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6123724356957944*fl[7])-0.3535533905932737*fl[6]-0.6123724356957944*fl[5]+0.6123724356957944*fl[4]-0.3535533905932737*fl[3]+0.3535533905932737*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  f_cl_quad = 0.6123724356957944*fc[7]-0.3535533905932737*fc[6]+0.6123724356957944*fc[5]-0.6123724356957944*fc[4]-0.3535533905932737*fc[3]+0.3535533905932737*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_cr_quad = (-0.6123724356957944*fc[7])-0.3535533905932737*fc[6]-0.6123724356957944*fc[5]+0.6123724356957944*fc[4]-0.3535533905932737*fc[3]+0.3535533905932737*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_rl_quad = 0.6123724356957944*fr[7]-0.3535533905932737*fr[6]+0.6123724356957944*fr[5]-0.6123724356957944*fr[4]-0.3535533905932737*fr[3]+0.3535533905932737*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  Fhat_l_nodal[2] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[2] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[2])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(0.8660254037844386*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.6123724356957944*fl[7]+0.3535533905932737*fl[6]+0.6123724356957944*fl[5]+0.6123724356957944*fl[4]+0.3535533905932737*fl[3]+0.3535533905932737*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  f_cl_quad = (-0.6123724356957944*fc[7])+0.3535533905932737*fc[6]-0.6123724356957944*fc[5]-0.6123724356957944*fc[4]+0.3535533905932737*fc[3]+0.3535533905932737*fc[2]-0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_cr_quad = 0.6123724356957944*fc[7]+0.3535533905932737*fc[6]+0.6123724356957944*fc[5]+0.6123724356957944*fc[4]+0.3535533905932737*fc[3]+0.3535533905932737*fc[2]+0.6123724356957944*fc[1]+0.3535533905932737*fc[0]; 
  f_rl_quad = (-0.6123724356957944*fr[7])+0.3535533905932737*fr[6]-0.6123724356957944*fr[5]-0.6123724356957944*fr[4]+0.3535533905932737*fr[3]+0.3535533905932737*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  Fhat_l_nodal[3] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[3] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  out[0] += ((-0.3535533905932737*Fhat_r_nodal[3])+0.3535533905932737*Fhat_l_nodal[3]-0.3535533905932737*Fhat_r_nodal[2]+0.3535533905932737*Fhat_l_nodal[2]-0.3535533905932737*Fhat_r_nodal[1]+0.3535533905932737*Fhat_l_nodal[1]-0.3535533905932737*Fhat_r_nodal[0]+0.3535533905932737*Fhat_l_nodal[0])*dx10; 
  out[1] += -0.6123724356957944*(Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0])*dx10; 
  out[2] += ((-0.3535533905932737*Fhat_r_nodal[3])+0.3535533905932737*Fhat_l_nodal[3]-0.3535533905932737*Fhat_r_nodal[2]+0.3535533905932737*(Fhat_l_nodal[2]+Fhat_r_nodal[1])-0.3535533905932737*Fhat_l_nodal[1]+0.3535533905932737*Fhat_r_nodal[0]-0.3535533905932737*Fhat_l_nodal[0])*dx10; 
  out[3] += ((-0.3535533905932737*Fhat_r_nodal[3])+0.3535533905932737*(Fhat_l_nodal[3]+Fhat_r_nodal[2])-0.3535533905932737*(Fhat_l_nodal[2]+Fhat_r_nodal[1])+0.3535533905932737*(Fhat_l_nodal[1]+Fhat_r_nodal[0])-0.3535533905932737*Fhat_l_nodal[0])*dx10; 
  out[4] += (0.6123724356957944*(Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0])-0.6123724356957944*(Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2]))*dx10; 
  out[5] += ((-0.6123724356957944*(Fhat_r_nodal[3]+Fhat_l_nodal[3]))+0.6123724356957944*(Fhat_r_nodal[2]+Fhat_l_nodal[2])-0.6123724356957944*(Fhat_r_nodal[1]+Fhat_l_nodal[1])+0.6123724356957944*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dx10; 
  out[6] += ((-0.3535533905932737*Fhat_r_nodal[3])+0.3535533905932737*(Fhat_l_nodal[3]+Fhat_r_nodal[2])-0.3535533905932737*Fhat_l_nodal[2]+0.3535533905932737*Fhat_r_nodal[1]-0.3535533905932737*(Fhat_l_nodal[1]+Fhat_r_nodal[0])+0.3535533905932737*Fhat_l_nodal[0])*dx10; 
  out[7] += ((-0.6123724356957944*(Fhat_r_nodal[3]+Fhat_l_nodal[3]))+0.6123724356957944*(Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.6123724356957944*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dx10; 

  return fabs(1.5*dx10*alpha_max);

} 
