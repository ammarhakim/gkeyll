#include <gkyl_nc_hamil_surf_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node2_x_1x1v_tensor_p1(const double *w, const double *dxv, 
   const double *poisson_tensor_conf, const double *hamil, const double *fl, 
   const double *fc, const double *fr, double* GKYL_RESTRICT alpha_max, double* GKYL_RESTRICT Ghat_l,
   double* GKYL_RESTRICT Ghat_r) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.612671688513684*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.612671688513684*hamil[2]+1.224744871391589*hamil[1])*dv10; 
  alpha_max[0] = fmax(alpha_max[0], fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6324956225107456*fl[5])-0.3651715179178389*fl[4]+0.5099715653772844*fl[3]+0.29443222054963*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.6324956225107456*fc[5]-0.3651715179178389*fc[4]-0.5099715653772844*fc[3]+0.29443222054963*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.6324956225107456*fc[5])-0.3651715179178389*fc[4]+0.5099715653772844*fc[3]+0.29443222054963*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.6324956225107456*fr[5]-0.3651715179178389*fr[4]-0.5099715653772844*fr[3]+0.29443222054963*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  double Ghat_l_nodal = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  double Ghat_r_nodal = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  Ghat_l[0] = 0.4611362613212576*Ghat_l_nodal; 
  Ghat_l[1] = 0.271546746793545*Ghat_l_nodal; 
  Ghat_l[2] = -0.3367876570272824*Ghat_l_nodal; 

  Ghat_r[0] = 0.4611362613212576*Ghat_r_nodal; 
  Ghat_r[1] = 0.271546746793545*Ghat_r_nodal; 
  Ghat_r[2] = -0.3367876570272824*Ghat_r_nodal; 

} 
