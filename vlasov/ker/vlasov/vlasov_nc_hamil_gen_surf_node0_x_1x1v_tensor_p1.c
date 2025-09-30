#include <gkyl_nc_hamil_surf_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_surf_node0_x_1x1v_tensor_p1(const double *w, const double *dxv, 
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
  alpha_l_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*dv10; 
  alpha_r_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*dv10; 
  alpha_max[0] = fmax(alpha_max[0], fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 1.185778957683234*fl[5]+0.6846098004178088*fl[4]-1.291704467391079*fl[3]-0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-1.185778957683234*fc[5])+0.6846098004178088*fc[4]+1.291704467391079*fc[3]-0.7457659219616816*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 1.185778957683234*fc[5]+0.6846098004178088*fc[4]-1.291704467391079*fc[3]-0.7457659219616816*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-1.185778957683234*fr[5])+0.6846098004178088*fr[4]+1.291704467391079*fr[3]-0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  double Ghat_l_nodal = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  double Ghat_r_nodal = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  Ghat_l[0] = 0.24597051986529*Ghat_l_nodal; 
  Ghat_l[1] = -0.3668728630454641*Ghat_l_nodal; 
  Ghat_l[2] = 0.3367876570272824*Ghat_l_nodal; 

  Ghat_r[0] = 0.24597051986529*Ghat_r_nodal; 
  Ghat_r[1] = -0.3668728630454641*Ghat_r_nodal; 
  Ghat_r[2] = 0.3367876570272824*Ghat_r_nodal; 

} 
