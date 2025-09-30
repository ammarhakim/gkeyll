#include <gkyl_nc_hamil_surf_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node0_x_1x1v_tensor_p1(const double *w, const double *dxv, const int edge,
   const double *poisson_tensor_conf, const double *hamil, const double *fl, 
   const double *fr, double* GKYL_RESTRICT Ghat) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 

  double alpha_quad = 0.0; 
  if (edge == -1) { 
    alpha_quad = (1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*dv10; 
  } else { 

    alpha_quad = (0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*(1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*dv10; 
  } 
 
  f_l_quad = 1.185778957683234*fl[5]+0.6846098004178088*fl[4]-1.291704467391079*fl[3]-0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_r_quad = (-1.185778957683234*fr[5])+0.6846098004178088*fr[4]+1.291704467391079*fr[3]-0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  double Ghat_nodal = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 
  Ghat[0] = 0.24597051986529*Ghat_nodal; 
  Ghat[1] = -0.3668728630454641*Ghat_nodal; 
  Ghat[2] = 0.3367876570272824*Ghat_nodal; 

} 
