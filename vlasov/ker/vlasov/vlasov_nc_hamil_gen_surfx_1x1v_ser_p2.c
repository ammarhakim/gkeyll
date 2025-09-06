#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  double Ghat_r[3] = {0.0}; 
  double Ghat_l[3] = {0.0}; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  double Ghat_l_nodal[4] = {0.0}; 
  double Ghat_r_nodal[4] = {0.0}; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = (1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_r_quad = (1.224744871391589*hamil[1]-4.084728180770504*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 1.185778957683234*fl[7]-1.667583296809123*fl[6]+0.6846098004178088*fl[5]+1.118033988749895*fl[4]-1.291704467391079*fl[3]-0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-1.185778957683234*fc[7])-1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]+1.291704467391079*fc[3]-0.7457659219616816*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 1.185778957683234*fc[7]-1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]-1.291704467391079*fc[3]-0.7457659219616816*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-1.185778957683234*fr[7])-1.667583296809123*fr[6]+0.6846098004178088*fr[5]+1.118033988749895*fr[4]+1.291704467391079*fr[3]-0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (1.224744871391589*hamil[1]-1.612671688513684*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_r_quad = (1.224744871391589*hamil[1]-1.612671688513684*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6324956225107456*fl[7])-0.6583704599151832*fl[6]-0.3651715179178389*fl[5]+1.118033988749895*fl[4]-0.5099715653772844*fl[3]-0.29443222054963*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.6324956225107456*fc[7]-0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]+0.5099715653772844*fc[3]-0.29443222054963*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.6324956225107456*fc[7])-0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]-0.5099715653772844*fc[3]-0.29443222054963*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.6324956225107456*fr[7]-0.6583704599151832*fr[6]-0.3651715179178389*fr[5]+1.118033988749895*fr[4]+0.5099715653772844*fr[3]-0.29443222054963*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (1.612671688513684*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_r_quad = (1.612671688513684*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6324956225107456*fl[7])+0.6583704599151832*fl[6]-0.3651715179178389*fl[5]+1.118033988749895*fl[4]+0.5099715653772844*fl[3]+0.29443222054963*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.6324956225107456*fc[7]+0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]-0.5099715653772844*fc[3]+0.29443222054963*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.6324956225107456*fc[7])+0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]+0.5099715653772844*fc[3]+0.29443222054963*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.6324956225107456*fr[7]+0.6583704599151832*fr[6]-0.3651715179178389*fr[5]+1.118033988749895*fr[4]-0.5099715653772844*fr[3]+0.29443222054963*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[2] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[2] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (4.084728180770504*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_r_quad = (4.084728180770504*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 1.185778957683234*fl[7]+1.667583296809123*fl[6]+0.6846098004178088*fl[5]+1.118033988749895*fl[4]+1.291704467391079*fl[3]+0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-1.185778957683234*fc[7])+1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]-1.291704467391079*fc[3]+0.7457659219616816*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 1.185778957683234*fc[7]+1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]+1.291704467391079*fc[3]+0.7457659219616816*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-1.185778957683234*fr[7])+1.667583296809123*fr[6]+0.6846098004178088*fr[5]+1.118033988749895*fr[4]-1.291704467391079*fr[3]+0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[3] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[3] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  Ghat_l[0] = 0.24597051986529*Ghat_l_nodal[3]+0.4611362613212576*Ghat_l_nodal[2]+0.4611362613212576*Ghat_l_nodal[1]+0.24597051986529*Ghat_l_nodal[0]; 
  Ghat_l[1] = 0.3668728630454641*Ghat_l_nodal[3]+0.271546746793545*Ghat_l_nodal[2]-0.271546746793545*Ghat_l_nodal[1]-0.3668728630454641*Ghat_l_nodal[0]; 
  Ghat_l[2] = 0.3367876570272824*Ghat_l_nodal[3]-0.3367876570272824*Ghat_l_nodal[2]-0.3367876570272824*Ghat_l_nodal[1]+0.3367876570272824*Ghat_l_nodal[0]; 

  Ghat_r[0] = 0.24597051986529*Ghat_r_nodal[3]+0.4611362613212576*Ghat_r_nodal[2]+0.4611362613212576*Ghat_r_nodal[1]+0.24597051986529*Ghat_r_nodal[0]; 
  Ghat_r[1] = 0.3668728630454641*Ghat_r_nodal[3]+0.271546746793545*Ghat_r_nodal[2]-0.271546746793545*Ghat_r_nodal[1]-0.3668728630454641*Ghat_r_nodal[0]; 
  Ghat_r[2] = 0.3367876570272824*Ghat_r_nodal[3]-0.3367876570272824*Ghat_r_nodal[2]-0.3367876570272824*Ghat_r_nodal[1]+0.3367876570272824*Ghat_r_nodal[0]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 

  return fabs(2.5*dv10*dx10*alpha_max);

} 
