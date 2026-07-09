#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_ho_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Fhat_l_nodal[4] = {0.0}; 
  double Fhat_r_nodal[4] = {0.0}; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = 5.002749890427371*hamil[7]+1.936491673103709*hamil[6]-2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = (-5.002749890427371*hamil[7])+1.936491673103709*hamil[6]-2.888338995926613*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 1.185778957683234*fl[7]-1.667583296809123*fl[6]+0.6846098004178088*fl[5]+1.118033988749895*fl[4]-1.291704467391079*fl[3]-0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-1.185778957683234*fc[7])-1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]+1.291704467391079*fc[3]-0.7457659219616816*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 1.185778957683234*fc[7]-1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]-1.291704467391079*fc[3]-0.7457659219616816*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-1.185778957683234*fr[7])-1.667583296809123*fr[6]+0.6846098004178088*fr[5]+1.118033988749895*fr[4]+1.291704467391079*fr[3]-0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 1.97511137974555*hamil[7]+1.936491673103709*hamil[6]-1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = (-1.97511137974555*hamil[7])+1.936491673103709*hamil[6]-1.140331086775586*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6324956225107456*fl[7])-0.6583704599151832*fl[6]-0.3651715179178389*fl[5]+1.118033988749895*fl[4]-0.5099715653772844*fl[3]-0.29443222054963*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.6324956225107456*fc[7]-0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]+0.5099715653772844*fc[3]-0.29443222054963*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.6324956225107456*fc[7])-0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]-0.5099715653772844*fc[3]-0.29443222054963*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.6324956225107456*fr[7]-0.6583704599151832*fr[6]-0.3651715179178389*fr[5]+1.118033988749895*fr[4]+0.5099715653772844*fr[3]-0.29443222054963*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-1.97511137974555*hamil[7])+1.936491673103709*hamil[6]+1.140331086775586*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = 1.97511137974555*hamil[7]+1.936491673103709*hamil[6]+1.140331086775586*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.6324956225107456*fl[7])+0.6583704599151832*fl[6]-0.3651715179178389*fl[5]+1.118033988749895*fl[4]+0.5099715653772844*fl[3]+0.29443222054963*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.6324956225107456*fc[7]+0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]-0.5099715653772844*fc[3]+0.29443222054963*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.6324956225107456*fc[7])+0.6583704599151832*fc[6]-0.3651715179178389*fc[5]+1.118033988749895*fc[4]+0.5099715653772844*fc[3]+0.29443222054963*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.6324956225107456*fr[7]+0.6583704599151832*fr[6]-0.3651715179178389*fr[5]+1.118033988749895*fr[4]-0.5099715653772844*fr[3]+0.29443222054963*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[2] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[2] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-5.002749890427371*hamil[7])+1.936491673103709*hamil[6]+2.888338995926613*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = 5.002749890427371*hamil[7]+1.936491673103709*hamil[6]+2.888338995926613*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 1.185778957683234*fl[7]+1.667583296809123*fl[6]+0.6846098004178088*fl[5]+1.118033988749895*fl[4]+1.291704467391079*fl[3]+0.7457659219616816*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-1.185778957683234*fc[7])+1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]-1.291704467391079*fc[3]+0.7457659219616816*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 1.185778957683234*fc[7]+1.667583296809123*fc[6]+0.6846098004178088*fc[5]+1.118033988749895*fc[4]+1.291704467391079*fc[3]+0.7457659219616816*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-1.185778957683234*fr[7])+1.667583296809123*fr[6]+0.6846098004178088*fr[5]+1.118033988749895*fr[4]-1.291704467391079*fr[3]+0.7457659219616816*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[3] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[3] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  out[0] += ((-0.1739274225687269*Fhat_r_nodal[3])+0.1739274225687269*Fhat_l_nodal[3]-0.3260725774312731*Fhat_r_nodal[2]+0.3260725774312731*Fhat_l_nodal[2]-0.3260725774312731*Fhat_r_nodal[1]+0.3260725774312731*Fhat_l_nodal[1]-0.1739274225687269*Fhat_r_nodal[0]+0.1739274225687269*Fhat_l_nodal[0])*dv10*dx10; 
  out[1] += ((-0.3012511327185366*(Fhat_r_nodal[3]+Fhat_l_nodal[3]))-0.5647742710659018*(Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.3012511327185366*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv10*dx10; 
  out[2] += ((-0.2594182892927712*Fhat_r_nodal[3])+0.2594182892927712*Fhat_l_nodal[3]-0.192012546066862*Fhat_r_nodal[2]+0.192012546066862*(Fhat_l_nodal[2]+Fhat_r_nodal[1])-0.192012546066862*Fhat_l_nodal[1]+0.2594182892927712*Fhat_r_nodal[0]-0.2594182892927712*Fhat_l_nodal[0])*dv10*dx10; 
  out[3] += ((-0.449325657467681*(Fhat_r_nodal[3]+Fhat_l_nodal[3]))-0.3325754854784645*(Fhat_r_nodal[2]+Fhat_l_nodal[2])+0.3325754854784645*(Fhat_r_nodal[1]+Fhat_l_nodal[1])+0.449325657467681*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv10*dx10; 
  out[4] += ((-0.3889135400150045*Fhat_r_nodal[3])+0.3889135400150045*Fhat_l_nodal[3]-0.7291204487348905*Fhat_r_nodal[2]+0.7291204487348905*Fhat_l_nodal[2]-0.7291204487348905*Fhat_r_nodal[1]+0.7291204487348905*Fhat_l_nodal[1]-0.3889135400150045*Fhat_r_nodal[0]+0.3889135400150045*Fhat_l_nodal[0])*dv10*dx10; 
  out[5] += ((-0.2381448361039205*Fhat_r_nodal[3])+0.2381448361039205*(Fhat_l_nodal[3]+Fhat_r_nodal[2])-0.2381448361039205*Fhat_l_nodal[2]+0.2381448361039205*Fhat_r_nodal[1]-0.2381448361039205*(Fhat_l_nodal[1]+Fhat_r_nodal[0])+0.2381448361039205*Fhat_l_nodal[0])*dv10*dx10; 
  out[6] += ((-0.5800769294653423*Fhat_r_nodal[3])+0.5800769294653423*Fhat_l_nodal[3]-0.4293531055383132*Fhat_r_nodal[2]+0.4293531055383132*(Fhat_l_nodal[2]+Fhat_r_nodal[1])-0.4293531055383132*Fhat_l_nodal[1]+0.5800769294653423*Fhat_r_nodal[0]-0.5800769294653423*Fhat_l_nodal[0])*dv10*dx10; 
  out[7] += ((-0.4124789556921534*(Fhat_r_nodal[3]+Fhat_l_nodal[3]))+0.4124789556921534*(Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.4124789556921534*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv10*dx10; 

  return fabs(2.5*dv10*dx10*alpha_max);

} 
