#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double ho_lax_flux_nodal_x_1x1v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *alpha_quad, const double *f_l, const double *f_r, 
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &Fhat_nodal[0]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = 1.185778957683234*f_l[7]-1.667583296809123*f_l[6]+0.6846098004178088*f_l[5]+1.118033988749895*f_l[4]-1.291704467391079*f_l[3]-0.7457659219616816*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = (-1.185778957683234*f_r[7])-1.667583296809123*f_r[6]+0.6846098004178088*f_r[5]+1.118033988749895*f_r[4]+1.291704467391079*f_r[3]-0.7457659219616816*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = (-0.6324956225107456*f_l[7])-0.6583704599151832*f_l[6]-0.3651715179178389*f_l[5]+1.118033988749895*f_l[4]-0.5099715653772844*f_l[3]-0.29443222054963*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 0.6324956225107456*f_r[7]-0.6583704599151832*f_r[6]-0.3651715179178389*f_r[5]+1.118033988749895*f_r[4]+0.5099715653772844*f_r[3]-0.29443222054963*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[2])); 
  f_l_quad = (-0.6324956225107456*f_l[7])+0.6583704599151832*f_l[6]-0.3651715179178389*f_l[5]+1.118033988749895*f_l[4]+0.5099715653772844*f_l[3]+0.29443222054963*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 0.6324956225107456*f_r[7]+0.6583704599151832*f_r[6]-0.3651715179178389*f_r[5]+1.118033988749895*f_r[4]-0.5099715653772844*f_r[3]+0.29443222054963*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[2] = 0.5*(alpha_quad[2]*(f_r_quad + f_l_quad) - fabs(alpha_quad[2])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[3])); 
  f_l_quad = 1.185778957683234*f_l[7]+1.667583296809123*f_l[6]+0.6846098004178088*f_l[5]+1.118033988749895*f_l[4]+1.291704467391079*f_l[3]+0.7457659219616816*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = (-1.185778957683234*f_r[7])+1.667583296809123*f_r[6]+0.6846098004178088*f_r[5]+1.118033988749895*f_r[4]-1.291704467391079*f_r[3]+0.7457659219616816*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[3] = 0.5*(alpha_quad[3]*(f_r_quad + f_l_quad) - fabs(alpha_quad[3])*(f_r_quad - f_l_quad)); 

  return 2.5*dx10*alpha_max;

} 
