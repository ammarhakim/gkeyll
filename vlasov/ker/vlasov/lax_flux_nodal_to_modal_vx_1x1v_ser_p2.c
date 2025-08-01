#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel, 
  const double *alpha_quad, const double *f_l, const double *f_r, 
  double *lax_nodal_quad, double* GKYL_RESTRICT vel_flux_surf) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *Ghat = &vel_flux_surf[0]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = -(1.6675832968091233*f_l[7])+1.1857789576832338*f_l[6]+1.118033988749895*f_l[5]+0.6846098004178088*f_l[4]-1.291704467391079*f_l[3]+0.8660254037844386*f_l[2]-0.7457659219616816*f_l[1]+0.5*f_l[0]; 
  f_r_quad = -(1.6675832968091233*f_r[7])-1.1857789576832338*f_r[6]+1.118033988749895*f_r[5]+0.6846098004178088*f_r[4]+1.291704467391079*f_r[3]-0.8660254037844386*f_r[2]-0.7457659219616816*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = -(0.6583704599151832*f_l[7])-0.6324956225107456*f_l[6]+1.118033988749895*f_l[5]-0.3651715179178389*f_l[4]-0.5099715653772844*f_l[3]+0.8660254037844386*f_l[2]-0.29443222054963*f_l[1]+0.5*f_l[0]; 
  f_r_quad = -(0.6583704599151832*f_r[7])+0.6324956225107456*f_r[6]+1.118033988749895*f_r[5]-0.3651715179178389*f_r[4]+0.5099715653772844*f_r[3]-0.8660254037844386*f_r[2]-0.29443222054963*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[2])); 
  f_l_quad = 0.6583704599151832*f_l[7]-0.6324956225107456*f_l[6]+1.118033988749895*f_l[5]-0.3651715179178389*f_l[4]+0.5099715653772844*f_l[3]+0.8660254037844386*f_l[2]+0.29443222054963*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 0.6583704599151832*f_r[7]+0.6324956225107456*f_r[6]+1.118033988749895*f_r[5]-0.3651715179178389*f_r[4]-0.5099715653772844*f_r[3]-0.8660254037844386*f_r[2]+0.29443222054963*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[2] = 0.5*(alpha_quad[2]*(f_r_quad + f_l_quad) - fabs(alpha_quad[2])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[3])); 
  f_l_quad = 1.6675832968091233*f_l[7]+1.1857789576832338*f_l[6]+1.118033988749895*f_l[5]+0.6846098004178088*f_l[4]+1.291704467391079*f_l[3]+0.8660254037844386*f_l[2]+0.7457659219616816*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 1.6675832968091233*f_r[7]-1.1857789576832338*f_r[6]+1.118033988749895*f_r[5]+0.6846098004178088*f_r[4]-1.291704467391079*f_r[3]-0.8660254037844386*f_r[2]+0.7457659219616816*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[3] = 0.5*(alpha_quad[3]*(f_r_quad + f_l_quad) - fabs(alpha_quad[3])*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.24597051986529*lax_nodal_quad[3]+0.46113626132125757*lax_nodal_quad[2]+0.46113626132125757*lax_nodal_quad[1]+0.24597051986529*lax_nodal_quad[0]; 
  Ghat[1] = 0.3668728630454641*lax_nodal_quad[3]+0.271546746793545*lax_nodal_quad[2]-0.271546746793545*lax_nodal_quad[1]-0.3668728630454641*lax_nodal_quad[0]; 
  Ghat[2] = 0.3367876570272824*lax_nodal_quad[3]-0.3367876570272824*lax_nodal_quad[2]-0.3367876570272824*lax_nodal_quad[1]+0.3367876570272824*lax_nodal_quad[0]; 

  return 2.5*dv10*alpha_max;

} 
