#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel_surf, 
  const double *alpha_quad, const double *f_l, const double *f_r, 
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[0]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = (-1.5*f_l[7])+0.7745966692414833*f_l[6]+1.118033988749895*f_l[5]+0.4472135954999579*f_l[4]-1.161895003862225*f_l[3]+0.8660254037844386*f_l[2]-0.6708203932499369*f_l[1]+0.5*f_l[0]; 
  f_r_quad = (-1.5*f_r[7])-0.7745966692414833*f_r[6]+1.118033988749895*f_r[5]+0.4472135954999579*f_r[4]+1.161895003862225*f_r[3]-0.8660254037844386*f_r[2]-0.6708203932499369*f_r[1]+0.5*f_r[0]; 
  out[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = (-0.9682458365518543*f_l[6])+1.118033988749895*f_l[5]-0.5590169943749475*f_l[4]+0.8660254037844386*f_l[2]+0.5*f_l[0]; 
  f_r_quad = 0.9682458365518543*f_r[6]+1.118033988749895*f_r[5]-0.5590169943749475*f_r[4]-0.8660254037844386*f_r[2]+0.5*f_r[0]; 
  out[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[2])); 
  f_l_quad = 1.5*f_l[7]+0.7745966692414833*f_l[6]+1.118033988749895*f_l[5]+0.4472135954999579*f_l[4]+1.161895003862225*f_l[3]+0.8660254037844386*f_l[2]+0.6708203932499369*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 1.5*f_r[7]-0.7745966692414833*f_r[6]+1.118033988749895*f_r[5]+0.4472135954999579*f_r[4]-1.161895003862225*f_r[3]-0.8660254037844386*f_r[2]+0.6708203932499369*f_r[1]+0.5*f_r[0]; 
  out[2] = 0.5*(alpha_quad[2]*(f_r_quad + f_l_quad) - fabs(alpha_quad[2])*(f_r_quad - f_l_quad)); 

  return 2.5*dv10*alpha_max;

} 
