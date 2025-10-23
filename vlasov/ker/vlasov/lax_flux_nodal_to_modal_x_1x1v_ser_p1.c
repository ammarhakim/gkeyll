#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x1v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *alpha_quad, const double *f_l, const double *f_r, 
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &Fhat_nodal[0]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = (-0.8660254037844386*f_l[3])-0.5*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 0.8660254037844386*f_r[3]-0.5*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = 0.8660254037844386*f_l[3]+0.5*f_l[2]+0.8660254037844386*f_l[1]+0.5*f_l[0]; 
  f_r_quad = (-0.8660254037844386*f_r[3])+0.5*f_r[2]-0.8660254037844386*f_r[1]+0.5*f_r[0]; 
  out[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  return 1.5*dx10*alpha_max;

} 
