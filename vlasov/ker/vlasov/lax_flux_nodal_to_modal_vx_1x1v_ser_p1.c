#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_vx_1x1v_ser_p1(const double *dxv, const double *alpha_quad, 
  const double *f_l, const double *f_r, double *lax_nodal_quad, 
  double* GKYL_RESTRICT vel_flux_surf) 
{ 
  const double dv10 = 2.0/dxv[1]; 

  double *Ghat = &vel_flux_surf[0]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = -(0.8660254037844386*f_l[3])+0.8660254037844386*f_l[2]-0.5*f_l[1]+0.5*f_l[0]; 
  f_r_quad = 0.8660254037844386*f_r[3]-0.8660254037844386*f_r[2]-0.5*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = 0.8660254037844386*f_l[3]+0.8660254037844386*f_l[2]+0.5*f_l[1]+0.5*f_l[0]; 
  f_r_quad = -(0.8660254037844386*f_r[3])-0.8660254037844386*f_r[2]+0.5*f_r[1]+0.5*f_r[0]; 
  lax_nodal_quad[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.7071067811865475*lax_nodal_quad[1]+0.7071067811865475*lax_nodal_quad[0]; 
  Ghat[1] = 0.7071067811865475*lax_nodal_quad[1]-0.7071067811865475*lax_nodal_quad[0]; 

  return 1.5*dv10*alpha_max;

} 
