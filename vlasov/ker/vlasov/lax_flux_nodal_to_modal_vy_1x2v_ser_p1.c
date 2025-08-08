#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *alpha_quad, const double *f_l, const double *f_r, 
  double *lax_nodal_quad, double* GKYL_RESTRICT vel_flux_surf) 
{ 
  double dv11 = 2.0/dxv[2]; 

  double *Ghat = &vel_flux_surf[4]; 
  double alpha_max = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_max = fmax(alpha_max, fabs(alpha_quad[0])); 
  f_l_quad = 0.6123724356957944*f_l[7]-0.6123724356957944*f_l[6]-0.6123724356957944*f_l[5]+0.3535533905932737*f_l[4]+0.6123724356957944*f_l[3]-0.3535533905932737*f_l[2]-0.3535533905932737*f_l[1]+0.3535533905932737*f_l[0]; 
  f_r_quad = -(0.6123724356957944*f_r[7])+0.6123724356957944*f_r[6]+0.6123724356957944*f_r[5]+0.3535533905932737*f_r[4]-0.6123724356957944*f_r[3]-0.3535533905932737*f_r[2]-0.3535533905932737*f_r[1]+0.3535533905932737*f_r[0]; 
  lax_nodal_quad[0] = 0.5*(alpha_quad[0]*(f_r_quad + f_l_quad) - fabs(alpha_quad[0])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[1])); 
  f_l_quad = -(0.6123724356957944*f_l[7])+0.6123724356957944*f_l[6]-0.6123724356957944*f_l[5]-0.3535533905932737*f_l[4]+0.6123724356957944*f_l[3]+0.3535533905932737*f_l[2]-0.3535533905932737*f_l[1]+0.3535533905932737*f_l[0]; 
  f_r_quad = 0.6123724356957944*f_r[7]-0.6123724356957944*f_r[6]+0.6123724356957944*f_r[5]-0.3535533905932737*f_r[4]-0.6123724356957944*f_r[3]+0.3535533905932737*f_r[2]-0.3535533905932737*f_r[1]+0.3535533905932737*f_r[0]; 
  lax_nodal_quad[1] = 0.5*(alpha_quad[1]*(f_r_quad + f_l_quad) - fabs(alpha_quad[1])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[2])); 
  f_l_quad = -(0.6123724356957944*f_l[7])-0.6123724356957944*f_l[6]+0.6123724356957944*f_l[5]-0.3535533905932737*f_l[4]+0.6123724356957944*f_l[3]-0.3535533905932737*f_l[2]+0.3535533905932737*f_l[1]+0.3535533905932737*f_l[0]; 
  f_r_quad = 0.6123724356957944*f_r[7]+0.6123724356957944*f_r[6]-0.6123724356957944*f_r[5]-0.3535533905932737*f_r[4]-0.6123724356957944*f_r[3]-0.3535533905932737*f_r[2]+0.3535533905932737*f_r[1]+0.3535533905932737*f_r[0]; 
  lax_nodal_quad[2] = 0.5*(alpha_quad[2]*(f_r_quad + f_l_quad) - fabs(alpha_quad[2])*(f_r_quad - f_l_quad)); 

  alpha_max = fmax(alpha_max, fabs(alpha_quad[3])); 
  f_l_quad = 0.6123724356957944*f_l[7]+0.6123724356957944*f_l[6]+0.6123724356957944*f_l[5]+0.3535533905932737*f_l[4]+0.6123724356957944*f_l[3]+0.3535533905932737*f_l[2]+0.3535533905932737*f_l[1]+0.3535533905932737*f_l[0]; 
  f_r_quad = -(0.6123724356957944*f_r[7])-0.6123724356957944*f_r[6]-0.6123724356957944*f_r[5]+0.3535533905932737*f_r[4]-0.6123724356957944*f_r[3]+0.3535533905932737*f_r[2]+0.3535533905932737*f_r[1]+0.3535533905932737*f_r[0]; 
  lax_nodal_quad[3] = 0.5*(alpha_quad[3]*(f_r_quad + f_l_quad) - fabs(alpha_quad[3])*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.5*lax_nodal_quad[3]+0.5*lax_nodal_quad[2]+0.5*lax_nodal_quad[1]+0.5*lax_nodal_quad[0]; 
  Ghat[1] = 0.5*lax_nodal_quad[3]+0.5*lax_nodal_quad[2]-0.5*lax_nodal_quad[1]-0.5*lax_nodal_quad[0]; 
  Ghat[2] = 0.5*lax_nodal_quad[3]-0.5*lax_nodal_quad[2]+0.5*lax_nodal_quad[1]-0.5*lax_nodal_quad[0]; 
  Ghat[3] = 0.5*lax_nodal_quad[3]-0.5*lax_nodal_quad[2]-0.5*lax_nodal_quad[1]+0.5*lax_nodal_quad[0]; 

  return 1.5*dv11*alpha_max;

} 
