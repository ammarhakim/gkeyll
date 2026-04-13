#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_y_3x_tensor_p1(const double *dxv, 
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double *out = &conf_flux_surf[32]; 
  double alpha_max = 0.0; 
  const double *em_no_J_l; 
  const double *em_no_J_r; 
  const double *flux_l_quad; 
  const double *flux_r_quad; 
  for (int i=0; i<6; ++i) {
    em_no_J_l = &field_no_J_con_l[i*8]; 
    em_no_J_r = &field_no_J_con_r[i*8]; 
    flux_l_quad = &flux_l[i*4]; 
    flux_r_quad = &flux_r[i*4]; 
    double U_l_quad = 0.0; 
    double U_r_quad = 0.0; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    U_l_quad = 0.6123724356957944*em_no_J_l[7]-0.6123724356957944*em_no_J_l[6]+0.3535533905932737*em_no_J_l[5]-0.6123724356957944*em_no_J_l[4]-0.3535533905932737*em_no_J_l[3]+0.6123724356957944*em_no_J_l[2]-0.3535533905932737*em_no_J_l[1]+0.3535533905932737*em_no_J_l[0]; 
    U_r_quad = -(0.6123724356957944*em_no_J_r[7])+0.6123724356957944*em_no_J_r[6]+0.3535533905932737*em_no_J_r[5]+0.6123724356957944*em_no_J_r[4]-0.3535533905932737*em_no_J_r[3]-0.6123724356957944*em_no_J_r[2]-0.3535533905932737*em_no_J_r[1]+0.3535533905932737*em_no_J_r[0]; 
    out[0+i*4] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*J_c[0]*(U_r_quad - U_l_quad)); 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    U_l_quad = -(0.6123724356957944*em_no_J_l[7])+0.6123724356957944*em_no_J_l[6]-0.3535533905932737*em_no_J_l[5]-0.6123724356957944*em_no_J_l[4]+0.3535533905932737*em_no_J_l[3]+0.6123724356957944*em_no_J_l[2]-0.3535533905932737*em_no_J_l[1]+0.3535533905932737*em_no_J_l[0]; 
    U_r_quad = 0.6123724356957944*em_no_J_r[7]-0.6123724356957944*em_no_J_r[6]-0.3535533905932737*em_no_J_r[5]+0.6123724356957944*em_no_J_r[4]+0.3535533905932737*em_no_J_r[3]-0.6123724356957944*em_no_J_r[2]-0.3535533905932737*em_no_J_r[1]+0.3535533905932737*em_no_J_r[0]; 
    out[1+i*4] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*J_c[1]*(U_r_quad - U_l_quad)); 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[2])); 
    U_l_quad = -(0.6123724356957944*em_no_J_l[7])-0.6123724356957944*em_no_J_l[6]-0.3535533905932737*em_no_J_l[5]+0.6123724356957944*em_no_J_l[4]-0.3535533905932737*em_no_J_l[3]+0.6123724356957944*em_no_J_l[2]+0.3535533905932737*em_no_J_l[1]+0.3535533905932737*em_no_J_l[0]; 
    U_r_quad = 0.6123724356957944*em_no_J_r[7]+0.6123724356957944*em_no_J_r[6]-0.3535533905932737*em_no_J_r[5]-0.6123724356957944*em_no_J_r[4]-0.3535533905932737*em_no_J_r[3]-0.6123724356957944*em_no_J_r[2]+0.3535533905932737*em_no_J_r[1]+0.3535533905932737*em_no_J_r[0]; 
    out[2+i*4] = 0.5*((flux_r_quad[2] + flux_l_quad[2]) - fabs(max_alpha_quad[2])*J_c[2]*(U_r_quad - U_l_quad)); 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[3])); 
    U_l_quad = 0.6123724356957944*em_no_J_l[7]+0.6123724356957944*em_no_J_l[6]+0.3535533905932737*em_no_J_l[5]+0.6123724356957944*em_no_J_l[4]+0.3535533905932737*em_no_J_l[3]+0.6123724356957944*em_no_J_l[2]+0.3535533905932737*em_no_J_l[1]+0.3535533905932737*em_no_J_l[0]; 
    U_r_quad = -(0.6123724356957944*em_no_J_r[7])-0.6123724356957944*em_no_J_r[6]+0.3535533905932737*em_no_J_r[5]-0.6123724356957944*em_no_J_r[4]+0.3535533905932737*em_no_J_r[3]-0.6123724356957944*em_no_J_r[2]+0.3535533905932737*em_no_J_r[1]+0.3535533905932737*em_no_J_r[0]; 
    out[3+i*4] = 0.5*((flux_r_quad[3] + flux_l_quad[3]) - fabs(max_alpha_quad[3])*J_c[3]*(U_r_quad - U_l_quad)); 
  }
  
  return 1.5*dx11*alpha_max;

} 
