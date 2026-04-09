#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_2x_ser_p2(const double *dxv, 
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &conf_flux_surf[0]; 
  double alpha_max = 0.0; 
  const double *em_no_J_l; 
  const double *em_no_J_r; 
  const double *flux_l_quad; 
  const double *flux_r_quad; 
  for (int i=0; i<6; ++i) {
    em_no_J_l = &field_no_J_con_l[i*8]; 
    em_no_J_r = &field_no_J_con_r[i*8]; 
    flux_l_quad = &flux_l[i*3]; 
    flux_r_quad = &flux_r[i*3]; 
    double U_l_quad = 0.0; 
    double U_r_quad = 0.0; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    U_l_quad = 0.7745966692414833*em_no_J_l[7]-1.5*em_no_J_l[6]+0.4472135954999579*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]-1.1618950038622249*em_no_J_l[3]-0.6708203932499369*em_no_J_l[2]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]; 
    U_r_quad = -(0.7745966692414833*em_no_J_r[7])-1.5*em_no_J_r[6]+0.4472135954999579*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]+1.1618950038622249*em_no_J_r[3]-0.6708203932499369*em_no_J_r[2]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]; 
    out[0+i*3] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*J_c[0]*(U_r_quad - U_l_quad)); 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    U_l_quad = -(0.9682458365518543*em_no_J_l[7])-0.5590169943749475*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]; 
    U_r_quad = 0.9682458365518543*em_no_J_r[7]-0.5590169943749475*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]; 
    out[1+i*3] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*J_c[1]*(U_r_quad - U_l_quad)); 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[2])); 
    U_l_quad = 0.7745966692414833*em_no_J_l[7]+1.5*em_no_J_l[6]+0.4472135954999579*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]+1.1618950038622249*em_no_J_l[3]+0.6708203932499369*em_no_J_l[2]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]; 
    U_r_quad = -(0.7745966692414833*em_no_J_r[7])+1.5*em_no_J_r[6]+0.4472135954999579*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]-1.1618950038622249*em_no_J_r[3]+0.6708203932499369*em_no_J_r[2]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]; 
    out[2+i*3] = 0.5*((flux_r_quad[2] + flux_l_quad[2]) - fabs(max_alpha_quad[2])*J_c[2]*(U_r_quad - U_l_quad)); 
  }
  
  return 2.5*dx10*alpha_max;

} 
