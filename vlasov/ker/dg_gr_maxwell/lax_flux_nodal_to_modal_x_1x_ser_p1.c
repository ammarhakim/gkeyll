#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_nodal_to_modal_x_1x_ser_p1(const double *dxv, 
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
    em_no_J_l = &field_no_J_con_l[i*2]; 
    em_no_J_r = &field_no_J_con_r[i*2]; 
    flux_l_quad = &flux_l[i*1]; 
    flux_r_quad = &flux_r[i*1]; 
    double U_l_quad = 0.0; 
    double U_r_quad = 0.0; 
    alpha_max = fmax(alpha_max, J_c[0]*fabs(max_alpha_quad[0])); 
    U_l_quad = 1.224744871391589*em_no_J_l[1]+0.7071067811865475*em_no_J_l[0]; 
    U_r_quad = 0.7071067811865475*em_no_J_r[0]-1.224744871391589*em_no_J_r[1]; 
    out[0+i*1] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*J_c[0]*(U_r_quad - U_l_quad)); 
  }
  
  return 1.5*dx10*alpha_max;

} 
