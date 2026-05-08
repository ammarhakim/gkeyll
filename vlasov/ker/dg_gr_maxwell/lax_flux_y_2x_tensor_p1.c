#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_y_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double *out = &conf_flux_surf[16]; 
  double alpha_max = 0.0; 
  const double *em_l; 
  const double *em_r; 
  const double *flux_l_quad; 
  const double *flux_r_quad; 
  for (int i=0; i<6; ++i) {
    em_l = &field_con_l[i*4]; 
    em_r = &field_con_r[i*4]; 
    flux_l_quad = &flux_l[i*2]; 
    flux_r_quad = &flux_r[i*2]; 
    double Q_l_quad = 0.0; 
    double Q_r_quad = 0.0; 
    // At the theta pole there is no penalty because: dQ = 0; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    Q_l_quad = -(0.8660254037844386*em_l[3])+0.8660254037844386*em_l[2]-0.5*em_l[1]+0.5*em_l[0]; 
    Q_r_quad = 0.8660254037844386*em_r[3]-0.8660254037844386*em_r[2]-0.5*em_r[1]+0.5*em_r[0]; 
    if (theta_pole) { 
      out[0+i*2] = 0.5*(flux_r_quad[0] + flux_l_quad[0]); 
    }
    else {
      out[0+i*2] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*(Q_r_quad - Q_l_quad)); 
    }
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    Q_l_quad = 0.8660254037844386*em_l[3]+0.8660254037844386*em_l[2]+0.5*em_l[1]+0.5*em_l[0]; 
    Q_r_quad = -(0.8660254037844386*em_r[3])-0.8660254037844386*em_r[2]+0.5*em_r[1]+0.5*em_r[0]; 
    if (theta_pole) { 
      out[1+i*2] = 0.5*(flux_r_quad[1] + flux_l_quad[1]); 
    }
    else {
      out[1+i*2] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*(Q_r_quad - Q_l_quad)); 
    }
  }
  
  return 1.5*dx11*alpha_max;

} 
