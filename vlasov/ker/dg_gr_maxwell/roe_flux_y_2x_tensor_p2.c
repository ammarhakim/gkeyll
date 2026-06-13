#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double roe_flux_y_2x_tensor_p2(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double *out = &conf_flux_surf[24]; 
  double alpha_max = 0.0; 
  const double *flux_l_quad; 
  const double *A_minus_dQ_quad; 
  for (int i=0; i<6; ++i) {
    flux_l_quad = &flux_l[i*3]; 
    A_minus_dQ_quad = &A_minus_dQ[i*3]; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    out[0+i*3] = flux_l_quad[0] + A_minus_dQ_quad[0]; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    out[1+i*3] = flux_l_quad[1] + A_minus_dQ_quad[1]; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[2])); 
    out[2+i*3] = flux_l_quad[2] + A_minus_dQ_quad[2]; 
  }
  
  return 2.5*dx11*alpha_max;

} 
