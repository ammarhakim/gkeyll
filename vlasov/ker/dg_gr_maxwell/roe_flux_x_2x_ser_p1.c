#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double roe_flux_x_2x_ser_p1(const double *dxv, 
      const double *flux_l, const double *A_minus_dQ, const double *max_alpha_quad, 
      double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &conf_flux_surf[0]; 
  double alpha_max = 0.0; 
  const double *flux_l_quad; 
  const double *A_minus_dQ_quad; 
  for (int i=0; i<6; ++i) {
    flux_l_quad = &flux_l[i*2]; 
    A_minus_dQ_quad = &A_minus_dQ[i*2]; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    out[0+i*2] = flux_l_quad[0] + A_minus_dQ_quad[0]; 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    out[1+i*2] = flux_l_quad[1] + A_minus_dQ_quad[1]; 
  }
  
  return 1.5*dx10*alpha_max;

} 
