#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p1.h> 
GKYL_CU_DH double lax_flux_nodal_y_2x2v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[8]; 
  double G_l[8]; 
  double G_r[8]; 
  for (int q = 0; q < 8; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 16; ++k) { 
    const int a = vst_2x2v_ser_p1_ph_x1_cmap[k]; 
    const double cfl = vst_2x2v_ser_p1_ph_x1_coefl[k]*f_l[k]; 
    const double cfr = vst_2x2v_ser_p1_ph_x1_coefr[k]*f_r[k]; 
    for (int j = 0; j < 4; ++j) { 
      G_l[j*2 + a] += vst_2x2v_ser_p1_ph_x1_V[j*8 + vst_2x2v_ser_p1_ph_x1_vlmap[k]]*cfl; 
      G_r[j*2 + a] += vst_2x2v_ser_p1_ph_x1_V[j*8 + vst_2x2v_ser_p1_ph_x1_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      const int n = i*4 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 2; ++a) { 
        f_l_quad += vst_2x2v_ser_p1_ph_x1_Cm[i*2 + a]*G_l[j*2 + a]; 
        f_r_quad += vst_2x2v_ser_p1_ph_x1_Cm[i*2 + a]*G_r[j*2 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 1.5*dx11*alpha_max;

} 
