#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x2v_ser_p2.h> 
GKYL_CU_DH double lax_flux_nodal_x_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &Fhat_nodal[0]; 
  double G_l[27]; 
  double G_r[27]; 
  for (int q = 0; q < 27; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 48; ++k) { 
    const int a = vst_2x2v_ser_p2_ph_x0_cmap[k]; 
    const double cfl = vst_2x2v_ser_p2_ph_x0_coefl[k]*f_l[k]; 
    const double cfr = vst_2x2v_ser_p2_ph_x0_coefr[k]*f_r[k]; 
    for (int j = 0; j < 9; ++j) { 
      G_l[j*3 + a] += vst_2x2v_ser_p2_ph_x0_V[j*20 + vst_2x2v_ser_p2_ph_x0_vlmap[k]]*cfl; 
      G_r[j*3 + a] += vst_2x2v_ser_p2_ph_x0_V[j*20 + vst_2x2v_ser_p2_ph_x0_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      const int n = i*9 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dx10*alpha_max;

} 
