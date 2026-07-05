#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double lax_flux_nodal_z_3x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx12 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[64]; 
  double G_l[32]; 
  double G_r[32]; 
  for (int q = 0; q < 32; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 64; ++k) { 
    const int a = vst_3x3v_ser_p1_ph_x2_cmap[k]; 
    const double cfl = vst_3x3v_ser_p1_ph_x2_coefl[k]*f_l[k]; 
    const double cfr = vst_3x3v_ser_p1_ph_x2_coefr[k]*f_r[k]; 
    for (int j = 0; j < 8; ++j) { 
      G_l[j*4 + a] += vst_3x3v_ser_p1_ph_x2_V[j*16 + vst_3x3v_ser_p1_ph_x2_vlmap[k]]*cfl; 
      G_r[j*4 + a] += vst_3x3v_ser_p1_ph_x2_V[j*16 + vst_3x3v_ser_p1_ph_x2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 8; ++j) { 
      const int n = i*8 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 4; ++a) { 
        f_l_quad += vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]*G_l[j*4 + a]; 
        f_r_quad += vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]*G_r[j*4 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 1.5*dx12*alpha_max;

} 
