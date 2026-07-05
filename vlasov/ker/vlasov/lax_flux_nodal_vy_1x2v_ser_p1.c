#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_ser_p1.h> 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv11 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[4]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  const double jac_tot = jacob_vel_surf_vx[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[3], jacob_vel_surf_r[3]); 
  double G_l[4]; 
  double G_r[4]; 
  for (int q = 0; q < 4; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 8; ++k) { 
    const int a = vst_1x2v_ser_p1_ph_v1_cmap[k]; 
    const double cfl = vst_1x2v_ser_p1_ph_v1_coefl[k]*f_l[k]; 
    const double cfr = vst_1x2v_ser_p1_ph_v1_coefr[k]*f_r[k]; 
    for (int j = 0; j < 2; ++j) { 
      G_l[j*2 + a] += vst_1x2v_ser_p1_ph_v1_V[j*4 + vst_1x2v_ser_p1_ph_v1_vlmap[k]]*cfl; 
      G_r[j*2 + a] += vst_1x2v_ser_p1_ph_v1_V[j*4 + vst_1x2v_ser_p1_ph_v1_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      const int n = i*2 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 2; ++a) { 
        f_l_quad += vst_1x2v_ser_p1_ph_v1_Cm[i*2 + a]*G_l[j*2 + a]; 
        f_r_quad += vst_1x2v_ser_p1_ph_v1_Cm[i*2 + a]*G_r[j*2 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_tot*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 1.5*dv11*alpha_max/jacob_vel_surf_min;

} 
