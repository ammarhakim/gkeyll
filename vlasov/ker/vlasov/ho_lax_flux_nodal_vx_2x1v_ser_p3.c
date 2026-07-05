#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x1v_ser_p3.h> 
GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  double G_l[12]; 
  double G_r[12]; 
  for (int q = 0; q < 12; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 32; ++k) { 
    const int a = vst_2x1v_ser_p3_ho_ph_v0_cmap[k]; 
    const double cfl = vst_2x1v_ser_p3_ho_ph_v0_coefl[k]*f_l[k]; 
    const double cfr = vst_2x1v_ser_p3_ho_ph_v0_coefr[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*12 + a] += vst_2x1v_ser_p3_ho_ph_v0_V[j*4 + vst_2x1v_ser_p3_ho_ph_v0_vlmap[k]]*cfl; 
      G_r[j*12 + a] += vst_2x1v_ser_p3_ho_ph_v0_V[j*4 + vst_2x1v_ser_p3_ho_ph_v0_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 12; ++a) { 
        f_l_quad += vst_2x1v_ser_p3_ho_ph_v0_Cm[i*12 + a]*G_l[j*12 + a]; 
        f_r_quad += vst_2x1v_ser_p3_ho_ph_v0_Cm[i*12 + a]*G_r[j*12 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 3.5*dv10*alpha_max/jacob_vel_surf_min;

} 
