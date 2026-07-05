#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv12 = 2.0/dxv[3]; 

  double *out = &Fhat_nodal[128]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf_r[4]; 
  const double jac_tot = jacob_vel_surf_vx[0]*jacob_vel_surf_vy[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[8], jacob_vel_surf_r[8]); 
  double G_l[48]; 
  double G_r[48]; 
  for (int q = 0; q < 48; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 48; ++k) { 
    const int a = vst_1x3v_ser_p2_ho_ph_v2_cmap[k]; 
    const double cfl = vst_1x3v_ser_p2_ho_ph_v2_coefl[k]*f_l[k]; 
    const double cfr = vst_1x3v_ser_p2_ho_ph_v2_coefr[k]*f_r[k]; 
    for (int j = 0; j < 16; ++j) { 
      G_l[j*3 + a] += vst_1x3v_ser_p2_ho_ph_v2_V[j*20 + vst_1x3v_ser_p2_ho_ph_v2_vlmap[k]]*cfl; 
      G_r[j*3 + a] += vst_1x3v_ser_p2_ho_ph_v2_V[j*20 + vst_1x3v_ser_p2_ho_ph_v2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      const int n = i*16 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_tot*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv12*alpha_max/jacob_vel_surf_min;

} 
