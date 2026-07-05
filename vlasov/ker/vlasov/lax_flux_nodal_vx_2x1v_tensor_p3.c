#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_2x1v_tensor_p3.h> 
GKYL_CU_DH double lax_flux_nodal_vx_2x1v_tensor_p3(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[0]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  double jac_nodal[25]; 
  jac_nodal[0] = 1; 
  jac_nodal[1] = 1; 
  jac_nodal[2] = 1; 
  jac_nodal[3] = 1; 
  jac_nodal[4] = 1; 
  jac_nodal[5] = 1; 
  jac_nodal[6] = 1; 
  jac_nodal[7] = 1; 
  jac_nodal[8] = 1; 
  jac_nodal[9] = 1; 
  jac_nodal[10] = 1; 
  jac_nodal[11] = 1; 
  jac_nodal[12] = 1; 
  jac_nodal[13] = 1; 
  jac_nodal[14] = 1; 
  jac_nodal[15] = 1; 
  jac_nodal[16] = 1; 
  jac_nodal[17] = 1; 
  jac_nodal[18] = 1; 
  jac_nodal[19] = 1; 
  jac_nodal[20] = 1; 
  jac_nodal[21] = 1; 
  jac_nodal[22] = 1; 
  jac_nodal[23] = 1; 
  jac_nodal[24] = 1; 
  double G_l[16]; 
  double G_r[16]; 
  for (int q = 0; q < 16; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 64; ++k) { 
    const int a = vst_2x1v_tensor_p3_ph_v0_cmap[k]; 
    const double cfl = vst_2x1v_tensor_p3_ph_v0_coefl[k]*f_l[k]; 
    const double cfr = vst_2x1v_tensor_p3_ph_v0_coefr[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*16 + a] += vst_2x1v_tensor_p3_ph_v0_V[j*4 + vst_2x1v_tensor_p3_ph_v0_vlmap[k]]*cfl; 
      G_r[j*16 + a] += vst_2x1v_tensor_p3_ph_v0_V[j*4 + vst_2x1v_tensor_p3_ph_v0_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 16; ++a) { 
        f_l_quad += vst_2x1v_tensor_p3_ph_v0_Cm[i*16 + a]*G_l[j*16 + a]; 
        f_r_quad += vst_2x1v_tensor_p3_ph_v0_Cm[i*16 + a]*G_r[j*16 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_nodal[n]*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 3.5*dv10*alpha_max/(0.07635866179581297*jacob_vel_surf_vx[4]-0.2679416522233886*jacob_vel_surf_vx[3]+0.5333333333333333*jacob_vel_surf_vx[2]-0.8931583920000739*jacob_vel_surf_vx[1]+1.5514080490943136*jacob_vel_surf_vx[0]);

} 
