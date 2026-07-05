#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_tables_1x2v_tensor_p2.h> 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv11 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[16]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf_r[4]; 
  double jac_nodal[16]; 
  jac_nodal[0] = 1.0*jacob_vel_surf_vx[0]; 
  jac_nodal[1] = 1.0*jacob_vel_surf_vx[1]; 
  jac_nodal[2] = 1.0*jacob_vel_surf_vx[2]; 
  jac_nodal[3] = 1.0*jacob_vel_surf_vx[3]; 
  jac_nodal[4] = 1.0*jacob_vel_surf_vx[0]; 
  jac_nodal[5] = 1.0*jacob_vel_surf_vx[1]; 
  jac_nodal[6] = 1.0*jacob_vel_surf_vx[2]; 
  jac_nodal[7] = 1.0*jacob_vel_surf_vx[3]; 
  jac_nodal[8] = 1.0*jacob_vel_surf_vx[0]; 
  jac_nodal[9] = 1.0*jacob_vel_surf_vx[1]; 
  jac_nodal[10] = 1.0*jacob_vel_surf_vx[2]; 
  jac_nodal[11] = 1.0*jacob_vel_surf_vx[3]; 
  jac_nodal[12] = 1.0*jacob_vel_surf_vx[0]; 
  jac_nodal[13] = 1.0*jacob_vel_surf_vx[1]; 
  jac_nodal[14] = 1.0*jacob_vel_surf_vx[2]; 
  jac_nodal[15] = 1.0*jacob_vel_surf_vx[3]; 
  double G_l[12]; 
  double G_r[12]; 
  for (int q = 0; q < 12; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 27; ++k) { 
    const int a = vst_1x2v_tensor_p2_ph_v1_cmap[k]; 
    const double cfl = vst_1x2v_tensor_p2_ph_v1_coefl[k]*f_l[k]; 
    const double cfr = vst_1x2v_tensor_p2_ph_v1_coefr[k]*f_r[k]; 
    for (int j = 0; j < 4; ++j) { 
      G_l[j*3 + a] += vst_1x2v_tensor_p2_ph_v1_V[j*9 + vst_1x2v_tensor_p2_ph_v1_vlmap[k]]*cfl; 
      G_r[j*3 + a] += vst_1x2v_tensor_p2_ph_v1_V[j*9 + vst_1x2v_tensor_p2_ph_v1_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      const int n = i*4 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += vst_1x2v_tensor_p2_ph_v1_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += vst_1x2v_tensor_p2_ph_v1_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_nodal[n]*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv11*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vy[3])+0.40076152031165013*jacob_vel_surf_vy[2]-0.8136324494869249*jacob_vel_surf_vy[1]+1.5267881254572662*jacob_vel_surf_vy[0]);

} 
