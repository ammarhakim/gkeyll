#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*4 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*1 + j (inner fastest). 
static const double lax_flux_nodal_vx_1x1v_tensor_p3_Cm[20] = { 
  1.0, 
  -2.26544961484666, 
  -2.814395652992812, 
  1.138707206919686, 
  1.0, 
  -1.3461732752642077, 
  0.2502930888902482, 
  -0.9485956902878822, 
  1.0, 
  0.0, 
  1.9230769230769231, 
  0.0, 
  1.0, 
  1.3461732752642077, 
  0.2502930888902482, 
  0.9485956902878822, 
  1.0, 
  2.26544961484666, 
  -2.814395652992812, 
  -1.138707206919686, 
}; 
static const double lax_flux_nodal_vx_1x1v_tensor_p3_Vl[4] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
  -2.272727272727273, 
}; 
static const double lax_flux_nodal_vx_1x1v_tensor_p3_Vr[4] = { 
  1.0, 
  -2.5, 
  -3.8461538461538463, 
  2.272727272727273, 
}; 
static const double lax_flux_nodal_vx_1x1v_tensor_p3_coef[16] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
  -0.58206528843421, 
  -0.58206528843421, 
  0.169, 
  -0.40326666115611387, 
  -0.40326666115611387, 
  0.33839976359329804, 
  0.33839976359329804, 
  0.6776, 
}; 
static const int lax_flux_nodal_vx_1x1v_tensor_p3_cmap[16] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
  3, 
  0, 
  2, 
  3, 
  1, 
  3, 
  2, 
  3, 
}; 
static const int lax_flux_nodal_vx_1x1v_tensor_p3_vlmap[16] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  0, 
  3, 
  2, 
  1, 
  3, 
  2, 
  3, 
  3, 
}; 
static const int lax_flux_nodal_vx_1x1v_tensor_p3_vrmap[16] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  0, 
  3, 
  2, 
  1, 
  3, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_tensor_p3(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[0]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  double jac_nodal[5]; 
  jac_nodal[0] = 1; 
  jac_nodal[1] = 1; 
  jac_nodal[2] = 1; 
  jac_nodal[3] = 1; 
  jac_nodal[4] = 1; 
  double G_l[4]; 
  double G_r[4]; 
  for (int q = 0; q < 4; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 16; ++k) { 
    const int a = lax_flux_nodal_vx_1x1v_tensor_p3_cmap[k]; 
    const double cfl = lax_flux_nodal_vx_1x1v_tensor_p3_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vx_1x1v_tensor_p3_coef[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*4 + a] += lax_flux_nodal_vx_1x1v_tensor_p3_Vl[j*4 + lax_flux_nodal_vx_1x1v_tensor_p3_vlmap[k]]*cfl; 
      G_r[j*4 + a] += lax_flux_nodal_vx_1x1v_tensor_p3_Vr[j*4 + lax_flux_nodal_vx_1x1v_tensor_p3_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 4; ++a) { 
        f_l_quad += lax_flux_nodal_vx_1x1v_tensor_p3_Cm[i*4 + a]*G_l[j*4 + a]; 
        f_r_quad += lax_flux_nodal_vx_1x1v_tensor_p3_Cm[i*4 + a]*G_r[j*4 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_nodal[n]*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 3.5*dv10*alpha_max/(0.07635866179581297*jacob_vel_surf_vx[4]-0.2679416522233886*jacob_vel_surf_vx[3]+0.5333333333333333*jacob_vel_surf_vx[2]-0.8931583920000739*jacob_vel_surf_vx[1]+1.5514080490943136*jacob_vel_surf_vx[0]);

} 
