#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*3 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*1 + j (inner fastest). 
static const double ho_lax_flux_nodal_vx_1x1v_ser_p2_Cm[12] = { 
  1.0, 
  -2.1528407789851314, 
  -2.355129310456592, 
  1.0, 
  -0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  2.1528407789851314, 
  -2.355129310456592, 
}; 
static const double ho_lax_flux_nodal_vx_1x1v_ser_p2_Vl[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
static const double ho_lax_flux_nodal_vx_1x1v_ser_p2_Vr[3] = { 
  1.0, 
  -2.5, 
  -3.8461538461538463, 
}; 
static const double ho_lax_flux_nodal_vx_1x1v_ser_p2_coef[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
static const int ho_lax_flux_nodal_vx_1x1v_ser_p2_cmap[8] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
}; 
static const int ho_lax_flux_nodal_vx_1x1v_ser_p2_vlmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
static const int ho_lax_flux_nodal_vx_1x1v_ser_p2_vrmap[8] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x1v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  double G_l[3]; 
  double G_r[3]; 
  for (int q = 0; q < 3; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 8; ++k) { 
    const int a = ho_lax_flux_nodal_vx_1x1v_ser_p2_cmap[k]; 
    const double cfl = ho_lax_flux_nodal_vx_1x1v_ser_p2_coef[k]*f_l[k]; 
    const double cfr = ho_lax_flux_nodal_vx_1x1v_ser_p2_coef[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*3 + a] += ho_lax_flux_nodal_vx_1x1v_ser_p2_Vl[j*3 + ho_lax_flux_nodal_vx_1x1v_ser_p2_vlmap[k]]*cfl; 
      G_r[j*3 + a] += ho_lax_flux_nodal_vx_1x1v_ser_p2_Vr[j*3 + ho_lax_flux_nodal_vx_1x1v_ser_p2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += ho_lax_flux_nodal_vx_1x1v_ser_p2_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += ho_lax_flux_nodal_vx_1x1v_ser_p2_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv10*alpha_max/jacob_vel_surf_min;

} 
