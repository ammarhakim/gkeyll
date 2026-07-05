#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*4 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*1 + j (inner fastest). 
static const double lax_flux_nodal_vx_1x1v_ser_p3_Cm[16] = { 
  1.0, 
  -2.1528407789851314, 
  -2.355129310456592, 
  0.6926067839891049, 
  1.0, 
  -0.8499526089621406, 
  1.2562282115554928, 
  -0.9357454538020444, 
  1.0, 
  0.8499526089621406, 
  1.2562282115554928, 
  0.9357454538020444, 
  1.0, 
  2.1528407789851314, 
  -2.355129310456592, 
  -0.6926067839891049, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p3_Vl[4] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
  -2.272727272727273, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p3_Vr[4] = { 
  1.0, 
  -2.5, 
  -3.8461538461538463, 
  2.272727272727273, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p3_coef[12] = { 
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
  -0.40326666115611387, 
  -0.40326666115611387, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p3_cmap[12] = { 
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
  3, 
  1, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p3_vlmap[12] = { 
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
  1, 
  3, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p3_vrmap[12] = { 
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
  1, 
  3, 
}; 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_ser_p3(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  double G_l[4]; 
  double G_r[4]; 
  for (int q = 0; q < 4; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 12; ++k) { 
    const int a = lax_flux_nodal_vx_1x1v_ser_p3_cmap[k]; 
    const double cfl = lax_flux_nodal_vx_1x1v_ser_p3_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vx_1x1v_ser_p3_coef[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*4 + a] += lax_flux_nodal_vx_1x1v_ser_p3_Vl[j*4 + lax_flux_nodal_vx_1x1v_ser_p3_vlmap[k]]*cfl; 
      G_r[j*4 + a] += lax_flux_nodal_vx_1x1v_ser_p3_Vr[j*4 + lax_flux_nodal_vx_1x1v_ser_p3_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 4; ++a) { 
        f_l_quad += lax_flux_nodal_vx_1x1v_ser_p3_Cm[i*4 + a]*G_l[j*4 + a]; 
        f_r_quad += lax_flux_nodal_vx_1x1v_ser_p3_Cm[i*4 + a]*G_r[j*4 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 3.5*dv10*alpha_max/jacob_vel_surf_min;

} 
