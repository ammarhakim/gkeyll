#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*1 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*3 + j (inner fastest). 
static const double lax_flux_nodal_x_1x1v_ser_p2_Cm[1] = { 
  1.0, 
}; 
static const double lax_flux_nodal_x_1x1v_ser_p2_Vl[24] = { 
  1.0, 
  2.5, 
  -1.9364916731037085, 
  -4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  7.448044896552726, 
  -3.8461538461538463, 
  1.0, 
  2.5, 
  0.0, 
  0.0, 
  -3.8461538461538463, 
  1.9230769230769231, 
  0.0, 
  4.8076923076923075, 
  1.0, 
  2.5, 
  1.9364916731037085, 
  4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  -7.448044896552726, 
  -3.8461538461538463, 
}; 
static const double lax_flux_nodal_x_1x1v_ser_p2_Vr[24] = { 
  1.0, 
  -2.5, 
  -1.9364916731037085, 
  4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  7.448044896552726, 
  3.8461538461538463, 
  1.0, 
  -2.5, 
  0.0, 
  0.0, 
  -3.8461538461538463, 
  1.9230769230769231, 
  0.0, 
  -4.8076923076923075, 
  1.0, 
  -2.5, 
  1.9364916731037085, 
  -4.841229182759272, 
  -3.8461538461538463, 
  -1.5384615384615385, 
  -7.448044896552726, 
  3.8461538461538463, 
}; 
static const double lax_flux_nodal_x_1x1v_ser_p2_coef[8] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
}; 
static const int lax_flux_nodal_x_1x1v_ser_p2_cmap[8] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
static const int lax_flux_nodal_x_1x1v_ser_p2_vlmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
static const int lax_flux_nodal_x_1x1v_ser_p2_vrmap[8] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_DH double lax_flux_nodal_x_1x1v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double *out = &Fhat_nodal[0]; 
  double G_l[3]; 
  double G_r[3]; 
  for (int q = 0; q < 3; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 8; ++k) { 
    const int a = lax_flux_nodal_x_1x1v_ser_p2_cmap[k]; 
    const double cfl = lax_flux_nodal_x_1x1v_ser_p2_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_x_1x1v_ser_p2_coef[k]*f_r[k]; 
    for (int j = 0; j < 3; ++j) { 
      G_l[j*1 + a] += lax_flux_nodal_x_1x1v_ser_p2_Vl[j*8 + lax_flux_nodal_x_1x1v_ser_p2_vlmap[k]]*cfl; 
      G_r[j*1 + a] += lax_flux_nodal_x_1x1v_ser_p2_Vr[j*8 + lax_flux_nodal_x_1x1v_ser_p2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      const int n = i*3 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 1; ++a) { 
        f_l_quad += lax_flux_nodal_x_1x1v_ser_p2_Cm[i*1 + a]*G_l[j*1 + a]; 
        f_r_quad += lax_flux_nodal_x_1x1v_ser_p2_Cm[i*1 + a]*G_r[j*1 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dx10*alpha_max;

} 
