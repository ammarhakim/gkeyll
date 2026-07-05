#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*3 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*3 + j (inner fastest). 
static const double lax_flux_nodal_vy_1x2v_ser_p2_Cm[9] = { 
  1.0, 
  -1.9364916731037085, 
  -1.5384615384615385, 
  1.0, 
  0.0, 
  1.9230769230769231, 
  1.0, 
  1.9364916731037085, 
  -1.5384615384615385, 
}; 
static const double lax_flux_nodal_vy_1x2v_ser_p2_Vl[24] = { 
  1.0, 
  -1.9364916731037085, 
  2.5, 
  -4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
  -3.8461538461538463, 
  7.448044896552726, 
  1.0, 
  0.0, 
  2.5, 
  0.0, 
  1.9230769230769231, 
  -3.8461538461538463, 
  4.8076923076923075, 
  0.0, 
  1.0, 
  1.9364916731037085, 
  2.5, 
  4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
  -3.8461538461538463, 
  -7.448044896552726, 
}; 
static const double lax_flux_nodal_vy_1x2v_ser_p2_Vr[24] = { 
  1.0, 
  -1.9364916731037085, 
  -2.5, 
  4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
  3.8461538461538463, 
  7.448044896552726, 
  1.0, 
  0.0, 
  -2.5, 
  0.0, 
  1.9230769230769231, 
  -3.8461538461538463, 
  -4.8076923076923075, 
  0.0, 
  1.0, 
  1.9364916731037085, 
  -2.5, 
  -4.841229182759272, 
  -1.5384615384615385, 
  -3.8461538461538463, 
  3.8461538461538463, 
  -7.448044896552726, 
}; 
static const double lax_flux_nodal_vy_1x2v_ser_p2_coef[20] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.20554804791094453, 
  -0.20554804791094453, 
  -0.20554804791094453, 
  0.11757550765359256, 
  -0.14240786495134314, 
  -0.14240786495134317, 
  -0.14240786495134314, 
  -0.14240786495134314, 
  -0.14240786495134317, 
  -0.14240786495134317, 
  -0.0986630629972534, 
  -0.0986630629972534, 
  -0.0986630629972534, 
}; 
static const int lax_flux_nodal_vy_1x2v_ser_p2_cmap[20] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  0, 
  1, 
  0, 
  2, 
  1, 
  1, 
}; 
static const int lax_flux_nodal_vy_1x2v_ser_p2_vlmap[20] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  0, 
  4, 
  5, 
  3, 
  1, 
  4, 
  2, 
  6, 
  5, 
  7, 
  3, 
  6, 
  7, 
}; 
static const int lax_flux_nodal_vy_1x2v_ser_p2_vrmap[20] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  0, 
  4, 
  5, 
  3, 
  1, 
  4, 
  2, 
  6, 
  5, 
  7, 
  3, 
  6, 
  7, 
}; 
GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv11 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[9]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  const double jac_tot = jacob_vel_surf_vx[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[4], jacob_vel_surf_r[4]); 
  double G_l[9]; 
  double G_r[9]; 
  for (int q = 0; q < 9; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 20; ++k) { 
    const int a = lax_flux_nodal_vy_1x2v_ser_p2_cmap[k]; 
    const double cfl = lax_flux_nodal_vy_1x2v_ser_p2_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vy_1x2v_ser_p2_coef[k]*f_r[k]; 
    for (int j = 0; j < 3; ++j) { 
      G_l[j*3 + a] += lax_flux_nodal_vy_1x2v_ser_p2_Vl[j*8 + lax_flux_nodal_vy_1x2v_ser_p2_vlmap[k]]*cfl; 
      G_r[j*3 + a] += lax_flux_nodal_vy_1x2v_ser_p2_Vr[j*8 + lax_flux_nodal_vy_1x2v_ser_p2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      const int n = i*3 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += lax_flux_nodal_vy_1x2v_ser_p2_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += lax_flux_nodal_vy_1x2v_ser_p2_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_tot*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv11*alpha_max/jacob_vel_surf_min;

} 
