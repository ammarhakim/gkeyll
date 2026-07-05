#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*3 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*4 + j (inner fastest). 
static const double ho_lax_flux_nodal_vy_1x2v_ser_p2_Cm[12] = { 
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
static const double ho_lax_flux_nodal_vy_1x2v_ser_p2_Vl[32] = { 
  1.0, 
  -2.1528407789851314, 
  2.5, 
  -5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  -5.887823276141479, 
  8.280156842250506, 
  1.0, 
  -0.8499526089621406, 
  2.5, 
  -2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  3.140570528888732, 
  3.269048496008233, 
  1.0, 
  0.8499526089621406, 
  2.5, 
  2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  3.140570528888732, 
  -3.269048496008233, 
  1.0, 
  2.1528407789851314, 
  2.5, 
  5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  -5.887823276141479, 
  -8.280156842250506, 
}; 
static const double ho_lax_flux_nodal_vy_1x2v_ser_p2_Vr[32] = { 
  1.0, 
  -2.1528407789851314, 
  -2.5, 
  5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  5.887823276141479, 
  8.280156842250506, 
  1.0, 
  -0.8499526089621406, 
  -2.5, 
  2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  -3.140570528888732, 
  3.269048496008233, 
  1.0, 
  0.8499526089621406, 
  -2.5, 
  -2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  -3.140570528888732, 
  -3.269048496008233, 
  1.0, 
  2.1528407789851314, 
  -2.5, 
  -5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  5.887823276141479, 
  -8.280156842250506, 
}; 
static const double ho_lax_flux_nodal_vy_1x2v_ser_p2_coef[20] = { 
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
static const int ho_lax_flux_nodal_vy_1x2v_ser_p2_cmap[20] = { 
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
static const int ho_lax_flux_nodal_vy_1x2v_ser_p2_vlmap[20] = { 
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
static const int ho_lax_flux_nodal_vy_1x2v_ser_p2_vrmap[20] = { 
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
GKYL_CU_DH double ho_lax_flux_nodal_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv11 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[16]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  const double jac_tot = jacob_vel_surf_vx[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[4], jacob_vel_surf_r[4]); 
  double G_l[12]; 
  double G_r[12]; 
  for (int q = 0; q < 12; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 20; ++k) { 
    const int a = ho_lax_flux_nodal_vy_1x2v_ser_p2_cmap[k]; 
    const double cfl = ho_lax_flux_nodal_vy_1x2v_ser_p2_coef[k]*f_l[k]; 
    const double cfr = ho_lax_flux_nodal_vy_1x2v_ser_p2_coef[k]*f_r[k]; 
    for (int j = 0; j < 4; ++j) { 
      G_l[j*3 + a] += ho_lax_flux_nodal_vy_1x2v_ser_p2_Vl[j*8 + ho_lax_flux_nodal_vy_1x2v_ser_p2_vlmap[k]]*cfl; 
      G_r[j*3 + a] += ho_lax_flux_nodal_vy_1x2v_ser_p2_Vr[j*8 + ho_lax_flux_nodal_vy_1x2v_ser_p2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      const int n = i*4 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += ho_lax_flux_nodal_vy_1x2v_ser_p2_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += ho_lax_flux_nodal_vy_1x2v_ser_p2_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_tot*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv11*alpha_max/jacob_vel_surf_min;

} 
