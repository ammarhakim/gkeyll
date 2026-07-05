#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*4 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*4 + j (inner fastest). 
static const double lax_flux_nodal_vx_2x3v_ser_p1_Cm[16] = { 
  1.0, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  2.0833333333333335, 
  1.0, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  -1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  1.4433756729740645, 
  2.0833333333333335, 
}; 
static const double lax_flux_nodal_vx_2x3v_ser_p1_Vl[32] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  -3.6084391824351614, 
  2.0833333333333335, 
  5.208333333333333, 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -3.6084391824351614, 
  3.6084391824351614, 
  -2.0833333333333335, 
  -5.208333333333333, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  -1.4433756729740645, 
  3.6084391824351614, 
  -3.6084391824351614, 
  -2.0833333333333335, 
  -5.208333333333333, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  1.4433756729740645, 
  3.6084391824351614, 
  3.6084391824351614, 
  2.0833333333333335, 
  5.208333333333333, 
}; 
static const double lax_flux_nodal_vx_2x3v_ser_p1_Vr[32] = { 
  1.0, 
  -2.5, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  3.6084391824351614, 
  3.6084391824351614, 
  2.0833333333333335, 
  -5.208333333333333, 
  1.0, 
  -2.5, 
  -1.4433756729740645, 
  1.4433756729740645, 
  3.6084391824351614, 
  -3.6084391824351614, 
  -2.0833333333333335, 
  5.208333333333333, 
  1.0, 
  -2.5, 
  1.4433756729740645, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  3.6084391824351614, 
  -2.0833333333333335, 
  5.208333333333333, 
  1.0, 
  -2.5, 
  1.4433756729740645, 
  1.4433756729740645, 
  -3.6084391824351614, 
  -3.6084391824351614, 
  2.0833333333333335, 
  -5.208333333333333, 
}; 
static const double lax_flux_nodal_vx_2x3v_ser_p1_coef[32] = { 
  0.1767766952966368, 
  0.12247448713915887, 
  0.12247448713915887, 
  0.12247448713915887, 
  0.12247448713915887, 
  0.12247448713915887, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.0848528137423857, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.058787753826796275, 
  0.040729350596345144, 
  0.040729350596345144, 
  0.040729350596345144, 
  0.040729350596345144, 
  0.040729350596345144, 
  0.028218121836862214, 
}; 
static const int lax_flux_nodal_vx_2x3v_ser_p1_cmap[32] = { 
  0, 
  1, 
  2, 
  0, 
  0, 
  0, 
  3, 
  1, 
  2, 
  1, 
  2, 
  0, 
  1, 
  2, 
  0, 
  0, 
  3, 
  3, 
  1, 
  2, 
  3, 
  1, 
  2, 
  1, 
  2, 
  0, 
  3, 
  3, 
  3, 
  1, 
  2, 
  3, 
}; 
static const int lax_flux_nodal_vx_2x3v_ser_p1_vlmap[32] = { 
  0, 
  0, 
  0, 
  1, 
  2, 
  3, 
  0, 
  1, 
  1, 
  2, 
  2, 
  4, 
  3, 
  3, 
  5, 
  6, 
  1, 
  2, 
  4, 
  4, 
  3, 
  5, 
  5, 
  6, 
  6, 
  7, 
  4, 
  5, 
  6, 
  7, 
  7, 
  7, 
}; 
static const int lax_flux_nodal_vx_2x3v_ser_p1_vrmap[32] = { 
  0, 
  0, 
  0, 
  1, 
  2, 
  3, 
  0, 
  1, 
  1, 
  2, 
  2, 
  4, 
  3, 
  3, 
  5, 
  6, 
  1, 
  2, 
  4, 
  4, 
  3, 
  5, 
  5, 
  6, 
  6, 
  7, 
  4, 
  5, 
  6, 
  7, 
  7, 
  7, 
}; 
GKYL_CU_DH double lax_flux_nodal_vx_2x3v_ser_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[2]; 

  double *out = &Fhat_nodal[0]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf_r[3]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf_r[6]; 
  const double jac_tot = jacob_vel_surf_vy[0]*jacob_vel_surf_vz[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  double G_l[16]; 
  double G_r[16]; 
  for (int q = 0; q < 16; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 32; ++k) { 
    const int a = lax_flux_nodal_vx_2x3v_ser_p1_cmap[k]; 
    const double cfl = lax_flux_nodal_vx_2x3v_ser_p1_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vx_2x3v_ser_p1_coef[k]*f_r[k]; 
    for (int j = 0; j < 4; ++j) { 
      G_l[j*4 + a] += lax_flux_nodal_vx_2x3v_ser_p1_Vl[j*8 + lax_flux_nodal_vx_2x3v_ser_p1_vlmap[k]]*cfl; 
      G_r[j*4 + a] += lax_flux_nodal_vx_2x3v_ser_p1_Vr[j*8 + lax_flux_nodal_vx_2x3v_ser_p1_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      const int n = i*4 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 4; ++a) { 
        f_l_quad += lax_flux_nodal_vx_2x3v_ser_p1_Cm[i*4 + a]*G_l[j*4 + a]; 
        f_r_quad += lax_flux_nodal_vx_2x3v_ser_p1_Cm[i*4 + a]*G_r[j*4 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_tot*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 1.5*dv10*alpha_max/jacob_vel_surf_min;

} 
