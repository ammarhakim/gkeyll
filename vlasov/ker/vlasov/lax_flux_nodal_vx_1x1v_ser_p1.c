#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*2 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*1 + j (inner fastest). 
static const double lax_flux_nodal_vx_1x1v_ser_p1_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p1_Vl[2] = { 
  1.0, 
  2.5, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p1_Vr[2] = { 
  1.0, 
  -2.5, 
}; 
static const double lax_flux_nodal_vx_1x1v_ser_p1_coef[4] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p1_cmap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p1_vlmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
static const int lax_flux_nodal_vx_1x1v_ser_p1_vrmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_DH double lax_flux_nodal_vx_1x1v_ser_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double *lax, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double dv10 = 2.0/dxv[1]; 

  double *out = &Fhat_nodal[0]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  double G_l[2]; 
  double G_r[2]; 
  for (int q = 0; q < 2; ++q) { G_l[q] = 0.0; G_r[q] = 0.0; } 
  for (int k = 0; k < 4; ++k) { 
    const int a = lax_flux_nodal_vx_1x1v_ser_p1_cmap[k]; 
    const double cfl = lax_flux_nodal_vx_1x1v_ser_p1_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vx_1x1v_ser_p1_coef[k]*f_r[k]; 
    for (int j = 0; j < 1; ++j) { 
      G_l[j*2 + a] += lax_flux_nodal_vx_1x1v_ser_p1_Vl[j*2 + lax_flux_nodal_vx_1x1v_ser_p1_vlmap[k]]*cfl; 
      G_r[j*2 + a] += lax_flux_nodal_vx_1x1v_ser_p1_Vr[j*2 + lax_flux_nodal_vx_1x1v_ser_p1_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      const int n = i*1 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 2; ++a) { 
        f_l_quad += lax_flux_nodal_vx_1x1v_ser_p1_Cm[i*2 + a]*G_l[j*2 + a]; 
        f_r_quad += lax_flux_nodal_vx_1x1v_ser_p1_Cm[i*2 + a]*G_r[j*2 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 1.5*dv10*alpha_max/jacob_vel_surf_min;

} 
