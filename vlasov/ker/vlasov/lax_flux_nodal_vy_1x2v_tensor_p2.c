#include <gkyl_vlasov_kernels.h> 
// Exact factorized evaluation of the phase basis at the surface tensor-product nodes: 
// phi_k = coef[k] * (outer shape cmap[k])(outer vars) * (inner shape)(normal, inner vars). 
// Cm[i*3 + a]: outer shape a at outer node i; Vl/Vr[j*NM + m]: inner shape m at the 
// +1/-1 edge and inner node j. Flat surface node index n = i*4 + j (inner fastest). 
static const double lax_flux_nodal_vy_1x2v_tensor_p2_Cm[12] = { 
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
static const double lax_flux_nodal_vy_1x2v_tensor_p2_Vl[36] = { 
  1.0, 
  -2.1528407789851314, 
  2.5, 
  -5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  -5.887823276141479, 
  8.280156842250506, 
  9.058189655602275, 
  1.0, 
  -0.8499526089621406, 
  2.5, 
  -2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  3.140570528888732, 
  3.269048496008233, 
  -4.8316469675211255, 
  1.0, 
  0.8499526089621406, 
  2.5, 
  2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  3.140570528888732, 
  -3.269048496008233, 
  -4.8316469675211255, 
  1.0, 
  2.1528407789851314, 
  2.5, 
  5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  -5.887823276141479, 
  -8.280156842250506, 
  9.058189655602275, 
}; 
static const double lax_flux_nodal_vy_1x2v_tensor_p2_Vr[36] = { 
  1.0, 
  -2.1528407789851314, 
  -2.5, 
  5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  5.887823276141479, 
  8.280156842250506, 
  9.058189655602275, 
  1.0, 
  -0.8499526089621406, 
  -2.5, 
  2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  -3.140570528888732, 
  3.269048496008233, 
  -4.8316469675211255, 
  1.0, 
  0.8499526089621406, 
  -2.5, 
  -2.1248815224053517, 
  1.2562282115554928, 
  -3.8461538461538463, 
  -3.140570528888732, 
  -3.269048496008233, 
  -4.8316469675211255, 
  1.0, 
  2.1528407789851314, 
  -2.5, 
  -5.382101947462829, 
  -2.355129310456592, 
  -3.8461538461538463, 
  5.887823276141479, 
  -8.280156842250506, 
  9.058189655602275, 
}; 
static const double lax_flux_nodal_vy_1x2v_tensor_p2_coef[27] = { 
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
  0.1195010460205265, 
  0.1195010460205265, 
  0.1195010460205265, 
  0.08279275330607137, 
  0.08279275330607137, 
  0.08279275330607137, 
  -0.06947524019389922, 
}; 
static const int lax_flux_nodal_vy_1x2v_tensor_p2_cmap[27] = { 
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
  2, 
  2, 
  0, 
  2, 
  2, 
  1, 
  2, 
}; 
static const int lax_flux_nodal_vy_1x2v_tensor_p2_vlmap[27] = { 
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
  4, 
  5, 
  8, 
  6, 
  7, 
  8, 
  8, 
}; 
static const int lax_flux_nodal_vy_1x2v_tensor_p2_vrmap[27] = { 
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
  4, 
  5, 
  8, 
  6, 
  7, 
  8, 
  8, 
}; 
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
    const int a = lax_flux_nodal_vy_1x2v_tensor_p2_cmap[k]; 
    const double cfl = lax_flux_nodal_vy_1x2v_tensor_p2_coef[k]*f_l[k]; 
    const double cfr = lax_flux_nodal_vy_1x2v_tensor_p2_coef[k]*f_r[k]; 
    for (int j = 0; j < 4; ++j) { 
      G_l[j*3 + a] += lax_flux_nodal_vy_1x2v_tensor_p2_Vl[j*9 + lax_flux_nodal_vy_1x2v_tensor_p2_vlmap[k]]*cfl; 
      G_r[j*3 + a] += lax_flux_nodal_vy_1x2v_tensor_p2_Vr[j*9 + lax_flux_nodal_vy_1x2v_tensor_p2_vrmap[k]]*cfr; 
    } 
  } 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      const int n = i*4 + j; 
      double f_l_quad = 0.0; 
      double f_r_quad = 0.0; 
      for (int a = 0; a < 3; ++a) { 
        f_l_quad += lax_flux_nodal_vy_1x2v_tensor_p2_Cm[i*3 + a]*G_l[j*3 + a]; 
        f_r_quad += lax_flux_nodal_vy_1x2v_tensor_p2_Cm[i*3 + a]*G_r[j*3 + a]; 
      } 
      alpha_max = fmax(alpha_max, fabs(alpha_quad[n])); 
      out[n] = 0.5*jac_nodal[n]*(alpha_quad[n]*(f_r_quad + f_l_quad) - fabs(alpha_quad[n])*(f_r_quad - f_l_quad)); 
    } 
  } 
  return 2.5*dv11*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vy[3])+0.40076152031165013*jacob_vel_surf_vy[2]-0.8136324494869249*jacob_vel_surf_vy[1]+1.5267881254572662*jacob_vel_surf_vy[0]);

} 
