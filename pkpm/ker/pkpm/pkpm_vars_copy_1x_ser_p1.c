#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p1(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf) 
{ 
  // count:     integer to indicate which matrix being fetched. 
  // x:         Input solution vector. 
  // prim:      Output volume expansion of primitive variables: 
  //            [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]. 
  // prim_surf: Output surface expansion of primitive variables 
  //            [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, Txx_xl, Txx_xr,  
  //             ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, Tyy_yl, Tyy_yr,  
  //             ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, Tzz_zl, Tzz_zr]  
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_pkpm_div_ppar = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_T_perp_over_m = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_T_perp_over_m_inv = gkyl_nmat_get(x, count+5); 
  struct gkyl_mat x_Txx = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_Tyy = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_Tzz = gkyl_nmat_get(x, count+8); 
  double *ux = &prim[0]; 
  double *uy = &prim[2]; 
  double *uz = &prim[4]; 
  double *p_force = &prim[6]; 
  double *T_perp_over_m = &prim[8]; 
  double *T_perp_over_m_inv = &prim[10]; 
  double *Txx = &prim[12]; 
  double *Tyy = &prim[14]; 
  double *Tzz = &prim[16]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  p_force[0] = gkyl_mat_get(&x_pkpm_div_ppar,0,0); 
  T_perp_over_m[0] = gkyl_mat_get(&x_T_perp_over_m,0,0); 
  T_perp_over_m_inv[0] = gkyl_mat_get(&x_T_perp_over_m_inv,0,0); 
  Txx[0] = 3.0*gkyl_mat_get(&x_Txx,0,0); 
  Tyy[0] = 3.0*gkyl_mat_get(&x_Tyy,0,0); 
  Tzz[0] = 3.0*gkyl_mat_get(&x_Tzz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  p_force[1] = gkyl_mat_get(&x_pkpm_div_ppar,1,0); 
  T_perp_over_m[1] = gkyl_mat_get(&x_T_perp_over_m,1,0); 
  T_perp_over_m_inv[1] = gkyl_mat_get(&x_T_perp_over_m_inv,1,0); 
  Txx[1] = 3.0*gkyl_mat_get(&x_Txx,1,0); 
  Tyy[1] = 3.0*gkyl_mat_get(&x_Tyy,1,0); 
  Tzz[1] = 3.0*gkyl_mat_get(&x_Tzz,1,0); 

  double *ux_xl = &prim_surf[0]; 
  double *ux_xr = &prim_surf[1]; 
  double *uy_xl = &prim_surf[2]; 
  double *uy_xr = &prim_surf[3]; 
  double *uz_xl = &prim_surf[4]; 
  double *uz_xr = &prim_surf[5]; 
  double *Txx_xl = &prim_surf[6]; 
  double *Txx_xr = &prim_surf[7]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  Txx_xl[0] = 0.7071067811865475*Txx[0]-1.224744871391589*Txx[1]; 
  Txx_xr[0] = 1.224744871391589*Txx[1]+0.7071067811865475*Txx[0]; 
 
} 
 
