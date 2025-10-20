#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p1(int count, struct gkyl_nmat *x, 
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
  double *uy = &prim[8]; 
  double *uz = &prim[16]; 
  double *p_force = &prim[24]; 
  double *T_perp_over_m = &prim[32]; 
  double *T_perp_over_m_inv = &prim[40]; 
  double *Txx = &prim[48]; 
  double *Tyy = &prim[56]; 
  double *Tzz = &prim[64]; 

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
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  p_force[2] = gkyl_mat_get(&x_pkpm_div_ppar,2,0); 
  T_perp_over_m[2] = gkyl_mat_get(&x_T_perp_over_m,2,0); 
  T_perp_over_m_inv[2] = gkyl_mat_get(&x_T_perp_over_m_inv,2,0); 
  Txx[2] = 3.0*gkyl_mat_get(&x_Txx,2,0); 
  Tyy[2] = 3.0*gkyl_mat_get(&x_Tyy,2,0); 
  Tzz[2] = 3.0*gkyl_mat_get(&x_Tzz,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  p_force[3] = gkyl_mat_get(&x_pkpm_div_ppar,3,0); 
  T_perp_over_m[3] = gkyl_mat_get(&x_T_perp_over_m,3,0); 
  T_perp_over_m_inv[3] = gkyl_mat_get(&x_T_perp_over_m_inv,3,0); 
  Txx[3] = 3.0*gkyl_mat_get(&x_Txx,3,0); 
  Tyy[3] = 3.0*gkyl_mat_get(&x_Tyy,3,0); 
  Tzz[3] = 3.0*gkyl_mat_get(&x_Tzz,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  p_force[4] = gkyl_mat_get(&x_pkpm_div_ppar,4,0); 
  T_perp_over_m[4] = gkyl_mat_get(&x_T_perp_over_m,4,0); 
  T_perp_over_m_inv[4] = gkyl_mat_get(&x_T_perp_over_m_inv,4,0); 
  Txx[4] = 3.0*gkyl_mat_get(&x_Txx,4,0); 
  Tyy[4] = 3.0*gkyl_mat_get(&x_Tyy,4,0); 
  Tzz[4] = 3.0*gkyl_mat_get(&x_Tzz,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  p_force[5] = gkyl_mat_get(&x_pkpm_div_ppar,5,0); 
  T_perp_over_m[5] = gkyl_mat_get(&x_T_perp_over_m,5,0); 
  T_perp_over_m_inv[5] = gkyl_mat_get(&x_T_perp_over_m_inv,5,0); 
  Txx[5] = 3.0*gkyl_mat_get(&x_Txx,5,0); 
  Tyy[5] = 3.0*gkyl_mat_get(&x_Tyy,5,0); 
  Tzz[5] = 3.0*gkyl_mat_get(&x_Tzz,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  p_force[6] = gkyl_mat_get(&x_pkpm_div_ppar,6,0); 
  T_perp_over_m[6] = gkyl_mat_get(&x_T_perp_over_m,6,0); 
  T_perp_over_m_inv[6] = gkyl_mat_get(&x_T_perp_over_m_inv,6,0); 
  Txx[6] = 3.0*gkyl_mat_get(&x_Txx,6,0); 
  Tyy[6] = 3.0*gkyl_mat_get(&x_Tyy,6,0); 
  Tzz[6] = 3.0*gkyl_mat_get(&x_Tzz,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 
  p_force[7] = gkyl_mat_get(&x_pkpm_div_ppar,7,0); 
  T_perp_over_m[7] = gkyl_mat_get(&x_T_perp_over_m,7,0); 
  T_perp_over_m_inv[7] = gkyl_mat_get(&x_T_perp_over_m_inv,7,0); 
  Txx[7] = 3.0*gkyl_mat_get(&x_Txx,7,0); 
  Tyy[7] = 3.0*gkyl_mat_get(&x_Tyy,7,0); 
  Tzz[7] = 3.0*gkyl_mat_get(&x_Tzz,7,0); 

  double *ux_xl = &prim_surf[0]; 
  double *ux_xr = &prim_surf[4]; 
  double *uy_xl = &prim_surf[8]; 
  double *uy_xr = &prim_surf[12]; 
  double *uz_xl = &prim_surf[16]; 
  double *uz_xr = &prim_surf[20]; 
  double *Txx_xl = &prim_surf[24]; 
  double *Txx_xr = &prim_surf[28]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xl[1] = 0.7071067811865475*ux[2]-1.224744871391589*ux[4]; 
  ux_xl[2] = 0.7071067811865475*ux[3]-1.224744871391589*ux[5]; 
  ux_xl[3] = 0.7071067811865475*ux[6]-1.224744871391589*ux[7]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xl[1] = 0.7071067811865475*uy[2]-1.224744871391589*uy[4]; 
  uy_xl[2] = 0.7071067811865475*uy[3]-1.224744871391589*uy[5]; 
  uy_xl[3] = 0.7071067811865475*uy[6]-1.224744871391589*uy[7]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xl[1] = 0.7071067811865475*uz[2]-1.224744871391589*uz[4]; 
  uz_xl[2] = 0.7071067811865475*uz[3]-1.224744871391589*uz[5]; 
  uz_xl[3] = 0.7071067811865475*uz[6]-1.224744871391589*uz[7]; 
  Txx_xl[0] = 0.7071067811865475*Txx[0]-1.224744871391589*Txx[1]; 
  Txx_xl[1] = 0.7071067811865475*Txx[2]-1.224744871391589*Txx[4]; 
  Txx_xl[2] = 0.7071067811865475*Txx[3]-1.224744871391589*Txx[5]; 
  Txx_xl[3] = 0.7071067811865475*Txx[6]-1.224744871391589*Txx[7]; 
 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xr[1] = 1.224744871391589*ux[4]+0.7071067811865475*ux[2]; 
  ux_xr[2] = 1.224744871391589*ux[5]+0.7071067811865475*ux[3]; 
  ux_xr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[6]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xr[1] = 1.224744871391589*uy[4]+0.7071067811865475*uy[2]; 
  uy_xr[2] = 1.224744871391589*uy[5]+0.7071067811865475*uy[3]; 
  uy_xr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[6]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xr[1] = 1.224744871391589*uz[4]+0.7071067811865475*uz[2]; 
  uz_xr[2] = 1.224744871391589*uz[5]+0.7071067811865475*uz[3]; 
  uz_xr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[6]; 
  Txx_xr[0] = 1.224744871391589*Txx[1]+0.7071067811865475*Txx[0]; 
  Txx_xr[1] = 1.224744871391589*Txx[4]+0.7071067811865475*Txx[2]; 
  Txx_xr[2] = 1.224744871391589*Txx[5]+0.7071067811865475*Txx[3]; 
  Txx_xr[3] = 1.224744871391589*Txx[7]+0.7071067811865475*Txx[6]; 
 
  double *ux_yl = &prim_surf[32]; 
  double *ux_yr = &prim_surf[36]; 
  double *uy_yl = &prim_surf[40]; 
  double *uy_yr = &prim_surf[44]; 
  double *uz_yl = &prim_surf[48]; 
  double *uz_yr = &prim_surf[52]; 
  double *Tyy_yl = &prim_surf[56]; 
  double *Tyy_yr = &prim_surf[60]; 
 
  ux_yl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[2]; 
  ux_yl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[4]; 
  ux_yl[2] = 0.7071067811865475*ux[3]-1.224744871391589*ux[6]; 
  ux_yl[3] = 0.7071067811865475*ux[5]-1.224744871391589*ux[7]; 
  uy_yl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[2]; 
  uy_yl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[4]; 
  uy_yl[2] = 0.7071067811865475*uy[3]-1.224744871391589*uy[6]; 
  uy_yl[3] = 0.7071067811865475*uy[5]-1.224744871391589*uy[7]; 
  uz_yl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[2]; 
  uz_yl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[4]; 
  uz_yl[2] = 0.7071067811865475*uz[3]-1.224744871391589*uz[6]; 
  uz_yl[3] = 0.7071067811865475*uz[5]-1.224744871391589*uz[7]; 
  Tyy_yl[0] = 0.7071067811865475*Tyy[0]-1.224744871391589*Tyy[2]; 
  Tyy_yl[1] = 0.7071067811865475*Tyy[1]-1.224744871391589*Tyy[4]; 
  Tyy_yl[2] = 0.7071067811865475*Tyy[3]-1.224744871391589*Tyy[6]; 
  Tyy_yl[3] = 0.7071067811865475*Tyy[5]-1.224744871391589*Tyy[7]; 
 
  ux_yr[0] = 1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yr[1] = 1.224744871391589*ux[4]+0.7071067811865475*ux[1]; 
  ux_yr[2] = 1.224744871391589*ux[6]+0.7071067811865475*ux[3]; 
  ux_yr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[5]; 
  uy_yr[0] = 1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yr[1] = 1.224744871391589*uy[4]+0.7071067811865475*uy[1]; 
  uy_yr[2] = 1.224744871391589*uy[6]+0.7071067811865475*uy[3]; 
  uy_yr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[5]; 
  uz_yr[0] = 1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yr[1] = 1.224744871391589*uz[4]+0.7071067811865475*uz[1]; 
  uz_yr[2] = 1.224744871391589*uz[6]+0.7071067811865475*uz[3]; 
  uz_yr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[5]; 
  Tyy_yr[0] = 1.224744871391589*Tyy[2]+0.7071067811865475*Tyy[0]; 
  Tyy_yr[1] = 1.224744871391589*Tyy[4]+0.7071067811865475*Tyy[1]; 
  Tyy_yr[2] = 1.224744871391589*Tyy[6]+0.7071067811865475*Tyy[3]; 
  Tyy_yr[3] = 1.224744871391589*Tyy[7]+0.7071067811865475*Tyy[5]; 
 
  double *ux_zl = &prim_surf[64]; 
  double *ux_zr = &prim_surf[68]; 
  double *uy_zl = &prim_surf[72]; 
  double *uy_zr = &prim_surf[76]; 
  double *uz_zl = &prim_surf[80]; 
  double *uz_zr = &prim_surf[84]; 
  double *Tzz_zl = &prim_surf[88]; 
  double *Tzz_zr = &prim_surf[92]; 
 
  ux_zl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[3]; 
  ux_zl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[5]; 
  ux_zl[2] = 0.7071067811865475*ux[2]-1.224744871391589*ux[6]; 
  ux_zl[3] = 0.7071067811865475*ux[4]-1.224744871391589*ux[7]; 
  uy_zl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[3]; 
  uy_zl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[5]; 
  uy_zl[2] = 0.7071067811865475*uy[2]-1.224744871391589*uy[6]; 
  uy_zl[3] = 0.7071067811865475*uy[4]-1.224744871391589*uy[7]; 
  uz_zl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[3]; 
  uz_zl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[5]; 
  uz_zl[2] = 0.7071067811865475*uz[2]-1.224744871391589*uz[6]; 
  uz_zl[3] = 0.7071067811865475*uz[4]-1.224744871391589*uz[7]; 
  Tzz_zl[0] = 0.7071067811865475*Tzz[0]-1.224744871391589*Tzz[3]; 
  Tzz_zl[1] = 0.7071067811865475*Tzz[1]-1.224744871391589*Tzz[5]; 
  Tzz_zl[2] = 0.7071067811865475*Tzz[2]-1.224744871391589*Tzz[6]; 
  Tzz_zl[3] = 0.7071067811865475*Tzz[4]-1.224744871391589*Tzz[7]; 
 
  ux_zr[0] = 1.224744871391589*ux[3]+0.7071067811865475*ux[0]; 
  ux_zr[1] = 1.224744871391589*ux[5]+0.7071067811865475*ux[1]; 
  ux_zr[2] = 1.224744871391589*ux[6]+0.7071067811865475*ux[2]; 
  ux_zr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[4]; 
  uy_zr[0] = 1.224744871391589*uy[3]+0.7071067811865475*uy[0]; 
  uy_zr[1] = 1.224744871391589*uy[5]+0.7071067811865475*uy[1]; 
  uy_zr[2] = 1.224744871391589*uy[6]+0.7071067811865475*uy[2]; 
  uy_zr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[4]; 
  uz_zr[0] = 1.224744871391589*uz[3]+0.7071067811865475*uz[0]; 
  uz_zr[1] = 1.224744871391589*uz[5]+0.7071067811865475*uz[1]; 
  uz_zr[2] = 1.224744871391589*uz[6]+0.7071067811865475*uz[2]; 
  uz_zr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[4]; 
  Tzz_zr[0] = 1.224744871391589*Tzz[3]+0.7071067811865475*Tzz[0]; 
  Tzz_zr[1] = 1.224744871391589*Tzz[5]+0.7071067811865475*Tzz[1]; 
  Tzz_zr[2] = 1.224744871391589*Tzz[6]+0.7071067811865475*Tzz[2]; 
  Tzz_zr[3] = 1.224744871391589*Tzz[7]+0.7071067811865475*Tzz[4]; 
 
} 
 
