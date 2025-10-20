#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_u_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf) 
{ 
  // count:  integer to indicate which matrix being fetched. 
  // x:      Input solution vector. 
  // u:      Output volume expansion of flow velocity: 
  //         [ux, uy, uz]. 
  // u_surf: Output surface expansion of flow velocity 
  //         [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //          ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr,  
  //          ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  double *ux = &u[0]; 
  double *uy = &u[9]; 
  double *uz = &u[18]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 
  ux[8] = gkyl_mat_get(&x_ux,8,0); 
  uy[8] = gkyl_mat_get(&x_uy,8,0); 
  uz[8] = gkyl_mat_get(&x_uz,8,0); 

  double *ux_xl = &u_surf[0]; 
  double *ux_xr = &u_surf[3]; 
  double *uy_xl = &u_surf[6]; 
  double *uy_xr = &u_surf[9]; 
  double *uz_xl = &u_surf[12]; 
  double *uz_xr = &u_surf[15]; 
 
  ux_xl[0] = 1.58113883008419*ux[4]-1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xl[1] = 1.58113883008419*ux[6]-1.224744871391589*ux[3]+0.7071067811865475*ux[2]; 
  ux_xl[2] = 1.58113883008419*ux[8]-1.224744871391589*ux[7]+0.7071067811865475*ux[5]; 
  uy_xl[0] = 1.58113883008419*uy[4]-1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xl[1] = 1.58113883008419*uy[6]-1.224744871391589*uy[3]+0.7071067811865475*uy[2]; 
  uy_xl[2] = 1.58113883008419*uy[8]-1.224744871391589*uy[7]+0.7071067811865475*uy[5]; 
  uz_xl[0] = 1.58113883008419*uz[4]-1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xl[1] = 1.58113883008419*uz[6]-1.224744871391589*uz[3]+0.7071067811865475*uz[2]; 
  uz_xl[2] = 1.58113883008419*uz[8]-1.224744871391589*uz[7]+0.7071067811865475*uz[5]; 
 
  ux_xr[0] = 1.58113883008419*ux[4]+1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xr[1] = 1.58113883008419*ux[6]+1.224744871391589*ux[3]+0.7071067811865475*ux[2]; 
  ux_xr[2] = 1.58113883008419*ux[8]+1.224744871391589*ux[7]+0.7071067811865475*ux[5]; 
  uy_xr[0] = 1.58113883008419*uy[4]+1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xr[1] = 1.58113883008419*uy[6]+1.224744871391589*uy[3]+0.7071067811865475*uy[2]; 
  uy_xr[2] = 1.58113883008419*uy[8]+1.224744871391589*uy[7]+0.7071067811865475*uy[5]; 
  uz_xr[0] = 1.58113883008419*uz[4]+1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xr[1] = 1.58113883008419*uz[6]+1.224744871391589*uz[3]+0.7071067811865475*uz[2]; 
  uz_xr[2] = 1.58113883008419*uz[8]+1.224744871391589*uz[7]+0.7071067811865475*uz[5]; 
 
  double *ux_yl = &u_surf[18]; 
  double *ux_yr = &u_surf[21]; 
  double *uy_yl = &u_surf[24]; 
  double *uy_yr = &u_surf[27]; 
  double *uz_yl = &u_surf[30]; 
  double *uz_yr = &u_surf[33]; 
 
  ux_yl[0] = 1.58113883008419*ux[5]-1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yl[1] = 1.58113883008419*ux[7]-1.224744871391589*ux[3]+0.7071067811865475*ux[1]; 
  ux_yl[2] = 1.58113883008419*ux[8]-1.224744871391589*ux[6]+0.7071067811865475*ux[4]; 
  uy_yl[0] = 1.58113883008419*uy[5]-1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yl[1] = 1.58113883008419*uy[7]-1.224744871391589*uy[3]+0.7071067811865475*uy[1]; 
  uy_yl[2] = 1.58113883008419*uy[8]-1.224744871391589*uy[6]+0.7071067811865475*uy[4]; 
  uz_yl[0] = 1.58113883008419*uz[5]-1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yl[1] = 1.58113883008419*uz[7]-1.224744871391589*uz[3]+0.7071067811865475*uz[1]; 
  uz_yl[2] = 1.58113883008419*uz[8]-1.224744871391589*uz[6]+0.7071067811865475*uz[4]; 
 
  ux_yr[0] = 1.58113883008419*ux[5]+1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yr[1] = 1.58113883008419*ux[7]+1.224744871391589*ux[3]+0.7071067811865475*ux[1]; 
  ux_yr[2] = 1.58113883008419*ux[8]+1.224744871391589*ux[6]+0.7071067811865475*ux[4]; 
  uy_yr[0] = 1.58113883008419*uy[5]+1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yr[1] = 1.58113883008419*uy[7]+1.224744871391589*uy[3]+0.7071067811865475*uy[1]; 
  uy_yr[2] = 1.58113883008419*uy[8]+1.224744871391589*uy[6]+0.7071067811865475*uy[4]; 
  uz_yr[0] = 1.58113883008419*uz[5]+1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yr[1] = 1.58113883008419*uz[7]+1.224744871391589*uz[3]+0.7071067811865475*uz[1]; 
  uz_yr[2] = 1.58113883008419*uz[8]+1.224744871391589*uz[6]+0.7071067811865475*uz[4]; 
 
} 
 
