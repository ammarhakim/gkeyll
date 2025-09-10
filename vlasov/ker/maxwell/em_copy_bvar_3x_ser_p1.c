#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_3x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_3x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:          Integer to indicate which matrix being fetched. 
  // x:              Input solution vector. 
  // em:             Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. 
  // bvar:           Output volume expansion of b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  struct gkyl_mat x_bxbx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxby = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_bxbz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_byby = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_bybz = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_bzbz = gkyl_nmat_get(x, count+5); 
  double *bx = &bvar[0]; 
  double *by = &bvar[8]; 
  double *bz = &bvar[16]; 
  double *bxbx = &bvar[24]; 
  double *bxby = &bvar[32]; 
  double *bxbz = &bvar[40]; 
  double *byby = &bvar[48]; 
  double *bybz = &bvar[56]; 
  double *bzbz = &bvar[64]; 
 
  bxbx[0] = gkyl_mat_get(&x_bxbx,0,0); 
  bxby[0] = gkyl_mat_get(&x_bxby,0,0); 
  bxbz[0] = gkyl_mat_get(&x_bxbz,0,0); 
  byby[0] = gkyl_mat_get(&x_byby,0,0); 
  bybz[0] = gkyl_mat_get(&x_bybz,0,0); 
  bzbz[0] = gkyl_mat_get(&x_bzbz,0,0); 
  bxbx[1] = gkyl_mat_get(&x_bxbx,1,0); 
  bxby[1] = gkyl_mat_get(&x_bxby,1,0); 
  bxbz[1] = gkyl_mat_get(&x_bxbz,1,0); 
  byby[1] = gkyl_mat_get(&x_byby,1,0); 
  bybz[1] = gkyl_mat_get(&x_bybz,1,0); 
  bzbz[1] = gkyl_mat_get(&x_bzbz,1,0); 
  bxbx[2] = gkyl_mat_get(&x_bxbx,2,0); 
  bxby[2] = gkyl_mat_get(&x_bxby,2,0); 
  bxbz[2] = gkyl_mat_get(&x_bxbz,2,0); 
  byby[2] = gkyl_mat_get(&x_byby,2,0); 
  bybz[2] = gkyl_mat_get(&x_bybz,2,0); 
  bzbz[2] = gkyl_mat_get(&x_bzbz,2,0); 
  bxbx[3] = gkyl_mat_get(&x_bxbx,3,0); 
  bxby[3] = gkyl_mat_get(&x_bxby,3,0); 
  bxbz[3] = gkyl_mat_get(&x_bxbz,3,0); 
  byby[3] = gkyl_mat_get(&x_byby,3,0); 
  bybz[3] = gkyl_mat_get(&x_bybz,3,0); 
  bzbz[3] = gkyl_mat_get(&x_bzbz,3,0); 
  bxbx[4] = gkyl_mat_get(&x_bxbx,4,0); 
  bxby[4] = gkyl_mat_get(&x_bxby,4,0); 
  bxbz[4] = gkyl_mat_get(&x_bxbz,4,0); 
  byby[4] = gkyl_mat_get(&x_byby,4,0); 
  bybz[4] = gkyl_mat_get(&x_bybz,4,0); 
  bzbz[4] = gkyl_mat_get(&x_bzbz,4,0); 
  bxbx[5] = gkyl_mat_get(&x_bxbx,5,0); 
  bxby[5] = gkyl_mat_get(&x_bxby,5,0); 
  bxbz[5] = gkyl_mat_get(&x_bxbz,5,0); 
  byby[5] = gkyl_mat_get(&x_byby,5,0); 
  bybz[5] = gkyl_mat_get(&x_bybz,5,0); 
  bzbz[5] = gkyl_mat_get(&x_bzbz,5,0); 
  bxbx[6] = gkyl_mat_get(&x_bxbx,6,0); 
  bxby[6] = gkyl_mat_get(&x_bxby,6,0); 
  bxbz[6] = gkyl_mat_get(&x_bxbz,6,0); 
  byby[6] = gkyl_mat_get(&x_byby,6,0); 
  bybz[6] = gkyl_mat_get(&x_bybz,6,0); 
  bzbz[6] = gkyl_mat_get(&x_bzbz,6,0); 
  bxbx[7] = gkyl_mat_get(&x_bxbx,7,0); 
  bxby[7] = gkyl_mat_get(&x_bxby,7,0); 
  bxbz[7] = gkyl_mat_get(&x_bxbz,7,0); 
  byby[7] = gkyl_mat_get(&x_byby,7,0); 
  bybz[7] = gkyl_mat_get(&x_bybz,7,0); 
  bzbz[7] = gkyl_mat_get(&x_bzbz,7,0); 
  const double *B_x = &em[24]; 
  const double *B_y = &em[32]; 
  const double *B_z = &em[40]; 
 
  int cell_avg = 0;
  if ((-0.3535533905932737*bxbx[7])+0.3535533905932737*bxbx[6]+0.3535533905932737*bxbx[5]+0.3535533905932737*bxbx[4]-0.3535533905932737*bxbx[3]-0.3535533905932737*bxbx[2]-0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*byby[7])+0.3535533905932737*byby[6]+0.3535533905932737*byby[5]+0.3535533905932737*byby[4]-0.3535533905932737*byby[3]-0.3535533905932737*byby[2]-0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bzbz[7])+0.3535533905932737*bzbz[6]+0.3535533905932737*bzbz[5]+0.3535533905932737*bzbz[4]-0.3535533905932737*bzbz[3]-0.3535533905932737*bzbz[2]-0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bxbx[7]-0.3535533905932737*bxbx[6]-0.3535533905932737*bxbx[5]+0.3535533905932737*bxbx[4]+0.3535533905932737*bxbx[3]-0.3535533905932737*bxbx[2]-0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*byby[7]-0.3535533905932737*byby[6]-0.3535533905932737*byby[5]+0.3535533905932737*byby[4]+0.3535533905932737*byby[3]-0.3535533905932737*byby[2]-0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bzbz[7]-0.3535533905932737*bzbz[6]-0.3535533905932737*bzbz[5]+0.3535533905932737*bzbz[4]+0.3535533905932737*bzbz[3]-0.3535533905932737*bzbz[2]-0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bxbx[7]-0.3535533905932737*bxbx[6]+0.3535533905932737*bxbx[5]-0.3535533905932737*bxbx[4]-0.3535533905932737*bxbx[3]+0.3535533905932737*bxbx[2]-0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*byby[7]-0.3535533905932737*byby[6]+0.3535533905932737*byby[5]-0.3535533905932737*byby[4]-0.3535533905932737*byby[3]+0.3535533905932737*byby[2]-0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bzbz[7]-0.3535533905932737*bzbz[6]+0.3535533905932737*bzbz[5]-0.3535533905932737*bzbz[4]-0.3535533905932737*bzbz[3]+0.3535533905932737*bzbz[2]-0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bxbx[7])+0.3535533905932737*bxbx[6]-0.3535533905932737*bxbx[5]-0.3535533905932737*bxbx[4]+0.3535533905932737*bxbx[3]+0.3535533905932737*bxbx[2]-0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*byby[7])+0.3535533905932737*byby[6]-0.3535533905932737*byby[5]-0.3535533905932737*byby[4]+0.3535533905932737*byby[3]+0.3535533905932737*byby[2]-0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bzbz[7])+0.3535533905932737*bzbz[6]-0.3535533905932737*bzbz[5]-0.3535533905932737*bzbz[4]+0.3535533905932737*bzbz[3]+0.3535533905932737*bzbz[2]-0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bxbx[7]+0.3535533905932737*bxbx[6]-0.3535533905932737*bxbx[5]-0.3535533905932737*bxbx[4]-0.3535533905932737*bxbx[3]-0.3535533905932737*bxbx[2]+0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*byby[7]+0.3535533905932737*byby[6]-0.3535533905932737*byby[5]-0.3535533905932737*byby[4]-0.3535533905932737*byby[3]-0.3535533905932737*byby[2]+0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bzbz[7]+0.3535533905932737*bzbz[6]-0.3535533905932737*bzbz[5]-0.3535533905932737*bzbz[4]-0.3535533905932737*bzbz[3]-0.3535533905932737*bzbz[2]+0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bxbx[7])-0.3535533905932737*bxbx[6]+0.3535533905932737*bxbx[5]-0.3535533905932737*bxbx[4]+0.3535533905932737*bxbx[3]-0.3535533905932737*bxbx[2]+0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*byby[7])-0.3535533905932737*byby[6]+0.3535533905932737*byby[5]-0.3535533905932737*byby[4]+0.3535533905932737*byby[3]-0.3535533905932737*byby[2]+0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bzbz[7])-0.3535533905932737*bzbz[6]+0.3535533905932737*bzbz[5]-0.3535533905932737*bzbz[4]+0.3535533905932737*bzbz[3]-0.3535533905932737*bzbz[2]+0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bxbx[7])-0.3535533905932737*bxbx[6]-0.3535533905932737*bxbx[5]+0.3535533905932737*bxbx[4]-0.3535533905932737*bxbx[3]+0.3535533905932737*bxbx[2]+0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*byby[7])-0.3535533905932737*byby[6]-0.3535533905932737*byby[5]+0.3535533905932737*byby[4]-0.3535533905932737*byby[3]+0.3535533905932737*byby[2]+0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.3535533905932737*bzbz[7])-0.3535533905932737*bzbz[6]-0.3535533905932737*bzbz[5]+0.3535533905932737*bzbz[4]-0.3535533905932737*bzbz[3]+0.3535533905932737*bzbz[2]+0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bxbx[7]+0.3535533905932737*bxbx[6]+0.3535533905932737*bxbx[5]+0.3535533905932737*bxbx[4]+0.3535533905932737*bxbx[3]+0.3535533905932737*bxbx[2]+0.3535533905932737*bxbx[1]+0.3535533905932737*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*byby[7]+0.3535533905932737*byby[6]+0.3535533905932737*byby[5]+0.3535533905932737*byby[4]+0.3535533905932737*byby[3]+0.3535533905932737*byby[2]+0.3535533905932737*byby[1]+0.3535533905932737*byby[0] < 0.0) cell_avg = 1; 
  if (0.3535533905932737*bzbz[7]+0.3535533905932737*bzbz[6]+0.3535533905932737*bzbz[5]+0.3535533905932737*bzbz[4]+0.3535533905932737*bzbz[3]+0.3535533905932737*bzbz[2]+0.3535533905932737*bzbz[1]+0.3535533905932737*bzbz[0] < 0.0) cell_avg = 1; 
  if (cell_avg || cell_avg_magB2[0]) { 
    bxbx[1] = 0.0; 
    bxby[1] = 0.0; 
    bxbz[1] = 0.0; 
    byby[1] = 0.0; 
    bybz[1] = 0.0; 
    bzbz[1] = 0.0; 
    bxbx[2] = 0.0; 
    bxby[2] = 0.0; 
    bxbz[2] = 0.0; 
    byby[2] = 0.0; 
    bybz[2] = 0.0; 
    bzbz[2] = 0.0; 
    bxbx[3] = 0.0; 
    bxby[3] = 0.0; 
    bxbz[3] = 0.0; 
    byby[3] = 0.0; 
    bybz[3] = 0.0; 
    bzbz[3] = 0.0; 
    bxbx[4] = 0.0; 
    bxby[4] = 0.0; 
    bxbz[4] = 0.0; 
    byby[4] = 0.0; 
    bybz[4] = 0.0; 
    bzbz[4] = 0.0; 
    bxbx[5] = 0.0; 
    bxby[5] = 0.0; 
    bxbz[5] = 0.0; 
    byby[5] = 0.0; 
    bybz[5] = 0.0; 
    bzbz[5] = 0.0; 
    bxbx[6] = 0.0; 
    bxby[6] = 0.0; 
    bxbz[6] = 0.0; 
    byby[6] = 0.0; 
    bybz[6] = 0.0; 
    bzbz[6] = 0.0; 
    bxbx[7] = 0.0; 
    bxby[7] = 0.0; 
    bxbz[7] = 0.0; 
    byby[7] = 0.0; 
    bybz[7] = 0.0; 
    bzbz[7] = 0.0; 
  // If bxbx, byby, or bzbz < 0.0 at the quadrature points, 
  // set cell_avg_magB2 to be true in case it was not true before. 
  cell_avg_magB2[0] = 1; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  ser_3x_p1_sqrt_with_sign(B_x, bxbx, bx); 
  ser_3x_p1_sqrt_with_sign(B_y, byby, by); 
  ser_3x_p1_sqrt_with_sign(B_z, bzbz, bz); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[4]; 
  double *bxbx_xl = &bvar_surf[8]; 
  double *bxbx_xr = &bvar_surf[12]; 
  double *bxby_xl = &bvar_surf[16]; 
  double *bxby_xr = &bvar_surf[20]; 
  double *bxbz_xl = &bvar_surf[24]; 
  double *bxbz_xr = &bvar_surf[28]; 
 
  bx_xl[0] = 0.7071067811865475*bx[0]-1.224744871391589*bx[1]; 
  bx_xl[1] = 0.7071067811865475*bx[2]-1.224744871391589*bx[4]; 
  bx_xl[2] = 0.7071067811865475*bx[3]-1.224744871391589*bx[5]; 
  bx_xl[3] = 0.7071067811865475*bx[6]-1.224744871391589*bx[7]; 
  bxbx_xl[0] = 0.7071067811865475*bxbx[0]-1.224744871391589*bxbx[1]; 
  bxbx_xl[1] = 0.7071067811865475*bxbx[2]-1.224744871391589*bxbx[4]; 
  bxbx_xl[2] = 0.7071067811865475*bxbx[3]-1.224744871391589*bxbx[5]; 
  bxbx_xl[3] = 0.7071067811865475*bxbx[6]-1.224744871391589*bxbx[7]; 
  bxby_xl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[1]; 
  bxby_xl[1] = 0.7071067811865475*bxby[2]-1.224744871391589*bxby[4]; 
  bxby_xl[2] = 0.7071067811865475*bxby[3]-1.224744871391589*bxby[5]; 
  bxby_xl[3] = 0.7071067811865475*bxby[6]-1.224744871391589*bxby[7]; 
  bxbz_xl[0] = 0.7071067811865475*bxbz[0]-1.224744871391589*bxbz[1]; 
  bxbz_xl[1] = 0.7071067811865475*bxbz[2]-1.224744871391589*bxbz[4]; 
  bxbz_xl[2] = 0.7071067811865475*bxbz[3]-1.224744871391589*bxbz[5]; 
  bxbz_xl[3] = 0.7071067811865475*bxbz[6]-1.224744871391589*bxbz[7]; 
 
  bx_xr[0] = 1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[1] = 1.224744871391589*bx[4]+0.7071067811865475*bx[2]; 
  bx_xr[2] = 1.224744871391589*bx[5]+0.7071067811865475*bx[3]; 
  bx_xr[3] = 1.224744871391589*bx[7]+0.7071067811865475*bx[6]; 
  bxbx_xr[0] = 1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xr[1] = 1.224744871391589*bxbx[4]+0.7071067811865475*bxbx[2]; 
  bxbx_xr[2] = 1.224744871391589*bxbx[5]+0.7071067811865475*bxbx[3]; 
  bxbx_xr[3] = 1.224744871391589*bxbx[7]+0.7071067811865475*bxbx[6]; 
  bxby_xr[0] = 1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xr[1] = 1.224744871391589*bxby[4]+0.7071067811865475*bxby[2]; 
  bxby_xr[2] = 1.224744871391589*bxby[5]+0.7071067811865475*bxby[3]; 
  bxby_xr[3] = 1.224744871391589*bxby[7]+0.7071067811865475*bxby[6]; 
  bxbz_xr[0] = 1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xr[1] = 1.224744871391589*bxbz[4]+0.7071067811865475*bxbz[2]; 
  bxbz_xr[2] = 1.224744871391589*bxbz[5]+0.7071067811865475*bxbz[3]; 
  bxbz_xr[3] = 1.224744871391589*bxbz[7]+0.7071067811865475*bxbz[6]; 
 
  double *by_yl = &bvar_surf[32]; 
  double *by_yr = &bvar_surf[36]; 
  double *bxby_yl = &bvar_surf[40]; 
  double *bxby_yr = &bvar_surf[44]; 
  double *byby_yl = &bvar_surf[48]; 
  double *byby_yr = &bvar_surf[52]; 
  double *bybz_yl = &bvar_surf[56]; 
  double *bybz_yr = &bvar_surf[60]; 
 
  by_yl[0] = 0.7071067811865475*by[0]-1.224744871391589*by[2]; 
  by_yl[1] = 0.7071067811865475*by[1]-1.224744871391589*by[4]; 
  by_yl[2] = 0.7071067811865475*by[3]-1.224744871391589*by[6]; 
  by_yl[3] = 0.7071067811865475*by[5]-1.224744871391589*by[7]; 
  bxby_yl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[2]; 
  bxby_yl[1] = 0.7071067811865475*bxby[1]-1.224744871391589*bxby[4]; 
  bxby_yl[2] = 0.7071067811865475*bxby[3]-1.224744871391589*bxby[6]; 
  bxby_yl[3] = 0.7071067811865475*bxby[5]-1.224744871391589*bxby[7]; 
  byby_yl[0] = 0.7071067811865475*byby[0]-1.224744871391589*byby[2]; 
  byby_yl[1] = 0.7071067811865475*byby[1]-1.224744871391589*byby[4]; 
  byby_yl[2] = 0.7071067811865475*byby[3]-1.224744871391589*byby[6]; 
  byby_yl[3] = 0.7071067811865475*byby[5]-1.224744871391589*byby[7]; 
  bybz_yl[0] = 0.7071067811865475*bybz[0]-1.224744871391589*bybz[2]; 
  bybz_yl[1] = 0.7071067811865475*bybz[1]-1.224744871391589*bybz[4]; 
  bybz_yl[2] = 0.7071067811865475*bybz[3]-1.224744871391589*bybz[6]; 
  bybz_yl[3] = 0.7071067811865475*bybz[5]-1.224744871391589*bybz[7]; 
 
  by_yr[0] = 1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yr[1] = 1.224744871391589*by[4]+0.7071067811865475*by[1]; 
  by_yr[2] = 1.224744871391589*by[6]+0.7071067811865475*by[3]; 
  by_yr[3] = 1.224744871391589*by[7]+0.7071067811865475*by[5]; 
  bxby_yr[0] = 1.224744871391589*bxby[2]+0.7071067811865475*bxby[0]; 
  bxby_yr[1] = 1.224744871391589*bxby[4]+0.7071067811865475*bxby[1]; 
  bxby_yr[2] = 1.224744871391589*bxby[6]+0.7071067811865475*bxby[3]; 
  bxby_yr[3] = 1.224744871391589*bxby[7]+0.7071067811865475*bxby[5]; 
  byby_yr[0] = 1.224744871391589*byby[2]+0.7071067811865475*byby[0]; 
  byby_yr[1] = 1.224744871391589*byby[4]+0.7071067811865475*byby[1]; 
  byby_yr[2] = 1.224744871391589*byby[6]+0.7071067811865475*byby[3]; 
  byby_yr[3] = 1.224744871391589*byby[7]+0.7071067811865475*byby[5]; 
  bybz_yr[0] = 1.224744871391589*bybz[2]+0.7071067811865475*bybz[0]; 
  bybz_yr[1] = 1.224744871391589*bybz[4]+0.7071067811865475*bybz[1]; 
  bybz_yr[2] = 1.224744871391589*bybz[6]+0.7071067811865475*bybz[3]; 
  bybz_yr[3] = 1.224744871391589*bybz[7]+0.7071067811865475*bybz[5]; 
 
  double *bz_zl = &bvar_surf[64]; 
  double *bz_zr = &bvar_surf[68]; 
  double *bxbz_zl = &bvar_surf[72]; 
  double *bxbz_zr = &bvar_surf[76]; 
  double *bybz_zl = &bvar_surf[80]; 
  double *bybz_zr = &bvar_surf[84]; 
  double *bzbz_zl = &bvar_surf[88]; 
  double *bzbz_zr = &bvar_surf[92]; 
 
  bz_zl[0] = 0.7071067811865475*bz[0]-1.224744871391589*bz[3]; 
  bz_zl[1] = 0.7071067811865475*bz[1]-1.224744871391589*bz[5]; 
  bz_zl[2] = 0.7071067811865475*bz[2]-1.224744871391589*bz[6]; 
  bz_zl[3] = 0.7071067811865475*bz[4]-1.224744871391589*bz[7]; 
  bxbz_zl[0] = 0.7071067811865475*bxbz[0]-1.224744871391589*bxbz[3]; 
  bxbz_zl[1] = 0.7071067811865475*bxbz[1]-1.224744871391589*bxbz[5]; 
  bxbz_zl[2] = 0.7071067811865475*bxbz[2]-1.224744871391589*bxbz[6]; 
  bxbz_zl[3] = 0.7071067811865475*bxbz[4]-1.224744871391589*bxbz[7]; 
  bybz_zl[0] = 0.7071067811865475*bybz[0]-1.224744871391589*bybz[3]; 
  bybz_zl[1] = 0.7071067811865475*bybz[1]-1.224744871391589*bybz[5]; 
  bybz_zl[2] = 0.7071067811865475*bybz[2]-1.224744871391589*bybz[6]; 
  bybz_zl[3] = 0.7071067811865475*bybz[4]-1.224744871391589*bybz[7]; 
  bzbz_zl[0] = 0.7071067811865475*bzbz[0]-1.224744871391589*bzbz[3]; 
  bzbz_zl[1] = 0.7071067811865475*bzbz[1]-1.224744871391589*bzbz[5]; 
  bzbz_zl[2] = 0.7071067811865475*bzbz[2]-1.224744871391589*bzbz[6]; 
  bzbz_zl[3] = 0.7071067811865475*bzbz[4]-1.224744871391589*bzbz[7]; 
 
  bz_zr[0] = 1.224744871391589*bz[3]+0.7071067811865475*bz[0]; 
  bz_zr[1] = 1.224744871391589*bz[5]+0.7071067811865475*bz[1]; 
  bz_zr[2] = 1.224744871391589*bz[6]+0.7071067811865475*bz[2]; 
  bz_zr[3] = 1.224744871391589*bz[7]+0.7071067811865475*bz[4]; 
  bxbz_zr[0] = 1.224744871391589*bxbz[3]+0.7071067811865475*bxbz[0]; 
  bxbz_zr[1] = 1.224744871391589*bxbz[5]+0.7071067811865475*bxbz[1]; 
  bxbz_zr[2] = 1.224744871391589*bxbz[6]+0.7071067811865475*bxbz[2]; 
  bxbz_zr[3] = 1.224744871391589*bxbz[7]+0.7071067811865475*bxbz[4]; 
  bybz_zr[0] = 1.224744871391589*bybz[3]+0.7071067811865475*bybz[0]; 
  bybz_zr[1] = 1.224744871391589*bybz[5]+0.7071067811865475*bybz[1]; 
  bybz_zr[2] = 1.224744871391589*bybz[6]+0.7071067811865475*bybz[2]; 
  bybz_zr[3] = 1.224744871391589*bybz[7]+0.7071067811865475*bybz[4]; 
  bzbz_zr[0] = 1.224744871391589*bzbz[3]+0.7071067811865475*bzbz[0]; 
  bzbz_zr[1] = 1.224744871391589*bzbz[5]+0.7071067811865475*bzbz[1]; 
  bzbz_zr[2] = 1.224744871391589*bzbz[6]+0.7071067811865475*bzbz[2]; 
  bzbz_zr[3] = 1.224744871391589*bzbz[7]+0.7071067811865475*bzbz[4]; 
 
} 
 
