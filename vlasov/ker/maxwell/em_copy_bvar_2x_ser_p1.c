#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_2x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
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
  double *by = &bvar[4]; 
  double *bz = &bvar[8]; 
  double *bxbx = &bvar[12]; 
  double *bxby = &bvar[16]; 
  double *bxbz = &bvar[20]; 
  double *byby = &bvar[24]; 
  double *bybz = &bvar[28]; 
  double *bzbz = &bvar[32]; 
 
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
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  int cell_avg = 0;
  if (0.5*bxbx[3]-0.5*bxbx[2]-0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.5*byby[3]-0.5*byby[2]-0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.5*bzbz[3]-0.5*bzbz[2]-0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[3])+0.5*bxbx[2]-0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[3])+0.5*byby[2]-0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[3])+0.5*bzbz[2]-0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[3])-0.5*bxbx[2]+0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[3])-0.5*byby[2]+0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[3])-0.5*bzbz[2]+0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.5*bxbx[3]+0.5*bxbx[2]+0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.5*byby[3]+0.5*byby[2]+0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.5*bzbz[3]+0.5*bzbz[2]+0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
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
  // If bxbx, byby, or bzbz < 0.0 at the quadrature points, 
  // set cell_avg_magB2 to be true in case it was not true before. 
  cell_avg_magB2[0] = 1; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  ser_2x_p1_sqrt_with_sign(B_x, bxbx, bx); 
  ser_2x_p1_sqrt_with_sign(B_y, byby, by); 
  ser_2x_p1_sqrt_with_sign(B_z, bzbz, bz); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[2]; 
  double *bxbx_xl = &bvar_surf[4]; 
  double *bxbx_xr = &bvar_surf[6]; 
  double *bxby_xl = &bvar_surf[8]; 
  double *bxby_xr = &bvar_surf[10]; 
  double *bxbz_xl = &bvar_surf[12]; 
  double *bxbz_xr = &bvar_surf[14]; 
 
  bx_xl[0] = 0.7071067811865475*bx[0]-1.224744871391589*bx[1]; 
  bx_xl[1] = 0.7071067811865475*bx[2]-1.224744871391589*bx[3]; 
  bxbx_xl[0] = 0.7071067811865475*bxbx[0]-1.224744871391589*bxbx[1]; 
  bxbx_xl[1] = 0.7071067811865475*bxbx[2]-1.224744871391589*bxbx[3]; 
  bxby_xl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[1]; 
  bxby_xl[1] = 0.7071067811865475*bxby[2]-1.224744871391589*bxby[3]; 
  bxbz_xl[0] = 0.7071067811865475*bxbz[0]-1.224744871391589*bxbz[1]; 
  bxbz_xl[1] = 0.7071067811865475*bxbz[2]-1.224744871391589*bxbz[3]; 
 
  bx_xr[0] = 1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[1] = 1.224744871391589*bx[3]+0.7071067811865475*bx[2]; 
  bxbx_xr[0] = 1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xr[1] = 1.224744871391589*bxbx[3]+0.7071067811865475*bxbx[2]; 
  bxby_xr[0] = 1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xr[1] = 1.224744871391589*bxby[3]+0.7071067811865475*bxby[2]; 
  bxbz_xr[0] = 1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xr[1] = 1.224744871391589*bxbz[3]+0.7071067811865475*bxbz[2]; 
 
  double *by_yl = &bvar_surf[16]; 
  double *by_yr = &bvar_surf[18]; 
  double *bxby_yl = &bvar_surf[20]; 
  double *bxby_yr = &bvar_surf[22]; 
  double *byby_yl = &bvar_surf[24]; 
  double *byby_yr = &bvar_surf[26]; 
  double *bybz_yl = &bvar_surf[28]; 
  double *bybz_yr = &bvar_surf[30]; 
 
  by_yl[0] = 0.7071067811865475*by[0]-1.224744871391589*by[2]; 
  by_yl[1] = 0.7071067811865475*by[1]-1.224744871391589*by[3]; 
  bxby_yl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[2]; 
  bxby_yl[1] = 0.7071067811865475*bxby[1]-1.224744871391589*bxby[3]; 
  byby_yl[0] = 0.7071067811865475*byby[0]-1.224744871391589*byby[2]; 
  byby_yl[1] = 0.7071067811865475*byby[1]-1.224744871391589*byby[3]; 
  bybz_yl[0] = 0.7071067811865475*bybz[0]-1.224744871391589*bybz[2]; 
  bybz_yl[1] = 0.7071067811865475*bybz[1]-1.224744871391589*bybz[3]; 
 
  by_yr[0] = 1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yr[1] = 1.224744871391589*by[3]+0.7071067811865475*by[1]; 
  bxby_yr[0] = 1.224744871391589*bxby[2]+0.7071067811865475*bxby[0]; 
  bxby_yr[1] = 1.224744871391589*bxby[3]+0.7071067811865475*bxby[1]; 
  byby_yr[0] = 1.224744871391589*byby[2]+0.7071067811865475*byby[0]; 
  byby_yr[1] = 1.224744871391589*byby[3]+0.7071067811865475*byby[1]; 
  bybz_yr[0] = 1.224744871391589*bybz[2]+0.7071067811865475*bybz[0]; 
  bybz_yr[1] = 1.224744871391589*bybz[3]+0.7071067811865475*bybz[1]; 
 
} 
 
