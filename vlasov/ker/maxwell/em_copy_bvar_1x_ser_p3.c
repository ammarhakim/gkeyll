#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p3_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_1x_ser_p3(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
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
  if ((-0.5701294036773671*bxbx[3])+0.9681844646844029*bxbx[2]-1.054672281193885*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5701294036773671*byby[3])+0.9681844646844029*byby[2]-1.054672281193885*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5701294036773671*bzbz[3])+0.9681844646844029*bzbz[2]-1.054672281193885*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.7702725556588816*bxbx[3]-0.5164305132317774*bxbx[2]-0.4163900395009129*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.7702725556588816*byby[3]-0.5164305132317774*byby[2]-0.4163900395009129*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg = 1; 
  if (0.7702725556588816*bzbz[3]-0.5164305132317774*bzbz[2]-0.4163900395009129*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.7702725556588816*bxbx[3])-0.5164305132317774*bxbx[2]+0.4163900395009129*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.7702725556588816*byby[3])-0.5164305132317774*byby[2]+0.4163900395009129*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.7702725556588816*bzbz[3])-0.5164305132317774*bzbz[2]+0.4163900395009129*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.5701294036773671*bxbx[3]+0.9681844646844029*bxbx[2]+1.054672281193885*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.5701294036773671*byby[3]+0.9681844646844029*byby[2]+1.054672281193885*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg = 1; 
  if (0.5701294036773671*bzbz[3]+0.9681844646844029*bzbz[2]+1.054672281193885*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg = 1; 
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
  ser_1x_p3_sqrt_with_sign(B_x, bxbx, bx); 
  ser_1x_p3_sqrt_with_sign(B_y, byby, by); 
  ser_1x_p3_sqrt_with_sign(B_z, bzbz, bz); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[1]; 
  double *bxbx_xl = &bvar_surf[2]; 
  double *bxbx_xr = &bvar_surf[3]; 
  double *bxby_xl = &bvar_surf[4]; 
  double *bxby_xr = &bvar_surf[5]; 
  double *bxbz_xl = &bvar_surf[6]; 
  double *bxbz_xr = &bvar_surf[7]; 
 
  bx_xl[0] = (-1.870828693386971*bx[3])+1.58113883008419*bx[2]-1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[0] = 1.870828693386971*bx[3]+1.58113883008419*bx[2]+1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bxbx_xl[0] = (-1.870828693386971*bxbx[3])+1.58113883008419*bxbx[2]-1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xr[0] = 1.870828693386971*bxbx[3]+1.58113883008419*bxbx[2]+1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxby_xl[0] = (-1.870828693386971*bxby[3])+1.58113883008419*bxby[2]-1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xr[0] = 1.870828693386971*bxby[3]+1.58113883008419*bxby[2]+1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxbz_xl[0] = (-1.870828693386971*bxbz[3])+1.58113883008419*bxbz[2]-1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xr[0] = 1.870828693386971*bxbz[3]+1.58113883008419*bxbz[2]+1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
 
} 
 
