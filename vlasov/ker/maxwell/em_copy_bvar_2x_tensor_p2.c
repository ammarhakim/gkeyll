#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
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
  double *by = &bvar[9]; 
  double *bz = &bvar[18]; 
  double *bxbx = &bvar[27]; 
  double *bxby = &bvar[36]; 
  double *bxbz = &bvar[45]; 
  double *byby = &bvar[54]; 
  double *bybz = &bvar[63]; 
  double *bzbz = &bvar[72]; 
 
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
  bxbx[8] = gkyl_mat_get(&x_bxbx,8,0); 
  bxby[8] = gkyl_mat_get(&x_bxby,8,0); 
  bxbz[8] = gkyl_mat_get(&x_bxbz,8,0); 
  byby[8] = gkyl_mat_get(&x_byby,8,0); 
  bybz[8] = gkyl_mat_get(&x_bybz,8,0); 
  bzbz[8] = gkyl_mat_get(&x_bzbz,8,0); 
  const double *B_x = &em[27]; 
  const double *B_y = &em[36]; 
  const double *B_z = &em[45]; 
 
  int cell_avg = 0;
  if (0.4*bxbx[8]-0.5999999999999995*bxbx[7]-0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]+0.9*bxbx[3]-0.6708203932499369*bxbx[2]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]-0.5999999999999995*byby[7]-0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]+0.9*byby[3]-0.6708203932499369*byby[2]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]-0.5999999999999995*bzbz[7]-0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]+0.9*bzbz[3]-0.6708203932499369*bzbz[2]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])+0.75*bxbx[7]-0.5590169943749475*bxbx[5]+0.4472135954999579*bxbx[4]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])+0.75*byby[7]-0.5590169943749475*byby[5]+0.4472135954999579*byby[4]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])+0.75*bzbz[7]-0.5590169943749475*bzbz[5]+0.4472135954999579*bzbz[4]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]-0.5999999999999995*bxbx[7]+0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]-0.9*bxbx[3]+0.6708203932499369*bxbx[2]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]-0.5999999999999995*byby[7]+0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]-0.9*byby[3]+0.6708203932499369*byby[2]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]-0.5999999999999995*bzbz[7]+0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]-0.9*bzbz[3]+0.6708203932499369*bzbz[2]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])+0.75*bxbx[6]+0.4472135954999579*bxbx[5]-0.5590169943749475*bxbx[4]-0.6708203932499369*bxbx[2]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])+0.75*byby[6]+0.4472135954999579*byby[5]-0.5590169943749475*byby[4]-0.6708203932499369*byby[2]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])+0.75*bzbz[6]+0.4472135954999579*bzbz[5]-0.5590169943749475*bzbz[4]-0.6708203932499369*bzbz[2]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.625*bxbx[8]-0.5590169943749475*bxbx[5]-0.5590169943749475*bxbx[4]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.625*byby[8]-0.5590169943749475*byby[5]-0.5590169943749475*byby[4]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.625*bzbz[8]-0.5590169943749475*bzbz[5]-0.5590169943749475*bzbz[4]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])-0.75*bxbx[6]+0.4472135954999579*bxbx[5]-0.5590169943749475*bxbx[4]+0.6708203932499369*bxbx[2]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])-0.75*byby[6]+0.4472135954999579*byby[5]-0.5590169943749475*byby[4]+0.6708203932499369*byby[2]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])-0.75*bzbz[6]+0.4472135954999579*bzbz[5]-0.5590169943749475*bzbz[4]+0.6708203932499369*bzbz[2]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]+0.5999999999999995*bxbx[7]-0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]-0.9*bxbx[3]-0.6708203932499369*bxbx[2]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]+0.5999999999999995*byby[7]-0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]-0.9*byby[3]-0.6708203932499369*byby[2]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]+0.5999999999999995*bzbz[7]-0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]-0.9*bzbz[3]-0.6708203932499369*bzbz[2]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])-0.75*bxbx[7]-0.5590169943749475*bxbx[5]+0.4472135954999579*bxbx[4]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])-0.75*byby[7]-0.5590169943749475*byby[5]+0.4472135954999579*byby[4]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])-0.75*bzbz[7]-0.5590169943749475*bzbz[5]+0.4472135954999579*bzbz[4]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]+0.5999999999999995*bxbx[7]+0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]+0.9*bxbx[3]+0.6708203932499369*bxbx[2]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]+0.5999999999999995*byby[7]+0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]+0.9*byby[3]+0.6708203932499369*byby[2]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]+0.5999999999999995*bzbz[7]+0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]+0.9*bzbz[3]+0.6708203932499369*bzbz[2]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
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
    bxbx[8] = 0.0; 
    bxby[8] = 0.0; 
    bxbz[8] = 0.0; 
    byby[8] = 0.0; 
    bybz[8] = 0.0; 
    bzbz[8] = 0.0; 
  // If bxbx, byby, or bzbz < 0.0 at the quadrature points, 
  // set cell_avg_magB2 to be true in case it was not true before. 
  cell_avg_magB2[0] = 1; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  tensor_2x_p2_sqrt_with_sign(B_x, bxbx, bx); 
  tensor_2x_p2_sqrt_with_sign(B_y, byby, by); 
  tensor_2x_p2_sqrt_with_sign(B_z, bzbz, bz); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[3]; 
  double *bxbx_xl = &bvar_surf[6]; 
  double *bxbx_xr = &bvar_surf[9]; 
  double *bxby_xl = &bvar_surf[12]; 
  double *bxby_xr = &bvar_surf[15]; 
  double *bxbz_xl = &bvar_surf[18]; 
  double *bxbz_xr = &bvar_surf[21]; 
 
  bx_xl[0] = 1.58113883008419*bx[4]-1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xl[1] = 1.58113883008419*bx[6]-1.224744871391589*bx[3]+0.7071067811865475*bx[2]; 
  bx_xl[2] = 1.58113883008419*bx[8]-1.224744871391589*bx[7]+0.7071067811865475*bx[5]; 
  bxbx_xl[0] = 1.58113883008419*bxbx[4]-1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xl[1] = 1.58113883008419*bxbx[6]-1.224744871391589*bxbx[3]+0.7071067811865475*bxbx[2]; 
  bxbx_xl[2] = 1.58113883008419*bxbx[8]-1.224744871391589*bxbx[7]+0.7071067811865475*bxbx[5]; 
  bxby_xl[0] = 1.58113883008419*bxby[4]-1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xl[1] = 1.58113883008419*bxby[6]-1.224744871391589*bxby[3]+0.7071067811865475*bxby[2]; 
  bxby_xl[2] = 1.58113883008419*bxby[8]-1.224744871391589*bxby[7]+0.7071067811865475*bxby[5]; 
  bxbz_xl[0] = 1.58113883008419*bxbz[4]-1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xl[1] = 1.58113883008419*bxbz[6]-1.224744871391589*bxbz[3]+0.7071067811865475*bxbz[2]; 
  bxbz_xl[2] = 1.58113883008419*bxbz[8]-1.224744871391589*bxbz[7]+0.7071067811865475*bxbz[5]; 
 
  bx_xr[0] = 1.58113883008419*bx[4]+1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[1] = 1.58113883008419*bx[6]+1.224744871391589*bx[3]+0.7071067811865475*bx[2]; 
  bx_xr[2] = 1.58113883008419*bx[8]+1.224744871391589*bx[7]+0.7071067811865475*bx[5]; 
  bxbx_xr[0] = 1.58113883008419*bxbx[4]+1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxbx_xr[1] = 1.58113883008419*bxbx[6]+1.224744871391589*bxbx[3]+0.7071067811865475*bxbx[2]; 
  bxbx_xr[2] = 1.58113883008419*bxbx[8]+1.224744871391589*bxbx[7]+0.7071067811865475*bxbx[5]; 
  bxby_xr[0] = 1.58113883008419*bxby[4]+1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxby_xr[1] = 1.58113883008419*bxby[6]+1.224744871391589*bxby[3]+0.7071067811865475*bxby[2]; 
  bxby_xr[2] = 1.58113883008419*bxby[8]+1.224744871391589*bxby[7]+0.7071067811865475*bxby[5]; 
  bxbz_xr[0] = 1.58113883008419*bxbz[4]+1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
  bxbz_xr[1] = 1.58113883008419*bxbz[6]+1.224744871391589*bxbz[3]+0.7071067811865475*bxbz[2]; 
  bxbz_xr[2] = 1.58113883008419*bxbz[8]+1.224744871391589*bxbz[7]+0.7071067811865475*bxbz[5]; 
 
  double *by_yl = &bvar_surf[24]; 
  double *by_yr = &bvar_surf[27]; 
  double *bxby_yl = &bvar_surf[30]; 
  double *bxby_yr = &bvar_surf[33]; 
  double *byby_yl = &bvar_surf[36]; 
  double *byby_yr = &bvar_surf[39]; 
  double *bybz_yl = &bvar_surf[42]; 
  double *bybz_yr = &bvar_surf[45]; 
 
  by_yl[0] = 1.58113883008419*by[5]-1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yl[1] = 1.58113883008419*by[7]-1.224744871391589*by[3]+0.7071067811865475*by[1]; 
  by_yl[2] = 1.58113883008419*by[8]-1.224744871391589*by[6]+0.7071067811865475*by[4]; 
  bxby_yl[0] = 1.58113883008419*bxby[5]-1.224744871391589*bxby[2]+0.7071067811865475*bxby[0]; 
  bxby_yl[1] = 1.58113883008419*bxby[7]-1.224744871391589*bxby[3]+0.7071067811865475*bxby[1]; 
  bxby_yl[2] = 1.58113883008419*bxby[8]-1.224744871391589*bxby[6]+0.7071067811865475*bxby[4]; 
  byby_yl[0] = 1.58113883008419*byby[5]-1.224744871391589*byby[2]+0.7071067811865475*byby[0]; 
  byby_yl[1] = 1.58113883008419*byby[7]-1.224744871391589*byby[3]+0.7071067811865475*byby[1]; 
  byby_yl[2] = 1.58113883008419*byby[8]-1.224744871391589*byby[6]+0.7071067811865475*byby[4]; 
  bybz_yl[0] = 1.58113883008419*bybz[5]-1.224744871391589*bybz[2]+0.7071067811865475*bybz[0]; 
  bybz_yl[1] = 1.58113883008419*bybz[7]-1.224744871391589*bybz[3]+0.7071067811865475*bybz[1]; 
  bybz_yl[2] = 1.58113883008419*bybz[8]-1.224744871391589*bybz[6]+0.7071067811865475*bybz[4]; 
 
  by_yr[0] = 1.58113883008419*by[5]+1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yr[1] = 1.58113883008419*by[7]+1.224744871391589*by[3]+0.7071067811865475*by[1]; 
  by_yr[2] = 1.58113883008419*by[8]+1.224744871391589*by[6]+0.7071067811865475*by[4]; 
  bxby_yr[0] = 1.58113883008419*bxby[5]+1.224744871391589*bxby[2]+0.7071067811865475*bxby[0]; 
  bxby_yr[1] = 1.58113883008419*bxby[7]+1.224744871391589*bxby[3]+0.7071067811865475*bxby[1]; 
  bxby_yr[2] = 1.58113883008419*bxby[8]+1.224744871391589*bxby[6]+0.7071067811865475*bxby[4]; 
  byby_yr[0] = 1.58113883008419*byby[5]+1.224744871391589*byby[2]+0.7071067811865475*byby[0]; 
  byby_yr[1] = 1.58113883008419*byby[7]+1.224744871391589*byby[3]+0.7071067811865475*byby[1]; 
  byby_yr[2] = 1.58113883008419*byby[8]+1.224744871391589*byby[6]+0.7071067811865475*byby[4]; 
  bybz_yr[0] = 1.58113883008419*bybz[5]+1.224744871391589*bybz[2]+0.7071067811865475*bybz[0]; 
  bybz_yr[1] = 1.58113883008419*bybz[7]+1.224744871391589*bybz[3]+0.7071067811865475*bybz[1]; 
  bybz_yr[2] = 1.58113883008419*bybz[8]+1.224744871391589*bybz[6]+0.7071067811865475*bybz[4]; 
 
} 
 
