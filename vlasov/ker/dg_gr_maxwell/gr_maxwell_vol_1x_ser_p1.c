#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *J_c, const double *fields_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx0 = 2.0/dx[0]; 
  
  double Dx_con_n[2] = {0.0};
  double Dy_con_n[2] = {0.0};
  double Dz_con_n[2] = {0.0};
  double Bx_con_n[2] = {0.0};
  double By_con_n[2] = {0.0};
  double Bz_con_n[2] = {0.0};
  double Ex_n[2] = {0.0};
  double Ey_n[2] = {0.0};
  double Ez_n[2] = {0.0};
  double Hx_n[2] = {0.0};
  double Hy_n[2] = {0.0};
  double Hz_n[2] = {0.0};
  
  const double *Dx = &fields_no_J[0]; 
  const double *Dy = &fields_no_J[2]; 
  const double *Dz = &fields_no_J[4]; 
  const double *Bx = &fields_no_J[6]; 
  const double *By = &fields_no_J[8]; 
  const double *Bz = &fields_no_J[10]; 
  
  double *outDx = &out[0]; 
  double *outDy = &out[2]; 
  double *outDz = &out[4]; 
  double *outBx = &out[6]; 
  double *outBy = &out[8]; 
  double *outBz = &out[10]; 
  
  Dx_con_n[0] = 0.7071067811865475*Dx[0]-0.7071067811865475*Dx[1];
  Dy_con_n[0] = 0.7071067811865475*Dy[0]-0.7071067811865475*Dy[1];
  Dz_con_n[0] = 0.7071067811865475*Dz[0]-0.7071067811865475*Dz[1];
  Bx_con_n[0] = 0.7071067811865475*Bx[0]-0.7071067811865475*Bx[1];
  By_con_n[0] = 0.7071067811865475*By[0]-0.7071067811865475*By[1];
  Bz_con_n[0] = 0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1];
  
  Dx_con_n[1] = 0.7071067811865475*Dx[1]+0.7071067811865475*Dx[0];
  Dy_con_n[1] = 0.7071067811865475*Dy[1]+0.7071067811865475*Dy[0];
  Dz_con_n[1] = 0.7071067811865475*Dz[1]+0.7071067811865475*Dz[0];
  Bx_con_n[1] = 0.7071067811865475*Bx[1]+0.7071067811865475*Bx[0];
  By_con_n[1] = 0.7071067811865475*By[1]+0.7071067811865475*By[0];
  Bz_con_n[1] = 0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[2]; 
  const double *h_xz_nodal = &h_ij_nodal[4]; 
  const double *h_yy_nodal = &h_ij_nodal[6]; 
  const double *h_yz_nodal = &h_ij_nodal[8]; 
  const double *h_zz_nodal = &h_ij_nodal[10]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[2]; 
  const double *shift_nodal_z = &shift_nodal[4]; 
  
  for (int i=0; i<2; ++i) {
    Hx_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_n[i] + h_xy_nodal[i]*By_con_n[i] + h_xz_nodal[i]*Bz_con_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_n[i] - shift_nodal_z[i]*Dy_con_n[i]); 
    Hy_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_n[i] + h_yy_nodal[i]*By_con_n[i] + h_yz_nodal[i]*Bz_con_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_n[i] + shift_nodal_z[i]*Dx_con_n[i]); 
    Hz_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_n[i] + h_yz_nodal[i]*By_con_n[i] + h_zz_nodal[i]*Bz_con_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_n[i] - shift_nodal_y[i]*Dx_con_n[i]); 
    Ex_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_n[i] + h_xy_nodal[i]*Dy_con_n[i] + h_xz_nodal[i]*Dz_con_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_n[i] - shift_nodal_z[i]*By_con_n[i]); 
    Ey_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_n[i] + h_yy_nodal[i]*Dy_con_n[i] + h_yz_nodal[i]*Dz_con_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_n[i] + shift_nodal_z[i]*Bx_con_n[i]); 
    Ez_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_n[i] + h_yz_nodal[i]*Dy_con_n[i] + h_zz_nodal[i]*Dz_con_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_n[i] - shift_nodal_y[i]*Bx_con_n[i]); 
  }
  

  outDy[1] += 1.224744871391589*Hz_n[1]*dx0+1.224744871391589*Hz_n[0]*dx0; 

  outDz[1] += -(1.224744871391589*Hy_n[1]*dx0)-1.224744871391589*Hy_n[0]*dx0; 


  outBy[1] += -(1.224744871391589*Ez_n[1]*dx0)-1.224744871391589*Ez_n[0]*dx0; 

  outBz[1] += 1.224744871391589*Ey_n[1]*dx0+1.224744871391589*Ey_n[0]*dx0; 

} 
