#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_1x_tensor_p1(const double *w, const double *dxv, 
      const int geom_edge, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *J_c, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad) 
{ 
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[1]; 
  const double *h_xz_nodal = &h_ij_nodal[2]; 
  const double *h_yy_nodal = &h_ij_nodal[3]; 
  const double *h_yz_nodal = &h_ij_nodal[4]; 
  const double *h_zz_nodal = &h_ij_nodal[5]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[1]; 
  const double *shift_nodal_z = &shift_nodal[2]; 
  
  double Dx_con_l_n[1] = {0.0};
  double Dy_con_l_n[1] = {0.0};
  double Dz_con_l_n[1] = {0.0};
  double Bx_con_l_n[1] = {0.0};
  double By_con_l_n[1] = {0.0};
  double Bz_con_l_n[1] = {0.0};
  double Ex_l_n[1] = {0.0};
  double Ey_l_n[1] = {0.0};
  double Ez_l_n[1] = {0.0};
  double Hx_l_n[1] = {0.0};
  double Hy_l_n[1] = {0.0};
  double Hz_l_n[1] = {0.0};
  
  double Dx_con_r_n[1] = {0.0};
  double Dy_con_r_n[1] = {0.0};
  double Dz_con_r_n[1] = {0.0};
  double Bx_con_r_n[1] = {0.0};
  double By_con_r_n[1] = {0.0};
  double Bz_con_r_n[1] = {0.0};
  double Ex_r_n[1] = {0.0};
  double Ey_r_n[1] = {0.0};
  double Ez_r_n[1] = {0.0};
  double Hx_r_n[1] = {0.0};
  double Hy_r_n[1] = {0.0};
  double Hz_r_n[1] = {0.0};
  
  const double *Dx_l = &field_no_J_con_l[0]; 
  const double *Dy_l = &field_no_J_con_l[2]; 
  const double *Dz_l = &field_no_J_con_l[4]; 
  const double *Bx_l = &field_no_J_con_l[6]; 
  const double *By_l = &field_no_J_con_l[8]; 
  const double *Bz_l = &field_no_J_con_l[10]; 
  
  Dx_con_l_n[0] = 1.224744871391589*Dx_l[1]+0.7071067811865475*Dx_l[0];
  Dy_con_l_n[0] = 1.224744871391589*Dy_l[1]+0.7071067811865475*Dy_l[0];
  Dz_con_l_n[0] = 1.224744871391589*Dz_l[1]+0.7071067811865475*Dz_l[0];
  Bx_con_l_n[0] = 1.224744871391589*Bx_l[1]+0.7071067811865475*Bx_l[0];
  By_con_l_n[0] = 1.224744871391589*By_l[1]+0.7071067811865475*By_l[0];
  Bz_con_l_n[0] = 1.224744871391589*Bz_l[1]+0.7071067811865475*Bz_l[0];
  
  for (int i=0; i<1; ++i) {
    Hx_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_l_n[i] + h_xy_nodal[i]*By_con_l_n[i] + h_xz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_l_n[i] - shift_nodal_z[i]*Dy_con_l_n[i]); 
    Hy_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_l_n[i] + h_yy_nodal[i]*By_con_l_n[i] + h_yz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_l_n[i] + shift_nodal_z[i]*Dx_con_l_n[i]); 
    Hz_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_l_n[i] + h_yz_nodal[i]*By_con_l_n[i] + h_zz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_l_n[i] - shift_nodal_y[i]*Dx_con_l_n[i]); 
    Ex_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_l_n[i] + h_xy_nodal[i]*Dy_con_l_n[i] + h_xz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_l_n[i] - shift_nodal_z[i]*By_con_l_n[i]); 
    Ey_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_l_n[i] + h_yy_nodal[i]*Dy_con_l_n[i] + h_yz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_l_n[i] + shift_nodal_z[i]*Bx_con_l_n[i]); 
    Ez_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_l_n[i] + h_yz_nodal[i]*Dy_con_l_n[i] + h_zz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_l_n[i] - shift_nodal_y[i]*Bx_con_l_n[i]); 
  }
  
  const double *Dx_r = &field_no_J_con_r[0]; 
  const double *Dy_r = &field_no_J_con_r[2]; 
  const double *Dz_r = &field_no_J_con_r[4]; 
  const double *Bx_r = &field_no_J_con_r[6]; 
  const double *By_r = &field_no_J_con_r[8]; 
  const double *Bz_r = &field_no_J_con_r[10]; 
  
  Dx_con_r_n[0] = 0.7071067811865475*Dx_r[0]-1.224744871391589*Dx_r[1];
  Dy_con_r_n[0] = 0.7071067811865475*Dy_r[0]-1.224744871391589*Dy_r[1];
  Dz_con_r_n[0] = 0.7071067811865475*Dz_r[0]-1.224744871391589*Dz_r[1];
  Bx_con_r_n[0] = 0.7071067811865475*Bx_r[0]-1.224744871391589*Bx_r[1];
  By_con_r_n[0] = 0.7071067811865475*By_r[0]-1.224744871391589*By_r[1];
  Bz_con_r_n[0] = 0.7071067811865475*Bz_r[0]-1.224744871391589*Bz_r[1];
  
  for (int i=0; i<1; ++i) {
    Hx_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_r_n[i] + h_xy_nodal[i]*By_con_r_n[i] + h_xz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_r_n[i] - shift_nodal_z[i]*Dy_con_r_n[i]); 
    Hy_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_r_n[i] + h_yy_nodal[i]*By_con_r_n[i] + h_yz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_r_n[i] + shift_nodal_z[i]*Dx_con_r_n[i]); 
    Hz_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_r_n[i] + h_yz_nodal[i]*By_con_r_n[i] + h_zz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_r_n[i] - shift_nodal_y[i]*Dx_con_r_n[i]); 
    Ex_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_r_n[i] + h_xy_nodal[i]*Dy_con_r_n[i] + h_xz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_r_n[i] - shift_nodal_z[i]*By_con_r_n[i]); 
    Ey_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_r_n[i] + h_yy_nodal[i]*Dy_con_r_n[i] + h_yz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_r_n[i] + shift_nodal_z[i]*Bx_con_r_n[i]); 
    Ez_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_r_n[i] + h_yz_nodal[i]*Dy_con_r_n[i] + h_zz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_r_n[i] - shift_nodal_y[i]*Bx_con_r_n[i]); 
  }
  
  double *flux_l_quad; 
  double *flux_r_quad; 
  flux_l_quad = &flux_l[0]; 
  flux_r_quad = &flux_r[0]; 
  flux_l_quad[0] = 0.0; 
  flux_r_quad[0] = 0.0; 
  flux_l_quad = &flux_l[1]; 
  flux_r_quad = &flux_r[1]; 
  flux_l_quad[0] = Hz_l_n[0]; 
  flux_r_quad[0] = Hz_r_n[0]; 
  flux_l_quad = &flux_l[2]; 
  flux_r_quad = &flux_r[2]; 
  flux_l_quad[0] = -(1.0*Hy_l_n[0]); 
  flux_r_quad[0] = -(1.0*Hy_r_n[0]); 
  flux_l_quad = &flux_l[3]; 
  flux_r_quad = &flux_r[3]; 
  flux_l_quad[0] = 0.0; 
  flux_r_quad[0] = 0.0; 
  flux_l_quad = &flux_l[4]; 
  flux_r_quad = &flux_r[4]; 
  flux_l_quad[0] = -(1.0*Ez_l_n[0]); 
  flux_r_quad[0] = -(1.0*Ez_r_n[0]); 
  flux_l_quad = &flux_l[5]; 
  flux_r_quad = &flux_r[5]; 
  flux_l_quad[0] = Ey_l_n[0]; 
  flux_r_quad[0] = Ey_r_n[0]; 
  
  for (int i=0; i<1; ++i) {
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( -shift_nodal_x[i] - lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i]) );
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( -shift_nodal_x[i] + lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i]) );
  }

} 
