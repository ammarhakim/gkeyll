#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const double *w, const double *dxv, 
      const int geom_edge, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *J_c, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad) 
{ 
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[3]; 
  const double *h_xz_nodal = &h_ij_nodal[6]; 
  const double *h_yy_nodal = &h_ij_nodal[9]; 
  const double *h_yz_nodal = &h_ij_nodal[12]; 
  const double *h_zz_nodal = &h_ij_nodal[15]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[3]; 
  const double *shift_nodal_z = &shift_nodal[6]; 
  
  double Dx_con_l_n[3] = {0.0};
  double Dy_con_l_n[3] = {0.0};
  double Dz_con_l_n[3] = {0.0};
  double Bx_con_l_n[3] = {0.0};
  double By_con_l_n[3] = {0.0};
  double Bz_con_l_n[3] = {0.0};
  double Ex_l_n[3] = {0.0};
  double Ey_l_n[3] = {0.0};
  double Ez_l_n[3] = {0.0};
  double Hx_l_n[3] = {0.0};
  double Hy_l_n[3] = {0.0};
  double Hz_l_n[3] = {0.0};
  
  double Dx_con_r_n[3] = {0.0};
  double Dy_con_r_n[3] = {0.0};
  double Dz_con_r_n[3] = {0.0};
  double Bx_con_r_n[3] = {0.0};
  double By_con_r_n[3] = {0.0};
  double Bz_con_r_n[3] = {0.0};
  double Ex_r_n[3] = {0.0};
  double Ey_r_n[3] = {0.0};
  double Ez_r_n[3] = {0.0};
  double Hx_r_n[3] = {0.0};
  double Hy_r_n[3] = {0.0};
  double Hz_r_n[3] = {0.0};
  
  const double *Dx_l = &field_no_J_con_l[0]; 
  const double *Dy_l = &field_no_J_con_l[8]; 
  const double *Dz_l = &field_no_J_con_l[16]; 
  const double *Bx_l = &field_no_J_con_l[24]; 
  const double *By_l = &field_no_J_con_l[32]; 
  const double *Bz_l = &field_no_J_con_l[40]; 
  
  Dx_con_l_n[0] = 0.7745966692414833*Dx_l[7]-1.5*Dx_l[6]+0.4472135954999579*Dx_l[5]+1.118033988749895*Dx_l[4]-1.1618950038622249*Dx_l[3]-0.6708203932499369*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[0] = 0.7745966692414833*Dy_l[7]-1.5*Dy_l[6]+0.4472135954999579*Dy_l[5]+1.118033988749895*Dy_l[4]-1.1618950038622249*Dy_l[3]-0.6708203932499369*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[0] = 0.7745966692414833*Dz_l[7]-1.5*Dz_l[6]+0.4472135954999579*Dz_l[5]+1.118033988749895*Dz_l[4]-1.1618950038622249*Dz_l[3]-0.6708203932499369*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[0] = 0.7745966692414833*Bx_l[7]-1.5*Bx_l[6]+0.4472135954999579*Bx_l[5]+1.118033988749895*Bx_l[4]-1.1618950038622249*Bx_l[3]-0.6708203932499369*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[0] = 0.7745966692414833*By_l[7]-1.5*By_l[6]+0.4472135954999579*By_l[5]+1.118033988749895*By_l[4]-1.1618950038622249*By_l[3]-0.6708203932499369*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[0] = 0.7745966692414833*Bz_l[7]-1.5*Bz_l[6]+0.4472135954999579*Bz_l[5]+1.118033988749895*Bz_l[4]-1.1618950038622249*Bz_l[3]-0.6708203932499369*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  
  Dx_con_l_n[1] = -(0.9682458365518543*Dx_l[7])-0.5590169943749475*Dx_l[5]+1.118033988749895*Dx_l[4]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[1] = -(0.9682458365518543*Dy_l[7])-0.5590169943749475*Dy_l[5]+1.118033988749895*Dy_l[4]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[1] = -(0.9682458365518543*Dz_l[7])-0.5590169943749475*Dz_l[5]+1.118033988749895*Dz_l[4]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[1] = -(0.9682458365518543*Bx_l[7])-0.5590169943749475*Bx_l[5]+1.118033988749895*Bx_l[4]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[1] = -(0.9682458365518543*By_l[7])-0.5590169943749475*By_l[5]+1.118033988749895*By_l[4]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[1] = -(0.9682458365518543*Bz_l[7])-0.5590169943749475*Bz_l[5]+1.118033988749895*Bz_l[4]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  
  Dx_con_l_n[2] = 0.7745966692414833*Dx_l[7]+1.5*Dx_l[6]+0.4472135954999579*Dx_l[5]+1.118033988749895*Dx_l[4]+1.1618950038622249*Dx_l[3]+0.6708203932499369*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[2] = 0.7745966692414833*Dy_l[7]+1.5*Dy_l[6]+0.4472135954999579*Dy_l[5]+1.118033988749895*Dy_l[4]+1.1618950038622249*Dy_l[3]+0.6708203932499369*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[2] = 0.7745966692414833*Dz_l[7]+1.5*Dz_l[6]+0.4472135954999579*Dz_l[5]+1.118033988749895*Dz_l[4]+1.1618950038622249*Dz_l[3]+0.6708203932499369*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[2] = 0.7745966692414833*Bx_l[7]+1.5*Bx_l[6]+0.4472135954999579*Bx_l[5]+1.118033988749895*Bx_l[4]+1.1618950038622249*Bx_l[3]+0.6708203932499369*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[2] = 0.7745966692414833*By_l[7]+1.5*By_l[6]+0.4472135954999579*By_l[5]+1.118033988749895*By_l[4]+1.1618950038622249*By_l[3]+0.6708203932499369*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[2] = 0.7745966692414833*Bz_l[7]+1.5*Bz_l[6]+0.4472135954999579*Bz_l[5]+1.118033988749895*Bz_l[4]+1.1618950038622249*Bz_l[3]+0.6708203932499369*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  
  for (int i=0; i<3; ++i) {
    Hx_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_l_n[i] + h_xy_nodal[i]*By_con_l_n[i] + h_xz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_l_n[i] - shift_nodal_z[i]*Dy_con_l_n[i]); 
    Hy_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_l_n[i] + h_yy_nodal[i]*By_con_l_n[i] + h_yz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_l_n[i] + shift_nodal_z[i]*Dx_con_l_n[i]); 
    Hz_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_l_n[i] + h_yz_nodal[i]*By_con_l_n[i] + h_zz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_l_n[i] - shift_nodal_y[i]*Dx_con_l_n[i]); 
    Ex_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_l_n[i] + h_xy_nodal[i]*Dy_con_l_n[i] + h_xz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_l_n[i] - shift_nodal_z[i]*By_con_l_n[i]); 
    Ey_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_l_n[i] + h_yy_nodal[i]*Dy_con_l_n[i] + h_yz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_l_n[i] + shift_nodal_z[i]*Bx_con_l_n[i]); 
    Ez_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_l_n[i] + h_yz_nodal[i]*Dy_con_l_n[i] + h_zz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_l_n[i] - shift_nodal_y[i]*Bx_con_l_n[i]); 
  }
  
  const double *Dx_r = &field_no_J_con_r[0]; 
  const double *Dy_r = &field_no_J_con_r[8]; 
  const double *Dz_r = &field_no_J_con_r[16]; 
  const double *Bx_r = &field_no_J_con_r[24]; 
  const double *By_r = &field_no_J_con_r[32]; 
  const double *Bz_r = &field_no_J_con_r[40]; 
  
  Dx_con_r_n[0] = -(0.7745966692414833*Dx_r[7])-1.5*Dx_r[6]+0.4472135954999579*Dx_r[5]+1.118033988749895*Dx_r[4]+1.1618950038622249*Dx_r[3]-0.6708203932499369*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[0] = -(0.7745966692414833*Dy_r[7])-1.5*Dy_r[6]+0.4472135954999579*Dy_r[5]+1.118033988749895*Dy_r[4]+1.1618950038622249*Dy_r[3]-0.6708203932499369*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[0] = -(0.7745966692414833*Dz_r[7])-1.5*Dz_r[6]+0.4472135954999579*Dz_r[5]+1.118033988749895*Dz_r[4]+1.1618950038622249*Dz_r[3]-0.6708203932499369*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[0] = -(0.7745966692414833*Bx_r[7])-1.5*Bx_r[6]+0.4472135954999579*Bx_r[5]+1.118033988749895*Bx_r[4]+1.1618950038622249*Bx_r[3]-0.6708203932499369*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[0] = -(0.7745966692414833*By_r[7])-1.5*By_r[6]+0.4472135954999579*By_r[5]+1.118033988749895*By_r[4]+1.1618950038622249*By_r[3]-0.6708203932499369*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[0] = -(0.7745966692414833*Bz_r[7])-1.5*Bz_r[6]+0.4472135954999579*Bz_r[5]+1.118033988749895*Bz_r[4]+1.1618950038622249*Bz_r[3]-0.6708203932499369*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  
  Dx_con_r_n[1] = 0.9682458365518543*Dx_r[7]-0.5590169943749475*Dx_r[5]+1.118033988749895*Dx_r[4]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[1] = 0.9682458365518543*Dy_r[7]-0.5590169943749475*Dy_r[5]+1.118033988749895*Dy_r[4]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[1] = 0.9682458365518543*Dz_r[7]-0.5590169943749475*Dz_r[5]+1.118033988749895*Dz_r[4]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[1] = 0.9682458365518543*Bx_r[7]-0.5590169943749475*Bx_r[5]+1.118033988749895*Bx_r[4]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[1] = 0.9682458365518543*By_r[7]-0.5590169943749475*By_r[5]+1.118033988749895*By_r[4]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[1] = 0.9682458365518543*Bz_r[7]-0.5590169943749475*Bz_r[5]+1.118033988749895*Bz_r[4]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  
  Dx_con_r_n[2] = -(0.7745966692414833*Dx_r[7])+1.5*Dx_r[6]+0.4472135954999579*Dx_r[5]+1.118033988749895*Dx_r[4]-1.1618950038622249*Dx_r[3]+0.6708203932499369*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[2] = -(0.7745966692414833*Dy_r[7])+1.5*Dy_r[6]+0.4472135954999579*Dy_r[5]+1.118033988749895*Dy_r[4]-1.1618950038622249*Dy_r[3]+0.6708203932499369*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[2] = -(0.7745966692414833*Dz_r[7])+1.5*Dz_r[6]+0.4472135954999579*Dz_r[5]+1.118033988749895*Dz_r[4]-1.1618950038622249*Dz_r[3]+0.6708203932499369*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[2] = -(0.7745966692414833*Bx_r[7])+1.5*Bx_r[6]+0.4472135954999579*Bx_r[5]+1.118033988749895*Bx_r[4]-1.1618950038622249*Bx_r[3]+0.6708203932499369*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[2] = -(0.7745966692414833*By_r[7])+1.5*By_r[6]+0.4472135954999579*By_r[5]+1.118033988749895*By_r[4]-1.1618950038622249*By_r[3]+0.6708203932499369*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[2] = -(0.7745966692414833*Bz_r[7])+1.5*Bz_r[6]+0.4472135954999579*Bz_r[5]+1.118033988749895*Bz_r[4]-1.1618950038622249*Bz_r[3]+0.6708203932499369*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  
  for (int i=0; i<3; ++i) {
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
  flux_l_quad[1] = 0.0; 
  flux_r_quad[1] = 0.0; 
  flux_l_quad[2] = 0.0; 
  flux_r_quad[2] = 0.0; 
  flux_l_quad = &flux_l[3]; 
  flux_r_quad = &flux_r[3]; 
  flux_l_quad[0] = 0.16666666666666666*(0.6*(5.0*Hz_l_n[2]-10.0*Hz_l_n[1]+5.0*Hz_l_n[0])-0.7745966692414833*(3.872983346207417*Hz_l_n[2]-3.872983346207417*Hz_l_n[0])+6.0*Hz_l_n[1]); 
  flux_r_quad[0] = 0.16666666666666666*(0.6*(5.0*Hz_r_n[2]-10.0*Hz_r_n[1]+5.0*Hz_r_n[0])-0.7745966692414833*(3.872983346207417*Hz_r_n[2]-3.872983346207417*Hz_r_n[0])+6.0*Hz_r_n[1]); 
  flux_l_quad[1] = Hz_l_n[1]; 
  flux_r_quad[1] = Hz_r_n[1]; 
  flux_l_quad[2] = 0.16666666666666666*(0.6*(5.0*Hz_l_n[2]-10.0*Hz_l_n[1]+5.0*Hz_l_n[0])+0.7745966692414833*(3.872983346207417*Hz_l_n[2]-3.872983346207417*Hz_l_n[0])+6.0*Hz_l_n[1]); 
  flux_r_quad[2] = 0.16666666666666666*(0.6*(5.0*Hz_r_n[2]-10.0*Hz_r_n[1]+5.0*Hz_r_n[0])+0.7745966692414833*(3.872983346207417*Hz_r_n[2]-3.872983346207417*Hz_r_n[0])+6.0*Hz_r_n[1]); 
  flux_l_quad = &flux_l[6]; 
  flux_r_quad = &flux_r[6]; 
  flux_l_quad[0] = -(0.16666666666666666*(0.6*(5.0*Hy_l_n[2]-10.0*Hy_l_n[1]+5.0*Hy_l_n[0])-0.7745966692414833*(3.872983346207417*Hy_l_n[2]-3.872983346207417*Hy_l_n[0])+6.0*Hy_l_n[1])); 
  flux_r_quad[0] = -(0.16666666666666666*(0.6*(5.0*Hy_r_n[2]-10.0*Hy_r_n[1]+5.0*Hy_r_n[0])-0.7745966692414833*(3.872983346207417*Hy_r_n[2]-3.872983346207417*Hy_r_n[0])+6.0*Hy_r_n[1])); 
  flux_l_quad[1] = -(1.0*Hy_l_n[1]); 
  flux_r_quad[1] = -(1.0*Hy_r_n[1]); 
  flux_l_quad[2] = -(0.16666666666666666*(0.6*(5.0*Hy_l_n[2]-10.0*Hy_l_n[1]+5.0*Hy_l_n[0])+0.7745966692414833*(3.872983346207417*Hy_l_n[2]-3.872983346207417*Hy_l_n[0])+6.0*Hy_l_n[1])); 
  flux_r_quad[2] = -(0.16666666666666666*(0.6*(5.0*Hy_r_n[2]-10.0*Hy_r_n[1]+5.0*Hy_r_n[0])+0.7745966692414833*(3.872983346207417*Hy_r_n[2]-3.872983346207417*Hy_r_n[0])+6.0*Hy_r_n[1])); 
  flux_l_quad = &flux_l[9]; 
  flux_r_quad = &flux_r[9]; 
  flux_l_quad[0] = 0.0; 
  flux_r_quad[0] = 0.0; 
  flux_l_quad[1] = 0.0; 
  flux_r_quad[1] = 0.0; 
  flux_l_quad[2] = 0.0; 
  flux_r_quad[2] = 0.0; 
  flux_l_quad = &flux_l[12]; 
  flux_r_quad = &flux_r[12]; 
  flux_l_quad[0] = -(0.16666666666666666*(0.6*(5.0*Ez_l_n[2]-10.0*Ez_l_n[1]+5.0*Ez_l_n[0])-0.7745966692414833*(3.872983346207417*Ez_l_n[2]-3.872983346207417*Ez_l_n[0])+6.0*Ez_l_n[1])); 
  flux_r_quad[0] = -(0.16666666666666666*(0.6*(5.0*Ez_r_n[2]-10.0*Ez_r_n[1]+5.0*Ez_r_n[0])-0.7745966692414833*(3.872983346207417*Ez_r_n[2]-3.872983346207417*Ez_r_n[0])+6.0*Ez_r_n[1])); 
  flux_l_quad[1] = -(1.0*Ez_l_n[1]); 
  flux_r_quad[1] = -(1.0*Ez_r_n[1]); 
  flux_l_quad[2] = -(0.16666666666666666*(0.6*(5.0*Ez_l_n[2]-10.0*Ez_l_n[1]+5.0*Ez_l_n[0])+0.7745966692414833*(3.872983346207417*Ez_l_n[2]-3.872983346207417*Ez_l_n[0])+6.0*Ez_l_n[1])); 
  flux_r_quad[2] = -(0.16666666666666666*(0.6*(5.0*Ez_r_n[2]-10.0*Ez_r_n[1]+5.0*Ez_r_n[0])+0.7745966692414833*(3.872983346207417*Ez_r_n[2]-3.872983346207417*Ez_r_n[0])+6.0*Ez_r_n[1])); 
  flux_l_quad = &flux_l[15]; 
  flux_r_quad = &flux_r[15]; 
  flux_l_quad[0] = 0.16666666666666666*(0.6*(5.0*Ey_l_n[2]-10.0*Ey_l_n[1]+5.0*Ey_l_n[0])-0.7745966692414833*(3.872983346207417*Ey_l_n[2]-3.872983346207417*Ey_l_n[0])+6.0*Ey_l_n[1]); 
  flux_r_quad[0] = 0.16666666666666666*(0.6*(5.0*Ey_r_n[2]-10.0*Ey_r_n[1]+5.0*Ey_r_n[0])-0.7745966692414833*(3.872983346207417*Ey_r_n[2]-3.872983346207417*Ey_r_n[0])+6.0*Ey_r_n[1]); 
  flux_l_quad[1] = Ey_l_n[1]; 
  flux_r_quad[1] = Ey_r_n[1]; 
  flux_l_quad[2] = 0.16666666666666666*(0.6*(5.0*Ey_l_n[2]-10.0*Ey_l_n[1]+5.0*Ey_l_n[0])+0.7745966692414833*(3.872983346207417*Ey_l_n[2]-3.872983346207417*Ey_l_n[0])+6.0*Ey_l_n[1]); 
  flux_r_quad[2] = 0.16666666666666666*(0.6*(5.0*Ey_r_n[2]-10.0*Ey_r_n[1]+5.0*Ey_r_n[0])+0.7745966692414833*(3.872983346207417*Ey_r_n[2]-3.872983346207417*Ey_r_n[0])+6.0*Ey_r_n[1]); 
  
  for (int i=0; i<3; ++i) {
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( -shift_nodal_x[i] - lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i]) );
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( -shift_nodal_x[i] + lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i]) );
  }

} 
