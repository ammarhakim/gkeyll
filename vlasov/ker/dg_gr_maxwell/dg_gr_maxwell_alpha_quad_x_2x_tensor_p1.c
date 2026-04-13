#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_tensor_p1(const double *w, const double *dxv, 
      const int theta_pole, const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *J_c, const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT A_plus_dQ, double* GKYL_RESTRICT A_minus_dQ,
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad) 
{ 
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[2]; 
  const double *h_xz_nodal = &h_ij_nodal[4]; 
  const double *h_yy_nodal = &h_ij_nodal[6]; 
  const double *h_yz_nodal = &h_ij_nodal[8]; 
  const double *h_zz_nodal = &h_ij_nodal[10]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[2]; 
  const double *shift_nodal_z = &shift_nodal[4]; 
  
  double Dx_con_l_n[2] = {0.0};
  double Dy_con_l_n[2] = {0.0};
  double Dz_con_l_n[2] = {0.0};
  double Bx_con_l_n[2] = {0.0};
  double By_con_l_n[2] = {0.0};
  double Bz_con_l_n[2] = {0.0};
  double Ex_l_n[2] = {0.0};
  double Ey_l_n[2] = {0.0};
  double Ez_l_n[2] = {0.0};
  double Hx_l_n[2] = {0.0};
  double Hy_l_n[2] = {0.0};
  double Hz_l_n[2] = {0.0};
  
  double Dx_con_r_n[2] = {0.0};
  double Dy_con_r_n[2] = {0.0};
  double Dz_con_r_n[2] = {0.0};
  double Bx_con_r_n[2] = {0.0};
  double By_con_r_n[2] = {0.0};
  double Bz_con_r_n[2] = {0.0};
  double Ex_r_n[2] = {0.0};
  double Ey_r_n[2] = {0.0};
  double Ez_r_n[2] = {0.0};
  double Hx_r_n[2] = {0.0};
  double Hy_r_n[2] = {0.0};
  double Hz_r_n[2] = {0.0};
  
  const double *Dx_l = &field_no_J_con_l[0]; 
  const double *Dy_l = &field_no_J_con_l[4]; 
  const double *Dz_l = &field_no_J_con_l[8]; 
  const double *Bx_l = &field_no_J_con_l[12]; 
  const double *By_l = &field_no_J_con_l[16]; 
  const double *Bz_l = &field_no_J_con_l[20]; 
  
  Dx_con_l_n[0] = -(0.8660254037844386*Dx_l[3])-0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[0] = -(0.8660254037844386*Dy_l[3])-0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[0] = -(0.8660254037844386*Dz_l[3])-0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[0] = -(0.8660254037844386*Bx_l[3])-0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[0] = -(0.8660254037844386*By_l[3])-0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[0] = -(0.8660254037844386*Bz_l[3])-0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  
  Dx_con_l_n[1] = 0.8660254037844386*Dx_l[3]+0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[1] = 0.8660254037844386*Dy_l[3]+0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[1] = 0.8660254037844386*Dz_l[3]+0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[1] = 0.8660254037844386*Bx_l[3]+0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[1] = 0.8660254037844386*By_l[3]+0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[1] = 0.8660254037844386*Bz_l[3]+0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  
  // If at a theta pole, zero out B^(theta), abd B^(theta) 
  if ( theta_pole ) { 
    Dy_con_l_n[0] = 0.0;
    By_con_l_n[0] = 0.0;
    Dy_con_l_n[1] = 0.0;
    By_con_l_n[1] = 0.0;
  } 
  
  for (int i=0; i<2; ++i) {
    Hx_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_l_n[i] + h_xy_nodal[i]*By_con_l_n[i] + h_xz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_l_n[i] - shift_nodal_z[i]*Dy_con_l_n[i]); 
    Hy_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_l_n[i] + h_yy_nodal[i]*By_con_l_n[i] + h_yz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_l_n[i] + shift_nodal_z[i]*Dx_con_l_n[i]); 
    Hz_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_l_n[i] + h_yz_nodal[i]*By_con_l_n[i] + h_zz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_l_n[i] - shift_nodal_y[i]*Dx_con_l_n[i]); 
    Ex_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_l_n[i] + h_xy_nodal[i]*Dy_con_l_n[i] + h_xz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_l_n[i] - shift_nodal_z[i]*By_con_l_n[i]); 
    Ey_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_l_n[i] + h_yy_nodal[i]*Dy_con_l_n[i] + h_yz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_l_n[i] + shift_nodal_z[i]*Bx_con_l_n[i]); 
    Ez_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_l_n[i] + h_yz_nodal[i]*Dy_con_l_n[i] + h_zz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_l_n[i] - shift_nodal_y[i]*Bx_con_l_n[i]); 
  }
  
  const double *Dx_r = &field_no_J_con_r[0]; 
  const double *Dy_r = &field_no_J_con_r[4]; 
  const double *Dz_r = &field_no_J_con_r[8]; 
  const double *Bx_r = &field_no_J_con_r[12]; 
  const double *By_r = &field_no_J_con_r[16]; 
  const double *Bz_r = &field_no_J_con_r[20]; 
  
  Dx_con_r_n[0] = 0.8660254037844386*Dx_r[3]-0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[0] = 0.8660254037844386*Dy_r[3]-0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[0] = 0.8660254037844386*Dz_r[3]-0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[0] = 0.8660254037844386*Bx_r[3]-0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[0] = 0.8660254037844386*By_r[3]-0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[0] = 0.8660254037844386*Bz_r[3]-0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  
  Dx_con_r_n[1] = -(0.8660254037844386*Dx_r[3])+0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[1] = -(0.8660254037844386*Dy_r[3])+0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[1] = -(0.8660254037844386*Dz_r[3])+0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[1] = -(0.8660254037844386*Bx_r[3])+0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[1] = -(0.8660254037844386*By_r[3])+0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[1] = -(0.8660254037844386*Bz_r[3])+0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  
  // If at a theta pole, zero out B^(theta), abd B^(theta) 
  if ( theta_pole ) { 
    Dy_con_r_n[0] = 0.0;
    By_con_r_n[0] = 0.0;
    Dy_con_r_n[1] = 0.0;
    By_con_r_n[1] = 0.0;
  } 
  
  for (int i=0; i<2; ++i) {
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
  flux_l_quad = &flux_l[2]; 
  flux_r_quad = &flux_r[2]; 
  flux_l_quad[0] = Hz_l_n[0]; 
  flux_r_quad[0] = Hz_r_n[0]; 
  flux_l_quad[1] = Hz_l_n[1]; 
  flux_r_quad[1] = Hz_r_n[1]; 
  flux_l_quad = &flux_l[4]; 
  flux_r_quad = &flux_r[4]; 
  flux_l_quad[0] = -(1.0*Hy_l_n[0]); 
  flux_r_quad[0] = -(1.0*Hy_r_n[0]); 
  flux_l_quad[1] = -(1.0*Hy_l_n[1]); 
  flux_r_quad[1] = -(1.0*Hy_r_n[1]); 
  flux_l_quad = &flux_l[6]; 
  flux_r_quad = &flux_r[6]; 
  flux_l_quad[0] = 0.0; 
  flux_r_quad[0] = 0.0; 
  flux_l_quad[1] = 0.0; 
  flux_r_quad[1] = 0.0; 
  flux_l_quad = &flux_l[8]; 
  flux_r_quad = &flux_r[8]; 
  flux_l_quad[0] = -(1.0*Ez_l_n[0]); 
  flux_r_quad[0] = -(1.0*Ez_r_n[0]); 
  flux_l_quad[1] = -(1.0*Ez_l_n[1]); 
  flux_r_quad[1] = -(1.0*Ez_r_n[1]); 
  flux_l_quad = &flux_l[10]; 
  flux_r_quad = &flux_r[10]; 
  flux_l_quad[0] = Ey_l_n[0]; 
  flux_r_quad[0] = Ey_r_n[0]; 
  flux_l_quad[1] = Ey_l_n[1]; 
  flux_r_quad[1] = Ey_r_n[1]; 
  
  double lambda_1[2] = {0.0};
  double lambda_2[2] = {0.0};
  double lambda_3[2] = {0.0};
  for (int i=0; i<2; ++i) {
    lambda_2[i] = -shift_nodal_x[i] + lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_3[i] = -shift_nodal_x[i] - lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_2[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_3[i] ));
  }
  // If at the theta pole, dU = dQ = 0, so A_plus_dQ = A_minus_dQ = 0
  if (theta_pole == 0) {    double dQ_n[6] = {0.0};
    double lambda_plus_n[6] = {0.0};
    double lambda_minus_n[6] = {0.0};
    double A_plus_dQ_n[6] = {0.0};
    double A_minus_dQ_n[6] = {0.0};
    for (int i=0; i<2; ++i) {
      dQ_n[0] = J_c[i]*(Dx_con_r_n[i] -  Dx_con_l_n[i]);
      dQ_n[1] = J_c[i]*(Dy_con_r_n[i] -  Dy_con_l_n[i]);
      dQ_n[2] = J_c[i]*(Dz_con_r_n[i] -  Dz_con_l_n[i]);
      dQ_n[3] = J_c[i]*(Bx_con_r_n[i] -  Bx_con_l_n[i]);
      dQ_n[4] = J_c[i]*(By_con_r_n[i] -  By_con_l_n[i]);
      dQ_n[5] = J_c[i]*(Bz_con_r_n[i] -  Bz_con_l_n[i]);
      lambda_plus_n[0] = (lambda_1[i] > 0.0) ? lambda_1[i] : 0.0;
      lambda_plus_n[1] = (lambda_1[i] > 0.0) ? lambda_1[i] : 0.0;
      lambda_plus_n[2] = (lambda_2[i] > 0.0) ? lambda_2[i] : 0.0;
      lambda_plus_n[3] = (lambda_2[i] > 0.0) ? lambda_2[i] : 0.0;
      lambda_plus_n[4] = (lambda_3[i] > 0.0) ? lambda_3[i] : 0.0;
      lambda_plus_n[5] = (lambda_3[i] > 0.0) ? lambda_3[i] : 0.0;
      lambda_minus_n[0] = (lambda_1[i] < 0.0) ? lambda_1[i] : 0.0;
      lambda_minus_n[1] = (lambda_1[i] < 0.0) ? lambda_1[i] : 0.0;
      lambda_minus_n[2] = (lambda_2[i] < 0.0) ? lambda_2[i] : 0.0;
      lambda_minus_n[3] = (lambda_2[i] < 0.0) ? lambda_2[i] : 0.0;
      lambda_minus_n[4] = (lambda_3[i] < 0.0) ? lambda_3[i] : 0.0;
      lambda_minus_n[5] = (lambda_3[i] < 0.0) ? lambda_3[i] : 0.0;
      A_dQ_x_calc(lapse_nodal[i], shift_nodal_x[i], shift_nodal_y[i], shift_nodal_z[i], 
                      h_xx_nodal[i], h_xy_nodal[i], h_xz_nodal[i], h_yy_nodal[i], h_yz_nodal[i], h_zz_nodal[i], 
                      J_c[i], lambda_plus_n, dQ_n, A_plus_dQ_n);
      A_plus_dQ[i + 0*2] = A_plus_dQ_n[0]; 
      A_plus_dQ[i + 1*2] = A_plus_dQ_n[1]; 
      A_plus_dQ[i + 2*2] = A_plus_dQ_n[2]; 
      A_plus_dQ[i + 3*2] = A_plus_dQ_n[3]; 
      A_plus_dQ[i + 4*2] = A_plus_dQ_n[4]; 
      A_plus_dQ[i + 5*2] = A_plus_dQ_n[5]; 
      A_dQ_x_calc(lapse_nodal[i], shift_nodal_x[i], shift_nodal_y[i], shift_nodal_z[i], 
                      h_xx_nodal[i], h_xy_nodal[i], h_xz_nodal[i], h_yy_nodal[i], h_yz_nodal[i], h_zz_nodal[i], 
                      J_c[i], lambda_minus_n, dQ_n, A_minus_dQ_n);
      A_minus_dQ[i + 0*2] = A_minus_dQ_n[0]; 
      A_minus_dQ[i + 1*2] = A_minus_dQ_n[1]; 
      A_minus_dQ[i + 2*2] = A_minus_dQ_n[2]; 
      A_minus_dQ[i + 3*2] = A_minus_dQ_n[3]; 
      A_minus_dQ[i + 4*2] = A_minus_dQ_n[4]; 
      A_minus_dQ[i + 5*2] = A_minus_dQ_n[5]; 
    }
  }

} 
