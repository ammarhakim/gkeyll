#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_vol_2x_ser_p1(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out) 
{ 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  
  double JDx_con_n[4] = {0.0};
  double JDy_con_n[4] = {0.0};
  double JDz_con_n[4] = {0.0};
  double JBx_con_n[4] = {0.0};
  double JBy_con_n[4] = {0.0};
  double JBz_con_n[4] = {0.0};
  double Ex_n[4] = {0.0};
  double Ey_n[4] = {0.0};
  double Ez_n[4] = {0.0};
  double Hx_n[4] = {0.0};
  double Hy_n[4] = {0.0};
  double Hz_n[4] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[4]; 
  const double *JDz = &fields[8]; 
  const double *JBx = &fields[12]; 
  const double *JBy = &fields[16]; 
  const double *JBz = &fields[20]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[4]; 
  double *outJDz = &out[8]; 
  double *outJBx = &out[12]; 
  double *outJBy = &out[16]; 
  double *outJBz = &out[20]; 
  
  JDx_con_n[0] = 0.5*JDx[3]-0.5*JDx[2]-0.5*JDx[1]+0.5*JDx[0];
  JDy_con_n[0] = 0.5*JDy[3]-0.5*JDy[2]-0.5*JDy[1]+0.5*JDy[0];
  JDz_con_n[0] = 0.5*JDz[3]-0.5*JDz[2]-0.5*JDz[1]+0.5*JDz[0];
  JBx_con_n[0] = 0.5*JBx[3]-0.5*JBx[2]-0.5*JBx[1]+0.5*JBx[0];
  JBy_con_n[0] = 0.5*JBy[3]-0.5*JBy[2]-0.5*JBy[1]+0.5*JBy[0];
  JBz_con_n[0] = 0.5*JBz[3]-0.5*JBz[2]-0.5*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[1] = -(0.5*JDx[3])+0.5*JDx[2]-0.5*JDx[1]+0.5*JDx[0];
  JDy_con_n[1] = -(0.5*JDy[3])+0.5*JDy[2]-0.5*JDy[1]+0.5*JDy[0];
  JDz_con_n[1] = -(0.5*JDz[3])+0.5*JDz[2]-0.5*JDz[1]+0.5*JDz[0];
  JBx_con_n[1] = -(0.5*JBx[3])+0.5*JBx[2]-0.5*JBx[1]+0.5*JBx[0];
  JBy_con_n[1] = -(0.5*JBy[3])+0.5*JBy[2]-0.5*JBy[1]+0.5*JBy[0];
  JBz_con_n[1] = -(0.5*JBz[3])+0.5*JBz[2]-0.5*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[2] = -(0.5*JDx[3])-0.5*JDx[2]+0.5*JDx[1]+0.5*JDx[0];
  JDy_con_n[2] = -(0.5*JDy[3])-0.5*JDy[2]+0.5*JDy[1]+0.5*JDy[0];
  JDz_con_n[2] = -(0.5*JDz[3])-0.5*JDz[2]+0.5*JDz[1]+0.5*JDz[0];
  JBx_con_n[2] = -(0.5*JBx[3])-0.5*JBx[2]+0.5*JBx[1]+0.5*JBx[0];
  JBy_con_n[2] = -(0.5*JBy[3])-0.5*JBy[2]+0.5*JBy[1]+0.5*JBy[0];
  JBz_con_n[2] = -(0.5*JBz[3])-0.5*JBz[2]+0.5*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[3] = 0.5*JDx[3]+0.5*JDx[2]+0.5*JDx[1]+0.5*JDx[0];
  JDy_con_n[3] = 0.5*JDy[3]+0.5*JDy[2]+0.5*JDy[1]+0.5*JDy[0];
  JDz_con_n[3] = 0.5*JDz[3]+0.5*JDz[2]+0.5*JDz[1]+0.5*JDz[0];
  JBx_con_n[3] = 0.5*JBx[3]+0.5*JBx[2]+0.5*JBx[1]+0.5*JBx[0];
  JBy_con_n[3] = 0.5*JBy[3]+0.5*JBy[2]+0.5*JBy[1]+0.5*JBy[0];
  JBz_con_n[3] = 0.5*JBz[3]+0.5*JBz[2]+0.5*JBz[1]+0.5*JBz[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[4]; 
  const double *h_xz_nodal = &h_ij_nodal[8]; 
  const double *h_yy_nodal = &h_ij_nodal[12]; 
  const double *h_yz_nodal = &h_ij_nodal[16]; 
  const double *h_zz_nodal = &h_ij_nodal[20]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[4]; 
  const double *shift_nodal_z = &shift_nodal[8]; 
  
  for (int i=0; i<4; ++i) {
    Hx_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JBx_con_n[i] + h_xy_nodal[i]*JBy_con_n[i] + h_xz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_y[i]*JDz_con_n[i] - shift_nodal_z[i]*JDy_con_n[i]); 
    Hy_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JBx_con_n[i] + h_yy_nodal[i]*JBy_con_n[i] + h_yz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( - shift_nodal_x[i]*JDz_con_n[i] + shift_nodal_z[i]*JDx_con_n[i]); 
    Hz_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JBx_con_n[i] + h_yz_nodal[i]*JBy_con_n[i] + h_zz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_x[i]*JDy_con_n[i] - shift_nodal_y[i]*JDx_con_n[i]); 
    Ex_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_n[i] + h_xy_nodal[i]*JDy_con_n[i] + h_xz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_n[i] - shift_nodal_z[i]*JBy_con_n[i]); 
    Ey_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_n[i] + h_yy_nodal[i]*JDy_con_n[i] + h_yz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_n[i] + shift_nodal_z[i]*JBx_con_n[i]); 
    Ez_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_n[i] + h_yz_nodal[i]*JDy_con_n[i] + h_zz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_n[i] - shift_nodal_y[i]*JBx_con_n[i]); 
  }
  
  outJDx[2] += -(0.8660254037844386*Hz_n[3]*dx1)-0.8660254037844386*Hz_n[2]*dx1-0.8660254037844386*Hz_n[1]*dx1-0.8660254037844386*Hz_n[0]*dx1; 
  outJDx[3] += -(0.8660254037844386*Hz_n[3]*dx1)-0.8660254037844386*Hz_n[2]*dx1+0.8660254037844386*Hz_n[1]*dx1+0.8660254037844386*Hz_n[0]*dx1; 

  outJDy[1] += 0.8660254037844386*Hz_n[3]*dx0+0.8660254037844386*Hz_n[2]*dx0+0.8660254037844386*Hz_n[1]*dx0+0.8660254037844386*Hz_n[0]*dx0; 
  outJDy[3] += 0.8660254037844386*Hz_n[3]*dx0-0.8660254037844386*Hz_n[2]*dx0+0.8660254037844386*Hz_n[1]*dx0-0.8660254037844386*Hz_n[0]*dx0; 

  outJDz[1] += -(0.8660254037844386*Hy_n[3]*dx0)-0.8660254037844386*Hy_n[2]*dx0-0.8660254037844386*Hy_n[1]*dx0-0.8660254037844386*Hy_n[0]*dx0; 
  outJDz[2] += 0.8660254037844386*Hx_n[3]*dx1+0.8660254037844386*Hx_n[2]*dx1+0.8660254037844386*Hx_n[1]*dx1+0.8660254037844386*Hx_n[0]*dx1; 
  outJDz[3] += 0.8660254037844386*Hx_n[3]*dx1+0.8660254037844386*Hx_n[2]*dx1-0.8660254037844386*Hx_n[1]*dx1-0.8660254037844386*Hx_n[0]*dx1-0.8660254037844386*Hy_n[3]*dx0+0.8660254037844386*Hy_n[2]*dx0-0.8660254037844386*Hy_n[1]*dx0+0.8660254037844386*Hy_n[0]*dx0; 

  outJBx[2] += 0.8660254037844386*Ez_n[3]*dx1+0.8660254037844386*Ez_n[2]*dx1+0.8660254037844386*Ez_n[1]*dx1+0.8660254037844386*Ez_n[0]*dx1; 
  outJBx[3] += 0.8660254037844386*Ez_n[3]*dx1+0.8660254037844386*Ez_n[2]*dx1-0.8660254037844386*Ez_n[1]*dx1-0.8660254037844386*Ez_n[0]*dx1; 

  outJBy[1] += -(0.8660254037844386*Ez_n[3]*dx0)-0.8660254037844386*Ez_n[2]*dx0-0.8660254037844386*Ez_n[1]*dx0-0.8660254037844386*Ez_n[0]*dx0; 
  outJBy[3] += -(0.8660254037844386*Ez_n[3]*dx0)+0.8660254037844386*Ez_n[2]*dx0-0.8660254037844386*Ez_n[1]*dx0+0.8660254037844386*Ez_n[0]*dx0; 

  outJBz[1] += 0.8660254037844386*Ey_n[3]*dx0+0.8660254037844386*Ey_n[2]*dx0+0.8660254037844386*Ey_n[1]*dx0+0.8660254037844386*Ey_n[0]*dx0; 
  outJBz[2] += -(0.8660254037844386*Ex_n[3]*dx1)-0.8660254037844386*Ex_n[2]*dx1-0.8660254037844386*Ex_n[1]*dx1-0.8660254037844386*Ex_n[0]*dx1; 
  outJBz[3] += -(0.8660254037844386*Ex_n[3]*dx1)-0.8660254037844386*Ex_n[2]*dx1+0.8660254037844386*Ex_n[1]*dx1+0.8660254037844386*Ex_n[0]*dx1+0.8660254037844386*Ey_n[3]*dx0-0.8660254037844386*Ey_n[2]*dx0+0.8660254037844386*Ey_n[1]*dx0-0.8660254037844386*Ey_n[0]*dx0; 

} 
