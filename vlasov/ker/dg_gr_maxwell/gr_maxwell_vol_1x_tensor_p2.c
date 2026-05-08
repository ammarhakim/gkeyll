#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p2(const double *w, const double *dx, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out) 
{ 
  double dx0 = 2.0/dx[0]; 
  
  double JDx_con_n[3] = {0.0};
  double JDy_con_n[3] = {0.0};
  double JDz_con_n[3] = {0.0};
  double JBx_con_n[3] = {0.0};
  double JBy_con_n[3] = {0.0};
  double JBz_con_n[3] = {0.0};
  double Ex_n[3] = {0.0};
  double Ey_n[3] = {0.0};
  double Ez_n[3] = {0.0};
  double Hx_n[3] = {0.0};
  double Hy_n[3] = {0.0};
  double Hz_n[3] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[3]; 
  const double *JDz = &fields[6]; 
  const double *JBx = &fields[9]; 
  const double *JBy = &fields[12]; 
  const double *JBz = &fields[15]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[3]; 
  double *outJDz = &out[6]; 
  double *outJBx = &out[9]; 
  double *outJBy = &out[12]; 
  double *outJBz = &out[15]; 
  
  JDx_con_n[0] = 0.6324555320336759*JDx[2]-0.9486832980505137*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[0] = 0.6324555320336759*JDy[2]-0.9486832980505137*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[0] = 0.6324555320336759*JDz[2]-0.9486832980505137*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[0] = 0.6324555320336759*JBx[2]-0.9486832980505137*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[0] = 0.6324555320336759*JBy[2]-0.9486832980505137*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[0] = 0.6324555320336759*JBz[2]-0.9486832980505137*JBz[1]+0.7071067811865475*JBz[0];
  
  JDx_con_n[1] = 0.7071067811865475*JDx[0]-0.7905694150420947*JDx[2];
  JDy_con_n[1] = 0.7071067811865475*JDy[0]-0.7905694150420947*JDy[2];
  JDz_con_n[1] = 0.7071067811865475*JDz[0]-0.7905694150420947*JDz[2];
  JBx_con_n[1] = 0.7071067811865475*JBx[0]-0.7905694150420947*JBx[2];
  JBy_con_n[1] = 0.7071067811865475*JBy[0]-0.7905694150420947*JBy[2];
  JBz_con_n[1] = 0.7071067811865475*JBz[0]-0.7905694150420947*JBz[2];
  
  JDx_con_n[2] = 0.6324555320336759*JDx[2]+0.9486832980505137*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[2] = 0.6324555320336759*JDy[2]+0.9486832980505137*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[2] = 0.6324555320336759*JDz[2]+0.9486832980505137*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[2] = 0.6324555320336759*JBx[2]+0.9486832980505137*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[2] = 0.6324555320336759*JBy[2]+0.9486832980505137*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[2] = 0.6324555320336759*JBz[2]+0.9486832980505137*JBz[1]+0.7071067811865475*JBz[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[3]; 
  const double *h_xz_nodal = &h_ij_nodal[6]; 
  const double *h_yy_nodal = &h_ij_nodal[9]; 
  const double *h_yz_nodal = &h_ij_nodal[12]; 
  const double *h_zz_nodal = &h_ij_nodal[15]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[3]; 
  const double *shift_nodal_z = &shift_nodal[6]; 
  
  for (int i=0; i<3; ++i) {
    Hx_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JBx_con_n[i] + h_xy_nodal[i]*JBy_con_n[i] + h_xz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_y[i]*JDz_con_n[i] - shift_nodal_z[i]*JDy_con_n[i]); 
    Hy_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JBx_con_n[i] + h_yy_nodal[i]*JBy_con_n[i] + h_yz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( - shift_nodal_x[i]*JDz_con_n[i] + shift_nodal_z[i]*JDx_con_n[i]); 
    Hz_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JBx_con_n[i] + h_yz_nodal[i]*JBy_con_n[i] + h_zz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_x[i]*JDy_con_n[i] - shift_nodal_y[i]*JDx_con_n[i]); 
    Ex_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_n[i] + h_xy_nodal[i]*JDy_con_n[i] + h_xz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_n[i] - shift_nodal_z[i]*JBy_con_n[i]); 
    Ey_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_n[i] + h_yy_nodal[i]*JDy_con_n[i] + h_yz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_n[i] + shift_nodal_z[i]*JBx_con_n[i]); 
    Ez_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_n[i] + h_yz_nodal[i]*JDy_con_n[i] + h_zz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_n[i] - shift_nodal_y[i]*JBx_con_n[i]); 
  }
  

  outJDy[1] += 0.6804138174397718*Hz_n[2]*dx0+1.0886621079036354*Hz_n[1]*dx0+0.6804138174397718*Hz_n[0]*dx0; 
  outJDy[2] += 2.041241452319315*Hz_n[2]*dx0-2.041241452319315*Hz_n[0]*dx0; 

  outJDz[1] += -(0.6804138174397718*Hy_n[2]*dx0)-1.0886621079036354*Hy_n[1]*dx0-0.6804138174397718*Hy_n[0]*dx0; 
  outJDz[2] += 2.041241452319315*Hy_n[0]*dx0-2.041241452319315*Hy_n[2]*dx0; 


  outJBy[1] += -(0.6804138174397718*Ez_n[2]*dx0)-1.0886621079036354*Ez_n[1]*dx0-0.6804138174397718*Ez_n[0]*dx0; 
  outJBy[2] += 2.041241452319315*Ez_n[0]*dx0-2.041241452319315*Ez_n[2]*dx0; 

  outJBz[1] += 0.6804138174397718*Ey_n[2]*dx0+1.0886621079036354*Ey_n[1]*dx0+0.6804138174397718*Ey_n[0]*dx0; 
  outJBz[2] += 2.041241452319315*Ey_n[2]*dx0-2.041241452319315*Ey_n[0]*dx0; 

} 
