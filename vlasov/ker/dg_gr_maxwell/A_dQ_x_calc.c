#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void A_dQ_x_calc(const double lapse, const double shift_x, const double shift_y, const double shift_z, 
      const double h_xx, const double h_xy, const double h_xz, const double h_yy, 
      const double h_yz, const double h_zz, const double J_c, const double *lambda, const double *dQ, double* GKYL_RESTRICT A_dQ) 
{ 
  A_dQ[0] = dQ[0]*lambda[0];
  A_dQ[1] = dQ[1]*(lambda[3]/2.0+lambda[5]/2.0)+(dQ[5]*sqrt(h_yy*h_zz)*(lambda[3]-lambda[5]))/(h_yy*2.0)-(dQ[3]*h_xz*lapse*(J_c*lambda[1]*shift_x*-2.0+J_c*lambda[3]*shift_x+J_c*lambda[5]*shift_x+lambda[3]*lapse*sqrt(h_yy*h_zz)-lambda[5]*lapse*sqrt(h_yy*h_zz)))/((J_c*J_c)*(shift_x*shift_x)*2.0-h_yy*h_zz*(lapse*lapse)*2.0);
  A_dQ[2] = dQ[2]*(lambda[2]/2.0+lambda[4]/2.0)-(dQ[4]*sqrt(h_yy*h_zz)*(lambda[2]-lambda[4]))/(h_zz*2.0)-(dQ[0]*h_xz*lapse*(h_yy*h_zz*lambda[0]*lapse*-2.0+h_yy*h_zz*lambda[2]*lapse+h_yy*h_zz*lambda[4]*lapse+J_c*lambda[2]*shift_x*sqrt(h_yy*h_zz)-J_c*lambda[4]*shift_x*sqrt(h_yy*h_zz)))/((J_c*J_c)*h_zz*(shift_x*shift_x)*2.0-h_yy*(h_zz*h_zz)*(lapse*lapse)*2.0);
  A_dQ[3] = dQ[3]*lambda[1];
  A_dQ[4] = dQ[4]*(lambda[2]/2.0+lambda[4]/2.0)-(dQ[2]*sqrt(h_yy*h_zz)*(lambda[2]-lambda[4]))/(h_yy*2.0)+(dQ[0]*h_xz*lapse*(J_c*lambda[0]*shift_x*-2.0+J_c*lambda[2]*shift_x+J_c*lambda[4]*shift_x+lambda[2]*lapse*sqrt(h_yy*h_zz)-lambda[4]*lapse*sqrt(h_yy*h_zz)))/((J_c*J_c)*(shift_x*shift_x)*2.0-h_yy*h_zz*(lapse*lapse)*2.0);
  A_dQ[5] = dQ[5]*(lambda[3]/2.0+lambda[5]/2.0)+(dQ[1]*sqrt(h_yy*h_zz)*(lambda[3]-lambda[5]))/(h_zz*2.0)-(dQ[3]*h_xz*lapse*(h_yy*h_zz*lambda[1]*lapse*-2.0+h_yy*h_zz*lambda[3]*lapse+h_yy*h_zz*lambda[5]*lapse+J_c*lambda[3]*shift_x*sqrt(h_yy*h_zz)-J_c*lambda[5]*shift_x*sqrt(h_yy*h_zz)))/((J_c*J_c)*h_zz*(shift_x*shift_x)*2.0-h_yy*(h_zz*h_zz)*(lapse*lapse)*2.0);
} 
