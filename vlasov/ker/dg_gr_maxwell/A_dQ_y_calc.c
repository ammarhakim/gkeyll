#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void A_dQ_y_calc(const double lapse, const double shift_x, const double shift_y, const double shift_z, 
      const double h_xx, const double h_xy, const double h_xz, const double h_yy, 
      const double h_yz, const double h_zz, const double J_c, const double *lambda, const double *dQ, double* GKYL_RESTRICT A_dQ) 
{ 
  A_dQ[0] = (dQ[0]*((h_xz*h_xz)*lambda[2]+(h_xz*h_xz)*lambda[4]-h_xx*h_zz*lambda[3]-h_xx*h_zz*lambda[5])*(-1.0/2.0))/(h_xx*h_zz-h_xz*h_xz)-(dQ[3]*h_xz*(lambda[2]-lambda[4])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0-(dQ[5]*h_zz*(lambda[3]-lambda[5])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0-(dQ[2]*h_xz*h_zz*(lambda[2]-lambda[3]+lambda[4]-lambda[5]))/(h_xx*h_zz*2.0-(h_xz*h_xz)*2.0)-(J_c*dQ[1]*shift_x*1.0/pow(h_xx*h_zz-h_xz*h_xz,3.0/2.0)*((h_xz*h_xz)*lambda[2]-(h_xz*h_xz)*lambda[4]-h_xx*h_zz*lambda[3]+h_xx*h_zz*lambda[5]))/(lapse*2.0)-(J_c*dQ[4]*h_xz*shift_x*(lambda[0]*-2.0+lambda[2]+lambda[4]))/(lapse*(h_xx*h_zz-h_xz*h_xz)*2.0);
  A_dQ[1] = dQ[1]*lambda[1];
  A_dQ[2] = (dQ[2]*((h_xz*h_xz)*lambda[3]+(h_xz*h_xz)*lambda[5]-h_xx*h_zz*lambda[2]-h_xx*h_zz*lambda[4])*(-1.0/2.0))/(h_xx*h_zz-h_xz*h_xz)+(dQ[3]*h_xx*(lambda[2]-lambda[4])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0+(dQ[5]*h_xz*(lambda[3]-lambda[5])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0+(dQ[0]*h_xx*h_xz*(lambda[2]-lambda[3]+lambda[4]-lambda[5]))/(h_xx*h_zz*2.0-(h_xz*h_xz)*2.0)+(J_c*dQ[4]*h_xx*shift_x*(lambda[0]*-2.0+lambda[2]+lambda[4]))/(lapse*(h_xx*h_zz-h_xz*h_xz)*2.0)+(J_c*dQ[1]*h_xx*h_xz*shift_x*1.0/pow(h_xx*h_zz-h_xz*h_xz,3.0/2.0)*(lambda[2]-lambda[3]-lambda[4]+lambda[5]))/(lapse*2.0);
  A_dQ[3] = dQ[3]*(lambda[2]/2.0+lambda[4]/2.0)+(dQ[0]*h_xz*(lambda[2]-lambda[4])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0+(dQ[2]*h_zz*(lambda[2]-lambda[4])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0+(J_c*dQ[4]*shift_x*(lambda[2]-lambda[4])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/(lapse*2.0)+(J_c*dQ[1]*h_xz*shift_x*(lambda[1]*-2.0+lambda[2]+lambda[4]))/(lapse*(h_xx*h_zz-h_xz*h_xz)*2.0);
  A_dQ[4] = dQ[4]*lambda[0];
  A_dQ[5] = dQ[5]*(lambda[3]/2.0+lambda[5]/2.0)-(dQ[0]*h_xx*(lambda[3]-lambda[5])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0-(dQ[2]*h_xz*(lambda[3]-lambda[5])*1.0/sqrt(h_xx*h_zz-h_xz*h_xz))/2.0-(J_c*dQ[1]*h_xx*shift_x*(lambda[1]*-2.0+lambda[3]+lambda[5]))/(lapse*(h_xx*h_zz-h_xz*h_xz)*2.0);
} 
