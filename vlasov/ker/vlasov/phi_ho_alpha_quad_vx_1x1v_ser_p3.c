#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_ho_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 

  double force_quad = 0.0;
  force_quad = jacob_cx_inv*(8.715656478198236*phi[3]-4.298388424359992*phi[2]+1.224744871391589*phi[1]);
  alpha_quad[0] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(1.2620965531989388*phi[3]-2.554184205050222*phi[2]+1.224744871391589*phi[1]);
  alpha_quad[1] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(1.224744871391589*phi[1]-2.806243040080455*phi[3]);
  alpha_quad[2] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(1.2620965531989388*phi[3]+2.554184205050222*phi[2]+1.224744871391589*phi[1]);
  alpha_quad[3] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(8.715656478198236*phi[3]+4.298388424359992*phi[2]+1.224744871391589*phi[1]);
  alpha_quad[4] -= dx10*force_quad;

} 
