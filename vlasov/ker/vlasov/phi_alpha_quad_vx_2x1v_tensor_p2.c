#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_2x1v_tensor_p2(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  double force_quad = 0.0;
  force_quad = (-2.323790007724448*phi[8])+0.7745966692414833*phi[7]+3.485685011586675*phi[6]-2.598076211353316*phi[4]-1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[0] -= dx10*force_quad;

  force_quad = 2.904737509655563*phi[8]-0.9682458365518543*phi[7]-2.598076211353316*phi[4]+0.8660254037844386*phi[1];
  alpha_quad[1] -= dx10*force_quad;

  force_quad = (-2.323790007724448*phi[8])+0.7745966692414833*phi[7]-3.485685011586675*phi[6]-2.598076211353316*phi[4]+1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[2] -= dx10*force_quad;

  force_quad = 0.7745966692414833*phi[7]-1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[3] -= dx10*force_quad;

  force_quad = 0.8660254037844386*phi[1]-0.9682458365518543*phi[7];
  alpha_quad[4] -= dx10*force_quad;

  force_quad = 0.7745966692414833*phi[7]+1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[5] -= dx10*force_quad;

  force_quad = 2.323790007724448*phi[8]+0.7745966692414833*phi[7]-3.485685011586675*phi[6]+2.598076211353316*phi[4]-1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[6] -= dx10*force_quad;

  force_quad = (-2.904737509655563*phi[8])-0.9682458365518543*phi[7]+2.598076211353316*phi[4]+0.8660254037844386*phi[1];
  alpha_quad[7] -= dx10*force_quad;

  force_quad = 2.323790007724448*phi[8]+0.7745966692414833*phi[7]+3.485685011586675*phi[6]+2.598076211353316*phi[4]+1.161895003862225*phi[3]+0.8660254037844386*phi[1];
  alpha_quad[8] -= dx10*force_quad;

} 
