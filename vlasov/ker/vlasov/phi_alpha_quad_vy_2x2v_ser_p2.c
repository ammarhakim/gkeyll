#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double force_quad = 0.0;
  force_quad = 3.485685011586675*phi[7]+0.7745966692414833*phi[6]-2.598076211353316*phi[5]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[0] -= dx11*force_quad;
  alpha_quad[1] -= dx11*force_quad;
  alpha_quad[2] -= dx11*force_quad;

  force_quad = 0.7745966692414833*phi[6]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[3] -= dx11*force_quad;
  alpha_quad[4] -= dx11*force_quad;
  alpha_quad[5] -= dx11*force_quad;

  force_quad = (-3.485685011586675*phi[7])+0.7745966692414833*phi[6]+2.598076211353316*phi[5]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[6] -= dx11*force_quad;
  alpha_quad[7] -= dx11*force_quad;
  alpha_quad[8] -= dx11*force_quad;

  force_quad = (-0.9682458365518543*phi[6])-2.598076211353316*phi[5]+0.8660254037844386*phi[2];
  alpha_quad[9] -= dx11*force_quad;
  alpha_quad[10] -= dx11*force_quad;
  alpha_quad[11] -= dx11*force_quad;

  force_quad = 0.8660254037844386*phi[2]-0.9682458365518543*phi[6];
  alpha_quad[12] -= dx11*force_quad;
  alpha_quad[13] -= dx11*force_quad;
  alpha_quad[14] -= dx11*force_quad;

  force_quad = (-0.9682458365518543*phi[6])+2.598076211353316*phi[5]+0.8660254037844386*phi[2];
  alpha_quad[15] -= dx11*force_quad;
  alpha_quad[16] -= dx11*force_quad;
  alpha_quad[17] -= dx11*force_quad;

  force_quad = (-3.485685011586675*phi[7])+0.7745966692414833*phi[6]-2.598076211353316*phi[5]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[18] -= dx11*force_quad;
  alpha_quad[19] -= dx11*force_quad;
  alpha_quad[20] -= dx11*force_quad;

  force_quad = 0.7745966692414833*phi[6]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[21] -= dx11*force_quad;
  alpha_quad[22] -= dx11*force_quad;
  alpha_quad[23] -= dx11*force_quad;

  force_quad = 3.485685011586675*phi[7]+0.7745966692414833*phi[6]+2.598076211353316*phi[5]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[24] -= dx11*force_quad;
  alpha_quad[25] -= dx11*force_quad;
  alpha_quad[26] -= dx11*force_quad;

} 
