#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 

  double force_quad = 0.0;
  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]+3.4856850115866753*phi[6]-2.5980762113533156*phi[4]-1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[0] -= dx10*force_quad;
  alpha_quad[1] -= dx10*force_quad;
  alpha_quad[2] -= dx10*force_quad;
  alpha_quad[3] -= dx10*force_quad;
  alpha_quad[4] -= dx10*force_quad;
  alpha_quad[5] -= dx10*force_quad;
  alpha_quad[6] -= dx10*force_quad;
  alpha_quad[7] -= dx10*force_quad;
  alpha_quad[8] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(-(0.9682458365518543*phi[7])-2.5980762113533156*phi[4]+0.8660254037844386*phi[1]);
  alpha_quad[9] -= dx10*force_quad;
  alpha_quad[10] -= dx10*force_quad;
  alpha_quad[11] -= dx10*force_quad;
  alpha_quad[12] -= dx10*force_quad;
  alpha_quad[13] -= dx10*force_quad;
  alpha_quad[14] -= dx10*force_quad;
  alpha_quad[15] -= dx10*force_quad;
  alpha_quad[16] -= dx10*force_quad;
  alpha_quad[17] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]-3.4856850115866753*phi[6]-2.5980762113533156*phi[4]+1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[18] -= dx10*force_quad;
  alpha_quad[19] -= dx10*force_quad;
  alpha_quad[20] -= dx10*force_quad;
  alpha_quad[21] -= dx10*force_quad;
  alpha_quad[22] -= dx10*force_quad;
  alpha_quad[23] -= dx10*force_quad;
  alpha_quad[24] -= dx10*force_quad;
  alpha_quad[25] -= dx10*force_quad;
  alpha_quad[26] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]-1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[27] -= dx10*force_quad;
  alpha_quad[28] -= dx10*force_quad;
  alpha_quad[29] -= dx10*force_quad;
  alpha_quad[30] -= dx10*force_quad;
  alpha_quad[31] -= dx10*force_quad;
  alpha_quad[32] -= dx10*force_quad;
  alpha_quad[33] -= dx10*force_quad;
  alpha_quad[34] -= dx10*force_quad;
  alpha_quad[35] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.8660254037844386*phi[1]-0.9682458365518543*phi[7]);
  alpha_quad[36] -= dx10*force_quad;
  alpha_quad[37] -= dx10*force_quad;
  alpha_quad[38] -= dx10*force_quad;
  alpha_quad[39] -= dx10*force_quad;
  alpha_quad[40] -= dx10*force_quad;
  alpha_quad[41] -= dx10*force_quad;
  alpha_quad[42] -= dx10*force_quad;
  alpha_quad[43] -= dx10*force_quad;
  alpha_quad[44] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]+1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[45] -= dx10*force_quad;
  alpha_quad[46] -= dx10*force_quad;
  alpha_quad[47] -= dx10*force_quad;
  alpha_quad[48] -= dx10*force_quad;
  alpha_quad[49] -= dx10*force_quad;
  alpha_quad[50] -= dx10*force_quad;
  alpha_quad[51] -= dx10*force_quad;
  alpha_quad[52] -= dx10*force_quad;
  alpha_quad[53] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]-3.4856850115866753*phi[6]+2.5980762113533156*phi[4]-1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[54] -= dx10*force_quad;
  alpha_quad[55] -= dx10*force_quad;
  alpha_quad[56] -= dx10*force_quad;
  alpha_quad[57] -= dx10*force_quad;
  alpha_quad[58] -= dx10*force_quad;
  alpha_quad[59] -= dx10*force_quad;
  alpha_quad[60] -= dx10*force_quad;
  alpha_quad[61] -= dx10*force_quad;
  alpha_quad[62] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(-(0.9682458365518543*phi[7])+2.5980762113533156*phi[4]+0.8660254037844386*phi[1]);
  alpha_quad[63] -= dx10*force_quad;
  alpha_quad[64] -= dx10*force_quad;
  alpha_quad[65] -= dx10*force_quad;
  alpha_quad[66] -= dx10*force_quad;
  alpha_quad[67] -= dx10*force_quad;
  alpha_quad[68] -= dx10*force_quad;
  alpha_quad[69] -= dx10*force_quad;
  alpha_quad[70] -= dx10*force_quad;
  alpha_quad[71] -= dx10*force_quad;

  force_quad = jacob_cx_inv*(0.7745966692414833*phi[7]+3.4856850115866753*phi[6]+2.5980762113533156*phi[4]+1.1618950038622249*phi[3]+0.8660254037844386*phi[1]);
  alpha_quad[72] -= dx10*force_quad;
  alpha_quad[73] -= dx10*force_quad;
  alpha_quad[74] -= dx10*force_quad;
  alpha_quad[75] -= dx10*force_quad;
  alpha_quad[76] -= dx10*force_quad;
  alpha_quad[77] -= dx10*force_quad;
  alpha_quad[78] -= dx10*force_quad;
  alpha_quad[79] -= dx10*force_quad;
  alpha_quad[80] -= dx10*force_quad;

} 
