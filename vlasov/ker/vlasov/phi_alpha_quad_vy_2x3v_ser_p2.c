#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vy_2x3v_ser_p2(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double force_quad = 0.0;
  force_quad = 3.485685011586675*phi[7]+0.7745966692414833*phi[6]-2.598076211353316*phi[5]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[0] -= dx11*force_quad;
  alpha_quad[1] -= dx11*force_quad;
  alpha_quad[2] -= dx11*force_quad;
  alpha_quad[3] -= dx11*force_quad;
  alpha_quad[4] -= dx11*force_quad;
  alpha_quad[5] -= dx11*force_quad;
  alpha_quad[6] -= dx11*force_quad;
  alpha_quad[7] -= dx11*force_quad;
  alpha_quad[8] -= dx11*force_quad;

  force_quad = 0.7745966692414833*phi[6]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[9] -= dx11*force_quad;
  alpha_quad[10] -= dx11*force_quad;
  alpha_quad[11] -= dx11*force_quad;
  alpha_quad[12] -= dx11*force_quad;
  alpha_quad[13] -= dx11*force_quad;
  alpha_quad[14] -= dx11*force_quad;
  alpha_quad[15] -= dx11*force_quad;
  alpha_quad[16] -= dx11*force_quad;
  alpha_quad[17] -= dx11*force_quad;

  force_quad = (-3.485685011586675*phi[7])+0.7745966692414833*phi[6]+2.598076211353316*phi[5]-1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[18] -= dx11*force_quad;
  alpha_quad[19] -= dx11*force_quad;
  alpha_quad[20] -= dx11*force_quad;
  alpha_quad[21] -= dx11*force_quad;
  alpha_quad[22] -= dx11*force_quad;
  alpha_quad[23] -= dx11*force_quad;
  alpha_quad[24] -= dx11*force_quad;
  alpha_quad[25] -= dx11*force_quad;
  alpha_quad[26] -= dx11*force_quad;

  force_quad = (-0.9682458365518543*phi[6])-2.598076211353316*phi[5]+0.8660254037844386*phi[2];
  alpha_quad[27] -= dx11*force_quad;
  alpha_quad[28] -= dx11*force_quad;
  alpha_quad[29] -= dx11*force_quad;
  alpha_quad[30] -= dx11*force_quad;
  alpha_quad[31] -= dx11*force_quad;
  alpha_quad[32] -= dx11*force_quad;
  alpha_quad[33] -= dx11*force_quad;
  alpha_quad[34] -= dx11*force_quad;
  alpha_quad[35] -= dx11*force_quad;

  force_quad = 0.8660254037844386*phi[2]-0.9682458365518543*phi[6];
  alpha_quad[36] -= dx11*force_quad;
  alpha_quad[37] -= dx11*force_quad;
  alpha_quad[38] -= dx11*force_quad;
  alpha_quad[39] -= dx11*force_quad;
  alpha_quad[40] -= dx11*force_quad;
  alpha_quad[41] -= dx11*force_quad;
  alpha_quad[42] -= dx11*force_quad;
  alpha_quad[43] -= dx11*force_quad;
  alpha_quad[44] -= dx11*force_quad;

  force_quad = (-0.9682458365518543*phi[6])+2.598076211353316*phi[5]+0.8660254037844386*phi[2];
  alpha_quad[45] -= dx11*force_quad;
  alpha_quad[46] -= dx11*force_quad;
  alpha_quad[47] -= dx11*force_quad;
  alpha_quad[48] -= dx11*force_quad;
  alpha_quad[49] -= dx11*force_quad;
  alpha_quad[50] -= dx11*force_quad;
  alpha_quad[51] -= dx11*force_quad;
  alpha_quad[52] -= dx11*force_quad;
  alpha_quad[53] -= dx11*force_quad;

  force_quad = (-3.485685011586675*phi[7])+0.7745966692414833*phi[6]-2.598076211353316*phi[5]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[54] -= dx11*force_quad;
  alpha_quad[55] -= dx11*force_quad;
  alpha_quad[56] -= dx11*force_quad;
  alpha_quad[57] -= dx11*force_quad;
  alpha_quad[58] -= dx11*force_quad;
  alpha_quad[59] -= dx11*force_quad;
  alpha_quad[60] -= dx11*force_quad;
  alpha_quad[61] -= dx11*force_quad;
  alpha_quad[62] -= dx11*force_quad;

  force_quad = 0.7745966692414833*phi[6]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[63] -= dx11*force_quad;
  alpha_quad[64] -= dx11*force_quad;
  alpha_quad[65] -= dx11*force_quad;
  alpha_quad[66] -= dx11*force_quad;
  alpha_quad[67] -= dx11*force_quad;
  alpha_quad[68] -= dx11*force_quad;
  alpha_quad[69] -= dx11*force_quad;
  alpha_quad[70] -= dx11*force_quad;
  alpha_quad[71] -= dx11*force_quad;

  force_quad = 3.485685011586675*phi[7]+0.7745966692414833*phi[6]+2.598076211353316*phi[5]+1.161895003862225*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[72] -= dx11*force_quad;
  alpha_quad[73] -= dx11*force_quad;
  alpha_quad[74] -= dx11*force_quad;
  alpha_quad[75] -= dx11*force_quad;
  alpha_quad[76] -= dx11*force_quad;
  alpha_quad[77] -= dx11*force_quad;
  alpha_quad[78] -= dx11*force_quad;
  alpha_quad[79] -= dx11*force_quad;
  alpha_quad[80] -= dx11*force_quad;

} 
