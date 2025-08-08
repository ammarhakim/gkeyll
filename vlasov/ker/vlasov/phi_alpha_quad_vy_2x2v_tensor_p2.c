#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  double force_quad = 0.0;
  force_quad = -(3.9547703670805854*phi[8])+4.308049588470176*phi[7]+1.1857789576832338*phi[6]-2.888338995926613*phi[5]-1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[0] -= dx11*force_quad;
  alpha_quad[1] -= dx11*force_quad;
  alpha_quad[2] -= dx11*force_quad;
  alpha_quad[3] -= dx11*force_quad;

  force_quad = -(1.5613636754553142*phi[8])+1.7008401285415227*phi[7]+1.1857789576832338*phi[6]-1.140331086775586*phi[5]-1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[4] -= dx11*force_quad;
  alpha_quad[5] -= dx11*force_quad;
  alpha_quad[6] -= dx11*force_quad;
  alpha_quad[7] -= dx11*force_quad;

  force_quad = 1.5613636754553142*phi[8]-1.7008401285415227*phi[7]+1.1857789576832338*phi[6]+1.140331086775586*phi[5]-1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[8] -= dx11*force_quad;
  alpha_quad[9] -= dx11*force_quad;
  alpha_quad[10] -= dx11*force_quad;
  alpha_quad[11] -= dx11*force_quad;

  force_quad = 3.9547703670805854*phi[8]-4.308049588470176*phi[7]+1.1857789576832338*phi[6]+2.888338995926613*phi[5]-1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[12] -= dx11*force_quad;
  alpha_quad[13] -= dx11*force_quad;
  alpha_quad[14] -= dx11*force_quad;
  alpha_quad[15] -= dx11*force_quad;

  force_quad = 2.1094782708076156*phi[8]+1.7008401285415227*phi[7]-0.6324956225107456*phi[6]-2.888338995926613*phi[5]-0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[16] -= dx11*force_quad;
  alpha_quad[17] -= dx11*force_quad;
  alpha_quad[18] -= dx11*force_quad;
  alpha_quad[19] -= dx11*force_quad;

  force_quad = 0.8328328677734793*phi[8]+0.6715004280822172*phi[7]-0.6324956225107456*phi[6]-1.140331086775586*phi[5]-0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[20] -= dx11*force_quad;
  alpha_quad[21] -= dx11*force_quad;
  alpha_quad[22] -= dx11*force_quad;
  alpha_quad[23] -= dx11*force_quad;

  force_quad = -(0.8328328677734793*phi[8])-0.6715004280822172*phi[7]-0.6324956225107456*phi[6]+1.140331086775586*phi[5]-0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[24] -= dx11*force_quad;
  alpha_quad[25] -= dx11*force_quad;
  alpha_quad[26] -= dx11*force_quad;
  alpha_quad[27] -= dx11*force_quad;

  force_quad = -(2.1094782708076156*phi[8])-1.7008401285415227*phi[7]-0.6324956225107456*phi[6]+2.888338995926613*phi[5]-0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[28] -= dx11*force_quad;
  alpha_quad[29] -= dx11*force_quad;
  alpha_quad[30] -= dx11*force_quad;
  alpha_quad[31] -= dx11*force_quad;

  force_quad = 2.1094782708076156*phi[8]-1.7008401285415227*phi[7]-0.6324956225107456*phi[6]-2.888338995926613*phi[5]+0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[32] -= dx11*force_quad;
  alpha_quad[33] -= dx11*force_quad;
  alpha_quad[34] -= dx11*force_quad;
  alpha_quad[35] -= dx11*force_quad;

  force_quad = 0.8328328677734793*phi[8]-0.6715004280822172*phi[7]-0.6324956225107456*phi[6]-1.140331086775586*phi[5]+0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[36] -= dx11*force_quad;
  alpha_quad[37] -= dx11*force_quad;
  alpha_quad[38] -= dx11*force_quad;
  alpha_quad[39] -= dx11*force_quad;

  force_quad = -(0.8328328677734793*phi[8])+0.6715004280822172*phi[7]-0.6324956225107456*phi[6]+1.140331086775586*phi[5]+0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[40] -= dx11*force_quad;
  alpha_quad[41] -= dx11*force_quad;
  alpha_quad[42] -= dx11*force_quad;
  alpha_quad[43] -= dx11*force_quad;

  force_quad = -(2.1094782708076156*phi[8])+1.7008401285415227*phi[7]-0.6324956225107456*phi[6]+2.888338995926613*phi[5]+0.5099715653772844*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[44] -= dx11*force_quad;
  alpha_quad[45] -= dx11*force_quad;
  alpha_quad[46] -= dx11*force_quad;
  alpha_quad[47] -= dx11*force_quad;

  force_quad = -(3.9547703670805854*phi[8])-4.308049588470176*phi[7]+1.1857789576832338*phi[6]-2.888338995926613*phi[5]+1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[48] -= dx11*force_quad;
  alpha_quad[49] -= dx11*force_quad;
  alpha_quad[50] -= dx11*force_quad;
  alpha_quad[51] -= dx11*force_quad;

  force_quad = -(1.5613636754553142*phi[8])-1.7008401285415227*phi[7]+1.1857789576832338*phi[6]-1.140331086775586*phi[5]+1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[52] -= dx11*force_quad;
  alpha_quad[53] -= dx11*force_quad;
  alpha_quad[54] -= dx11*force_quad;
  alpha_quad[55] -= dx11*force_quad;

  force_quad = 1.5613636754553142*phi[8]+1.7008401285415227*phi[7]+1.1857789576832338*phi[6]+1.140331086775586*phi[5]+1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[56] -= dx11*force_quad;
  alpha_quad[57] -= dx11*force_quad;
  alpha_quad[58] -= dx11*force_quad;
  alpha_quad[59] -= dx11*force_quad;

  force_quad = 3.9547703670805854*phi[8]+4.308049588470176*phi[7]+1.1857789576832338*phi[6]+2.888338995926613*phi[5]+1.291704467391079*phi[3]+0.8660254037844386*phi[2];
  alpha_quad[60] -= dx11*force_quad;
  alpha_quad[61] -= dx11*force_quad;
  alpha_quad[62] -= dx11*force_quad;
  alpha_quad[63] -= dx11*force_quad;

} 
