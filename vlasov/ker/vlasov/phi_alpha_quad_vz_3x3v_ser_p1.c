#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void phi_alpha_quad_vz_3x3v_ser_p1(const double *dxv, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx12 = 2.0/dxv[2]; 

  double force_quad = 0.0;
  force_quad = 0.6123724356957944*phi[7]-0.6123724356957944*phi[6]-0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[0] -= dx12*force_quad;
  alpha_quad[1] -= dx12*force_quad;
  alpha_quad[2] -= dx12*force_quad;
  alpha_quad[3] -= dx12*force_quad;

  force_quad = 0.6123724356957944*phi[7]-0.6123724356957944*phi[6]-0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[4] -= dx12*force_quad;
  alpha_quad[5] -= dx12*force_quad;
  alpha_quad[6] -= dx12*force_quad;
  alpha_quad[7] -= dx12*force_quad;

  force_quad = -(0.6123724356957944*phi[7])+0.6123724356957944*phi[6]-0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[8] -= dx12*force_quad;
  alpha_quad[9] -= dx12*force_quad;
  alpha_quad[10] -= dx12*force_quad;
  alpha_quad[11] -= dx12*force_quad;

  force_quad = -(0.6123724356957944*phi[7])+0.6123724356957944*phi[6]-0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[12] -= dx12*force_quad;
  alpha_quad[13] -= dx12*force_quad;
  alpha_quad[14] -= dx12*force_quad;
  alpha_quad[15] -= dx12*force_quad;

  force_quad = -(0.6123724356957944*phi[7])-0.6123724356957944*phi[6]+0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[16] -= dx12*force_quad;
  alpha_quad[17] -= dx12*force_quad;
  alpha_quad[18] -= dx12*force_quad;
  alpha_quad[19] -= dx12*force_quad;

  force_quad = -(0.6123724356957944*phi[7])-0.6123724356957944*phi[6]+0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[20] -= dx12*force_quad;
  alpha_quad[21] -= dx12*force_quad;
  alpha_quad[22] -= dx12*force_quad;
  alpha_quad[23] -= dx12*force_quad;

  force_quad = 0.6123724356957944*phi[7]+0.6123724356957944*phi[6]+0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[24] -= dx12*force_quad;
  alpha_quad[25] -= dx12*force_quad;
  alpha_quad[26] -= dx12*force_quad;
  alpha_quad[27] -= dx12*force_quad;

  force_quad = 0.6123724356957944*phi[7]+0.6123724356957944*phi[6]+0.6123724356957944*phi[5]+0.6123724356957944*phi[3];
  alpha_quad[28] -= dx12*force_quad;
  alpha_quad[29] -= dx12*force_quad;
  alpha_quad[30] -= dx12*force_quad;
  alpha_quad[31] -= dx12*force_quad;

} 
