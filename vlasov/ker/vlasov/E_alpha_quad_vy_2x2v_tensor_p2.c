#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void E_alpha_quad_vy_2x2v_tensor_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double *Ey = &qmem[9]; 

  double force_quad = 0.0;
  force_quad = 0.9373811576562242*Ey[8]-1.02111731798518*Ey[7]-1.02111731798518*Ey[6]+0.6846098004178088*Ey[5]+0.6846098004178088*Ey[4]+1.1123336207187138*Ey[3]-0.7457659219616816*Ey[2]-0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[0] += force_quad;
  alpha_quad[1] += force_quad;
  alpha_quad[2] += force_quad;
  alpha_quad[3] += force_quad;

  force_quad = -(0.5*Ey[8])+0.5446649474682878*Ey[7]-0.4031423674941089*Ey[6]-0.3651715179178389*Ey[5]+0.6846098004178088*Ey[4]+0.43915503282683993*Ey[3]-0.29443222054963*Ey[2]-0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[4] += force_quad;
  alpha_quad[5] += force_quad;
  alpha_quad[6] += force_quad;
  alpha_quad[7] += force_quad;

  force_quad = -(0.5*Ey[8])+0.5446649474682878*Ey[7]+0.4031423674941089*Ey[6]-0.3651715179178389*Ey[5]+0.6846098004178088*Ey[4]-0.43915503282683993*Ey[3]+0.29443222054963*Ey[2]-0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[8] += force_quad;
  alpha_quad[9] += force_quad;
  alpha_quad[10] += force_quad;
  alpha_quad[11] += force_quad;

  force_quad = 0.9373811576562242*Ey[8]-1.02111731798518*Ey[7]+1.02111731798518*Ey[6]+0.6846098004178088*Ey[5]+0.6846098004178088*Ey[4]-1.1123336207187138*Ey[3]+0.7457659219616816*Ey[2]-0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[12] += force_quad;
  alpha_quad[13] += force_quad;
  alpha_quad[14] += force_quad;
  alpha_quad[15] += force_quad;

  force_quad = -(0.5*Ey[8])-0.40314236749410887*Ey[7]+0.5446649474682878*Ey[6]+0.6846098004178088*Ey[5]-0.3651715179178389*Ey[4]+0.43915503282683993*Ey[3]-0.7457659219616816*Ey[2]-0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[16] += force_quad;
  alpha_quad[17] += force_quad;
  alpha_quad[18] += force_quad;
  alpha_quad[19] += force_quad;

  force_quad = 0.26670047499683713*Ey[8]+0.2150365218040566*Ey[7]+0.21503652180405658*Ey[6]-0.3651715179178389*Ey[5]-0.3651715179178389*Ey[4]+0.1733806649955719*Ey[3]-0.29443222054963*Ey[2]-0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[20] += force_quad;
  alpha_quad[21] += force_quad;
  alpha_quad[22] += force_quad;
  alpha_quad[23] += force_quad;

  force_quad = 0.26670047499683713*Ey[8]+0.2150365218040566*Ey[7]-0.21503652180405658*Ey[6]-0.3651715179178389*Ey[5]-0.3651715179178389*Ey[4]-0.1733806649955719*Ey[3]+0.29443222054963*Ey[2]-0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[24] += force_quad;
  alpha_quad[25] += force_quad;
  alpha_quad[26] += force_quad;
  alpha_quad[27] += force_quad;

  force_quad = -(0.5*Ey[8])-0.40314236749410887*Ey[7]-0.5446649474682878*Ey[6]+0.6846098004178088*Ey[5]-0.3651715179178389*Ey[4]-0.43915503282683993*Ey[3]+0.7457659219616816*Ey[2]-0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[28] += force_quad;
  alpha_quad[29] += force_quad;
  alpha_quad[30] += force_quad;
  alpha_quad[31] += force_quad;

  force_quad = -(0.5*Ey[8])+0.40314236749410887*Ey[7]+0.5446649474682878*Ey[6]+0.6846098004178088*Ey[5]-0.3651715179178389*Ey[4]-0.43915503282683993*Ey[3]-0.7457659219616816*Ey[2]+0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[32] += force_quad;
  alpha_quad[33] += force_quad;
  alpha_quad[34] += force_quad;
  alpha_quad[35] += force_quad;

  force_quad = 0.26670047499683713*Ey[8]-0.2150365218040566*Ey[7]+0.21503652180405658*Ey[6]-0.3651715179178389*Ey[5]-0.3651715179178389*Ey[4]-0.1733806649955719*Ey[3]-0.29443222054963*Ey[2]+0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[36] += force_quad;
  alpha_quad[37] += force_quad;
  alpha_quad[38] += force_quad;
  alpha_quad[39] += force_quad;

  force_quad = 0.26670047499683713*Ey[8]-0.2150365218040566*Ey[7]-0.21503652180405658*Ey[6]-0.3651715179178389*Ey[5]-0.3651715179178389*Ey[4]+0.1733806649955719*Ey[3]+0.29443222054963*Ey[2]+0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[40] += force_quad;
  alpha_quad[41] += force_quad;
  alpha_quad[42] += force_quad;
  alpha_quad[43] += force_quad;

  force_quad = -(0.5*Ey[8])+0.40314236749410887*Ey[7]-0.5446649474682878*Ey[6]+0.6846098004178088*Ey[5]-0.3651715179178389*Ey[4]+0.43915503282683993*Ey[3]+0.7457659219616816*Ey[2]+0.29443222054963*Ey[1]+0.5*Ey[0];
  alpha_quad[44] += force_quad;
  alpha_quad[45] += force_quad;
  alpha_quad[46] += force_quad;
  alpha_quad[47] += force_quad;

  force_quad = 0.9373811576562242*Ey[8]+1.02111731798518*Ey[7]-1.02111731798518*Ey[6]+0.6846098004178088*Ey[5]+0.6846098004178088*Ey[4]-1.1123336207187138*Ey[3]-0.7457659219616816*Ey[2]+0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[48] += force_quad;
  alpha_quad[49] += force_quad;
  alpha_quad[50] += force_quad;
  alpha_quad[51] += force_quad;

  force_quad = -(0.5*Ey[8])-0.5446649474682878*Ey[7]-0.4031423674941089*Ey[6]-0.3651715179178389*Ey[5]+0.6846098004178088*Ey[4]-0.43915503282683993*Ey[3]-0.29443222054963*Ey[2]+0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[52] += force_quad;
  alpha_quad[53] += force_quad;
  alpha_quad[54] += force_quad;
  alpha_quad[55] += force_quad;

  force_quad = -(0.5*Ey[8])-0.5446649474682878*Ey[7]+0.4031423674941089*Ey[6]-0.3651715179178389*Ey[5]+0.6846098004178088*Ey[4]+0.43915503282683993*Ey[3]+0.29443222054963*Ey[2]+0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[56] += force_quad;
  alpha_quad[57] += force_quad;
  alpha_quad[58] += force_quad;
  alpha_quad[59] += force_quad;

  force_quad = 0.9373811576562242*Ey[8]+1.02111731798518*Ey[7]+1.02111731798518*Ey[6]+0.6846098004178088*Ey[5]+0.6846098004178088*Ey[4]+1.1123336207187138*Ey[3]+0.7457659219616816*Ey[2]+0.7457659219616816*Ey[1]+0.5*Ey[0];
  alpha_quad[60] += force_quad;
  alpha_quad[61] += force_quad;
  alpha_quad[62] += force_quad;
  alpha_quad[63] += force_quad;

} 
