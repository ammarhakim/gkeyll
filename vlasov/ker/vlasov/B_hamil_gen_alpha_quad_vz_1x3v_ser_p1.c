#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_gen_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  double Bx_quad = 0.0;
  double By_quad = 0.0;
  Bx_quad = 0.7071067811865475*Bx[0]-0.7071067811865475*Bx[1];
  By_quad = 0.7071067811865475*By[0]-0.7071067811865475*By[1];
  alpha_quad[0] += dv10*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[11]-0.75*hamil[10]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[1] += dv10*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[11]-0.75*hamil[10]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[2] += dv10*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[11]-0.75*hamil[10]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[3] += dv10*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[11]-0.75*hamil[10]+0.4330127018922193*hamil[7]-0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;

  Bx_quad = 0.7071067811865475*Bx[1]+0.7071067811865475*Bx[0];
  By_quad = 0.7071067811865475*By[1]+0.7071067811865475*By[0];
  alpha_quad[4] += dv10*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[11]-0.75*hamil[10]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[5] += dv10*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[11]-0.75*hamil[10]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[6] += dv10*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[12]-0.4330127018922193*hamil[11]-0.75*hamil[9]-0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[11]-0.75*hamil[10]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;
  alpha_quad[7] += dv10*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[12]+0.4330127018922193*hamil[11]-0.75*hamil[9]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*By_quad - dv11*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[11]-0.75*hamil[10]+0.4330127018922193*hamil[7]+0.4330127018922193*hamil[6]+0.4330127018922193*hamil[3])*Bx_quad;

} 
