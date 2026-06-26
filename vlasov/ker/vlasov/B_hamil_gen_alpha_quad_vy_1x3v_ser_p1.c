#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_gen_alpha_quad_vy_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[6]; 
  const double *Bz = &qmem[10]; 
  double Bx_quad = 0.0;
  double Bz_quad = 0.0;
  Bx_quad = 0.7071067811865475*Bx[0]-0.7071067811865475*Bx[1];
  Bz_quad = 0.7071067811865475*Bz[0]-0.7071067811865475*Bz[1];
  alpha_quad[0] += dv12*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[12]+0.75*hamil[11]-0.4330127018922193*hamil[9]-0.75*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[1] += dv12*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[12]+0.75*hamil[11]+0.4330127018922193*hamil[9]-0.75*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[2] += dv12*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(-(0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[12]+0.75*hamil[11]-0.4330127018922193*hamil[9]-0.75*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[3] += dv12*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[12]+0.75*hamil[11]+0.4330127018922193*hamil[9]-0.75*hamil[7]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;

  Bx_quad = 0.7071067811865475*Bx[1]+0.7071067811865475*Bx[0];
  Bz_quad = 0.7071067811865475*Bz[1]+0.7071067811865475*Bz[0];
  alpha_quad[4] += dv12*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[12]-0.75*hamil[11]-0.4330127018922193*hamil[9]-0.75*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[5] += dv12*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[12]-0.75*hamil[11]+0.4330127018922193*hamil[9]-0.75*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[6] += dv12*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[12]-0.75*hamil[11]-0.4330127018922193*hamil[9]-0.75*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;
  alpha_quad[7] += dv12*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bx_quad - dv10*(-(0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[12]-0.75*hamil[11]+0.4330127018922193*hamil[9]-0.75*hamil[7]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*Bz_quad;

} 
