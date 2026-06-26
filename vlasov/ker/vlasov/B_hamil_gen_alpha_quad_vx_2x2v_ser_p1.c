#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void B_hamil_gen_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *jacob_vel_surf, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double *Bz = &qmem[20]; 
  double Bz_quad = 0.0;
  Bz_quad = 0.5*Bz[3]-0.5*Bz[2]-0.5*Bz[1]+0.5*Bz[0];
  alpha_quad[0] += 2.0/dxv[3]*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;
  alpha_quad[1] += 2.0/dxv[3]*(-(0.75*hamil[15])+0.75*hamil[14]+0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;

  Bz_quad = -(0.5*Bz[3])+0.5*Bz[2]-0.5*Bz[1]+0.5*Bz[0];
  alpha_quad[2] += 2.0/dxv[3]*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;
  alpha_quad[3] += 2.0/dxv[3]*(0.75*hamil[15]-0.75*hamil[14]+0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;

  Bz_quad = -(0.5*Bz[3])-0.5*Bz[2]+0.5*Bz[1]+0.5*Bz[0];
  alpha_quad[4] += 2.0/dxv[3]*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;
  alpha_quad[5] += 2.0/dxv[3]*(0.75*hamil[15]+0.75*hamil[14]-0.75*hamil[13]-0.4330127018922193*hamil[12]-0.75*hamil[10]-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;

  Bz_quad = 0.5*Bz[3]+0.5*Bz[2]+0.5*Bz[1]+0.5*Bz[0];
  alpha_quad[6] += 2.0/dxv[3]*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;
  alpha_quad[7] += 2.0/dxv[3]*(-(0.75*hamil[15])-0.75*hamil[14]-0.75*hamil[13]+0.4330127018922193*hamil[12]-0.75*hamil[10]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4])*Bz_quad;

} 
