#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += -0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv10; 
  out[2] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[4] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv10; 
  out[5] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv10; 
  out[6] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[7] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[8] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv10; 
  out[9] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[10] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv10; 
  out[11] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[12] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[13] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv10; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv10; 

  } else { 

  out[0] += 0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[2] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4])+0.25*(Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[4] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv10; 
  out[5] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv10; 
  out[6] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[7] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[8] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv10; 
  out[9] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[10] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv10; 
  out[11] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[12] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[13] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv10; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv10; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv10; 

  } 
  return 0.0;

} 
