#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[3]; 

  const double *Fhat_nodal = &flux[8]; 
  if (edge == -1) { 

  out[0] += -0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[1] += (0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11; 
  out[2] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[3] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv11; 
  out[4] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[5] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[6] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv11; 
  out[7] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv11; 
  out[8] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11; 
  out[9] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[10] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[11] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv11; 
  out[12] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[13] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv11; 

  } else { 

  out[0] += 0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[1] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[2] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4])+0.25*(Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[3] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv11; 
  out[4] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[5] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[6] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv11; 
  out[7] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv11; 
  out[8] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11; 
  out[9] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[10] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[11] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv11; 
  out[12] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[13] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv11; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv11; 

  } 
  return 0.0;

} 
