#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[2]; 

  const double *Fhat_nodal = &flux[4]; 
  if (edge == -1) { 

  out[0] += -0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[1] += (0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0])-0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]))*dv11; 
  out[2] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*Fhat_nodal[2]-0.3535533905932737*Fhat_nodal[1]+0.3535533905932737*Fhat_nodal[0])*dv11; 
  out[3] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[4] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])-0.3535533905932737*Fhat_nodal[0])*dv11; 
  out[5] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv11; 
  out[6] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv11; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv11; 

  } else { 

  out[0] += 0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[1] += (0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2])-0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11; 
  out[2] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*Fhat_nodal[2]+0.3535533905932737*Fhat_nodal[1]-0.3535533905932737*Fhat_nodal[0])*dv11; 
  out[3] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11; 
  out[4] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])+0.3535533905932737*Fhat_nodal[0])*dv11; 
  out[5] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv11; 
  out[6] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv11; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv11; 

  } 
  return 0.0;

} 
