#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += -0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0])-0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10; 
  out[2] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*Fhat_nodal[2]-0.3535533905932737*Fhat_nodal[1]+0.3535533905932737*Fhat_nodal[0])*dv10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10; 
  out[5] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])-0.3535533905932737*Fhat_nodal[0])*dv10; 
  out[6] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv10; 

  } else { 

  out[0] += 0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2])-0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10; 
  out[2] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*Fhat_nodal[2]+0.3535533905932737*Fhat_nodal[1]-0.3535533905932737*Fhat_nodal[0])*dv10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10; 
  out[5] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])+0.3535533905932737*Fhat_nodal[0])*dv10; 
  out[6] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv10; 

  } 
  return 0.0;

} 
