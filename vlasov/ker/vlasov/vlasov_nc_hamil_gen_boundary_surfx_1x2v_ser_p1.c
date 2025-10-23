#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += -0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0])-0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]))*dx10; 
  out[3] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*Fhat_nodal[2]-0.3535533905932737*Fhat_nodal[1]+0.3535533905932737*Fhat_nodal[0])*dx10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dx10; 
  out[5] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dx10; 
  out[6] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])-0.3535533905932737*Fhat_nodal[0])*dx10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dx10; 

  } else { 

  out[0] += 0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2])-0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0]))*dx10; 
  out[3] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*Fhat_nodal[2]+0.3535533905932737*Fhat_nodal[1]-0.3535533905932737*Fhat_nodal[0])*dx10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dx10; 
  out[5] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dx10; 
  out[6] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])+0.3535533905932737*Fhat_nodal[0])*dx10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dx10; 

  } 
  return 0.0;

} 
