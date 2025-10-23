#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x2v_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  out[0] += ((-0.3535533905932737*Fhat_r_nodal_c[3])+0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*Fhat_l_nodal_c[2]-0.3535533905932737*Fhat_r_nodal_c[1]+0.3535533905932737*Fhat_l_nodal_c[1]-0.3535533905932737*Fhat_r_nodal_c[0]+0.3535533905932737*Fhat_l_nodal_c[0])*dx10; 
  out[1] += -0.6123724356957944*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0])*dx10; 
  out[2] += ((-0.3535533905932737*Fhat_r_nodal_c[3])+0.3535533905932737*Fhat_l_nodal_c[3]-0.3535533905932737*Fhat_r_nodal_c[2]+0.3535533905932737*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1])-0.3535533905932737*Fhat_l_nodal_c[1]+0.3535533905932737*Fhat_r_nodal_c[0]-0.3535533905932737*Fhat_l_nodal_c[0])*dx10; 
  out[3] += ((-0.3535533905932737*Fhat_r_nodal_c[3])+0.3535533905932737*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2])-0.3535533905932737*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1])+0.3535533905932737*(Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0])-0.3535533905932737*Fhat_l_nodal_c[0])*dx10; 
  out[4] += (0.6123724356957944*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0])-0.6123724356957944*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]))*dx10; 
  out[5] += ((-0.6123724356957944*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]))+0.6123724356957944*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])-0.6123724356957944*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.6123724356957944*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dx10; 
  out[6] += ((-0.3535533905932737*Fhat_r_nodal_c[3])+0.3535533905932737*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2])-0.3535533905932737*Fhat_l_nodal_c[2]+0.3535533905932737*Fhat_r_nodal_c[1]-0.3535533905932737*(Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0])+0.3535533905932737*Fhat_l_nodal_c[0])*dx10; 
  out[7] += ((-0.6123724356957944*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]))+0.6123724356957944*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.6123724356957944*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dx10; 

  return 0.0;

} 
