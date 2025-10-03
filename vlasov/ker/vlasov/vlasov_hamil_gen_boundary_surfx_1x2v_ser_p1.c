#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Fhat_nodal[4]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = 0.6123724356957944*fskin[7]+0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = (-0.6123724356957944*fedge[7])+0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = (-0.6123724356957944*fskin[7])-0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = 0.6123724356957944*fedge[7]-0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-1.060660171779821*hamil[7])-0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = (-0.6123724356957944*fskin[7])-0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = 0.6123724356957944*fedge[7]-0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.060660171779821*hamil[7]+0.6123724356957944*hamil[6]+1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = 0.6123724356957944*fskin[7]+0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]+0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  f_r_quad = (-0.6123724356957944*fedge[7])+0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]-0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  Fhat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += -0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10*dx10; 
  out[1] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10*dx10; 
  out[2] += (0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0])-0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10*dx10; 
  out[3] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*Fhat_nodal[2]-0.3535533905932737*Fhat_nodal[1]+0.3535533905932737*Fhat_nodal[0])*dv10*dx10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10*dx10; 
  out[5] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv10*dx10; 
  out[6] += ((-0.3535533905932737*Fhat_nodal[3])+0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])-0.3535533905932737*Fhat_nodal[0])*dv10*dx10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv10*dx10; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = 0.6123724356957944*fedge[7]+0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = (-0.6123724356957944*fskin[7])+0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = (-0.6123724356957944*fedge[7])-0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]-0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]-0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = 0.6123724356957944*fskin[7]-0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]+0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]-0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.060660171779821*hamil[7]-0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = (-0.6123724356957944*fedge[7])-0.3535533905932737*fedge[6]-0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]-0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = 0.6123724356957944*fskin[7]-0.3535533905932737*fskin[6]+0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]-0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-1.060660171779821*hamil[7])+0.6123724356957944*hamil[6]-1.060660171779821*hamil[4]+0.6123724356957944*hamil[2]; 
  f_l_quad = 0.6123724356957944*fedge[7]+0.3535533905932737*fedge[6]+0.6123724356957944*fedge[5]+0.6123724356957944*fedge[4]+0.3535533905932737*fedge[3]+0.3535533905932737*fedge[2]+0.6123724356957944*fedge[1]+0.3535533905932737*fedge[0]; 
  f_r_quad = (-0.6123724356957944*fskin[7])+0.3535533905932737*fskin[6]-0.6123724356957944*fskin[5]-0.6123724356957944*fskin[4]+0.3535533905932737*fskin[3]+0.3535533905932737*fskin[2]-0.6123724356957944*fskin[1]+0.3535533905932737*fskin[0]; 
  Fhat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += 0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10*dx10; 
  out[1] += -0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv10*dx10; 
  out[2] += (0.3535533905932737*(Fhat_nodal[3]+Fhat_nodal[2])-0.3535533905932737*(Fhat_nodal[1]+Fhat_nodal[0]))*dv10*dx10; 
  out[3] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*Fhat_nodal[2]+0.3535533905932737*Fhat_nodal[1]-0.3535533905932737*Fhat_nodal[0])*dv10*dx10; 
  out[4] += (0.6123724356957944*(Fhat_nodal[1]+Fhat_nodal[0])-0.6123724356957944*(Fhat_nodal[3]+Fhat_nodal[2]))*dv10*dx10; 
  out[5] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*Fhat_nodal[2]-0.6123724356957944*Fhat_nodal[1]+0.6123724356957944*Fhat_nodal[0])*dv10*dx10; 
  out[6] += (0.3535533905932737*Fhat_nodal[3]-0.3535533905932737*(Fhat_nodal[2]+Fhat_nodal[1])+0.3535533905932737*Fhat_nodal[0])*dv10*dx10; 
  out[7] += ((-0.6123724356957944*Fhat_nodal[3])+0.6123724356957944*(Fhat_nodal[2]+Fhat_nodal[1])-0.6123724356957944*Fhat_nodal[0])*dv10*dx10; 

  } 
  return 0.0;

} 
