#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 

  double Fhat_nodal[8]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 0.75*hamil[15]-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fskin[15])+0.4330127018922193*fskin[14]-0.25*fskin[13]+0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]+0.25*fskin[10]-0.4330127018922193*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]+0.25*fskin[6]-0.4330127018922193*fskin[5]-0.25*fskin[4]-0.25*fskin[3]+0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = 0.4330127018922193*fedge[15]-0.4330127018922193*fedge[14]-0.25*fedge[13]-0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]+0.25*fedge[10]+0.4330127018922193*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]+0.25*fedge[6]+0.4330127018922193*fedge[5]-0.25*fedge[4]-0.25*fedge[3]-0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fskin[15]-0.4330127018922193*fskin[14]+0.25*fskin[13]-0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]-0.25*fskin[10]+0.4330127018922193*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]+0.25*fskin[6]-0.4330127018922193*fskin[5]+0.25*fskin[4]-0.25*fskin[3]+0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = (-0.4330127018922193*fedge[15])+0.4330127018922193*fedge[14]+0.25*fedge[13]+0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]-0.25*fedge[10]-0.4330127018922193*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]+0.25*fedge[6]+0.4330127018922193*fedge[5]+0.25*fedge[4]-0.25*fedge[3]-0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fskin[15]-0.4330127018922193*fskin[14]+0.25*fskin[13]+0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]-0.25*fskin[10]-0.4330127018922193*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]-0.25*fskin[6]-0.4330127018922193*fskin[5]-0.25*fskin[4]+0.25*fskin[3]+0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = (-0.4330127018922193*fedge[15])+0.4330127018922193*fedge[14]+0.25*fedge[13]-0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]-0.25*fedge[10]+0.4330127018922193*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]-0.25*fedge[6]+0.4330127018922193*fedge[5]-0.25*fedge[4]+0.25*fedge[3]-0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fskin[15])+0.4330127018922193*fskin[14]-0.25*fskin[13]-0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]+0.25*fskin[10]+0.4330127018922193*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]-0.25*fskin[6]-0.4330127018922193*fskin[5]+0.25*fskin[4]+0.25*fskin[3]+0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = 0.4330127018922193*fedge[15]-0.4330127018922193*fedge[14]-0.25*fedge[13]+0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]+0.25*fedge[10]-0.4330127018922193*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]-0.25*fedge[6]+0.4330127018922193*fedge[5]+0.25*fedge[4]+0.25*fedge[3]-0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fskin[15]+0.4330127018922193*fskin[14]+0.25*fskin[13]-0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]+0.25*fskin[10]-0.4330127018922193*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]-0.25*fskin[6]+0.4330127018922193*fskin[5]-0.25*fskin[4]-0.25*fskin[3]+0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = (-0.4330127018922193*fedge[15])-0.4330127018922193*fedge[14]+0.25*fedge[13]+0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]+0.25*fedge[10]+0.4330127018922193*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]-0.25*fedge[6]-0.4330127018922193*fedge[5]-0.25*fedge[4]-0.25*fedge[3]-0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[4] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fskin[15])-0.4330127018922193*fskin[14]-0.25*fskin[13]+0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]-0.25*fskin[10]+0.4330127018922193*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]-0.25*fskin[6]+0.4330127018922193*fskin[5]+0.25*fskin[4]-0.25*fskin[3]+0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = 0.4330127018922193*fedge[15]+0.4330127018922193*fedge[14]-0.25*fedge[13]-0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]-0.25*fedge[10]-0.4330127018922193*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]-0.25*fedge[6]-0.4330127018922193*fedge[5]+0.25*fedge[4]-0.25*fedge[3]-0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[5] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fskin[15])-0.4330127018922193*fskin[14]-0.25*fskin[13]-0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]-0.25*fskin[10]-0.4330127018922193*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]+0.25*fskin[6]+0.4330127018922193*fskin[5]-0.25*fskin[4]+0.25*fskin[3]+0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = 0.4330127018922193*fedge[15]+0.4330127018922193*fedge[14]-0.25*fedge[13]+0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]-0.25*fedge[10]+0.4330127018922193*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]+0.25*fedge[6]-0.4330127018922193*fedge[5]-0.25*fedge[4]+0.25*fedge[3]-0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[6] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fskin[15]+0.4330127018922193*fskin[14]+0.25*fskin[13]+0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]+0.25*fskin[10]+0.4330127018922193*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]+0.25*fskin[6]+0.4330127018922193*fskin[5]+0.25*fskin[4]+0.25*fskin[3]+0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  f_r_quad = (-0.4330127018922193*fedge[15])-0.4330127018922193*fedge[14]+0.25*fedge[13]-0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]+0.25*fedge[10]-0.4330127018922193*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]+0.25*fedge[6]-0.4330127018922193*fedge[5]+0.25*fedge[4]+0.25*fedge[3]-0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  Fhat_nodal[7] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += -0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11*dx11; 
  out[1] += (0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11*dx11; 
  out[2] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11*dx11; 
  out[3] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[4] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv11*dx11; 
  out[5] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11*dx11; 
  out[6] += ((-0.25*(Fhat_nodal[7]+Fhat_nodal[6]))+0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[7] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[8] += ((-0.25*Fhat_nodal[7])+0.25*Fhat_nodal[6]-0.25*Fhat_nodal[5]+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv11*dx11; 
  out[9] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[10] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv11*dx11; 
  out[11] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[12] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[13] += ((-0.25*Fhat_nodal[7])+0.25*(Fhat_nodal[6]+Fhat_nodal[5])-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv11*dx11; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (-0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fedge[15])+0.4330127018922193*fedge[14]-0.25*fedge[13]+0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]+0.25*fedge[10]-0.4330127018922193*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]+0.25*fedge[6]-0.4330127018922193*fedge[5]-0.25*fedge[4]-0.25*fedge[3]+0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = 0.4330127018922193*fskin[15]-0.4330127018922193*fskin[14]-0.25*fskin[13]-0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]+0.25*fskin[10]+0.4330127018922193*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]+0.25*fskin[6]+0.4330127018922193*fskin[5]-0.25*fskin[4]-0.25*fskin[3]-0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fedge[15]-0.4330127018922193*fedge[14]+0.25*fedge[13]-0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]-0.25*fedge[10]+0.4330127018922193*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]+0.25*fedge[6]-0.4330127018922193*fedge[5]+0.25*fedge[4]-0.25*fedge[3]+0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = (-0.4330127018922193*fskin[15])+0.4330127018922193*fskin[14]+0.25*fskin[13]+0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]-0.25*fskin[10]-0.4330127018922193*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]+0.25*fskin[6]+0.4330127018922193*fskin[5]+0.25*fskin[4]-0.25*fskin[3]-0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fedge[15]-0.4330127018922193*fedge[14]+0.25*fedge[13]+0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]-0.25*fedge[10]-0.4330127018922193*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]-0.25*fedge[6]-0.4330127018922193*fedge[5]-0.25*fedge[4]+0.25*fedge[3]+0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = (-0.4330127018922193*fskin[15])+0.4330127018922193*fskin[14]+0.25*fskin[13]-0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]-0.25*fskin[10]+0.4330127018922193*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]-0.25*fskin[6]+0.4330127018922193*fskin[5]-0.25*fskin[4]+0.25*fskin[3]-0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fedge[15])+0.4330127018922193*fedge[14]-0.25*fedge[13]-0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]+0.25*fedge[10]+0.4330127018922193*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]-0.25*fedge[6]-0.4330127018922193*fedge[5]+0.25*fedge[4]+0.25*fedge[3]+0.4330127018922193*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = 0.4330127018922193*fskin[15]-0.4330127018922193*fskin[14]-0.25*fskin[13]+0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]+0.25*fskin[10]-0.4330127018922193*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]-0.25*fskin[6]+0.4330127018922193*fskin[5]+0.25*fskin[4]+0.25*fskin[3]-0.4330127018922193*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[3] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fedge[15]+0.4330127018922193*fedge[14]+0.25*fedge[13]-0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]+0.25*fedge[10]-0.4330127018922193*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]-0.25*fedge[6]+0.4330127018922193*fedge[5]-0.25*fedge[4]-0.25*fedge[3]+0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = (-0.4330127018922193*fskin[15])-0.4330127018922193*fskin[14]+0.25*fskin[13]+0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]+0.25*fskin[10]+0.4330127018922193*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]-0.25*fskin[6]-0.4330127018922193*fskin[5]-0.25*fskin[4]-0.25*fskin[3]-0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[4] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fedge[15])-0.4330127018922193*fedge[14]-0.25*fedge[13]+0.4330127018922193*fedge[12]-0.4330127018922193*fedge[11]-0.25*fedge[10]+0.4330127018922193*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]-0.25*fedge[6]+0.4330127018922193*fedge[5]+0.25*fedge[4]-0.25*fedge[3]+0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = 0.4330127018922193*fskin[15]+0.4330127018922193*fskin[14]-0.25*fskin[13]-0.4330127018922193*fskin[12]+0.4330127018922193*fskin[11]-0.25*fskin[10]-0.4330127018922193*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]-0.25*fskin[6]-0.4330127018922193*fskin[5]+0.25*fskin[4]-0.25*fskin[3]-0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[5] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = (-0.4330127018922193*fedge[15])-0.4330127018922193*fedge[14]-0.25*fedge[13]-0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]-0.25*fedge[10]-0.4330127018922193*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]+0.25*fedge[6]+0.4330127018922193*fedge[5]-0.25*fedge[4]+0.25*fedge[3]+0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = 0.4330127018922193*fskin[15]+0.4330127018922193*fskin[14]-0.25*fskin[13]+0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]-0.25*fskin[10]+0.4330127018922193*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]+0.25*fskin[6]-0.4330127018922193*fskin[5]-0.25*fskin[4]+0.25*fskin[3]-0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[6] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  f_l_quad = 0.4330127018922193*fedge[15]+0.4330127018922193*fedge[14]+0.25*fedge[13]+0.4330127018922193*fedge[12]+0.4330127018922193*fedge[11]+0.25*fedge[10]+0.4330127018922193*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]+0.25*fedge[6]+0.4330127018922193*fedge[5]+0.25*fedge[4]+0.25*fedge[3]+0.4330127018922193*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  f_r_quad = (-0.4330127018922193*fskin[15])-0.4330127018922193*fskin[14]+0.25*fskin[13]-0.4330127018922193*fskin[12]-0.4330127018922193*fskin[11]+0.25*fskin[10]-0.4330127018922193*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]+0.25*fskin[6]-0.4330127018922193*fskin[5]+0.25*fskin[4]+0.25*fskin[3]-0.4330127018922193*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  Fhat_nodal[7] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += 0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11*dx11; 
  out[1] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4])-0.25*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[2] += -0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])*dv11*dx11; 
  out[3] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4])+0.25*(Fhat_nodal[3]+Fhat_nodal[2])-0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[4] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*Fhat_nodal[4]+0.25*Fhat_nodal[3]-0.25*Fhat_nodal[2]+0.25*Fhat_nodal[1]-0.25*Fhat_nodal[0])*dv11*dx11; 
  out[5] += (0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2]+Fhat_nodal[1]+Fhat_nodal[0])-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]+Fhat_nodal[5]+Fhat_nodal[4]))*dv11*dx11; 
  out[6] += (0.25*(Fhat_nodal[7]+Fhat_nodal[6])-0.25*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])+0.25*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[7] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4])-0.4330127018922193*(Fhat_nodal[3]+Fhat_nodal[2])+0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[8] += (0.25*Fhat_nodal[7]-0.25*Fhat_nodal[6]+0.25*Fhat_nodal[5]-0.25*(Fhat_nodal[4]+Fhat_nodal[3])+0.25*Fhat_nodal[2]-0.25*Fhat_nodal[1]+0.25*Fhat_nodal[0])*dv11*dx11; 
  out[9] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*Fhat_nodal[4]-0.4330127018922193*Fhat_nodal[3]+0.4330127018922193*Fhat_nodal[2]-0.4330127018922193*Fhat_nodal[1]+0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[10] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*(Fhat_nodal[4]+Fhat_nodal[3])-0.25*(Fhat_nodal[2]+Fhat_nodal[1])+0.25*Fhat_nodal[0])*dv11*dx11; 
  out[11] += ((-0.4330127018922193*(Fhat_nodal[7]+Fhat_nodal[6]))+0.4330127018922193*(Fhat_nodal[5]+Fhat_nodal[4]+Fhat_nodal[3]+Fhat_nodal[2])-0.4330127018922193*(Fhat_nodal[1]+Fhat_nodal[0]))*dv11*dx11; 
  out[12] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*Fhat_nodal[6]-0.4330127018922193*Fhat_nodal[5]+0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])-0.4330127018922193*Fhat_nodal[2]+0.4330127018922193*Fhat_nodal[1]-0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[13] += (0.25*Fhat_nodal[7]-0.25*(Fhat_nodal[6]+Fhat_nodal[5])+0.25*Fhat_nodal[4]-0.25*Fhat_nodal[3]+0.25*(Fhat_nodal[2]+Fhat_nodal[1])-0.25*Fhat_nodal[0])*dv11*dx11; 
  out[14] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*(Fhat_nodal[4]+Fhat_nodal[3])+0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])-0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 
  out[15] += ((-0.4330127018922193*Fhat_nodal[7])+0.4330127018922193*(Fhat_nodal[6]+Fhat_nodal[5])-0.4330127018922193*Fhat_nodal[4]+0.4330127018922193*Fhat_nodal[3]-0.4330127018922193*(Fhat_nodal[2]+Fhat_nodal[1])+0.4330127018922193*Fhat_nodal[0])*dv11*dx11; 

  } 
  return 0.0;

} 
