#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 

  double Fhat_l_nodal[8] = {0.0}; 
  double Fhat_r_nodal[8] = {0.0}; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = (-0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = 0.75*hamil[15]-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.4330127018922193*fl[15])+0.4330127018922193*fl[14]-0.25*fl[13]+0.4330127018922193*fl[12]+0.4330127018922193*fl[11]+0.25*fl[10]-0.4330127018922193*fl[9]+0.25*fl[8]-0.4330127018922193*fl[7]+0.25*fl[6]-0.4330127018922193*fl[5]-0.25*fl[4]-0.25*fl[3]+0.4330127018922193*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]-0.25*fc[13]-0.4330127018922193*fc[12]-0.4330127018922193*fc[11]+0.25*fc[10]+0.4330127018922193*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]+0.25*fc[6]+0.4330127018922193*fc[5]-0.25*fc[4]-0.25*fc[3]-0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.25*fc[13]+0.4330127018922193*fc[12]+0.4330127018922193*fc[11]+0.25*fc[10]-0.4330127018922193*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]+0.25*fc[6]-0.4330127018922193*fc[5]-0.25*fc[4]-0.25*fc[3]+0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = 0.4330127018922193*fr[15]-0.4330127018922193*fr[14]-0.25*fr[13]-0.4330127018922193*fr[12]-0.4330127018922193*fr[11]+0.25*fr[10]+0.4330127018922193*fr[9]+0.25*fr[8]+0.4330127018922193*fr[7]+0.25*fr[6]+0.4330127018922193*fr[5]-0.25*fr[4]-0.25*fr[3]-0.4330127018922193*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-0.75*hamil[15])+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = 0.75*hamil[15]-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.4330127018922193*fl[15]-0.4330127018922193*fl[14]+0.25*fl[13]-0.4330127018922193*fl[12]+0.4330127018922193*fl[11]-0.25*fl[10]+0.4330127018922193*fl[9]-0.25*fl[8]-0.4330127018922193*fl[7]+0.25*fl[6]-0.4330127018922193*fl[5]+0.25*fl[4]-0.25*fl[3]+0.4330127018922193*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]+0.25*fc[13]+0.4330127018922193*fc[12]-0.4330127018922193*fc[11]-0.25*fc[10]-0.4330127018922193*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]+0.25*fc[6]+0.4330127018922193*fc[5]+0.25*fc[4]-0.25*fc[3]-0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.25*fc[13]-0.4330127018922193*fc[12]+0.4330127018922193*fc[11]-0.25*fc[10]+0.4330127018922193*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]+0.25*fc[6]-0.4330127018922193*fc[5]+0.25*fc[4]-0.25*fc[3]+0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = (-0.4330127018922193*fr[15])+0.4330127018922193*fr[14]+0.25*fr[13]+0.4330127018922193*fr[12]-0.4330127018922193*fr[11]-0.25*fr[10]-0.4330127018922193*fr[9]-0.25*fr[8]+0.4330127018922193*fr[7]+0.25*fr[6]+0.4330127018922193*fr[5]+0.25*fr[4]-0.25*fr[3]-0.4330127018922193*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = (-0.75*hamil[15])+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.4330127018922193*fl[15]-0.4330127018922193*fl[14]+0.25*fl[13]+0.4330127018922193*fl[12]-0.4330127018922193*fl[11]-0.25*fl[10]-0.4330127018922193*fl[9]+0.25*fl[8]+0.4330127018922193*fl[7]-0.25*fl[6]-0.4330127018922193*fl[5]-0.25*fl[4]+0.25*fl[3]+0.4330127018922193*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]+0.25*fc[13]-0.4330127018922193*fc[12]+0.4330127018922193*fc[11]-0.25*fc[10]+0.4330127018922193*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]-0.25*fc[6]+0.4330127018922193*fc[5]-0.25*fc[4]+0.25*fc[3]-0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.25*fc[13]+0.4330127018922193*fc[12]-0.4330127018922193*fc[11]-0.25*fc[10]-0.4330127018922193*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]-0.25*fc[6]-0.4330127018922193*fc[5]-0.25*fc[4]+0.25*fc[3]+0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = (-0.4330127018922193*fr[15])+0.4330127018922193*fr[14]+0.25*fr[13]-0.4330127018922193*fr[12]+0.4330127018922193*fr[11]-0.25*fr[10]+0.4330127018922193*fr[9]+0.25*fr[8]-0.4330127018922193*fr[7]-0.25*fr[6]+0.4330127018922193*fr[5]-0.25*fr[4]+0.25*fr[3]-0.4330127018922193*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[2] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[2] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 0.75*hamil[15]-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = (-0.75*hamil[15])+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.4330127018922193*fl[15])+0.4330127018922193*fl[14]-0.25*fl[13]-0.4330127018922193*fl[12]-0.4330127018922193*fl[11]+0.25*fl[10]+0.4330127018922193*fl[9]-0.25*fl[8]+0.4330127018922193*fl[7]-0.25*fl[6]-0.4330127018922193*fl[5]+0.25*fl[4]+0.25*fl[3]+0.4330127018922193*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]-0.25*fc[13]+0.4330127018922193*fc[12]+0.4330127018922193*fc[11]+0.25*fc[10]-0.4330127018922193*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]-0.25*fc[6]+0.4330127018922193*fc[5]+0.25*fc[4]+0.25*fc[3]-0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.25*fc[13]-0.4330127018922193*fc[12]-0.4330127018922193*fc[11]+0.25*fc[10]+0.4330127018922193*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]-0.25*fc[6]-0.4330127018922193*fc[5]+0.25*fc[4]+0.25*fc[3]+0.4330127018922193*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = 0.4330127018922193*fr[15]-0.4330127018922193*fr[14]-0.25*fr[13]+0.4330127018922193*fr[12]+0.4330127018922193*fr[11]+0.25*fr[10]-0.4330127018922193*fr[9]-0.25*fr[8]-0.4330127018922193*fr[7]-0.25*fr[6]+0.4330127018922193*fr[5]+0.25*fr[4]+0.25*fr[3]-0.4330127018922193*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[3] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[3] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = (-0.75*hamil[15])-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.4330127018922193*fl[15]+0.4330127018922193*fl[14]+0.25*fl[13]-0.4330127018922193*fl[12]-0.4330127018922193*fl[11]+0.25*fl[10]-0.4330127018922193*fl[9]-0.25*fl[8]-0.4330127018922193*fl[7]-0.25*fl[6]+0.4330127018922193*fl[5]-0.25*fl[4]-0.25*fl[3]+0.4330127018922193*fl[2]+0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = (-0.4330127018922193*fc[15])-0.4330127018922193*fc[14]+0.25*fc[13]+0.4330127018922193*fc[12]+0.4330127018922193*fc[11]+0.25*fc[10]+0.4330127018922193*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]-0.25*fc[6]-0.4330127018922193*fc[5]-0.25*fc[4]-0.25*fc[3]-0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = 0.4330127018922193*fc[15]+0.4330127018922193*fc[14]+0.25*fc[13]-0.4330127018922193*fc[12]-0.4330127018922193*fc[11]+0.25*fc[10]-0.4330127018922193*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]-0.25*fc[6]+0.4330127018922193*fc[5]-0.25*fc[4]-0.25*fc[3]+0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = (-0.4330127018922193*fr[15])-0.4330127018922193*fr[14]+0.25*fr[13]+0.4330127018922193*fr[12]+0.4330127018922193*fr[11]+0.25*fr[10]+0.4330127018922193*fr[9]-0.25*fr[8]+0.4330127018922193*fr[7]-0.25*fr[6]-0.4330127018922193*fr[5]-0.25*fr[4]-0.25*fr[3]-0.4330127018922193*fr[2]+0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[4] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[4] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 0.75*hamil[15]+0.75*hamil[14]-0.4330127018922193*hamil[13]-0.75*hamil[12]-0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = (-0.75*hamil[15])-0.75*hamil[14]-0.4330127018922193*hamil[13]+0.75*hamil[12]-0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.4330127018922193*fl[15])-0.4330127018922193*fl[14]-0.25*fl[13]+0.4330127018922193*fl[12]-0.4330127018922193*fl[11]-0.25*fl[10]+0.4330127018922193*fl[9]+0.25*fl[8]-0.4330127018922193*fl[7]-0.25*fl[6]+0.4330127018922193*fl[5]+0.25*fl[4]-0.25*fl[3]+0.4330127018922193*fl[2]+0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = 0.4330127018922193*fc[15]+0.4330127018922193*fc[14]-0.25*fc[13]-0.4330127018922193*fc[12]+0.4330127018922193*fc[11]-0.25*fc[10]-0.4330127018922193*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]-0.25*fc[6]-0.4330127018922193*fc[5]+0.25*fc[4]-0.25*fc[3]-0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = (-0.4330127018922193*fc[15])-0.4330127018922193*fc[14]-0.25*fc[13]+0.4330127018922193*fc[12]-0.4330127018922193*fc[11]-0.25*fc[10]+0.4330127018922193*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]-0.25*fc[6]+0.4330127018922193*fc[5]+0.25*fc[4]-0.25*fc[3]+0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = 0.4330127018922193*fr[15]+0.4330127018922193*fr[14]-0.25*fr[13]-0.4330127018922193*fr[12]+0.4330127018922193*fr[11]-0.25*fr[10]-0.4330127018922193*fr[9]+0.25*fr[8]+0.4330127018922193*fr[7]-0.25*fr[6]-0.4330127018922193*fr[5]+0.25*fr[4]-0.25*fr[3]-0.4330127018922193*fr[2]+0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[5] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[5] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = 0.75*hamil[15]+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.4330127018922193*fl[15])-0.4330127018922193*fl[14]-0.25*fl[13]-0.4330127018922193*fl[12]+0.4330127018922193*fl[11]-0.25*fl[10]-0.4330127018922193*fl[9]-0.25*fl[8]+0.4330127018922193*fl[7]+0.25*fl[6]+0.4330127018922193*fl[5]-0.25*fl[4]+0.25*fl[3]+0.4330127018922193*fl[2]+0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = 0.4330127018922193*fc[15]+0.4330127018922193*fc[14]-0.25*fc[13]+0.4330127018922193*fc[12]-0.4330127018922193*fc[11]-0.25*fc[10]+0.4330127018922193*fc[9]-0.25*fc[8]-0.4330127018922193*fc[7]+0.25*fc[6]-0.4330127018922193*fc[5]-0.25*fc[4]+0.25*fc[3]-0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = (-0.4330127018922193*fc[15])-0.4330127018922193*fc[14]-0.25*fc[13]-0.4330127018922193*fc[12]+0.4330127018922193*fc[11]-0.25*fc[10]-0.4330127018922193*fc[9]-0.25*fc[8]+0.4330127018922193*fc[7]+0.25*fc[6]+0.4330127018922193*fc[5]-0.25*fc[4]+0.25*fc[3]+0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = 0.4330127018922193*fr[15]+0.4330127018922193*fr[14]-0.25*fr[13]+0.4330127018922193*fr[12]-0.4330127018922193*fr[11]-0.25*fr[10]+0.4330127018922193*fr[9]-0.25*fr[8]-0.4330127018922193*fr[7]+0.25*fr[6]-0.4330127018922193*fr[5]-0.25*fr[4]+0.25*fr[3]-0.4330127018922193*fr[2]+0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[6] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[6] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-0.75*hamil[15])-0.75*hamil[14]+0.4330127018922193*hamil[13]-0.75*hamil[12]+0.4330127018922193*hamil[10]-0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_r_quad = 0.75*hamil[15]+0.75*hamil[14]+0.4330127018922193*hamil[13]+0.75*hamil[12]+0.4330127018922193*hamil[10]+0.75*hamil[9]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[4]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.4330127018922193*fl[15]+0.4330127018922193*fl[14]+0.25*fl[13]+0.4330127018922193*fl[12]+0.4330127018922193*fl[11]+0.25*fl[10]+0.4330127018922193*fl[9]+0.25*fl[8]+0.4330127018922193*fl[7]+0.25*fl[6]+0.4330127018922193*fl[5]+0.25*fl[4]+0.25*fl[3]+0.4330127018922193*fl[2]+0.25*fl[1]+0.25*fl[0]; 
  f_cl_quad = (-0.4330127018922193*fc[15])-0.4330127018922193*fc[14]+0.25*fc[13]-0.4330127018922193*fc[12]-0.4330127018922193*fc[11]+0.25*fc[10]-0.4330127018922193*fc[9]+0.25*fc[8]-0.4330127018922193*fc[7]+0.25*fc[6]-0.4330127018922193*fc[5]+0.25*fc[4]+0.25*fc[3]-0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_cr_quad = 0.4330127018922193*fc[15]+0.4330127018922193*fc[14]+0.25*fc[13]+0.4330127018922193*fc[12]+0.4330127018922193*fc[11]+0.25*fc[10]+0.4330127018922193*fc[9]+0.25*fc[8]+0.4330127018922193*fc[7]+0.25*fc[6]+0.4330127018922193*fc[5]+0.25*fc[4]+0.25*fc[3]+0.4330127018922193*fc[2]+0.25*fc[1]+0.25*fc[0]; 
  f_rl_quad = (-0.4330127018922193*fr[15])-0.4330127018922193*fr[14]+0.25*fr[13]-0.4330127018922193*fr[12]-0.4330127018922193*fr[11]+0.25*fr[10]-0.4330127018922193*fr[9]+0.25*fr[8]-0.4330127018922193*fr[7]+0.25*fr[6]-0.4330127018922193*fr[5]+0.25*fr[4]+0.25*fr[3]-0.4330127018922193*fr[2]+0.25*fr[1]+0.25*fr[0]; 
  Fhat_l_nodal[7] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[7] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  out[0] += ((-0.25*Fhat_r_nodal[7])+0.25*Fhat_l_nodal[7]-0.25*Fhat_r_nodal[6]+0.25*Fhat_l_nodal[6]-0.25*Fhat_r_nodal[5]+0.25*Fhat_l_nodal[5]-0.25*Fhat_r_nodal[4]+0.25*Fhat_l_nodal[4]-0.25*Fhat_r_nodal[3]+0.25*Fhat_l_nodal[3]-0.25*Fhat_r_nodal[2]+0.25*Fhat_l_nodal[2]-0.25*Fhat_r_nodal[1]+0.25*Fhat_l_nodal[1]-0.25*Fhat_r_nodal[0]+0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[1] += ((-0.25*Fhat_r_nodal[7])+0.25*Fhat_l_nodal[7]-0.25*Fhat_r_nodal[6]+0.25*Fhat_l_nodal[6]-0.25*Fhat_r_nodal[5]+0.25*Fhat_l_nodal[5]-0.25*Fhat_r_nodal[4]+0.25*(Fhat_l_nodal[4]+Fhat_r_nodal[3])-0.25*Fhat_l_nodal[3]+0.25*Fhat_r_nodal[2]-0.25*Fhat_l_nodal[2]+0.25*Fhat_r_nodal[1]-0.25*Fhat_l_nodal[1]+0.25*Fhat_r_nodal[0]-0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[2] += -0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]+Fhat_r_nodal[6]+Fhat_l_nodal[6]+Fhat_r_nodal[5]+Fhat_l_nodal[5]+Fhat_r_nodal[4]+Fhat_l_nodal[4]+Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0])*dv11*dx11; 
  out[3] += ((-0.25*Fhat_r_nodal[7])+0.25*Fhat_l_nodal[7]-0.25*Fhat_r_nodal[6]+0.25*(Fhat_l_nodal[6]+Fhat_r_nodal[5])-0.25*Fhat_l_nodal[5]+0.25*Fhat_r_nodal[4]-0.25*(Fhat_l_nodal[4]+Fhat_r_nodal[3])+0.25*Fhat_l_nodal[3]-0.25*Fhat_r_nodal[2]+0.25*(Fhat_l_nodal[2]+Fhat_r_nodal[1])-0.25*Fhat_l_nodal[1]+0.25*Fhat_r_nodal[0]-0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[4] += ((-0.25*Fhat_r_nodal[7])+0.25*(Fhat_l_nodal[7]+Fhat_r_nodal[6])-0.25*(Fhat_l_nodal[6]+Fhat_r_nodal[5])+0.25*(Fhat_l_nodal[5]+Fhat_r_nodal[4])-0.25*(Fhat_l_nodal[4]+Fhat_r_nodal[3])+0.25*(Fhat_l_nodal[3]+Fhat_r_nodal[2])-0.25*(Fhat_l_nodal[2]+Fhat_r_nodal[1])+0.25*(Fhat_l_nodal[1]+Fhat_r_nodal[0])-0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[5] += (0.4330127018922193*(Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0])-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]+Fhat_r_nodal[6]+Fhat_l_nodal[6]+Fhat_r_nodal[5]+Fhat_l_nodal[5]+Fhat_r_nodal[4]+Fhat_l_nodal[4]))*dv11*dx11; 
  out[6] += ((-0.25*Fhat_r_nodal[7])+0.25*Fhat_l_nodal[7]-0.25*Fhat_r_nodal[6]+0.25*(Fhat_l_nodal[6]+Fhat_r_nodal[5])-0.25*Fhat_l_nodal[5]+0.25*Fhat_r_nodal[4]-0.25*Fhat_l_nodal[4]+0.25*Fhat_r_nodal[3]-0.25*Fhat_l_nodal[3]+0.25*Fhat_r_nodal[2]-0.25*(Fhat_l_nodal[2]+Fhat_r_nodal[1])+0.25*Fhat_l_nodal[1]-0.25*Fhat_r_nodal[0]+0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[7] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]+Fhat_r_nodal[6]+Fhat_l_nodal[6]))+0.4330127018922193*(Fhat_r_nodal[5]+Fhat_l_nodal[5]+Fhat_r_nodal[4]+Fhat_l_nodal[4])-0.4330127018922193*(Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2])+0.4330127018922193*(Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 
  out[8] += ((-0.25*Fhat_r_nodal[7])+0.25*(Fhat_l_nodal[7]+Fhat_r_nodal[6])-0.25*(Fhat_l_nodal[6]+Fhat_r_nodal[5])+0.25*(Fhat_l_nodal[5]+Fhat_r_nodal[4])-0.25*Fhat_l_nodal[4]+0.25*Fhat_r_nodal[3]-0.25*(Fhat_l_nodal[3]+Fhat_r_nodal[2])+0.25*(Fhat_l_nodal[2]+Fhat_r_nodal[1])-0.25*(Fhat_l_nodal[1]+Fhat_r_nodal[0])+0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[9] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]))+0.4330127018922193*(Fhat_r_nodal[6]+Fhat_l_nodal[6])-0.4330127018922193*(Fhat_r_nodal[5]+Fhat_l_nodal[5])+0.4330127018922193*(Fhat_r_nodal[4]+Fhat_l_nodal[4])-0.4330127018922193*(Fhat_r_nodal[3]+Fhat_l_nodal[3])+0.4330127018922193*(Fhat_r_nodal[2]+Fhat_l_nodal[2])-0.4330127018922193*(Fhat_r_nodal[1]+Fhat_l_nodal[1])+0.4330127018922193*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 
  out[10] += ((-0.25*Fhat_r_nodal[7])+0.25*(Fhat_l_nodal[7]+Fhat_r_nodal[6])-0.25*Fhat_l_nodal[6]+0.25*Fhat_r_nodal[5]-0.25*(Fhat_l_nodal[5]+Fhat_r_nodal[4])+0.25*Fhat_l_nodal[4]-0.25*Fhat_r_nodal[3]+0.25*(Fhat_l_nodal[3]+Fhat_r_nodal[2])-0.25*Fhat_l_nodal[2]+0.25*Fhat_r_nodal[1]-0.25*(Fhat_l_nodal[1]+Fhat_r_nodal[0])+0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[11] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]+Fhat_r_nodal[6]+Fhat_l_nodal[6]))+0.4330127018922193*(Fhat_r_nodal[5]+Fhat_l_nodal[5]+Fhat_r_nodal[4]+Fhat_l_nodal[4]+Fhat_r_nodal[3]+Fhat_l_nodal[3]+Fhat_r_nodal[2]+Fhat_l_nodal[2])-0.4330127018922193*(Fhat_r_nodal[1]+Fhat_l_nodal[1]+Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 
  out[12] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]))+0.4330127018922193*(Fhat_r_nodal[6]+Fhat_l_nodal[6])-0.4330127018922193*(Fhat_r_nodal[5]+Fhat_l_nodal[5])+0.4330127018922193*(Fhat_r_nodal[4]+Fhat_l_nodal[4]+Fhat_r_nodal[3]+Fhat_l_nodal[3])-0.4330127018922193*(Fhat_r_nodal[2]+Fhat_l_nodal[2])+0.4330127018922193*(Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.4330127018922193*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 
  out[13] += ((-0.25*Fhat_r_nodal[7])+0.25*(Fhat_l_nodal[7]+Fhat_r_nodal[6])-0.25*Fhat_l_nodal[6]+0.25*Fhat_r_nodal[5]-0.25*(Fhat_l_nodal[5]+Fhat_r_nodal[4])+0.25*(Fhat_l_nodal[4]+Fhat_r_nodal[3])-0.25*(Fhat_l_nodal[3]+Fhat_r_nodal[2])+0.25*Fhat_l_nodal[2]-0.25*Fhat_r_nodal[1]+0.25*(Fhat_l_nodal[1]+Fhat_r_nodal[0])-0.25*Fhat_l_nodal[0])*dv11*dx11; 
  out[14] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]))+0.4330127018922193*(Fhat_r_nodal[6]+Fhat_l_nodal[6]+Fhat_r_nodal[5]+Fhat_l_nodal[5])-0.4330127018922193*(Fhat_r_nodal[4]+Fhat_l_nodal[4]+Fhat_r_nodal[3]+Fhat_l_nodal[3])+0.4330127018922193*(Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.4330127018922193*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 
  out[15] += ((-0.4330127018922193*(Fhat_r_nodal[7]+Fhat_l_nodal[7]))+0.4330127018922193*(Fhat_r_nodal[6]+Fhat_l_nodal[6]+Fhat_r_nodal[5]+Fhat_l_nodal[5])-0.4330127018922193*(Fhat_r_nodal[4]+Fhat_l_nodal[4])+0.4330127018922193*(Fhat_r_nodal[3]+Fhat_l_nodal[3])-0.4330127018922193*(Fhat_r_nodal[2]+Fhat_l_nodal[2]+Fhat_r_nodal[1]+Fhat_l_nodal[1])+0.4330127018922193*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv11*dx11; 

  return fabs(1.5*dv11*dx11*alpha_max);

} 
