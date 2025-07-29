#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  double Ghat_r[2]; 
  double Ghat_l[2]; 
  double Ghat_l_nodal[2]; 
  double Ghat_r_nodal[2]; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = 0.8660254037844386*hamil[2]-1.5*hamil[3]; 
  alpha_r_quad = 1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = -(0.8660254037844386*fl[3])-0.5*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.8660254037844386*fc[3]-0.5*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = -(0.8660254037844386*fl[3])-0.5*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_rl_quad = 0.8660254037844386*fr[3]-0.5*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 0.8660254037844386*hamil[2]-1.5*hamil[3]; 
  alpha_r_quad = 1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.8660254037844386*fl[3]+0.5*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = -(0.8660254037844386*fc[3])+0.5*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 0.8660254037844386*fl[3]+0.5*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_rl_quad = -(0.8660254037844386*fr[3])+0.5*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Ghat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Ghat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  Ghat_l[0] = 0.7071067811865475*Ghat_l_nodal[1]+0.7071067811865475*Ghat_l_nodal[0]; 
  Ghat_l[1] = 0.7071067811865475*Ghat_l_nodal[1]-0.7071067811865475*Ghat_l_nodal[0]; 

  Ghat_r[0] = 0.7071067811865475*Ghat_r_nodal[1]+0.7071067811865475*Ghat_r_nodal[0]; 
  Ghat_r[1] = 0.7071067811865475*Ghat_r_nodal[1]-0.7071067811865475*Ghat_r_nodal[0]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10*dx10; 
  out[1] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10*dx10); 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10*dx10; 
  out[3] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10*dx10); 

  return 1.5*dv10*dx10*alpha_max;

} 
