#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Fhat_l_nodal[3] = {0.0}; 
  double Fhat_r_nodal[3] = {0.0}; 
  double alpha_l_quad = 0.0; 
  double alpha_r_quad = 0.0; 
  double alpha_max = 0.0; 
  double f_lr_quad = 0.0; 
  double f_cl_quad = 0.0; 
  double f_cr_quad = 0.0; 
  double f_rl_quad = 0.0; 
  alpha_l_quad = 4.499999999999999*hamil[7]+1.936491673103709*hamil[6]-2.598076211353316*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = (-4.499999999999999*hamil[7])+1.936491673103709*hamil[6]-2.598076211353316*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.7745966692414833*fl[7]-1.5*fl[6]+0.4472135954999579*fl[5]+1.118033988749895*fl[4]-1.161895003862225*fl[3]-0.6708203932499369*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-0.7745966692414833*fc[7])-1.5*fc[6]+0.4472135954999579*fc[5]+1.118033988749895*fc[4]+1.161895003862225*fc[3]-0.6708203932499369*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 0.7745966692414833*fc[7]-1.5*fc[6]+0.4472135954999579*fc[5]+1.118033988749895*fc[4]-1.161895003862225*fc[3]-0.6708203932499369*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-0.7745966692414833*fr[7])-1.5*fr[6]+0.4472135954999579*fr[5]+1.118033988749895*fr[4]+1.161895003862225*fr[3]-0.6708203932499369*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[0] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[0] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = 1.936491673103709*hamil[6]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = 1.936491673103709*hamil[6]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = (-0.9682458365518543*fl[7])-0.5590169943749475*fl[5]+1.118033988749895*fl[4]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = 0.9682458365518543*fc[7]-0.5590169943749475*fc[5]+1.118033988749895*fc[4]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = (-0.9682458365518543*fc[7])-0.5590169943749475*fc[5]+1.118033988749895*fc[4]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = 0.9682458365518543*fr[7]-0.5590169943749475*fr[5]+1.118033988749895*fr[4]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[1] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[1] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  alpha_l_quad = (-4.499999999999999*hamil[7])+1.936491673103709*hamil[6]+2.598076211353316*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_r_quad = 4.499999999999999*hamil[7]+1.936491673103709*hamil[6]+2.598076211353316*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  alpha_max = fmax(alpha_max, fmax(fabs(alpha_l_quad), fabs(alpha_r_quad))); 
  f_lr_quad = 0.7745966692414833*fl[7]+1.5*fl[6]+0.4472135954999579*fl[5]+1.118033988749895*fl[4]+1.161895003862225*fl[3]+0.6708203932499369*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  f_cl_quad = (-0.7745966692414833*fc[7])+1.5*fc[6]+0.4472135954999579*fc[5]+1.118033988749895*fc[4]-1.161895003862225*fc[3]+0.6708203932499369*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_cr_quad = 0.7745966692414833*fc[7]+1.5*fc[6]+0.4472135954999579*fc[5]+1.118033988749895*fc[4]+1.161895003862225*fc[3]+0.6708203932499369*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  f_rl_quad = (-0.7745966692414833*fr[7])+1.5*fr[6]+0.4472135954999579*fr[5]+1.118033988749895*fr[4]-1.161895003862225*fr[3]+0.6708203932499369*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  Fhat_l_nodal[2] = 0.5*(alpha_l_quad*(f_cl_quad + f_lr_quad) - fabs(alpha_l_quad)*(f_cl_quad - f_lr_quad)); 
  Fhat_r_nodal[2] = 0.5*(alpha_r_quad*(f_rl_quad + f_cr_quad) - fabs(alpha_r_quad)*(f_rl_quad - f_cr_quad)); 

  out[0] += ((-0.2777777777777778*Fhat_r_nodal[2])+0.2777777777777778*Fhat_l_nodal[2]-0.4444444444444444*Fhat_r_nodal[1]+0.4444444444444444*Fhat_l_nodal[1]-0.2777777777777778*Fhat_r_nodal[0]+0.2777777777777778*Fhat_l_nodal[0])*dv10*dx10; 
  out[1] += ((-0.4811252243246882*(Fhat_r_nodal[2]+Fhat_l_nodal[2]))-0.7698003589195012*(Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.4811252243246882*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv10*dx10; 
  out[2] += ((-0.3726779962499649*Fhat_r_nodal[2])+0.3726779962499649*(Fhat_l_nodal[2]+Fhat_r_nodal[0])-0.3726779962499649*Fhat_l_nodal[0])*dv10*dx10; 
  out[3] += (0.6454972243679029*(Fhat_r_nodal[0]+Fhat_l_nodal[0])-0.6454972243679029*(Fhat_r_nodal[2]+Fhat_l_nodal[2]))*dv10*dx10; 
  out[4] += ((-0.6211299937499417*Fhat_r_nodal[2])+0.6211299937499417*Fhat_l_nodal[2]-0.9938079899999065*Fhat_r_nodal[1]+0.9938079899999065*Fhat_l_nodal[1]-0.6211299937499417*Fhat_r_nodal[0]+0.6211299937499417*Fhat_l_nodal[0])*dv10*dx10; 
  out[5] += ((-0.2484519974999766*Fhat_r_nodal[2])+0.2484519974999766*Fhat_l_nodal[2]+0.4969039949999532*Fhat_r_nodal[1]-0.4969039949999532*Fhat_l_nodal[1]-0.2484519974999766*Fhat_r_nodal[0]+0.2484519974999766*Fhat_l_nodal[0])*dv10*dx10; 
  out[6] += ((-0.8333333333333335*Fhat_r_nodal[2])+0.8333333333333335*(Fhat_l_nodal[2]+Fhat_r_nodal[0])-0.8333333333333335*Fhat_l_nodal[0])*dv10*dx10; 
  out[7] += ((-0.4303314829119352*(Fhat_r_nodal[2]+Fhat_l_nodal[2]))+0.8606629658238704*(Fhat_r_nodal[1]+Fhat_l_nodal[1])-0.4303314829119352*(Fhat_r_nodal[0]+Fhat_l_nodal[0]))*dv10*dx10; 

  return fabs(2.5*dv10*dx10*alpha_max);

} 
