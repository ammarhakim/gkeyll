#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Fhat_nodal[3]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = (-4.5*hamil[7])+1.936491673103709*hamil[6]-2.598076211353316*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = 0.7745966692414833*fskin[7]-1.5*fskin[6]+0.4472135954999579*fskin[5]+1.118033988749895*fskin[4]-1.161895003862225*fskin[3]-0.6708203932499369*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.7745966692414833*fedge[7])-1.5*fedge[6]+0.4472135954999579*fedge[5]+1.118033988749895*fedge[4]+1.161895003862225*fedge[3]-0.6708203932499369*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.936491673103709*hamil[6]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = (-0.9682458365518543*fskin[7])-0.5590169943749475*fskin[5]+1.118033988749895*fskin[4]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.9682458365518543*fedge[7]-0.5590169943749475*fedge[5]+1.118033988749895*fedge[4]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 4.5*hamil[7]+1.936491673103709*hamil[6]+2.598076211353316*hamil[5]+1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = 0.7745966692414833*fskin[7]+1.5*fskin[6]+0.4472135954999579*fskin[5]+1.118033988749895*fskin[4]+1.161895003862225*fskin[3]+0.6708203932499369*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = (-0.7745966692414833*fedge[7])+1.5*fedge[6]+0.4472135954999579*fedge[5]+1.118033988749895*fedge[4]-1.161895003862225*fedge[3]+0.6708203932499369*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += ((-0.2777777777777778*Fhat_nodal[2])-0.4444444444444444*Fhat_nodal[1]-0.2777777777777778*Fhat_nodal[0])*dv10*dx10; 
  out[1] += ((-0.4811252243246882*Fhat_nodal[2])-0.7698003589195012*Fhat_nodal[1]-0.4811252243246882*Fhat_nodal[0])*dv10*dx10; 
  out[2] += (0.3726779962499649*Fhat_nodal[0]-0.3726779962499649*Fhat_nodal[2])*dv10*dx10; 
  out[3] += (0.6454972243679029*Fhat_nodal[0]-0.6454972243679029*Fhat_nodal[2])*dv10*dx10; 
  out[4] += ((-0.6211299937499417*Fhat_nodal[2])-0.9938079899999065*Fhat_nodal[1]-0.6211299937499417*Fhat_nodal[0])*dv10*dx10; 
  out[5] += ((-0.2484519974999766*Fhat_nodal[2])+0.4969039949999532*Fhat_nodal[1]-0.2484519974999766*Fhat_nodal[0])*dv10*dx10; 
  out[6] += (0.8333333333333335*Fhat_nodal[0]-0.8333333333333335*Fhat_nodal[2])*dv10*dx10; 
  out[7] += ((-0.4303314829119352*Fhat_nodal[2])+0.8606629658238704*Fhat_nodal[1]-0.4303314829119352*Fhat_nodal[0])*dv10*dx10; 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 4.5*hamil[7]+1.936491673103709*hamil[6]-2.598076211353316*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = 0.7745966692414833*fedge[7]-1.5*fedge[6]+0.4472135954999579*fedge[5]+1.118033988749895*fedge[4]-1.161895003862225*fedge[3]-0.6708203932499369*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.7745966692414833*fskin[7])-1.5*fskin[6]+0.4472135954999579*fskin[5]+1.118033988749895*fskin[4]+1.161895003862225*fskin[3]-0.6708203932499369*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Fhat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.936491673103709*hamil[6]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = (-0.9682458365518543*fedge[7])-0.5590169943749475*fedge[5]+1.118033988749895*fedge[4]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.9682458365518543*fskin[7]-0.5590169943749475*fskin[5]+1.118033988749895*fskin[4]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Fhat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = (-4.5*hamil[7])+1.936491673103709*hamil[6]+2.598076211353316*hamil[5]-1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = 0.7745966692414833*fedge[7]+1.5*fedge[6]+0.4472135954999579*fedge[5]+1.118033988749895*fedge[4]+1.161895003862225*fedge[3]+0.6708203932499369*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = (-0.7745966692414833*fskin[7])+1.5*fskin[6]+0.4472135954999579*fskin[5]+1.118033988749895*fskin[4]-1.161895003862225*fskin[3]+0.6708203932499369*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Fhat_nodal[2] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  out[0] += (0.2777777777777778*Fhat_nodal[2]+0.4444444444444444*Fhat_nodal[1]+0.2777777777777778*Fhat_nodal[0])*dv10*dx10; 
  out[1] += ((-0.4811252243246882*Fhat_nodal[2])-0.7698003589195012*Fhat_nodal[1]-0.4811252243246882*Fhat_nodal[0])*dv10*dx10; 
  out[2] += (0.3726779962499649*Fhat_nodal[2]-0.3726779962499649*Fhat_nodal[0])*dv10*dx10; 
  out[3] += (0.6454972243679029*Fhat_nodal[0]-0.6454972243679029*Fhat_nodal[2])*dv10*dx10; 
  out[4] += (0.6211299937499417*Fhat_nodal[2]+0.9938079899999065*Fhat_nodal[1]+0.6211299937499417*Fhat_nodal[0])*dv10*dx10; 
  out[5] += (0.2484519974999766*Fhat_nodal[2]-0.4969039949999532*Fhat_nodal[1]+0.2484519974999766*Fhat_nodal[0])*dv10*dx10; 
  out[6] += (0.8333333333333335*Fhat_nodal[2]-0.8333333333333335*Fhat_nodal[0])*dv10*dx10; 
  out[7] += ((-0.4303314829119352*Fhat_nodal[2])+0.8606629658238704*Fhat_nodal[1]-0.4303314829119352*Fhat_nodal[0])*dv10*dx10; 

  } 
  return 0.0;

} 
