#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p2(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  out[0] += ((-0.2777777777777778*Fhat_r_nodal_c[2])+0.2777777777777778*Fhat_l_nodal_c[2]-0.4444444444444444*Fhat_r_nodal_c[1]+0.4444444444444444*Fhat_l_nodal_c[1]-0.2777777777777778*Fhat_r_nodal_c[0]+0.2777777777777778*Fhat_l_nodal_c[0])*dv10; 
  out[1] += ((-0.3726779962499649*Fhat_r_nodal_c[2])+0.3726779962499649*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[0])-0.3726779962499649*Fhat_l_nodal_c[0])*dv10; 
  out[2] += ((-0.4811252243246882*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]))-0.7698003589195012*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.4811252243246882*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[3] += (0.6454972243679029*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0])-0.6454972243679029*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]))*dv10; 
  out[4] += ((-0.2484519974999766*Fhat_r_nodal_c[2])+0.2484519974999766*Fhat_l_nodal_c[2]+0.4969039949999532*Fhat_r_nodal_c[1]-0.4969039949999532*Fhat_l_nodal_c[1]-0.2484519974999766*Fhat_r_nodal_c[0]+0.2484519974999766*Fhat_l_nodal_c[0])*dv10; 
  out[5] += ((-0.6211299937499417*Fhat_r_nodal_c[2])+0.6211299937499417*Fhat_l_nodal_c[2]-0.9938079899999065*Fhat_r_nodal_c[1]+0.9938079899999065*Fhat_l_nodal_c[1]-0.6211299937499417*Fhat_r_nodal_c[0]+0.6211299937499417*Fhat_l_nodal_c[0])*dv10; 
  out[6] += ((-0.4303314829119352*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]))+0.8606629658238704*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.4303314829119352*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[7] += ((-0.8333333333333335*Fhat_r_nodal_c[2])+0.8333333333333335*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[0])-0.8333333333333335*Fhat_l_nodal_c[0])*dv10; 
  out[8] += ((-0.5555555555555556*Fhat_r_nodal_c[2])+0.5555555555555556*Fhat_l_nodal_c[2]+1.111111111111111*Fhat_r_nodal_c[1]-1.111111111111111*Fhat_l_nodal_c[1]-0.5555555555555556*Fhat_r_nodal_c[0]+0.5555555555555556*Fhat_l_nodal_c[0])*dv10; 

  return 0.0;

} 
