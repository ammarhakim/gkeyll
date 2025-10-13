#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += ((-0.2777777777777778*Fhat_nodal[2])-0.4444444444444444*Fhat_nodal[1]-0.2777777777777778*Fhat_nodal[0])*dx10; 
  out[1] += ((-0.4811252243246882*Fhat_nodal[2])-0.7698003589195012*Fhat_nodal[1]-0.4811252243246882*Fhat_nodal[0])*dx10; 
  out[2] += (0.3726779962499649*Fhat_nodal[0]-0.3726779962499649*Fhat_nodal[2])*dx10; 
  out[3] += (0.6454972243679029*Fhat_nodal[0]-0.6454972243679029*Fhat_nodal[2])*dx10; 
  out[4] += ((-0.6211299937499417*Fhat_nodal[2])-0.9938079899999065*Fhat_nodal[1]-0.6211299937499417*Fhat_nodal[0])*dx10; 
  out[5] += ((-0.2484519974999766*Fhat_nodal[2])+0.4969039949999532*Fhat_nodal[1]-0.2484519974999766*Fhat_nodal[0])*dx10; 
  out[6] += (0.8333333333333335*Fhat_nodal[0]-0.8333333333333335*Fhat_nodal[2])*dx10; 
  out[7] += ((-0.4303314829119352*Fhat_nodal[2])+0.8606629658238704*Fhat_nodal[1]-0.4303314829119352*Fhat_nodal[0])*dx10; 

  } else { 

  out[0] += (0.2777777777777778*Fhat_nodal[2]+0.4444444444444444*Fhat_nodal[1]+0.2777777777777778*Fhat_nodal[0])*dx10; 
  out[1] += ((-0.4811252243246882*Fhat_nodal[2])-0.7698003589195012*Fhat_nodal[1]-0.4811252243246882*Fhat_nodal[0])*dx10; 
  out[2] += (0.3726779962499649*Fhat_nodal[2]-0.3726779962499649*Fhat_nodal[0])*dx10; 
  out[3] += (0.6454972243679029*Fhat_nodal[0]-0.6454972243679029*Fhat_nodal[2])*dx10; 
  out[4] += (0.6211299937499417*Fhat_nodal[2]+0.9938079899999065*Fhat_nodal[1]+0.6211299937499417*Fhat_nodal[0])*dx10; 
  out[5] += (0.2484519974999766*Fhat_nodal[2]-0.4969039949999532*Fhat_nodal[1]+0.2484519974999766*Fhat_nodal[0])*dx10; 
  out[6] += (0.8333333333333335*Fhat_nodal[2]-0.8333333333333335*Fhat_nodal[0])*dx10; 
  out[7] += ((-0.4303314829119352*Fhat_nodal[2])+0.8606629658238704*Fhat_nodal[1]-0.4303314829119352*Fhat_nodal[0])*dx10; 

  } 
  return 0.0;

} 
