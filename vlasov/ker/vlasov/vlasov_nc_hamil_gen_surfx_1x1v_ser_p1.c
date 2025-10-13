#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  out[0] += ((-0.5*Fhat_r_nodal_c[1])+0.5*Fhat_l_nodal_c[1]-0.5*Fhat_r_nodal_c[0]+0.5*Fhat_l_nodal_c[0])*dx10; 
  out[1] += -0.8660254037844386*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0])*dx10; 
  out[2] += ((-0.5*Fhat_r_nodal_c[1])+0.5*(Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0])-0.5*Fhat_l_nodal_c[0])*dx10; 
  out[3] += (0.8660254037844386*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0])-0.8660254037844386*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1]))*dx10; 

  return 0.0;

} 
