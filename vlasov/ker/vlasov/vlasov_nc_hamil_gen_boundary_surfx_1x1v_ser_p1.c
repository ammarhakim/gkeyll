#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += -0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.5*Fhat_nodal[0]-0.5*Fhat_nodal[1])*dx10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dx10; 

  } else { 

  out[0] += 0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[1] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dx10; 
  out[2] += (0.5*Fhat_nodal[1]-0.5*Fhat_nodal[0])*dx10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dx10; 

  } 
  return 0.0;

} 
