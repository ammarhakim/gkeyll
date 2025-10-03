#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += -0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.5*Fhat_nodal[0]-0.5*Fhat_nodal[1])*dv10; 
  out[2] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dv10; 

  } else { 

  out[0] += 0.5*(Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[1] += (0.5*Fhat_nodal[1]-0.5*Fhat_nodal[0])*dv10; 
  out[2] += -0.8660254037844386*(Fhat_nodal[1]+Fhat_nodal[0])*dv10; 
  out[3] += (0.8660254037844386*Fhat_nodal[0]-0.8660254037844386*Fhat_nodal[1])*dv10; 

  } 
  return 0.0;

} 
