#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_lorentz_1v_ser_p2(const double *w, const double *dxv, double* GKYL_RESTRICT gamma, double* GKYL_RESTRICT gamma_inv) 
{ 
  // w:   Cell-center coordinates.
  // dxv: Cell spacing.
  // gamma:  Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv: Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
 
  const double w0 = w[0]; 
  const double dv0 = dxv[0]; 
  const double w0_sq = w0*w0, dv0_sq = dv0*dv0; 
  double p0_sq[3] = {0.0};
  p0_sq[0] = 1.414213562373095*w0_sq+0.1178511301977579*dv0_sq; 
  p0_sq[1] = 0.8164965809277261*dv0*w0; 
  p0_sq[2] = 0.105409255338946*dv0_sq; 

  double gamma_nodal[3] = {0.0};
  double gamma_inv_nodal[3] = {0.0};

  gamma_nodal[0] = sqrt(1.0 + 1.58113883008419*p0_sq[2]-1.224744871391589*p0_sq[1]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[0] = 1.0/gamma_nodal[0];
  gamma_nodal[1] = sqrt(1.0 + 0.7071067811865475*p0_sq[0]-0.7905694150420947*p0_sq[2]);
  gamma_inv_nodal[1] = 1.0/gamma_nodal[1];
  gamma_nodal[2] = sqrt(1.0 + 1.58113883008419*p0_sq[2]+1.224744871391589*p0_sq[1]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[2] = 1.0/gamma_nodal[2];

  gamma[0] = 0.2357022603955158*gamma_nodal[2]+0.9428090415820636*gamma_nodal[1]+0.2357022603955158*gamma_nodal[0]; 
  gamma[1] = 0.408248290463863*gamma_nodal[2]-0.408248290463863*gamma_nodal[0]; 
  gamma[2] = 0.210818510677892*gamma_nodal[2]-0.421637021355784*gamma_nodal[1]+0.210818510677892*gamma_nodal[0]; 

  gamma_inv[0] = 0.2357022603955158*gamma_inv_nodal[2]+0.9428090415820636*gamma_inv_nodal[1]+0.2357022603955158*gamma_inv_nodal[0]; 
  gamma_inv[1] = 0.408248290463863*gamma_inv_nodal[2]-0.408248290463863*gamma_inv_nodal[0]; 
  gamma_inv[2] = 0.210818510677892*gamma_inv_nodal[2]-0.421637021355784*gamma_inv_nodal[1]+0.210818510677892*gamma_inv_nodal[0]; 
} 
