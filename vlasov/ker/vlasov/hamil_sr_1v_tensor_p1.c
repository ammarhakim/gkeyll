#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_1v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Utilized by relativistic simulations. 
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double hamil_nodal[3] = {0.0};
  double hamil_inv_nodal[3] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(wx1, 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];

  hamil[0] = 0.2357022603955158*hamil_nodal[2]+0.9428090415820636*hamil_nodal[1]+0.2357022603955158*hamil_nodal[0]; 
  hamil[1] = 0.408248290463863*hamil_nodal[2]-0.408248290463863*hamil_nodal[0]; 
  hamil[2] = 0.210818510677892*hamil_nodal[2]-0.421637021355784*hamil_nodal[1]+0.210818510677892*hamil_nodal[0]; 

  hamil_inv[0] = 0.2357022603955158*hamil_inv_nodal[2]+0.9428090415820636*hamil_inv_nodal[1]+0.2357022603955158*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.408248290463863*hamil_inv_nodal[2]-0.408248290463863*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.210818510677892*hamil_inv_nodal[2]-0.421637021355784*hamil_inv_nodal[1]+0.210818510677892*hamil_inv_nodal[0]; 

} 
