#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_1v_tensor_p2(const double *w, const double *dxv, double* GKYL_RESTRICT hamil) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // hamil: Particle Hamiltonian.
 
  const double wx1 = w[0], dv1 = dxv[0]; 
  double hamil_nodal[3] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0));
  hamil_nodal[1] = sqrt(1.0 + pow(wx1, 2.0));
  hamil_nodal[2] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0));

  hamil[0] = 0.2357022603955158*hamil_nodal[2]+0.9428090415820636*hamil_nodal[1]+0.2357022603955158*hamil_nodal[0]; 
  hamil[1] = 0.408248290463863*hamil_nodal[2]-0.408248290463863*hamil_nodal[0]; 
  hamil[2] = 0.21081851067789195*hamil_nodal[2]-0.421637021355784*hamil_nodal[1]+0.21081851067789195*hamil_nodal[0]; 

} 
