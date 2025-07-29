#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_1v_tensor_p1(const double *w, const double *dxv, double* GKYL_RESTRICT hamil) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // hamil: Particle Hamiltonian.
 
  const double wx1 = w[0], dv1 = dxv[0]; 
  double hamil_nodal[2] = {0.0};
  hamil_nodal[0] = 0.5*(pow(wx1-0.5*dv1, 2.0));
  hamil_nodal[1] = 0.5*(pow(wx1+0.5*dv1, 2.0));

  hamil[0] = 0.7071067811865475*hamil_nodal[1]+0.7071067811865475*hamil_nodal[0]; 
  hamil[1] = 0.408248290463863*hamil_nodal[1]-0.408248290463863*hamil_nodal[0]; 

} 
