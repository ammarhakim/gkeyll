#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_1v_ser_p3(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Utilized by relativistic simulations. 
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double hamil_nodal[4] = {0.0};
  double hamil_inv_nodal[4] = {0.0};
  hamil_nodal[0] = 0.5*(pow(wx1-0.5*dv1, 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = 0.5*(pow(wx1-0.22360679774997896*dv1, 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = 0.5*(pow(wx1+0.22360679774997896*dv1, 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = 0.5*(pow(wx1+0.5*dv1, 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];

  hamil[0] = 0.11785113019775789*hamil_nodal[3]+0.5892556509887895*hamil_nodal[2]+0.5892556509887895*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[1] = 0.20412414523193148*hamil_nodal[3]+0.4564354645876384*hamil_nodal[2]-0.4564354645876384*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[2] = 0.2635231383473649*hamil_nodal[3]-0.2635231383473649*hamil_nodal[2]-0.2635231383473649*hamil_nodal[1]+0.2635231383473649*hamil_nodal[0]; 
  hamil[3] = 0.13363062095621214*hamil_nodal[3]-0.29880715233359834*hamil_nodal[2]+0.29880715233359834*hamil_nodal[1]-0.13363062095621214*hamil_nodal[0]; 

  hamil_inv[0] = 0.11785113019775789*hamil_inv_nodal[3]+0.5892556509887895*hamil_inv_nodal[2]+0.5892556509887895*hamil_inv_nodal[1]+0.11785113019775789*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.20412414523193148*hamil_inv_nodal[3]+0.4564354645876384*hamil_inv_nodal[2]-0.4564354645876384*hamil_inv_nodal[1]-0.20412414523193148*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.2635231383473649*hamil_inv_nodal[3]-0.2635231383473649*hamil_inv_nodal[2]-0.2635231383473649*hamil_inv_nodal[1]+0.2635231383473649*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.13363062095621214*hamil_inv_nodal[3]-0.29880715233359834*hamil_inv_nodal[2]+0.29880715233359834*hamil_inv_nodal[1]-0.13363062095621214*hamil_inv_nodal[0]; 

} 
