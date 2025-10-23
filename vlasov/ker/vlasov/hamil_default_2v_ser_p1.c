#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_2v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Utilized by relativistic simulations. 
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double wx2 = w[1], dv2 = dxv[1]; 
  double hamil_nodal[4] = {0.0};
  double hamil_inv_nodal[4] = {0.0};
  hamil_nodal[0] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];

  hamil[0] = 0.5*hamil_nodal[3]+0.5*hamil_nodal[2]+0.5*hamil_nodal[1]+0.5*hamil_nodal[0]; 
  hamil[1] = 0.2886751345948129*hamil_nodal[3]-0.2886751345948129*hamil_nodal[2]+0.2886751345948129*hamil_nodal[1]-0.2886751345948129*hamil_nodal[0]; 
  hamil[2] = 0.2886751345948129*hamil_nodal[3]+0.2886751345948129*hamil_nodal[2]-0.2886751345948129*hamil_nodal[1]-0.2886751345948129*hamil_nodal[0]; 
  hamil[3] = 0.1666666666666667*hamil_nodal[3]-0.1666666666666667*hamil_nodal[2]-0.1666666666666667*hamil_nodal[1]+0.1666666666666667*hamil_nodal[0]; 

  hamil_inv[0] = 0.5*hamil_inv_nodal[3]+0.5*hamil_inv_nodal[2]+0.5*hamil_inv_nodal[1]+0.5*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.2886751345948129*hamil_inv_nodal[3]-0.2886751345948129*hamil_inv_nodal[2]+0.2886751345948129*hamil_inv_nodal[1]-0.2886751345948129*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.2886751345948129*hamil_inv_nodal[3]+0.2886751345948129*hamil_inv_nodal[2]-0.2886751345948129*hamil_inv_nodal[1]-0.2886751345948129*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.1666666666666667*hamil_inv_nodal[3]-0.1666666666666667*hamil_inv_nodal[2]-0.1666666666666667*hamil_inv_nodal[1]+0.1666666666666667*hamil_inv_nodal[0]; 

} 
