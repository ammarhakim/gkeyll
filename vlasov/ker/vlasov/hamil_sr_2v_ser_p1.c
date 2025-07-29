#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_2v_ser_p1(const double *w, const double *dxv, double* GKYL_RESTRICT hamil) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // hamil: Particle Hamiltonian.
 
  const double wx1 = w[0], dv1 = dxv[0]; 
  const double wx2 = w[1], dv2 = dxv[1]; 
  double hamil_nodal[4] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_nodal[1] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_nodal[2] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_nodal[3] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));

  hamil[0] = 0.5*hamil_nodal[3]+0.5*hamil_nodal[2]+0.5*hamil_nodal[1]+0.5*hamil_nodal[0]; 
  hamil[1] = 0.2886751345948129*hamil_nodal[3]-0.2886751345948129*hamil_nodal[2]+0.2886751345948129*hamil_nodal[1]-0.2886751345948129*hamil_nodal[0]; 
  hamil[2] = 0.2886751345948129*hamil_nodal[3]+0.2886751345948129*hamil_nodal[2]-0.2886751345948129*hamil_nodal[1]-0.2886751345948129*hamil_nodal[0]; 
  hamil[3] = 0.16666666666666666*hamil_nodal[3]-0.16666666666666666*hamil_nodal[2]-0.16666666666666666*hamil_nodal[1]+0.16666666666666666*hamil_nodal[0]; 

} 
