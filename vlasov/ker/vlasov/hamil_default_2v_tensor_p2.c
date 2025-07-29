#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_2v_tensor_p2(const double *w, const double *dxv, double* GKYL_RESTRICT hamil) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // hamil: Particle Hamiltonian.
 
  const double wx1 = w[0], dv1 = dxv[0]; 
  const double wx2 = w[1], dv2 = dxv[1]; 
  double hamil_nodal[9] = {0.0};
  hamil_nodal[0] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_nodal[1] = 0.5*(pow(wx1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_nodal[2] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_nodal[3] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2, 2.0));
  hamil_nodal[4] = 0.5*(pow(wx1, 2.0) + pow(wx2, 2.0));
  hamil_nodal[5] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2, 2.0));
  hamil_nodal[6] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_nodal[7] = 0.5*(pow(wx1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_nodal[8] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));

  hamil[0] = 0.05555555555555555*hamil_nodal[8]+0.2222222222222222*hamil_nodal[7]+0.05555555555555555*hamil_nodal[6]+0.2222222222222222*hamil_nodal[5]+0.8888888888888888*hamil_nodal[4]+0.2222222222222222*hamil_nodal[3]+0.05555555555555555*hamil_nodal[2]+0.2222222222222222*hamil_nodal[1]+0.05555555555555555*hamil_nodal[0]; 
  hamil[1] = 0.09622504486493764*hamil_nodal[8]-0.09622504486493764*hamil_nodal[6]+0.3849001794597506*hamil_nodal[5]-0.3849001794597506*hamil_nodal[3]+0.09622504486493764*hamil_nodal[2]-0.09622504486493764*hamil_nodal[0]; 
  hamil[2] = 0.09622504486493764*hamil_nodal[8]+0.3849001794597506*hamil_nodal[7]+0.09622504486493764*hamil_nodal[6]-0.09622504486493764*hamil_nodal[2]-0.3849001794597506*hamil_nodal[1]-0.09622504486493764*hamil_nodal[0]; 
  hamil[3] = 0.16666666666666666*hamil_nodal[8]-0.16666666666666666*hamil_nodal[6]-0.16666666666666666*hamil_nodal[2]+0.16666666666666666*hamil_nodal[0]; 
  hamil[4] = 0.04969039949999532*hamil_nodal[8]-0.09938079899999064*hamil_nodal[7]+0.04969039949999532*hamil_nodal[6]+0.19876159799998128*hamil_nodal[5]-0.39752319599996255*hamil_nodal[4]+0.19876159799998128*hamil_nodal[3]+0.04969039949999532*hamil_nodal[2]-0.09938079899999064*hamil_nodal[1]+0.04969039949999532*hamil_nodal[0]; 
  hamil[5] = 0.04969039949999532*hamil_nodal[8]+0.19876159799998128*hamil_nodal[7]+0.04969039949999532*hamil_nodal[6]-0.09938079899999064*hamil_nodal[5]-0.39752319599996255*hamil_nodal[4]-0.09938079899999064*hamil_nodal[3]+0.04969039949999532*hamil_nodal[2]+0.19876159799998128*hamil_nodal[1]+0.04969039949999532*hamil_nodal[0]; 
  hamil[6] = 0.08606629658238703*hamil_nodal[8]-0.17213259316477406*hamil_nodal[7]+0.08606629658238703*hamil_nodal[6]-0.08606629658238703*hamil_nodal[2]+0.17213259316477406*hamil_nodal[1]-0.08606629658238703*hamil_nodal[0]; 
  hamil[7] = 0.08606629658238703*hamil_nodal[8]-0.08606629658238703*hamil_nodal[6]-0.17213259316477406*hamil_nodal[5]+0.17213259316477406*hamil_nodal[3]+0.08606629658238703*hamil_nodal[2]-0.08606629658238703*hamil_nodal[0]; 
  hamil[8] = 0.044444444444444446*hamil_nodal[8]-0.08888888888888889*hamil_nodal[7]+0.044444444444444446*hamil_nodal[6]-0.08888888888888889*hamil_nodal[5]+0.17777777777777778*hamil_nodal[4]-0.08888888888888889*hamil_nodal[3]+0.044444444444444446*hamil_nodal[2]-0.08888888888888889*hamil_nodal[1]+0.044444444444444446*hamil_nodal[0]; 

} 
