#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_2v_tensor_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Utilized by relativistic simulations. 
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double wx2 = w[1], dv2 = dxv[1]; 
  double hamil_nodal[9] = {0.0};
  double hamil_inv_nodal[9] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(wx1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2, 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];
  hamil_nodal[4] = sqrt(1.0 + pow(wx1, 2.0) + pow(wx2, 2.0));
  hamil_inv_nodal[4] = 1.0/hamil_nodal[4];
  hamil_nodal[5] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2, 2.0));
  hamil_inv_nodal[5] = 1.0/hamil_nodal[5];
  hamil_nodal[6] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_inv_nodal[6] = 1.0/hamil_nodal[6];
  hamil_nodal[7] = sqrt(1.0 + pow(wx1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_inv_nodal[7] = 1.0/hamil_nodal[7];
  hamil_nodal[8] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0));
  hamil_inv_nodal[8] = 1.0/hamil_nodal[8];

  hamil[0] = 0.05555555555555555*hamil_nodal[8]+0.2222222222222222*hamil_nodal[7]+0.05555555555555555*hamil_nodal[6]+0.2222222222222222*hamil_nodal[5]+0.8888888888888888*hamil_nodal[4]+0.2222222222222222*hamil_nodal[3]+0.05555555555555555*hamil_nodal[2]+0.2222222222222222*hamil_nodal[1]+0.05555555555555555*hamil_nodal[0]; 
  hamil[1] = 0.09622504486493764*hamil_nodal[8]-0.09622504486493764*hamil_nodal[6]+0.3849001794597506*hamil_nodal[5]-0.3849001794597506*hamil_nodal[3]+0.09622504486493764*hamil_nodal[2]-0.09622504486493764*hamil_nodal[0]; 
  hamil[2] = 0.09622504486493764*hamil_nodal[8]+0.3849001794597506*hamil_nodal[7]+0.09622504486493764*hamil_nodal[6]-0.09622504486493764*hamil_nodal[2]-0.3849001794597506*hamil_nodal[1]-0.09622504486493764*hamil_nodal[0]; 
  hamil[3] = 0.1666666666666667*hamil_nodal[8]-0.1666666666666667*hamil_nodal[6]-0.1666666666666667*hamil_nodal[2]+0.1666666666666667*hamil_nodal[0]; 
  hamil[4] = 0.04969039949999532*hamil_nodal[8]-0.09938079899999064*hamil_nodal[7]+0.04969039949999532*hamil_nodal[6]+0.1987615979999813*hamil_nodal[5]-0.3975231959999626*hamil_nodal[4]+0.1987615979999813*hamil_nodal[3]+0.04969039949999532*hamil_nodal[2]-0.09938079899999064*hamil_nodal[1]+0.04969039949999532*hamil_nodal[0]; 
  hamil[5] = 0.04969039949999532*hamil_nodal[8]+0.1987615979999813*hamil_nodal[7]+0.04969039949999532*hamil_nodal[6]-0.09938079899999064*hamil_nodal[5]-0.3975231959999626*hamil_nodal[4]-0.09938079899999064*hamil_nodal[3]+0.04969039949999532*hamil_nodal[2]+0.1987615979999813*hamil_nodal[1]+0.04969039949999532*hamil_nodal[0]; 
  hamil[6] = 0.08606629658238703*hamil_nodal[8]-0.1721325931647741*hamil_nodal[7]+0.08606629658238703*hamil_nodal[6]-0.08606629658238703*hamil_nodal[2]+0.1721325931647741*hamil_nodal[1]-0.08606629658238703*hamil_nodal[0]; 
  hamil[7] = 0.08606629658238703*hamil_nodal[8]-0.08606629658238703*hamil_nodal[6]-0.1721325931647741*hamil_nodal[5]+0.1721325931647741*hamil_nodal[3]+0.08606629658238703*hamil_nodal[2]-0.08606629658238703*hamil_nodal[0]; 
  hamil[8] = 0.04444444444444445*hamil_nodal[8]-0.08888888888888889*hamil_nodal[7]+0.04444444444444445*hamil_nodal[6]-0.08888888888888889*hamil_nodal[5]+0.1777777777777778*hamil_nodal[4]-0.08888888888888889*hamil_nodal[3]+0.04444444444444445*hamil_nodal[2]-0.08888888888888889*hamil_nodal[1]+0.04444444444444445*hamil_nodal[0]; 

  hamil_inv[0] = 0.05555555555555555*hamil_inv_nodal[8]+0.2222222222222222*hamil_inv_nodal[7]+0.05555555555555555*hamil_inv_nodal[6]+0.2222222222222222*hamil_inv_nodal[5]+0.8888888888888888*hamil_inv_nodal[4]+0.2222222222222222*hamil_inv_nodal[3]+0.05555555555555555*hamil_inv_nodal[2]+0.2222222222222222*hamil_inv_nodal[1]+0.05555555555555555*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.09622504486493764*hamil_inv_nodal[8]-0.09622504486493764*hamil_inv_nodal[6]+0.3849001794597506*hamil_inv_nodal[5]-0.3849001794597506*hamil_inv_nodal[3]+0.09622504486493764*hamil_inv_nodal[2]-0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.09622504486493764*hamil_inv_nodal[8]+0.3849001794597506*hamil_inv_nodal[7]+0.09622504486493764*hamil_inv_nodal[6]-0.09622504486493764*hamil_inv_nodal[2]-0.3849001794597506*hamil_inv_nodal[1]-0.09622504486493764*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.1666666666666667*hamil_inv_nodal[8]-0.1666666666666667*hamil_inv_nodal[6]-0.1666666666666667*hamil_inv_nodal[2]+0.1666666666666667*hamil_inv_nodal[0]; 
  hamil_inv[4] = 0.04969039949999532*hamil_inv_nodal[8]-0.09938079899999064*hamil_inv_nodal[7]+0.04969039949999532*hamil_inv_nodal[6]+0.1987615979999813*hamil_inv_nodal[5]-0.3975231959999626*hamil_inv_nodal[4]+0.1987615979999813*hamil_inv_nodal[3]+0.04969039949999532*hamil_inv_nodal[2]-0.09938079899999064*hamil_inv_nodal[1]+0.04969039949999532*hamil_inv_nodal[0]; 
  hamil_inv[5] = 0.04969039949999532*hamil_inv_nodal[8]+0.1987615979999813*hamil_inv_nodal[7]+0.04969039949999532*hamil_inv_nodal[6]-0.09938079899999064*hamil_inv_nodal[5]-0.3975231959999626*hamil_inv_nodal[4]-0.09938079899999064*hamil_inv_nodal[3]+0.04969039949999532*hamil_inv_nodal[2]+0.1987615979999813*hamil_inv_nodal[1]+0.04969039949999532*hamil_inv_nodal[0]; 
  hamil_inv[6] = 0.08606629658238703*hamil_inv_nodal[8]-0.1721325931647741*hamil_inv_nodal[7]+0.08606629658238703*hamil_inv_nodal[6]-0.08606629658238703*hamil_inv_nodal[2]+0.1721325931647741*hamil_inv_nodal[1]-0.08606629658238703*hamil_inv_nodal[0]; 
  hamil_inv[7] = 0.08606629658238703*hamil_inv_nodal[8]-0.08606629658238703*hamil_inv_nodal[6]-0.1721325931647741*hamil_inv_nodal[5]+0.1721325931647741*hamil_inv_nodal[3]+0.08606629658238703*hamil_inv_nodal[2]-0.08606629658238703*hamil_inv_nodal[0]; 
  hamil_inv[8] = 0.04444444444444445*hamil_inv_nodal[8]-0.08888888888888889*hamil_inv_nodal[7]+0.04444444444444445*hamil_inv_nodal[6]-0.08888888888888889*hamil_inv_nodal[5]+0.1777777777777778*hamil_inv_nodal[4]-0.08888888888888889*hamil_inv_nodal[3]+0.04444444444444445*hamil_inv_nodal[2]-0.08888888888888889*hamil_inv_nodal[1]+0.04444444444444445*hamil_inv_nodal[0]; 

} 
