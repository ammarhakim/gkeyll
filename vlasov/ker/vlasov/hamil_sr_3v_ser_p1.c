#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_sr_3v_ser_p1(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT hamil, double* GKYL_RESTRICT hamil_inv) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // hamil: Particle Hamiltonian.
  // hamil_inv: Inverse particle Hamiltonian. Utilized by relativistic simulations. 
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double wx2 = w[1], dv2 = dxv[1]; 
  double wx3 = w[2], dv3 = dxv[2]; 
  double hamil_nodal[8] = {0.0};
  double hamil_inv_nodal[8] = {0.0};
  hamil_nodal[0] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_inv_nodal[0] = 1.0/hamil_nodal[0];
  hamil_nodal[1] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_inv_nodal[1] = 1.0/hamil_nodal[1];
  hamil_nodal[2] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_inv_nodal[2] = 1.0/hamil_nodal[2];
  hamil_nodal[3] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_inv_nodal[3] = 1.0/hamil_nodal[3];
  hamil_nodal[4] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_inv_nodal[4] = 1.0/hamil_nodal[4];
  hamil_nodal[5] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_inv_nodal[5] = 1.0/hamil_nodal[5];
  hamil_nodal[6] = sqrt(1.0 + pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_inv_nodal[6] = 1.0/hamil_nodal[6];
  hamil_nodal[7] = sqrt(1.0 + pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_inv_nodal[7] = 1.0/hamil_nodal[7];

  hamil[0] = 0.3535533905932737*hamil_nodal[7]+0.3535533905932737*hamil_nodal[6]+0.3535533905932737*hamil_nodal[5]+0.3535533905932737*hamil_nodal[4]+0.3535533905932737*hamil_nodal[3]+0.3535533905932737*hamil_nodal[2]+0.3535533905932737*hamil_nodal[1]+0.3535533905932737*hamil_nodal[0]; 
  hamil[1] = 0.20412414523193148*hamil_nodal[7]-0.20412414523193148*hamil_nodal[6]+0.20412414523193148*hamil_nodal[5]-0.20412414523193148*hamil_nodal[4]+0.20412414523193148*hamil_nodal[3]-0.20412414523193148*hamil_nodal[2]+0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[2] = 0.20412414523193148*hamil_nodal[7]+0.20412414523193148*hamil_nodal[6]-0.20412414523193148*hamil_nodal[5]-0.20412414523193148*hamil_nodal[4]+0.20412414523193148*hamil_nodal[3]+0.20412414523193148*hamil_nodal[2]-0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[3] = 0.20412414523193148*hamil_nodal[7]+0.20412414523193148*hamil_nodal[6]+0.20412414523193148*hamil_nodal[5]+0.20412414523193148*hamil_nodal[4]-0.20412414523193148*hamil_nodal[3]-0.20412414523193148*hamil_nodal[2]-0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[4] = 0.11785113019775789*hamil_nodal[7]-0.11785113019775789*hamil_nodal[6]-0.11785113019775789*hamil_nodal[5]+0.11785113019775789*hamil_nodal[4]+0.11785113019775789*hamil_nodal[3]-0.11785113019775789*hamil_nodal[2]-0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[5] = 0.11785113019775789*hamil_nodal[7]-0.11785113019775789*hamil_nodal[6]+0.11785113019775789*hamil_nodal[5]-0.11785113019775789*hamil_nodal[4]-0.11785113019775789*hamil_nodal[3]+0.11785113019775789*hamil_nodal[2]-0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[6] = 0.11785113019775789*hamil_nodal[7]+0.11785113019775789*hamil_nodal[6]-0.11785113019775789*hamil_nodal[5]-0.11785113019775789*hamil_nodal[4]-0.11785113019775789*hamil_nodal[3]-0.11785113019775789*hamil_nodal[2]+0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[7] = 0.06804138174397717*hamil_nodal[7]-0.06804138174397717*hamil_nodal[6]-0.06804138174397717*hamil_nodal[5]+0.06804138174397717*hamil_nodal[4]-0.06804138174397717*hamil_nodal[3]+0.06804138174397717*hamil_nodal[2]+0.06804138174397717*hamil_nodal[1]-0.06804138174397717*hamil_nodal[0]; 

  hamil_inv[0] = 0.3535533905932737*hamil_inv_nodal[7]+0.3535533905932737*hamil_inv_nodal[6]+0.3535533905932737*hamil_inv_nodal[5]+0.3535533905932737*hamil_inv_nodal[4]+0.3535533905932737*hamil_inv_nodal[3]+0.3535533905932737*hamil_inv_nodal[2]+0.3535533905932737*hamil_inv_nodal[1]+0.3535533905932737*hamil_inv_nodal[0]; 
  hamil_inv[1] = 0.20412414523193148*hamil_inv_nodal[7]-0.20412414523193148*hamil_inv_nodal[6]+0.20412414523193148*hamil_inv_nodal[5]-0.20412414523193148*hamil_inv_nodal[4]+0.20412414523193148*hamil_inv_nodal[3]-0.20412414523193148*hamil_inv_nodal[2]+0.20412414523193148*hamil_inv_nodal[1]-0.20412414523193148*hamil_inv_nodal[0]; 
  hamil_inv[2] = 0.20412414523193148*hamil_inv_nodal[7]+0.20412414523193148*hamil_inv_nodal[6]-0.20412414523193148*hamil_inv_nodal[5]-0.20412414523193148*hamil_inv_nodal[4]+0.20412414523193148*hamil_inv_nodal[3]+0.20412414523193148*hamil_inv_nodal[2]-0.20412414523193148*hamil_inv_nodal[1]-0.20412414523193148*hamil_inv_nodal[0]; 
  hamil_inv[3] = 0.20412414523193148*hamil_inv_nodal[7]+0.20412414523193148*hamil_inv_nodal[6]+0.20412414523193148*hamil_inv_nodal[5]+0.20412414523193148*hamil_inv_nodal[4]-0.20412414523193148*hamil_inv_nodal[3]-0.20412414523193148*hamil_inv_nodal[2]-0.20412414523193148*hamil_inv_nodal[1]-0.20412414523193148*hamil_inv_nodal[0]; 
  hamil_inv[4] = 0.11785113019775789*hamil_inv_nodal[7]-0.11785113019775789*hamil_inv_nodal[6]-0.11785113019775789*hamil_inv_nodal[5]+0.11785113019775789*hamil_inv_nodal[4]+0.11785113019775789*hamil_inv_nodal[3]-0.11785113019775789*hamil_inv_nodal[2]-0.11785113019775789*hamil_inv_nodal[1]+0.11785113019775789*hamil_inv_nodal[0]; 
  hamil_inv[5] = 0.11785113019775789*hamil_inv_nodal[7]-0.11785113019775789*hamil_inv_nodal[6]+0.11785113019775789*hamil_inv_nodal[5]-0.11785113019775789*hamil_inv_nodal[4]-0.11785113019775789*hamil_inv_nodal[3]+0.11785113019775789*hamil_inv_nodal[2]-0.11785113019775789*hamil_inv_nodal[1]+0.11785113019775789*hamil_inv_nodal[0]; 
  hamil_inv[6] = 0.11785113019775789*hamil_inv_nodal[7]+0.11785113019775789*hamil_inv_nodal[6]-0.11785113019775789*hamil_inv_nodal[5]-0.11785113019775789*hamil_inv_nodal[4]-0.11785113019775789*hamil_inv_nodal[3]-0.11785113019775789*hamil_inv_nodal[2]+0.11785113019775789*hamil_inv_nodal[1]+0.11785113019775789*hamil_inv_nodal[0]; 
  hamil_inv[7] = 0.06804138174397717*hamil_inv_nodal[7]-0.06804138174397717*hamil_inv_nodal[6]-0.06804138174397717*hamil_inv_nodal[5]+0.06804138174397717*hamil_inv_nodal[4]-0.06804138174397717*hamil_inv_nodal[3]+0.06804138174397717*hamil_inv_nodal[2]+0.06804138174397717*hamil_inv_nodal[1]-0.06804138174397717*hamil_inv_nodal[0]; 

} 
