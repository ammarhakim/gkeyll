#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_default_3v_ser_p1(const double *w, const double *dxv, double* GKYL_RESTRICT hamil) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // hamil: Particle Hamiltonian.
 
  const double wx1 = w[0], dv1 = dxv[0]; 
  const double wx2 = w[1], dv2 = dxv[1]; 
  const double wx3 = w[2], dv3 = dxv[2]; 
  double hamil_nodal[8] = {0.0};
  hamil_nodal[0] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_nodal[1] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_nodal[2] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_nodal[3] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3-0.5*dv3, 2.0));
  hamil_nodal[4] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_nodal[5] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2-0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_nodal[6] = 0.5*(pow(wx1-0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));
  hamil_nodal[7] = 0.5*(pow(wx1+0.5*dv1, 2.0) + pow(wx2+0.5*dv2, 2.0) + pow(wx3+0.5*dv3, 2.0));

  hamil[0] = 0.3535533905932737*hamil_nodal[7]+0.3535533905932737*hamil_nodal[6]+0.3535533905932737*hamil_nodal[5]+0.3535533905932737*hamil_nodal[4]+0.3535533905932737*hamil_nodal[3]+0.3535533905932737*hamil_nodal[2]+0.3535533905932737*hamil_nodal[1]+0.3535533905932737*hamil_nodal[0]; 
  hamil[1] = 0.20412414523193148*hamil_nodal[7]-0.20412414523193148*hamil_nodal[6]+0.20412414523193148*hamil_nodal[5]-0.20412414523193148*hamil_nodal[4]+0.20412414523193148*hamil_nodal[3]-0.20412414523193148*hamil_nodal[2]+0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[2] = 0.20412414523193148*hamil_nodal[7]+0.20412414523193148*hamil_nodal[6]-0.20412414523193148*hamil_nodal[5]-0.20412414523193148*hamil_nodal[4]+0.20412414523193148*hamil_nodal[3]+0.20412414523193148*hamil_nodal[2]-0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[3] = 0.20412414523193148*hamil_nodal[7]+0.20412414523193148*hamil_nodal[6]+0.20412414523193148*hamil_nodal[5]+0.20412414523193148*hamil_nodal[4]-0.20412414523193148*hamil_nodal[3]-0.20412414523193148*hamil_nodal[2]-0.20412414523193148*hamil_nodal[1]-0.20412414523193148*hamil_nodal[0]; 
  hamil[4] = 0.11785113019775789*hamil_nodal[7]-0.11785113019775789*hamil_nodal[6]-0.11785113019775789*hamil_nodal[5]+0.11785113019775789*hamil_nodal[4]+0.11785113019775789*hamil_nodal[3]-0.11785113019775789*hamil_nodal[2]-0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[5] = 0.11785113019775789*hamil_nodal[7]-0.11785113019775789*hamil_nodal[6]+0.11785113019775789*hamil_nodal[5]-0.11785113019775789*hamil_nodal[4]-0.11785113019775789*hamil_nodal[3]+0.11785113019775789*hamil_nodal[2]-0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[6] = 0.11785113019775789*hamil_nodal[7]+0.11785113019775789*hamil_nodal[6]-0.11785113019775789*hamil_nodal[5]-0.11785113019775789*hamil_nodal[4]-0.11785113019775789*hamil_nodal[3]-0.11785113019775789*hamil_nodal[2]+0.11785113019775789*hamil_nodal[1]+0.11785113019775789*hamil_nodal[0]; 
  hamil[7] = 0.06804138174397717*hamil_nodal[7]-0.06804138174397717*hamil_nodal[6]-0.06804138174397717*hamil_nodal[5]+0.06804138174397717*hamil_nodal[4]-0.06804138174397717*hamil_nodal[3]+0.06804138174397717*hamil_nodal[2]+0.06804138174397717*hamil_nodal[1]-0.06804138174397717*hamil_nodal[0]; 

} 
