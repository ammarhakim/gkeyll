#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_sol_stencil_2x_ser_p1(const double *sol_nodal_global, long perpOff, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: map between local nodes and global nodes.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.5*sol_nodal_global[perpOff+globalIdxs[3]]+0.5*sol_nodal_global[perpOff+globalIdxs[2]]+0.5*sol_nodal_global[perpOff+globalIdxs[1]]+0.5*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[1] = 0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[3]]-0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[2]]+0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[1]]-0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[2] = 0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[3]]+0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[2]]-0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[1]]-0.28867513459481287*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[3] = 0.16666666666666666*sol_nodal_global[perpOff+globalIdxs[3]]-0.16666666666666666*sol_nodal_global[perpOff+globalIdxs[2]]-0.16666666666666666*sol_nodal_global[perpOff+globalIdxs[1]]+0.16666666666666666*sol_nodal_global[perpOff+globalIdxs[0]];

}
GKYL_CU_DH void fem_poisson_perp_sol_stencil_3x_ser_p1(const double *sol_nodal_global, long perpOff, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // perpOff: memory offset due to other perpendicular planes (perp index * global number of nodes).
  // globalIdxs: map between local nodes and global nodes.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[7]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[6]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[5]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[4]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[3]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[2]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[1]]+0.3535533905932737*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[1] = 0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[7]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[6]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[5]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[4]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[3]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[2]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[1]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[2] = 0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[7]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[6]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[5]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[4]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[3]]+0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[2]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[1]]-0.20412414523193145*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[3] = 0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[7]]+0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[6]]+0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[5]]+0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[4]]-0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[3]]-0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[2]]-0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[1]]-0.2041241452319315*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[4] = 0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[7]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[6]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[5]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[4]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[3]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[2]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[1]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[5] = 0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[7]]-0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[6]]+0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[5]]-0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[4]]-0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[3]]+0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[2]]-0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[1]]+0.11785113019775792*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[6] = 0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[7]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[6]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[5]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[4]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[3]]-0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[2]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[1]]+0.11785113019775789*sol_nodal_global[perpOff+globalIdxs[0]];
  sol_modal_local[7] = 0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[7]]-0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[6]]-0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[5]]+0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[4]]-0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[3]]+0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[2]]+0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[1]]-0.06804138174397716*sol_nodal_global[perpOff+globalIdxs[0]];

}
