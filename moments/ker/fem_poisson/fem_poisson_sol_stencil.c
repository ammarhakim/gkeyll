#include <gkyl_fem_poisson_kernels.h> 
 
GKYL_CU_DH void fem_poisson_sol_stencil_1x_ser_p1(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.7071067811865476*sol_nodal_global[globalIdxs[1]]+0.7071067811865476*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = 0.408248290463863*sol_nodal_global[globalIdxs[1]]-0.408248290463863*sol_nodal_global[globalIdxs[0]];

}
GKYL_CU_DH void fem_poisson_sol_stencil_1x_ser_p2(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.2357022603955158*sol_nodal_global[globalIdxs[2]]+0.9428090415820637*sol_nodal_global[globalIdxs[1]]+0.2357022603955158*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = 0.408248290463863*sol_nodal_global[globalIdxs[2]]-0.408248290463863*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[2] = 0.210818510677892*sol_nodal_global[globalIdxs[2]]-0.421637021355784*sol_nodal_global[globalIdxs[1]]+0.210818510677892*sol_nodal_global[globalIdxs[0]];

}
GKYL_CU_DH void fem_poisson_sol_stencil_2x_ser_p1(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.5*sol_nodal_global[globalIdxs[3]]+0.5*sol_nodal_global[globalIdxs[2]]+0.5*sol_nodal_global[globalIdxs[1]]+0.5*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = 0.2886751345948129*sol_nodal_global[globalIdxs[3]]-0.2886751345948129*sol_nodal_global[globalIdxs[2]]+0.2886751345948129*sol_nodal_global[globalIdxs[1]]-0.2886751345948129*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[2] = 0.2886751345948129*sol_nodal_global[globalIdxs[3]]+0.2886751345948129*sol_nodal_global[globalIdxs[2]]-0.2886751345948129*sol_nodal_global[globalIdxs[1]]-0.2886751345948129*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[3] = 0.1666666666666667*sol_nodal_global[globalIdxs[3]]-0.1666666666666667*sol_nodal_global[globalIdxs[2]]-0.1666666666666667*sol_nodal_global[globalIdxs[1]]+0.1666666666666667*sol_nodal_global[globalIdxs[0]];

}
GKYL_CU_DH void fem_poisson_sol_stencil_2x_ser_p2(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = (-0.1666666666666667*sol_nodal_global[globalIdxs[7]])+0.6666666666666666*sol_nodal_global[globalIdxs[6]]-0.1666666666666667*sol_nodal_global[globalIdxs[5]]+0.6666666666666666*sol_nodal_global[globalIdxs[4]]+0.6666666666666666*sol_nodal_global[globalIdxs[3]]-0.1666666666666667*sol_nodal_global[globalIdxs[2]]+0.6666666666666666*sol_nodal_global[globalIdxs[1]]-0.1666666666666667*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = 0.09622504486493762*sol_nodal_global[globalIdxs[7]]-0.09622504486493762*sol_nodal_global[globalIdxs[5]]+0.3849001794597505*sol_nodal_global[globalIdxs[4]]-0.3849001794597505*sol_nodal_global[globalIdxs[3]]+0.09622504486493762*sol_nodal_global[globalIdxs[2]]-0.09622504486493762*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[2] = 0.09622504486493762*sol_nodal_global[globalIdxs[7]]+0.3849001794597505*sol_nodal_global[globalIdxs[6]]+0.09622504486493762*sol_nodal_global[globalIdxs[5]]-0.09622504486493762*sol_nodal_global[globalIdxs[2]]-0.3849001794597505*sol_nodal_global[globalIdxs[1]]-0.09622504486493762*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[3] = 0.1666666666666667*sol_nodal_global[globalIdxs[7]]-0.1666666666666667*sol_nodal_global[globalIdxs[5]]-0.1666666666666667*sol_nodal_global[globalIdxs[2]]+0.1666666666666667*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[4] = 0.149071198499986*sol_nodal_global[globalIdxs[7]]-0.298142396999972*sol_nodal_global[globalIdxs[6]]+0.149071198499986*sol_nodal_global[globalIdxs[5]]+0.149071198499986*sol_nodal_global[globalIdxs[2]]-0.298142396999972*sol_nodal_global[globalIdxs[1]]+0.149071198499986*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[5] = 0.149071198499986*sol_nodal_global[globalIdxs[7]]+0.149071198499986*sol_nodal_global[globalIdxs[5]]-0.298142396999972*sol_nodal_global[globalIdxs[4]]-0.298142396999972*sol_nodal_global[globalIdxs[3]]+0.149071198499986*sol_nodal_global[globalIdxs[2]]+0.149071198499986*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[6] = 0.08606629658238704*sol_nodal_global[globalIdxs[7]]-0.1721325931647741*sol_nodal_global[globalIdxs[6]]+0.08606629658238704*sol_nodal_global[globalIdxs[5]]-0.08606629658238704*sol_nodal_global[globalIdxs[2]]+0.1721325931647741*sol_nodal_global[globalIdxs[1]]-0.08606629658238704*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[7] = 0.08606629658238704*sol_nodal_global[globalIdxs[7]]-0.08606629658238704*sol_nodal_global[globalIdxs[5]]-0.1721325931647741*sol_nodal_global[globalIdxs[4]]+0.1721325931647741*sol_nodal_global[globalIdxs[3]]+0.08606629658238704*sol_nodal_global[globalIdxs[2]]-0.08606629658238704*sol_nodal_global[globalIdxs[0]];

}
