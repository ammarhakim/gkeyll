#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_compton_1v_ser_p3(const double *w, const double *dxv, const double *vmap, 
  double t_cool, double p0, double* GKYL_RESTRICT rad) 
{ 
  // w:   Cell-center coordinates of velocity grid.
  // dxv: Cell spacing of velocity grid.
  // vmap: Velocity-space map for nonuniform meshes.
  // t_cool: Cooling time.
  // p0: Momentum to relax to.
  // rad: Radiation field (-(p - p0)/tcool for Compton, -sign(p)*p^4 = -|p|*p^3 for curvature).
 
  double wx1 = w[0], dv1 = dxv[0]; 
  double *rad_vx = &rad[0];
  double rad_nodal_vx[4] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(wx1-0.5*dv1 - p0);

  rad_nodal_vx[1] = -1.0/t_cool*(wx1-0.22360679774997896*dv1 - p0);

  rad_nodal_vx[2] = -1.0/t_cool*(wx1+0.22360679774997896*dv1 - p0);

  rad_nodal_vx[3] = -1.0/t_cool*(wx1+0.5*dv1 - p0);

  rad_vx[0] = 0.11785113019775789*rad_nodal_vx[3]+0.5892556509887895*rad_nodal_vx[2]+0.5892556509887895*rad_nodal_vx[1]+0.11785113019775789*rad_nodal_vx[0]; 
  rad_vx[1] = 0.20412414523193148*rad_nodal_vx[3]+0.4564354645876384*rad_nodal_vx[2]-0.4564354645876384*rad_nodal_vx[1]-0.20412414523193148*rad_nodal_vx[0]; 
  rad_vx[2] = 0.2635231383473649*rad_nodal_vx[3]-0.2635231383473649*rad_nodal_vx[2]-0.2635231383473649*rad_nodal_vx[1]+0.2635231383473649*rad_nodal_vx[0]; 
  rad_vx[3] = 0.13363062095621214*rad_nodal_vx[3]-0.29880715233359834*rad_nodal_vx[2]+0.29880715233359834*rad_nodal_vx[1]-0.13363062095621214*rad_nodal_vx[0]; 

} 
