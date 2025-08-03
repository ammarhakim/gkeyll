#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_curvature_1v_ser_p1(const double *w, const double *dxv, const double *vmap, 
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
  double rad_nodal_vx[2] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(pow(wx1-0.5*dv1, 3.0)*fabs(wx1-0.5*dv1));

  rad_nodal_vx[1] = -1.0/t_cool*(pow(wx1+0.5*dv1, 3.0)*fabs(wx1+0.5*dv1));

  rad_vx[0] = 0.7071067811865475*rad_nodal_vx[1]+0.7071067811865475*rad_nodal_vx[0]; 
  rad_vx[1] = 0.408248290463863*rad_nodal_vx[1]-0.408248290463863*rad_nodal_vx[0]; 

} 
