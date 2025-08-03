#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void rad_curvature_2v_ser_p1(const double *w, const double *dxv, const double *vmap, 
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
  double wx2 = w[1], dv2 = dxv[1]; 
  double *rad_vy = &rad[4];
  double rad_nodal_vy[4] = {0.0};

  rad_nodal_vx[0] = -1.0/t_cool*(pow(wx1-0.5*dv1, 3.0)*fabs(wx1-0.5*dv1));
  rad_nodal_vy[0] = -1.0/t_cool*(pow(wx2-0.5*dv2, 3.0)*fabs(wx2-0.5*dv2));

  rad_nodal_vx[1] = -1.0/t_cool*(pow(wx1+0.5*dv1, 3.0)*fabs(wx1+0.5*dv1));
  rad_nodal_vy[1] = -1.0/t_cool*(pow(wx2-0.5*dv2, 3.0)*fabs(wx2-0.5*dv2));

  rad_nodal_vx[2] = -1.0/t_cool*(pow(wx1-0.5*dv1, 3.0)*fabs(wx1-0.5*dv1));
  rad_nodal_vy[2] = -1.0/t_cool*(pow(wx2+0.5*dv2, 3.0)*fabs(wx2+0.5*dv2));

  rad_nodal_vx[3] = -1.0/t_cool*(pow(wx1+0.5*dv1, 3.0)*fabs(wx1+0.5*dv1));
  rad_nodal_vy[3] = -1.0/t_cool*(pow(wx2+0.5*dv2, 3.0)*fabs(wx2+0.5*dv2));

  rad_vx[0] = 0.5*rad_nodal_vx[3]+0.5*rad_nodal_vx[2]+0.5*rad_nodal_vx[1]+0.5*rad_nodal_vx[0]; 
  rad_vx[1] = 0.2886751345948129*rad_nodal_vx[3]-0.2886751345948129*rad_nodal_vx[2]+0.2886751345948129*rad_nodal_vx[1]-0.2886751345948129*rad_nodal_vx[0]; 
  rad_vx[2] = 0.2886751345948129*rad_nodal_vx[3]+0.2886751345948129*rad_nodal_vx[2]-0.2886751345948129*rad_nodal_vx[1]-0.2886751345948129*rad_nodal_vx[0]; 
  rad_vx[3] = 0.16666666666666666*rad_nodal_vx[3]-0.16666666666666666*rad_nodal_vx[2]-0.16666666666666666*rad_nodal_vx[1]+0.16666666666666666*rad_nodal_vx[0]; 

  rad_vy[0] = 0.5*rad_nodal_vy[3]+0.5*rad_nodal_vy[2]+0.5*rad_nodal_vy[1]+0.5*rad_nodal_vy[0]; 
  rad_vy[1] = 0.2886751345948129*rad_nodal_vy[3]-0.2886751345948129*rad_nodal_vy[2]+0.2886751345948129*rad_nodal_vy[1]-0.2886751345948129*rad_nodal_vy[0]; 
  rad_vy[2] = 0.2886751345948129*rad_nodal_vy[3]+0.2886751345948129*rad_nodal_vy[2]-0.2886751345948129*rad_nodal_vy[1]-0.2886751345948129*rad_nodal_vy[0]; 
  rad_vy[3] = 0.16666666666666666*rad_nodal_vy[3]-0.16666666666666666*rad_nodal_vy[2]-0.16666666666666666*rad_nodal_vy[1]+0.16666666666666666*rad_nodal_vy[0]; 

} 
