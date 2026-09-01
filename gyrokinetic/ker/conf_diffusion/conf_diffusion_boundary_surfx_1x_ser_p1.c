#include <gkyl_conf_diffusion_kernels.h>

GKYL_CU_DH double conf_diffusion_boundary_surfxx_1x_ser_p1(const double *wskin, const double *dxskin,
  const double *Kedge, const double *Kskin, int edge,
  const double *fedge, const double *fskin, double* GKYL_RESTRICT out)
{
  // Zero normal flux at the physical face; edge is the adjacent interior cell.
  // edge: -1 for the lower physical boundary, +1 for the upper boundary.

  const double rdx20 = 2.0/dxskin[0];
  const double rdx2sq0 = pow(2.0/dxskin[0],2.0);

  if (edge == -1) {
    out[0] += 0.015625*((21.213203435596427*Kedge[1]-21.213203435596427*Kskin[1]-12.24744871391589*Kedge[0]-12.24744871391589*Kskin[0])*fedge[1]+(21.213203435596427*Kedge[1]-21.213203435596427*Kskin[1]-12.24744871391589*Kedge[0]-12.24744871391589*Kskin[0])*fskin[1]+(22.0454076850486*fskin[0]-22.0454076850486*fedge[0])*Kedge[1]+(22.0454076850486*fedge[0]-22.0454076850486*fskin[0])*Kskin[1]+(12.727922061357857*Kedge[0]+12.727922061357857*Kskin[0])*fedge[0]+(-(12.727922061357857*Kedge[0])-12.727922061357857*Kskin[0])*fskin[0])*rdx2sq0;
    out[1] += 0.015625*((36.74234614174767*Kedge[1]-36.74234614174767*Kskin[1]-21.213203435596427*Kedge[0]-21.213203435596427*Kskin[0])*fedge[1]+(36.74234614174767*Kedge[1]-36.74234614174767*Kskin[1]-21.213203435596427*Kedge[0]-21.213203435596427*Kskin[0])*fskin[1]+(38.18376618407357*fskin[0]-38.18376618407357*fedge[0])*Kedge[1]+(38.18376618407357*fedge[0]-38.18376618407357*fskin[0])*Kskin[1]+(22.0454076850486*Kedge[0]+22.0454076850486*Kskin[0])*fedge[0]+(-(22.0454076850486*Kedge[0])-22.0454076850486*Kskin[0])*fskin[0])*rdx2sq0;
  }
  else {
    out[0] += 0.015625*((21.213203435596427*Kedge[1]-21.213203435596427*Kskin[1]+12.24744871391589*Kedge[0]+12.24744871391589*Kskin[0])*fedge[1]+(21.213203435596427*Kedge[1]-21.213203435596427*Kskin[1]+12.24744871391589*Kedge[0]+12.24744871391589*Kskin[0])*fskin[1]+(22.0454076850486*fedge[0]-22.0454076850486*fskin[0])*Kedge[1]+(22.0454076850486*fskin[0]-22.0454076850486*fedge[0])*Kskin[1]+(12.727922061357857*Kedge[0]+12.727922061357857*Kskin[0])*fedge[0]+(-(12.727922061357857*Kedge[0])-12.727922061357857*Kskin[0])*fskin[0])*rdx2sq0;
    out[1] += -0.015625*((36.74234614174767*Kedge[1]-36.74234614174767*Kskin[1]+21.213203435596427*Kedge[0]+21.213203435596427*Kskin[0])*fedge[1]+(36.74234614174767*Kedge[1]-36.74234614174767*Kskin[1]+21.213203435596427*Kedge[0]+21.213203435596427*Kskin[0])*fskin[1]+(38.18376618407357*fedge[0]-38.18376618407357*fskin[0])*Kedge[1]+(38.18376618407357*fskin[0]-38.18376618407357*fedge[0])*Kskin[1]+(22.0454076850486*Kedge[0]+22.0454076850486*Kskin[0])*fedge[0]+(-(22.0454076850486*Kedge[0])-22.0454076850486*Kskin[0])*fskin[0])*rdx2sq0;
  }

  return 0.;
}
