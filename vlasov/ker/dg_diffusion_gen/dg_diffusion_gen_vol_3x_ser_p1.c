#include <gkyl_dg_diffusion_gen_kernels.h>

GKYL_CU_DH double
dg_diffusion_gen_vol_3x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* qIn, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jxx = 4/dx[0]/dx[0];
  const double Jxy = 4/dx[0]/dx[1];
  const double Jxz = 4/dx[0]/dx[2];
  const double Jyy = 4/dx[1]/dx[1];
  const double Jyz = 4/dx[1]/dx[2];
  const double Jzz = 4/dx[2]/dx[2];

  const double* Dxx = &Dij[0];
  const double* Dxy = &Dij[8];
  const double* Dxz = &Dij[16];
  const double* Dyy = &Dij[24];
  const double* Dyz = &Dij[32];
  const double* Dzz = &Dij[40];

  return Jxx*Dxx[0] + Jxy*Dxy[0] + Jxz*Dxz[0] + Jyy*Dyy[0] + Jyz*Dyz[0] + Jzz*Dzz[0];
}
