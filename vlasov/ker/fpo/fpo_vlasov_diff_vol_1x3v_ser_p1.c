#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH double
fpo_vlasov_diff_vol_1x3v_ser_p1(const double* w, const double* dx,
  const double* g, const double* f, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // g: 
  // f: 
  // out: Incremented output

  const double Jvxvx = 4/dx[1]/dx[1];

  out[0] +=  Jvxvx*(0.0);
  out[1] +=  Jvxvx*(0.0);
  out[2] +=  Jvxvx*(0.0);
  out[3] +=  Jvxvx*(0.0);
  out[4] +=  Jvxvx*(0.0);
  out[5] +=  Jvxvx*(0.0);
  out[6] +=  Jvxvx*(0.0);
  out[7] +=  Jvxvx*(0.0);
  out[8] +=  Jvxvx*(0.0);
  out[9] +=  Jvxvx*(0.0);
  out[10] +=  Jvxvx*(0.0);
  out[11] +=  Jvxvx*(0.0);
  out[12] +=  Jvxvx*(0.0);
  out[13] +=  Jvxvx*(0.0);
  out[14] +=  Jvxvx*(0.0);
  out[15] +=  Jvxvx*(0.0);

  const double Jvxvy = 4/dx[1]/dx[2];

  out[0] +=  Jvxvy*(0.0);
  out[1] +=  Jvxvy*(0.0);
  out[2] +=  Jvxvy*(0.0);
  out[3] +=  Jvxvy*(0.0);
  out[4] +=  Jvxvy*(0.0);
  out[5] +=  Jvxvy*(0.0);
  out[6] +=  Jvxvy*(0.0);
  out[7] +=  Jvxvy*(2.25*f[8]*g[15]+2.25*f[4]*g[14]+2.25*f[1]*g[11]+2.25*f[0]*g[7]);
  out[8] +=  Jvxvy*(0.0);
  out[9] +=  Jvxvy*(0.0);
  out[10] +=  Jvxvy*(0.0);
  out[11] +=  Jvxvy*(2.25*f[4]*g[15]+2.25*f[8]*g[14]+2.25*f[0]*g[11]+2.25*f[1]*g[7]);
  out[12] +=  Jvxvy*(0.0);
  out[13] +=  Jvxvy*(0.0);
  out[14] +=  Jvxvy*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[8]*g[11]+2.25*f[4]*g[7]);
  out[15] +=  Jvxvy*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[4]*g[11]+2.25*g[7]*f[8]);

  const double Jvxvz = 4/dx[1]/dx[3];

  out[0] +=  Jvxvz*(0.0);
  out[1] +=  Jvxvz*(0.0);
  out[2] +=  Jvxvz*(0.0);
  out[3] +=  Jvxvz*(0.0);
  out[4] +=  Jvxvz*(0.0);
  out[5] +=  Jvxvz*(0.0);
  out[6] +=  Jvxvz*(0.0);
  out[7] +=  Jvxvz*(0.0);
  out[8] +=  Jvxvz*(0.0);
  out[9] +=  Jvxvz*(2.25*f[6]*g[15]+2.25*f[3]*g[14]+2.25*f[1]*g[12]+2.25*f[0]*g[9]);
  out[10] +=  Jvxvz*(0.0);
  out[11] +=  Jvxvz*(0.0);
  out[12] +=  Jvxvz*(2.25*f[3]*g[15]+2.25*f[6]*g[14]+2.25*f[0]*g[12]+2.25*f[1]*g[9]);
  out[13] +=  Jvxvz*(0.0);
  out[14] +=  Jvxvz*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[6]*g[12]+2.25*f[3]*g[9]);
  out[15] +=  Jvxvz*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[3]*g[12]+2.25*f[6]*g[9]);

  const double Jvyvx = 4/dx[2]/dx[1];

  out[0] +=  Jvyvx*(0.0);
  out[1] +=  Jvyvx*(0.0);
  out[2] +=  Jvyvx*(0.0);
  out[3] +=  Jvyvx*(0.0);
  out[4] +=  Jvyvx*(0.0);
  out[5] +=  Jvyvx*(0.0);
  out[6] +=  Jvyvx*(0.0);
  out[7] +=  Jvyvx*(2.25*f[8]*g[15]+2.25*f[4]*g[14]+2.25*f[1]*g[11]+2.25*f[0]*g[7]);
  out[8] +=  Jvyvx*(0.0);
  out[9] +=  Jvyvx*(0.0);
  out[10] +=  Jvyvx*(0.0);
  out[11] +=  Jvyvx*(2.25*f[4]*g[15]+2.25*f[8]*g[14]+2.25*f[0]*g[11]+2.25*f[1]*g[7]);
  out[12] +=  Jvyvx*(0.0);
  out[13] +=  Jvyvx*(0.0);
  out[14] +=  Jvyvx*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[8]*g[11]+2.25*f[4]*g[7]);
  out[15] +=  Jvyvx*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[4]*g[11]+2.25*g[7]*f[8]);

  const double Jvyvy = 4/dx[2]/dx[2];

  out[0] +=  Jvyvy*(0.0);
  out[1] +=  Jvyvy*(0.0);
  out[2] +=  Jvyvy*(0.0);
  out[3] +=  Jvyvy*(0.0);
  out[4] +=  Jvyvy*(0.0);
  out[5] +=  Jvyvy*(0.0);
  out[6] +=  Jvyvy*(0.0);
  out[7] +=  Jvyvy*(0.0);
  out[8] +=  Jvyvy*(0.0);
  out[9] +=  Jvyvy*(0.0);
  out[10] +=  Jvyvy*(0.0);
  out[11] +=  Jvyvy*(0.0);
  out[12] +=  Jvyvy*(0.0);
  out[13] +=  Jvyvy*(0.0);
  out[14] +=  Jvyvy*(0.0);
  out[15] +=  Jvyvy*(0.0);

  const double Jvyvz = 4/dx[2]/dx[3];

  out[0] +=  Jvyvz*(0.0);
  out[1] +=  Jvyvz*(0.0);
  out[2] +=  Jvyvz*(0.0);
  out[3] +=  Jvyvz*(0.0);
  out[4] +=  Jvyvz*(0.0);
  out[5] +=  Jvyvz*(0.0);
  out[6] +=  Jvyvz*(0.0);
  out[7] +=  Jvyvz*(0.0);
  out[8] +=  Jvyvz*(0.0);
  out[9] +=  Jvyvz*(0.0);
  out[10] +=  Jvyvz*(2.25*f[5]*g[15]+2.25*f[2]*g[14]+2.25*f[1]*g[13]+2.25*f[0]*g[10]);
  out[11] +=  Jvyvz*(0.0);
  out[12] +=  Jvyvz*(0.0);
  out[13] +=  Jvyvz*(2.25*f[2]*g[15]+2.25*f[5]*g[14]+2.25*f[0]*g[13]+2.25*f[1]*g[10]);
  out[14] +=  Jvyvz*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[5]*g[13]+2.25*f[2]*g[10]);
  out[15] +=  Jvyvz*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[2]*g[13]+2.25*f[5]*g[10]);

  const double Jvzvx = 4/dx[3]/dx[1];

  out[0] +=  Jvzvx*(0.0);
  out[1] +=  Jvzvx*(0.0);
  out[2] +=  Jvzvx*(0.0);
  out[3] +=  Jvzvx*(0.0);
  out[4] +=  Jvzvx*(0.0);
  out[5] +=  Jvzvx*(0.0);
  out[6] +=  Jvzvx*(0.0);
  out[7] +=  Jvzvx*(0.0);
  out[8] +=  Jvzvx*(0.0);
  out[9] +=  Jvzvx*(2.25*f[6]*g[15]+2.25*f[3]*g[14]+2.25*f[1]*g[12]+2.25*f[0]*g[9]);
  out[10] +=  Jvzvx*(0.0);
  out[11] +=  Jvzvx*(0.0);
  out[12] +=  Jvzvx*(2.25*f[3]*g[15]+2.25*f[6]*g[14]+2.25*f[0]*g[12]+2.25*f[1]*g[9]);
  out[13] +=  Jvzvx*(0.0);
  out[14] +=  Jvzvx*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[6]*g[12]+2.25*f[3]*g[9]);
  out[15] +=  Jvzvx*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[3]*g[12]+2.25*f[6]*g[9]);

  const double Jvzvy = 4/dx[3]/dx[2];

  out[0] +=  Jvzvy*(0.0);
  out[1] +=  Jvzvy*(0.0);
  out[2] +=  Jvzvy*(0.0);
  out[3] +=  Jvzvy*(0.0);
  out[4] +=  Jvzvy*(0.0);
  out[5] +=  Jvzvy*(0.0);
  out[6] +=  Jvzvy*(0.0);
  out[7] +=  Jvzvy*(0.0);
  out[8] +=  Jvzvy*(0.0);
  out[9] +=  Jvzvy*(0.0);
  out[10] +=  Jvzvy*(2.25*f[5]*g[15]+2.25*f[2]*g[14]+2.25*f[1]*g[13]+2.25*f[0]*g[10]);
  out[11] +=  Jvzvy*(0.0);
  out[12] +=  Jvzvy*(0.0);
  out[13] +=  Jvzvy*(2.25*f[2]*g[15]+2.25*f[5]*g[14]+2.25*f[0]*g[13]+2.25*f[1]*g[10]);
  out[14] +=  Jvzvy*(2.25*f[1]*g[15]+2.25*f[0]*g[14]+2.25*f[5]*g[13]+2.25*f[2]*g[10]);
  out[15] +=  Jvzvy*(2.25*f[0]*g[15]+2.25*f[1]*g[14]+2.25*f[2]*g[13]+2.25*f[5]*g[10]);

  const double Jvzvz = 4/dx[3]/dx[3];

  out[0] +=  Jvzvz*(0.0);
  out[1] +=  Jvzvz*(0.0);
  out[2] +=  Jvzvz*(0.0);
  out[3] +=  Jvzvz*(0.0);
  out[4] +=  Jvzvz*(0.0);
  out[5] +=  Jvzvz*(0.0);
  out[6] +=  Jvzvz*(0.0);
  out[7] +=  Jvzvz*(0.0);
  out[8] +=  Jvzvz*(0.0);
  out[9] +=  Jvzvz*(0.0);
  out[10] +=  Jvzvz*(0.0);
  out[11] +=  Jvzvz*(0.0);
  out[12] +=  Jvzvz*(0.0);
  out[13] +=  Jvzvz*(0.0);
  out[14] +=  Jvzvz*(0.0);
  out[15] +=  Jvzvz*(0.0);

  return Jvxvx*(0.0) + Jvxvy*(0.75*g[7]) + Jvxvz*(0.75*g[9]) + Jvyvx*(0.75*g[7]) + Jvyvy*(0.0) + Jvyvz*(0.75*g[10]) + Jvzvx*(0.75*g[9]) + Jvzvy*(0.75*g[10]) + Jvzvz*(0.0);
}
