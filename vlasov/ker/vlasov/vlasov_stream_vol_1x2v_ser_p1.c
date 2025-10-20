#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w1Ddx0  = w[1]/dxv[0]; 
  double dv1Ddx0 = dxv[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w1Ddx0+f[2]*dv1Ddx0; 
  out[4] += 3.464101615137754*f[2]*w1Ddx0+(0.8944271909999159*f[8]+f[0])*dv1Ddx0; 
  out[5] += 3.464101615137754*f[3]*w1Ddx0+f[6]*dv1Ddx0; 
  out[7] += 3.464101615137754*f[6]*w1Ddx0+(0.8944271909999161*f[10]+f[3])*dv1Ddx0; 
  out[9] += 3.464101615137755*f[8]*w1Ddx0+0.8944271909999161*f[2]*dv1Ddx0; 
  out[11] += 3.464101615137755*f[10]*w1Ddx0+0.8944271909999159*f[6]*dv1Ddx0; 
  out[13] += 3.464101615137755*f[12]*w1Ddx0+f[14]*dv1Ddx0; 
  out[15] += 3.464101615137755*f[14]*w1Ddx0+f[12]*dv1Ddx0; 

  return 3.0*(fabs(w1Ddx0)+0.5*dv1Ddx0);
} 
