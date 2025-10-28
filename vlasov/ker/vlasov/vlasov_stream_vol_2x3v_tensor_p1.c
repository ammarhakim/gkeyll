#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_2x3v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w2Ddx0  = w[2]/dxv[0]; 
  double dv2Ddx0 = dxv[2]/dxv[0]; 
  double w3Ddx1  = w[3]/dxv[1]; 
  double dv3Ddx1 = dxv[3]/dxv[1]; 

  out[1] += 3.464101615137754*f[0]*w2Ddx0+f[3]*dv2Ddx0; 
  out[2] += 3.464101615137754*f[0]*w3Ddx1+f[4]*dv3Ddx1; 
  out[6] += 3.464101615137754*f[1]*w3Ddx1+3.464101615137754*f[2]*w2Ddx0+f[9]*dv3Ddx1+f[8]*dv2Ddx0; 
  out[7] += 3.464101615137754*f[3]*w2Ddx0+f[0]*dv2Ddx0; 
  out[8] += 3.464101615137754*f[3]*w3Ddx1+f[11]*dv3Ddx1; 
  out[9] += 3.464101615137754*f[4]*w2Ddx0+f[11]*dv2Ddx0; 
  out[10] += 3.464101615137754*f[4]*w3Ddx1+f[0]*dv3Ddx1; 
  out[12] += 3.464101615137754*f[5]*w2Ddx0+f[14]*dv2Ddx0; 
  out[13] += 3.464101615137754*f[5]*w3Ddx1+f[15]*dv3Ddx1; 
  out[16] += 3.464101615137754*f[7]*w3Ddx1+3.464101615137754*f[8]*w2Ddx0+f[18]*dv3Ddx1+f[2]*dv2Ddx0; 
  out[17] += 3.464101615137754*f[9]*w3Ddx1+3.464101615137754*f[10]*w2Ddx0+f[1]*dv3Ddx1+f[19]*dv2Ddx0; 
  out[18] += 3.464101615137754*f[11]*w2Ddx0+f[4]*dv2Ddx0; 
  out[19] += 3.464101615137754*f[11]*w3Ddx1+f[3]*dv3Ddx1; 
  out[20] += 3.464101615137754*f[12]*w3Ddx1+3.464101615137754*f[13]*w2Ddx0+f[23]*dv3Ddx1+f[22]*dv2Ddx0; 
  out[21] += 3.464101615137754*f[14]*w2Ddx0+f[5]*dv2Ddx0; 
  out[22] += 3.464101615137754*f[14]*w3Ddx1+f[25]*dv3Ddx1; 
  out[23] += 3.464101615137754*f[15]*w2Ddx0+f[25]*dv2Ddx0; 
  out[24] += 3.464101615137754*f[15]*w3Ddx1+f[5]*dv3Ddx1; 
  out[26] += 3.464101615137754*f[18]*w3Ddx1+3.464101615137754*f[19]*w2Ddx0+f[7]*dv3Ddx1+f[10]*dv2Ddx0; 
  out[27] += 3.464101615137754*f[21]*w3Ddx1+3.464101615137754*f[22]*w2Ddx0+f[29]*dv3Ddx1+f[13]*dv2Ddx0; 
  out[28] += 3.464101615137754*f[23]*w3Ddx1+3.464101615137754*f[24]*w2Ddx0+f[12]*dv3Ddx1+f[30]*dv2Ddx0; 
  out[29] += 3.464101615137754*f[25]*w2Ddx0+f[15]*dv2Ddx0; 
  out[30] += 3.464101615137754*f[25]*w3Ddx1+f[14]*dv3Ddx1; 
  out[31] += 3.464101615137754*f[29]*w3Ddx1+3.464101615137754*f[30]*w2Ddx0+f[21]*dv3Ddx1+f[24]*dv2Ddx0; 

  return 3.0*(fabs(w2Ddx0)+0.5*dv2Ddx0+fabs(w3Ddx1)+0.5*dv3Ddx1);
} 
