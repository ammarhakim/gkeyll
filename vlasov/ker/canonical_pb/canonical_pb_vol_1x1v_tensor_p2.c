#include <gkyl_canonical_pb_kernels.h> 
double canonical_pb_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *hamil, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 0.4*(33.54101966249684*f[6]*hamil[8]+33.54101966249684*f[3]*hamil[7]+15.0*f[4]*hamil[6]+33.54101966249685*f[2]*hamil[5]+15.0*f[1]*hamil[3]+15.0*f[0]*hamil[2])*dxdvInv0; 
  out[2] += -0.4*(33.54101966249684*f[7]*hamil[8]+15.0*f[5]*hamil[7]+33.54101966249684*f[3]*hamil[6]+33.54101966249685*f[1]*hamil[4]+15.0*f[2]*hamil[3]+15.0*f[0]*hamil[1])*dxdvInv0; 
  out[3] += -2.0*(6.708203932499369*f[5]*hamil[8]-6.708203932499369*f[4]*hamil[8]-3.0*f[7]*hamil[7]-6.708203932499369*f[1]*hamil[7]+3.0*f[6]*hamil[6]+6.708203932499369*f[2]*hamil[6]-6.0*f[5]*hamil[5]-6.708203932499369*f[0]*hamil[5]+6.0*f[4]*hamil[4]+6.708203932499369*f[0]*hamil[4]-3.0*f[2]*hamil[2]+3.0*f[1]*hamil[1])*dxdvInv0; 
  out[4] += 0.4*(67.0820393249937*f[3]*hamil[8]+67.0820393249937*f[6]*hamil[7]+75.00000000000001*f[2]*hamil[7]+30.0*f[1]*hamil[6]+75.0*f[3]*hamil[5]+30.0*hamil[3]*f[4]+33.54101966249685*f[0]*hamil[3]+33.54101966249685*f[1]*hamil[2])*dxdvInv0; 
  out[5] += -0.4*(67.0820393249937*f[3]*hamil[8]+30.0*f[2]*hamil[7]+67.0820393249937*hamil[6]*f[7]+75.00000000000001*f[1]*hamil[6]+30.0*hamil[3]*f[5]+75.0*f[3]*hamil[4]+33.54101966249685*f[0]*hamil[3]+33.54101966249685*hamil[1]*f[2])*dxdvInv0; 
  out[6] += 0.4*(30.0*f[7]*hamil[8]+67.08203932499369*f[1]*hamil[8]+45.0*hamil[7]*f[8]+67.0820393249937*f[5]*hamil[7]+67.0820393249937*f[4]*hamil[7]+75.0*f[0]*hamil[7]+67.0820393249937*hamil[5]*f[7]+15.0*hamil[3]*f[6]+75.00000000000001*f[1]*hamil[5]-30.0*f[1]*hamil[4]-15.0*hamil[1]*f[4]+33.54101966249684*f[2]*hamil[3]+33.54101966249684*hamil[2]*f[3])*dxdvInv0; 
  out[7] += -0.4*(30.0*f[6]*hamil[8]+67.08203932499369*f[2]*hamil[8]+45.0*hamil[6]*f[8]+15.0*hamil[3]*f[7]+67.0820393249937*f[5]*hamil[6]+67.0820393249937*f[4]*hamil[6]+75.0*f[0]*hamil[6]+67.0820393249937*hamil[4]*f[6]-30.0*f[2]*hamil[5]-15.0*hamil[2]*f[5]+75.00000000000001*f[2]*hamil[4]+33.54101966249684*f[1]*hamil[3]+33.54101966249684*hamil[1]*f[3])*dxdvInv0; 
  out[8] += 2.0*(6.0*f[6]*hamil[7]+13.41640786499874*f[2]*hamil[7]-6.0*hamil[6]*f[7]+6.708203932499369*hamil[2]*f[7]-13.41640786499874*f[1]*hamil[6]-6.708203932499369*hamil[1]*f[6]+13.41640786499874*f[3]*hamil[5]+6.708203932499369*hamil[3]*f[5]-13.41640786499874*f[3]*hamil[4]-6.708203932499369*hamil[3]*f[4])*dxdvInv0; 
  return 0.; 
} 
