#include <gkyl_mom_canonical_pb_kernels.h> 
GKYL_CU_DH void canonical_pb_M1i_from_H_1x2v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double rdvx2 = 2.0/dxv[1];
  const double rdvy2 = 2.0/dxv[2];

  out[0] += (1.224744871391589*f[13]*hamil[15]*rdvx2+1.224744871391589*f[12]*hamil[14]*rdvx2+2.7386127875258306*f[7]*hamil[11]*rdvx2+2.7386127875258306*f[6]*hamil[10]*rdvx2+2.7386127875258306*f[4]*hamil[9]*rdvx2+2.7386127875258306*f[2]*hamil[8]*rdvx2+1.224744871391589*f[5]*hamil[7]*rdvx2+1.224744871391589*f[3]*hamil[6]*rdvx2+1.224744871391589*f[1]*hamil[4]*rdvx2+1.224744871391589*f[0]*hamil[2]*rdvx2)*volFact; 
  out[1] += (1.224744871391589*f[12]*hamil[15]*rdvx2+1.224744871391589*f[13]*hamil[14]*rdvx2+2.7386127875258306*f[6]*hamil[11]*rdvx2+2.7386127875258306*f[7]*hamil[10]*rdvx2+2.7386127875258306*f[2]*hamil[9]*rdvx2+2.7386127875258306*f[4]*hamil[8]*rdvx2+1.224744871391589*f[3]*hamil[7]*rdvx2+1.224744871391589*f[5]*hamil[6]*rdvx2+1.224744871391589*f[0]*hamil[4]*rdvx2+1.224744871391589*f[1]*hamil[2]*rdvx2)*volFact; 
  out[2] += (2.7386127875258306*f[7]*hamil[15]*rdvy2+2.7386127875258306*f[6]*hamil[14]*rdvy2+2.7386127875258306*f[5]*hamil[13]*rdvy2+2.7386127875258306*f[3]*hamil[12]*rdvy2+1.224744871391589*f[9]*hamil[11]*rdvy2+1.224744871391589*f[8]*hamil[10]*rdvy2+1.224744871391589*f[4]*hamil[7]*rdvy2+1.224744871391589*f[2]*hamil[6]*rdvy2+1.224744871391589*f[1]*hamil[5]*rdvy2+1.224744871391589*f[0]*hamil[3]*rdvy2)*volFact; 
  out[3] += (2.7386127875258306*f[6]*hamil[15]*rdvy2+2.7386127875258306*f[7]*hamil[14]*rdvy2+2.7386127875258306*f[3]*hamil[13]*rdvy2+2.7386127875258306*f[5]*hamil[12]*rdvy2+1.224744871391589*f[8]*hamil[11]*rdvy2+1.224744871391589*f[9]*hamil[10]*rdvy2+1.224744871391589*f[2]*hamil[7]*rdvy2+1.224744871391589*f[4]*hamil[6]*rdvy2+1.224744871391589*f[0]*hamil[5]*rdvy2+1.224744871391589*f[1]*hamil[3]*rdvy2)*volFact; 
} 
GKYL_CU_DH void canonical_pb_MEnergy_1x2v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 

  out[0] += (0.7071067811865475*f[15]*hamil[15]+0.7071067811865475*f[14]*hamil[14]+0.7071067811865475*f[13]*hamil[13]+0.7071067811865475*f[12]*hamil[12]+0.7071067811865475*f[11]*hamil[11]+0.7071067811865475*f[10]*hamil[10]+0.7071067811865475*f[9]*hamil[9]+0.7071067811865475*f[8]*hamil[8]+0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[14]*hamil[15]+0.7071067811865475*hamil[14]*f[15]+0.7071067811865475*f[12]*hamil[13]+0.7071067811865475*hamil[12]*f[13]+0.7071067811865475*f[10]*hamil[11]+0.7071067811865475*hamil[10]*f[11]+0.7071067811865475*f[8]*hamil[9]+0.7071067811865475*hamil[8]*f[9]+0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void canonical_pb_int_five_moments_1x2v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv1 = dxv[1]; 
  const double dv2 = dxv[2]; 
  const double rdvx2 = 2.0/dxv[1];
  const double rdvy2 = 2.0/dxv[2];

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688774*f[13]*hamil[15]*rdvx2+1.7320508075688774*f[12]*hamil[14]*rdvx2+3.872983346207417*f[7]*hamil[11]*rdvx2+3.872983346207417*f[6]*hamil[10]*rdvx2+3.872983346207417*f[4]*hamil[9]*rdvx2+3.872983346207417*f[2]*hamil[8]*rdvx2+1.7320508075688772*f[5]*hamil[7]*rdvx2+1.7320508075688772*f[3]*hamil[6]*rdvx2+1.7320508075688772*f[1]*hamil[4]*rdvx2+1.7320508075688772*f[0]*hamil[2]*rdvx2)*volFact; 
  out[2] += (3.872983346207417*f[7]*hamil[15]*rdvy2+3.872983346207417*f[6]*hamil[14]*rdvy2+3.872983346207417*f[5]*hamil[13]*rdvy2+3.872983346207417*f[3]*hamil[12]*rdvy2+1.7320508075688774*f[9]*hamil[11]*rdvy2+1.7320508075688774*f[8]*hamil[10]*rdvy2+1.7320508075688772*f[4]*hamil[7]*rdvy2+1.7320508075688772*f[2]*hamil[6]*rdvy2+1.7320508075688772*f[1]*hamil[5]*rdvy2+1.7320508075688772*f[0]*hamil[3]*rdvy2)*volFact; 
  out[3] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
