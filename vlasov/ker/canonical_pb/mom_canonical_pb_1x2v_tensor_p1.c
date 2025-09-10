#include <gkyl_mom_canonical_pb_kernels.h> 
GKYL_CU_DH void canonical_pb_M1i_from_H_1x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double rdvx2 = 2.0/dxv[1];
  const double rdvy2 = 2.0/dxv[2];

  out[0] += (1.224744871391589*f[5]*hamil[7]*rdvx2+1.224744871391589*f[3]*hamil[6]*rdvx2+1.224744871391589*f[1]*hamil[4]*rdvx2+1.224744871391589*f[0]*hamil[2]*rdvx2)*volFact; 
  out[1] += (1.224744871391589*f[3]*hamil[7]*rdvx2+1.224744871391589*f[5]*hamil[6]*rdvx2+1.224744871391589*f[0]*hamil[4]*rdvx2+1.224744871391589*f[1]*hamil[2]*rdvx2)*volFact; 
  out[2] += (1.224744871391589*f[4]*hamil[7]*rdvy2+1.224744871391589*f[2]*hamil[6]*rdvy2+1.224744871391589*f[1]*hamil[5]*rdvy2+1.224744871391589*f[0]*hamil[3]*rdvy2)*volFact; 
  out[3] += (1.224744871391589*f[2]*hamil[7]*rdvy2+1.224744871391589*f[4]*hamil[6]*rdvy2+1.224744871391589*f[0]*hamil[5]*rdvy2+1.224744871391589*f[1]*hamil[3]*rdvy2)*volFact; 
} 
GKYL_CU_DH void canonical_pb_MEnergy_1x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 

  out[0] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[6]*hamil[7]+0.7071067811865475*hamil[6]*f[7]+0.7071067811865475*f[3]*hamil[5]+0.7071067811865475*hamil[3]*f[5]+0.7071067811865475*f[2]*hamil[4]+0.7071067811865475*hamil[2]*f[4]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void canonical_pb_int_five_moments_1x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double dv1 = dxv[1]; 
  const double dv2 = dxv[2]; 
  const double rdvx2 = 2.0/dxv[1];
  const double rdvy2 = 2.0/dxv[2];

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[5]*hamil[7]*rdvx2+1.7320508075688772*f[3]*hamil[6]*rdvx2+1.7320508075688772*f[1]*hamil[4]*rdvx2+1.7320508075688772*f[0]*hamil[2]*rdvx2)*volFact; 
  out[2] += (1.7320508075688772*f[4]*hamil[7]*rdvy2+1.7320508075688772*f[2]*hamil[6]*rdvy2+1.7320508075688772*f[1]*hamil[5]*rdvy2+1.7320508075688772*f[0]*hamil[3]*rdvy2)*volFact; 
  out[3] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
