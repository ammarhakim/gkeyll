#include <gkyl_mom_canonical_pb_kernels.h> 
GKYL_CU_DH void canonical_pb_M1i_from_H_1x1v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double rdvx2 = 2.0/dxv[1];

  out[0] += (1.224744871391589*f[1]*hamil[3]*rdvx2+1.224744871391589*f[0]*hamil[2]*rdvx2)*volFact; 
  out[1] += (1.224744871391589*f[0]*hamil[3]*rdvx2+1.224744871391589*f[1]*hamil[2]*rdvx2)*volFact; 
} 
GKYL_CU_DH void canonical_pb_MEnergy_1x1v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 

  out[0] += (0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
} 
GKYL_CU_DH void canonical_pb_int_five_moments_1x1v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double dv1 = dxv[1]; 
  const double rdvx2 = 2.0/dxv[1];

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[1]*hamil[3]*rdvx2+1.7320508075688772*f[0]*hamil[2]*rdvx2)*volFact; 
  out[2] += (f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
