#include <gkyl_mom_canonical_pb_kernels.h> 
GKYL_CU_DH void canonical_pb_M1i_from_H_2x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double rdvx2 = 2.0/dxv[2];
  const double rdvy2 = 2.0/dxv[3];

  out[0] += (0.8660254037844386*f[12]*hamil[15]*rdvx2+0.8660254037844386*f[9]*hamil[14]*rdvx2+0.8660254037844386*f[8]*hamil[13]*rdvx2+0.8660254037844386*f[5]*hamil[11]*rdvx2+0.8660254037844386*f[4]*hamil[10]*rdvx2+0.8660254037844386*f[2]*hamil[7]*rdvx2+0.8660254037844386*f[1]*hamil[6]*rdvx2+0.8660254037844386*f[0]*hamil[3]*rdvx2)*volFact; 
  out[1] += (0.8660254037844386*f[9]*hamil[15]*rdvx2+0.8660254037844386*f[12]*hamil[14]*rdvx2+0.8660254037844386*f[4]*hamil[13]*rdvx2+0.8660254037844386*f[2]*hamil[11]*rdvx2+0.8660254037844386*f[8]*hamil[10]*rdvx2+0.8660254037844386*f[5]*hamil[7]*rdvx2+0.8660254037844386*f[0]*hamil[6]*rdvx2+0.8660254037844386*f[1]*hamil[3]*rdvx2)*volFact; 
  out[2] += (0.8660254037844386*f[8]*hamil[15]*rdvx2+0.8660254037844386*f[4]*hamil[14]*rdvx2+0.8660254037844386*f[12]*hamil[13]*rdvx2+0.8660254037844386*f[1]*hamil[11]*rdvx2+0.8660254037844386*f[9]*hamil[10]*rdvx2+0.8660254037844386*f[0]*hamil[7]*rdvx2+0.8660254037844386*f[5]*hamil[6]*rdvx2+0.8660254037844386*f[2]*hamil[3]*rdvx2)*volFact; 
  out[3] += (0.8660254037844386*f[4]*hamil[15]*rdvx2+0.8660254037844386*f[8]*hamil[14]*rdvx2+0.8660254037844386*f[9]*hamil[13]*rdvx2+0.8660254037844386*hamil[10]*f[12]*rdvx2+0.8660254037844386*f[0]*hamil[11]*rdvx2+0.8660254037844386*f[1]*hamil[7]*rdvx2+0.8660254037844386*f[2]*hamil[6]*rdvx2+0.8660254037844386*hamil[3]*f[5]*rdvx2)*volFact; 
  out[4] += (0.8660254037844386*f[11]*hamil[15]*rdvy2+0.8660254037844386*f[7]*hamil[14]*rdvy2+0.8660254037844386*f[6]*hamil[13]*rdvy2+0.8660254037844386*f[5]*hamil[12]*rdvy2+0.8660254037844386*f[3]*hamil[10]*rdvy2+0.8660254037844386*f[2]*hamil[9]*rdvy2+0.8660254037844386*f[1]*hamil[8]*rdvy2+0.8660254037844386*f[0]*hamil[4]*rdvy2)*volFact; 
  out[5] += (0.8660254037844386*f[7]*hamil[15]*rdvy2+0.8660254037844386*f[11]*hamil[14]*rdvy2+0.8660254037844386*f[3]*hamil[13]*rdvy2+0.8660254037844386*f[2]*hamil[12]*rdvy2+0.8660254037844386*f[6]*hamil[10]*rdvy2+0.8660254037844386*f[5]*hamil[9]*rdvy2+0.8660254037844386*f[0]*hamil[8]*rdvy2+0.8660254037844386*f[1]*hamil[4]*rdvy2)*volFact; 
  out[6] += (0.8660254037844386*f[6]*hamil[15]*rdvy2+0.8660254037844386*f[3]*hamil[14]*rdvy2+0.8660254037844386*f[11]*hamil[13]*rdvy2+0.8660254037844386*f[1]*hamil[12]*rdvy2+0.8660254037844386*f[7]*hamil[10]*rdvy2+0.8660254037844386*f[0]*hamil[9]*rdvy2+0.8660254037844386*f[5]*hamil[8]*rdvy2+0.8660254037844386*f[2]*hamil[4]*rdvy2)*volFact; 
  out[7] += (0.8660254037844386*f[3]*hamil[15]*rdvy2+0.8660254037844386*f[6]*hamil[14]*rdvy2+0.8660254037844386*f[7]*hamil[13]*rdvy2+0.8660254037844386*f[0]*hamil[12]*rdvy2+0.8660254037844386*hamil[10]*f[11]*rdvy2+0.8660254037844386*f[1]*hamil[9]*rdvy2+0.8660254037844386*f[2]*hamil[8]*rdvy2+0.8660254037844386*hamil[4]*f[5]*rdvy2)*volFact; 
} 
GKYL_CU_DH void canonical_pb_MEnergy_2x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 

  out[0] += (0.5*f[15]*hamil[15]+0.5*f[14]*hamil[14]+0.5*f[13]*hamil[13]+0.5*f[12]*hamil[12]+0.5*f[11]*hamil[11]+0.5*f[10]*hamil[10]+0.5*f[9]*hamil[9]+0.5*f[8]*hamil[8]+0.5*f[7]*hamil[7]+0.5*f[6]*hamil[6]+0.5*f[5]*hamil[5]+0.5*f[4]*hamil[4]+0.5*f[3]*hamil[3]+0.5*f[2]*hamil[2]+0.5*f[1]*hamil[1]+0.5*f[0]*hamil[0])*volFact; 
  out[1] += (0.5*f[14]*hamil[15]+0.5*hamil[14]*f[15]+0.5*f[10]*hamil[13]+0.5*hamil[10]*f[13]+0.5*f[9]*hamil[12]+0.5*hamil[9]*f[12]+0.5*f[7]*hamil[11]+0.5*hamil[7]*f[11]+0.5*f[4]*hamil[8]+0.5*hamil[4]*f[8]+0.5*f[3]*hamil[6]+0.5*hamil[3]*f[6]+0.5*f[2]*hamil[5]+0.5*hamil[2]*f[5]+0.5*f[0]*hamil[1]+0.5*hamil[0]*f[1])*volFact; 
  out[2] += (0.5*f[13]*hamil[15]+0.5*hamil[13]*f[15]+0.5*f[10]*hamil[14]+0.5*hamil[10]*f[14]+0.5*f[8]*hamil[12]+0.5*hamil[8]*f[12]+0.5*f[6]*hamil[11]+0.5*hamil[6]*f[11]+0.5*f[4]*hamil[9]+0.5*hamil[4]*f[9]+0.5*f[3]*hamil[7]+0.5*hamil[3]*f[7]+0.5*f[1]*hamil[5]+0.5*hamil[1]*f[5]+0.5*f[0]*hamil[2]+0.5*hamil[0]*f[2])*volFact; 
  out[3] += (0.5*f[10]*hamil[15]+0.5*hamil[10]*f[15]+0.5*f[13]*hamil[14]+0.5*hamil[13]*f[14]+0.5*f[4]*hamil[12]+0.5*hamil[4]*f[12]+0.5*f[3]*hamil[11]+0.5*hamil[3]*f[11]+0.5*f[8]*hamil[9]+0.5*hamil[8]*f[9]+0.5*f[6]*hamil[7]+0.5*hamil[6]*f[7]+0.5*f[0]*hamil[5]+0.5*hamil[0]*f[5]+0.5*f[1]*hamil[2]+0.5*hamil[1]*f[2])*volFact; 
} 
GKYL_CU_DH void canonical_pb_int_five_moments_2x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double dv1 = dxv[2]; 
  const double dv2 = dxv[3]; 
  const double rdvx2 = 2.0/dxv[2];
  const double rdvy2 = 2.0/dxv[3];

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.7320508075688772*f[12]*hamil[15]*rdvx2+1.7320508075688772*f[9]*hamil[14]*rdvx2+1.7320508075688772*f[8]*hamil[13]*rdvx2+1.7320508075688772*f[5]*hamil[11]*rdvx2+1.7320508075688772*f[4]*hamil[10]*rdvx2+1.7320508075688772*f[2]*hamil[7]*rdvx2+1.7320508075688772*f[1]*hamil[6]*rdvx2+1.7320508075688772*f[0]*hamil[3]*rdvx2)*volFact; 
  out[2] += (1.7320508075688772*f[11]*hamil[15]*rdvy2+1.7320508075688772*f[7]*hamil[14]*rdvy2+1.7320508075688772*f[6]*hamil[13]*rdvy2+1.7320508075688772*f[5]*hamil[12]*rdvy2+1.7320508075688772*f[3]*hamil[10]*rdvy2+1.7320508075688772*f[2]*hamil[9]*rdvy2+1.7320508075688772*f[1]*hamil[8]*rdvy2+1.7320508075688772*f[0]*hamil[4]*rdvy2)*volFact; 
  out[3] += (f[15]*hamil[15]+f[14]*hamil[14]+f[13]*hamil[13]+f[12]*hamil[12]+f[11]*hamil[11]+f[10]*hamil[10]+f[9]*hamil[9]+f[8]*hamil[8]+f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
