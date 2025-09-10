#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 

  double alphaDrag[32]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.0*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.0*nuUSum[1]-2.0*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -0.5773502691896258*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[4] = -0.5773502691896258*dxv[1]*nuSum[1]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[16] = (2.0*nuUSum[2]-2.0*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[17] = (2.0*nuUSum[3]-2.0*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[19] = -0.5773502691896258*nuSum[0]*dxv[2]*rdvy2; 
  alphaDrag[21] = -0.5773502691896258*nuSum[1]*dxv[2]*rdvy2; 

  out[2] += 0.6123724356957944*(alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[5]*alphaDrag[21]+f[3]*alphaDrag[19]+f[1]*alphaDrag[17]+f[0]*alphaDrag[16]); 
  out[4] += 0.6123724356957944*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[3]*alphaDrag[21]+f[5]*alphaDrag[19]+f[0]*alphaDrag[17]+f[1]*alphaDrag[16]); 
  out[6] += 0.6123724356957944*(f[7]*alphaDrag[21]+f[6]*alphaDrag[19]+f[4]*alphaDrag[17]+f[2]*alphaDrag[16]+alphaDrag[4]*f[7]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[6]*alphaDrag[21]+f[7]*alphaDrag[19]+f[2]*alphaDrag[17]+f[4]*alphaDrag[16]+alphaDrag[2]*f[7]+alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[8] += 1.224744871391589*(alphaDrag[4]*f[9]+alphaDrag[2]*f[8])+1.369306393762915*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[2]*f[9]+alphaDrag[4]*f[8])+1.369306393762915*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += 0.6123724356957944*(f[11]*alphaDrag[21]+f[10]*alphaDrag[19]+f[9]*alphaDrag[17]+f[8]*alphaDrag[16])+1.224744871391589*(alphaDrag[4]*f[11]+alphaDrag[2]*f[10])+1.369306393762915*(alphaDrag[1]*f[7]+alphaDrag[0]*f[6]+alphaDrag[4]*f[5]+alphaDrag[2]*f[3]); 
  out[11] += 0.6123724356957944*(f[10]*alphaDrag[21]+f[11]*alphaDrag[19]+f[8]*alphaDrag[17]+f[9]*alphaDrag[16])+1.224744871391589*(alphaDrag[2]*f[11]+alphaDrag[4]*f[10])+1.369306393762915*(alphaDrag[0]*f[7]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]+f[3]*alphaDrag[4]); 
  out[12] += (1.224744871391589*f[13]+1.369306393762915*f[1])*alphaDrag[21]+1.224744871391589*f[12]*alphaDrag[19]+1.369306393762915*(f[0]*alphaDrag[19]+f[5]*alphaDrag[17]+f[3]*alphaDrag[16]); 
  out[13] += (1.224744871391589*f[12]+1.369306393762915*f[0])*alphaDrag[21]+1.224744871391589*f[13]*alphaDrag[19]+1.369306393762915*(f[1]*alphaDrag[19]+f[3]*alphaDrag[17]+f[5]*alphaDrag[16]); 
  out[14] += (1.224744871391589*f[15]+1.369306393762915*f[4])*alphaDrag[21]+1.224744871391589*f[14]*alphaDrag[19]+1.369306393762915*(f[2]*alphaDrag[19]+f[7]*alphaDrag[17]+f[6]*alphaDrag[16])+0.6123724356957944*(alphaDrag[4]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[12]); 
  out[15] += (1.224744871391589*f[14]+1.369306393762915*f[2])*alphaDrag[21]+1.224744871391589*f[15]*alphaDrag[19]+1.369306393762915*(f[4]*alphaDrag[19]+f[6]*alphaDrag[17]+f[7]*alphaDrag[16])+0.6123724356957944*(alphaDrag[2]*f[15]+alphaDrag[4]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[12]); 

  return fabs(0.8838834764831842*alphaDrag[0])+fabs(0.8838834764831842*alphaDrag[16]); 

} 
