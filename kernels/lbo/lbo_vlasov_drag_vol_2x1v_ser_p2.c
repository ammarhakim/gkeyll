#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[2]; 

  double alphaDrag[20]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.414213562373095*nuUSum[0]-1.414213562373095*nuSum[0]*w[2])*rdvx2; 
  alphaDrag[1] = (1.414213562373095*nuUSum[1]-1.414213562373095*nuSum[1]*w[2])*rdvx2; 
  alphaDrag[2] = (1.414213562373095*nuUSum[2]-1.414213562373095*nuSum[2]*w[2])*rdvx2; 
  alphaDrag[3] = -0.408248290463863*nuSum[0]*dxv[2]*rdvx2; 
  alphaDrag[4] = (1.414213562373095*nuUSum[3]-1.414213562373095*w[2]*nuSum[3])*rdvx2; 
  alphaDrag[5] = -0.408248290463863*nuSum[1]*dxv[2]*rdvx2; 
  alphaDrag[6] = -0.408248290463863*dxv[2]*nuSum[2]*rdvx2; 
  alphaDrag[7] = (1.414213562373095*nuUSum[4]-1.414213562373095*w[2]*nuSum[4])*rdvx2; 
  alphaDrag[8] = (1.414213562373095*nuUSum[5]-1.414213562373095*w[2]*nuSum[5])*rdvx2; 
  alphaDrag[10] = -0.408248290463863*dxv[2]*nuSum[3]*rdvx2; 
  alphaDrag[11] = (1.414213562373095*nuUSum[6]-1.414213562373095*w[2]*nuSum[6])*rdvx2; 
  alphaDrag[12] = (1.414213562373095*nuUSum[7]-1.414213562373095*w[2]*nuSum[7])*rdvx2; 
  alphaDrag[13] = -0.408248290463863*dxv[2]*nuSum[4]*rdvx2; 
  alphaDrag[14] = -0.408248290463863*dxv[2]*nuSum[5]*rdvx2; 
  alphaDrag[17] = -0.408248290463863*dxv[2]*nuSum[6]*rdvx2; 
  alphaDrag[18] = -0.408248290463863*dxv[2]*nuSum[7]*rdvx2; 

  out[3] += 0.6123724356957944*(alphaDrag[18]*f[18]+alphaDrag[17]*f[17]+alphaDrag[14]*f[14]+alphaDrag[13]*f[13]+alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.6123724356957944*(alphaDrag[14]*f[18]+f[14]*alphaDrag[18])+0.5477225575051661*(alphaDrag[10]*f[17]+f[10]*alphaDrag[17]+alphaDrag[5]*f[13]+f[5]*alphaDrag[13])+0.6123724356957944*(alphaDrag[8]*f[12]+f[8]*alphaDrag[12])+0.5477225575051661*(alphaDrag[4]*f[11]+f[4]*alphaDrag[11])+0.6123724356957944*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10])+0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.5477225575051661*(alphaDrag[10]*f[18]+f[10]*alphaDrag[18])+0.6123724356957944*(alphaDrag[13]*f[17]+f[13]*alphaDrag[17])+0.5477225575051661*(alphaDrag[6]*f[14]+f[6]*alphaDrag[14]+alphaDrag[4]*f[12]+f[4]*alphaDrag[12])+0.6123724356957944*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[5]*f[10]+f[5]*alphaDrag[10])+0.5477225575051661*(alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+0.6123724356957944*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*alphaDrag[10]*f[19]+1.369306393762915*(alphaDrag[12]*f[18]+f[12]*alphaDrag[18]+alphaDrag[11]*f[17]+f[11]*alphaDrag[17])+1.224744871391589*(alphaDrag[6]*f[16]+alphaDrag[5]*f[15])+1.369306393762915*(alphaDrag[8]*f[14]+f[8]*alphaDrag[14]+alphaDrag[7]*f[13]+f[7]*alphaDrag[13]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10])+1.224744871391589*alphaDrag[3]*f[9]+1.369306393762915*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[10] += (0.4898979485566357*alphaDrag[17]+0.5477225575051661*alphaDrag[6])*f[18]+0.4898979485566357*f[17]*alphaDrag[18]+0.5477225575051661*(f[6]*alphaDrag[18]+alphaDrag[5]*f[17]+f[5]*alphaDrag[17]+alphaDrag[10]*f[14]+f[10]*alphaDrag[14]+alphaDrag[10]*f[13]+f[10]*alphaDrag[13])+(0.4898979485566357*alphaDrag[11]+0.5477225575051661*alphaDrag[2])*f[12]+0.4898979485566357*f[11]*alphaDrag[12]+0.5477225575051661*(f[2]*alphaDrag[12]+alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.6123724356957944*(alphaDrag[3]*f[10]+f[3]*alphaDrag[10])+0.5477225575051661*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8]+alphaDrag[4]*f[7]+f[4]*alphaDrag[7])+0.6123724356957944*(alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.5477225575051661*alphaDrag[18]*f[18]+0.3912303982179757*alphaDrag[17]*f[17]+0.6123724356957944*(alphaDrag[6]*f[17]+f[6]*alphaDrag[17])+0.3912303982179757*alphaDrag[13]*f[13]+0.6123724356957944*(alphaDrag[3]*f[13]+f[3]*alphaDrag[13])+0.5477225575051661*alphaDrag[12]*f[12]+0.3912303982179757*alphaDrag[11]*f[11]+0.6123724356957944*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11])+0.5477225575051661*alphaDrag[10]*f[10]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*(alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[1]*f[1]); 
  out[14] += 0.3912303982179757*alphaDrag[18]*f[18]+0.6123724356957944*(alphaDrag[5]*f[18]+f[5]*alphaDrag[18])+0.5477225575051661*alphaDrag[17]*f[17]+0.3912303982179757*alphaDrag[14]*f[14]+0.6123724356957944*(alphaDrag[3]*f[14]+f[3]*alphaDrag[14])+0.3912303982179757*alphaDrag[12]*f[12]+0.6123724356957944*(alphaDrag[1]*f[12]+f[1]*alphaDrag[12])+0.5477225575051661*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10])+0.3912303982179757*alphaDrag[8]*f[8]+0.6123724356957944*(alphaDrag[0]*f[8]+f[0]*alphaDrag[8])+0.5477225575051661*(alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[2]*f[2]); 
  out[15] += (1.095445115010332*alphaDrag[17]+1.224744871391589*alphaDrag[6])*f[19]+1.369306393762915*(alphaDrag[8]*f[18]+f[8]*alphaDrag[18])+1.224744871391589*(alphaDrag[4]*f[17]+f[4]*alphaDrag[17]+alphaDrag[10]*f[16])+(1.095445115010332*alphaDrag[13]+1.224744871391589*alphaDrag[3])*f[15]+1.369306393762915*(alphaDrag[12]*f[14]+f[12]*alphaDrag[14])+1.224744871391589*(alphaDrag[1]*f[13]+f[1]*alphaDrag[13]+alphaDrag[10]*f[11]+f[10]*alphaDrag[11])+1.369306393762915*(alphaDrag[2]*f[10]+f[2]*alphaDrag[10])+1.224744871391589*(alphaDrag[5]*(f[9]+f[7])+f[5]*alphaDrag[7])+1.369306393762915*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[16] += 1.095445115010332*alphaDrag[18]*f[19]+1.224744871391589*(alphaDrag[5]*f[19]+alphaDrag[4]*f[18]+f[4]*alphaDrag[18])+1.369306393762915*(alphaDrag[7]*f[17]+f[7]*alphaDrag[17])+1.095445115010332*alphaDrag[14]*f[16]+1.224744871391589*(alphaDrag[3]*f[16]+alphaDrag[10]*f[15]+alphaDrag[2]*f[14]+f[2]*alphaDrag[14])+1.369306393762915*(alphaDrag[11]*f[13]+f[11]*alphaDrag[13])+1.224744871391589*(alphaDrag[10]*f[12]+f[10]*alphaDrag[12])+1.369306393762915*(alphaDrag[1]*f[10]+f[1]*alphaDrag[10])+1.224744871391589*(alphaDrag[6]*(f[9]+f[8])+f[6]*alphaDrag[8])+1.369306393762915*(alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[17] += 0.4898979485566357*(alphaDrag[10]*f[18]+f[10]*alphaDrag[18])+(0.5477225575051661*alphaDrag[14]+0.3912303982179757*alphaDrag[13]+0.6123724356957944*alphaDrag[3])*f[17]+(0.5477225575051661*f[14]+0.3912303982179757*f[13])*alphaDrag[17]+0.6123724356957944*(f[3]*alphaDrag[17]+alphaDrag[6]*f[13]+f[6]*alphaDrag[13])+0.4898979485566357*(alphaDrag[4]*f[12]+f[4]*alphaDrag[12])+(0.5477225575051661*alphaDrag[8]+0.3912303982179757*alphaDrag[7]+0.6123724356957944*alphaDrag[0])*f[11]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alphaDrag[11]+0.5477225575051661*(alphaDrag[5]*f[10]+f[5]*alphaDrag[10])+0.6123724356957944*(alphaDrag[2]*f[7]+f[2]*alphaDrag[7])+0.5477225575051661*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[18] += (0.3912303982179757*alphaDrag[14]+0.5477225575051661*alphaDrag[13]+0.6123724356957944*alphaDrag[3])*f[18]+(0.3912303982179757*f[14]+0.5477225575051661*f[13]+0.6123724356957944*f[3])*alphaDrag[18]+0.4898979485566357*(alphaDrag[10]*f[17]+f[10]*alphaDrag[17])+0.6123724356957944*(alphaDrag[5]*f[14]+f[5]*alphaDrag[14])+(0.3912303982179757*alphaDrag[8]+0.5477225575051661*alphaDrag[7]+0.6123724356957944*alphaDrag[0])*f[12]+(0.3912303982179757*f[8]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alphaDrag[12]+0.4898979485566357*(alphaDrag[4]*f[11]+f[4]*alphaDrag[11])+0.5477225575051661*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10])+0.6123724356957944*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+0.5477225575051661*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]); 
  out[19] += (1.095445115010332*(alphaDrag[14]+alphaDrag[13])+1.224744871391589*alphaDrag[3])*f[19]+(1.095445115010332*alphaDrag[11]+1.224744871391589*alphaDrag[2])*f[18]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*alphaDrag[18]+(1.095445115010332*alphaDrag[12]+1.224744871391589*alphaDrag[1])*f[17]+1.095445115010332*(f[15]+f[12])*alphaDrag[17]+1.224744871391589*(f[1]*alphaDrag[17]+alphaDrag[5]*f[16]+alphaDrag[6]*f[15]+alphaDrag[4]*f[14]+f[4]*alphaDrag[14]+alphaDrag[4]*f[13]+f[4]*alphaDrag[13]+alphaDrag[6]*f[12]+f[6]*alphaDrag[12]+alphaDrag[5]*f[11]+f[5]*alphaDrag[11])+(1.224744871391589*(alphaDrag[8]+alphaDrag[7])+1.369306393762915*alphaDrag[0])*f[10]+1.224744871391589*(f[9]+f[8]+f[7])*alphaDrag[10]+1.369306393762915*(f[0]*alphaDrag[10]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*(alphaDrag[8]+alphaDrag[7])); 

} 
