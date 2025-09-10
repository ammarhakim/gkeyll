#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[2]; 
  const double rdvy2 = 2.0/dxv[3]; 

  double alphaDrag[64]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.0*nuSum[0]*w[2])*rdvx2; 
  alphaDrag[1] = (2.0*nuUSum[1]-2.0*nuSum[1]*w[2])*rdvx2; 
  alphaDrag[2] = (2.0*nuUSum[2]-2.0*nuSum[2]*w[2])*rdvx2; 
  alphaDrag[3] = -0.5773502691896258*nuSum[0]*dxv[2]*rdvx2; 
  alphaDrag[5] = (2.0*nuUSum[3]-2.0*w[2]*nuSum[3])*rdvx2; 
  alphaDrag[6] = -0.5773502691896258*nuSum[1]*dxv[2]*rdvx2; 
  alphaDrag[7] = -0.5773502691896258*dxv[2]*nuSum[2]*rdvx2; 
  alphaDrag[11] = -0.5773502691896258*dxv[2]*nuSum[3]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[32] = (2.0*nuUSum[4]-2.0*nuSum[0]*w[3])*rdvy2; 
  alphaDrag[33] = (2.0*nuUSum[5]-2.0*nuSum[1]*w[3])*rdvy2; 
  alphaDrag[34] = (2.0*nuUSum[6]-2.0*nuSum[2]*w[3])*rdvy2; 
  alphaDrag[36] = -0.5773502691896258*nuSum[0]*dxv[3]*rdvy2; 
  alphaDrag[37] = (2.0*nuUSum[7]-2.0*nuSum[3]*w[3])*rdvy2; 
  alphaDrag[40] = -0.5773502691896258*nuSum[1]*dxv[3]*rdvy2; 
  alphaDrag[41] = -0.5773502691896258*nuSum[2]*dxv[3]*rdvy2; 
  alphaDrag[44] = -0.5773502691896258*dxv[3]*nuSum[3]*rdvy2; 

  out[3] += 0.4330127018922193*(alphaDrag[11]*f[11]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alphaDrag[44]+f[9]*alphaDrag[41]+f[8]*alphaDrag[40]+f[5]*alphaDrag[37]+f[4]*alphaDrag[36]+f[2]*alphaDrag[34]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[6] += 0.4330127018922193*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.4330127018922193*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.4330127018922193*(f[9]*alphaDrag[44]+f[12]*alphaDrag[41]+f[4]*alphaDrag[40]+f[2]*alphaDrag[37]+f[8]*alphaDrag[36]+f[5]*alphaDrag[34]+f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[9] += 0.4330127018922193*(f[8]*alphaDrag[44]+f[4]*alphaDrag[41]+f[12]*alphaDrag[40]+f[1]*alphaDrag[37]+f[9]*alphaDrag[36]+f[0]*alphaDrag[34]+f[5]*alphaDrag[33]+f[2]*alphaDrag[32]); 
  out[10] += 0.4330127018922193*(f[15]*alphaDrag[44]+f[14]*alphaDrag[41]+f[13]*alphaDrag[40]+f[11]*alphaDrag[37]+f[10]*alphaDrag[36]+f[7]*alphaDrag[34]+f[6]*alphaDrag[33]+f[3]*alphaDrag[32]+alphaDrag[11]*f[15]+alphaDrag[7]*f[14]+alphaDrag[6]*f[13]+alphaDrag[5]*f[12]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[12] += 0.4330127018922193*(f[4]*alphaDrag[44]+f[8]*alphaDrag[41]+f[9]*alphaDrag[40]+f[0]*alphaDrag[37]+f[12]*alphaDrag[36]+f[1]*alphaDrag[34]+f[2]*alphaDrag[33]+f[5]*alphaDrag[32]); 
  out[13] += 0.4330127018922193*(f[14]*alphaDrag[44]+f[15]*alphaDrag[41]+f[10]*alphaDrag[40]+f[7]*alphaDrag[37]+f[13]*alphaDrag[36]+f[11]*alphaDrag[34]+f[3]*alphaDrag[33]+f[6]*alphaDrag[32]+alphaDrag[7]*f[15]+alphaDrag[11]*f[14]+alphaDrag[3]*f[13]+alphaDrag[2]*f[12]+alphaDrag[6]*f[10]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[14] += 0.4330127018922193*(f[13]*alphaDrag[44]+f[10]*alphaDrag[41]+f[15]*alphaDrag[40]+f[6]*alphaDrag[37]+f[14]*alphaDrag[36]+f[3]*alphaDrag[34]+f[11]*alphaDrag[33]+f[7]*alphaDrag[32]+alphaDrag[6]*f[15]+alphaDrag[3]*f[14]+alphaDrag[11]*f[13]+alphaDrag[1]*f[12]+alphaDrag[7]*f[10]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[15] += 0.4330127018922193*(f[10]*alphaDrag[44]+f[13]*alphaDrag[41]+f[14]*alphaDrag[40]+f[3]*alphaDrag[37]+f[15]*alphaDrag[36]+f[6]*alphaDrag[34]+f[7]*alphaDrag[33]+f[11]*alphaDrag[32]+alphaDrag[3]*f[15]+alphaDrag[6]*f[14]+alphaDrag[7]*f[13]+alphaDrag[0]*f[12]+f[10]*alphaDrag[11]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[16] += 0.8660254037844386*(alphaDrag[11]*f[20]+alphaDrag[7]*f[18]+alphaDrag[6]*f[17]+alphaDrag[3]*f[16])+0.9682458365518543*(alphaDrag[5]*f[11]+f[5]*alphaDrag[11]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[17] += 0.8660254037844386*(alphaDrag[7]*f[20]+alphaDrag[11]*f[18]+alphaDrag[3]*f[17]+alphaDrag[6]*f[16])+0.9682458365518543*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11]+alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[18] += 0.8660254037844386*(alphaDrag[6]*f[20]+alphaDrag[3]*f[18]+alphaDrag[11]*f[17]+alphaDrag[7]*f[16])+0.9682458365518543*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[19] += 0.4330127018922193*(f[23]*alphaDrag[44]+f[22]*alphaDrag[41]+f[21]*alphaDrag[40]+f[20]*alphaDrag[37]+f[19]*alphaDrag[36]+f[18]*alphaDrag[34]+f[17]*alphaDrag[33]+f[16]*alphaDrag[32])+0.8660254037844386*(alphaDrag[11]*f[23]+alphaDrag[7]*f[22]+alphaDrag[6]*f[21]+alphaDrag[3]*f[19])+0.9682458365518543*(alphaDrag[5]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[11]*f[12]+alphaDrag[0]*f[10]+alphaDrag[7]*f[9]+alphaDrag[6]*f[8]+alphaDrag[3]*f[4]); 
  out[20] += 0.8660254037844386*(alphaDrag[3]*f[20]+alphaDrag[6]*f[18]+alphaDrag[7]*f[17]+alphaDrag[11]*f[16])+0.9682458365518543*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]); 
  out[21] += 0.4330127018922193*(f[22]*alphaDrag[44]+f[23]*alphaDrag[41]+f[19]*alphaDrag[40]+f[18]*alphaDrag[37]+f[21]*alphaDrag[36]+f[20]*alphaDrag[34]+f[16]*alphaDrag[33]+f[17]*alphaDrag[32])+0.8660254037844386*(alphaDrag[7]*f[23]+alphaDrag[11]*f[22]+alphaDrag[3]*f[21]+alphaDrag[6]*f[19])+0.9682458365518543*(alphaDrag[2]*f[15]+alphaDrag[5]*f[14]+alphaDrag[0]*f[13]+alphaDrag[7]*f[12]+f[9]*alphaDrag[11]+alphaDrag[1]*f[10]+alphaDrag[3]*f[8]+f[4]*alphaDrag[6]); 
  out[22] += 0.4330127018922193*(f[21]*alphaDrag[44]+f[19]*alphaDrag[41]+f[23]*alphaDrag[40]+f[17]*alphaDrag[37]+f[22]*alphaDrag[36]+f[16]*alphaDrag[34]+f[20]*alphaDrag[33]+f[18]*alphaDrag[32])+0.8660254037844386*(alphaDrag[6]*f[23]+alphaDrag[3]*f[22]+alphaDrag[11]*f[21]+alphaDrag[7]*f[19])+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[14]+alphaDrag[5]*f[13]+alphaDrag[6]*f[12]+f[8]*alphaDrag[11]+alphaDrag[2]*f[10]+alphaDrag[3]*f[9]+f[4]*alphaDrag[7]); 
  out[23] += 0.4330127018922193*(f[19]*alphaDrag[44]+f[21]*alphaDrag[41]+f[22]*alphaDrag[40]+f[16]*alphaDrag[37]+f[23]*alphaDrag[36]+f[17]*alphaDrag[34]+f[18]*alphaDrag[33]+f[20]*alphaDrag[32])+0.8660254037844386*(alphaDrag[3]*f[23]+alphaDrag[6]*f[22]+alphaDrag[7]*f[21]+alphaDrag[11]*f[19])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[14]+alphaDrag[2]*f[13]+alphaDrag[3]*f[12]+f[4]*alphaDrag[11]+alphaDrag[5]*f[10]+alphaDrag[6]*f[9]+alphaDrag[7]*f[8]); 
  out[24] += (0.8660254037844386*f[28]+0.9682458365518543*f[5])*alphaDrag[44]+(0.8660254037844386*f[26]+0.9682458365518543*f[2])*alphaDrag[41]+0.8660254037844386*f[25]*alphaDrag[40]+0.9682458365518543*(f[1]*alphaDrag[40]+f[12]*alphaDrag[37])+0.8660254037844386*f[24]*alphaDrag[36]+0.9682458365518543*(f[0]*alphaDrag[36]+f[9]*alphaDrag[34]+f[8]*alphaDrag[33]+f[4]*alphaDrag[32]); 
  out[25] += (0.8660254037844386*f[26]+0.9682458365518543*f[2])*alphaDrag[44]+(0.8660254037844386*f[28]+0.9682458365518543*f[5])*alphaDrag[41]+0.8660254037844386*f[24]*alphaDrag[40]+0.9682458365518543*(f[0]*alphaDrag[40]+f[9]*alphaDrag[37])+0.8660254037844386*f[25]*alphaDrag[36]+0.9682458365518543*(f[1]*alphaDrag[36]+f[12]*alphaDrag[34]+f[4]*alphaDrag[33]+f[8]*alphaDrag[32]); 
  out[26] += (0.8660254037844386*f[25]+0.9682458365518543*f[1])*alphaDrag[44]+(0.8660254037844386*f[24]+0.9682458365518543*f[0])*alphaDrag[41]+0.8660254037844386*f[28]*alphaDrag[40]+0.9682458365518543*(f[5]*alphaDrag[40]+f[8]*alphaDrag[37])+0.8660254037844386*f[26]*alphaDrag[36]+0.9682458365518543*(f[2]*alphaDrag[36]+f[4]*alphaDrag[34]+f[12]*alphaDrag[33]+f[9]*alphaDrag[32]); 
  out[27] += (0.8660254037844386*f[31]+0.9682458365518543*f[11])*alphaDrag[44]+(0.8660254037844386*f[30]+0.9682458365518543*f[7])*alphaDrag[41]+0.8660254037844386*f[29]*alphaDrag[40]+0.9682458365518543*(f[6]*alphaDrag[40]+f[15]*alphaDrag[37])+0.8660254037844386*f[27]*alphaDrag[36]+0.9682458365518543*(f[3]*alphaDrag[36]+f[14]*alphaDrag[34]+f[13]*alphaDrag[33]+f[10]*alphaDrag[32])+0.4330127018922193*(alphaDrag[11]*f[31]+alphaDrag[7]*f[30]+alphaDrag[6]*f[29]+alphaDrag[5]*f[28]+alphaDrag[3]*f[27]+alphaDrag[2]*f[26]+alphaDrag[1]*f[25]+alphaDrag[0]*f[24]); 
  out[28] += (0.8660254037844386*f[24]+0.9682458365518543*f[0])*alphaDrag[44]+(0.8660254037844386*f[25]+0.9682458365518543*f[1])*alphaDrag[41]+0.8660254037844386*f[26]*alphaDrag[40]+0.9682458365518543*(f[2]*alphaDrag[40]+f[4]*alphaDrag[37])+0.8660254037844386*f[28]*alphaDrag[36]+0.9682458365518543*(f[5]*alphaDrag[36]+f[8]*alphaDrag[34]+f[9]*alphaDrag[33]+f[12]*alphaDrag[32]); 
  out[29] += (0.8660254037844386*f[30]+0.9682458365518543*f[7])*alphaDrag[44]+(0.8660254037844386*f[31]+0.9682458365518543*f[11])*alphaDrag[41]+0.8660254037844386*f[27]*alphaDrag[40]+0.9682458365518543*(f[3]*alphaDrag[40]+f[14]*alphaDrag[37])+0.8660254037844386*f[29]*alphaDrag[36]+0.9682458365518543*(f[6]*alphaDrag[36]+f[15]*alphaDrag[34]+f[10]*alphaDrag[33]+f[13]*alphaDrag[32])+0.4330127018922193*(alphaDrag[7]*f[31]+alphaDrag[11]*f[30]+alphaDrag[3]*f[29]+alphaDrag[2]*f[28]+alphaDrag[6]*f[27]+alphaDrag[5]*f[26]+alphaDrag[0]*f[25]+alphaDrag[1]*f[24]); 
  out[30] += (0.8660254037844386*f[29]+0.9682458365518543*f[6])*alphaDrag[44]+(0.8660254037844386*f[27]+0.9682458365518543*f[3])*alphaDrag[41]+0.8660254037844386*f[31]*alphaDrag[40]+0.9682458365518543*(f[11]*alphaDrag[40]+f[13]*alphaDrag[37])+0.8660254037844386*f[30]*alphaDrag[36]+0.9682458365518543*(f[7]*alphaDrag[36]+f[10]*alphaDrag[34]+f[15]*alphaDrag[33]+f[14]*alphaDrag[32])+0.4330127018922193*(alphaDrag[6]*f[31]+alphaDrag[3]*f[30]+alphaDrag[11]*f[29]+alphaDrag[1]*f[28]+alphaDrag[7]*f[27]+alphaDrag[0]*f[26]+alphaDrag[5]*f[25]+alphaDrag[2]*f[24]); 
  out[31] += (0.8660254037844386*f[27]+0.9682458365518543*f[3])*alphaDrag[44]+(0.8660254037844386*f[29]+0.9682458365518543*f[6])*alphaDrag[41]+0.8660254037844386*f[30]*alphaDrag[40]+0.9682458365518543*(f[7]*alphaDrag[40]+f[10]*alphaDrag[37])+0.8660254037844386*f[31]*alphaDrag[36]+0.9682458365518543*(f[11]*alphaDrag[36]+f[13]*alphaDrag[34]+f[14]*alphaDrag[33]+f[15]*alphaDrag[32])+0.4330127018922193*(alphaDrag[3]*f[31]+alphaDrag[6]*f[30]+alphaDrag[7]*f[29]+alphaDrag[0]*f[28]+alphaDrag[11]*f[27]+alphaDrag[1]*f[26]+alphaDrag[2]*f[25]+alphaDrag[5]*f[24]); 

  return fabs(0.625*alphaDrag[0])+fabs(0.625*alphaDrag[32]); 

} 
