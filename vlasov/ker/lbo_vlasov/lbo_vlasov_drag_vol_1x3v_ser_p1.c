#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvz2 = 2.0/dxv[3]; 

  double alphaDrag[120]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-2.828427124746191*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.828427124746191*nuUSum[1]-2.828427124746191*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[5] = -0.8164965809277261*dxv[1]*nuSum[1]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[40] = (2.828427124746191*nuUSum[2]-2.828427124746191*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[41] = (2.828427124746191*nuUSum[3]-2.828427124746191*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[43] = -0.8164965809277261*nuSum[0]*dxv[2]*rdvy2; 
  alphaDrag[46] = -0.8164965809277261*nuSum[1]*dxv[2]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[80] = (2.828427124746191*nuUSum[4]-2.828427124746191*nuSum[0]*w[3])*rdvz2; 
  alphaDrag[81] = (2.828427124746191*nuUSum[5]-2.828427124746191*nuSum[1]*w[3])*rdvz2; 
  alphaDrag[84] = -0.8164965809277261*nuSum[0]*dxv[3]*rdvz2; 
  alphaDrag[88] = -0.8164965809277261*nuSum[1]*dxv[3]*rdvz2; 

  out[2] += 0.4330127018922193*(alphaDrag[5]*f[5]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[6]*alphaDrag[46]+f[3]*alphaDrag[43]+f[1]*alphaDrag[41]+f[0]*alphaDrag[40]); 
  out[4] += 0.4330127018922193*(f[8]*alphaDrag[88]+f[4]*alphaDrag[84]+f[1]*alphaDrag[81]+f[0]*alphaDrag[80]); 
  out[5] += 0.4330127018922193*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.4330127018922193*(f[3]*alphaDrag[46]+f[6]*alphaDrag[43]+f[0]*alphaDrag[41]+f[1]*alphaDrag[40]); 
  out[7] += 0.4330127018922193*(f[11]*alphaDrag[46]+f[7]*alphaDrag[43]+f[5]*alphaDrag[41]+f[2]*alphaDrag[40]+alphaDrag[5]*f[11]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[4]*alphaDrag[88]+f[8]*alphaDrag[84]+f[0]*alphaDrag[81]+f[1]*alphaDrag[80]); 
  out[9] += 0.4330127018922193*(f[12]*alphaDrag[88]+f[9]*alphaDrag[84]+f[5]*alphaDrag[81]+f[2]*alphaDrag[80]+alphaDrag[5]*f[12]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[13]*alphaDrag[88]+f[10]*alphaDrag[84]+f[6]*alphaDrag[81]+f[3]*alphaDrag[80]+f[13]*alphaDrag[46]+f[10]*alphaDrag[43]+f[8]*alphaDrag[41]+f[4]*alphaDrag[40]); 
  out[11] += 0.4330127018922193*(f[7]*alphaDrag[46]+f[11]*alphaDrag[43]+f[2]*alphaDrag[41]+f[5]*alphaDrag[40]+alphaDrag[2]*f[11]+alphaDrag[5]*f[7]+alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[12] += 0.4330127018922193*(f[9]*alphaDrag[88]+f[12]*alphaDrag[84]+f[2]*alphaDrag[81]+f[5]*alphaDrag[80]+alphaDrag[2]*f[12]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[13] += 0.4330127018922193*(f[10]*alphaDrag[88]+f[13]*alphaDrag[84]+f[3]*alphaDrag[81]+f[6]*alphaDrag[80]+f[10]*alphaDrag[46]+f[13]*alphaDrag[43]+f[4]*alphaDrag[41]+f[8]*alphaDrag[40]); 
  out[14] += 0.4330127018922193*(f[15]*alphaDrag[88]+f[14]*alphaDrag[84]+f[11]*alphaDrag[81]+f[7]*alphaDrag[80]+f[15]*alphaDrag[46]+f[14]*alphaDrag[43]+f[12]*alphaDrag[41]+f[9]*alphaDrag[40]+alphaDrag[5]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[10]); 
  out[15] += 0.4330127018922193*(f[14]*alphaDrag[88]+f[15]*alphaDrag[84]+f[7]*alphaDrag[81]+f[11]*alphaDrag[80]+f[14]*alphaDrag[46]+f[15]*alphaDrag[43]+f[9]*alphaDrag[41]+f[12]*alphaDrag[40]+alphaDrag[2]*f[15]+alphaDrag[5]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[10]); 
  out[16] += 0.8660254037844386*(alphaDrag[5]*f[17]+alphaDrag[2]*f[16])+0.9682458365518543*(alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[17] += 0.8660254037844386*(alphaDrag[2]*f[17]+alphaDrag[5]*f[16])+0.9682458365518543*(alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[18] += 0.4330127018922193*(f[20]*alphaDrag[46]+f[18]*alphaDrag[43]+f[17]*alphaDrag[41]+f[16]*alphaDrag[40])+0.8660254037844386*(alphaDrag[5]*f[20]+alphaDrag[2]*f[18])+0.9682458365518543*(alphaDrag[1]*f[11]+alphaDrag[0]*f[7]+alphaDrag[5]*f[6]+alphaDrag[2]*f[3]); 
  out[19] += 0.4330127018922193*(f[21]*alphaDrag[88]+f[19]*alphaDrag[84]+f[17]*alphaDrag[81]+f[16]*alphaDrag[80])+0.8660254037844386*(alphaDrag[5]*f[21]+alphaDrag[2]*f[19])+0.9682458365518543*(alphaDrag[1]*f[12]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[20] += 0.4330127018922193*(f[18]*alphaDrag[46]+f[20]*alphaDrag[43]+f[16]*alphaDrag[41]+f[17]*alphaDrag[40])+0.8660254037844386*(alphaDrag[2]*f[20]+alphaDrag[5]*f[18])+0.9682458365518543*(alphaDrag[0]*f[11]+alphaDrag[1]*f[7]+alphaDrag[2]*f[6]+f[3]*alphaDrag[5]); 
  out[21] += 0.4330127018922193*(f[19]*alphaDrag[88]+f[21]*alphaDrag[84]+f[16]*alphaDrag[81]+f[17]*alphaDrag[80])+0.8660254037844386*(alphaDrag[2]*f[21]+alphaDrag[5]*f[19])+0.9682458365518543*(alphaDrag[0]*f[12]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[22] += 0.4330127018922193*(f[23]*alphaDrag[88]+f[22]*alphaDrag[84]+f[20]*alphaDrag[81]+f[18]*alphaDrag[80]+f[23]*alphaDrag[46]+f[22]*alphaDrag[43]+f[21]*alphaDrag[41]+f[19]*alphaDrag[40])+0.8660254037844386*(alphaDrag[5]*f[23]+alphaDrag[2]*f[22])+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[14]+alphaDrag[5]*f[13]+alphaDrag[2]*f[10]); 
  out[23] += 0.4330127018922193*(f[22]*alphaDrag[88]+f[23]*alphaDrag[84]+f[18]*alphaDrag[81]+f[20]*alphaDrag[80]+f[22]*alphaDrag[46]+f[23]*alphaDrag[43]+f[19]*alphaDrag[41]+f[21]*alphaDrag[40])+0.8660254037844386*(alphaDrag[2]*f[23]+alphaDrag[5]*f[22])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[14]+alphaDrag[2]*f[13]+alphaDrag[5]*f[10]); 
  out[24] += (0.8660254037844386*f[25]+0.9682458365518543*f[1])*alphaDrag[46]+0.8660254037844386*f[24]*alphaDrag[43]+0.9682458365518543*(f[0]*alphaDrag[43]+f[6]*alphaDrag[41]+f[3]*alphaDrag[40]); 
  out[25] += (0.8660254037844386*f[24]+0.9682458365518543*f[0])*alphaDrag[46]+0.8660254037844386*f[25]*alphaDrag[43]+0.9682458365518543*(f[1]*alphaDrag[43]+f[3]*alphaDrag[41]+f[6]*alphaDrag[40]); 
  out[26] += (0.8660254037844386*f[28]+0.9682458365518543*f[5])*alphaDrag[46]+0.8660254037844386*f[26]*alphaDrag[43]+0.9682458365518543*(f[2]*alphaDrag[43]+f[11]*alphaDrag[41]+f[7]*alphaDrag[40])+0.4330127018922193*(alphaDrag[5]*f[28]+alphaDrag[2]*f[26]+alphaDrag[1]*f[25]+alphaDrag[0]*f[24]); 
  out[27] += 0.4330127018922193*(f[29]*alphaDrag[88]+f[27]*alphaDrag[84]+f[25]*alphaDrag[81]+f[24]*alphaDrag[80])+(0.8660254037844386*f[29]+0.9682458365518543*f[8])*alphaDrag[46]+0.8660254037844386*f[27]*alphaDrag[43]+0.9682458365518543*(f[4]*alphaDrag[43]+f[13]*alphaDrag[41]+f[10]*alphaDrag[40]); 
  out[28] += (0.8660254037844386*f[26]+0.9682458365518543*f[2])*alphaDrag[46]+0.8660254037844386*f[28]*alphaDrag[43]+0.9682458365518543*(f[5]*alphaDrag[43]+f[7]*alphaDrag[41]+f[11]*alphaDrag[40])+0.4330127018922193*(alphaDrag[2]*f[28]+alphaDrag[5]*f[26]+alphaDrag[0]*f[25]+alphaDrag[1]*f[24]); 
  out[29] += 0.4330127018922193*(f[27]*alphaDrag[88]+f[29]*alphaDrag[84]+f[24]*alphaDrag[81]+f[25]*alphaDrag[80])+(0.8660254037844386*f[27]+0.9682458365518543*f[4])*alphaDrag[46]+0.8660254037844386*f[29]*alphaDrag[43]+0.9682458365518543*(f[8]*alphaDrag[43]+f[10]*alphaDrag[41]+f[13]*alphaDrag[40]); 
  out[30] += 0.4330127018922193*(f[31]*alphaDrag[88]+f[30]*alphaDrag[84]+f[28]*alphaDrag[81]+f[26]*alphaDrag[80])+(0.8660254037844386*f[31]+0.9682458365518543*f[12])*alphaDrag[46]+0.8660254037844386*f[30]*alphaDrag[43]+0.9682458365518543*(f[9]*alphaDrag[43]+f[15]*alphaDrag[41]+f[14]*alphaDrag[40])+0.4330127018922193*(alphaDrag[5]*f[31]+alphaDrag[2]*f[30]+alphaDrag[1]*f[29]+alphaDrag[0]*f[27]); 
  out[31] += 0.4330127018922193*(f[30]*alphaDrag[88]+f[31]*alphaDrag[84]+f[26]*alphaDrag[81]+f[28]*alphaDrag[80])+(0.8660254037844386*f[30]+0.9682458365518543*f[9])*alphaDrag[46]+0.8660254037844386*f[31]*alphaDrag[43]+0.9682458365518543*(f[12]*alphaDrag[43]+f[14]*alphaDrag[41]+f[15]*alphaDrag[40])+0.4330127018922193*(alphaDrag[2]*f[31]+alphaDrag[5]*f[30]+alphaDrag[0]*f[29]+alphaDrag[1]*f[27]); 
  out[32] += (0.8660254037844386*f[33]+0.9682458365518543*f[1])*alphaDrag[88]+0.8660254037844386*f[32]*alphaDrag[84]+0.9682458365518543*(f[0]*alphaDrag[84]+f[8]*alphaDrag[81]+f[4]*alphaDrag[80]); 
  out[33] += (0.8660254037844386*f[32]+0.9682458365518543*f[0])*alphaDrag[88]+0.8660254037844386*f[33]*alphaDrag[84]+0.9682458365518543*(f[1]*alphaDrag[84]+f[4]*alphaDrag[81]+f[8]*alphaDrag[80]); 
  out[34] += (0.8660254037844386*f[36]+0.9682458365518543*f[5])*alphaDrag[88]+0.8660254037844386*f[34]*alphaDrag[84]+0.9682458365518543*(f[2]*alphaDrag[84]+f[12]*alphaDrag[81]+f[9]*alphaDrag[80])+0.4330127018922193*(alphaDrag[5]*f[36]+alphaDrag[2]*f[34]+alphaDrag[1]*f[33]+alphaDrag[0]*f[32]); 
  out[35] += (0.8660254037844386*f[37]+0.9682458365518543*f[6])*alphaDrag[88]+0.8660254037844386*f[35]*alphaDrag[84]+0.9682458365518543*(f[3]*alphaDrag[84]+f[13]*alphaDrag[81]+f[10]*alphaDrag[80])+0.4330127018922193*(f[37]*alphaDrag[46]+f[35]*alphaDrag[43]+f[33]*alphaDrag[41]+f[32]*alphaDrag[40]); 
  out[36] += (0.8660254037844386*f[34]+0.9682458365518543*f[2])*alphaDrag[88]+0.8660254037844386*f[36]*alphaDrag[84]+0.9682458365518543*(f[5]*alphaDrag[84]+f[9]*alphaDrag[81]+f[12]*alphaDrag[80])+0.4330127018922193*(alphaDrag[2]*f[36]+alphaDrag[5]*f[34]+alphaDrag[0]*f[33]+alphaDrag[1]*f[32]); 
  out[37] += (0.8660254037844386*f[35]+0.9682458365518543*f[3])*alphaDrag[88]+0.8660254037844386*f[37]*alphaDrag[84]+0.9682458365518543*(f[6]*alphaDrag[84]+f[10]*alphaDrag[81]+f[13]*alphaDrag[80])+0.4330127018922193*(f[35]*alphaDrag[46]+f[37]*alphaDrag[43]+f[32]*alphaDrag[41]+f[33]*alphaDrag[40]); 
  out[38] += (0.8660254037844386*f[39]+0.9682458365518543*f[11])*alphaDrag[88]+0.8660254037844386*f[38]*alphaDrag[84]+0.9682458365518543*(f[7]*alphaDrag[84]+f[15]*alphaDrag[81]+f[14]*alphaDrag[80])+0.4330127018922193*(f[39]*alphaDrag[46]+f[38]*alphaDrag[43]+f[36]*alphaDrag[41]+f[34]*alphaDrag[40]+alphaDrag[5]*f[39]+alphaDrag[2]*f[38]+alphaDrag[1]*f[37]+alphaDrag[0]*f[35]); 
  out[39] += (0.8660254037844386*f[38]+0.9682458365518543*f[7])*alphaDrag[88]+0.8660254037844386*f[39]*alphaDrag[84]+0.9682458365518543*(f[11]*alphaDrag[84]+f[14]*alphaDrag[81]+f[15]*alphaDrag[80])+0.4330127018922193*(f[38]*alphaDrag[46]+f[39]*alphaDrag[43]+f[34]*alphaDrag[41]+f[36]*alphaDrag[40]+alphaDrag[2]*f[39]+alphaDrag[5]*f[38]+alphaDrag[0]*f[37]+alphaDrag[1]*f[35]); 

  return fabs(0.625*alphaDrag[0])+fabs(0.625*alphaDrag[40])+fabs(0.625*alphaDrag[80]); 

} 
