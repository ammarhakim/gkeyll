#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double
fpo_vlasov_drag_vol_1x3v_ser_p1(const double* w, const double* dxv,
  const double* h, const double* f, double* GKYL_RESTRICT out) 
{
  // w[4]: Cell-center coordinates
  // dxv[4]: Cell spacing
  // h: Input Rosenbluth potential
  // f: Input distribution function
  // out: Incremented output
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvz2 = 2.0/dxv[3]; 

  double alphaDrag[48]; 
  // Expand rdv2*(dh/dvx) in phase basis.
  alphaDrag[0] = 1.732050807568877*h[2]*rdvx2; 
  alphaDrag[1] = 1.732050807568877*h[5]*rdvx2; 
  alphaDrag[3] = 1.732050807568877*h[7]*rdvx2; 
  alphaDrag[4] = 1.732050807568877*h[9]*rdvx2; 
  alphaDrag[6] = 1.732050807568877*h[11]*rdvx2; 
  alphaDrag[8] = 1.732050807568877*h[12]*rdvx2; 
  alphaDrag[10] = 1.732050807568877*h[14]*rdvx2; 
  alphaDrag[13] = 1.732050807568877*h[15]*rdvx2; 

  // Expand rdv2*(dh/dvy) in phase basis.
  alphaDrag[16] = 1.732050807568877*h[3]*rdvy2; 
  alphaDrag[17] = 1.732050807568877*h[6]*rdvy2; 
  alphaDrag[18] = 1.732050807568877*h[7]*rdvy2; 
  alphaDrag[20] = 1.732050807568877*h[10]*rdvy2; 
  alphaDrag[21] = 1.732050807568877*h[11]*rdvy2; 
  alphaDrag[24] = 1.732050807568877*h[13]*rdvy2; 
  alphaDrag[25] = 1.732050807568877*h[14]*rdvy2; 
  alphaDrag[28] = 1.732050807568877*h[15]*rdvy2; 

  // Expand rdv2*(dh/dvz) in phase basis.
  alphaDrag[32] = 1.732050807568877*h[4]*rdvz2; 
  alphaDrag[33] = 1.732050807568877*h[8]*rdvz2; 
  alphaDrag[34] = 1.732050807568877*h[9]*rdvz2; 
  alphaDrag[35] = 1.732050807568877*h[10]*rdvz2; 
  alphaDrag[37] = 1.732050807568877*h[12]*rdvz2; 
  alphaDrag[38] = 1.732050807568877*h[13]*rdvz2; 
  alphaDrag[39] = 1.732050807568877*h[14]*rdvz2; 
  alphaDrag[43] = 1.732050807568877*h[15]*rdvz2; 

  out[2] += 0.4330127018922193*(alphaDrag[13]*f[13]+alphaDrag[10]*f[10]+alphaDrag[8]*f[8]+alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[12]*alphaDrag[28]+f[9]*alphaDrag[25]+f[8]*alphaDrag[24]+f[5]*alphaDrag[21]+f[4]*alphaDrag[20]+f[2]*alphaDrag[18]+f[1]*alphaDrag[17]+f[0]*alphaDrag[16]); 
  out[4] += 0.4330127018922193*(f[11]*alphaDrag[43]+f[7]*alphaDrag[39]+f[6]*alphaDrag[38]+f[5]*alphaDrag[37]+f[3]*alphaDrag[35]+f[2]*alphaDrag[34]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[5] += 0.4330127018922193*(alphaDrag[10]*f[13]+f[10]*alphaDrag[13]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.4330127018922193*(f[9]*alphaDrag[28]+f[12]*alphaDrag[25]+f[4]*alphaDrag[24]+f[2]*alphaDrag[21]+f[8]*alphaDrag[20]+f[5]*alphaDrag[18]+f[0]*alphaDrag[17]+f[1]*alphaDrag[16]); 
  out[7] += 0.4330127018922193*(f[8]*alphaDrag[28]+f[4]*alphaDrag[25]+f[12]*alphaDrag[24]+f[1]*alphaDrag[21]+f[9]*alphaDrag[20]+f[0]*alphaDrag[18]+f[5]*alphaDrag[17]+f[2]*alphaDrag[16]+alphaDrag[8]*f[13]+f[8]*alphaDrag[13]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[8] += 0.4330127018922193*(f[7]*alphaDrag[43]+f[11]*alphaDrag[39]+f[3]*alphaDrag[38]+f[2]*alphaDrag[37]+f[6]*alphaDrag[35]+f[5]*alphaDrag[34]+f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[9] += 0.4330127018922193*(f[6]*alphaDrag[43]+f[3]*alphaDrag[39]+f[11]*alphaDrag[38]+f[1]*alphaDrag[37]+f[7]*alphaDrag[35]+f[0]*alphaDrag[34]+f[5]*alphaDrag[33]+f[2]*alphaDrag[32]+alphaDrag[6]*f[13]+f[6]*alphaDrag[13]+alphaDrag[3]*f[10]+f[3]*alphaDrag[10]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]); 
  out[10] += 0.4330127018922193*(f[5]*alphaDrag[43]+f[2]*alphaDrag[39]+f[1]*alphaDrag[38]+f[11]*alphaDrag[37]+f[0]*alphaDrag[35]+f[7]*alphaDrag[34]+f[6]*alphaDrag[33]+f[3]*alphaDrag[32]+f[5]*alphaDrag[28]+f[2]*alphaDrag[25]+f[1]*alphaDrag[24]+f[12]*alphaDrag[21]+f[0]*alphaDrag[20]+f[9]*alphaDrag[18]+f[8]*alphaDrag[17]+f[4]*alphaDrag[16]); 
  out[11] += 0.4330127018922193*(f[4]*alphaDrag[28]+f[8]*alphaDrag[25]+f[9]*alphaDrag[24]+f[0]*alphaDrag[21]+f[12]*alphaDrag[20]+f[1]*alphaDrag[18]+f[2]*alphaDrag[17]+f[5]*alphaDrag[16]+alphaDrag[4]*f[13]+f[4]*alphaDrag[13]+alphaDrag[8]*f[10]+f[8]*alphaDrag[10]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[12] += 0.4330127018922193*(f[3]*alphaDrag[43]+f[6]*alphaDrag[39]+f[7]*alphaDrag[38]+f[0]*alphaDrag[37]+f[11]*alphaDrag[35]+f[1]*alphaDrag[34]+f[2]*alphaDrag[33]+f[5]*alphaDrag[32]+alphaDrag[3]*f[13]+f[3]*alphaDrag[13]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[0]*f[8]+f[0]*alphaDrag[8]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[13] += 0.4330127018922193*(f[2]*alphaDrag[43]+f[5]*alphaDrag[39]+f[0]*alphaDrag[38]+f[7]*alphaDrag[37]+f[1]*alphaDrag[35]+f[11]*alphaDrag[34]+f[3]*alphaDrag[33]+f[6]*alphaDrag[32]+f[2]*alphaDrag[28]+f[5]*alphaDrag[25]+f[0]*alphaDrag[24]+f[9]*alphaDrag[21]+f[1]*alphaDrag[20]+f[12]*alphaDrag[18]+f[4]*alphaDrag[17]+f[8]*alphaDrag[16]); 
  out[14] += 0.4330127018922193*(f[1]*alphaDrag[43]+f[0]*alphaDrag[39]+f[5]*alphaDrag[38]+f[6]*alphaDrag[37]+f[2]*alphaDrag[35]+f[3]*alphaDrag[34]+f[11]*alphaDrag[33]+f[7]*alphaDrag[32]+f[1]*alphaDrag[28]+f[0]*alphaDrag[25]+f[5]*alphaDrag[24]+f[8]*alphaDrag[21]+f[2]*alphaDrag[20]+f[4]*alphaDrag[18]+f[12]*alphaDrag[17]+f[9]*alphaDrag[16]+alphaDrag[1]*f[13]+f[1]*alphaDrag[13]+alphaDrag[0]*f[10]+f[0]*alphaDrag[10]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 
  out[15] += 0.4330127018922193*(f[0]*alphaDrag[43]+f[1]*alphaDrag[39]+f[2]*alphaDrag[38]+f[3]*alphaDrag[37]+f[5]*alphaDrag[35]+f[6]*alphaDrag[34]+f[7]*alphaDrag[33]+f[11]*alphaDrag[32]+f[0]*alphaDrag[28]+f[1]*alphaDrag[25]+f[2]*alphaDrag[24]+f[4]*alphaDrag[21]+f[5]*alphaDrag[20]+f[8]*alphaDrag[18]+f[9]*alphaDrag[17]+f[12]*alphaDrag[16]+alphaDrag[0]*f[13]+f[0]*alphaDrag[13]+alphaDrag[1]*f[10]+f[1]*alphaDrag[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[4]*f[6]+f[4]*alphaDrag[6]); 

  return fabs(0.125*alphaDrag[0])+fabs(0.125*alphaDrag[16])+fabs(0.125*alphaDrag[32]); 

} 
