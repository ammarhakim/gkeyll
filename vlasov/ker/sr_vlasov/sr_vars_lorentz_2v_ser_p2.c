#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_lorentz_2v_ser_p2(const double *w, const double *dxv, double* GKYL_RESTRICT gamma, double* GKYL_RESTRICT gamma_inv) 
{ 
  // w:   Cell-center coordinates.
  // dxv: Cell spacing.
  // gamma:  Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv: Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
 
  const double w0 = w[0]; 
  const double dv0 = dxv[0]; 
  const double w0_sq = w0*w0, dv0_sq = dv0*dv0; 
  double p0_sq[3] = {0.0};
  const double w1 = w[1]; 
  const double dv1 = dxv[1]; 
  const double w1_sq = w1*w1, dv1_sq = dv1*dv1; 
  double p1_sq[3] = {0.0};
  p0_sq[0] = 1.414213562373095*w0_sq+0.1178511301977579*dv0_sq; 
  p0_sq[1] = 0.8164965809277261*dv0*w0; 
  p0_sq[2] = 0.105409255338946*dv0_sq; 

  p1_sq[0] = 1.414213562373095*w1_sq+0.1178511301977579*dv1_sq; 
  p1_sq[1] = 0.8164965809277261*dv1*w1; 
  p1_sq[2] = 0.105409255338946*dv1_sq; 

  double gamma_nodal[8] = {0.0};
  double gamma_inv_nodal[8] = {0.0};

  gamma_nodal[0] = sqrt(1.0 + 1.58113883008419*p1_sq[2]+1.58113883008419*p0_sq[2]-1.224744871391589*p1_sq[1]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[0] = 1.0/gamma_nodal[0];
  gamma_nodal[1] = sqrt(1.0 + 1.58113883008419*p1_sq[2]-0.7905694150420947*p0_sq[2]-1.224744871391589*p1_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[1] = 1.0/gamma_nodal[1];
  gamma_nodal[2] = sqrt(1.0 + 1.58113883008419*p1_sq[2]+1.58113883008419*p0_sq[2]-1.224744871391589*p1_sq[1]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[2] = 1.0/gamma_nodal[2];
  gamma_nodal[3] = sqrt(1.0 + (-0.7905694150420947*p1_sq[2])+1.58113883008419*p0_sq[2]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[3] = 1.0/gamma_nodal[3];
  gamma_nodal[4] = sqrt(1.0 + (-0.7905694150420947*p1_sq[2])+1.58113883008419*p0_sq[2]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[4] = 1.0/gamma_nodal[4];
  gamma_nodal[5] = sqrt(1.0 + 1.58113883008419*p1_sq[2]+1.58113883008419*p0_sq[2]+1.224744871391589*p1_sq[1]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[5] = 1.0/gamma_nodal[5];
  gamma_nodal[6] = sqrt(1.0 + 1.58113883008419*p1_sq[2]-0.7905694150420947*p0_sq[2]+1.224744871391589*p1_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[6] = 1.0/gamma_nodal[6];
  gamma_nodal[7] = sqrt(1.0 + 1.58113883008419*p1_sq[2]+1.58113883008419*p0_sq[2]+1.224744871391589*p1_sq[1]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[7] = 1.0/gamma_nodal[7];

  gamma[0] = (-0.1666666666666667*gamma_nodal[7])+0.6666666666666666*gamma_nodal[6]-0.1666666666666667*gamma_nodal[5]+0.6666666666666666*gamma_nodal[4]+0.6666666666666666*gamma_nodal[3]-0.1666666666666667*gamma_nodal[2]+0.6666666666666666*gamma_nodal[1]-0.1666666666666667*gamma_nodal[0]; 
  gamma[1] = 0.09622504486493764*gamma_nodal[7]-0.09622504486493764*gamma_nodal[5]+0.3849001794597506*gamma_nodal[4]-0.3849001794597506*gamma_nodal[3]+0.09622504486493764*gamma_nodal[2]-0.09622504486493764*gamma_nodal[0]; 
  gamma[2] = 0.09622504486493764*gamma_nodal[7]+0.3849001794597506*gamma_nodal[6]+0.09622504486493764*gamma_nodal[5]-0.09622504486493764*gamma_nodal[2]-0.3849001794597506*gamma_nodal[1]-0.09622504486493764*gamma_nodal[0]; 
  gamma[3] = 0.1666666666666667*gamma_nodal[7]-0.1666666666666667*gamma_nodal[5]-0.1666666666666667*gamma_nodal[2]+0.1666666666666667*gamma_nodal[0]; 
  gamma[4] = 0.149071198499986*gamma_nodal[7]-0.2981423969999719*gamma_nodal[6]+0.149071198499986*gamma_nodal[5]+0.149071198499986*gamma_nodal[2]-0.2981423969999719*gamma_nodal[1]+0.149071198499986*gamma_nodal[0]; 
  gamma[5] = 0.149071198499986*gamma_nodal[7]+0.149071198499986*gamma_nodal[5]-0.2981423969999719*gamma_nodal[4]-0.2981423969999719*gamma_nodal[3]+0.149071198499986*gamma_nodal[2]+0.149071198499986*gamma_nodal[0]; 
  gamma[6] = 0.08606629658238703*gamma_nodal[7]-0.1721325931647741*gamma_nodal[6]+0.08606629658238703*gamma_nodal[5]-0.08606629658238703*gamma_nodal[2]+0.1721325931647741*gamma_nodal[1]-0.08606629658238703*gamma_nodal[0]; 
  gamma[7] = 0.08606629658238703*gamma_nodal[7]-0.08606629658238703*gamma_nodal[5]-0.1721325931647741*gamma_nodal[4]+0.1721325931647741*gamma_nodal[3]+0.08606629658238703*gamma_nodal[2]-0.08606629658238703*gamma_nodal[0]; 

  gamma_inv[0] = (-0.1666666666666667*gamma_inv_nodal[7])+0.6666666666666666*gamma_inv_nodal[6]-0.1666666666666667*gamma_inv_nodal[5]+0.6666666666666666*gamma_inv_nodal[4]+0.6666666666666666*gamma_inv_nodal[3]-0.1666666666666667*gamma_inv_nodal[2]+0.6666666666666666*gamma_inv_nodal[1]-0.1666666666666667*gamma_inv_nodal[0]; 
  gamma_inv[1] = 0.09622504486493764*gamma_inv_nodal[7]-0.09622504486493764*gamma_inv_nodal[5]+0.3849001794597506*gamma_inv_nodal[4]-0.3849001794597506*gamma_inv_nodal[3]+0.09622504486493764*gamma_inv_nodal[2]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[2] = 0.09622504486493764*gamma_inv_nodal[7]+0.3849001794597506*gamma_inv_nodal[6]+0.09622504486493764*gamma_inv_nodal[5]-0.09622504486493764*gamma_inv_nodal[2]-0.3849001794597506*gamma_inv_nodal[1]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[3] = 0.1666666666666667*gamma_inv_nodal[7]-0.1666666666666667*gamma_inv_nodal[5]-0.1666666666666667*gamma_inv_nodal[2]+0.1666666666666667*gamma_inv_nodal[0]; 
  gamma_inv[4] = 0.149071198499986*gamma_inv_nodal[7]-0.2981423969999719*gamma_inv_nodal[6]+0.149071198499986*gamma_inv_nodal[5]+0.149071198499986*gamma_inv_nodal[2]-0.2981423969999719*gamma_inv_nodal[1]+0.149071198499986*gamma_inv_nodal[0]; 
  gamma_inv[5] = 0.149071198499986*gamma_inv_nodal[7]+0.149071198499986*gamma_inv_nodal[5]-0.2981423969999719*gamma_inv_nodal[4]-0.2981423969999719*gamma_inv_nodal[3]+0.149071198499986*gamma_inv_nodal[2]+0.149071198499986*gamma_inv_nodal[0]; 
  gamma_inv[6] = 0.08606629658238703*gamma_inv_nodal[7]-0.1721325931647741*gamma_inv_nodal[6]+0.08606629658238703*gamma_inv_nodal[5]-0.08606629658238703*gamma_inv_nodal[2]+0.1721325931647741*gamma_inv_nodal[1]-0.08606629658238703*gamma_inv_nodal[0]; 
  gamma_inv[7] = 0.08606629658238703*gamma_inv_nodal[7]-0.08606629658238703*gamma_inv_nodal[5]-0.1721325931647741*gamma_inv_nodal[4]+0.1721325931647741*gamma_inv_nodal[3]+0.08606629658238703*gamma_inv_nodal[2]-0.08606629658238703*gamma_inv_nodal[0]; 
} 
