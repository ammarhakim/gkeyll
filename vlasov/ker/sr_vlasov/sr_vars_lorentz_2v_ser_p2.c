#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_lorentz_2v_ser_p2(const double *w, const double *dxv, double* GKYL_RESTRICT gamma_inv) 
{ 
  // w:   Cell-center coordinates.
  // dxv: Cell spacing.
  // gamma_inv: Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
 
  const double w0 = w[0]; 
  const double dv0 = dxv[0]; 
  const double w0_sq = w0*w0, dv0_sq = dv0*dv0; 
  double p0_sq[3] = {0.0};
  const double w1 = w[1]; 
  const double dv1 = dxv[1]; 
  const double w1_sq = w1*w1, dv1_sq = dv1*dv1; 
  double p1_sq[3] = {0.0};
  p0_sq[0] = 1.4142135623730951*w0_sq+0.11785113019775789*dv0_sq; 
  p0_sq[1] = 0.8164965809277261*dv0*w0; 
  p0_sq[2] = 0.10540925533894596*dv0_sq; 

  p1_sq[0] = 1.4142135623730951*w1_sq+0.11785113019775789*dv1_sq; 
  p1_sq[1] = 0.8164965809277261*dv1*w1; 
  p1_sq[2] = 0.10540925533894596*dv1_sq; 

  double gamma_nodal[8] = {0.0};
  double gamma_inv_nodal[8] = {0.0};

  gamma_nodal[0] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]+1.5811388300841895*p0_sq[2]-1.224744871391589*p1_sq[1]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[0] = 1.0/gamma_nodal[0];
  gamma_nodal[1] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]-0.7905694150420947*p0_sq[2]-1.224744871391589*p1_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[1] = 1.0/gamma_nodal[1];
  gamma_nodal[2] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]+1.5811388300841895*p0_sq[2]-1.224744871391589*p1_sq[1]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[2] = 1.0/gamma_nodal[2];
  gamma_nodal[3] = sqrt(1.0 + -(0.7905694150420947*p1_sq[2])+1.5811388300841895*p0_sq[2]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[3] = 1.0/gamma_nodal[3];
  gamma_nodal[4] = sqrt(1.0 + -(0.7905694150420947*p1_sq[2])+1.5811388300841895*p0_sq[2]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[4] = 1.0/gamma_nodal[4];
  gamma_nodal[5] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]+1.5811388300841895*p0_sq[2]+1.224744871391589*p1_sq[1]-1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[5] = 1.0/gamma_nodal[5];
  gamma_nodal[6] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]-0.7905694150420947*p0_sq[2]+1.224744871391589*p1_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[6] = 1.0/gamma_nodal[6];
  gamma_nodal[7] = sqrt(1.0 + 1.5811388300841895*p1_sq[2]+1.5811388300841895*p0_sq[2]+1.224744871391589*p1_sq[1]+1.224744871391589*p0_sq[1]+0.7071067811865475*p1_sq[0]+0.7071067811865475*p0_sq[0]);
  gamma_inv_nodal[7] = 1.0/gamma_nodal[7];

  gamma_inv[0] = -(0.16666666666666666*gamma_inv_nodal[7])+0.6666666666666666*gamma_inv_nodal[6]-0.16666666666666666*gamma_inv_nodal[5]+0.6666666666666666*gamma_inv_nodal[4]+0.6666666666666666*gamma_inv_nodal[3]-0.16666666666666666*gamma_inv_nodal[2]+0.6666666666666666*gamma_inv_nodal[1]-0.16666666666666666*gamma_inv_nodal[0]; 
  gamma_inv[1] = 0.09622504486493764*gamma_inv_nodal[7]-0.09622504486493764*gamma_inv_nodal[5]+0.3849001794597506*gamma_inv_nodal[4]-0.3849001794597506*gamma_inv_nodal[3]+0.09622504486493764*gamma_inv_nodal[2]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[2] = 0.09622504486493764*gamma_inv_nodal[7]+0.3849001794597506*gamma_inv_nodal[6]+0.09622504486493764*gamma_inv_nodal[5]-0.09622504486493764*gamma_inv_nodal[2]-0.3849001794597506*gamma_inv_nodal[1]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[3] = 0.16666666666666666*gamma_inv_nodal[7]-0.16666666666666666*gamma_inv_nodal[5]-0.16666666666666666*gamma_inv_nodal[2]+0.16666666666666666*gamma_inv_nodal[0]; 
  gamma_inv[4] = 0.14907119849998596*gamma_inv_nodal[7]-0.2981423969999719*gamma_inv_nodal[6]+0.14907119849998596*gamma_inv_nodal[5]+0.14907119849998596*gamma_inv_nodal[2]-0.2981423969999719*gamma_inv_nodal[1]+0.14907119849998596*gamma_inv_nodal[0]; 
  gamma_inv[5] = 0.14907119849998596*gamma_inv_nodal[7]+0.14907119849998596*gamma_inv_nodal[5]-0.2981423969999719*gamma_inv_nodal[4]-0.2981423969999719*gamma_inv_nodal[3]+0.14907119849998596*gamma_inv_nodal[2]+0.14907119849998596*gamma_inv_nodal[0]; 
  gamma_inv[6] = 0.08606629658238703*gamma_inv_nodal[7]-0.17213259316477406*gamma_inv_nodal[6]+0.08606629658238703*gamma_inv_nodal[5]-0.08606629658238703*gamma_inv_nodal[2]+0.17213259316477406*gamma_inv_nodal[1]-0.08606629658238703*gamma_inv_nodal[0]; 
  gamma_inv[7] = 0.08606629658238703*gamma_inv_nodal[7]-0.08606629658238703*gamma_inv_nodal[5]-0.17213259316477406*gamma_inv_nodal[4]+0.17213259316477406*gamma_inv_nodal[3]+0.08606629658238703*gamma_inv_nodal[2]-0.08606629658238703*gamma_inv_nodal[0]; 
} 
