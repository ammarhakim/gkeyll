#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_2x2v_ser_p1(const double *w, const double *dxv, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[2]*dxv[3]/4; 
 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double *V_1 = &u_i[4]; 
  const double *V_1_sq = &u_i_sq[4]; 
 
  double temp[32] = {0.0}; 
  double temp_sq[32] = {0.0}; 
  double p_fac[32] = {0.0}; 
  temp[0] = 2.0*V_1[0]*wx2+2.0*V_0[0]*wx1; 
  temp[1] = 2.0*V_1[1]*wx2+2.0*V_0[1]*wx1; 
  temp[2] = 2.0*V_1[2]*wx2+2.0*V_0[2]*wx1; 
  temp[3] = 0.5773502691896258*V_0[0]*dv1; 
  temp[4] = 0.5773502691896258*V_1[0]*dv2; 
  temp[5] = 2.0*V_1[3]*wx2+2.0*V_0[3]*wx1; 
  temp[6] = 0.5773502691896258*V_0[1]*dv1; 
  temp[7] = 0.5773502691896258*V_0[2]*dv1; 
  temp[8] = 0.5773502691896258*V_1[1]*dv2; 
  temp[9] = 0.5773502691896258*V_1[2]*dv2; 
  temp[11] = 0.5773502691896258*V_0[3]*dv1; 
  temp[12] = 0.5773502691896258*V_1[3]*dv2; 

  temp_sq[0] = 2.0*V_1_sq[0]*wx2_sq+2.0*V_0[3]*V_1[3]*wx1*wx2+2.0*V_0[2]*V_1[2]*wx1*wx2+2.0*V_0[1]*V_1[1]*wx1*wx2+2.0*V_0[0]*V_1[0]*wx1*wx2+2.0*V_0_sq[0]*wx1_sq+0.1666666666666667*V_1_sq[0]*dv2_sq+0.1666666666666667*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 2.0*V_1_sq[1]*wx2_sq+2.0*V_0[2]*V_1[3]*wx1*wx2+2.0*V_1[2]*V_0[3]*wx1*wx2+2.0*V_0[0]*V_1[1]*wx1*wx2+2.0*V_1[0]*V_0[1]*wx1*wx2+2.0*V_0_sq[1]*wx1_sq+0.1666666666666667*V_1_sq[1]*dv2_sq+0.1666666666666667*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 2.0*V_1_sq[2]*wx2_sq+2.0*V_0[1]*V_1[3]*wx1*wx2+2.0*V_1[1]*V_0[3]*wx1*wx2+2.0*V_0[0]*V_1[2]*wx1*wx2+2.0*V_1[0]*V_0[2]*wx1*wx2+2.0*V_0_sq[2]*wx1_sq+0.1666666666666667*V_1_sq[2]*dv2_sq+0.1666666666666667*V_0_sq[2]*dv1_sq; 
  temp_sq[3] = 0.5773502691896258*V_0[3]*V_1[3]*dv1*wx2+0.5773502691896258*V_0[2]*V_1[2]*dv1*wx2+0.5773502691896258*V_0[1]*V_1[1]*dv1*wx2+0.5773502691896258*V_0[0]*V_1[0]*dv1*wx2+1.154700538379252*V_0_sq[0]*dv1*wx1; 
  temp_sq[4] = 1.154700538379252*V_1_sq[0]*dv2*wx2+0.5773502691896258*V_0[3]*V_1[3]*dv2*wx1+0.5773502691896258*V_0[2]*V_1[2]*dv2*wx1+0.5773502691896258*V_0[1]*V_1[1]*dv2*wx1+0.5773502691896258*V_0[0]*V_1[0]*dv2*wx1; 
  temp_sq[5] = 2.0*V_1_sq[3]*wx2_sq+2.0*V_0[0]*V_1[3]*wx1*wx2+2.0*V_1[0]*V_0[3]*wx1*wx2+2.0*V_0[1]*V_1[2]*wx1*wx2+2.0*V_1[1]*V_0[2]*wx1*wx2+2.0*V_0_sq[3]*wx1_sq+0.1666666666666667*V_1_sq[3]*dv2_sq+0.1666666666666667*V_0_sq[3]*dv1_sq; 
  temp_sq[6] = 0.5773502691896258*V_0[2]*V_1[3]*dv1*wx2+0.5773502691896258*V_1[2]*V_0[3]*dv1*wx2+0.5773502691896258*V_0[0]*V_1[1]*dv1*wx2+0.5773502691896258*V_1[0]*V_0[1]*dv1*wx2+1.154700538379252*V_0_sq[1]*dv1*wx1; 
  temp_sq[7] = 0.5773502691896258*V_0[1]*V_1[3]*dv1*wx2+0.5773502691896258*V_1[1]*V_0[3]*dv1*wx2+0.5773502691896258*V_0[0]*V_1[2]*dv1*wx2+0.5773502691896258*V_1[0]*V_0[2]*dv1*wx2+1.154700538379252*V_0_sq[2]*dv1*wx1; 
  temp_sq[8] = 1.154700538379252*V_1_sq[1]*dv2*wx2+0.5773502691896258*V_0[2]*V_1[3]*dv2*wx1+0.5773502691896258*V_1[2]*V_0[3]*dv2*wx1+0.5773502691896258*V_0[0]*V_1[1]*dv2*wx1+0.5773502691896258*V_1[0]*V_0[1]*dv2*wx1; 
  temp_sq[9] = 1.154700538379252*V_1_sq[2]*dv2*wx2+0.5773502691896258*V_0[1]*V_1[3]*dv2*wx1+0.5773502691896258*V_1[1]*V_0[3]*dv2*wx1+0.5773502691896258*V_0[0]*V_1[2]*dv2*wx1+0.5773502691896258*V_1[0]*V_0[2]*dv2*wx1; 
  temp_sq[10] = 0.1666666666666667*V_0[3]*V_1[3]*dv1*dv2+0.1666666666666667*V_0[2]*V_1[2]*dv1*dv2+0.1666666666666667*V_0[1]*V_1[1]*dv1*dv2+0.1666666666666667*V_0[0]*V_1[0]*dv1*dv2; 
  temp_sq[11] = 0.5773502691896258*V_0[0]*V_1[3]*dv1*wx2+0.5773502691896258*V_1[0]*V_0[3]*dv1*wx2+0.5773502691896258*V_0[1]*V_1[2]*dv1*wx2+0.5773502691896258*V_1[1]*V_0[2]*dv1*wx2+1.154700538379252*V_0_sq[3]*dv1*wx1; 
  temp_sq[12] = 1.154700538379252*V_1_sq[3]*dv2*wx2+0.5773502691896258*V_0[0]*V_1[3]*dv2*wx1+0.5773502691896258*V_1[0]*V_0[3]*dv2*wx1+0.5773502691896258*V_0[1]*V_1[2]*dv2*wx1+0.5773502691896258*V_1[1]*V_0[2]*dv2*wx1; 
  temp_sq[13] = 0.1666666666666667*V_0[2]*V_1[3]*dv1*dv2+0.1666666666666667*V_1[2]*V_0[3]*dv1*dv2+0.1666666666666667*V_0[0]*V_1[1]*dv1*dv2+0.1666666666666667*V_1[0]*V_0[1]*dv1*dv2; 
  temp_sq[14] = 0.1666666666666667*V_0[1]*V_1[3]*dv1*dv2+0.1666666666666667*V_1[1]*V_0[3]*dv1*dv2+0.1666666666666667*V_0[0]*V_1[2]*dv1*dv2+0.1666666666666667*V_1[0]*V_0[2]*dv1*dv2; 
  temp_sq[15] = 0.1666666666666667*V_0[0]*V_1[3]*dv1*dv2+0.1666666666666667*V_1[0]*V_0[3]*dv1*dv2+0.1666666666666667*V_0[1]*V_1[2]*dv1*dv2+0.1666666666666667*V_1[1]*V_0[2]*dv1*dv2; 
  temp_sq[16] = 0.149071198499986*V_0_sq[0]*dv1_sq; 
  temp_sq[17] = 0.149071198499986*V_0_sq[1]*dv1_sq; 
  temp_sq[18] = 0.149071198499986*V_0_sq[2]*dv1_sq; 
  temp_sq[20] = 0.149071198499986*V_0_sq[3]*dv1_sq; 
  temp_sq[24] = 0.149071198499986*V_1_sq[0]*dv2_sq; 
  temp_sq[25] = 0.149071198499986*V_1_sq[1]*dv2_sq; 
  temp_sq[26] = 0.149071198499986*V_1_sq[2]*dv2_sq; 
  temp_sq[28] = 0.149071198499986*V_1_sq[3]*dv2_sq; 

  p_fac[0] = 0.5*gamma_inv[5]*temp_sq[24]+0.5*gamma_inv[4]*temp_sq[16]+0.5*gamma_inv[3]*temp_sq[10]-1.0*GammaV[3]*temp[5]+0.5*gamma_inv[2]*temp_sq[4]+0.5*gamma_inv[1]*temp_sq[3]-1.0*GammaV[2]*temp[2]-1.0*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.5*gamma_inv[0]*temp_sq[0]-1.0*GammaV[0]*temp[0]-2.0*gamma_inv[0]; 
  p_fac[1] = 0.5000000000000001*gamma_inv[5]*temp_sq[25]+0.5000000000000001*gamma_inv[4]*temp_sq[17]+0.5*gamma_inv[3]*temp_sq[13]+0.5*gamma_inv[2]*temp_sq[8]+0.5*gamma_inv[1]*temp_sq[6]-1.0*GammaV[2]*temp[5]-1.0*temp[2]*GammaV[3]+0.5*gamma_inv[0]*temp_sq[1]-1.0*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.0*temp[0]*GammaV[1]; 
  p_fac[2] = 0.5000000000000001*gamma_inv[5]*temp_sq[26]+0.5000000000000001*gamma_inv[4]*temp_sq[18]+0.5*gamma_inv[3]*temp_sq[14]+0.5*gamma_inv[2]*temp_sq[9]+0.5*gamma_inv[1]*temp_sq[7]-1.0*GammaV[1]*temp[5]-1.0*temp[1]*GammaV[3]+0.5*gamma_inv[0]*temp_sq[2]-1.0*GammaV[0]*temp[2]+gamma[0]*GammaV_sq[2]-1.0*temp[0]*GammaV[2]; 
  p_fac[3] = 0.5000000000000001*gamma_inv[7]*temp_sq[24]+0.4472135954999579*gamma_inv[1]*temp_sq[16]-1.0*GammaV[3]*temp[11]+0.447213595499958*gamma_inv[6]*temp_sq[10]+0.5*gamma_inv[2]*temp_sq[10]-1.0*GammaV[2]*temp[7]-1.0*GammaV[1]*temp[6]+0.5*gamma_inv[3]*temp_sq[4]+0.4472135954999579*temp_sq[3]*gamma_inv[4]+0.5*gamma_inv[0]*temp_sq[3]-1.0*GammaV[0]*temp[3]+GammaV_sq[0]*gamma[1]+0.5*temp_sq[0]*gamma_inv[1]-2.0*gamma_inv[1]; 
  p_fac[4] = 0.4472135954999579*gamma_inv[2]*temp_sq[24]+0.5000000000000001*gamma_inv[6]*temp_sq[16]-1.0*GammaV[3]*temp[12]+0.447213595499958*gamma_inv[7]*temp_sq[10]+0.5*gamma_inv[1]*temp_sq[10]-1.0*GammaV[2]*temp[9]-1.0*GammaV[1]*temp[8]+0.4472135954999579*temp_sq[4]*gamma_inv[5]+0.5*gamma_inv[0]*temp_sq[4]-1.0*GammaV[0]*temp[4]+0.5*gamma_inv[3]*temp_sq[3]+GammaV_sq[0]*gamma[2]+0.5*temp_sq[0]*gamma_inv[2]-2.0*gamma_inv[2]; 
  p_fac[5] = 0.5*gamma_inv[5]*temp_sq[28]+0.5*gamma_inv[4]*temp_sq[20]+0.5*gamma_inv[3]*temp_sq[15]+0.5*gamma_inv[2]*temp_sq[12]+0.5*gamma_inv[1]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[5]-1.0*GammaV[0]*temp[5]+gamma[0]*GammaV_sq[3]-1.0*temp[0]*GammaV[3]-1.0*GammaV[1]*temp[2]-1.0*temp[1]*GammaV[2]; 
  p_fac[6] = 0.5*gamma_inv[7]*temp_sq[25]+0.447213595499958*gamma_inv[1]*temp_sq[17]+0.447213595499958*gamma_inv[6]*temp_sq[13]+0.5*gamma_inv[2]*temp_sq[13]-1.0*GammaV[2]*temp[11]+0.5*gamma_inv[3]*temp_sq[8]-1.0*GammaV[3]*temp[7]+0.4472135954999579*gamma_inv[4]*temp_sq[6]+0.5*gamma_inv[0]*temp_sq[6]-1.0*GammaV[0]*temp[6]-1.0*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[1]+0.5*gamma_inv[1]*temp_sq[1]; 
  p_fac[7] = 0.5*gamma_inv[7]*temp_sq[26]+0.447213595499958*gamma_inv[1]*temp_sq[18]+0.447213595499958*gamma_inv[6]*temp_sq[14]+0.5*gamma_inv[2]*temp_sq[14]-1.0*GammaV[1]*temp[11]+0.5*gamma_inv[3]*temp_sq[9]+0.4472135954999579*gamma_inv[4]*temp_sq[7]+0.5*gamma_inv[0]*temp_sq[7]-1.0*GammaV[0]*temp[7]-1.0*GammaV[3]*temp[6]-1.0*GammaV[2]*temp[3]+0.5*gamma_inv[1]*temp_sq[2]+gamma[1]*GammaV_sq[2]; 
  p_fac[8] = 0.447213595499958*gamma_inv[2]*temp_sq[25]+0.5*gamma_inv[6]*temp_sq[17]+0.447213595499958*gamma_inv[7]*temp_sq[13]+0.5*gamma_inv[1]*temp_sq[13]-1.0*GammaV[2]*temp[12]-1.0*GammaV[3]*temp[9]+0.4472135954999579*gamma_inv[5]*temp_sq[8]+0.5*gamma_inv[0]*temp_sq[8]-1.0*GammaV[0]*temp[8]+0.5*gamma_inv[3]*temp_sq[6]-1.0*GammaV[1]*temp[4]+GammaV_sq[1]*gamma[2]+0.5*temp_sq[1]*gamma_inv[2]; 
  p_fac[9] = 0.447213595499958*gamma_inv[2]*temp_sq[26]+0.5*gamma_inv[6]*temp_sq[18]+0.447213595499958*gamma_inv[7]*temp_sq[14]+0.5*gamma_inv[1]*temp_sq[14]-1.0*GammaV[1]*temp[12]+0.4472135954999579*gamma_inv[5]*temp_sq[9]+0.5*gamma_inv[0]*temp_sq[9]-1.0*GammaV[0]*temp[9]-1.0*GammaV[3]*temp[8]+0.5*gamma_inv[3]*temp_sq[7]-1.0*GammaV[2]*temp[4]+GammaV_sq[2]*gamma[2]+0.5*gamma_inv[2]*temp_sq[2]; 
  p_fac[10] = 0.4472135954999579*gamma_inv[3]*temp_sq[24]+0.4472135954999579*gamma_inv[3]*temp_sq[16]+0.4472135954999579*gamma_inv[5]*temp_sq[10]+0.4472135954999579*gamma_inv[4]*temp_sq[10]+0.5*gamma_inv[0]*temp_sq[10]+0.447213595499958*temp_sq[4]*gamma_inv[7]+0.447213595499958*temp_sq[3]*gamma_inv[6]+0.5*gamma_inv[1]*temp_sq[4]+GammaV_sq[0]*gamma[3]+0.5*gamma_inv[2]*temp_sq[3]+0.5*temp_sq[0]*gamma_inv[3]-2.0*gamma_inv[3]; 
  p_fac[11] = 0.5000000000000001*gamma_inv[7]*temp_sq[28]+0.4472135954999579*gamma_inv[1]*temp_sq[20]+0.447213595499958*gamma_inv[6]*temp_sq[15]+0.5*gamma_inv[2]*temp_sq[15]+0.5*gamma_inv[3]*temp_sq[12]+0.4472135954999579*gamma_inv[4]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[11]-1.0*GammaV[0]*temp[11]-1.0*GammaV[1]*temp[7]-1.0*GammaV[2]*temp[6]+0.5*gamma_inv[1]*temp_sq[5]-1.0*GammaV[3]*temp[3]+gamma[1]*GammaV_sq[3]; 
  p_fac[12] = 0.4472135954999579*gamma_inv[2]*temp_sq[28]+0.5000000000000001*gamma_inv[6]*temp_sq[20]+0.447213595499958*gamma_inv[7]*temp_sq[15]+0.5*gamma_inv[1]*temp_sq[15]+0.4472135954999579*gamma_inv[5]*temp_sq[12]+0.5*gamma_inv[0]*temp_sq[12]-1.0*GammaV[0]*temp[12]+0.5*gamma_inv[3]*temp_sq[11]-1.0*GammaV[1]*temp[9]-1.0*GammaV[2]*temp[8]+0.5*gamma_inv[2]*temp_sq[5]-1.0*GammaV[3]*temp[4]+gamma[2]*GammaV_sq[3]; 
  p_fac[13] = 0.447213595499958*gamma_inv[3]*temp_sq[25]+0.447213595499958*gamma_inv[3]*temp_sq[17]+0.4472135954999579*gamma_inv[5]*temp_sq[13]+0.4472135954999579*gamma_inv[4]*temp_sq[13]+0.5*gamma_inv[0]*temp_sq[13]+0.447213595499958*gamma_inv[7]*temp_sq[8]+0.5*gamma_inv[1]*temp_sq[8]+0.447213595499958*gamma_inv[6]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[6]+GammaV_sq[1]*gamma[3]+0.5*temp_sq[1]*gamma_inv[3]; 
  p_fac[14] = 0.447213595499958*gamma_inv[3]*temp_sq[26]+0.447213595499958*gamma_inv[3]*temp_sq[18]+0.4472135954999579*gamma_inv[5]*temp_sq[14]+0.4472135954999579*gamma_inv[4]*temp_sq[14]+0.5*gamma_inv[0]*temp_sq[14]+0.447213595499958*gamma_inv[7]*temp_sq[9]+0.5*gamma_inv[1]*temp_sq[9]+0.447213595499958*gamma_inv[6]*temp_sq[7]+0.5*gamma_inv[2]*temp_sq[7]+GammaV_sq[2]*gamma[3]+0.5*temp_sq[2]*gamma_inv[3]; 
  p_fac[15] = 0.4472135954999579*gamma_inv[3]*temp_sq[28]+0.4472135954999579*gamma_inv[3]*temp_sq[20]+0.4472135954999579*gamma_inv[5]*temp_sq[15]+0.4472135954999579*gamma_inv[4]*temp_sq[15]+0.5*gamma_inv[0]*temp_sq[15]+0.447213595499958*gamma_inv[7]*temp_sq[12]+0.5*gamma_inv[1]*temp_sq[12]+0.447213595499958*gamma_inv[6]*temp_sq[11]+0.5*gamma_inv[2]*temp_sq[11]+0.5*gamma_inv[3]*temp_sq[5]+GammaV_sq[3]*gamma[3]; 
  p_fac[16] = 0.31943828249997*gamma_inv[4]*temp_sq[16]+0.5*gamma_inv[0]*temp_sq[16]+0.4472135954999579*gamma_inv[3]*temp_sq[10]+0.5000000000000001*temp_sq[4]*gamma_inv[6]+GammaV_sq[0]*gamma[4]+0.5*temp_sq[0]*gamma_inv[4]-2.0*gamma_inv[4]+0.4472135954999579*gamma_inv[1]*temp_sq[3]; 
  p_fac[17] = 0.31943828249997*gamma_inv[4]*temp_sq[17]+0.5*gamma_inv[0]*temp_sq[17]+0.447213595499958*gamma_inv[3]*temp_sq[13]+0.5*gamma_inv[6]*temp_sq[8]+0.447213595499958*gamma_inv[1]*temp_sq[6]+1.0*GammaV_sq[1]*gamma[4]+0.5000000000000001*temp_sq[1]*gamma_inv[4]; 
  p_fac[18] = 0.31943828249997*gamma_inv[4]*temp_sq[18]+0.5*gamma_inv[0]*temp_sq[18]+0.447213595499958*gamma_inv[3]*temp_sq[14]+0.5*gamma_inv[6]*temp_sq[9]+0.447213595499958*gamma_inv[1]*temp_sq[7]+1.0*GammaV_sq[2]*gamma[4]+0.5000000000000001*temp_sq[2]*gamma_inv[4]; 
  p_fac[19] = 0.4472135954999579*gamma_inv[6]*temp_sq[24]+0.31943828249997*gamma_inv[6]*temp_sq[16]+0.5000000000000001*gamma_inv[2]*temp_sq[16]+0.4*gamma_inv[7]*temp_sq[10]+0.447213595499958*gamma_inv[1]*temp_sq[10]+GammaV_sq[0]*gamma[6]+0.5*temp_sq[0]*gamma_inv[6]-2.0*gamma_inv[6]+0.5000000000000001*gamma_inv[4]*temp_sq[4]+0.447213595499958*gamma_inv[3]*temp_sq[3]; 
  p_fac[20] = 0.31943828249997*gamma_inv[4]*temp_sq[20]+0.5*gamma_inv[0]*temp_sq[20]+0.4472135954999579*gamma_inv[3]*temp_sq[15]+0.5000000000000001*gamma_inv[6]*temp_sq[12]+0.4472135954999579*gamma_inv[1]*temp_sq[11]+0.5*gamma_inv[4]*temp_sq[5]+GammaV_sq[3]*gamma[4]; 
  p_fac[21] = 0.4472135954999579*gamma_inv[6]*temp_sq[25]+0.31943828249997*gamma_inv[6]*temp_sq[17]+0.5000000000000001*gamma_inv[2]*temp_sq[17]+0.4*gamma_inv[7]*temp_sq[13]+0.4472135954999579*gamma_inv[1]*temp_sq[13]+0.5*gamma_inv[4]*temp_sq[8]+1.0*GammaV_sq[1]*gamma[6]+0.4472135954999579*gamma_inv[3]*temp_sq[6]+0.5000000000000001*temp_sq[1]*gamma_inv[6]; 
  p_fac[22] = 0.4472135954999579*gamma_inv[6]*temp_sq[26]+0.31943828249997*gamma_inv[6]*temp_sq[18]+0.5000000000000001*gamma_inv[2]*temp_sq[18]+0.4*gamma_inv[7]*temp_sq[14]+0.4472135954999579*gamma_inv[1]*temp_sq[14]+0.5*gamma_inv[4]*temp_sq[9]+0.4472135954999579*gamma_inv[3]*temp_sq[7]+1.0*GammaV_sq[2]*gamma[6]+0.5000000000000001*temp_sq[2]*gamma_inv[6]; 
  p_fac[23] = 0.4472135954999579*gamma_inv[6]*temp_sq[28]+0.31943828249997*gamma_inv[6]*temp_sq[20]+0.5000000000000001*gamma_inv[2]*temp_sq[20]+0.4*gamma_inv[7]*temp_sq[15]+0.447213595499958*gamma_inv[1]*temp_sq[15]+0.5000000000000001*gamma_inv[4]*temp_sq[12]+0.447213595499958*gamma_inv[3]*temp_sq[11]+GammaV_sq[3]*gamma[6]+0.5*temp_sq[5]*gamma_inv[6]; 
  p_fac[24] = 0.31943828249997*gamma_inv[5]*temp_sq[24]+0.5*gamma_inv[0]*temp_sq[24]+0.4472135954999579*gamma_inv[3]*temp_sq[10]+0.5000000000000001*temp_sq[3]*gamma_inv[7]+GammaV_sq[0]*gamma[5]+0.5*temp_sq[0]*gamma_inv[5]-2.0*gamma_inv[5]+0.4472135954999579*gamma_inv[2]*temp_sq[4]; 
  p_fac[25] = 0.31943828249997*gamma_inv[5]*temp_sq[25]+0.5*gamma_inv[0]*temp_sq[25]+0.447213595499958*gamma_inv[3]*temp_sq[13]+0.447213595499958*gamma_inv[2]*temp_sq[8]+0.5*temp_sq[6]*gamma_inv[7]+1.0*GammaV_sq[1]*gamma[5]+0.5000000000000001*temp_sq[1]*gamma_inv[5]; 
  p_fac[26] = 0.31943828249997*gamma_inv[5]*temp_sq[26]+0.5*gamma_inv[0]*temp_sq[26]+0.447213595499958*gamma_inv[3]*temp_sq[14]+0.447213595499958*gamma_inv[2]*temp_sq[9]+0.5*gamma_inv[7]*temp_sq[7]+1.0*GammaV_sq[2]*gamma[5]+0.5000000000000001*temp_sq[2]*gamma_inv[5]; 
  p_fac[27] = 0.31943828249997*gamma_inv[7]*temp_sq[24]+0.5000000000000001*gamma_inv[1]*temp_sq[24]+0.4472135954999579*gamma_inv[7]*temp_sq[16]+0.4*gamma_inv[6]*temp_sq[10]+0.447213595499958*gamma_inv[2]*temp_sq[10]+GammaV_sq[0]*gamma[7]+0.5*temp_sq[0]*gamma_inv[7]-2.0*gamma_inv[7]+0.5000000000000001*temp_sq[3]*gamma_inv[5]+0.447213595499958*gamma_inv[3]*temp_sq[4]; 
  p_fac[28] = 0.31943828249997*gamma_inv[5]*temp_sq[28]+0.5*gamma_inv[0]*temp_sq[28]+0.4472135954999579*gamma_inv[3]*temp_sq[15]+0.4472135954999579*gamma_inv[2]*temp_sq[12]+0.5000000000000001*gamma_inv[7]*temp_sq[11]+GammaV_sq[3]*gamma[5]+0.5*gamma_inv[5]*temp_sq[5]; 
  p_fac[29] = 0.31943828249997*gamma_inv[7]*temp_sq[25]+0.5000000000000001*gamma_inv[1]*temp_sq[25]+0.4472135954999579*gamma_inv[7]*temp_sq[17]+0.4*gamma_inv[6]*temp_sq[13]+0.4472135954999579*gamma_inv[2]*temp_sq[13]+0.4472135954999579*gamma_inv[3]*temp_sq[8]+1.0*GammaV_sq[1]*gamma[7]+0.5000000000000001*temp_sq[1]*gamma_inv[7]+0.5*gamma_inv[5]*temp_sq[6]; 
  p_fac[30] = 0.31943828249997*gamma_inv[7]*temp_sq[26]+0.5000000000000001*gamma_inv[1]*temp_sq[26]+0.4472135954999579*gamma_inv[7]*temp_sq[18]+0.4*gamma_inv[6]*temp_sq[14]+0.4472135954999579*gamma_inv[2]*temp_sq[14]+0.4472135954999579*gamma_inv[3]*temp_sq[9]+1.0*GammaV_sq[2]*gamma[7]+0.5*gamma_inv[5]*temp_sq[7]+0.5000000000000001*temp_sq[2]*gamma_inv[7]; 
  p_fac[31] = 0.31943828249997*gamma_inv[7]*temp_sq[28]+0.5000000000000001*gamma_inv[1]*temp_sq[28]+0.4472135954999579*gamma_inv[7]*temp_sq[20]+0.4*gamma_inv[6]*temp_sq[15]+0.447213595499958*gamma_inv[2]*temp_sq[15]+0.447213595499958*gamma_inv[3]*temp_sq[12]+0.5000000000000001*gamma_inv[5]*temp_sq[11]+GammaV_sq[3]*gamma[7]+0.5*temp_sq[5]*gamma_inv[7]; 

  sr_pressure[0] += (0.25*f[31]*p_fac[31]+0.25*f[30]*p_fac[30]+0.25*f[29]*p_fac[29]+0.25*f[28]*p_fac[28]+0.25*f[27]*p_fac[27]+0.25*f[26]*p_fac[26]+0.25*f[25]*p_fac[25]+0.25*f[24]*p_fac[24]+0.25*f[23]*p_fac[23]+0.25*f[22]*p_fac[22]+0.25*f[21]*p_fac[21]+0.25*f[20]*p_fac[20]+0.25*f[19]*p_fac[19]+0.25*f[18]*p_fac[18]+0.25*f[17]*p_fac[17]+0.25*f[16]*p_fac[16]+0.25*f[15]*p_fac[15]+0.25*f[14]*p_fac[14]+0.25*f[13]*p_fac[13]+0.25*f[12]*p_fac[12]+0.25*f[11]*p_fac[11]+0.25*f[10]*p_fac[10]+0.25*f[9]*p_fac[9]+0.25*f[8]*p_fac[8]+0.25*f[7]*p_fac[7]+0.25*f[6]*p_fac[6]+0.25*f[5]*p_fac[5]+0.25*f[4]*p_fac[4]+0.25*f[3]*p_fac[3]+0.25*f[2]*p_fac[2]+0.25*f[1]*p_fac[1]+0.25*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.2500000000000001*f[30]*p_fac[31]+0.2500000000000001*p_fac[30]*f[31]+0.2500000000000001*f[27]*p_fac[29]+0.2500000000000001*p_fac[27]*f[29]+0.2500000000000001*f[26]*p_fac[28]+0.2500000000000001*p_fac[26]*f[28]+0.2500000000000001*f[24]*p_fac[25]+0.2500000000000001*p_fac[24]*f[25]+0.2500000000000001*f[22]*p_fac[23]+0.2500000000000001*p_fac[22]*f[23]+0.2500000000000001*f[19]*p_fac[21]+0.2500000000000001*p_fac[19]*f[21]+0.2500000000000001*f[18]*p_fac[20]+0.2500000000000001*p_fac[18]*f[20]+0.2500000000000001*f[16]*p_fac[17]+0.2500000000000001*p_fac[16]*f[17]+0.25*f[14]*p_fac[15]+0.25*p_fac[14]*f[15]+0.25*f[10]*p_fac[13]+0.25*p_fac[10]*f[13]+0.25*f[9]*p_fac[12]+0.25*p_fac[9]*f[12]+0.25*f[7]*p_fac[11]+0.25*p_fac[7]*f[11]+0.25*f[4]*p_fac[8]+0.25*p_fac[4]*f[8]+0.25*f[3]*p_fac[6]+0.25*p_fac[3]*f[6]+0.25*f[2]*p_fac[5]+0.25*p_fac[2]*f[5]+0.25*f[0]*p_fac[1]+0.25*p_fac[0]*f[1])*volFact; 
  sr_pressure[2] += (0.2500000000000001*f[29]*p_fac[31]+0.2500000000000001*p_fac[29]*f[31]+0.2500000000000001*f[27]*p_fac[30]+0.2500000000000001*p_fac[27]*f[30]+0.2500000000000001*f[25]*p_fac[28]+0.2500000000000001*p_fac[25]*f[28]+0.2500000000000001*f[24]*p_fac[26]+0.2500000000000001*p_fac[24]*f[26]+0.2500000000000001*f[21]*p_fac[23]+0.2500000000000001*p_fac[21]*f[23]+0.2500000000000001*f[19]*p_fac[22]+0.2500000000000001*p_fac[19]*f[22]+0.2500000000000001*f[17]*p_fac[20]+0.2500000000000001*p_fac[17]*f[20]+0.2500000000000001*f[16]*p_fac[18]+0.2500000000000001*p_fac[16]*f[18]+0.25*f[13]*p_fac[15]+0.25*p_fac[13]*f[15]+0.25*f[10]*p_fac[14]+0.25*p_fac[10]*f[14]+0.25*f[8]*p_fac[12]+0.25*p_fac[8]*f[12]+0.25*f[6]*p_fac[11]+0.25*p_fac[6]*f[11]+0.25*f[4]*p_fac[9]+0.25*p_fac[4]*f[9]+0.25*f[3]*p_fac[7]+0.25*p_fac[3]*f[7]+0.25*f[1]*p_fac[5]+0.25*p_fac[1]*f[5]+0.25*f[0]*p_fac[2]+0.25*p_fac[0]*f[2])*volFact; 
  sr_pressure[3] += (0.25*f[27]*p_fac[31]+0.25*p_fac[27]*f[31]+0.25*f[29]*p_fac[30]+0.25*p_fac[29]*f[30]+0.25*f[24]*p_fac[28]+0.25*p_fac[24]*f[28]+0.25*f[25]*p_fac[26]+0.25*p_fac[25]*f[26]+0.25*f[19]*p_fac[23]+0.25*p_fac[19]*f[23]+0.25*f[21]*p_fac[22]+0.25*p_fac[21]*f[22]+0.25*f[16]*p_fac[20]+0.25*p_fac[16]*f[20]+0.25*f[17]*p_fac[18]+0.25*p_fac[17]*f[18]+0.25*f[10]*p_fac[15]+0.25*p_fac[10]*f[15]+0.25*f[13]*p_fac[14]+0.25*p_fac[13]*f[14]+0.25*f[4]*p_fac[12]+0.25*p_fac[4]*f[12]+0.25*f[3]*p_fac[11]+0.25*p_fac[3]*f[11]+0.25*f[8]*p_fac[9]+0.25*p_fac[8]*f[9]+0.25*f[6]*p_fac[7]+0.25*p_fac[6]*f[7]+0.25*f[0]*p_fac[5]+0.25*p_fac[0]*f[5]+0.25*f[1]*p_fac[2]+0.25*p_fac[1]*f[2])*volFact; 
} 
