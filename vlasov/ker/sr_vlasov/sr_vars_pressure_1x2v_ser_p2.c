#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[1]*dxv[2]/4; 
 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double *V_1 = &u_i[3]; 
  const double *V_1_sq = &u_i_sq[3]; 
 
  double temp[20] = {0.0}; 
  double temp_sq[20] = {0.0}; 
  double p_fac[20] = {0.0}; 
  temp[0] = 2.0*V_1[0]*wx2+2.0*V_0[0]*wx1; 
  temp[1] = 2.0*V_1[1]*wx2+2.0*V_0[1]*wx1; 
  temp[2] = 0.5773502691896258*V_0[0]*dv1; 
  temp[3] = 0.5773502691896258*V_1[0]*dv2; 
  temp[4] = 0.5773502691896258*V_0[1]*dv1; 
  temp[5] = 0.5773502691896258*V_1[1]*dv2; 
  temp[7] = 2.0*V_1[2]*wx2+2.0*V_0[2]*wx1; 
  temp[11] = 0.5773502691896257*V_0[2]*dv1; 
  temp[13] = 0.5773502691896257*V_1[2]*dv2; 

  temp_sq[0] = 2.0*V_1_sq[0]*wx2_sq+2.828427124746191*V_0[2]*V_1[2]*wx1*wx2+2.828427124746191*V_0[1]*V_1[1]*wx1*wx2+2.828427124746191*V_0[0]*V_1[0]*wx1*wx2+2.0*V_0_sq[0]*wx1_sq+0.1666666666666667*V_1_sq[0]*dv2_sq+0.1666666666666667*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 2.0*V_1_sq[1]*wx2_sq+2.529822128134704*V_0[1]*V_1[2]*wx1*wx2+2.529822128134704*V_1[1]*V_0[2]*wx1*wx2+2.828427124746191*V_0[0]*V_1[1]*wx1*wx2+2.828427124746191*V_1[0]*V_0[1]*wx1*wx2+2.0*V_0_sq[1]*wx1_sq+0.1666666666666667*V_1_sq[1]*dv2_sq+0.1666666666666667*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 0.8164965809277261*V_0[2]*V_1[2]*dv1*wx2+0.8164965809277261*V_0[1]*V_1[1]*dv1*wx2+0.8164965809277261*V_0[0]*V_1[0]*dv1*wx2+1.154700538379252*V_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 1.154700538379252*V_1_sq[0]*dv2*wx2+0.8164965809277261*V_0[2]*V_1[2]*dv2*wx1+0.8164965809277261*V_0[1]*V_1[1]*dv2*wx1+0.8164965809277261*V_0[0]*V_1[0]*dv2*wx1; 
  temp_sq[4] = 0.7302967433402218*V_0[1]*V_1[2]*dv1*wx2+0.7302967433402218*V_1[1]*V_0[2]*dv1*wx2+0.8164965809277261*V_0[0]*V_1[1]*dv1*wx2+0.8164965809277261*V_1[0]*V_0[1]*dv1*wx2+1.154700538379252*V_0_sq[1]*dv1*wx1; 
  temp_sq[5] = 1.154700538379252*V_1_sq[1]*dv2*wx2+0.7302967433402218*V_0[1]*V_1[2]*dv2*wx1+0.7302967433402218*V_1[1]*V_0[2]*dv2*wx1+0.8164965809277261*V_0[0]*V_1[1]*dv2*wx1+0.8164965809277261*V_1[0]*V_0[1]*dv2*wx1; 
  temp_sq[6] = 0.2357022603955158*V_0[2]*V_1[2]*dv1*dv2+0.2357022603955158*V_0[1]*V_1[1]*dv1*dv2+0.2357022603955158*V_0[0]*V_1[0]*dv1*dv2; 
  temp_sq[7] = 2.0*V_1_sq[2]*wx2_sq+1.807015805810503*V_0[2]*V_1[2]*wx1*wx2+2.828427124746191*V_0[0]*V_1[2]*wx1*wx2+2.828427124746191*V_1[0]*V_0[2]*wx1*wx2+2.529822128134704*V_0[1]*V_1[1]*wx1*wx2+2.0*V_0_sq[2]*wx1_sq+0.1666666666666667*V_1_sq[2]*dv2_sq+0.1666666666666667*V_0_sq[2]*dv1_sq; 
  temp_sq[8] = 0.149071198499986*V_0_sq[0]*dv1_sq; 
  temp_sq[9] = 0.149071198499986*V_1_sq[0]*dv2_sq; 
  temp_sq[10] = 0.210818510677892*V_0[1]*V_1[2]*dv1*dv2+0.210818510677892*V_1[1]*V_0[2]*dv1*dv2+0.2357022603955158*V_0[0]*V_1[1]*dv1*dv2+0.2357022603955158*V_1[0]*V_0[1]*dv1*dv2; 
  temp_sq[11] = 0.5216405309573011*V_0[2]*V_1[2]*dv1*wx2+0.816496580927726*V_0[0]*V_1[2]*dv1*wx2+0.816496580927726*V_1[0]*V_0[2]*dv1*wx2+0.7302967433402215*V_0[1]*V_1[1]*dv1*wx2+1.154700538379251*V_0_sq[2]*dv1*wx1; 
  temp_sq[12] = 0.149071198499986*V_0_sq[1]*dv1_sq; 
  temp_sq[13] = 1.154700538379251*V_1_sq[2]*dv2*wx2+0.5216405309573011*V_0[2]*V_1[2]*dv2*wx1+0.816496580927726*V_0[0]*V_1[2]*dv2*wx1+0.816496580927726*V_1[0]*V_0[2]*dv2*wx1+0.7302967433402215*V_0[1]*V_1[1]*dv2*wx1; 
  temp_sq[15] = 0.149071198499986*V_1_sq[1]*dv2_sq; 
  temp_sq[17] = 0.1505846504842085*V_0[2]*V_1[2]*dv1*dv2+0.2357022603955158*V_0[0]*V_1[2]*dv1*dv2+0.2357022603955158*V_1[0]*V_0[2]*dv1*dv2+0.210818510677892*V_0[1]*V_1[1]*dv1*dv2; 

  p_fac[0] = 0.5*gamma_inv[5]*temp_sq[9]+0.5*gamma_inv[4]*temp_sq[8]-1.414213562373095*GammaV[2]*temp[7]+0.5*gamma_inv[3]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[3]+0.5*gamma_inv[1]*temp_sq[2]-1.414213562373095*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.5*gamma_inv[0]*temp_sq[0]-1.414213562373095*GammaV[0]*temp[0]-1.414213562373095*gamma_inv[0]; 
  p_fac[1] = 0.5000000000000001*gamma_inv[5]*temp_sq[15]+0.5000000000000001*gamma_inv[4]*temp_sq[12]+0.5*gamma_inv[3]*temp_sq[10]-1.264911064067352*GammaV[1]*temp[7]+0.5*gamma_inv[2]*temp_sq[5]+0.5*gamma_inv[1]*temp_sq[4]-1.264911064067352*temp[1]*GammaV[2]+0.5*gamma_inv[0]*temp_sq[1]-1.414213562373095*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.414213562373095*temp[0]*GammaV[1]; 
  p_fac[2] = (-1.414213562373095*GammaV[2]*temp[11])+0.5000000000000001*gamma_inv[7]*temp_sq[9]+0.4472135954999579*gamma_inv[1]*temp_sq[8]+0.447213595499958*gamma_inv[6]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[6]-1.414213562373095*GammaV[1]*temp[4]+0.4472135954999579*temp_sq[2]*gamma_inv[4]+0.5*gamma_inv[3]*temp_sq[3]+0.5*gamma_inv[0]*temp_sq[2]-1.414213562373095*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.5*temp_sq[0]*gamma_inv[1]-1.414213562373095*gamma_inv[1]; 
  p_fac[3] = (-1.414213562373095*GammaV[2]*temp[13])+0.4472135954999579*gamma_inv[2]*temp_sq[9]+0.5000000000000001*gamma_inv[6]*temp_sq[8]+0.447213595499958*temp_sq[6]*gamma_inv[7]+0.5*gamma_inv[1]*temp_sq[6]-1.414213562373095*GammaV[1]*temp[5]+0.4472135954999579*temp_sq[3]*gamma_inv[5]+0.5*gamma_inv[0]*temp_sq[3]-1.414213562373095*GammaV[0]*temp[3]+0.5*temp_sq[2]*gamma_inv[3]+GammaV_sq[0]*gamma[2]+0.5*temp_sq[0]*gamma_inv[2]-1.414213562373095*gamma_inv[2]; 
  p_fac[4] = 0.5*gamma_inv[7]*temp_sq[15]+0.447213595499958*gamma_inv[1]*temp_sq[12]-1.264911064067352*GammaV[1]*temp[11]+0.447213595499958*gamma_inv[6]*temp_sq[10]+0.5*gamma_inv[2]*temp_sq[10]+0.5*gamma_inv[3]*temp_sq[5]+0.4472135954999579*gamma_inv[4]*temp_sq[4]+0.5*gamma_inv[0]*temp_sq[4]-1.264911064067352*GammaV[2]*temp[4]-1.414213562373095*GammaV[0]*temp[4]-1.414213562373095*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.5*gamma_inv[1]*temp_sq[1]; 
  p_fac[5] = 0.447213595499958*gamma_inv[2]*temp_sq[15]-1.264911064067352*GammaV[1]*temp[13]+0.5*gamma_inv[6]*temp_sq[12]+0.447213595499958*gamma_inv[7]*temp_sq[10]+0.5*gamma_inv[1]*temp_sq[10]+0.4472135954999579*gamma_inv[5]*temp_sq[5]+0.5*gamma_inv[0]*temp_sq[5]-1.264911064067352*GammaV[2]*temp[5]-1.414213562373095*GammaV[0]*temp[5]+0.5*gamma_inv[3]*temp_sq[4]-1.414213562373095*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[2]+0.5*temp_sq[1]*gamma_inv[2]; 
  p_fac[6] = 0.4472135954999579*gamma_inv[3]*temp_sq[9]+0.4472135954999579*gamma_inv[3]*temp_sq[8]+0.447213595499958*temp_sq[3]*gamma_inv[7]+0.4472135954999579*gamma_inv[5]*temp_sq[6]+0.4472135954999579*gamma_inv[4]*temp_sq[6]+0.5*gamma_inv[0]*temp_sq[6]+0.447213595499958*temp_sq[2]*gamma_inv[6]+GammaV_sq[0]*gamma[3]+0.5*gamma_inv[1]*temp_sq[3]+0.5*temp_sq[0]*gamma_inv[3]-1.414213562373095*gamma_inv[3]+0.5*gamma_inv[2]*temp_sq[2]; 
  p_fac[7] = 0.5*gamma_inv[3]*temp_sq[17]+0.5000000000000001*gamma_inv[2]*temp_sq[13]+0.5000000000000001*gamma_inv[1]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[7]-0.9035079029052515*GammaV[2]*temp[7]-1.414213562373095*GammaV[0]*temp[7]+gamma[0]*GammaV_sq[2]-1.414213562373095*temp[0]*GammaV[2]-1.264911064067352*GammaV[1]*temp[1]; 
  p_fac[8] = 0.31943828249997*gamma_inv[4]*temp_sq[8]+0.5*gamma_inv[0]*temp_sq[8]+0.4472135954999579*gamma_inv[3]*temp_sq[6]+0.5000000000000001*temp_sq[3]*gamma_inv[6]+GammaV_sq[0]*gamma[4]+0.5*temp_sq[0]*gamma_inv[4]-1.414213562373095*gamma_inv[4]+0.4472135954999579*gamma_inv[1]*temp_sq[2]; 
  p_fac[9] = 0.31943828249997*gamma_inv[5]*temp_sq[9]+0.5*gamma_inv[0]*temp_sq[9]+0.5000000000000001*temp_sq[2]*gamma_inv[7]+0.4472135954999579*gamma_inv[3]*temp_sq[6]+GammaV_sq[0]*gamma[5]+0.5*temp_sq[0]*gamma_inv[5]-1.414213562373095*gamma_inv[5]+0.4472135954999579*gamma_inv[2]*temp_sq[3]; 
  p_fac[10] = 0.447213595499958*gamma_inv[3]*temp_sq[15]+0.447213595499958*gamma_inv[3]*temp_sq[12]+0.4472135954999579*gamma_inv[5]*temp_sq[10]+0.4472135954999579*gamma_inv[4]*temp_sq[10]+0.5*gamma_inv[0]*temp_sq[10]+0.447213595499958*temp_sq[5]*gamma_inv[7]+0.447213595499958*temp_sq[4]*gamma_inv[6]+0.5*gamma_inv[1]*temp_sq[5]+0.5*gamma_inv[2]*temp_sq[4]+GammaV_sq[1]*gamma[3]+0.5*temp_sq[1]*gamma_inv[3]; 
  p_fac[11] = 0.4472135954999579*gamma_inv[6]*temp_sq[17]+0.5000000000000001*gamma_inv[2]*temp_sq[17]+0.5*gamma_inv[3]*temp_sq[13]+0.4472135954999579*gamma_inv[4]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[11]-0.9035079029052515*GammaV[2]*temp[11]-1.414213562373095*GammaV[0]*temp[11]+0.5000000000000001*gamma_inv[1]*temp_sq[7]-1.264911064067352*GammaV[1]*temp[4]-1.414213562373095*GammaV[2]*temp[2]+1.0*gamma[1]*GammaV_sq[2]; 
  p_fac[12] = 0.31943828249997*gamma_inv[4]*temp_sq[12]+0.5*gamma_inv[0]*temp_sq[12]+0.447213595499958*gamma_inv[3]*temp_sq[10]+0.5*temp_sq[5]*gamma_inv[6]+1.0*GammaV_sq[1]*gamma[4]+0.447213595499958*gamma_inv[1]*temp_sq[4]+0.5000000000000001*temp_sq[1]*gamma_inv[4]; 
  p_fac[13] = 0.4472135954999579*gamma_inv[7]*temp_sq[17]+0.5000000000000001*gamma_inv[1]*temp_sq[17]+0.4472135954999579*gamma_inv[5]*temp_sq[13]+0.5*gamma_inv[0]*temp_sq[13]-0.9035079029052515*GammaV[2]*temp[13]-1.414213562373095*GammaV[0]*temp[13]+0.5*gamma_inv[3]*temp_sq[11]+0.5000000000000001*gamma_inv[2]*temp_sq[7]-1.264911064067352*GammaV[1]*temp[5]-1.414213562373095*GammaV[2]*temp[3]+1.0*GammaV_sq[2]*gamma[2]; 
  p_fac[14] = 0.4472135954999579*gamma_inv[6]*temp_sq[9]+0.31943828249997*gamma_inv[6]*temp_sq[8]+0.5000000000000001*gamma_inv[2]*temp_sq[8]+0.4*temp_sq[6]*gamma_inv[7]+GammaV_sq[0]*gamma[6]+0.447213595499958*gamma_inv[1]*temp_sq[6]+0.5*temp_sq[0]*gamma_inv[6]-1.414213562373095*gamma_inv[6]+0.5000000000000001*temp_sq[3]*gamma_inv[4]+0.447213595499958*temp_sq[2]*gamma_inv[3]; 
  p_fac[15] = 0.31943828249997*gamma_inv[5]*temp_sq[15]+0.5*gamma_inv[0]*temp_sq[15]+0.447213595499958*gamma_inv[3]*temp_sq[10]+0.5*temp_sq[4]*gamma_inv[7]+1.0*GammaV_sq[1]*gamma[5]+0.447213595499958*gamma_inv[2]*temp_sq[5]+0.5000000000000001*temp_sq[1]*gamma_inv[5]; 
  p_fac[16] = 0.31943828249997*gamma_inv[7]*temp_sq[9]+0.5000000000000001*gamma_inv[1]*temp_sq[9]+0.4472135954999579*gamma_inv[7]*temp_sq[8]+GammaV_sq[0]*gamma[7]+0.5*temp_sq[0]*gamma_inv[7]-1.414213562373095*gamma_inv[7]+0.4*gamma_inv[6]*temp_sq[6]+0.447213595499958*gamma_inv[2]*temp_sq[6]+0.5000000000000001*temp_sq[2]*gamma_inv[5]+0.447213595499958*gamma_inv[3]*temp_sq[3]; 
  p_fac[17] = 0.4472135954999579*gamma_inv[5]*temp_sq[17]+0.4472135954999579*gamma_inv[4]*temp_sq[17]+0.5*gamma_inv[0]*temp_sq[17]+0.4472135954999579*gamma_inv[7]*temp_sq[13]+0.5000000000000001*gamma_inv[1]*temp_sq[13]+0.4472135954999579*gamma_inv[6]*temp_sq[11]+0.5000000000000001*gamma_inv[2]*temp_sq[11]+0.5*gamma_inv[3]*temp_sq[7]+GammaV_sq[2]*gamma[3]; 
  p_fac[18] = 0.4472135954999579*gamma_inv[6]*temp_sq[15]+0.31943828249997*gamma_inv[6]*temp_sq[12]+0.5000000000000001*gamma_inv[2]*temp_sq[12]+0.4*gamma_inv[7]*temp_sq[10]+0.4472135954999579*gamma_inv[1]*temp_sq[10]+1.0*GammaV_sq[1]*gamma[6]+0.5000000000000001*temp_sq[1]*gamma_inv[6]+0.5*gamma_inv[4]*temp_sq[5]+0.4472135954999579*gamma_inv[3]*temp_sq[4]; 
  p_fac[19] = 0.31943828249997*gamma_inv[7]*temp_sq[15]+0.5000000000000001*gamma_inv[1]*temp_sq[15]+0.4472135954999579*gamma_inv[7]*temp_sq[12]+0.4*gamma_inv[6]*temp_sq[10]+0.4472135954999579*gamma_inv[2]*temp_sq[10]+1.0*GammaV_sq[1]*gamma[7]+0.5000000000000001*temp_sq[1]*gamma_inv[7]+0.4472135954999579*gamma_inv[3]*temp_sq[5]+0.5*temp_sq[4]*gamma_inv[5]; 

  sr_pressure[0] += (0.3535533905932737*f[19]*p_fac[19]+0.3535533905932737*f[18]*p_fac[18]+0.3535533905932737*f[17]*p_fac[17]+0.3535533905932737*f[16]*p_fac[16]+0.3535533905932737*f[15]*p_fac[15]+0.3535533905932737*f[14]*p_fac[14]+0.3535533905932737*f[13]*p_fac[13]+0.3535533905932737*f[12]*p_fac[12]+0.3535533905932737*f[11]*p_fac[11]+0.3535533905932737*f[10]*p_fac[10]+0.3535533905932737*f[9]*p_fac[9]+0.3535533905932737*f[8]*p_fac[8]+0.3535533905932737*f[7]*p_fac[7]+0.3535533905932737*f[6]*p_fac[6]+0.3535533905932737*f[5]*p_fac[5]+0.3535533905932737*f[4]*p_fac[4]+0.3535533905932737*f[3]*p_fac[3]+0.3535533905932737*f[2]*p_fac[2]+0.3535533905932737*f[1]*p_fac[1]+0.3535533905932737*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.3535533905932737*f[16]*p_fac[19]+0.3535533905932737*p_fac[16]*f[19]+0.3535533905932737*f[14]*p_fac[18]+0.3535533905932737*p_fac[14]*f[18]+0.3162277660168379*f[10]*p_fac[17]+0.3162277660168379*p_fac[10]*f[17]+0.3535533905932737*f[9]*p_fac[15]+0.3535533905932737*p_fac[9]*f[15]+0.3162277660168379*f[5]*p_fac[13]+0.3162277660168379*p_fac[5]*f[13]+0.3535533905932737*f[8]*p_fac[12]+0.3535533905932737*p_fac[8]*f[12]+0.3162277660168379*f[4]*p_fac[11]+0.3162277660168379*p_fac[4]*f[11]+0.3535533905932737*f[6]*p_fac[10]+0.3535533905932737*p_fac[6]*f[10]+0.3162277660168379*f[1]*p_fac[7]+0.3162277660168379*p_fac[1]*f[7]+0.3535533905932737*f[3]*p_fac[5]+0.3535533905932737*p_fac[3]*f[5]+0.3535533905932737*f[2]*p_fac[4]+0.3535533905932737*p_fac[2]*f[4]+0.3535533905932737*f[0]*p_fac[1]+0.3535533905932737*p_fac[0]*f[1])*volFact; 
  sr_pressure[2] += (0.3162277660168379*f[19]*p_fac[19]+0.3162277660168379*f[18]*p_fac[18]+0.2258769757263128*f[17]*p_fac[17]+0.3535533905932737*f[6]*p_fac[17]+0.3535533905932737*p_fac[6]*f[17]+0.3162277660168379*f[15]*p_fac[15]+0.2258769757263128*f[13]*p_fac[13]+0.3535533905932737*f[3]*p_fac[13]+0.3535533905932737*p_fac[3]*f[13]+0.3162277660168379*f[12]*p_fac[12]+0.2258769757263128*f[11]*p_fac[11]+0.3535533905932737*f[2]*p_fac[11]+0.3535533905932737*p_fac[2]*f[11]+0.3162277660168379*f[10]*p_fac[10]+0.2258769757263128*f[7]*p_fac[7]+0.3535533905932737*f[0]*p_fac[7]+0.3535533905932737*p_fac[0]*f[7]+0.3162277660168379*f[5]*p_fac[5]+0.3162277660168379*f[4]*p_fac[4]+0.3162277660168379*f[1]*p_fac[1])*volFact; 
} 
