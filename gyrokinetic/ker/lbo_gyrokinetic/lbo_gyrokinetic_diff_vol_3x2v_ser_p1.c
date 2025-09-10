#include <gkyl_lbo_gyrokinetic_kernels.h> 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // dxv[5]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[8];

  double rdvVmapPrimeSq4[2]; 
  rdvVmapPrimeSq4[0] = 4.0/(dxv[3]*dxv[3]*vmap_prime[0]*vmap_prime[0]); 
  rdvVmapPrimeSq4[1] = 4.0/(dxv[4]*dxv[4]*vmap_prime[1]*vmap_prime[1]); 

  // Expand nuVtSqSum/vpar'^2 in conf basis.
  double facDiffVpar[8] = {0.};
  facDiffVpar[0] = nuVtSqSum[0]*rdvVmapPrimeSq4[0]; 
  facDiffVpar[1] = rdvVmapPrimeSq4[0]*nuVtSqSum[1]; 
  facDiffVpar[2] = rdvVmapPrimeSq4[0]*nuVtSqSum[2]; 
  facDiffVpar[3] = rdvVmapPrimeSq4[0]*nuVtSqSum[3]; 
  facDiffVpar[4] = rdvVmapPrimeSq4[0]*nuVtSqSum[4]; 
  facDiffVpar[5] = rdvVmapPrimeSq4[0]*nuVtSqSum[5]; 
  facDiffVpar[6] = rdvVmapPrimeSq4[0]*nuVtSqSum[6]; 
  facDiffVpar[7] = rdvVmapPrimeSq4[0]*nuVtSqSum[7]; 

  // Expand 2*m*nuVtSqSum/bmag/mu'^2 in conf basis.
  double facDiffMu[8] = {0.};
  facDiffMu[0] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[5]*nuVtSqSum[5]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiffMu[1] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[3]*nuVtSqSum[5]+nuVtSqSum[3]*bmag_inv[5]+bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 
  facDiffMu[2] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[5]*nuVtSqSum[7]+bmag_inv[3]*nuVtSqSum[6]+bmag_inv[1]*nuVtSqSum[4]+bmag_inv[0]*nuVtSqSum[2])*m_; 
  facDiffMu[3] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[1]*nuVtSqSum[5]+nuVtSqSum[1]*bmag_inv[5]+bmag_inv[0]*nuVtSqSum[3]+nuVtSqSum[0]*bmag_inv[3])*m_; 
  facDiffMu[4] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[3]*nuVtSqSum[7]+bmag_inv[5]*nuVtSqSum[6]+bmag_inv[0]*nuVtSqSum[4]+bmag_inv[1]*nuVtSqSum[2])*m_; 
  facDiffMu[5] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[0]*nuVtSqSum[5]+nuVtSqSum[0]*bmag_inv[5]+bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3])*m_; 
  facDiffMu[6] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[1]*nuVtSqSum[7]+bmag_inv[0]*nuVtSqSum[6]+nuVtSqSum[4]*bmag_inv[5]+nuVtSqSum[2]*bmag_inv[3])*m_; 
  facDiffMu[7] = 0.7071067811865475*rdvVmapPrimeSq4[1]*(bmag_inv[0]*nuVtSqSum[7]+bmag_inv[1]*nuVtSqSum[6]+nuVtSqSum[2]*bmag_inv[5]+bmag_inv[3]*nuVtSqSum[4])*m_; 

  out[5] += 0.75*vmap[3]*facDiffMu[7]*fin[16]+0.75*vmap[3]*facDiffMu[6]*fin[8]+0.75*vmap[3]*facDiffMu[5]*fin[7]+0.75*vmap[3]*facDiffMu[4]*fin[6]+0.75*facDiffMu[3]*fin[3]*vmap[3]+0.75*facDiffMu[2]*fin[2]*vmap[3]+0.75*facDiffMu[1]*fin[1]*vmap[3]+0.75*facDiffMu[0]*fin[0]*vmap[3]; 
  out[12] += 0.75*vmap[3]*facDiffMu[6]*fin[16]+0.75*vmap[3]*facDiffMu[7]*fin[8]+0.75*facDiffMu[3]*vmap[3]*fin[7]+0.75*facDiffMu[2]*vmap[3]*fin[6]+0.75*fin[3]*vmap[3]*facDiffMu[5]+0.75*fin[2]*vmap[3]*facDiffMu[4]+0.75*facDiffMu[0]*fin[1]*vmap[3]+0.75*fin[0]*facDiffMu[1]*vmap[3]; 
  out[13] += 0.75*vmap[3]*facDiffMu[5]*fin[16]+0.75*facDiffMu[3]*vmap[3]*fin[8]+0.75*vmap[3]*facDiffMu[7]*fin[7]+0.75*facDiffMu[1]*vmap[3]*fin[6]+0.75*fin[3]*vmap[3]*facDiffMu[6]+0.75*fin[1]*vmap[3]*facDiffMu[4]+0.75*facDiffMu[0]*fin[2]*vmap[3]+0.75*fin[0]*facDiffMu[2]*vmap[3]; 
  out[14] += 0.75*vmap[3]*facDiffMu[4]*fin[16]+0.75*facDiffMu[2]*vmap[3]*fin[8]+0.75*facDiffMu[1]*vmap[3]*fin[7]+0.75*vmap[3]*fin[6]*facDiffMu[7]+0.75*fin[2]*vmap[3]*facDiffMu[6]+0.75*fin[1]*vmap[3]*facDiffMu[5]+0.75*facDiffMu[0]*fin[3]*vmap[3]+0.75*fin[0]*facDiffMu[3]*vmap[3]; 
  out[15] += 0.75*vmap[3]*facDiffMu[7]*fin[26]+0.75*vmap[3]*facDiffMu[6]*fin[19]+0.75*vmap[3]*facDiffMu[5]*fin[18]+0.75*vmap[3]*facDiffMu[4]*fin[17]+0.75*facDiffMu[3]*vmap[3]*fin[11]+0.75*facDiffMu[2]*vmap[3]*fin[10]+0.75*facDiffMu[1]*vmap[3]*fin[9]+0.75*facDiffMu[0]*vmap[3]*fin[4]; 
  out[20] += 0.75*facDiffMu[3]*vmap[3]*fin[16]+0.75*vmap[3]*facDiffMu[5]*fin[8]+0.75*vmap[3]*facDiffMu[6]*fin[7]+0.75*fin[3]*vmap[3]*facDiffMu[7]+0.75*facDiffMu[0]*vmap[3]*fin[6]+0.75*fin[0]*vmap[3]*facDiffMu[4]+0.75*facDiffMu[1]*fin[2]*vmap[3]+0.75*fin[1]*facDiffMu[2]*vmap[3]; 
  out[21] += 0.75*facDiffMu[2]*vmap[3]*fin[16]+0.75*vmap[3]*facDiffMu[4]*fin[8]+0.75*facDiffMu[0]*vmap[3]*fin[7]+0.75*fin[2]*vmap[3]*facDiffMu[7]+0.75*vmap[3]*facDiffMu[6]*fin[6]+0.75*fin[0]*vmap[3]*facDiffMu[5]+0.75*facDiffMu[1]*fin[3]*vmap[3]+0.75*fin[1]*facDiffMu[3]*vmap[3]; 
  out[22] += 0.75*facDiffMu[1]*vmap[3]*fin[16]+0.75*facDiffMu[0]*vmap[3]*fin[8]+0.75*vmap[3]*facDiffMu[4]*fin[7]+0.75*fin[1]*vmap[3]*facDiffMu[7]+0.75*vmap[3]*facDiffMu[5]*fin[6]+0.75*fin[0]*vmap[3]*facDiffMu[6]+0.75*facDiffMu[2]*fin[3]*vmap[3]+0.75*fin[2]*facDiffMu[3]*vmap[3]; 
  out[23] += 0.75*vmap[3]*facDiffMu[6]*fin[26]+0.75*vmap[3]*facDiffMu[7]*fin[19]+0.75*facDiffMu[3]*vmap[3]*fin[18]+0.75*facDiffMu[2]*vmap[3]*fin[17]+0.75*vmap[3]*facDiffMu[5]*fin[11]+0.75*vmap[3]*facDiffMu[4]*fin[10]+0.75*facDiffMu[0]*vmap[3]*fin[9]+0.75*facDiffMu[1]*vmap[3]*fin[4]; 
  out[24] += 0.75*vmap[3]*facDiffMu[5]*fin[26]+0.75*facDiffMu[3]*vmap[3]*fin[19]+0.75*vmap[3]*facDiffMu[7]*fin[18]+0.75*facDiffMu[1]*vmap[3]*fin[17]+0.75*vmap[3]*facDiffMu[6]*fin[11]+0.75*facDiffMu[0]*vmap[3]*fin[10]+0.75*vmap[3]*facDiffMu[4]*fin[9]+0.75*facDiffMu[2]*vmap[3]*fin[4]; 
  out[25] += 0.75*vmap[3]*facDiffMu[4]*fin[26]+0.75*facDiffMu[2]*vmap[3]*fin[19]+0.75*facDiffMu[1]*vmap[3]*fin[18]+0.75*vmap[3]*facDiffMu[7]*fin[17]+0.75*facDiffMu[0]*vmap[3]*fin[11]+0.75*vmap[3]*facDiffMu[6]*fin[10]+0.75*vmap[3]*facDiffMu[5]*fin[9]+0.75*facDiffMu[3]*vmap[3]*fin[4]; 
  out[27] += 0.75*facDiffMu[0]*vmap[3]*fin[16]+0.75*facDiffMu[1]*vmap[3]*fin[8]+0.75*facDiffMu[2]*vmap[3]*fin[7]+0.75*fin[0]*vmap[3]*facDiffMu[7]+0.75*facDiffMu[3]*vmap[3]*fin[6]+0.75*fin[1]*vmap[3]*facDiffMu[6]+0.75*fin[2]*vmap[3]*facDiffMu[5]+0.75*fin[3]*vmap[3]*facDiffMu[4]; 
  out[28] += 0.75*facDiffMu[3]*vmap[3]*fin[26]+0.75*vmap[3]*facDiffMu[5]*fin[19]+0.75*vmap[3]*facDiffMu[6]*fin[18]+0.75*facDiffMu[0]*vmap[3]*fin[17]+0.75*vmap[3]*facDiffMu[7]*fin[11]+0.75*facDiffMu[1]*vmap[3]*fin[10]+0.75*facDiffMu[2]*vmap[3]*fin[9]+0.75*vmap[3]*facDiffMu[4]*fin[4]; 
  out[29] += 0.75*facDiffMu[2]*vmap[3]*fin[26]+0.75*vmap[3]*facDiffMu[4]*fin[19]+0.75*facDiffMu[0]*vmap[3]*fin[18]+0.75*vmap[3]*facDiffMu[6]*fin[17]+0.75*facDiffMu[1]*vmap[3]*fin[11]+0.75*vmap[3]*facDiffMu[7]*fin[10]+0.75*facDiffMu[3]*vmap[3]*fin[9]+0.75*vmap[3]*fin[4]*facDiffMu[5]; 
  out[30] += 0.75*facDiffMu[1]*vmap[3]*fin[26]+0.75*facDiffMu[0]*vmap[3]*fin[19]+0.75*vmap[3]*facDiffMu[4]*fin[18]+0.75*vmap[3]*facDiffMu[5]*fin[17]+0.75*facDiffMu[2]*vmap[3]*fin[11]+0.75*facDiffMu[3]*vmap[3]*fin[10]+0.75*vmap[3]*facDiffMu[7]*fin[9]+0.75*vmap[3]*fin[4]*facDiffMu[6]; 
  out[31] += 0.75*facDiffMu[0]*vmap[3]*fin[26]+0.75*facDiffMu[1]*vmap[3]*fin[19]+0.75*facDiffMu[2]*vmap[3]*fin[18]+0.75*facDiffMu[3]*vmap[3]*fin[17]+0.75*vmap[3]*facDiffMu[4]*fin[11]+0.75*vmap[3]*facDiffMu[5]*fin[10]+0.75*vmap[3]*facDiffMu[6]*fin[9]+0.75*vmap[3]*fin[4]*facDiffMu[7]; 
  out[32] += 2.371708245126284*facDiffVpar[7]*fin[16]+2.371708245126284*facDiffVpar[6]*fin[8]+2.371708245126284*facDiffVpar[5]*fin[7]+2.371708245126284*facDiffVpar[4]*fin[6]+2.371708245126284*facDiffVpar[3]*fin[3]+2.371708245126284*facDiffVpar[2]*fin[2]+2.371708245126284*facDiffVpar[1]*fin[1]+2.371708245126284*facDiffVpar[0]*fin[0]; 
  out[33] += 2.371708245126284*facDiffVpar[6]*fin[16]+2.371708245126284*facDiffVpar[7]*fin[8]+2.371708245126284*facDiffVpar[3]*fin[7]+2.371708245126284*facDiffVpar[2]*fin[6]+2.371708245126284*fin[3]*facDiffVpar[5]+2.371708245126284*fin[2]*facDiffVpar[4]+2.371708245126284*facDiffVpar[0]*fin[1]+2.371708245126284*fin[0]*facDiffVpar[1]; 
  out[34] += 2.371708245126284*facDiffVpar[5]*fin[16]+2.371708245126284*facDiffVpar[3]*fin[8]+2.371708245126284*facDiffVpar[7]*fin[7]+2.371708245126284*facDiffVpar[1]*fin[6]+2.371708245126284*fin[3]*facDiffVpar[6]+2.371708245126284*fin[1]*facDiffVpar[4]+2.371708245126284*facDiffVpar[0]*fin[2]+2.371708245126284*fin[0]*facDiffVpar[2]; 
  out[35] += 2.371708245126284*facDiffVpar[4]*fin[16]+2.371708245126284*facDiffVpar[2]*fin[8]+2.371708245126284*facDiffVpar[1]*fin[7]+2.371708245126284*fin[6]*facDiffVpar[7]+2.371708245126284*fin[2]*facDiffVpar[6]+2.371708245126284*fin[1]*facDiffVpar[5]+2.371708245126284*facDiffVpar[0]*fin[3]+2.371708245126284*fin[0]*facDiffVpar[3]; 
  out[36] += 0.75*vmap[3]*facDiffMu[7]*fin[43]+0.75*vmap[3]*facDiffMu[6]*fin[39]+0.75*vmap[3]*facDiffMu[5]*fin[38]+0.75*vmap[3]*facDiffMu[4]*fin[37]+0.75*facDiffMu[3]*vmap[3]*fin[35]+0.75*facDiffMu[2]*vmap[3]*fin[34]+0.75*facDiffMu[1]*vmap[3]*fin[33]+0.75*facDiffMu[0]*vmap[3]*fin[32]+2.371708245126284*facDiffVpar[7]*fin[27]+2.371708245126284*facDiffVpar[6]*fin[22]+2.371708245126284*facDiffVpar[5]*fin[21]+2.371708245126284*facDiffVpar[4]*fin[20]+2.371708245126284*facDiffVpar[3]*fin[14]+2.371708245126284*facDiffVpar[2]*fin[13]+2.371708245126284*facDiffVpar[1]*fin[12]+2.371708245126284*facDiffVpar[0]*fin[5]; 
  out[37] += 2.371708245126284*facDiffVpar[3]*fin[16]+2.371708245126284*facDiffVpar[5]*fin[8]+2.371708245126284*facDiffVpar[6]*fin[7]+2.371708245126284*fin[3]*facDiffVpar[7]+2.371708245126284*facDiffVpar[0]*fin[6]+2.371708245126284*fin[0]*facDiffVpar[4]+2.371708245126284*facDiffVpar[1]*fin[2]+2.371708245126284*fin[1]*facDiffVpar[2]; 
  out[38] += 2.371708245126284*facDiffVpar[2]*fin[16]+2.371708245126284*facDiffVpar[4]*fin[8]+2.371708245126284*facDiffVpar[0]*fin[7]+2.371708245126284*fin[2]*facDiffVpar[7]+2.371708245126284*facDiffVpar[6]*fin[6]+2.371708245126284*fin[0]*facDiffVpar[5]+2.371708245126284*facDiffVpar[1]*fin[3]+2.371708245126284*fin[1]*facDiffVpar[3]; 
  out[39] += 2.371708245126284*facDiffVpar[1]*fin[16]+2.371708245126284*facDiffVpar[0]*fin[8]+2.371708245126284*facDiffVpar[4]*fin[7]+2.371708245126284*fin[1]*facDiffVpar[7]+2.371708245126284*facDiffVpar[5]*fin[6]+2.371708245126284*fin[0]*facDiffVpar[6]+2.371708245126284*facDiffVpar[2]*fin[3]+2.371708245126284*fin[2]*facDiffVpar[3]; 
  out[40] += 0.75*vmap[3]*facDiffMu[6]*fin[43]+0.75*vmap[3]*facDiffMu[7]*fin[39]+0.75*facDiffMu[3]*vmap[3]*fin[38]+0.75*facDiffMu[2]*vmap[3]*fin[37]+0.75*vmap[3]*facDiffMu[5]*fin[35]+0.75*vmap[3]*facDiffMu[4]*fin[34]+0.75*facDiffMu[0]*vmap[3]*fin[33]+0.75*facDiffMu[1]*vmap[3]*fin[32]+2.371708245126284*facDiffVpar[6]*fin[27]+2.371708245126284*facDiffVpar[7]*fin[22]+2.371708245126284*facDiffVpar[3]*fin[21]+2.371708245126284*facDiffVpar[2]*fin[20]+2.371708245126284*facDiffVpar[5]*fin[14]+2.371708245126284*facDiffVpar[4]*fin[13]+2.371708245126284*facDiffVpar[0]*fin[12]+2.371708245126284*facDiffVpar[1]*fin[5]; 
  out[41] += 0.75*vmap[3]*facDiffMu[5]*fin[43]+0.75*facDiffMu[3]*vmap[3]*fin[39]+0.75*vmap[3]*facDiffMu[7]*fin[38]+0.75*facDiffMu[1]*vmap[3]*fin[37]+0.75*vmap[3]*facDiffMu[6]*fin[35]+0.75*facDiffMu[0]*vmap[3]*fin[34]+0.75*vmap[3]*facDiffMu[4]*fin[33]+0.75*facDiffMu[2]*vmap[3]*fin[32]+2.371708245126284*facDiffVpar[5]*fin[27]+2.371708245126284*facDiffVpar[3]*fin[22]+2.371708245126284*facDiffVpar[7]*fin[21]+2.371708245126284*facDiffVpar[1]*fin[20]+2.371708245126284*facDiffVpar[6]*fin[14]+2.371708245126284*facDiffVpar[0]*fin[13]+2.371708245126284*facDiffVpar[4]*fin[12]+2.371708245126284*facDiffVpar[2]*fin[5]; 
  out[42] += 0.75*vmap[3]*facDiffMu[4]*fin[43]+0.75*facDiffMu[2]*vmap[3]*fin[39]+0.75*facDiffMu[1]*vmap[3]*fin[38]+0.75*vmap[3]*facDiffMu[7]*fin[37]+0.75*facDiffMu[0]*vmap[3]*fin[35]+0.75*vmap[3]*facDiffMu[6]*fin[34]+0.75*vmap[3]*facDiffMu[5]*fin[33]+0.75*facDiffMu[3]*vmap[3]*fin[32]+2.371708245126284*facDiffVpar[4]*fin[27]+2.371708245126284*facDiffVpar[2]*fin[22]+2.371708245126284*facDiffVpar[1]*fin[21]+2.371708245126284*facDiffVpar[7]*fin[20]+2.371708245126284*facDiffVpar[0]*fin[14]+2.371708245126284*facDiffVpar[6]*fin[13]+2.371708245126284*facDiffVpar[5]*fin[12]+2.371708245126284*facDiffVpar[3]*fin[5]; 
  out[43] += 2.371708245126284*facDiffVpar[0]*fin[16]+2.371708245126284*facDiffVpar[1]*fin[8]+2.371708245126284*facDiffVpar[2]*fin[7]+2.371708245126284*fin[0]*facDiffVpar[7]+2.371708245126284*facDiffVpar[3]*fin[6]+2.371708245126284*fin[1]*facDiffVpar[6]+2.371708245126284*fin[2]*facDiffVpar[5]+2.371708245126284*fin[3]*facDiffVpar[4]; 
  out[44] += 0.75*facDiffMu[3]*vmap[3]*fin[43]+0.75*vmap[3]*facDiffMu[5]*fin[39]+0.75*vmap[3]*facDiffMu[6]*fin[38]+0.75*facDiffMu[0]*vmap[3]*fin[37]+0.75*vmap[3]*facDiffMu[7]*fin[35]+0.75*facDiffMu[1]*vmap[3]*fin[34]+0.75*facDiffMu[2]*vmap[3]*fin[33]+0.75*vmap[3]*facDiffMu[4]*fin[32]+2.371708245126284*facDiffVpar[3]*fin[27]+2.371708245126284*facDiffVpar[5]*fin[22]+2.371708245126284*facDiffVpar[6]*fin[21]+2.371708245126284*facDiffVpar[0]*fin[20]+2.371708245126284*facDiffVpar[7]*fin[14]+2.371708245126284*facDiffVpar[1]*fin[13]+2.371708245126284*facDiffVpar[2]*fin[12]+2.371708245126284*facDiffVpar[4]*fin[5]; 
  out[45] += 0.75*facDiffMu[2]*vmap[3]*fin[43]+0.75*vmap[3]*facDiffMu[4]*fin[39]+0.75*facDiffMu[0]*vmap[3]*fin[38]+0.75*vmap[3]*facDiffMu[6]*fin[37]+0.75*facDiffMu[1]*vmap[3]*fin[35]+0.75*vmap[3]*facDiffMu[7]*fin[34]+0.75*facDiffMu[3]*vmap[3]*fin[33]+0.75*vmap[3]*facDiffMu[5]*fin[32]+2.371708245126284*facDiffVpar[2]*fin[27]+2.371708245126284*facDiffVpar[4]*fin[22]+2.371708245126284*facDiffVpar[0]*fin[21]+2.371708245126284*facDiffVpar[6]*fin[20]+2.371708245126284*facDiffVpar[1]*fin[14]+2.371708245126284*facDiffVpar[7]*fin[13]+2.371708245126284*facDiffVpar[3]*fin[12]+2.371708245126284*facDiffVpar[5]*fin[5]; 
  out[46] += 0.75*facDiffMu[1]*vmap[3]*fin[43]+0.75*facDiffMu[0]*vmap[3]*fin[39]+0.75*vmap[3]*facDiffMu[4]*fin[38]+0.75*vmap[3]*facDiffMu[5]*fin[37]+0.75*facDiffMu[2]*vmap[3]*fin[35]+0.75*facDiffMu[3]*vmap[3]*fin[34]+0.75*vmap[3]*facDiffMu[7]*fin[33]+0.75*vmap[3]*facDiffMu[6]*fin[32]+2.371708245126284*facDiffVpar[1]*fin[27]+2.371708245126284*facDiffVpar[0]*fin[22]+2.371708245126284*facDiffVpar[4]*fin[21]+2.371708245126284*facDiffVpar[5]*fin[20]+2.371708245126284*facDiffVpar[2]*fin[14]+2.371708245126284*facDiffVpar[3]*fin[13]+2.371708245126284*facDiffVpar[7]*fin[12]+2.371708245126284*fin[5]*facDiffVpar[6]; 
  out[47] += 0.75*facDiffMu[0]*vmap[3]*fin[43]+0.75*facDiffMu[1]*vmap[3]*fin[39]+0.75*facDiffMu[2]*vmap[3]*fin[38]+0.75*facDiffMu[3]*vmap[3]*fin[37]+0.75*vmap[3]*facDiffMu[4]*fin[35]+0.75*vmap[3]*facDiffMu[5]*fin[34]+0.75*vmap[3]*facDiffMu[6]*fin[33]+0.75*vmap[3]*facDiffMu[7]*fin[32]+2.371708245126284*facDiffVpar[0]*fin[27]+2.371708245126284*facDiffVpar[1]*fin[22]+2.371708245126284*facDiffVpar[2]*fin[21]+2.371708245126284*facDiffVpar[3]*fin[20]+2.371708245126284*facDiffVpar[4]*fin[14]+2.371708245126284*facDiffVpar[5]*fin[13]+2.371708245126284*facDiffVpar[6]*fin[12]+2.371708245126284*fin[5]*facDiffVpar[7]; 

  return fabs(facDiffMu[0]*vmap[2]+3.181980515339463*facDiffVpar[0]); 

} 
