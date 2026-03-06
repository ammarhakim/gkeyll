#include <gkyl_mom_weighted_gyrokinetic_kernels.h> 
GKYL_CU_DH void int_mom_gyrokinetic_M0_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M1_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (vmap[1]*wgt[1]*f[3]+wgt[0]*vmap[1]*f[2]+vmap[0]*f[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2PAR_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.632455532033676*vmap1R2*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M2_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.632455532033676*vmap1R2*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_M0M1M2_WEIGHTX_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.4142135623730951*f[1]*wgt[1]+1.4142135623730951*f[0]*wgt[0])*volFact; 
  out[1] += (vmap[1]*wgt[1]*f[3]+wgt[0]*vmap[1]*f[2]+vmap[0]*f[1]*wgt[1]+f[0]*vmap[0]*wgt[0])*volFact; 
  out[2] += (0.632455532033676*vmap1R2*wgt[1]*f[5]+0.6324555320336759*wgt[0]*vmap1R2*f[4]+1.4142135623730951*vmap[0]*vmap[1]*wgt[1]*f[3]+1.4142135623730951*vmap[0]*wgt[0]*vmap[1]*f[2]+0.7071067811865475*f[1]*vmap1R2*wgt[1]+0.7071067811865475*vmap0R2*f[1]*wgt[1]+0.7071067811865475*f[0]*wgt[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2*wgt[0])*volFact; 
} 

GKYL_CU_DH void int_mom_gyrokinetic_HAMILTONIAN_NOWEIGHT_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, const double *wgt, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[2] += (1.4142135623730951*f[1]*phi[1]*q_+1.4142135623730951*f[0]*phi[0]*q_)*volFact+(0.4472135954999579*vmap1R2*f[4]+vmap[0]*vmap[1]*f[2]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact; 
} 

