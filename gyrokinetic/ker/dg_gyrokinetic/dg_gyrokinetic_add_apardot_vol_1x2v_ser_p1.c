#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_1x2v_ser_p1(const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[12] = {0.}; 
  alphavpar[0] = -((1.0*apardot[0]*q_)/m_); 
  alphavpar[1] = -((1.0*apardot[1]*q_)/m_); 
  alphavpar[2] = -((1.0*apardot[2]*q_)/m_); 
  alphavpar[3] = -((1.0*apardot[3]*q_)/m_); 
  alphavpar[4] = -((1.0*apardot[4]*q_)/m_); 
  alphavpar[5] = -((1.0*apardot[5]*q_)/m_); 
  alphavpar[6] = -((1.0*apardot[6]*q_)/m_); 
  alphavpar[7] = -((1.0*apardot[7]*q_)/m_); 
  alphavpar[8] = -((1.0*apardot[8]*q_)/m_); 
  alphavpar[9] = -((1.0*apardot[9]*q_)/m_); 
  alphavpar[10] = -((1.0*apardot[10]*q_)/m_); 
  alphavpar[11] = -((1.0*apardot[11]*q_)/m_); 


  out[2] += 0.6123724356957944*(alphavpar[11]*fin[11]+alphavpar[10]*fin[10]+alphavpar[9]*fin[9]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.6123724356957944*(alphavpar[10]*fin[11]+fin[10]*alphavpar[11]+alphavpar[8]*fin[9]+fin[8]*alphavpar[9]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphavpar[2]*fin[4]+fin[2]*alphavpar[4]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[6] += 0.6123724356957944*(alphavpar[9]*fin[11]+fin[9]*alphavpar[11]+alphavpar[8]*fin[10]+fin[8]*alphavpar[10]+alphavpar[4]*fin[7]+fin[4]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*(alphavpar[8]*fin[11]+fin[8]*alphavpar[11]+alphavpar[9]*fin[10]+fin[9]*alphavpar[10]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[8] += 1.224744871391589*(alphavpar[7]*fin[11]+fin[7]*alphavpar[11]+alphavpar[6]*fin[10]+fin[6]*alphavpar[10]+alphavpar[4]*fin[9]+fin[4]*alphavpar[9]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8])+1.369306393762915*(alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[9] += 1.224744871391589*(alphavpar[6]*fin[11]+fin[6]*alphavpar[11]+alphavpar[7]*fin[10]+fin[7]*alphavpar[10]+alphavpar[2]*fin[9]+fin[2]*alphavpar[9]+alphavpar[4]*fin[8]+fin[4]*alphavpar[8])+1.369306393762915*(alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[10] += 1.224744871391589*(alphavpar[4]*fin[11]+fin[4]*alphavpar[11]+alphavpar[2]*fin[10]+fin[2]*alphavpar[10]+alphavpar[7]*fin[9]+fin[7]*alphavpar[9]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8])+1.369306393762915*(alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[11] += 1.224744871391589*(alphavpar[2]*fin[11]+fin[2]*alphavpar[11]+alphavpar[4]*fin[10]+fin[4]*alphavpar[10]+alphavpar[6]*fin[9]+fin[6]*alphavpar[9]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8])+1.369306393762915*(alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4]); 

  return 0.; 
} 
