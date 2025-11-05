#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_apardot_vol_1x1v_ser_p1(const double q_, const double m_, const double *apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // apardot: time derivative of parallel component of magnetic vector potential.
  // fin: Distribution function.
  // out: output increment.

  double alphavpar[6] = {0.}; 
  alphavpar[0] = -((1.0*apardot[0]*q_)/m_); 
  alphavpar[1] = -((1.0*apardot[1]*q_)/m_); 
  alphavpar[2] = -((1.0*apardot[2]*q_)/m_); 
  alphavpar[3] = -((1.0*apardot[3]*q_)/m_); 
  alphavpar[4] = -((1.0*apardot[4]*q_)/m_); 
  alphavpar[5] = -((1.0*apardot[5]*q_)/m_); 


  out[2] += 0.8660254037844386*(alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.8660254037844387*(alphavpar[4]*fin[5]+fin[4]*alphavpar[5])+0.8660254037844386*(alphavpar[2]*fin[3]+fin[2]*alphavpar[3]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.7320508075688774*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5])+1.7320508075688772*(alphavpar[2]*fin[4]+fin[2]*alphavpar[4])+1.9364916731037085*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[5] += 1.7320508075688772*(alphavpar[2]*fin[5]+fin[2]*alphavpar[5])+1.7320508075688774*(alphavpar[3]*fin[4]+fin[3]*alphavpar[4])+1.9364916731037085*(alphavpar[0]*fin[3]+fin[0]*alphavpar[3]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 

  return 0.; 
} 
