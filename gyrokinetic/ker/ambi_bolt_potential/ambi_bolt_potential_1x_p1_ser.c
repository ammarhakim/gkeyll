#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // cmag:        Clebsch function in definition of magnetic field.
  // jacobtotInv: reciprocal of the phase-space and conf-space Jacobians (1/(J*B)).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0Ion:       ion density.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  // Particle number density evaluate at the sheath entrance
  out[0] = m0Ion[0]-1.7320508075688772*m0Ion[1]; 

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = -(0.7329037678543797*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e))+0.42314218766081724*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)-0.24430125595145993*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)-0.24430125595145993*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)-0.24430125595145993*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e)+0.14104739588693907*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 

  double m0JacIonB[1];
  m0JacIonB[0] = 0.7071067811865475*m0JacIon[0]-1.224744871391589*m0JacIon[1]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[2] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // cmag:        Clebsch function in definition of magnetic field.
  // jacobtotInv: reciprocal of the phase-space and conf-space Jacobians (1/(J*B)).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0Ion:       ion density.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  // Particle number density evaluate at the sheath entrance
  out[0] = 1.7320508075688772*m0Ion[1]+m0Ion[0]; 

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.3535533905932737*GammaJac_i[0]*sheathDirDx-0.6123724356957944*GammaJac_i[1]*sheathDirDx; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = 0.7329037678543797*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.24430125595145993*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)+0.24430125595145993*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)+0.24430125595145993*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e)+0.14104739588693907*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[2] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi) 
{ 
  // q_e: electron change.
  // T_e: electron temperature.
  // m0Ion: ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi: electrostatic potential in domain volume.

  double phi_qp[2];
  phi_qp[0] = -((1.0*log(m0Ion[0]/(sheathvals[0]-1.0*sheathvals[1])-(1.0*m0Ion[1])/(sheathvals[0]-1.0*sheathvals[1]))*T_e)/q_e)-0.7071067811865475*sheathvals[3]+0.7071067811865475*sheathvals[2]; 
  phi_qp[1] = -((1.0*log(m0Ion[1]/(sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.7071067811865475*sheathvals[3]+0.7071067811865475*sheathvals[2]; 

  phi[0] = 0.7071067811865475*(phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.7071067811865475*phi_qp[1]-0.7071067811865475*phi_qp[0]; 
}

