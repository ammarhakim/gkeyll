#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
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
  out[0] = 2.23606797749979*m0Ion[2]-1.7320508075688772*m0Ion[1]+m0Ion[0]; 

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.7905694150420947*GammaJac_i[2]*sheathDirDx+0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = 1.5769578262626003*cmag[2]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)-1.2215062797572998*cmag[1]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.7052369794346953*cmag[0]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)-1.2215062797572998*jacobtotInv[1]*cmag[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.7052369794346953*jacobtotInv[0]*cmag[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.9461746957575601*cmag[1]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)-0.5462742152960396*cmag[0]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)-0.5462742152960396*jacobtotInv[0]*cmag[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.31539156525252*cmag[0]*jacobtotInv[0]*m0JacIon[2]*sqrt(T_e/m_e)-1.2215062797572998*m0JacIon[1]*cmag[2]*jacobtotInv[2]*sqrt(T_e/m_e)+0.7052369794346953*m0JacIon[0]*cmag[2]*jacobtotInv[2]*sqrt(T_e/m_e)+0.9461746957575601*cmag[1]*m0JacIon[1]*jacobtotInv[2]*sqrt(T_e/m_e)-0.5462742152960396*cmag[0]*m0JacIon[1]*jacobtotInv[2]*sqrt(T_e/m_e)-0.5462742152960396*m0JacIon[0]*cmag[1]*jacobtotInv[2]*sqrt(T_e/m_e)+0.31539156525252*cmag[0]*m0JacIon[0]*jacobtotInv[2]*sqrt(T_e/m_e)+0.9461746957575601*jacobtotInv[1]*m0JacIon[1]*cmag[2]*sqrt(T_e/m_e)-0.5462742152960396*jacobtotInv[0]*m0JacIon[1]*cmag[2]*sqrt(T_e/m_e)-0.5462742152960396*m0JacIon[0]*jacobtotInv[1]*cmag[2]*sqrt(T_e/m_e)+0.31539156525252*jacobtotInv[0]*m0JacIon[0]*cmag[2]*sqrt(T_e/m_e)-0.7329037678543797*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)-0.24430125595145993*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)-0.24430125595145993*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)-0.24430125595145993*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e)+0.14104739588693907*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.5811388300841895*m0JacIon[2]-1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[3] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
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
  out[0] = 2.23606797749979*m0Ion[2]+1.7320508075688772*m0Ion[1]+m0Ion[0]; 

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.7905694150420947*GammaJac_i[2]*sheathDirDx-0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = 1.5769578262626003*cmag[2]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)+1.2215062797572998*cmag[1]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.7052369794346953*cmag[0]*jacobtotInv[2]*m0JacIon[2]*sqrt(T_e/m_e)+1.2215062797572998*jacobtotInv[1]*cmag[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.7052369794346953*jacobtotInv[0]*cmag[2]*m0JacIon[2]*sqrt(T_e/m_e)+0.9461746957575601*cmag[1]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.5462742152960396*cmag[0]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.5462742152960396*jacobtotInv[0]*cmag[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.31539156525252*cmag[0]*jacobtotInv[0]*m0JacIon[2]*sqrt(T_e/m_e)+1.2215062797572998*m0JacIon[1]*cmag[2]*jacobtotInv[2]*sqrt(T_e/m_e)+0.7052369794346953*m0JacIon[0]*cmag[2]*jacobtotInv[2]*sqrt(T_e/m_e)+0.9461746957575601*cmag[1]*m0JacIon[1]*jacobtotInv[2]*sqrt(T_e/m_e)+0.5462742152960396*cmag[0]*m0JacIon[1]*jacobtotInv[2]*sqrt(T_e/m_e)+0.5462742152960396*m0JacIon[0]*cmag[1]*jacobtotInv[2]*sqrt(T_e/m_e)+0.31539156525252*cmag[0]*m0JacIon[0]*jacobtotInv[2]*sqrt(T_e/m_e)+0.9461746957575601*jacobtotInv[1]*m0JacIon[1]*cmag[2]*sqrt(T_e/m_e)+0.5462742152960396*jacobtotInv[0]*m0JacIon[1]*cmag[2]*sqrt(T_e/m_e)+0.5462742152960396*m0JacIon[0]*jacobtotInv[1]*cmag[2]*sqrt(T_e/m_e)+0.31539156525252*jacobtotInv[0]*m0JacIon[0]*cmag[2]*sqrt(T_e/m_e)+0.7329037678543797*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.24430125595145993*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.42314218766081724*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)+0.24430125595145993*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)+0.24430125595145993*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e)+0.14104739588693907*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.5811388300841895*m0JacIon[2]+1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[3] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi) 
{ 
  // q_e: electron change.
  // T_e: electron temperature.
  // m0Ion: ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi: electrostatic potential in domain volume.

  double phi_qp[3];
  phi_qp[0] = -((1.0*log(-((3.0*m0Ion[1])/(2.0*sheathvals[2]-3.0*sheathvals[1]+2.23606797749979*sheathvals[0]))+(1.4142135623730951*m0Ion[2])/(1.4142135623730951*sheathvals[2]-2.1213203435596424*sheathvals[1]+1.5811388300841895*sheathvals[0])+m0Ion[0]/(0.8944271909999159*sheathvals[2]-1.3416407864998738*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.6324555320336759*sheathvals[5]-0.9486832980505137*sheathvals[4]+0.7071067811865475*sheathvals[3]; 
  phi_qp[1] = -((1.0*log(m0Ion[0]/(sheathvals[0]-1.118033988749895*sheathvals[2])-(2.23606797749979*m0Ion[2])/(2.0*sheathvals[0]-2.23606797749979*sheathvals[2]))*T_e)/q_e)-0.7905694150420947*sheathvals[5]+0.7071067811865475*sheathvals[3]; 
  phi_qp[2] = -((1.0*log((3.0*m0Ion[1])/(2.0*sheathvals[2]+3.0*sheathvals[1]+2.23606797749979*sheathvals[0])+(1.4142135623730951*m0Ion[2])/(1.4142135623730951*sheathvals[2]+2.1213203435596424*sheathvals[1]+1.5811388300841895*sheathvals[0])+m0Ion[0]/(0.8944271909999159*sheathvals[2]+1.3416407864998738*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.6324555320336759*sheathvals[5]+0.9486832980505137*sheathvals[4]+0.7071067811865475*sheathvals[3]; 

  phi[0] = 0.39283710065919303*phi_qp[2]+0.6285393610547091*phi_qp[1]+0.39283710065919303*phi_qp[0]; 
  phi[1] = 0.5270462766947298*phi_qp[2]-0.5270462766947298*phi_qp[0]; 
  phi[2] = 0.35136418446315326*phi_qp[2]-0.7027283689263066*phi_qp[1]+0.35136418446315326*phi_qp[0]; 
}

