#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
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
  out[0] = m0Ion[0]-1.7320508075688772*m0Ion[2]; 
  out[1] = m0Ion[1]-1.7320508075688772*m0Ion[3]; 

  double GammaJacIonB[2];
  GammaJacIonB[0] = 0.6123724356957944*GammaJac_i[2]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 
  GammaJacIonB[1] = 0.6123724356957944*GammaJac_i[3]*sheathDirDx+0.3535533905932737*GammaJac_i[1]*sheathDirDx; 

  double x3HalfMomJacElcB[2];
  x3HalfMomJacElcB[0] = -(0.12215062797572997*cmag[0]*jacobtotInv[1]*m0JacIon[3]*sqrt(T_e/m_e))-0.12215062797572997*jacobtotInv[0]*cmag[1]*m0JacIon[3]*sqrt(T_e/m_e)-0.12215062797572997*cmag[1]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)-0.12215062797572997*cmag[0]*jacobtotInv[0]*m0JacIon[2]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 
  x3HalfMomJacElcB[1] = -(0.2198711303563139*cmag[1]*jacobtotInv[1]*m0JacIon[3]*sqrt(T_e/m_e))-0.12215062797572997*cmag[0]*jacobtotInv[0]*m0JacIon[3]*sqrt(T_e/m_e)-0.12215062797572997*cmag[0]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)-0.12215062797572997*jacobtotInv[0]*cmag[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.12694265629824517*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)+0.07052369794346953*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e); 

  double m0JacIonB[2];
  m0JacIonB[0] = 0.7071067811865475*m0JacIon[0]-1.224744871391589*m0JacIon[2]; 
  m0JacIonB[1] = 0.7071067811865475*m0JacIon[1]-1.224744871391589*m0JacIon[3]; 

  double phiS_qp[2];
  if ((isfinite(0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1])) && (0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1]>0.) && (0.7071067811865475*m0JacIonB[0]-0.7071067811865475*m0JacIonB[1]>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/(x3HalfMomJacElcB[0]-1.0*x3HalfMomJacElcB[1])-(1.0*GammaJacIonB[1])/(x3HalfMomJacElcB[0]-1.0*x3HalfMomJacElcB[1]))*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0])) && (0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0]>0.) && (0.7071067811865475*m0JacIonB[1]+0.7071067811865475*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (log(GammaJacIonB[1]/(x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }

  // Sheath potential
  out[4] = phiS_qp[1]+phiS_qp[0]; 
  out[5] = phiS_qp[1]-1.0*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out) 
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
  out[0] = 1.7320508075688772*m0Ion[2]+m0Ion[0]; 
  out[1] = 1.7320508075688772*m0Ion[3]+m0Ion[1]; 

  double GammaJacIonB[2];
  GammaJacIonB[0] = 0.3535533905932737*GammaJac_i[0]*sheathDirDx-0.6123724356957944*GammaJac_i[2]*sheathDirDx; 
  GammaJacIonB[1] = 0.3535533905932737*GammaJac_i[1]*sheathDirDx-0.6123724356957944*GammaJac_i[3]*sheathDirDx; 

  double x3HalfMomJacElcB[2];
  x3HalfMomJacElcB[0] = 0.12215062797572997*cmag[0]*jacobtotInv[1]*m0JacIon[3]*sqrt(T_e/m_e)+0.12215062797572997*jacobtotInv[0]*cmag[1]*m0JacIon[3]*sqrt(T_e/m_e)+0.12215062797572997*cmag[1]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.12215062797572997*cmag[0]*jacobtotInv[0]*m0JacIon[2]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*jacobtotInv[0]*cmag[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*m0JacIon[0]*cmag[1]*jacobtotInv[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[0]*m0JacIon[0]*sqrt(T_e/m_e); 
  x3HalfMomJacElcB[1] = 0.2198711303563139*cmag[1]*jacobtotInv[1]*m0JacIon[3]*sqrt(T_e/m_e)+0.12215062797572997*cmag[0]*jacobtotInv[0]*m0JacIon[3]*sqrt(T_e/m_e)+0.12215062797572997*cmag[0]*jacobtotInv[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.12215062797572997*jacobtotInv[0]*cmag[1]*m0JacIon[2]*sqrt(T_e/m_e)+0.12694265629824517*cmag[1]*jacobtotInv[1]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*jacobtotInv[0]*m0JacIon[1]*sqrt(T_e/m_e)+0.07052369794346953*cmag[0]*m0JacIon[0]*jacobtotInv[1]*sqrt(T_e/m_e)+0.07052369794346953*jacobtotInv[0]*m0JacIon[0]*cmag[1]*sqrt(T_e/m_e); 

  double m0JacIonB[2];
  m0JacIonB[0] = 1.224744871391589*m0JacIon[2]+0.7071067811865475*m0JacIon[0]; 
  m0JacIonB[1] = 1.224744871391589*m0JacIon[3]+0.7071067811865475*m0JacIon[1]; 

  double phiS_qp[2];
  if ((isfinite(0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1])) && (0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1]>0.) && (0.7071067811865475*m0JacIonB[0]-0.7071067811865475*m0JacIonB[1]>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/(x3HalfMomJacElcB[0]-1.0*x3HalfMomJacElcB[1])-(1.0*GammaJacIonB[1])/(x3HalfMomJacElcB[0]-1.0*x3HalfMomJacElcB[1]))*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0])) && (0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0]>0.) && (0.7071067811865475*m0JacIonB[1]+0.7071067811865475*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (log(GammaJacIonB[1]/(x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }

  // Sheath potential
  out[4] = phiS_qp[1]+phiS_qp[0]; 
  out[5] = phiS_qp[1]-1.0*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_ser_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi) 
{ 
  // q_e: electron change.
  // T_e: electron temperature.
  // m0Ion: ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi: electrostatic potential in domain volume.

  double phi_qp[4];
  phi_qp[0] = -((1.0*log(m0Ion[3]/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.5*sheathvals[7]-0.5*sheathvals[6]-0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[1] = -((1.0*log(-((1.0*m0Ion[3])/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))+m0Ion[2]/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.5*sheathvals[7]+0.5*sheathvals[6]-0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[2] = -((1.0*log(-((1.0*m0Ion[3])/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))-(1.0*m0Ion[2])/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.5*sheathvals[7]-0.5*sheathvals[6]+0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[3] = -((1.0*log(m0Ion[3]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[2]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.5*sheathvals[7]+0.5*sheathvals[6]+0.5*sheathvals[5]+0.5*sheathvals[4]; 

  phi[0] = 0.5*(phi_qp[3]+phi_qp[2]+phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.5*(phi_qp[3]+phi_qp[2])-0.5*(phi_qp[1]+phi_qp[0]); 
  phi[2] = 0.5*phi_qp[3]-0.5*phi_qp[2]+0.5*phi_qp[1]-0.5*phi_qp[0]; 
  phi[3] = 0.5*phi_qp[3]-0.5*(phi_qp[2]+phi_qp[1])+0.5*phi_qp[0]; 
}

