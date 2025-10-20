#pragma once 
#include <math.h> 
#include <gkyl_util.h> 

GKYL_CU_DH static inline void 
iz_react_rate_1x_ser_p1(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 
 
  double n_neut0 = 0.7071067811865476*n_neut[0]; 
  double vt_sq_neut0 = 0.7071067811865476*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.7071067811865476*vt_sq_elc[0]; 
  double T0 = (0.7071067811865476*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
  }
} 

GKYL_CU_DH static inline void 
iz_react_rate_1x_ser_p2(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 

  double n_neut0 = 0.7071067811865476*n_neut[0]; 
  double vt_sq_neut0 = 0.7071067811865476*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.7071067811865476*vt_sq_elc[0]; 
  double T0 = (0.7071067811865476*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
  }
}

GKYL_CU_DH static inline void 
iz_react_rate_2x_ser_p1(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 

  double n_neut0 = 0.5*n_neut[0]; 
  double vt_sq_neut0 = 0.5*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.5*vt_sq_elc[0]; 
  double T0 = (0.5*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
  }
}

GKYL_CU_DH static inline void 
iz_react_rate_2x_ser_p2(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 

  double n_neut0 = 0.5*n_neut[0]; 
  double vt_sq_neut0 = 0.5*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.5*vt_sq_elc[0]; 
  double T0 = (0.5*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
  }
}

GKYL_CU_DH static inline void 
iz_react_rate_3x_ser_p1(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 

  double n_neut0 = 0.3535533905932738*n_neut[0]; 
  double vt_sq_neut0 = 0.3535533905932738*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.3535533905932738*vt_sq_elc[0]; 
  double T0 = (0.3535533905932738*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
  }
}

GKYL_CU_DH static inline void 
iz_react_rate_3x_ser_p2(double elem_charge, double mass, double E, double A, double K, double P, double X, const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz) 
{ 
  // elem_charge : elementary charge (J - eV conversion factor). 
  // mass :        mass of electron
  // E :           Voronov ionization energy. 
  // A :           Voronov constant. 
  // K :           Voronov constant. 
  // P :           Voronov constant. 
  // X :           Voronov constant. 
  // n_neut :      neutral density. 
  // vt_sq_neut :  neutral squared thermal speed, sqrt(T/m). 
  // vt_sq_elc :   electron squared thermal speed, sqrt(T/m). 
  // coef_iz :     ionization reaction rate. 

  double n_neut0 = 0.3535533905932738*n_neut[0]; 
  double vt_sq_neut0 = 0.3535533905932738*vt_sq_neut[0]; 
  double vt_sq_elc0 = 0.3535533905932738*vt_sq_elc[0]; 
  double T0 = (0.3535533905932738*vt_sq_elc[0]*mass)/elem_charge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || n_neut0 <= 0.0 || vt_sq_neut0 <= 0.0 || vt_sq_elc0 <= 0.0) { 
    coef_iz[0] = 0.0;
  } 
  else {
    coef_iz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
  }
}
