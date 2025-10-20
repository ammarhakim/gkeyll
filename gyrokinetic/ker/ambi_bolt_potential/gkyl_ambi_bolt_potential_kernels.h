// Gkyl ------------------------------------------------------------------------
//
// Header file for Ambipolar Boltzmann electron potential solver.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
#pragma once

#include <gkyl_util.h>
#include <math.h>


EXTERN_C_BEG

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_ser_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_3x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_3x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_3x_ser_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 


  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_ser_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_3x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_3x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_3x_ser_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 


  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_tensor_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_tensor_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_3x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_3x_tensor_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_3x_tensor_p1(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 


  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_tensor_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_tensor_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 

  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_3x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_3x_tensor_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out); 
  GKYL_CU_DH void ambi_bolt_potential_phi_calc_3x_tensor_p2(double q_e, double T_e, const double *m0Ion, const double *sheathvals, double *phi); 



EXTERN_C_END

