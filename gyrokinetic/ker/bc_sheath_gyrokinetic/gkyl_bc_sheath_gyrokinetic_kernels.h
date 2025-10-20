#pragma once

#include <gkyl_util.h>
#include <math.h>

// approximation for inverse Langevin function 
GKYL_CU_DH
static inline double invL(double x) {
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

EXTERN_C_BEG

GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 
GKYL_CU_DH void bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1(const double *vmap, const double q2Dm, const double *phi, const double *phiWall, const double *f, double *fRefl); 

EXTERN_C_END

