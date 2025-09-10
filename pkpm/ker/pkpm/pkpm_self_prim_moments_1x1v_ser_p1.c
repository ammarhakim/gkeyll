#include <gkyl_prim_lbo_pkpm_kernels.h> 
 
GKYL_CU_DH void pkpm_self_prim_moments_1x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // vlasov_pkpm_moms:     [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model. 
  // boundary_corrections: boundary corrections to vtSq. 
  // nu:                   collision frequency. 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[2]; 
  const double *p_perp = &vlasov_pkpm_moms[4]; 
  const double *M1 = &vlasov_pkpm_moms[6]; 
  // If a corner value is below zero, use cell average.
  bool cellAvg = false;
  if (-(0.5*(2.4494897427831783*rho[1]-1.4142135623730951*rho[0])) < 0) cellAvg = true; 
  if (-(0.5*(2.4494897427831783*p_parallel[1]-1.4142135623730951*p_parallel[0])) < 0) cellAvg = true; 
  if (-(0.5*(2.4494897427831783*p_perp[1]-1.4142135623730951*p_perp[0])) < 0) cellAvg = true; 
  if (0.5*(2.4494897427831783*rho[1]+1.4142135623730951*rho[0]) < 0) cellAvg = true; 
  if (0.5*(2.4494897427831783*p_parallel[1]+1.4142135623730951*p_parallel[0]) < 0) cellAvg = true; 
  if (0.5*(2.4494897427831783*p_perp[1]+1.4142135623730951*p_perp[0]) < 0) cellAvg = true; 
 
  double m0r[2] = {0.0}; 
  double m1r[2] = {0.0}; 
  double cMr[2] = {0.0}; 
  double cEr[2] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = rho[0]; 
    m0r[1] = 0.0; 
    m1r[0] = M1[0]; 
    m1r[1] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    gkyl_mat_set(rhs,0,0,M1[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = 0.0; 
    gkyl_mat_set(rhs,2,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,3,0,0.0); 
  } else { 
    m0r[0] = rho[0]; 
    m0r[1] = rho[1]; 
    m1r[0] = M1[0]; 
    m1r[1] = M1[1]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    gkyl_mat_set(rhs,0,0,M1[0]); 
    gkyl_mat_set(rhs,1,0,M1[1]); 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = boundary_corrections[3]; 
    gkyl_mat_set(rhs,2,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,3,0,2.0*p_perp[1]+p_parallel[1]); 
  } 
 
  // ....... Block from weak multiply of u (correction to M1) and rho  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to u (correction to M1).......... // 
  gkyl_mat_set(A,0,2,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,3,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,2,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,3,-(0.7071067811865475*cMr[0])); 
 
  // ....... Block from weak multiply of u (correction to M1) and M1  .......... // 
  gkyl_mat_set(A,2,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,2,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,1,0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,2,2,2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,2,3,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,2,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,3,2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
