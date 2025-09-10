#include <gkyl_prim_lbo_pkpm_kernels.h> 
 
GKYL_CU_DH void pkpm_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // vlasov_pkpm_moms:     [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model. 
  // boundary_corrections: boundary corrections to vtSq. 
  // nu:                   collision frequency. 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 
  const double *M1 = &vlasov_pkpm_moms[9]; 
  // If a corner value is below zero, use cell average.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*rho[2]-1.7320508075688772*rho[1]+rho[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*p_parallel[2]-1.7320508075688772*p_parallel[1]+p_parallel[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*p_perp[2]-1.7320508075688772*p_perp[1]+p_perp[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*rho[2]+1.7320508075688772*rho[1]+rho[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*p_parallel[2]+1.7320508075688772*p_parallel[1]+p_parallel[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*p_perp[2]+1.7320508075688772*p_perp[1]+p_perp[0]) < 0) cellAvg = true; 
 
  double m0r[3] = {0.0}; 
  double m1r[3] = {0.0}; 
  double cMr[3] = {0.0}; 
  double cEr[3] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = rho[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = M1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    gkyl_mat_set(rhs,0,0,M1[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cEr[0] = boundary_corrections[3]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,3,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
  } else { 
    m0r[0] = rho[0]; 
    m0r[1] = rho[1]; 
    m0r[2] = rho[2]; 
    m1r[0] = M1[0]; 
    m1r[1] = M1[1]; 
    m1r[2] = M1[2]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    gkyl_mat_set(rhs,0,0,M1[0]); 
    gkyl_mat_set(rhs,1,0,M1[1]); 
    gkyl_mat_set(rhs,2,0,M1[2]); 
    cEr[0] = boundary_corrections[3]; 
    cEr[1] = boundary_corrections[4]; 
    cEr[2] = boundary_corrections[5]; 
    gkyl_mat_set(rhs,3,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,4,0,2.0*p_perp[1]+p_parallel[1]); 
    gkyl_mat_set(rhs,5,0,2.0*p_perp[2]+p_parallel[2]); 
  } 
 
  // ....... Block from weak multiply of u (correction to M1) and rho  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,2,0.45175395145262565*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to u (correction to M1).......... // 
  gkyl_mat_set(A,0,3,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,4,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,0,5,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,1,3,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,4,-(0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,5,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,3,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,2,4,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,5,-(0.45175395145262565*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of u (correction to M1) and M1  .......... // 
  gkyl_mat_set(A,3,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,3,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,4,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,4,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,4,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,5,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,2,0.45175395145262565*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,3,3,2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,3,4,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,5,2.1213203435596424*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,4,3,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,4,4,1.8973665961010278*m0r[2]-0.6324555320336759*cEr[2]+2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,4,5,1.8973665961010278*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,3,2.1213203435596424*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,5,4,1.8973665961010278*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,5,1.355261854357877*m0r[2]-0.45175395145262565*cEr[2]+2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
