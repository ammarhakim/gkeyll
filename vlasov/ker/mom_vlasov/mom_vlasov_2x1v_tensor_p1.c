#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_M0_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M2ij_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  out[0] += (0.421637021355784*f[7]*vmap_vx3_sq+0.7071067811865475*f[0]*vmap_vx3_sq+0.45175395145262565*f[7]*vmap_vx2_sq+0.7071067811865475*f[0]*vmap_vx2_sq+0.6324555320336759*f[7]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx1_sq+0.7071067811865475*f[0]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vx[3]*f[7]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[7]+1.2421180068162374*vmap_vx[2]*f[3]*vmap_vx[3]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[3]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[3])*volFact; 
  out[1] += (0.4216370213557841*f[9]*vmap_vx3_sq+0.7071067811865475*f[1]*vmap_vx3_sq+0.4517539514526257*f[9]*vmap_vx2_sq+0.7071067811865475*f[1]*vmap_vx2_sq+0.632455532033676*f[9]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx1_sq+0.7071067811865475*f[1]*vmap_vx0_sq+1.2421180068162376*vmap_vx[1]*vmap_vx[3]*f[9]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[9]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[5]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[5]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[5])*volFact; 
  out[2] += (0.4216370213557841*f[10]*vmap_vx3_sq+0.7071067811865475*f[2]*vmap_vx3_sq+0.4517539514526257*f[10]*vmap_vx2_sq+0.7071067811865475*f[2]*vmap_vx2_sq+0.632455532033676*f[10]*vmap_vx1_sq+0.7071067811865475*f[2]*vmap_vx1_sq+0.7071067811865475*f[2]*vmap_vx0_sq+1.2421180068162376*vmap_vx[1]*vmap_vx[3]*f[10]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[10]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[6]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[6]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[6])*volFact; 
  out[3] += (0.421637021355784*f[11]*vmap_vx3_sq+0.7071067811865475*f[4]*vmap_vx3_sq+0.45175395145262565*f[11]*vmap_vx2_sq+0.7071067811865475*f[4]*vmap_vx2_sq+0.6324555320336759*f[11]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx1_sq+0.7071067811865475*f[4]*vmap_vx0_sq+1.2421180068162374*vmap_vx[1]*vmap_vx[3]*f[11]+1.4142135623730951*vmap_vx[0]*vmap_vx[2]*f[11]+1.2421180068162374*vmap_vx[2]*vmap_vx[3]*f[8]+1.264911064067352*vmap_vx[1]*vmap_vx[2]*f[8]+1.4142135623730951*vmap_vx[0]*vmap_vx[1]*f[8])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M3ijk_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
  const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  const double *vmap_vx = &vmap[0]; 
  const double vmap_vx0_sq = vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx0_cu = vmap_vx[0]*vmap_vx[0]*vmap_vx[0]; 
  const double vmap_vx1_sq = vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx1_cu = vmap_vx[1]*vmap_vx[1]*vmap_vx[1]; 
  const double vmap_vx2_sq = vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx2_cu = vmap_vx[2]*vmap_vx[2]*vmap_vx[2]; 
  const double vmap_vx3_sq = vmap_vx[3]*vmap_vx[3]; 
  const double vmap_vx3_cu = vmap_vx[3]*vmap_vx[3]*vmap_vx[3]; 
  out[0] += (2.772727272727273*vmap_vx[2]*f[7]*vmap_vx3_sq+0.8944271909999159*vmap_vx[0]*f[7]*vmap_vx3_sq+2.3*vmap_vx[1]*f[3]*vmap_vx3_sq+0.8944271909999159*f[0]*vmap_vx[2]*vmap_vx3_sq+1.5*f[0]*vmap_vx[0]*vmap_vx3_sq+0.499917348540637*f[3]*vmap_vx3_cu+0.9583148474999099*vmap_vx[0]*f[7]*vmap_vx2_sq+1.963961012123931*f[3]*vmap_vx[3]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[3]*vmap_vx2_sq+1.5*f[0]*vmap_vx[0]*vmap_vx2_sq+1.0714285714285714*f[7]*vmap_vx2_cu+0.31943828249996997*f[0]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[7]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[7]*vmap_vx1_sq+1.1783766072743584*f[3]*vmap_vx[3]*vmap_vx1_sq+1.3416407864998738*f[0]*vmap_vx[2]*vmap_vx1_sq+1.5*f[0]*vmap_vx[0]*vmap_vx1_sq+0.9*f[3]*vmap_vx1_cu+1.5*vmap_vx[2]*f[7]*vmap_vx0_sq+1.5*vmap_vx[1]*f[3]*vmap_vx0_sq+0.5*f[0]*vmap_vx0_cu+3.927922024247862*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[7]+2.6349301969610384*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[7]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*f[3]*vmap_vx[3]+2.6349301969610384*f[0]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[3])*volFact; 
  out[1] += (2.7727272727272734*vmap_vx[2]*f[9]*vmap_vx3_sq+0.8944271909999161*vmap_vx[0]*f[9]*vmap_vx3_sq+2.3*vmap_vx[1]*f[5]*vmap_vx3_sq+0.8944271909999159*f[1]*vmap_vx[2]*vmap_vx3_sq+1.5*vmap_vx[0]*f[1]*vmap_vx3_sq+0.499917348540637*f[5]*vmap_vx3_cu+0.9583148474999098*vmap_vx[0]*f[9]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[5]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[5]*vmap_vx2_sq+1.5*vmap_vx[0]*f[1]*vmap_vx2_sq+1.0714285714285714*f[9]*vmap_vx2_cu+0.31943828249996997*f[1]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[9]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[9]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[5]*vmap_vx1_sq+1.3416407864998738*f[1]*vmap_vx[2]*vmap_vx1_sq+1.5*vmap_vx[0]*f[1]*vmap_vx1_sq+0.9*f[5]*vmap_vx1_cu+1.5*vmap_vx[2]*f[9]*vmap_vx0_sq+1.5*vmap_vx[1]*f[5]*vmap_vx0_sq+0.5*f[1]*vmap_vx0_cu+3.927922024247863*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[9]+2.6349301969610397*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[9]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[5]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[5]+2.6349301969610384*f[1]*vmap_vx[1]*vmap_vx[2]*vmap_vx[3])*volFact; 
  out[2] += (2.7727272727272734*vmap_vx[2]*f[10]*vmap_vx3_sq+0.8944271909999161*vmap_vx[0]*f[10]*vmap_vx3_sq+2.3*vmap_vx[1]*f[6]*vmap_vx3_sq+0.8944271909999159*f[2]*vmap_vx[2]*vmap_vx3_sq+1.5*vmap_vx[0]*f[2]*vmap_vx3_sq+0.499917348540637*f[6]*vmap_vx3_cu+0.9583148474999098*vmap_vx[0]*f[10]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[6]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[6]*vmap_vx2_sq+1.5*vmap_vx[0]*f[2]*vmap_vx2_sq+1.0714285714285714*f[10]*vmap_vx2_cu+0.31943828249996997*f[2]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[10]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[10]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[6]*vmap_vx1_sq+1.3416407864998738*f[2]*vmap_vx[2]*vmap_vx1_sq+1.5*vmap_vx[0]*f[2]*vmap_vx1_sq+0.9*f[6]*vmap_vx1_cu+1.5*vmap_vx[2]*f[10]*vmap_vx0_sq+1.5*vmap_vx[1]*f[6]*vmap_vx0_sq+0.5*f[2]*vmap_vx0_cu+3.927922024247863*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[10]+2.6349301969610397*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[10]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[6]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[6]+2.6349301969610384*vmap_vx[1]*f[2]*vmap_vx[2]*vmap_vx[3])*volFact; 
  out[3] += (2.772727272727273*vmap_vx[2]*f[11]*vmap_vx3_sq+0.8944271909999159*vmap_vx[0]*f[11]*vmap_vx3_sq+2.3*vmap_vx[1]*f[8]*vmap_vx3_sq+0.8944271909999159*vmap_vx[2]*f[4]*vmap_vx3_sq+1.5*vmap_vx[0]*f[4]*vmap_vx3_sq+0.499917348540637*f[8]*vmap_vx3_cu+0.9583148474999099*vmap_vx[0]*f[11]*vmap_vx2_sq+1.963961012123931*vmap_vx[3]*f[8]*vmap_vx2_sq+2.357142857142857*vmap_vx[1]*f[8]*vmap_vx2_sq+1.5*vmap_vx[0]*f[4]*vmap_vx2_sq+1.0714285714285714*f[11]*vmap_vx2_cu+0.31943828249996997*f[4]*vmap_vx2_cu+2.357142857142857*vmap_vx[2]*f[11]*vmap_vx1_sq+1.3416407864998738*vmap_vx[0]*f[11]*vmap_vx1_sq+1.1783766072743584*vmap_vx[3]*f[8]*vmap_vx1_sq+1.3416407864998738*vmap_vx[2]*f[4]*vmap_vx1_sq+1.5*vmap_vx[0]*f[4]*vmap_vx1_sq+0.9*f[8]*vmap_vx1_cu+1.5*vmap_vx[2]*f[11]*vmap_vx0_sq+1.5*vmap_vx[1]*f[8]*vmap_vx0_sq+0.5*f[4]*vmap_vx0_cu+3.927922024247862*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[11]+2.6349301969610384*vmap_vx[0]*vmap_vx[1]*vmap_vx[3]*f[11]+2.6349301969610384*vmap_vx[0]*vmap_vx[2]*vmap_vx[3]*f[8]+2.6832815729997477*vmap_vx[0]*vmap_vx[1]*vmap_vx[2]*f[8]+2.6349301969610384*vmap_vx[1]*vmap_vx[2]*vmap_vx[3]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += ((1.2909944487358058*hamil[2]*f[7]+0.4303314829119353*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]+0.6454972243679029*hamil[1]*f[3]+1.4433756729740645*f[0]*hamil[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[7])+0.4303314829119353*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]-0.6454972243679029*hamil[1]*f[3]-1.4433756729740645*f[0]*hamil[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[0]*hamil[1]-0.8606629658238706*hamil[1]*f[7])*dv10*volFact)/jacob_vx[1]; 
  out[1] += ((1.2909944487358056*hamil[2]*f[9]+0.4303314829119352*hamil[1]*f[9]+1.9364916731037085*hamil[2]*f[5]+0.6454972243679029*hamil[1]*f[5]+1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[9])+0.4303314829119352*hamil[1]*f[9]+1.9364916731037085*hamil[2]*f[5]-0.6454972243679029*hamil[1]*f[5]-1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[1]*hamil[1]-0.8606629658238704*hamil[1]*f[9])*dv10*volFact)/jacob_vx[1]; 
  out[2] += ((1.2909944487358056*hamil[2]*f[10]+0.4303314829119352*hamil[1]*f[10]+1.9364916731037085*hamil[2]*f[6]+0.6454972243679029*hamil[1]*f[6]+1.4433756729740645*f[2]*hamil[2]+0.48112522432468824*hamil[1]*f[2])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[10])+0.4303314829119352*hamil[1]*f[10]+1.9364916731037085*hamil[2]*f[6]-0.6454972243679029*hamil[1]*f[6]-1.4433756729740645*f[2]*hamil[2]+0.48112522432468824*hamil[1]*f[2])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[2]-0.8606629658238704*hamil[1]*f[10])*dv10*volFact)/jacob_vx[1]; 
  out[3] += ((1.2909944487358058*hamil[2]*f[11]+0.4303314829119353*hamil[1]*f[11]+1.9364916731037085*hamil[2]*f[8]+0.6454972243679029*hamil[1]*f[8]+1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[11])+0.4303314829119353*hamil[1]*f[11]+1.9364916731037085*hamil[2]*f[8]-0.6454972243679029*hamil[1]*f[8]-1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[4]-0.8606629658238706*hamil[1]*f[11])*dv10*volFact)/jacob_vx[1]; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  out[0] += (hamil[2]*f[7]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[1] += (1.0000000000000002*hamil[2]*f[9]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[2] += (1.0000000000000002*hamil[2]*f[10]+hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[3] += (hamil[2]*f[11]+hamil[1]*f[8]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[2]/2; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  out[4] += ((1.2909944487358058*hamil[2]*f[7]+0.4303314829119353*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]+0.6454972243679029*hamil[1]*f[3]+1.4433756729740645*f[0]*hamil[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[7])+0.4303314829119353*hamil[1]*f[7]+1.9364916731037085*hamil[2]*f[3]-0.6454972243679029*hamil[1]*f[3]-1.4433756729740645*f[0]*hamil[2]+0.48112522432468824*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[0]*hamil[1]-0.8606629658238706*hamil[1]*f[7])*dv10*volFact)/jacob_vx[1]; 
  out[5] += ((1.2909944487358056*hamil[2]*f[9]+0.4303314829119352*hamil[1]*f[9]+1.9364916731037085*hamil[2]*f[5]+0.6454972243679029*hamil[1]*f[5]+1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[9])+0.4303314829119352*hamil[1]*f[9]+1.9364916731037085*hamil[2]*f[5]-0.6454972243679029*hamil[1]*f[5]-1.4433756729740645*f[1]*hamil[2]+0.48112522432468824*f[1]*hamil[1])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*f[1]*hamil[1]-0.8606629658238704*hamil[1]*f[9])*dv10*volFact)/jacob_vx[1]; 
  out[6] += ((1.2909944487358056*hamil[2]*f[10]+0.4303314829119352*hamil[1]*f[10]+1.9364916731037085*hamil[2]*f[6]+0.6454972243679029*hamil[1]*f[6]+1.4433756729740645*f[2]*hamil[2]+0.48112522432468824*hamil[1]*f[2])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358056*hamil[2]*f[10])+0.4303314829119352*hamil[1]*f[10]+1.9364916731037085*hamil[2]*f[6]-0.6454972243679029*hamil[1]*f[6]-1.4433756729740645*f[2]*hamil[2]+0.48112522432468824*hamil[1]*f[2])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[2]-0.8606629658238704*hamil[1]*f[10])*dv10*volFact)/jacob_vx[1]; 
  out[7] += ((1.2909944487358058*hamil[2]*f[11]+0.4303314829119353*hamil[1]*f[11]+1.9364916731037085*hamil[2]*f[8]+0.6454972243679029*hamil[1]*f[8]+1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[2]+((-(1.2909944487358058*hamil[2]*f[11])+0.4303314829119353*hamil[1]*f[11]+1.9364916731037085*hamil[2]*f[8]-0.6454972243679029*hamil[1]*f[8]-1.4433756729740645*hamil[2]*f[4]+0.48112522432468824*hamil[1]*f[4])*dv10*volFact)/jacob_vx[0]+((0.7698003589195012*hamil[1]*f[4]-0.8606629658238706*hamil[1]*f[11])*dv10*volFact)/jacob_vx[1]; 
  out[8] += (hamil[2]*f[7]+hamil[1]*f[3]+f[0]*hamil[0])*volFact; 
  out[9] += (1.0000000000000002*hamil[2]*f[9]+hamil[1]*f[5]+hamil[0]*f[1])*volFact; 
  out[10] += (1.0000000000000002*hamil[2]*f[10]+hamil[1]*f[6]+hamil[0]*f[2])*volFact; 
  out[11] += (hamil[2]*f[11]+hamil[1]*f[8]+hamil[0]*f[4])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *jacob_vel, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += ((2.5819888974716116*hamil[2]*f[7]+0.8606629658238706*hamil[1]*f[7]+3.872983346207417*hamil[2]*f[3]+1.2909944487358058*hamil[1]*f[3]+2.886751345948129*f[0]*hamil[2]+0.9622504486493765*f[0]*hamil[1])*dv10*volFact)/jacob_vx[2]+((-(2.5819888974716116*hamil[2]*f[7])+0.8606629658238706*hamil[1]*f[7]+3.872983346207417*hamil[2]*f[3]-1.2909944487358058*hamil[1]*f[3]-2.886751345948129*f[0]*hamil[2]+0.9622504486493765*f[0]*hamil[1])*dv10*volFact)/jacob_vx[0]+((1.5396007178390023*f[0]*hamil[1]-1.7213259316477412*hamil[1]*f[7])*dv10*volFact)/jacob_vx[1]; 
  out[2] += (2.0*hamil[2]*f[7]+2.0*hamil[1]*f[3]+2.0*f[0]*hamil[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_upper_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7200822998230954*px[3]-0.7071067811865475*px[1]+0.7071067811865475*px[0]) > v_thresh) { 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])-0.7200822998230954*px[3]) > v_thresh) { 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
GKYL_CU_DH void mom_vlasov_M0_lower_2x1v_tensor_p1(const double *w, const double *dxv, const int *idx, 
    const double *vmap, const double *jacob_vel, double v_thresh, double f_thresh, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double *px = &vmap[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_nodes[8] = {0.0};
  double fquad = 0.0;
  if ((0.7200822998230954*px[3]-0.7071067811865475*px[1]+0.7071067811865475*px[0]) < -v_thresh) { 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*(f[6]+f[5]+f[4])-0.3535533905932737*(f[3]+f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[0] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*f[6]+0.3535533905932737*f[5]-0.3535533905932737*(f[4]+f[3])+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[2] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6])-0.3535533905932737*(f[5]+f[4]+f[3]+f[2])+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[4] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]+f[5]))+0.3535533905932737*f[4]-0.3535533905932737*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
    f_nodes[6] = fquad > jacob_vx[0]*f_thresh ? fquad : 0.0; 
  } 
  if ((0.7071067811865475*(px[1]+px[0])-0.7200822998230954*px[3]) < -v_thresh) { 
    fquad = 0.3535533905932737*f[8]-0.3535533905932737*(f[6]+f[5])+0.3535533905932737*(f[4]+f[3])-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
    f_nodes[1] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*f[8])+0.3535533905932737*f[6]-0.3535533905932737*(f[5]+f[4])+0.3535533905932737*(f[3]+f[2])-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
    f_nodes[3] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = -(0.3535533905932737*(f[8]+f[6]))+0.3535533905932737*f[5]-0.3535533905932737*f[4]+0.3535533905932737*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
    f_nodes[5] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
    fquad = 0.3535533905932737*(f[8]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0]); 
    f_nodes[7] = fquad > jacob_vx[1]*f_thresh ? fquad : 0.0; 
  } 
  out[0] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
  out[1] += (0.5*f_nodes[7]+0.5*f_nodes[6]+0.5*f_nodes[5]+0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[2] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]+0.5*f_nodes[3]+0.5*f_nodes[2]-0.5*f_nodes[1]-0.5*f_nodes[0])*volFact; 
  out[3] += (0.5*f_nodes[7]+0.5*f_nodes[6]-0.5*f_nodes[5]-0.5*f_nodes[4]-0.5*f_nodes[3]-0.5*f_nodes[2]+0.5*f_nodes[1]+0.5*f_nodes[0])*volFact; 
} 
