#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_diagy_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],6.);

  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[2] = 19.062339907114627*coeff[1]*qSkin[12]-4.921875*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 19.062339907114634*coeff[1]*qSkin[20]-4.921875*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 19.062339907114634*coeff[1]*qSkin[22]-4.921875*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 19.062339907114634*coeff[1]*qSkin[26]-4.921875*coeff[1]*qSkin[9]; 
  boundSurf_incr[12] = 19.062339907114627*coeff[1]*qSkin[2]-73.828125*coeff[1]*qSkin[12]; 
  boundSurf_incr[15] = 19.062339907114627*coeff[1]*qSkin[33]-4.921875*coeff[1]*qSkin[15]; 
  boundSurf_incr[16] = 19.062339907114627*coeff[1]*qSkin[36]-4.921875*coeff[1]*qSkin[16]; 
  boundSurf_incr[18] = 19.062339907114627*coeff[1]*qSkin[38]-4.921875*coeff[1]*qSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*qSkin[19]); 
  boundSurf_incr[20] = 19.062339907114634*coeff[1]*qSkin[5]-73.828125*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 19.062339907114634*coeff[1]*qSkin[7]-73.828125*coeff[1]*qSkin[22]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*qSkin[24]); 
  boundSurf_incr[26] = 19.062339907114634*coeff[1]*qSkin[9]-73.828125*coeff[1]*qSkin[26]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*qSkin[29]); 
  boundSurf_incr[31] = 19.062339907114634*coeff[1]*qSkin[45]-4.921875*coeff[1]*qSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*qSkin[32]); 
  boundSurf_incr[33] = 19.062339907114627*coeff[1]*qSkin[15]-73.828125*coeff[1]*qSkin[33]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*qSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*qSkin[35]); 
  boundSurf_incr[36] = 19.062339907114627*coeff[1]*qSkin[16]-73.828125*coeff[1]*qSkin[36]; 
  boundSurf_incr[38] = 19.062339907114627*coeff[1]*qSkin[18]-73.828125*coeff[1]*qSkin[38]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*qSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*qSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*qSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*qSkin[44]); 
  boundSurf_incr[45] = 19.062339907114634*coeff[1]*qSkin[31]-73.828125*coeff[1]*qSkin[45]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*qSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*qSkin[47]); 

  } else { 

  boundSurf_incr[2] = -(19.062339907114627*coeff[1]*qSkin[12])-4.921875*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = -(19.062339907114634*coeff[1]*qSkin[20])-4.921875*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = -(19.062339907114634*coeff[1]*qSkin[22])-4.921875*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = -(19.062339907114634*coeff[1]*qSkin[26])-4.921875*coeff[1]*qSkin[9]; 
  boundSurf_incr[12] = -(73.828125*coeff[1]*qSkin[12])-19.062339907114627*coeff[1]*qSkin[2]; 
  boundSurf_incr[15] = -(19.062339907114627*coeff[1]*qSkin[33])-4.921875*coeff[1]*qSkin[15]; 
  boundSurf_incr[16] = -(19.062339907114627*coeff[1]*qSkin[36])-4.921875*coeff[1]*qSkin[16]; 
  boundSurf_incr[18] = -(19.062339907114627*coeff[1]*qSkin[38])-4.921875*coeff[1]*qSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*qSkin[19]); 
  boundSurf_incr[20] = -(73.828125*coeff[1]*qSkin[20])-19.062339907114634*coeff[1]*qSkin[5]; 
  boundSurf_incr[22] = -(73.828125*coeff[1]*qSkin[22])-19.062339907114634*coeff[1]*qSkin[7]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*qSkin[24]); 
  boundSurf_incr[26] = -(73.828125*coeff[1]*qSkin[26])-19.062339907114634*coeff[1]*qSkin[9]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*qSkin[29]); 
  boundSurf_incr[31] = -(19.062339907114634*coeff[1]*qSkin[45])-4.921875*coeff[1]*qSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*qSkin[32]); 
  boundSurf_incr[33] = -(73.828125*coeff[1]*qSkin[33])-19.062339907114627*coeff[1]*qSkin[15]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*qSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*qSkin[35]); 
  boundSurf_incr[36] = -(73.828125*coeff[1]*qSkin[36])-19.062339907114627*coeff[1]*qSkin[16]; 
  boundSurf_incr[38] = -(73.828125*coeff[1]*qSkin[38])-19.062339907114627*coeff[1]*qSkin[18]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*qSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*qSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*qSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*qSkin[44]); 
  boundSurf_incr[45] = -(73.828125*coeff[1]*qSkin[45])-19.062339907114634*coeff[1]*qSkin[31]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*qSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*qSkin[47]); 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 
  out[12] += boundSurf_incr[12]*rdx2Sq; 
  out[13] += boundSurf_incr[13]*rdx2Sq; 
  out[14] += boundSurf_incr[14]*rdx2Sq; 
  out[15] += boundSurf_incr[15]*rdx2Sq; 
  out[16] += boundSurf_incr[16]*rdx2Sq; 
  out[17] += boundSurf_incr[17]*rdx2Sq; 
  out[18] += boundSurf_incr[18]*rdx2Sq; 
  out[19] += boundSurf_incr[19]*rdx2Sq; 
  out[20] += boundSurf_incr[20]*rdx2Sq; 
  out[21] += boundSurf_incr[21]*rdx2Sq; 
  out[22] += boundSurf_incr[22]*rdx2Sq; 
  out[23] += boundSurf_incr[23]*rdx2Sq; 
  out[24] += boundSurf_incr[24]*rdx2Sq; 
  out[25] += boundSurf_incr[25]*rdx2Sq; 
  out[26] += boundSurf_incr[26]*rdx2Sq; 
  out[27] += boundSurf_incr[27]*rdx2Sq; 
  out[28] += boundSurf_incr[28]*rdx2Sq; 
  out[29] += boundSurf_incr[29]*rdx2Sq; 
  out[30] += boundSurf_incr[30]*rdx2Sq; 
  out[31] += boundSurf_incr[31]*rdx2Sq; 
  out[32] += boundSurf_incr[32]*rdx2Sq; 
  out[33] += boundSurf_incr[33]*rdx2Sq; 
  out[34] += boundSurf_incr[34]*rdx2Sq; 
  out[35] += boundSurf_incr[35]*rdx2Sq; 
  out[36] += boundSurf_incr[36]*rdx2Sq; 
  out[37] += boundSurf_incr[37]*rdx2Sq; 
  out[38] += boundSurf_incr[38]*rdx2Sq; 
  out[39] += boundSurf_incr[39]*rdx2Sq; 
  out[40] += boundSurf_incr[40]*rdx2Sq; 
  out[41] += boundSurf_incr[41]*rdx2Sq; 
  out[42] += boundSurf_incr[42]*rdx2Sq; 
  out[43] += boundSurf_incr[43]*rdx2Sq; 
  out[44] += boundSurf_incr[44]*rdx2Sq; 
  out[45] += boundSurf_incr[45]*rdx2Sq; 
  out[46] += boundSurf_incr[46]*rdx2Sq; 
  out[47] += boundSurf_incr[47]*rdx2Sq; 

  return 0.;
}

