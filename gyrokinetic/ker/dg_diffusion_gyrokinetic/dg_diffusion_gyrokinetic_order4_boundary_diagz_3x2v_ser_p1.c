#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_diagz_3x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[2],4.);

  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[3] = 1.5*coeff[2]*qSkin[3]; 
  boundSurf_incr[7] = 1.5*coeff[2]*qSkin[7]; 
  boundSurf_incr[8] = 1.5*coeff[2]*qSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[2]*qSkin[11]; 
  boundSurf_incr[14] = 1.5*coeff[2]*qSkin[14]; 
  boundSurf_incr[16] = 1.5*coeff[2]*qSkin[16]; 
  boundSurf_incr[18] = 1.5*coeff[2]*qSkin[18]; 
  boundSurf_incr[19] = 1.5*coeff[2]*qSkin[19]; 
  boundSurf_incr[21] = 1.5*coeff[2]*qSkin[21]; 
  boundSurf_incr[22] = 1.5*coeff[2]*qSkin[22]; 
  boundSurf_incr[25] = 1.5*coeff[2]*qSkin[25]; 
  boundSurf_incr[26] = 1.5*coeff[2]*qSkin[26]; 
  boundSurf_incr[27] = 1.5*coeff[2]*qSkin[27]; 
  boundSurf_incr[29] = 1.5*coeff[2]*qSkin[29]; 
  boundSurf_incr[30] = 1.5*coeff[2]*qSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[2]*qSkin[31]; 
  boundSurf_incr[35] = 1.5*coeff[2]*qSkin[35]; 
  boundSurf_incr[38] = 1.5*coeff[2]*qSkin[38]; 
  boundSurf_incr[39] = 1.5*coeff[2]*qSkin[39]; 
  boundSurf_incr[42] = 1.5*coeff[2]*qSkin[42]; 
  boundSurf_incr[43] = 1.5*coeff[2]*qSkin[43]; 
  boundSurf_incr[45] = 1.5*coeff[2]*qSkin[45]; 
  boundSurf_incr[46] = 1.5*coeff[2]*qSkin[46]; 
  boundSurf_incr[47] = 1.5*coeff[2]*qSkin[47]; 

  } else { 

  boundSurf_incr[3] = 1.5*coeff[2]*qSkin[3]; 
  boundSurf_incr[7] = 1.5*coeff[2]*qSkin[7]; 
  boundSurf_incr[8] = 1.5*coeff[2]*qSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[2]*qSkin[11]; 
  boundSurf_incr[14] = 1.5*coeff[2]*qSkin[14]; 
  boundSurf_incr[16] = 1.5*coeff[2]*qSkin[16]; 
  boundSurf_incr[18] = 1.5*coeff[2]*qSkin[18]; 
  boundSurf_incr[19] = 1.5*coeff[2]*qSkin[19]; 
  boundSurf_incr[21] = 1.5*coeff[2]*qSkin[21]; 
  boundSurf_incr[22] = 1.5*coeff[2]*qSkin[22]; 
  boundSurf_incr[25] = 1.5*coeff[2]*qSkin[25]; 
  boundSurf_incr[26] = 1.5*coeff[2]*qSkin[26]; 
  boundSurf_incr[27] = 1.5*coeff[2]*qSkin[27]; 
  boundSurf_incr[29] = 1.5*coeff[2]*qSkin[29]; 
  boundSurf_incr[30] = 1.5*coeff[2]*qSkin[30]; 
  boundSurf_incr[31] = 1.5*coeff[2]*qSkin[31]; 
  boundSurf_incr[35] = 1.5*coeff[2]*qSkin[35]; 
  boundSurf_incr[38] = 1.5*coeff[2]*qSkin[38]; 
  boundSurf_incr[39] = 1.5*coeff[2]*qSkin[39]; 
  boundSurf_incr[42] = 1.5*coeff[2]*qSkin[42]; 
  boundSurf_incr[43] = 1.5*coeff[2]*qSkin[43]; 
  boundSurf_incr[45] = 1.5*coeff[2]*qSkin[45]; 
  boundSurf_incr[46] = 1.5*coeff[2]*qSkin[46]; 
  boundSurf_incr[47] = 1.5*coeff[2]*qSkin[47]; 

  }

  out[0] += -(1.0*boundSurf_incr[0]*rdx2Sq); 
  out[1] += -(1.0*boundSurf_incr[1]*rdx2Sq); 
  out[2] += -(1.0*boundSurf_incr[2]*rdx2Sq); 
  out[3] += -(1.0*boundSurf_incr[3]*rdx2Sq); 
  out[4] += -(1.0*boundSurf_incr[4]*rdx2Sq); 
  out[5] += -(1.0*boundSurf_incr[5]*rdx2Sq); 
  out[6] += -(1.0*boundSurf_incr[6]*rdx2Sq); 
  out[7] += -(1.0*boundSurf_incr[7]*rdx2Sq); 
  out[8] += -(1.0*boundSurf_incr[8]*rdx2Sq); 
  out[9] += -(1.0*boundSurf_incr[9]*rdx2Sq); 
  out[10] += -(1.0*boundSurf_incr[10]*rdx2Sq); 
  out[11] += -(1.0*boundSurf_incr[11]*rdx2Sq); 
  out[12] += -(1.0*boundSurf_incr[12]*rdx2Sq); 
  out[13] += -(1.0*boundSurf_incr[13]*rdx2Sq); 
  out[14] += -(1.0*boundSurf_incr[14]*rdx2Sq); 
  out[15] += -(1.0*boundSurf_incr[15]*rdx2Sq); 
  out[16] += -(1.0*boundSurf_incr[16]*rdx2Sq); 
  out[17] += -(1.0*boundSurf_incr[17]*rdx2Sq); 
  out[18] += -(1.0*boundSurf_incr[18]*rdx2Sq); 
  out[19] += -(1.0*boundSurf_incr[19]*rdx2Sq); 
  out[20] += -(1.0*boundSurf_incr[20]*rdx2Sq); 
  out[21] += -(1.0*boundSurf_incr[21]*rdx2Sq); 
  out[22] += -(1.0*boundSurf_incr[22]*rdx2Sq); 
  out[23] += -(1.0*boundSurf_incr[23]*rdx2Sq); 
  out[24] += -(1.0*boundSurf_incr[24]*rdx2Sq); 
  out[25] += -(1.0*boundSurf_incr[25]*rdx2Sq); 
  out[26] += -(1.0*boundSurf_incr[26]*rdx2Sq); 
  out[27] += -(1.0*boundSurf_incr[27]*rdx2Sq); 
  out[28] += -(1.0*boundSurf_incr[28]*rdx2Sq); 
  out[29] += -(1.0*boundSurf_incr[29]*rdx2Sq); 
  out[30] += -(1.0*boundSurf_incr[30]*rdx2Sq); 
  out[31] += -(1.0*boundSurf_incr[31]*rdx2Sq); 
  out[32] += -(1.0*boundSurf_incr[32]*rdx2Sq); 
  out[33] += -(1.0*boundSurf_incr[33]*rdx2Sq); 
  out[34] += -(1.0*boundSurf_incr[34]*rdx2Sq); 
  out[35] += -(1.0*boundSurf_incr[35]*rdx2Sq); 
  out[36] += -(1.0*boundSurf_incr[36]*rdx2Sq); 
  out[37] += -(1.0*boundSurf_incr[37]*rdx2Sq); 
  out[38] += -(1.0*boundSurf_incr[38]*rdx2Sq); 
  out[39] += -(1.0*boundSurf_incr[39]*rdx2Sq); 
  out[40] += -(1.0*boundSurf_incr[40]*rdx2Sq); 
  out[41] += -(1.0*boundSurf_incr[41]*rdx2Sq); 
  out[42] += -(1.0*boundSurf_incr[42]*rdx2Sq); 
  out[43] += -(1.0*boundSurf_incr[43]*rdx2Sq); 
  out[44] += -(1.0*boundSurf_incr[44]*rdx2Sq); 
  out[45] += -(1.0*boundSurf_incr[45]*rdx2Sq); 
  out[46] += -(1.0*boundSurf_incr[46]*rdx2Sq); 
  out[47] += -(1.0*boundSurf_incr[47]*rdx2Sq); 

  return 0.;
}

