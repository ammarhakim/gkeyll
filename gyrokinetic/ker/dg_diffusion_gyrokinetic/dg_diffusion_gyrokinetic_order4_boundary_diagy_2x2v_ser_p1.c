#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_diagy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],4.);

  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[2] = 1.5*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 1.5*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[1]*qSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[1]*qSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[1]*qSkin[12]; 
  boundSurf_incr[14] = 1.5*coeff[1]*qSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[1]*qSkin[15]; 
  boundSurf_incr[18] = 1.5*coeff[1]*qSkin[18]; 
  boundSurf_incr[20] = 1.5*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*qSkin[22]; 
  boundSurf_incr[23] = 1.5*coeff[1]*qSkin[23]; 

  } else { 

  boundSurf_incr[2] = 1.5*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 1.5*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 1.5*coeff[1]*qSkin[9]; 
  boundSurf_incr[11] = 1.5*coeff[1]*qSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[1]*qSkin[12]; 
  boundSurf_incr[14] = 1.5*coeff[1]*qSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[1]*qSkin[15]; 
  boundSurf_incr[18] = 1.5*coeff[1]*qSkin[18]; 
  boundSurf_incr[20] = 1.5*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 1.5*coeff[1]*qSkin[22]; 
  boundSurf_incr[23] = 1.5*coeff[1]*qSkin[23]; 

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

  return 0.;
}

