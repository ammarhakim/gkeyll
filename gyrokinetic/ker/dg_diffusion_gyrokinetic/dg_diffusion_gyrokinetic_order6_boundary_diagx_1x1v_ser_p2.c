#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_diagx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],6.);

  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 19.062339907114627*coeff[0]*qSkin[4]-4.921875*coeff[0]*qSkin[1]; 
  boundSurf_incr[3] = 19.062339907114634*coeff[0]*qSkin[6]-4.921875*coeff[0]*qSkin[3]; 
  boundSurf_incr[4] = 19.062339907114627*coeff[0]*qSkin[1]-73.828125*coeff[0]*qSkin[4]; 
  boundSurf_incr[6] = 19.062339907114634*coeff[0]*qSkin[3]-73.828125*coeff[0]*qSkin[6]; 
  boundSurf_incr[7] = -(4.921875*coeff[0]*qSkin[7]); 

  } else { 

  boundSurf_incr[1] = -(19.062339907114627*coeff[0]*qSkin[4])-4.921875*coeff[0]*qSkin[1]; 
  boundSurf_incr[3] = -(19.062339907114634*coeff[0]*qSkin[6])-4.921875*coeff[0]*qSkin[3]; 
  boundSurf_incr[4] = -(73.828125*coeff[0]*qSkin[4])-19.062339907114627*coeff[0]*qSkin[1]; 
  boundSurf_incr[6] = -(73.828125*coeff[0]*qSkin[6])-19.062339907114634*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = -(4.921875*coeff[0]*qSkin[7]); 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 

  return 0.;
}

