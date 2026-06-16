#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*qSkin[4]-1.25*coeff[0]*qSkin[1]+0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = 0.9682458365518543*coeff[0]*qSkin[6]-1.25*coeff[0]*qSkin[3]+0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = -(3.75*coeff[0]*qSkin[4])+4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = -(3.75*coeff[0]*qSkin[6])+4.841229182759271*coeff[0]*qSkin[3]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = 0.8660254037844387*coeff[0]*qSkin[5]-1.25*coeff[0]*qSkin[7]; 

  } else { 

  boundSurf_incr[1] = -(0.9682458365518543*coeff[0]*qSkin[4])-1.25*coeff[0]*qSkin[1]-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = -(0.9682458365518543*coeff[0]*qSkin[6])-1.25*coeff[0]*qSkin[3]-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = -(3.75*coeff[0]*qSkin[4])-4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = -(3.75*coeff[0]*qSkin[6])-4.841229182759271*coeff[0]*qSkin[3]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = -(1.25*coeff[0]*qSkin[7])-0.8660254037844387*coeff[0]*qSkin[5]; 

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

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x1v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fGhost[8];
  fGhost[0] = 0.7071067811865476*(jacobgeo_inv[2]*qGhost[4]+jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qGhost[4]+qGhost[1]*jacobgeo_inv[2])+7.0710678118654755*(jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1])); 
  fGhost[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qGhost[6]+21.213203435596427*(jacobgeo_inv[1]*qGhost[3]+jacobgeo_inv[0]*qGhost[2])); 
  fGhost[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qGhost[6]+(18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qGhost[3]+21.213203435596427*jacobgeo_inv[1]*qGhost[2]); 
  fGhost[4] = 0.014285714285714285*((31.622776601683807*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qGhost[4]+49.49747468305833*qGhost[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qGhost[1]); 
  fGhost[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qGhost[7]+21.213203435596427*jacobgeo_inv[0]*qGhost[5]); 
  fGhost[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.49242404917499*jacobgeo_inv[0])*qGhost[6]+132.81566172707196*jacobgeo_inv[1]*qGhost[3]+148.49242404917499*jacobgeo_inv[2]*qGhost[2]); 
  fGhost[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qGhost[7]+21.21320343559643*jacobgeo_inv[1]*qGhost[5]); 

  double fSkin[8];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[2]*qSkin[4]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qSkin[4]+qSkin[1]*jacobgeo_inv[2])+7.0710678118654755*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qSkin[6]+21.213203435596427*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qSkin[6]+(18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qSkin[3]+21.213203435596427*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.014285714285714285*((31.622776601683807*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qSkin[4]+49.49747468305833*qSkin[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[7]+21.213203435596427*jacobgeo_inv[0]*qSkin[5]); 
  fSkin[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.49242404917499*jacobgeo_inv[0])*qSkin[6]+132.81566172707196*jacobgeo_inv[1]*qSkin[3]+148.49242404917499*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qSkin[7]+21.21320343559643*jacobgeo_inv[1]*qSkin[5]); 

  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 1.5309310892394856*coeff[2]*fSkin[4]-1.185854122563142*coeff[1]*fSkin[4]+0.6846531968814573*coeff[0]*fSkin[4]-1.9764235376052366*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 1.5309310892394858*coeff[2]*fSkin[6]-1.1858541225631418*coeff[1]*fSkin[6]+0.6846531968814574*coeff[0]*fSkin[6]-1.9764235376052366*coeff[2]*fSkin[3]+1.5309310892394856*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]+1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = -(5.929270612815711*coeff[2]*fSkin[4])+4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.4232659844072875*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(5.929270612815711*coeff[2]*fSkin[6])+4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]+7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]+3.4232659844072866*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = -(1.9764235376052366*coeff[2]*fSkin[7])+1.5309310892394856*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]+1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]+0.6123724356957944*coeff[0]*fSkin[5]; 

  } else { 

  boundSurf_incr[1] = -(1.5309310892394856*coeff[2]*fSkin[4])-1.185854122563142*coeff[1]*fSkin[4]-0.6846531968814573*coeff[0]*fSkin[4]-1.9764235376052366*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(1.5309310892394858*coeff[2]*fSkin[6])-1.1858541225631418*coeff[1]*fSkin[6]-0.6846531968814574*coeff[0]*fSkin[6]-1.9764235376052366*coeff[2]*fSkin[3]-1.5309310892394856*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]-1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = -(5.929270612815711*coeff[2]*fSkin[4])-4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.4232659844072875*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(5.929270612815711*coeff[2]*fSkin[6])-4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]-7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]-3.4232659844072866*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = -(1.9764235376052366*coeff[2]*fSkin[7])-1.5309310892394856*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]-1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]-0.6123724356957944*coeff[0]*fSkin[5]; 

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

