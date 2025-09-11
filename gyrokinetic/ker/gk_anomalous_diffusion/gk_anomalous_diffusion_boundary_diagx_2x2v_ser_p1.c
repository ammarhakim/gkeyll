#include <gkyl_gk_anomalous_diffusion_kernels.h>

GKYL_CU_DH double gk_anomalous_diffusion_boundary_diagx_2x2v_ser_p1(const double *wSkin, const double *dxSkin, const double *nuSkin, const double *nuGhost, const double *jacobgeo_inv, int edge, const double *JfSkin, const double *JfGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // nuSkin/nuGhost: Diffusivity in skin and ghost cells.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // JfSkin/JfGhost: distribution times conf-space Jacobian in skin and ghost cells.
  // out: Incremented output.

  const double rdx2sq = pow(2./dxSkin[0],2.0);

  double fSkin[24];
  fSkin[0] = 0.5*(jacobgeo_inv[3]*JfSkin[5]+JfSkin[2]*jacobgeo_inv[2]+JfSkin[1]*jacobgeo_inv[1]+JfSkin[0]*jacobgeo_inv[0]); 
  fSkin[1] = 0.5*(jacobgeo_inv[2]*JfSkin[5]+JfSkin[2]*jacobgeo_inv[3]+JfSkin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*JfSkin[1]); 
  fSkin[2] = 0.5*(jacobgeo_inv[1]*JfSkin[5]+JfSkin[1]*jacobgeo_inv[3]+JfSkin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*JfSkin[2]); 
  fSkin[3] = 0.5*(jacobgeo_inv[3]*JfSkin[11]+jacobgeo_inv[2]*JfSkin[7]+jacobgeo_inv[1]*JfSkin[6]+jacobgeo_inv[0]*JfSkin[3]); 
  fSkin[4] = 0.5*(jacobgeo_inv[3]*JfSkin[12]+jacobgeo_inv[2]*JfSkin[9]+jacobgeo_inv[1]*JfSkin[8]+jacobgeo_inv[0]*JfSkin[4]); 
  fSkin[5] = 0.5*(jacobgeo_inv[0]*JfSkin[5]+JfSkin[0]*jacobgeo_inv[3]+JfSkin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*JfSkin[2]); 
  fSkin[6] = 0.5*(jacobgeo_inv[2]*JfSkin[11]+jacobgeo_inv[3]*JfSkin[7]+jacobgeo_inv[0]*JfSkin[6]+jacobgeo_inv[1]*JfSkin[3]); 
  fSkin[7] = 0.5*(jacobgeo_inv[1]*JfSkin[11]+jacobgeo_inv[0]*JfSkin[7]+jacobgeo_inv[3]*JfSkin[6]+jacobgeo_inv[2]*JfSkin[3]); 
  fSkin[8] = 0.5*(jacobgeo_inv[2]*JfSkin[12]+jacobgeo_inv[3]*JfSkin[9]+jacobgeo_inv[0]*JfSkin[8]+jacobgeo_inv[1]*JfSkin[4]); 
  fSkin[9] = 0.5*(jacobgeo_inv[1]*JfSkin[12]+jacobgeo_inv[0]*JfSkin[9]+jacobgeo_inv[3]*JfSkin[8]+jacobgeo_inv[2]*JfSkin[4]); 
  fSkin[10] = 0.5*(jacobgeo_inv[3]*JfSkin[15]+jacobgeo_inv[2]*JfSkin[14]+jacobgeo_inv[1]*JfSkin[13]+jacobgeo_inv[0]*JfSkin[10]); 
  fSkin[11] = 0.5*(jacobgeo_inv[0]*JfSkin[11]+jacobgeo_inv[1]*JfSkin[7]+jacobgeo_inv[2]*JfSkin[6]+JfSkin[3]*jacobgeo_inv[3]); 
  fSkin[12] = 0.5*(jacobgeo_inv[0]*JfSkin[12]+jacobgeo_inv[1]*JfSkin[9]+jacobgeo_inv[2]*JfSkin[8]+jacobgeo_inv[3]*JfSkin[4]); 
  fSkin[13] = 0.5*(jacobgeo_inv[2]*JfSkin[15]+jacobgeo_inv[3]*JfSkin[14]+jacobgeo_inv[0]*JfSkin[13]+jacobgeo_inv[1]*JfSkin[10]); 
  fSkin[14] = 0.5*(jacobgeo_inv[1]*JfSkin[15]+jacobgeo_inv[0]*JfSkin[14]+jacobgeo_inv[3]*JfSkin[13]+jacobgeo_inv[2]*JfSkin[10]); 
  fSkin[15] = 0.5*(jacobgeo_inv[0]*JfSkin[15]+jacobgeo_inv[1]*JfSkin[14]+jacobgeo_inv[2]*JfSkin[13]+jacobgeo_inv[3]*JfSkin[10]); 
  fSkin[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*JfSkin[20]+15.000000000000002*(jacobgeo_inv[2]*JfSkin[18]+jacobgeo_inv[1]*JfSkin[17])+15.0*jacobgeo_inv[0]*JfSkin[16]); 
  fSkin[17] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*JfSkin[20]+15.0*(jacobgeo_inv[3]*JfSkin[18]+jacobgeo_inv[0]*JfSkin[17])+15.000000000000002*jacobgeo_inv[1]*JfSkin[16]); 
  fSkin[18] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*JfSkin[20]+15.0*(jacobgeo_inv[0]*JfSkin[18]+jacobgeo_inv[3]*JfSkin[17])+15.000000000000002*jacobgeo_inv[2]*JfSkin[16]); 
  fSkin[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*JfSkin[23]+15.000000000000002*(jacobgeo_inv[2]*JfSkin[22]+jacobgeo_inv[1]*JfSkin[21])+15.0*jacobgeo_inv[0]*JfSkin[19]); 
  fSkin[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*JfSkin[20]+15.000000000000002*(jacobgeo_inv[1]*JfSkin[18]+jacobgeo_inv[2]*JfSkin[17])+15.0*jacobgeo_inv[3]*JfSkin[16]); 
  fSkin[21] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[2]*JfSkin[23]+15.0*(jacobgeo_inv[3]*JfSkin[22]+jacobgeo_inv[0]*JfSkin[21])+15.000000000000002*jacobgeo_inv[1]*JfSkin[19]); 
  fSkin[22] = 0.03333333333333333*(15.000000000000002*jacobgeo_inv[1]*JfSkin[23]+15.0*(jacobgeo_inv[0]*JfSkin[22]+jacobgeo_inv[3]*JfSkin[21])+15.000000000000002*jacobgeo_inv[2]*JfSkin[19]); 
  fSkin[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*JfSkin[23]+15.000000000000002*(jacobgeo_inv[1]*JfSkin[22]+jacobgeo_inv[2]*JfSkin[21])+15.0*jacobgeo_inv[3]*JfSkin[19]); 

  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[0] = -(0.75*nuSkin[3]*fSkin[5])+0.4330127018922193*nuSkin[2]*fSkin[5]-0.75*fSkin[1]*nuSkin[1]+0.4330127018922193*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = 1.2990381056766578*nuSkin[3]*fSkin[5]-0.75*nuSkin[2]*fSkin[5]+1.2990381056766578*fSkin[1]*nuSkin[1]-0.75*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = -(0.75*nuSkin[1]*fSkin[5])+0.4330127018922193*nuSkin[0]*fSkin[5]-0.75*fSkin[1]*nuSkin[3]+0.4330127018922193*fSkin[1]*nuSkin[2]; 
  boundSurf_incr[3] = -(0.75*nuSkin[3]*fSkin[11])+0.4330127018922193*nuSkin[2]*fSkin[11]-0.75*nuSkin[1]*fSkin[6]+0.4330127018922193*nuSkin[0]*fSkin[6]; 
  boundSurf_incr[4] = -(0.75*nuSkin[3]*fSkin[12])+0.4330127018922193*nuSkin[2]*fSkin[12]-0.75*nuSkin[1]*fSkin[8]+0.4330127018922193*nuSkin[0]*fSkin[8]; 
  boundSurf_incr[5] = 1.2990381056766578*nuSkin[1]*fSkin[5]-0.75*nuSkin[0]*fSkin[5]+1.2990381056766578*fSkin[1]*nuSkin[3]-0.75*fSkin[1]*nuSkin[2]; 
  boundSurf_incr[6] = 1.2990381056766578*nuSkin[3]*fSkin[11]-0.75*nuSkin[2]*fSkin[11]+1.2990381056766578*nuSkin[1]*fSkin[6]-0.75*nuSkin[0]*fSkin[6]; 
  boundSurf_incr[7] = -(0.75*nuSkin[1]*fSkin[11])+0.4330127018922193*nuSkin[0]*fSkin[11]-0.75*nuSkin[3]*fSkin[6]+0.4330127018922193*nuSkin[2]*fSkin[6]; 
  boundSurf_incr[8] = 1.2990381056766578*nuSkin[3]*fSkin[12]-0.75*nuSkin[2]*fSkin[12]+1.2990381056766578*nuSkin[1]*fSkin[8]-0.75*nuSkin[0]*fSkin[8]; 
  boundSurf_incr[9] = -(0.75*nuSkin[1]*fSkin[12])+0.4330127018922193*nuSkin[0]*fSkin[12]-0.75*nuSkin[3]*fSkin[8]+0.4330127018922193*nuSkin[2]*fSkin[8]; 
  boundSurf_incr[10] = -(0.75*nuSkin[3]*fSkin[15])+0.4330127018922193*nuSkin[2]*fSkin[15]-0.75*nuSkin[1]*fSkin[13]+0.4330127018922193*nuSkin[0]*fSkin[13]; 
  boundSurf_incr[11] = 1.2990381056766578*nuSkin[1]*fSkin[11]-0.75*nuSkin[0]*fSkin[11]+1.2990381056766578*nuSkin[3]*fSkin[6]-0.75*nuSkin[2]*fSkin[6]; 
  boundSurf_incr[12] = 1.2990381056766578*nuSkin[1]*fSkin[12]-0.75*nuSkin[0]*fSkin[12]+1.2990381056766578*nuSkin[3]*fSkin[8]-0.75*nuSkin[2]*fSkin[8]; 
  boundSurf_incr[13] = 1.2990381056766578*nuSkin[3]*fSkin[15]-0.75*nuSkin[2]*fSkin[15]+1.2990381056766578*nuSkin[1]*fSkin[13]-0.75*nuSkin[0]*fSkin[13]; 
  boundSurf_incr[14] = -(0.75*nuSkin[1]*fSkin[15])+0.4330127018922193*nuSkin[0]*fSkin[15]-0.75*nuSkin[3]*fSkin[13]+0.4330127018922193*nuSkin[2]*fSkin[13]; 
  boundSurf_incr[15] = 1.2990381056766578*nuSkin[1]*fSkin[15]-0.75*nuSkin[0]*fSkin[15]+1.2990381056766578*nuSkin[3]*fSkin[13]-0.75*nuSkin[2]*fSkin[13]; 
  boundSurf_incr[16] = -(0.75*nuSkin[3]*fSkin[20])+0.4330127018922193*nuSkin[2]*fSkin[20]-0.75*nuSkin[1]*fSkin[17]+0.43301270189221935*nuSkin[0]*fSkin[17]; 
  boundSurf_incr[17] = 1.299038105676658*nuSkin[3]*fSkin[20]-0.75*nuSkin[2]*fSkin[20]+1.2990381056766578*nuSkin[1]*fSkin[17]-0.75*nuSkin[0]*fSkin[17]; 
  boundSurf_incr[18] = -(0.75*nuSkin[1]*fSkin[20])+0.43301270189221935*nuSkin[0]*fSkin[20]-0.75*nuSkin[3]*fSkin[17]+0.4330127018922193*nuSkin[2]*fSkin[17]; 
  boundSurf_incr[19] = -(0.75*nuSkin[3]*fSkin[23])+0.4330127018922193*nuSkin[2]*fSkin[23]-0.75*nuSkin[1]*fSkin[21]+0.43301270189221935*nuSkin[0]*fSkin[21]; 
  boundSurf_incr[20] = 1.2990381056766578*nuSkin[1]*fSkin[20]-0.75*nuSkin[0]*fSkin[20]+1.299038105676658*nuSkin[3]*fSkin[17]-0.75*nuSkin[2]*fSkin[17]; 
  boundSurf_incr[21] = 1.299038105676658*nuSkin[3]*fSkin[23]-0.75*nuSkin[2]*fSkin[23]+1.2990381056766578*nuSkin[1]*fSkin[21]-0.75*nuSkin[0]*fSkin[21]; 
  boundSurf_incr[22] = -(0.75*nuSkin[1]*fSkin[23])+0.43301270189221935*nuSkin[0]*fSkin[23]-0.75*nuSkin[3]*fSkin[21]+0.4330127018922193*nuSkin[2]*fSkin[21]; 
  boundSurf_incr[23] = 1.2990381056766578*nuSkin[1]*fSkin[23]-0.75*nuSkin[0]*fSkin[23]+1.299038105676658*nuSkin[3]*fSkin[21]-0.75*nuSkin[2]*fSkin[21]; 

  } else { 

  boundSurf_incr[0] = -(0.75*nuSkin[3]*fSkin[5])-0.4330127018922193*nuSkin[2]*fSkin[5]-0.75*fSkin[1]*nuSkin[1]-0.4330127018922193*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[1] = -(1.2990381056766578*nuSkin[3]*fSkin[5])-0.75*nuSkin[2]*fSkin[5]-1.2990381056766578*fSkin[1]*nuSkin[1]-0.75*nuSkin[0]*fSkin[1]; 
  boundSurf_incr[2] = -(0.75*nuSkin[1]*fSkin[5])-0.4330127018922193*nuSkin[0]*fSkin[5]-0.75*fSkin[1]*nuSkin[3]-0.4330127018922193*fSkin[1]*nuSkin[2]; 
  boundSurf_incr[3] = -(0.75*nuSkin[3]*fSkin[11])-0.4330127018922193*nuSkin[2]*fSkin[11]-0.75*nuSkin[1]*fSkin[6]-0.4330127018922193*nuSkin[0]*fSkin[6]; 
  boundSurf_incr[4] = -(0.75*nuSkin[3]*fSkin[12])-0.4330127018922193*nuSkin[2]*fSkin[12]-0.75*nuSkin[1]*fSkin[8]-0.4330127018922193*nuSkin[0]*fSkin[8]; 
  boundSurf_incr[5] = -(1.2990381056766578*nuSkin[1]*fSkin[5])-0.75*nuSkin[0]*fSkin[5]-1.2990381056766578*fSkin[1]*nuSkin[3]-0.75*fSkin[1]*nuSkin[2]; 
  boundSurf_incr[6] = -(1.2990381056766578*nuSkin[3]*fSkin[11])-0.75*nuSkin[2]*fSkin[11]-1.2990381056766578*nuSkin[1]*fSkin[6]-0.75*nuSkin[0]*fSkin[6]; 
  boundSurf_incr[7] = -(0.75*nuSkin[1]*fSkin[11])-0.4330127018922193*nuSkin[0]*fSkin[11]-0.75*nuSkin[3]*fSkin[6]-0.4330127018922193*nuSkin[2]*fSkin[6]; 
  boundSurf_incr[8] = -(1.2990381056766578*nuSkin[3]*fSkin[12])-0.75*nuSkin[2]*fSkin[12]-1.2990381056766578*nuSkin[1]*fSkin[8]-0.75*nuSkin[0]*fSkin[8]; 
  boundSurf_incr[9] = -(0.75*nuSkin[1]*fSkin[12])-0.4330127018922193*nuSkin[0]*fSkin[12]-0.75*nuSkin[3]*fSkin[8]-0.4330127018922193*nuSkin[2]*fSkin[8]; 
  boundSurf_incr[10] = -(0.75*nuSkin[3]*fSkin[15])-0.4330127018922193*nuSkin[2]*fSkin[15]-0.75*nuSkin[1]*fSkin[13]-0.4330127018922193*nuSkin[0]*fSkin[13]; 
  boundSurf_incr[11] = -(1.2990381056766578*nuSkin[1]*fSkin[11])-0.75*nuSkin[0]*fSkin[11]-1.2990381056766578*nuSkin[3]*fSkin[6]-0.75*nuSkin[2]*fSkin[6]; 
  boundSurf_incr[12] = -(1.2990381056766578*nuSkin[1]*fSkin[12])-0.75*nuSkin[0]*fSkin[12]-1.2990381056766578*nuSkin[3]*fSkin[8]-0.75*nuSkin[2]*fSkin[8]; 
  boundSurf_incr[13] = -(1.2990381056766578*nuSkin[3]*fSkin[15])-0.75*nuSkin[2]*fSkin[15]-1.2990381056766578*nuSkin[1]*fSkin[13]-0.75*nuSkin[0]*fSkin[13]; 
  boundSurf_incr[14] = -(0.75*nuSkin[1]*fSkin[15])-0.4330127018922193*nuSkin[0]*fSkin[15]-0.75*nuSkin[3]*fSkin[13]-0.4330127018922193*nuSkin[2]*fSkin[13]; 
  boundSurf_incr[15] = -(1.2990381056766578*nuSkin[1]*fSkin[15])-0.75*nuSkin[0]*fSkin[15]-1.2990381056766578*nuSkin[3]*fSkin[13]-0.75*nuSkin[2]*fSkin[13]; 
  boundSurf_incr[16] = -(0.75*nuSkin[3]*fSkin[20])-0.4330127018922193*nuSkin[2]*fSkin[20]-0.75*nuSkin[1]*fSkin[17]-0.43301270189221935*nuSkin[0]*fSkin[17]; 
  boundSurf_incr[17] = -(1.299038105676658*nuSkin[3]*fSkin[20])-0.75*nuSkin[2]*fSkin[20]-1.2990381056766578*nuSkin[1]*fSkin[17]-0.75*nuSkin[0]*fSkin[17]; 
  boundSurf_incr[18] = -(0.75*nuSkin[1]*fSkin[20])-0.43301270189221935*nuSkin[0]*fSkin[20]-0.75*nuSkin[3]*fSkin[17]-0.4330127018922193*nuSkin[2]*fSkin[17]; 
  boundSurf_incr[19] = -(0.75*nuSkin[3]*fSkin[23])-0.4330127018922193*nuSkin[2]*fSkin[23]-0.75*nuSkin[1]*fSkin[21]-0.43301270189221935*nuSkin[0]*fSkin[21]; 
  boundSurf_incr[20] = -(1.2990381056766578*nuSkin[1]*fSkin[20])-0.75*nuSkin[0]*fSkin[20]-1.299038105676658*nuSkin[3]*fSkin[17]-0.75*nuSkin[2]*fSkin[17]; 
  boundSurf_incr[21] = -(1.299038105676658*nuSkin[3]*fSkin[23])-0.75*nuSkin[2]*fSkin[23]-1.2990381056766578*nuSkin[1]*fSkin[21]-0.75*nuSkin[0]*fSkin[21]; 
  boundSurf_incr[22] = -(0.75*nuSkin[1]*fSkin[23])-0.43301270189221935*nuSkin[0]*fSkin[23]-0.75*nuSkin[3]*fSkin[21]-0.4330127018922193*nuSkin[2]*fSkin[21]; 
  boundSurf_incr[23] = -(1.2990381056766578*nuSkin[1]*fSkin[23])-0.75*nuSkin[0]*fSkin[23]-1.299038105676658*nuSkin[3]*fSkin[21]-0.75*nuSkin[2]*fSkin[21]; 

  }

  out[0] += boundSurf_incr[0]*rdx2sq; 
  out[1] += boundSurf_incr[1]*rdx2sq; 
  out[2] += boundSurf_incr[2]*rdx2sq; 
  out[3] += boundSurf_incr[3]*rdx2sq; 
  out[4] += boundSurf_incr[4]*rdx2sq; 
  out[5] += boundSurf_incr[5]*rdx2sq; 
  out[6] += boundSurf_incr[6]*rdx2sq; 
  out[7] += boundSurf_incr[7]*rdx2sq; 
  out[8] += boundSurf_incr[8]*rdx2sq; 
  out[9] += boundSurf_incr[9]*rdx2sq; 
  out[10] += boundSurf_incr[10]*rdx2sq; 
  out[11] += boundSurf_incr[11]*rdx2sq; 
  out[12] += boundSurf_incr[12]*rdx2sq; 
  out[13] += boundSurf_incr[13]*rdx2sq; 
  out[14] += boundSurf_incr[14]*rdx2sq; 
  out[15] += boundSurf_incr[15]*rdx2sq; 
  out[16] += boundSurf_incr[16]*rdx2sq; 
  out[17] += boundSurf_incr[17]*rdx2sq; 
  out[18] += boundSurf_incr[18]*rdx2sq; 
  out[19] += boundSurf_incr[19]*rdx2sq; 
  out[20] += boundSurf_incr[20]*rdx2sq; 
  out[21] += boundSurf_incr[21]*rdx2sq; 
  out[22] += boundSurf_incr[22]*rdx2sq; 
  out[23] += boundSurf_incr[23]*rdx2sq; 

  return 0.;
}

