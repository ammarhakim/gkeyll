#include <gkyl_gk_numerical_diffusion_kernels.h>

static GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffx(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[0],2); 
  out[1] += -(1.060660171779821*nu[6]*fin[16]*rdx2fac)-1.060660171779821*nu[3]*fin[7]*rdx2fac-1.060660171779821*nu[2]*fin[6]*rdx2fac-1.060660171779821*nu[0]*fin[1]*rdx2fac; 
  out[6] += -(1.060660171779821*nu[3]*fin[16]*rdx2fac)-1.060660171779821*nu[6]*fin[7]*rdx2fac-1.060660171779821*nu[0]*fin[6]*rdx2fac-1.060660171779821*fin[1]*nu[2]*rdx2fac; 
  out[7] += -(1.060660171779821*nu[2]*fin[16]*rdx2fac)-1.060660171779821*nu[0]*fin[7]*rdx2fac-1.060660171779821*fin[6]*nu[6]*rdx2fac-1.060660171779821*fin[1]*nu[3]*rdx2fac; 
  out[9] += -(1.060660171779821*nu[6]*fin[26]*rdx2fac)-1.060660171779821*nu[3]*fin[18]*rdx2fac-1.060660171779821*nu[2]*fin[17]*rdx2fac-1.060660171779821*nu[0]*fin[9]*rdx2fac; 
  out[12] += -(1.060660171779821*nu[6]*fin[27]*rdx2fac)-1.060660171779821*nu[3]*fin[21]*rdx2fac-1.060660171779821*nu[2]*fin[20]*rdx2fac-1.060660171779821*nu[0]*fin[12]*rdx2fac; 
  out[16] += -(1.060660171779821*nu[0]*fin[16]*rdx2fac)-1.060660171779821*nu[2]*fin[7]*rdx2fac-1.060660171779821*fin[1]*nu[6]*rdx2fac-1.060660171779821*nu[3]*fin[6]*rdx2fac; 
  out[17] += -(1.060660171779821*nu[3]*fin[26]*rdx2fac)-1.060660171779821*nu[6]*fin[18]*rdx2fac-1.060660171779821*nu[0]*fin[17]*rdx2fac-1.060660171779821*nu[2]*fin[9]*rdx2fac; 
  out[18] += -(1.060660171779821*nu[2]*fin[26]*rdx2fac)-1.060660171779821*nu[0]*fin[18]*rdx2fac-1.060660171779821*nu[6]*fin[17]*rdx2fac-1.060660171779821*nu[3]*fin[9]*rdx2fac; 
  out[20] += -(1.060660171779821*nu[3]*fin[27]*rdx2fac)-1.060660171779821*nu[6]*fin[21]*rdx2fac-1.060660171779821*nu[0]*fin[20]*rdx2fac-1.060660171779821*nu[2]*fin[12]*rdx2fac; 
  out[21] += -(1.060660171779821*nu[2]*fin[27]*rdx2fac)-1.060660171779821*nu[0]*fin[21]*rdx2fac-1.060660171779821*nu[6]*fin[20]*rdx2fac-1.060660171779821*nu[3]*fin[12]*rdx2fac; 
  out[23] += -(1.060660171779821*nu[6]*fin[31]*rdx2fac)-1.060660171779821*nu[3]*fin[29]*rdx2fac-1.060660171779821*nu[2]*fin[28]*rdx2fac-1.060660171779821*nu[0]*fin[23]*rdx2fac; 
  out[26] += -(1.060660171779821*nu[0]*fin[26]*rdx2fac)-1.060660171779821*nu[2]*fin[18]*rdx2fac-1.060660171779821*nu[3]*fin[17]*rdx2fac-1.060660171779821*nu[6]*fin[9]*rdx2fac; 
  out[27] += -(1.060660171779821*nu[0]*fin[27]*rdx2fac)-1.060660171779821*nu[2]*fin[21]*rdx2fac-1.060660171779821*nu[3]*fin[20]*rdx2fac-1.060660171779821*nu[6]*fin[12]*rdx2fac; 
  out[28] += -(1.060660171779821*nu[3]*fin[31]*rdx2fac)-1.060660171779821*nu[6]*fin[29]*rdx2fac-1.060660171779821*nu[0]*fin[28]*rdx2fac-1.060660171779821*nu[2]*fin[23]*rdx2fac; 
  out[29] += -(1.060660171779821*nu[2]*fin[31]*rdx2fac)-1.060660171779821*nu[0]*fin[29]*rdx2fac-1.060660171779821*nu[6]*fin[28]*rdx2fac-1.060660171779821*nu[3]*fin[23]*rdx2fac; 
  out[31] += -(1.060660171779821*nu[0]*fin[31]*rdx2fac)-1.060660171779821*nu[2]*fin[29]*rdx2fac-1.060660171779821*nu[3]*fin[28]*rdx2fac-1.060660171779821*nu[6]*fin[23]*rdx2fac; 
  out[33] += -(1.060660171779821*nu[6]*fin[43]*rdx2fac)-1.060660171779821*nu[3]*fin[38]*rdx2fac-1.060660171779821*nu[2]*fin[37]*rdx2fac-1.060660171779821*nu[0]*fin[33]*rdx2fac; 
  out[37] += -(1.060660171779821*nu[3]*fin[43]*rdx2fac)-1.060660171779821*nu[6]*fin[38]*rdx2fac-1.060660171779821*nu[0]*fin[37]*rdx2fac-1.060660171779821*nu[2]*fin[33]*rdx2fac; 
  out[38] += -(1.060660171779821*nu[2]*fin[43]*rdx2fac)-1.060660171779821*nu[0]*fin[38]*rdx2fac-1.060660171779821*nu[6]*fin[37]*rdx2fac-1.060660171779821*nu[3]*fin[33]*rdx2fac; 
  out[40] += -(1.060660171779821*nu[6]*fin[47]*rdx2fac)-1.060660171779821*nu[3]*fin[45]*rdx2fac-1.060660171779821*nu[2]*fin[44]*rdx2fac-1.060660171779821*nu[0]*fin[40]*rdx2fac; 
  out[43] += -(1.060660171779821*nu[0]*fin[43]*rdx2fac)-1.060660171779821*nu[2]*fin[38]*rdx2fac-1.060660171779821*nu[3]*fin[37]*rdx2fac-1.060660171779821*nu[6]*fin[33]*rdx2fac; 
  out[44] += -(1.060660171779821*nu[3]*fin[47]*rdx2fac)-1.060660171779821*nu[6]*fin[45]*rdx2fac-1.060660171779821*nu[0]*fin[44]*rdx2fac-1.060660171779821*nu[2]*fin[40]*rdx2fac; 
  out[45] += -(1.060660171779821*nu[2]*fin[47]*rdx2fac)-1.060660171779821*nu[0]*fin[45]*rdx2fac-1.060660171779821*nu[6]*fin[44]*rdx2fac-1.060660171779821*nu[3]*fin[40]*rdx2fac; 
  out[47] += -(1.060660171779821*nu[0]*fin[47]*rdx2fac)-1.060660171779821*nu[2]*fin[45]*rdx2fac-1.060660171779821*nu[3]*fin[44]*rdx2fac-1.060660171779821*nu[6]*fin[40]*rdx2fac; 

  return 1.4142135623730951*nu[0]*rdx2fac; 
}

static GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffy(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[1],2); 
  out[2] += -(1.060660171779821*nu[13]*fin[16]*rdx2fac)-1.060660171779821*fin[8]*nu[11]*rdx2fac-1.060660171779821*fin[6]*nu[9]*rdx2fac-1.060660171779821*fin[2]*nu[8]*rdx2fac; 
  out[6] += -(1.060660171779821*nu[11]*fin[16]*rdx2fac)-1.060660171779821*fin[8]*nu[13]*rdx2fac-1.060660171779821*fin[2]*nu[9]*rdx2fac-1.060660171779821*fin[6]*nu[8]*rdx2fac; 
  out[8] += -(1.060660171779821*nu[9]*fin[16]*rdx2fac)-1.060660171779821*fin[6]*nu[13]*rdx2fac-1.060660171779821*fin[2]*nu[11]*rdx2fac-1.060660171779821*fin[8]*nu[8]*rdx2fac; 
  out[10] += -(1.060660171779821*nu[13]*fin[26]*rdx2fac)-1.060660171779821*nu[11]*fin[19]*rdx2fac-1.060660171779821*nu[9]*fin[17]*rdx2fac-1.060660171779821*nu[8]*fin[10]*rdx2fac; 
  out[13] += -(1.060660171779821*nu[13]*fin[27]*rdx2fac)-1.060660171779821*nu[11]*fin[22]*rdx2fac-1.060660171779821*nu[9]*fin[20]*rdx2fac-1.060660171779821*nu[8]*fin[13]*rdx2fac; 
  out[16] += -(1.060660171779821*nu[8]*fin[16]*rdx2fac)-1.060660171779821*fin[2]*nu[13]*rdx2fac-1.060660171779821*fin[6]*nu[11]*rdx2fac-1.060660171779821*fin[8]*nu[9]*rdx2fac; 
  out[17] += -(1.060660171779821*nu[11]*fin[26]*rdx2fac)-1.060660171779821*nu[13]*fin[19]*rdx2fac-1.060660171779821*nu[8]*fin[17]*rdx2fac-1.060660171779821*nu[9]*fin[10]*rdx2fac; 
  out[19] += -(1.060660171779821*nu[9]*fin[26]*rdx2fac)-1.060660171779821*nu[8]*fin[19]*rdx2fac-1.060660171779821*nu[13]*fin[17]*rdx2fac-1.060660171779821*fin[10]*nu[11]*rdx2fac; 
  out[20] += -(1.060660171779821*nu[11]*fin[27]*rdx2fac)-1.060660171779821*nu[13]*fin[22]*rdx2fac-1.060660171779821*nu[8]*fin[20]*rdx2fac-1.060660171779821*nu[9]*fin[13]*rdx2fac; 
  out[22] += -(1.060660171779821*nu[9]*fin[27]*rdx2fac)-1.060660171779821*nu[8]*fin[22]*rdx2fac-1.060660171779821*nu[13]*fin[20]*rdx2fac-1.060660171779821*nu[11]*fin[13]*rdx2fac; 
  out[24] += -(1.060660171779821*nu[13]*fin[31]*rdx2fac)-1.060660171779821*nu[11]*fin[30]*rdx2fac-1.060660171779821*nu[9]*fin[28]*rdx2fac-1.060660171779821*nu[8]*fin[24]*rdx2fac; 
  out[26] += -(1.060660171779821*nu[8]*fin[26]*rdx2fac)-1.060660171779821*nu[9]*fin[19]*rdx2fac-1.060660171779821*nu[11]*fin[17]*rdx2fac-1.060660171779821*fin[10]*nu[13]*rdx2fac; 
  out[27] += -(1.060660171779821*nu[8]*fin[27]*rdx2fac)-1.060660171779821*nu[9]*fin[22]*rdx2fac-1.060660171779821*nu[11]*fin[20]*rdx2fac-1.060660171779821*fin[13]*nu[13]*rdx2fac; 
  out[28] += -(1.060660171779821*nu[11]*fin[31]*rdx2fac)-1.060660171779821*nu[13]*fin[30]*rdx2fac-1.060660171779821*nu[8]*fin[28]*rdx2fac-1.060660171779821*nu[9]*fin[24]*rdx2fac; 
  out[30] += -(1.060660171779821*nu[9]*fin[31]*rdx2fac)-1.060660171779821*nu[8]*fin[30]*rdx2fac-1.060660171779821*nu[13]*fin[28]*rdx2fac-1.060660171779821*nu[11]*fin[24]*rdx2fac; 
  out[31] += -(1.060660171779821*nu[8]*fin[31]*rdx2fac)-1.060660171779821*nu[9]*fin[30]*rdx2fac-1.060660171779821*nu[11]*fin[28]*rdx2fac-1.060660171779821*nu[13]*fin[24]*rdx2fac; 
  out[34] += -(1.060660171779821*nu[13]*fin[43]*rdx2fac)-1.060660171779821*nu[11]*fin[39]*rdx2fac-1.060660171779821*nu[9]*fin[37]*rdx2fac-1.060660171779821*nu[8]*fin[34]*rdx2fac; 
  out[37] += -(1.060660171779821*nu[11]*fin[43]*rdx2fac)-1.060660171779821*nu[13]*fin[39]*rdx2fac-1.060660171779821*nu[8]*fin[37]*rdx2fac-1.060660171779821*nu[9]*fin[34]*rdx2fac; 
  out[39] += -(1.060660171779821*nu[9]*fin[43]*rdx2fac)-1.060660171779821*nu[8]*fin[39]*rdx2fac-1.060660171779821*nu[13]*fin[37]*rdx2fac-1.060660171779821*nu[11]*fin[34]*rdx2fac; 
  out[41] += -(1.060660171779821*nu[13]*fin[47]*rdx2fac)-1.060660171779821*nu[11]*fin[46]*rdx2fac-1.060660171779821*nu[9]*fin[44]*rdx2fac-1.060660171779821*nu[8]*fin[41]*rdx2fac; 
  out[43] += -(1.060660171779821*nu[8]*fin[43]*rdx2fac)-1.060660171779821*nu[9]*fin[39]*rdx2fac-1.060660171779821*nu[11]*fin[37]*rdx2fac-1.060660171779821*nu[13]*fin[34]*rdx2fac; 
  out[44] += -(1.060660171779821*nu[11]*fin[47]*rdx2fac)-1.060660171779821*nu[13]*fin[46]*rdx2fac-1.060660171779821*nu[8]*fin[44]*rdx2fac-1.060660171779821*nu[9]*fin[41]*rdx2fac; 
  out[46] += -(1.060660171779821*nu[9]*fin[47]*rdx2fac)-1.060660171779821*nu[8]*fin[46]*rdx2fac-1.060660171779821*nu[13]*fin[44]*rdx2fac-1.060660171779821*nu[11]*fin[41]*rdx2fac; 
  out[47] += -(1.060660171779821*nu[8]*fin[47]*rdx2fac)-1.060660171779821*nu[9]*fin[46]*rdx2fac-1.060660171779821*nu[11]*fin[44]*rdx2fac-1.060660171779821*nu[13]*fin[41]*rdx2fac; 

  return 1.4142135623730951*nu[8]*rdx2fac; 
}

static GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffz(const double *w, const double *dx, const double *nu, const double *fin, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity.
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // fin: Input field times conf-space Jacobian
  // out: Incremented output

  const double rdx2fac = pow(2.0/dx[2],2); 
  out[3] += -(1.060660171779821*fin[16]*nu[20]*rdx2fac)-1.060660171779821*fin[8]*nu[18]*rdx2fac-1.060660171779821*fin[7]*nu[17]*rdx2fac-1.060660171779821*fin[3]*nu[16]*rdx2fac; 
  out[7] += -(1.060660171779821*fin[8]*nu[20]*rdx2fac)-1.060660171779821*fin[16]*nu[18]*rdx2fac-1.060660171779821*fin[3]*nu[17]*rdx2fac-1.060660171779821*fin[7]*nu[16]*rdx2fac; 
  out[8] += -(1.060660171779821*fin[7]*nu[20]*rdx2fac)-1.060660171779821*fin[3]*nu[18]*rdx2fac-1.060660171779821*fin[16]*nu[17]*rdx2fac-1.060660171779821*fin[8]*nu[16]*rdx2fac; 
  out[11] += -(1.060660171779821*nu[20]*fin[26]*rdx2fac)-1.060660171779821*nu[18]*fin[19]*rdx2fac-1.060660171779821*nu[17]*fin[18]*rdx2fac-1.060660171779821*fin[11]*nu[16]*rdx2fac; 
  out[14] += -(1.060660171779821*nu[20]*fin[27]*rdx2fac)-1.060660171779821*nu[18]*fin[22]*rdx2fac-1.060660171779821*nu[17]*fin[21]*rdx2fac-1.060660171779821*fin[14]*nu[16]*rdx2fac; 
  out[16] += -(1.060660171779821*fin[3]*nu[20]*rdx2fac)-1.060660171779821*fin[7]*nu[18]*rdx2fac-1.060660171779821*fin[8]*nu[17]*rdx2fac-1.060660171779821*fin[16]*nu[16]*rdx2fac; 
  out[18] += -(1.060660171779821*nu[18]*fin[26]*rdx2fac)-1.060660171779821*fin[19]*nu[20]*rdx2fac-1.060660171779821*nu[16]*fin[18]*rdx2fac-1.060660171779821*fin[11]*nu[17]*rdx2fac; 
  out[19] += -(1.060660171779821*nu[17]*fin[26]*rdx2fac)-1.060660171779821*fin[18]*nu[20]*rdx2fac-1.060660171779821*nu[16]*fin[19]*rdx2fac-1.060660171779821*fin[11]*nu[18]*rdx2fac; 
  out[21] += -(1.060660171779821*nu[18]*fin[27]*rdx2fac)-1.060660171779821*nu[20]*fin[22]*rdx2fac-1.060660171779821*nu[16]*fin[21]*rdx2fac-1.060660171779821*fin[14]*nu[17]*rdx2fac; 
  out[22] += -(1.060660171779821*nu[17]*fin[27]*rdx2fac)-1.060660171779821*nu[16]*fin[22]*rdx2fac-1.060660171779821*nu[20]*fin[21]*rdx2fac-1.060660171779821*fin[14]*nu[18]*rdx2fac; 
  out[25] += -(1.060660171779821*nu[20]*fin[31]*rdx2fac)-1.060660171779821*nu[18]*fin[30]*rdx2fac-1.060660171779821*nu[17]*fin[29]*rdx2fac-1.060660171779821*nu[16]*fin[25]*rdx2fac; 
  out[26] += -(1.060660171779821*nu[16]*fin[26]*rdx2fac)-1.060660171779821*fin[11]*nu[20]*rdx2fac-1.060660171779821*nu[17]*fin[19]*rdx2fac-1.060660171779821*fin[18]*nu[18]*rdx2fac; 
  out[27] += -(1.060660171779821*nu[16]*fin[27]*rdx2fac)-1.060660171779821*nu[17]*fin[22]*rdx2fac-1.060660171779821*nu[18]*fin[21]*rdx2fac-1.060660171779821*fin[14]*nu[20]*rdx2fac; 
  out[29] += -(1.060660171779821*nu[18]*fin[31]*rdx2fac)-1.060660171779821*nu[20]*fin[30]*rdx2fac-1.060660171779821*nu[16]*fin[29]*rdx2fac-1.060660171779821*nu[17]*fin[25]*rdx2fac; 
  out[30] += -(1.060660171779821*nu[17]*fin[31]*rdx2fac)-1.060660171779821*nu[16]*fin[30]*rdx2fac-1.060660171779821*nu[20]*fin[29]*rdx2fac-1.060660171779821*nu[18]*fin[25]*rdx2fac; 
  out[31] += -(1.060660171779821*nu[16]*fin[31]*rdx2fac)-1.060660171779821*nu[17]*fin[30]*rdx2fac-1.060660171779821*nu[18]*fin[29]*rdx2fac-1.060660171779821*nu[20]*fin[25]*rdx2fac; 
  out[35] += -(1.060660171779821*nu[20]*fin[43]*rdx2fac)-1.060660171779821*nu[18]*fin[39]*rdx2fac-1.060660171779821*nu[17]*fin[38]*rdx2fac-1.060660171779821*nu[16]*fin[35]*rdx2fac; 
  out[38] += -(1.060660171779821*nu[18]*fin[43]*rdx2fac)-1.060660171779821*nu[20]*fin[39]*rdx2fac-1.060660171779821*nu[16]*fin[38]*rdx2fac-1.060660171779821*nu[17]*fin[35]*rdx2fac; 
  out[39] += -(1.060660171779821*nu[17]*fin[43]*rdx2fac)-1.060660171779821*nu[16]*fin[39]*rdx2fac-1.060660171779821*nu[20]*fin[38]*rdx2fac-1.060660171779821*nu[18]*fin[35]*rdx2fac; 
  out[42] += -(1.060660171779821*nu[20]*fin[47]*rdx2fac)-1.060660171779821*nu[18]*fin[46]*rdx2fac-1.060660171779821*nu[17]*fin[45]*rdx2fac-1.060660171779821*nu[16]*fin[42]*rdx2fac; 
  out[43] += -(1.060660171779821*nu[16]*fin[43]*rdx2fac)-1.060660171779821*nu[17]*fin[39]*rdx2fac-1.060660171779821*nu[18]*fin[38]*rdx2fac-1.060660171779821*nu[20]*fin[35]*rdx2fac; 
  out[45] += -(1.060660171779821*nu[18]*fin[47]*rdx2fac)-1.060660171779821*nu[20]*fin[46]*rdx2fac-1.060660171779821*nu[16]*fin[45]*rdx2fac-1.060660171779821*nu[17]*fin[42]*rdx2fac; 
  out[46] += -(1.060660171779821*nu[17]*fin[47]*rdx2fac)-1.060660171779821*nu[16]*fin[46]*rdx2fac-1.060660171779821*nu[20]*fin[45]*rdx2fac-1.060660171779821*nu[18]*fin[42]*rdx2fac; 
  out[47] += -(1.060660171779821*nu[16]*fin[47]*rdx2fac)-1.060660171779821*nu[17]*fin[46]*rdx2fac-1.060660171779821*nu[18]*fin[45]*rdx2fac-1.060660171779821*nu[20]*fin[42]*rdx2fac; 

  return 1.4142135623730951*nu[16]*rdx2fac; 
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsx(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxy(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffy(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxyz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffy(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffx(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsy(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffy(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsyz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffy(w, dx, nu, fin, out);
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

GKYL_CU_DH double gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsz(const double *w, const double *dx, const double *nu, const double *jacobgeo_inv, const double *Jfin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // nu: Diffusivity
  // jacobgeo_inv: reciprocal of the configuration space Jacobian.
  // Jfin: Input field times conf-space Jacobian
  // out: Incremented output

  double fin[48];
  fin[0] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[16]+jacobgeo_inv[6]*Jfin[8]+jacobgeo_inv[5]*Jfin[7]+jacobgeo_inv[4]*Jfin[6]+Jfin[3]*jacobgeo_inv[3]+Jfin[2]*jacobgeo_inv[2]+Jfin[1]*jacobgeo_inv[1]+Jfin[0]*jacobgeo_inv[0]); 
  fin[1] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[16]+jacobgeo_inv[7]*Jfin[8]+jacobgeo_inv[3]*Jfin[7]+jacobgeo_inv[2]*Jfin[6]+Jfin[3]*jacobgeo_inv[5]+Jfin[2]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*Jfin[1]); 
  fin[2] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[16]+jacobgeo_inv[3]*Jfin[8]+Jfin[7]*jacobgeo_inv[7]+Jfin[3]*jacobgeo_inv[6]+jacobgeo_inv[1]*Jfin[6]+Jfin[1]*jacobgeo_inv[4]+Jfin[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*Jfin[2]); 
  fin[3] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[16]+jacobgeo_inv[2]*Jfin[8]+Jfin[6]*jacobgeo_inv[7]+jacobgeo_inv[1]*Jfin[7]+Jfin[2]*jacobgeo_inv[6]+Jfin[1]*jacobgeo_inv[5]+Jfin[0]*jacobgeo_inv[3]+jacobgeo_inv[0]*Jfin[3]); 
  fin[4] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[26]+jacobgeo_inv[6]*Jfin[19]+jacobgeo_inv[5]*Jfin[18]+jacobgeo_inv[4]*Jfin[17]+jacobgeo_inv[3]*Jfin[11]+jacobgeo_inv[2]*Jfin[10]+jacobgeo_inv[1]*Jfin[9]+jacobgeo_inv[0]*Jfin[4]); 
  fin[5] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[27]+jacobgeo_inv[6]*Jfin[22]+jacobgeo_inv[5]*Jfin[21]+jacobgeo_inv[4]*Jfin[20]+jacobgeo_inv[3]*Jfin[14]+jacobgeo_inv[2]*Jfin[13]+jacobgeo_inv[1]*Jfin[12]+jacobgeo_inv[0]*Jfin[5]); 
  fin[6] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[16]+jacobgeo_inv[5]*Jfin[8]+Jfin[3]*jacobgeo_inv[7]+jacobgeo_inv[6]*Jfin[7]+jacobgeo_inv[0]*Jfin[6]+Jfin[0]*jacobgeo_inv[4]+Jfin[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*Jfin[2]); 
  fin[7] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[16]+jacobgeo_inv[4]*Jfin[8]+Jfin[2]*jacobgeo_inv[7]+jacobgeo_inv[0]*Jfin[7]+Jfin[6]*jacobgeo_inv[6]+Jfin[0]*jacobgeo_inv[5]+Jfin[1]*jacobgeo_inv[3]+jacobgeo_inv[1]*Jfin[3]); 
  fin[8] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[16]+jacobgeo_inv[0]*Jfin[8]+Jfin[1]*jacobgeo_inv[7]+jacobgeo_inv[4]*Jfin[7]+Jfin[0]*jacobgeo_inv[6]+jacobgeo_inv[5]*Jfin[6]+Jfin[2]*jacobgeo_inv[3]+jacobgeo_inv[2]*Jfin[3]); 
  fin[9] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[26]+jacobgeo_inv[7]*Jfin[19]+jacobgeo_inv[3]*Jfin[18]+jacobgeo_inv[2]*Jfin[17]+jacobgeo_inv[5]*Jfin[11]+jacobgeo_inv[4]*Jfin[10]+jacobgeo_inv[0]*Jfin[9]+jacobgeo_inv[1]*Jfin[4]); 
  fin[10] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[26]+jacobgeo_inv[3]*Jfin[19]+jacobgeo_inv[7]*Jfin[18]+jacobgeo_inv[1]*Jfin[17]+jacobgeo_inv[6]*Jfin[11]+jacobgeo_inv[0]*Jfin[10]+jacobgeo_inv[4]*Jfin[9]+jacobgeo_inv[2]*Jfin[4]); 
  fin[11] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[26]+jacobgeo_inv[2]*Jfin[19]+jacobgeo_inv[1]*Jfin[18]+jacobgeo_inv[7]*Jfin[17]+jacobgeo_inv[0]*Jfin[11]+jacobgeo_inv[6]*Jfin[10]+jacobgeo_inv[5]*Jfin[9]+jacobgeo_inv[3]*Jfin[4]); 
  fin[12] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[27]+jacobgeo_inv[7]*Jfin[22]+jacobgeo_inv[3]*Jfin[21]+jacobgeo_inv[2]*Jfin[20]+jacobgeo_inv[5]*Jfin[14]+jacobgeo_inv[4]*Jfin[13]+jacobgeo_inv[0]*Jfin[12]+jacobgeo_inv[1]*Jfin[5]); 
  fin[13] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[27]+jacobgeo_inv[3]*Jfin[22]+jacobgeo_inv[7]*Jfin[21]+jacobgeo_inv[1]*Jfin[20]+jacobgeo_inv[6]*Jfin[14]+jacobgeo_inv[0]*Jfin[13]+jacobgeo_inv[4]*Jfin[12]+jacobgeo_inv[2]*Jfin[5]); 
  fin[14] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[27]+jacobgeo_inv[2]*Jfin[22]+jacobgeo_inv[1]*Jfin[21]+jacobgeo_inv[7]*Jfin[20]+jacobgeo_inv[0]*Jfin[14]+jacobgeo_inv[6]*Jfin[13]+jacobgeo_inv[5]*Jfin[12]+jacobgeo_inv[3]*Jfin[5]); 
  fin[15] = 0.3535533905932737*(jacobgeo_inv[7]*Jfin[31]+jacobgeo_inv[6]*Jfin[30]+jacobgeo_inv[5]*Jfin[29]+jacobgeo_inv[4]*Jfin[28]+jacobgeo_inv[3]*Jfin[25]+jacobgeo_inv[2]*Jfin[24]+jacobgeo_inv[1]*Jfin[23]+jacobgeo_inv[0]*Jfin[15]); 
  fin[16] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[16]+jacobgeo_inv[1]*Jfin[8]+Jfin[0]*jacobgeo_inv[7]+jacobgeo_inv[2]*Jfin[7]+Jfin[1]*jacobgeo_inv[6]+jacobgeo_inv[3]*Jfin[6]+Jfin[2]*jacobgeo_inv[5]+Jfin[3]*jacobgeo_inv[4]); 
  fin[17] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[26]+jacobgeo_inv[5]*Jfin[19]+jacobgeo_inv[6]*Jfin[18]+jacobgeo_inv[0]*Jfin[17]+jacobgeo_inv[7]*Jfin[11]+jacobgeo_inv[1]*Jfin[10]+jacobgeo_inv[2]*Jfin[9]+Jfin[4]*jacobgeo_inv[4]); 
  fin[18] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[26]+jacobgeo_inv[4]*Jfin[19]+jacobgeo_inv[0]*Jfin[18]+jacobgeo_inv[6]*Jfin[17]+jacobgeo_inv[1]*Jfin[11]+jacobgeo_inv[7]*Jfin[10]+jacobgeo_inv[3]*Jfin[9]+Jfin[4]*jacobgeo_inv[5]); 
  fin[19] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[26]+jacobgeo_inv[0]*Jfin[19]+jacobgeo_inv[4]*Jfin[18]+jacobgeo_inv[5]*Jfin[17]+jacobgeo_inv[2]*Jfin[11]+jacobgeo_inv[3]*Jfin[10]+jacobgeo_inv[7]*Jfin[9]+Jfin[4]*jacobgeo_inv[6]); 
  fin[20] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[27]+jacobgeo_inv[5]*Jfin[22]+jacobgeo_inv[6]*Jfin[21]+jacobgeo_inv[0]*Jfin[20]+jacobgeo_inv[7]*Jfin[14]+jacobgeo_inv[1]*Jfin[13]+jacobgeo_inv[2]*Jfin[12]+jacobgeo_inv[4]*Jfin[5]); 
  fin[21] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[27]+jacobgeo_inv[4]*Jfin[22]+jacobgeo_inv[0]*Jfin[21]+jacobgeo_inv[6]*Jfin[20]+jacobgeo_inv[1]*Jfin[14]+jacobgeo_inv[7]*Jfin[13]+jacobgeo_inv[3]*Jfin[12]+Jfin[5]*jacobgeo_inv[5]); 
  fin[22] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[27]+jacobgeo_inv[0]*Jfin[22]+jacobgeo_inv[4]*Jfin[21]+jacobgeo_inv[5]*Jfin[20]+jacobgeo_inv[2]*Jfin[14]+jacobgeo_inv[3]*Jfin[13]+jacobgeo_inv[7]*Jfin[12]+Jfin[5]*jacobgeo_inv[6]); 
  fin[23] = 0.3535533905932737*(jacobgeo_inv[6]*Jfin[31]+jacobgeo_inv[7]*Jfin[30]+jacobgeo_inv[3]*Jfin[29]+jacobgeo_inv[2]*Jfin[28]+jacobgeo_inv[5]*Jfin[25]+jacobgeo_inv[4]*Jfin[24]+jacobgeo_inv[0]*Jfin[23]+jacobgeo_inv[1]*Jfin[15]); 
  fin[24] = 0.3535533905932737*(jacobgeo_inv[5]*Jfin[31]+jacobgeo_inv[3]*Jfin[30]+jacobgeo_inv[7]*Jfin[29]+jacobgeo_inv[1]*Jfin[28]+jacobgeo_inv[6]*Jfin[25]+jacobgeo_inv[0]*Jfin[24]+jacobgeo_inv[4]*Jfin[23]+jacobgeo_inv[2]*Jfin[15]); 
  fin[25] = 0.3535533905932737*(jacobgeo_inv[4]*Jfin[31]+jacobgeo_inv[2]*Jfin[30]+jacobgeo_inv[1]*Jfin[29]+jacobgeo_inv[7]*Jfin[28]+jacobgeo_inv[0]*Jfin[25]+jacobgeo_inv[6]*Jfin[24]+jacobgeo_inv[5]*Jfin[23]+jacobgeo_inv[3]*Jfin[15]); 
  fin[26] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[26]+jacobgeo_inv[1]*Jfin[19]+jacobgeo_inv[2]*Jfin[18]+jacobgeo_inv[3]*Jfin[17]+jacobgeo_inv[4]*Jfin[11]+jacobgeo_inv[5]*Jfin[10]+jacobgeo_inv[6]*Jfin[9]+Jfin[4]*jacobgeo_inv[7]); 
  fin[27] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[27]+jacobgeo_inv[1]*Jfin[22]+jacobgeo_inv[2]*Jfin[21]+jacobgeo_inv[3]*Jfin[20]+jacobgeo_inv[4]*Jfin[14]+jacobgeo_inv[5]*Jfin[13]+jacobgeo_inv[6]*Jfin[12]+Jfin[5]*jacobgeo_inv[7]); 
  fin[28] = 0.3535533905932737*(jacobgeo_inv[3]*Jfin[31]+jacobgeo_inv[5]*Jfin[30]+jacobgeo_inv[6]*Jfin[29]+jacobgeo_inv[0]*Jfin[28]+jacobgeo_inv[7]*Jfin[25]+jacobgeo_inv[1]*Jfin[24]+jacobgeo_inv[2]*Jfin[23]+jacobgeo_inv[4]*Jfin[15]); 
  fin[29] = 0.3535533905932737*(jacobgeo_inv[2]*Jfin[31]+jacobgeo_inv[4]*Jfin[30]+jacobgeo_inv[0]*Jfin[29]+jacobgeo_inv[6]*Jfin[28]+jacobgeo_inv[1]*Jfin[25]+jacobgeo_inv[7]*Jfin[24]+jacobgeo_inv[3]*Jfin[23]+jacobgeo_inv[5]*Jfin[15]); 
  fin[30] = 0.3535533905932737*(jacobgeo_inv[1]*Jfin[31]+jacobgeo_inv[0]*Jfin[30]+jacobgeo_inv[4]*Jfin[29]+jacobgeo_inv[5]*Jfin[28]+jacobgeo_inv[2]*Jfin[25]+jacobgeo_inv[3]*Jfin[24]+jacobgeo_inv[7]*Jfin[23]+jacobgeo_inv[6]*Jfin[15]); 
  fin[31] = 0.3535533905932737*(jacobgeo_inv[0]*Jfin[31]+jacobgeo_inv[1]*Jfin[30]+jacobgeo_inv[2]*Jfin[29]+jacobgeo_inv[3]*Jfin[28]+jacobgeo_inv[4]*Jfin[25]+jacobgeo_inv[5]*Jfin[24]+jacobgeo_inv[6]*Jfin[23]+jacobgeo_inv[7]*Jfin[15]); 
  fin[32] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[43]+15.0*(jacobgeo_inv[6]*Jfin[39]+jacobgeo_inv[5]*Jfin[38]+jacobgeo_inv[4]*Jfin[37])+15.000000000000002*(jacobgeo_inv[3]*Jfin[35]+jacobgeo_inv[2]*Jfin[34]+jacobgeo_inv[1]*Jfin[33])+15.0*jacobgeo_inv[0]*Jfin[32]); 
  fin[33] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[43]+15.000000000000002*(jacobgeo_inv[7]*Jfin[39]+jacobgeo_inv[3]*Jfin[38]+jacobgeo_inv[2]*Jfin[37])+15.0*(jacobgeo_inv[5]*Jfin[35]+jacobgeo_inv[4]*Jfin[34]+jacobgeo_inv[0]*Jfin[33])+15.000000000000002*jacobgeo_inv[1]*Jfin[32]); 
  fin[34] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[43]+15.000000000000002*(jacobgeo_inv[3]*Jfin[39]+jacobgeo_inv[7]*Jfin[38]+jacobgeo_inv[1]*Jfin[37])+15.0*(jacobgeo_inv[6]*Jfin[35]+jacobgeo_inv[0]*Jfin[34]+jacobgeo_inv[4]*Jfin[33])+15.000000000000002*jacobgeo_inv[2]*Jfin[32]); 
  fin[35] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[43]+15.000000000000002*(jacobgeo_inv[2]*Jfin[39]+jacobgeo_inv[1]*Jfin[38]+jacobgeo_inv[7]*Jfin[37])+15.0*(jacobgeo_inv[0]*Jfin[35]+jacobgeo_inv[6]*Jfin[34]+jacobgeo_inv[5]*Jfin[33])+15.000000000000002*jacobgeo_inv[3]*Jfin[32]); 
  fin[36] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[7]*Jfin[47]+15.0*(jacobgeo_inv[6]*Jfin[46]+jacobgeo_inv[5]*Jfin[45]+jacobgeo_inv[4]*Jfin[44])+15.000000000000002*(jacobgeo_inv[3]*Jfin[42]+jacobgeo_inv[2]*Jfin[41]+jacobgeo_inv[1]*Jfin[40])+15.0*jacobgeo_inv[0]*Jfin[36]); 
  fin[37] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[43]+15.0*(jacobgeo_inv[5]*Jfin[39]+jacobgeo_inv[6]*Jfin[38]+jacobgeo_inv[0]*Jfin[37])+15.000000000000002*(jacobgeo_inv[7]*Jfin[35]+jacobgeo_inv[1]*Jfin[34]+jacobgeo_inv[2]*Jfin[33])+15.0*jacobgeo_inv[4]*Jfin[32]); 
  fin[38] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[43]+15.0*(jacobgeo_inv[4]*Jfin[39]+jacobgeo_inv[0]*Jfin[38]+jacobgeo_inv[6]*Jfin[37])+15.000000000000002*(jacobgeo_inv[1]*Jfin[35]+jacobgeo_inv[7]*Jfin[34]+jacobgeo_inv[3]*Jfin[33])+15.0*jacobgeo_inv[5]*Jfin[32]); 
  fin[39] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[43]+15.0*(jacobgeo_inv[0]*Jfin[39]+jacobgeo_inv[4]*Jfin[38]+jacobgeo_inv[5]*Jfin[37])+15.000000000000002*(jacobgeo_inv[2]*Jfin[35]+jacobgeo_inv[3]*Jfin[34]+jacobgeo_inv[7]*Jfin[33])+15.0*jacobgeo_inv[6]*Jfin[32]); 
  fin[40] = 0.023570226039551577*(15.0*jacobgeo_inv[6]*Jfin[47]+15.000000000000002*(jacobgeo_inv[7]*Jfin[46]+jacobgeo_inv[3]*Jfin[45]+jacobgeo_inv[2]*Jfin[44])+15.0*(jacobgeo_inv[5]*Jfin[42]+jacobgeo_inv[4]*Jfin[41]+jacobgeo_inv[0]*Jfin[40])+15.000000000000002*jacobgeo_inv[1]*Jfin[36]); 
  fin[41] = 0.023570226039551577*(15.0*jacobgeo_inv[5]*Jfin[47]+15.000000000000002*(jacobgeo_inv[3]*Jfin[46]+jacobgeo_inv[7]*Jfin[45]+jacobgeo_inv[1]*Jfin[44])+15.0*(jacobgeo_inv[6]*Jfin[42]+jacobgeo_inv[0]*Jfin[41]+jacobgeo_inv[4]*Jfin[40])+15.000000000000002*jacobgeo_inv[2]*Jfin[36]); 
  fin[42] = 0.023570226039551577*(15.0*jacobgeo_inv[4]*Jfin[47]+15.000000000000002*(jacobgeo_inv[2]*Jfin[46]+jacobgeo_inv[1]*Jfin[45]+jacobgeo_inv[7]*Jfin[44])+15.0*(jacobgeo_inv[0]*Jfin[42]+jacobgeo_inv[6]*Jfin[41]+jacobgeo_inv[5]*Jfin[40])+15.000000000000002*jacobgeo_inv[3]*Jfin[36]); 
  fin[43] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[43]+15.000000000000002*(jacobgeo_inv[1]*Jfin[39]+jacobgeo_inv[2]*Jfin[38]+jacobgeo_inv[3]*Jfin[37])+15.0*(jacobgeo_inv[4]*Jfin[35]+jacobgeo_inv[5]*Jfin[34]+jacobgeo_inv[6]*Jfin[33])+15.000000000000002*jacobgeo_inv[7]*Jfin[32]); 
  fin[44] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[3]*Jfin[47]+15.0*(jacobgeo_inv[5]*Jfin[46]+jacobgeo_inv[6]*Jfin[45]+jacobgeo_inv[0]*Jfin[44])+15.000000000000002*(jacobgeo_inv[7]*Jfin[42]+jacobgeo_inv[1]*Jfin[41]+jacobgeo_inv[2]*Jfin[40])+15.0*jacobgeo_inv[4]*Jfin[36]); 
  fin[45] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[2]*Jfin[47]+15.0*(jacobgeo_inv[4]*Jfin[46]+jacobgeo_inv[0]*Jfin[45]+jacobgeo_inv[6]*Jfin[44])+15.000000000000002*(jacobgeo_inv[1]*Jfin[42]+jacobgeo_inv[7]*Jfin[41]+jacobgeo_inv[3]*Jfin[40])+15.0*jacobgeo_inv[5]*Jfin[36]); 
  fin[46] = 0.023570226039551577*(15.000000000000002*jacobgeo_inv[1]*Jfin[47]+15.0*(jacobgeo_inv[0]*Jfin[46]+jacobgeo_inv[4]*Jfin[45]+jacobgeo_inv[5]*Jfin[44])+15.000000000000002*(jacobgeo_inv[2]*Jfin[42]+jacobgeo_inv[3]*Jfin[41]+jacobgeo_inv[7]*Jfin[40])+15.0*jacobgeo_inv[6]*Jfin[36]); 
  fin[47] = 0.023570226039551577*(15.0*jacobgeo_inv[0]*Jfin[47]+15.000000000000002*(jacobgeo_inv[1]*Jfin[46]+jacobgeo_inv[2]*Jfin[45]+jacobgeo_inv[3]*Jfin[44])+15.0*(jacobgeo_inv[4]*Jfin[42]+jacobgeo_inv[5]*Jfin[41]+jacobgeo_inv[6]*Jfin[40])+15.000000000000002*jacobgeo_inv[7]*Jfin[36]); 

  double cflFreq = 0.;
  
  cflFreq += gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffz(w, dx, nu, fin, out);

  return cflFreq;
}

