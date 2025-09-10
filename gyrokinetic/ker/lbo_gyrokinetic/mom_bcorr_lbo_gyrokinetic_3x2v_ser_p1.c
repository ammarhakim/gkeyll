#include <gkyl_mom_bcorr_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_3x2v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direction. 
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // _m:        species mass. 
  // fIn:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  dS = (3.141592653589793*dxv[4])/_m; 
 
  if (edge == GKYL_VX_LOWER) {
 
    out[0] += (((-2.23606797749979*fIn[32])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS)/vmap_prime[0]; 
    out[1] += (((-2.23606797749979*fIn[33])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS)/vmap_prime[0]; 
    out[2] += (((-2.23606797749979*fIn[34])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS)/vmap_prime[0]; 
    out[3] += (((-2.23606797749979*fIn[35])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS)/vmap_prime[0]; 
    out[4] += (((-2.23606797749979*fIn[37])+1.732050807568877*fIn[17]-1.0*fIn[6])*dS)/vmap_prime[0]; 
    out[5] += (((-2.23606797749979*fIn[38])+1.732050807568877*fIn[18]-1.0*fIn[7])*dS)/vmap_prime[0]; 
    out[6] += (((-2.23606797749979*fIn[39])+1.732050807568877*fIn[19]-1.0*fIn[8])*dS)/vmap_prime[0]; 
    out[7] += (((-2.23606797749979*fIn[43])+1.732050807568877*fIn[26]-1.0*fIn[16])*dS)/vmap_prime[0]; 
    out[8] += (vBoundary[0]*((-2.23606797749979*fIn[32])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS)/vmap_prime[0]; 
    out[9] += (vBoundary[0]*((-2.23606797749979*fIn[33])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS)/vmap_prime[0]; 
    out[10] += (vBoundary[0]*((-2.23606797749979*fIn[34])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS)/vmap_prime[0]; 
    out[11] += (vBoundary[0]*((-2.23606797749979*fIn[35])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS)/vmap_prime[0]; 
    out[12] += (vBoundary[0]*((-2.23606797749979*fIn[37])+1.732050807568877*fIn[17]-1.0*fIn[6])*dS)/vmap_prime[0]; 
    out[13] += (vBoundary[0]*((-2.23606797749979*fIn[38])+1.732050807568877*fIn[18]-1.0*fIn[7])*dS)/vmap_prime[0]; 
    out[14] += (vBoundary[0]*((-2.23606797749979*fIn[39])+1.732050807568877*fIn[19]-1.0*fIn[8])*dS)/vmap_prime[0]; 
    out[15] += (vBoundary[0]*((-2.23606797749979*fIn[43])+1.732050807568877*fIn[26]-1.0*fIn[16])*dS)/vmap_prime[0]; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    out[0] += ((2.23606797749979*fIn[32]+1.732050807568877*fIn[4]+fIn[0])*dS)/vmap_prime[0]; 
    out[1] += ((2.23606797749979*fIn[33]+1.732050807568877*fIn[9]+fIn[1])*dS)/vmap_prime[0]; 
    out[2] += ((2.23606797749979*fIn[34]+1.732050807568877*fIn[10]+fIn[2])*dS)/vmap_prime[0]; 
    out[3] += ((2.23606797749979*fIn[35]+1.732050807568877*fIn[11]+fIn[3])*dS)/vmap_prime[0]; 
    out[4] += ((2.23606797749979*fIn[37]+1.732050807568877*fIn[17]+fIn[6])*dS)/vmap_prime[0]; 
    out[5] += ((2.23606797749979*fIn[38]+1.732050807568877*fIn[18]+fIn[7])*dS)/vmap_prime[0]; 
    out[6] += ((2.23606797749979*fIn[39]+1.732050807568877*fIn[19]+fIn[8])*dS)/vmap_prime[0]; 
    out[7] += ((2.23606797749979*fIn[43]+1.732050807568877*fIn[26]+fIn[16])*dS)/vmap_prime[0]; 
    out[8] += (vBoundary[2]*(2.23606797749979*fIn[32]+1.732050807568877*fIn[4]+fIn[0])*dS)/vmap_prime[0]; 
    out[9] += (vBoundary[2]*(2.23606797749979*fIn[33]+1.732050807568877*fIn[9]+fIn[1])*dS)/vmap_prime[0]; 
    out[10] += (vBoundary[2]*(2.23606797749979*fIn[34]+1.732050807568877*fIn[10]+fIn[2])*dS)/vmap_prime[0]; 
    out[11] += (vBoundary[2]*(2.23606797749979*fIn[35]+1.732050807568877*fIn[11]+fIn[3])*dS)/vmap_prime[0]; 
    out[12] += (vBoundary[2]*(2.23606797749979*fIn[37]+1.732050807568877*fIn[17]+fIn[6])*dS)/vmap_prime[0]; 
    out[13] += (vBoundary[2]*(2.23606797749979*fIn[38]+1.732050807568877*fIn[18]+fIn[7])*dS)/vmap_prime[0]; 
    out[14] += (vBoundary[2]*(2.23606797749979*fIn[39]+1.732050807568877*fIn[19]+fIn[8])*dS)/vmap_prime[0]; 
    out[15] += (vBoundary[2]*(2.23606797749979*fIn[43]+1.732050807568877*fIn[26]+fIn[16])*dS)/vmap_prime[0]; 
 
  }
 
  dS = (6.283185307179586*dxv[3])/_m; 
 
  if (edge == GKYL_VY_LOWER) {
 
    out[8] += (vBoundary[1]*(1.732050807568877*fIn[5]-1.0*fIn[0])*dS)/vmap_prime[1]; 
    out[9] += (vBoundary[1]*(1.732050807568877*fIn[12]-1.0*fIn[1])*dS)/vmap_prime[1]; 
    out[10] += (vBoundary[1]*(1.732050807568877*fIn[13]-1.0*fIn[2])*dS)/vmap_prime[1]; 
    out[11] += (vBoundary[1]*(1.732050807568877*fIn[14]-1.0*fIn[3])*dS)/vmap_prime[1]; 
    out[12] += (vBoundary[1]*(1.732050807568877*fIn[20]-1.0*fIn[6])*dS)/vmap_prime[1]; 
    out[13] += (vBoundary[1]*(1.732050807568877*fIn[21]-1.0*fIn[7])*dS)/vmap_prime[1]; 
    out[14] += (vBoundary[1]*(1.732050807568877*fIn[22]-1.0*fIn[8])*dS)/vmap_prime[1]; 
    out[15] += (vBoundary[1]*(1.732050807568877*fIn[27]-1.0*fIn[16])*dS)/vmap_prime[1]; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    out[8] += (vBoundary[3]*(1.732050807568877*fIn[5]+fIn[0])*dS)/vmap_prime[1]; 
    out[9] += (vBoundary[3]*(1.732050807568877*fIn[12]+fIn[1])*dS)/vmap_prime[1]; 
    out[10] += (vBoundary[3]*(1.732050807568877*fIn[13]+fIn[2])*dS)/vmap_prime[1]; 
    out[11] += (vBoundary[3]*(1.732050807568877*fIn[14]+fIn[3])*dS)/vmap_prime[1]; 
    out[12] += (vBoundary[3]*(1.732050807568877*fIn[20]+fIn[6])*dS)/vmap_prime[1]; 
    out[13] += (vBoundary[3]*(1.732050807568877*fIn[21]+fIn[7])*dS)/vmap_prime[1]; 
    out[14] += (vBoundary[3]*(1.732050807568877*fIn[22]+fIn[8])*dS)/vmap_prime[1]; 
    out[15] += (vBoundary[3]*(1.732050807568877*fIn[27]+fIn[16])*dS)/vmap_prime[1]; 
 
  }
 
} 
 
