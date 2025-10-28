#include <gkyl_mom_bcorr_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_2x2v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *vmap_prime, double _m, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direction. 
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // _m:        species mass. 
  // fIn:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  dS = (3.141592653589793*dxv[3])/_m; 
 
  if (edge == GKYL_VX_LOWER) {
 
    out[0] += (((-2.23606797749979*fIn[16])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS)/vmap_prime[0]; 
    out[1] += (((-2.23606797749979*fIn[17])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS)/vmap_prime[0]; 
    out[2] += (((-2.23606797749979*fIn[18])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS)/vmap_prime[0]; 
    out[3] += (((-2.23606797749979*fIn[20])+1.732050807568877*fIn[11]-1.0*fIn[5])*dS)/vmap_prime[0]; 
    out[4] += (vBoundary[0]*((-2.23606797749979*fIn[16])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS)/vmap_prime[0]; 
    out[5] += (vBoundary[0]*((-2.23606797749979*fIn[17])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS)/vmap_prime[0]; 
    out[6] += (vBoundary[0]*((-2.23606797749979*fIn[18])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS)/vmap_prime[0]; 
    out[7] += (vBoundary[0]*((-2.23606797749979*fIn[20])+1.732050807568877*fIn[11]-1.0*fIn[5])*dS)/vmap_prime[0]; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    out[0] += ((2.23606797749979*fIn[16]+1.732050807568877*fIn[3]+fIn[0])*dS)/vmap_prime[0]; 
    out[1] += ((2.23606797749979*fIn[17]+1.732050807568877*fIn[6]+fIn[1])*dS)/vmap_prime[0]; 
    out[2] += ((2.23606797749979*fIn[18]+1.732050807568877*fIn[7]+fIn[2])*dS)/vmap_prime[0]; 
    out[3] += ((2.23606797749979*fIn[20]+1.732050807568877*fIn[11]+fIn[5])*dS)/vmap_prime[0]; 
    out[4] += (vBoundary[2]*(2.23606797749979*fIn[16]+1.732050807568877*fIn[3]+fIn[0])*dS)/vmap_prime[0]; 
    out[5] += (vBoundary[2]*(2.23606797749979*fIn[17]+1.732050807568877*fIn[6]+fIn[1])*dS)/vmap_prime[0]; 
    out[6] += (vBoundary[2]*(2.23606797749979*fIn[18]+1.732050807568877*fIn[7]+fIn[2])*dS)/vmap_prime[0]; 
    out[7] += (vBoundary[2]*(2.23606797749979*fIn[20]+1.732050807568877*fIn[11]+fIn[5])*dS)/vmap_prime[0]; 
 
  }
 
  dS = (6.283185307179586*dxv[2])/_m; 
 
  if (edge == GKYL_VY_LOWER) {
 
    out[4] += (vBoundary[1]*(1.732050807568877*fIn[4]-1.0*fIn[0])*dS)/vmap_prime[1]; 
    out[5] += (vBoundary[1]*(1.732050807568877*fIn[8]-1.0*fIn[1])*dS)/vmap_prime[1]; 
    out[6] += (vBoundary[1]*(1.732050807568877*fIn[9]-1.0*fIn[2])*dS)/vmap_prime[1]; 
    out[7] += (vBoundary[1]*(1.732050807568877*fIn[12]-1.0*fIn[5])*dS)/vmap_prime[1]; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    out[4] += (vBoundary[3]*(1.732050807568877*fIn[4]+fIn[0])*dS)/vmap_prime[1]; 
    out[5] += (vBoundary[3]*(1.732050807568877*fIn[8]+fIn[1])*dS)/vmap_prime[1]; 
    out[6] += (vBoundary[3]*(1.732050807568877*fIn[9]+fIn[2])*dS)/vmap_prime[1]; 
    out[7] += (vBoundary[3]*(1.732050807568877*fIn[12]+fIn[5])*dS)/vmap_prime[1]; 
 
  }
 
} 
 
