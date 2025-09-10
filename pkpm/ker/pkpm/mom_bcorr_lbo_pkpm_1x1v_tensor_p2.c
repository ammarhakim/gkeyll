#include <gkyl_mom_bcorr_lbo_pkpm_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_pkpm_1x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[9]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f^(vmax)_(vmin) and vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*mass; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*mass; 
    out[2] += ((-1.58113883008419*fIn[8])+1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*mass; 
    out[3] += out[0]*vBoundary[0]; 
    out[4] += vBoundary[0]*out[1]; 
    out[5] += vBoundary[0]*out[2]; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*mass; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*mass; 
    out[2] += (1.58113883008419*fIn[8]+1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*mass; 
    out[3] += out[0]*vBoundary[1]; 
    out[4] += out[1]*vBoundary[1]; 
    out[5] += vBoundary[1]*out[2]; 
  }
 
} 
 
