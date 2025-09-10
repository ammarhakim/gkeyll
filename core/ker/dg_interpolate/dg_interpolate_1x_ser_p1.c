#include <gkyl_dg_interpolate_kernels.h> 
 
GKYL_CU_DH void dg_interpolate_1x_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
{ 
  // wDo: cell center of donor cell.
  // wTar: cell center of target cell.
  // dxDo: cell length of donor cell.
  // dxTar: cell length of target cell.
  // fldDo: donor field.
  // fldTar: target field in cells pointed to by the stencil.

  double eLo = fmax(-1.0,1.0-(2.0*(wTar[0]-1.0*wDo[0]+0.5*dxTar[0]+0.5*dxDo[0]))/dxTar[0]);
  double eUp = fmin( 1.0,(2.0*(-(1.0*wTar[0])+wDo[0]+0.5*dxTar[0]+0.5*dxDo[0]))/dxTar[0]-1.0);

#ifdef __CUDA_ARCH__
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  atomicAdd(&fldTar[0], (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2); 
#else
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);

  fldTar[0] += (dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]+0.5*fldDo[0]*eUp-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]+0.4330127018922193*fldDo[0]*eUpR2-0.4330127018922193*fldDo[0]*eLoR2; 
#endif

}

