#include <gkyl_dg_interpolate_kernels.h> 
 
GKYL_CU_DH void dg_interpolate_1x_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar) 
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
  const double wDo0R2 = pow(wDo[0],2);
  const double wTar0R2 = pow(wTar[0],2);
  const double dxDo0R2 = pow(dxDo[0],2);
  const double dxTar0R2 = pow(dxTar[0],2);
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eLoR4 = pow(eLo,4);
  const double eLoR5 = pow(eLo,5);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);
  const double eUpR4 = pow(eUp,4);
  const double eUpR5 = pow(eUp,5);

  atomicAdd(&fldTar[0], (dxTar0R2*(0.5590169943749475*fldDo[2]*eUpR3-0.5590169943749475*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(3.3541019662496847*wTar[0]*fldDo[2]*eUpR2-3.3541019662496847*wDo[0]*fldDo[2]*eUpR2-3.3541019662496847*wTar[0]*fldDo[2]*eLoR2+3.3541019662496847*wDo[0]*fldDo[2]*eLoR2))/dxDo0R2+(dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(6.708203932499369*wTar0R2*fldDo[2]*eUp-13.416407864998739*wDo[0]*wTar[0]*fldDo[2]*eUp+6.708203932499369*wDo0R2*fldDo[2]*eUp-6.708203932499369*wTar0R2*fldDo[2]*eLo+13.416407864998739*wDo[0]*wTar[0]*fldDo[2]*eLo-6.708203932499369*wDo0R2*fldDo[2]*eLo)/dxDo0R2+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]-0.5590169943749475*fldDo[2]*eUp+0.5*fldDo[0]*eUp+0.5590169943749475*fldDo[2]*eLo-0.5*fldDo[0]*eLo); 
  atomicAdd(&fldTar[1], (dxTar0R2*(0.7261843774138906*fldDo[2]*eUpR4-0.7261843774138906*fldDo[2]*eLoR4))/dxDo0R2+(dxTar[0]*(3.872983346207417*wTar[0]*fldDo[2]*eUpR3-3.872983346207417*wDo[0]*fldDo[2]*eUpR3-3.872983346207417*wTar[0]*fldDo[2]*eLoR3+3.872983346207417*wDo[0]*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(5.809475019311125*wTar0R2*fldDo[2]*eUpR2-11.61895003862225*wDo[0]*wTar[0]*fldDo[2]*eUpR2+5.809475019311125*wDo0R2*fldDo[2]*eUpR2-5.809475019311125*wTar0R2*fldDo[2]*eLoR2+11.61895003862225*wDo[0]*wTar[0]*fldDo[2]*eLoR2-5.809475019311125*wDo0R2*fldDo[2]*eLoR2)/dxDo0R2+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]-0.4841229182759271*fldDo[2]*eUpR2+0.4330127018922193*fldDo[0]*eUpR2+0.4841229182759271*fldDo[2]*eLoR2-0.4330127018922193*fldDo[0]*eLoR2); 
  atomicAdd(&fldTar[2], (dxTar0R2*(1.125*fldDo[2]*eUpR5-0.625*fldDo[2]*eUpR3-1.125*fldDo[2]*eLoR5+0.625*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(5.625*wTar[0]*fldDo[2]*eUpR4-5.625*wDo[0]*fldDo[2]*eUpR4-3.75*wTar[0]*fldDo[2]*eUpR2+3.75*wDo[0]*fldDo[2]*eUpR2-5.625*wTar[0]*fldDo[2]*eLoR4+5.625*wDo[0]*fldDo[2]*eLoR4+3.75*wTar[0]*fldDo[2]*eLoR2-3.75*wDo[0]*fldDo[2]*eLoR2))/dxDo0R2+(dxTar[0]*(0.7261843774138906*fldDo[1]*eUpR4-0.4841229182759271*fldDo[1]*eUpR2-0.7261843774138906*fldDo[1]*eLoR4+0.4841229182759271*fldDo[1]*eLoR2))/dxDo[0]+(7.5*wTar0R2*fldDo[2]*eUpR3-15.0*wDo[0]*wTar[0]*fldDo[2]*eUpR3+7.5*wDo0R2*fldDo[2]*eUpR3-7.5*wTar0R2*fldDo[2]*eUp+15.0*wDo[0]*wTar[0]*fldDo[2]*eUp-7.5*wDo0R2*fldDo[2]*eUp-7.5*wTar0R2*fldDo[2]*eLoR3+15.0*wDo[0]*wTar[0]*fldDo[2]*eLoR3-7.5*wDo0R2*fldDo[2]*eLoR3+7.5*wTar0R2*fldDo[2]*eLo-15.0*wDo[0]*wTar[0]*fldDo[2]*eLo+7.5*wDo0R2*fldDo[2]*eLo)/dxDo0R2+(1.9364916731037085*wTar[0]*fldDo[1]*eUpR3-1.9364916731037085*wDo[0]*fldDo[1]*eUpR3-1.9364916731037085*wTar[0]*fldDo[1]*eUp+1.9364916731037085*wDo[0]*fldDo[1]*eUp-1.9364916731037085*wTar[0]*fldDo[1]*eLoR3+1.9364916731037085*wDo[0]*fldDo[1]*eLoR3+1.9364916731037085*wTar[0]*fldDo[1]*eLo-1.9364916731037085*wDo[0]*fldDo[1]*eLo)/dxDo[0]-0.625*fldDo[2]*eUpR3+0.5590169943749475*fldDo[0]*eUpR3+0.625*fldDo[2]*eUp-0.5590169943749475*fldDo[0]*eUp+0.625*fldDo[2]*eLoR3-0.5590169943749475*fldDo[0]*eLoR3-0.625*fldDo[2]*eLo+0.5590169943749475*fldDo[0]*eLo); 
#else
  const double wDo0R2 = pow(wDo[0],2);
  const double wTar0R2 = pow(wTar[0],2);
  const double dxDo0R2 = pow(dxDo[0],2);
  const double dxTar0R2 = pow(dxTar[0],2);
  const double eLoR2 = pow(eLo,2);
  const double eLoR3 = pow(eLo,3);
  const double eLoR4 = pow(eLo,4);
  const double eLoR5 = pow(eLo,5);
  const double eUpR2 = pow(eUp,2);
  const double eUpR3 = pow(eUp,3);
  const double eUpR4 = pow(eUp,4);
  const double eUpR5 = pow(eUp,5);

  fldTar[0] += (dxTar0R2*(0.5590169943749475*fldDo[2]*eUpR3-0.5590169943749475*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(3.3541019662496847*wTar[0]*fldDo[2]*eUpR2-3.3541019662496847*wDo[0]*fldDo[2]*eUpR2-3.3541019662496847*wTar[0]*fldDo[2]*eLoR2+3.3541019662496847*wDo[0]*fldDo[2]*eLoR2))/dxDo0R2+(dxTar[0]*(0.4330127018922193*fldDo[1]*eUpR2-0.4330127018922193*fldDo[1]*eLoR2))/dxDo[0]+(6.708203932499369*wTar0R2*fldDo[2]*eUp-13.416407864998739*wDo[0]*wTar[0]*fldDo[2]*eUp+6.708203932499369*wDo0R2*fldDo[2]*eUp-6.708203932499369*wTar0R2*fldDo[2]*eLo+13.416407864998739*wDo[0]*wTar[0]*fldDo[2]*eLo-6.708203932499369*wDo0R2*fldDo[2]*eLo)/dxDo0R2+(1.7320508075688772*wTar[0]*fldDo[1]*eUp-1.7320508075688772*wDo[0]*fldDo[1]*eUp-1.7320508075688772*wTar[0]*fldDo[1]*eLo+1.7320508075688772*wDo[0]*fldDo[1]*eLo)/dxDo[0]-0.5590169943749475*fldDo[2]*eUp+0.5*fldDo[0]*eUp+0.5590169943749475*fldDo[2]*eLo-0.5*fldDo[0]*eLo; 
  fldTar[1] += (dxTar0R2*(0.7261843774138906*fldDo[2]*eUpR4-0.7261843774138906*fldDo[2]*eLoR4))/dxDo0R2+(dxTar[0]*(3.872983346207417*wTar[0]*fldDo[2]*eUpR3-3.872983346207417*wDo[0]*fldDo[2]*eUpR3-3.872983346207417*wTar[0]*fldDo[2]*eLoR3+3.872983346207417*wDo[0]*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(0.5*fldDo[1]*eUpR3-0.5*fldDo[1]*eLoR3))/dxDo[0]+(5.809475019311125*wTar0R2*fldDo[2]*eUpR2-11.61895003862225*wDo[0]*wTar[0]*fldDo[2]*eUpR2+5.809475019311125*wDo0R2*fldDo[2]*eUpR2-5.809475019311125*wTar0R2*fldDo[2]*eLoR2+11.61895003862225*wDo[0]*wTar[0]*fldDo[2]*eLoR2-5.809475019311125*wDo0R2*fldDo[2]*eLoR2)/dxDo0R2+(1.5*wTar[0]*fldDo[1]*eUpR2-1.5*wDo[0]*fldDo[1]*eUpR2-1.5*wTar[0]*fldDo[1]*eLoR2+1.5*wDo[0]*fldDo[1]*eLoR2)/dxDo[0]-0.4841229182759271*fldDo[2]*eUpR2+0.4330127018922193*fldDo[0]*eUpR2+0.4841229182759271*fldDo[2]*eLoR2-0.4330127018922193*fldDo[0]*eLoR2; 
  fldTar[2] += (dxTar0R2*(1.125*fldDo[2]*eUpR5-0.625*fldDo[2]*eUpR3-1.125*fldDo[2]*eLoR5+0.625*fldDo[2]*eLoR3))/dxDo0R2+(dxTar[0]*(5.625*wTar[0]*fldDo[2]*eUpR4-5.625*wDo[0]*fldDo[2]*eUpR4-3.75*wTar[0]*fldDo[2]*eUpR2+3.75*wDo[0]*fldDo[2]*eUpR2-5.625*wTar[0]*fldDo[2]*eLoR4+5.625*wDo[0]*fldDo[2]*eLoR4+3.75*wTar[0]*fldDo[2]*eLoR2-3.75*wDo[0]*fldDo[2]*eLoR2))/dxDo0R2+(dxTar[0]*(0.7261843774138906*fldDo[1]*eUpR4-0.4841229182759271*fldDo[1]*eUpR2-0.7261843774138906*fldDo[1]*eLoR4+0.4841229182759271*fldDo[1]*eLoR2))/dxDo[0]+(7.5*wTar0R2*fldDo[2]*eUpR3-15.0*wDo[0]*wTar[0]*fldDo[2]*eUpR3+7.5*wDo0R2*fldDo[2]*eUpR3-7.5*wTar0R2*fldDo[2]*eUp+15.0*wDo[0]*wTar[0]*fldDo[2]*eUp-7.5*wDo0R2*fldDo[2]*eUp-7.5*wTar0R2*fldDo[2]*eLoR3+15.0*wDo[0]*wTar[0]*fldDo[2]*eLoR3-7.5*wDo0R2*fldDo[2]*eLoR3+7.5*wTar0R2*fldDo[2]*eLo-15.0*wDo[0]*wTar[0]*fldDo[2]*eLo+7.5*wDo0R2*fldDo[2]*eLo)/dxDo0R2+(1.9364916731037085*wTar[0]*fldDo[1]*eUpR3-1.9364916731037085*wDo[0]*fldDo[1]*eUpR3-1.9364916731037085*wTar[0]*fldDo[1]*eUp+1.9364916731037085*wDo[0]*fldDo[1]*eUp-1.9364916731037085*wTar[0]*fldDo[1]*eLoR3+1.9364916731037085*wDo[0]*fldDo[1]*eLoR3+1.9364916731037085*wTar[0]*fldDo[1]*eLo-1.9364916731037085*wDo[0]*fldDo[1]*eLo)/dxDo[0]-0.625*fldDo[2]*eUpR3+0.5590169943749475*fldDo[0]*eUpR3+0.625*fldDo[2]*eUp-0.5590169943749475*fldDo[0]*eUp+0.625*fldDo[2]*eLoR3-0.5590169943749475*fldDo[0]*eLoR3-0.625*fldDo[2]*eLo+0.5590169943749475*fldDo[0]*eLo; 
#endif

}

