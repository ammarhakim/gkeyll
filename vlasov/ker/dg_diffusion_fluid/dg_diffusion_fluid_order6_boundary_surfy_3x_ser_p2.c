#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[8])+35.21807064562169*coeff[1]*fEdge[8]-34.09975027401226*coeff[1]*fSkin[2]-34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[12])+35.21807064562168*coeff[1]*fEdge[12]-34.09975027401226*coeff[1]*fSkin[4]-34.09975027401226*coeff[1]*fEdge[4]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(70.53065765632411*coeff[1]*fSkin[8])+51.46831774920949*coeff[1]*fEdge[8]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]-34.09975027401226*fSkin[0]*coeff[1]+34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[14])+35.21807064562168*coeff[1]*fEdge[14]-34.09975027401226*coeff[1]*fSkin[6]-34.09975027401226*coeff[1]*fEdge[6]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(70.53065765632414*coeff[1]*fSkin[12])+51.468317749209504*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[4]-56.6015625*coeff[1]*fEdge[4]-34.09975027401226*coeff[1]*fSkin[1]+34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[1]*fSkin[18])+35.21807064562169*coeff[1]*fEdge[18]-34.09975027401226*coeff[1]*fSkin[10]-34.09975027401226*coeff[1]*fEdge[10]-19.6875*coeff[1]*fSkin[5]+19.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = -(70.53065765632414*coeff[1]*fSkin[14])+51.468317749209504*coeff[1]*fEdge[14]-61.5234375*coeff[1]*fSkin[6]-56.6015625*coeff[1]*fEdge[6]-34.09975027401226*coeff[1]*fSkin[3]+34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = -(34.09975027401227*coeff[1]*fSkin[11])-34.09975027401227*coeff[1]*fEdge[11]-19.6875*coeff[1]*fSkin[7]+19.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[1]*fSkin[8])-3.1640625*coeff[1]*fEdge[8]-31.316701275974005*coeff[1]*fSkin[2]-12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = -(34.09975027401227*coeff[1]*fSkin[16])-34.09975027401227*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[9]+19.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = -(70.53065765632411*coeff[1]*fSkin[18])+51.46831774920949*coeff[1]*fEdge[18]-61.5234375*coeff[1]*fSkin[10]-56.6015625*coeff[1]*fEdge[10]-34.09975027401226*coeff[1]*fSkin[5]+34.09975027401226*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = -(61.5234375*coeff[1]*fSkin[11])-56.6015625*coeff[1]*fEdge[11]-34.09975027401227*coeff[1]*fSkin[7]+34.09975027401227*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]-31.316701275974033*coeff[1]*fSkin[4]-12.2543613688594*coeff[1]*fEdge[4]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = -(34.09975027401227*coeff[1]*fSkin[17])-34.09975027401227*coeff[1]*fEdge[17]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = -(70.6640625*coeff[1]*fSkin[14])-3.1640625*coeff[1]*fEdge[14]-31.316701275974033*coeff[1]*fSkin[6]-12.2543613688594*coeff[1]*fEdge[6]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = -(34.09975027401227*coeff[1]*fSkin[19])-34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[15]+19.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(61.5234375*coeff[1]*fSkin[16])-56.6015625*coeff[1]*fEdge[16]-34.09975027401227*coeff[1]*fSkin[9]+34.09975027401227*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = -(61.5234375*coeff[1]*fSkin[17])-56.6015625*coeff[1]*fEdge[17]-34.09975027401227*coeff[1]*fSkin[13]+34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = -(70.6640625*coeff[1]*fSkin[18])-3.1640625*coeff[1]*fEdge[18]-31.316701275974005*coeff[1]*fSkin[10]-12.2543613688594*coeff[1]*fEdge[10]-12.57788237343632*coeff[1]*fSkin[5]+12.57788237343632*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*fSkin[19])-56.6015625*coeff[1]*fEdge[19]-34.09975027401227*coeff[1]*fSkin[15]+34.09975027401227*coeff[1]*fEdge[15]; 

  boundSurf_incr[2] = 19.062339907114627*coeff[1]*fSkin[8]-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 19.062339907114634*coeff[1]*fSkin[12]-4.921875*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 19.062339907114634*coeff[1]*fSkin[14]-4.921875*coeff[1]*fSkin[6]; 
  boundSurf_incr[8] = 19.062339907114627*coeff[1]*fSkin[2]-73.828125*coeff[1]*fSkin[8]; 
  boundSurf_incr[10] = 19.062339907114627*coeff[1]*fSkin[18]-4.921875*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = -(4.921875*coeff[1]*fSkin[11]); 
  boundSurf_incr[12] = 19.062339907114634*coeff[1]*fSkin[4]-73.828125*coeff[1]*fSkin[12]; 
  boundSurf_incr[14] = 19.062339907114634*coeff[1]*fSkin[6]-73.828125*coeff[1]*fSkin[14]; 
  boundSurf_incr[16] = -(4.921875*coeff[1]*fSkin[16]); 
  boundSurf_incr[17] = -(4.921875*coeff[1]*fSkin[17]); 
  boundSurf_incr[18] = 19.062339907114627*coeff[1]*fSkin[10]-73.828125*coeff[1]*fSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*fSkin[19]); 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[8])+35.21807064562169*coeff[1]*fEdge[8]+34.09975027401226*coeff[1]*fSkin[2]+34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[12])+35.21807064562168*coeff[1]*fEdge[12]+34.09975027401226*coeff[1]*fSkin[4]+34.09975027401226*coeff[1]*fEdge[4]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 70.53065765632411*coeff[1]*fSkin[8]-51.46831774920949*coeff[1]*fEdge[8]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]+34.09975027401226*fSkin[0]*coeff[1]-34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[14])+35.21807064562168*coeff[1]*fEdge[14]+34.09975027401226*coeff[1]*fSkin[6]+34.09975027401226*coeff[1]*fEdge[6]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 70.53065765632414*coeff[1]*fSkin[12]-51.468317749209504*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[4]-56.6015625*coeff[1]*fEdge[4]+34.09975027401226*coeff[1]*fSkin[1]-34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[1]*fSkin[18])+35.21807064562169*coeff[1]*fEdge[18]+34.09975027401226*coeff[1]*fSkin[10]+34.09975027401226*coeff[1]*fEdge[10]-19.6875*coeff[1]*fSkin[5]+19.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = 70.53065765632414*coeff[1]*fSkin[14]-51.468317749209504*coeff[1]*fEdge[14]-61.5234375*coeff[1]*fSkin[6]-56.6015625*coeff[1]*fEdge[6]+34.09975027401226*coeff[1]*fSkin[3]-34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = 34.09975027401227*coeff[1]*fSkin[11]+34.09975027401227*coeff[1]*fEdge[11]-19.6875*coeff[1]*fSkin[7]+19.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[1]*fSkin[8])-3.1640625*coeff[1]*fEdge[8]+31.316701275974005*coeff[1]*fSkin[2]+12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = 34.09975027401227*coeff[1]*fSkin[16]+34.09975027401227*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[9]+19.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = 70.53065765632411*coeff[1]*fSkin[18]-51.46831774920949*coeff[1]*fEdge[18]-61.5234375*coeff[1]*fSkin[10]-56.6015625*coeff[1]*fEdge[10]+34.09975027401226*coeff[1]*fSkin[5]-34.09975027401226*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = -(61.5234375*coeff[1]*fSkin[11])-56.6015625*coeff[1]*fEdge[11]+34.09975027401227*coeff[1]*fSkin[7]-34.09975027401227*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]+31.316701275974033*coeff[1]*fSkin[4]+12.2543613688594*coeff[1]*fEdge[4]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = 34.09975027401227*coeff[1]*fSkin[17]+34.09975027401227*coeff[1]*fEdge[17]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = -(70.6640625*coeff[1]*fSkin[14])-3.1640625*coeff[1]*fEdge[14]+31.316701275974033*coeff[1]*fSkin[6]+12.2543613688594*coeff[1]*fEdge[6]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = 34.09975027401227*coeff[1]*fSkin[19]+34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[15]+19.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(61.5234375*coeff[1]*fSkin[16])-56.6015625*coeff[1]*fEdge[16]+34.09975027401227*coeff[1]*fSkin[9]-34.09975027401227*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = -(61.5234375*coeff[1]*fSkin[17])-56.6015625*coeff[1]*fEdge[17]+34.09975027401227*coeff[1]*fSkin[13]-34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = -(70.6640625*coeff[1]*fSkin[18])-3.1640625*coeff[1]*fEdge[18]+31.316701275974005*coeff[1]*fSkin[10]+12.2543613688594*coeff[1]*fEdge[10]-12.57788237343632*coeff[1]*fSkin[5]+12.57788237343632*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*fSkin[19])-56.6015625*coeff[1]*fEdge[19]+34.09975027401227*coeff[1]*fSkin[15]-34.09975027401227*coeff[1]*fEdge[15]; 

  boundSurf_incr[2] = -(19.062339907114627*coeff[1]*fSkin[8])-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = -(19.062339907114634*coeff[1]*fSkin[12])-4.921875*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = -(19.062339907114634*coeff[1]*fSkin[14])-4.921875*coeff[1]*fSkin[6]; 
  boundSurf_incr[8] = -(73.828125*coeff[1]*fSkin[8])-19.062339907114627*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = -(19.062339907114627*coeff[1]*fSkin[18])-4.921875*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = -(4.921875*coeff[1]*fSkin[11]); 
  boundSurf_incr[12] = -(73.828125*coeff[1]*fSkin[12])-19.062339907114634*coeff[1]*fSkin[4]; 
  boundSurf_incr[14] = -(73.828125*coeff[1]*fSkin[14])-19.062339907114634*coeff[1]*fSkin[6]; 
  boundSurf_incr[16] = -(4.921875*coeff[1]*fSkin[16]); 
  boundSurf_incr[17] = -(4.921875*coeff[1]*fSkin[17]); 
  boundSurf_incr[18] = -(73.828125*coeff[1]*fSkin[18])-19.062339907114627*coeff[1]*fSkin[10]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*fSkin[19]); 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  return 0.;
}

