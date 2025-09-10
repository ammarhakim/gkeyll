#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_3x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  double vol_incr[27] = {0.0}; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[8])+35.21807064562169*coeff[1]*fEdge[8]-34.09975027401226*coeff[1]*fSkin[2]-34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[12])+35.21807064562168*coeff[1]*fEdge[12]-34.09975027401226*coeff[1]*fSkin[4]-34.09975027401226*coeff[1]*fEdge[4]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(70.53065765632411*coeff[1]*fSkin[8])+51.46831774920949*coeff[1]*fEdge[8]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]-34.09975027401226*fSkin[0]*coeff[1]+34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[14])+35.21807064562168*coeff[1]*fEdge[14]-34.09975027401226*coeff[1]*fSkin[6]-34.09975027401226*coeff[1]*fEdge[6]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(70.53065765632414*coeff[1]*fSkin[12])+51.468317749209504*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[4]-56.6015625*coeff[1]*fEdge[4]-34.09975027401226*coeff[1]*fSkin[1]+34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[1]*fSkin[18])+35.21807064562169*coeff[1]*fEdge[18]-34.09975027401226*coeff[1]*fSkin[10]-34.09975027401226*coeff[1]*fEdge[10]-19.6875*coeff[1]*fSkin[5]+19.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = -(70.53065765632414*coeff[1]*fSkin[14])+51.468317749209504*coeff[1]*fEdge[14]-61.5234375*coeff[1]*fSkin[6]-56.6015625*coeff[1]*fEdge[6]-34.09975027401226*coeff[1]*fSkin[3]+34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = -(35.21807064562169*coeff[1]*fSkin[20])+35.21807064562169*coeff[1]*fEdge[20]-34.09975027401227*coeff[1]*fSkin[11]-34.09975027401227*coeff[1]*fEdge[11]-19.6875*coeff[1]*fSkin[7]+19.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[1]*fSkin[8])-3.1640625*coeff[1]*fEdge[8]-31.316701275974005*coeff[1]*fSkin[2]-12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = -(35.21807064562169*coeff[1]*fSkin[22])+35.21807064562169*coeff[1]*fEdge[22]-34.09975027401227*coeff[1]*fSkin[16]-34.09975027401227*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[9]+19.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = -(70.53065765632411*coeff[1]*fSkin[18])+51.46831774920949*coeff[1]*fEdge[18]-61.5234375*coeff[1]*fSkin[10]-56.6015625*coeff[1]*fEdge[10]-34.09975027401226*coeff[1]*fSkin[5]+34.09975027401226*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = -(70.53065765632414*coeff[1]*fSkin[20])+51.468317749209504*coeff[1]*fEdge[20]-61.5234375*coeff[1]*fSkin[11]-56.6015625*coeff[1]*fEdge[11]-34.09975027401227*coeff[1]*fSkin[7]+34.09975027401227*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]-31.316701275974033*coeff[1]*fSkin[4]-12.2543613688594*coeff[1]*fEdge[4]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = -(35.21807064562168*coeff[1]*fSkin[23])+35.21807064562168*coeff[1]*fEdge[23]-34.09975027401227*coeff[1]*fSkin[17]-34.09975027401227*coeff[1]*fEdge[17]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = -(70.6640625*coeff[1]*fSkin[14])-3.1640625*coeff[1]*fEdge[14]-31.316701275974033*coeff[1]*fSkin[6]-12.2543613688594*coeff[1]*fEdge[6]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = -(35.21807064562168*coeff[1]*fSkin[25])+35.21807064562168*coeff[1]*fEdge[25]-34.09975027401227*coeff[1]*fSkin[19]-34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[15]+19.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(70.53065765632414*coeff[1]*fSkin[22])+51.468317749209504*coeff[1]*fEdge[22]-61.5234375*coeff[1]*fSkin[16]-56.6015625*coeff[1]*fEdge[16]-34.09975027401227*coeff[1]*fSkin[9]+34.09975027401227*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = -(70.53065765632411*coeff[1]*fSkin[23])+51.46831774920949*coeff[1]*fEdge[23]-61.5234375*coeff[1]*fSkin[17]-56.6015625*coeff[1]*fEdge[17]-34.09975027401227*coeff[1]*fSkin[13]+34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = -(70.6640625*coeff[1]*fSkin[18])-3.1640625*coeff[1]*fEdge[18]-31.316701275974005*coeff[1]*fSkin[10]-12.2543613688594*coeff[1]*fEdge[10]-12.57788237343632*coeff[1]*fSkin[5]+12.57788237343632*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = -(70.53065765632411*coeff[1]*fSkin[25])+51.46831774920949*coeff[1]*fEdge[25]-61.5234375*coeff[1]*fSkin[19]-56.6015625*coeff[1]*fEdge[19]-34.09975027401227*coeff[1]*fSkin[15]+34.09975027401227*coeff[1]*fEdge[15]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*fSkin[20])-3.1640625*coeff[1]*fEdge[20]-31.316701275974033*coeff[1]*fSkin[11]-12.2543613688594*coeff[1]*fEdge[11]-12.57788237343632*coeff[1]*fSkin[7]+12.57788237343632*coeff[1]*fEdge[7]; 
  edgeSurf_incr[21] = -(35.21807064562169*coeff[1]*fSkin[26])+35.21807064562169*coeff[1]*fEdge[26]-34.09975027401226*coeff[1]*fSkin[24]-34.09975027401226*coeff[1]*fEdge[24]-19.6875*coeff[1]*fSkin[21]+19.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*fSkin[22])-3.1640625*coeff[1]*fEdge[22]-31.316701275974033*coeff[1]*fSkin[16]-12.2543613688594*coeff[1]*fEdge[16]-12.57788237343632*coeff[1]*fSkin[9]+12.57788237343632*coeff[1]*fEdge[9]; 
  edgeSurf_incr[23] = -(70.6640625*coeff[1]*fSkin[23])-3.1640625*coeff[1]*fEdge[23]-31.316701275974005*coeff[1]*fSkin[17]-12.2543613688594*coeff[1]*fEdge[17]-12.577882373436315*coeff[1]*fSkin[13]+12.577882373436315*coeff[1]*fEdge[13]; 
  edgeSurf_incr[24] = -(70.53065765632411*coeff[1]*fSkin[26])+51.46831774920949*coeff[1]*fEdge[26]-61.5234375*coeff[1]*fSkin[24]-56.6015625*coeff[1]*fEdge[24]-34.09975027401226*coeff[1]*fSkin[21]+34.09975027401226*coeff[1]*fEdge[21]; 
  edgeSurf_incr[25] = -(70.6640625*coeff[1]*fSkin[25])-3.1640625*coeff[1]*fEdge[25]-31.316701275974005*coeff[1]*fSkin[19]-12.2543613688594*coeff[1]*fEdge[19]-12.577882373436315*coeff[1]*fSkin[15]+12.577882373436315*coeff[1]*fEdge[15]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*fSkin[26])-3.1640625*coeff[1]*fEdge[26]-31.316701275974005*coeff[1]*fSkin[24]-12.2543613688594*coeff[1]*fEdge[24]-12.57788237343632*coeff[1]*fSkin[21]+12.57788237343632*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = 19.062339907114627*coeff[1]*fSkin[8]-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 19.062339907114634*coeff[1]*fSkin[12]-4.921875*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 19.062339907114634*coeff[1]*fSkin[14]-4.921875*coeff[1]*fSkin[6]; 
  boundSurf_incr[8] = 19.062339907114627*coeff[1]*fSkin[2]-73.828125*coeff[1]*fSkin[8]; 
  boundSurf_incr[10] = 19.062339907114627*coeff[1]*fSkin[18]-4.921875*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = 19.062339907114634*coeff[1]*fSkin[20]-4.921875*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = 19.062339907114634*coeff[1]*fSkin[4]-73.828125*coeff[1]*fSkin[12]; 
  boundSurf_incr[14] = 19.062339907114634*coeff[1]*fSkin[6]-73.828125*coeff[1]*fSkin[14]; 
  boundSurf_incr[16] = 19.062339907114634*coeff[1]*fSkin[22]-4.921875*coeff[1]*fSkin[16]; 
  boundSurf_incr[17] = 19.062339907114627*coeff[1]*fSkin[23]-4.921875*coeff[1]*fSkin[17]; 
  boundSurf_incr[18] = 19.062339907114627*coeff[1]*fSkin[10]-73.828125*coeff[1]*fSkin[18]; 
  boundSurf_incr[19] = 19.062339907114627*coeff[1]*fSkin[25]-4.921875*coeff[1]*fSkin[19]; 
  boundSurf_incr[20] = 19.062339907114634*coeff[1]*fSkin[11]-73.828125*coeff[1]*fSkin[20]; 
  boundSurf_incr[22] = 19.062339907114634*coeff[1]*fSkin[16]-73.828125*coeff[1]*fSkin[22]; 
  boundSurf_incr[23] = 19.062339907114627*coeff[1]*fSkin[17]-73.828125*coeff[1]*fSkin[23]; 
  boundSurf_incr[24] = 19.062339907114627*coeff[1]*fSkin[26]-4.921875*coeff[1]*fSkin[24]; 
  boundSurf_incr[25] = 19.062339907114627*coeff[1]*fSkin[19]-73.828125*coeff[1]*fSkin[25]; 
  boundSurf_incr[26] = 19.062339907114627*coeff[1]*fSkin[24]-73.828125*coeff[1]*fSkin[26]; 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[8])+35.21807064562169*coeff[1]*fEdge[8]+34.09975027401226*coeff[1]*fSkin[2]+34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[12])+35.21807064562168*coeff[1]*fEdge[12]+34.09975027401226*coeff[1]*fSkin[4]+34.09975027401226*coeff[1]*fEdge[4]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 70.53065765632411*coeff[1]*fSkin[8]-51.46831774920949*coeff[1]*fEdge[8]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]+34.09975027401226*fSkin[0]*coeff[1]-34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[14])+35.21807064562168*coeff[1]*fEdge[14]+34.09975027401226*coeff[1]*fSkin[6]+34.09975027401226*coeff[1]*fEdge[6]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 70.53065765632414*coeff[1]*fSkin[12]-51.468317749209504*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[4]-56.6015625*coeff[1]*fEdge[4]+34.09975027401226*coeff[1]*fSkin[1]-34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = -(35.21807064562169*coeff[1]*fSkin[18])+35.21807064562169*coeff[1]*fEdge[18]+34.09975027401226*coeff[1]*fSkin[10]+34.09975027401226*coeff[1]*fEdge[10]-19.6875*coeff[1]*fSkin[5]+19.6875*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = 70.53065765632414*coeff[1]*fSkin[14]-51.468317749209504*coeff[1]*fEdge[14]-61.5234375*coeff[1]*fSkin[6]-56.6015625*coeff[1]*fEdge[6]+34.09975027401226*coeff[1]*fSkin[3]-34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = -(35.21807064562169*coeff[1]*fSkin[20])+35.21807064562169*coeff[1]*fEdge[20]+34.09975027401227*coeff[1]*fSkin[11]+34.09975027401227*coeff[1]*fEdge[11]-19.6875*coeff[1]*fSkin[7]+19.6875*coeff[1]*fEdge[7]; 
  edgeSurf_incr[8] = -(70.6640625*coeff[1]*fSkin[8])-3.1640625*coeff[1]*fEdge[8]+31.316701275974005*coeff[1]*fSkin[2]+12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[9] = -(35.21807064562169*coeff[1]*fSkin[22])+35.21807064562169*coeff[1]*fEdge[22]+34.09975027401227*coeff[1]*fSkin[16]+34.09975027401227*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[9]+19.6875*coeff[1]*fEdge[9]; 
  edgeSurf_incr[10] = 70.53065765632411*coeff[1]*fSkin[18]-51.46831774920949*coeff[1]*fEdge[18]-61.5234375*coeff[1]*fSkin[10]-56.6015625*coeff[1]*fEdge[10]+34.09975027401226*coeff[1]*fSkin[5]-34.09975027401226*coeff[1]*fEdge[5]; 
  edgeSurf_incr[11] = 70.53065765632414*coeff[1]*fSkin[20]-51.468317749209504*coeff[1]*fEdge[20]-61.5234375*coeff[1]*fSkin[11]-56.6015625*coeff[1]*fEdge[11]+34.09975027401227*coeff[1]*fSkin[7]-34.09975027401227*coeff[1]*fEdge[7]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]+31.316701275974033*coeff[1]*fSkin[4]+12.2543613688594*coeff[1]*fEdge[4]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = -(35.21807064562168*coeff[1]*fSkin[23])+35.21807064562168*coeff[1]*fEdge[23]+34.09975027401227*coeff[1]*fSkin[17]+34.09975027401227*coeff[1]*fEdge[17]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = -(70.6640625*coeff[1]*fSkin[14])-3.1640625*coeff[1]*fEdge[14]+31.316701275974033*coeff[1]*fSkin[6]+12.2543613688594*coeff[1]*fEdge[6]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[15] = -(35.21807064562168*coeff[1]*fSkin[25])+35.21807064562168*coeff[1]*fEdge[25]+34.09975027401227*coeff[1]*fSkin[19]+34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[15]+19.6875*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = 70.53065765632414*coeff[1]*fSkin[22]-51.468317749209504*coeff[1]*fEdge[22]-61.5234375*coeff[1]*fSkin[16]-56.6015625*coeff[1]*fEdge[16]+34.09975027401227*coeff[1]*fSkin[9]-34.09975027401227*coeff[1]*fEdge[9]; 
  edgeSurf_incr[17] = 70.53065765632411*coeff[1]*fSkin[23]-51.46831774920949*coeff[1]*fEdge[23]-61.5234375*coeff[1]*fSkin[17]-56.6015625*coeff[1]*fEdge[17]+34.09975027401227*coeff[1]*fSkin[13]-34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[18] = -(70.6640625*coeff[1]*fSkin[18])-3.1640625*coeff[1]*fEdge[18]+31.316701275974005*coeff[1]*fSkin[10]+12.2543613688594*coeff[1]*fEdge[10]-12.57788237343632*coeff[1]*fSkin[5]+12.57788237343632*coeff[1]*fEdge[5]; 
  edgeSurf_incr[19] = 70.53065765632411*coeff[1]*fSkin[25]-51.46831774920949*coeff[1]*fEdge[25]-61.5234375*coeff[1]*fSkin[19]-56.6015625*coeff[1]*fEdge[19]+34.09975027401227*coeff[1]*fSkin[15]-34.09975027401227*coeff[1]*fEdge[15]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*fSkin[20])-3.1640625*coeff[1]*fEdge[20]+31.316701275974033*coeff[1]*fSkin[11]+12.2543613688594*coeff[1]*fEdge[11]-12.57788237343632*coeff[1]*fSkin[7]+12.57788237343632*coeff[1]*fEdge[7]; 
  edgeSurf_incr[21] = -(35.21807064562169*coeff[1]*fSkin[26])+35.21807064562169*coeff[1]*fEdge[26]+34.09975027401226*coeff[1]*fSkin[24]+34.09975027401226*coeff[1]*fEdge[24]-19.6875*coeff[1]*fSkin[21]+19.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*fSkin[22])-3.1640625*coeff[1]*fEdge[22]+31.316701275974033*coeff[1]*fSkin[16]+12.2543613688594*coeff[1]*fEdge[16]-12.57788237343632*coeff[1]*fSkin[9]+12.57788237343632*coeff[1]*fEdge[9]; 
  edgeSurf_incr[23] = -(70.6640625*coeff[1]*fSkin[23])-3.1640625*coeff[1]*fEdge[23]+31.316701275974005*coeff[1]*fSkin[17]+12.2543613688594*coeff[1]*fEdge[17]-12.577882373436315*coeff[1]*fSkin[13]+12.577882373436315*coeff[1]*fEdge[13]; 
  edgeSurf_incr[24] = 70.53065765632411*coeff[1]*fSkin[26]-51.46831774920949*coeff[1]*fEdge[26]-61.5234375*coeff[1]*fSkin[24]-56.6015625*coeff[1]*fEdge[24]+34.09975027401226*coeff[1]*fSkin[21]-34.09975027401226*coeff[1]*fEdge[21]; 
  edgeSurf_incr[25] = -(70.6640625*coeff[1]*fSkin[25])-3.1640625*coeff[1]*fEdge[25]+31.316701275974005*coeff[1]*fSkin[19]+12.2543613688594*coeff[1]*fEdge[19]-12.577882373436315*coeff[1]*fSkin[15]+12.577882373436315*coeff[1]*fEdge[15]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*fSkin[26])-3.1640625*coeff[1]*fEdge[26]+31.316701275974005*coeff[1]*fSkin[24]+12.2543613688594*coeff[1]*fEdge[24]-12.57788237343632*coeff[1]*fSkin[21]+12.57788237343632*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = -(19.062339907114627*coeff[1]*fSkin[8])-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = -(19.062339907114634*coeff[1]*fSkin[12])-4.921875*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = -(19.062339907114634*coeff[1]*fSkin[14])-4.921875*coeff[1]*fSkin[6]; 
  boundSurf_incr[8] = -(73.828125*coeff[1]*fSkin[8])-19.062339907114627*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = -(19.062339907114627*coeff[1]*fSkin[18])-4.921875*coeff[1]*fSkin[10]; 
  boundSurf_incr[11] = -(19.062339907114634*coeff[1]*fSkin[20])-4.921875*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = -(73.828125*coeff[1]*fSkin[12])-19.062339907114634*coeff[1]*fSkin[4]; 
  boundSurf_incr[14] = -(73.828125*coeff[1]*fSkin[14])-19.062339907114634*coeff[1]*fSkin[6]; 
  boundSurf_incr[16] = -(19.062339907114634*coeff[1]*fSkin[22])-4.921875*coeff[1]*fSkin[16]; 
  boundSurf_incr[17] = -(19.062339907114627*coeff[1]*fSkin[23])-4.921875*coeff[1]*fSkin[17]; 
  boundSurf_incr[18] = -(73.828125*coeff[1]*fSkin[18])-19.062339907114627*coeff[1]*fSkin[10]; 
  boundSurf_incr[19] = -(19.062339907114627*coeff[1]*fSkin[25])-4.921875*coeff[1]*fSkin[19]; 
  boundSurf_incr[20] = -(73.828125*coeff[1]*fSkin[20])-19.062339907114634*coeff[1]*fSkin[11]; 
  boundSurf_incr[22] = -(73.828125*coeff[1]*fSkin[22])-19.062339907114634*coeff[1]*fSkin[16]; 
  boundSurf_incr[23] = -(73.828125*coeff[1]*fSkin[23])-19.062339907114627*coeff[1]*fSkin[17]; 
  boundSurf_incr[24] = -(19.062339907114627*coeff[1]*fSkin[26])-4.921875*coeff[1]*fSkin[24]; 
  boundSurf_incr[25] = -(73.828125*coeff[1]*fSkin[25])-19.062339907114627*coeff[1]*fSkin[19]; 
  boundSurf_incr[26] = -(73.828125*coeff[1]*fSkin[26])-19.062339907114627*coeff[1]*fSkin[24]; 

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
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 

  return 0.;
}

