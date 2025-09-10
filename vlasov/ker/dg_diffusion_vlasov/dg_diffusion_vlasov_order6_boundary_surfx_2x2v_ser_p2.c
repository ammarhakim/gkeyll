#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order6_boundary_surfx_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[11])+35.21807064562169*coeff[0]*fEdge[11]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(70.53065765632411*coeff[0]*fSkin[11])+51.46831774920949*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(35.21807064562168*coeff[0]*fSkin[19])+35.21807064562168*coeff[0]*fEdge[19]-34.09975027401226*coeff[0]*fSkin[5]-34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[0]*fSkin[21])+35.21807064562168*coeff[0]*fEdge[21]-34.09975027401226*coeff[0]*fSkin[6]-34.09975027401226*coeff[0]*fEdge[6]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[0]*fSkin[25])+35.21807064562168*coeff[0]*fEdge[25]-34.09975027401226*coeff[0]*fSkin[8]-34.09975027401226*coeff[0]*fEdge[8]-19.6875*coeff[0]*fSkin[4]+19.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(70.53065765632414*coeff[0]*fSkin[19])+51.468317749209504*coeff[0]*fEdge[19]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = -(70.53065765632414*coeff[0]*fSkin[21])+51.468317749209504*coeff[0]*fEdge[21]-61.5234375*coeff[0]*fSkin[6]-56.6015625*coeff[0]*fEdge[6]-34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = -(35.21807064562169*coeff[0]*fSkin[32])+35.21807064562169*coeff[0]*fEdge[32]-34.09975027401226*coeff[0]*fSkin[15]-34.09975027401226*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[7]+19.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = -(70.53065765632414*coeff[0]*fSkin[25])+51.468317749209504*coeff[0]*fEdge[25]-61.5234375*coeff[0]*fSkin[8]-56.6015625*coeff[0]*fEdge[8]-34.09975027401226*coeff[0]*fSkin[4]+34.09975027401226*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = -(35.21807064562169*coeff[0]*fSkin[35])+35.21807064562169*coeff[0]*fEdge[35]-34.09975027401226*coeff[0]*fSkin[16]-34.09975027401226*coeff[0]*fEdge[16]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[0]*fSkin[37])+35.21807064562169*coeff[0]*fEdge[37]-34.09975027401226*coeff[0]*fSkin[17]-34.09975027401226*coeff[0]*fEdge[17]-19.6875*coeff[0]*fSkin[10]+19.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]-31.316701275974005*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = -(34.09975027401227*coeff[0]*fSkin[20])-34.09975027401227*coeff[0]*fEdge[20]-19.6875*coeff[0]*fSkin[12]+19.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = -(34.09975027401227*coeff[0]*fSkin[23])-34.09975027401227*coeff[0]*fEdge[23]-19.6875*coeff[0]*fSkin[13]+19.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = -(34.09975027401227*coeff[0]*fSkin[28])-34.09975027401227*coeff[0]*fEdge[28]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = -(70.53065765632411*coeff[0]*fSkin[32])+51.46831774920949*coeff[0]*fEdge[32]-61.5234375*coeff[0]*fSkin[15]-56.6015625*coeff[0]*fEdge[15]-34.09975027401226*coeff[0]*fSkin[7]+34.09975027401226*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = -(70.53065765632411*coeff[0]*fSkin[35])+51.46831774920949*coeff[0]*fEdge[35]-61.5234375*coeff[0]*fSkin[16]-56.6015625*coeff[0]*fEdge[16]-34.09975027401226*coeff[0]*fSkin[9]+34.09975027401226*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = -(70.53065765632411*coeff[0]*fSkin[37])+51.46831774920949*coeff[0]*fEdge[37]-61.5234375*coeff[0]*fSkin[17]-56.6015625*coeff[0]*fEdge[17]-34.09975027401226*coeff[0]*fSkin[10]+34.09975027401226*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = -(35.21807064562168*coeff[0]*fSkin[44])+35.21807064562168*coeff[0]*fEdge[44]-34.09975027401226*coeff[0]*fSkin[31]-34.09975027401226*coeff[0]*fEdge[31]-19.6875*coeff[0]*fSkin[18]+19.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = -(70.6640625*coeff[0]*fSkin[19])-3.1640625*coeff[0]*fEdge[19]-31.316701275974033*coeff[0]*fSkin[5]-12.2543613688594*coeff[0]*fEdge[5]-12.577882373436315*coeff[0]*fSkin[2]+12.577882373436315*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = -(61.5234375*coeff[0]*fSkin[20])-56.6015625*coeff[0]*fEdge[20]-34.09975027401227*coeff[0]*fSkin[12]+34.09975027401227*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = -(70.6640625*coeff[0]*fSkin[21])-3.1640625*coeff[0]*fEdge[21]-31.316701275974033*coeff[0]*fSkin[6]-12.2543613688594*coeff[0]*fEdge[6]-12.577882373436315*coeff[0]*fSkin[3]+12.577882373436315*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = -(34.09975027401227*coeff[0]*fSkin[33])-34.09975027401227*coeff[0]*fEdge[33]-19.6875*coeff[0]*fSkin[22]+19.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(61.5234375*coeff[0]*fSkin[23])-56.6015625*coeff[0]*fEdge[23]-34.09975027401227*coeff[0]*fSkin[13]+34.09975027401227*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = -(34.09975027401227*coeff[0]*fSkin[34])-34.09975027401227*coeff[0]*fEdge[34]-19.6875*coeff[0]*fSkin[24]+19.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = -(70.6640625*coeff[0]*fSkin[25])-3.1640625*coeff[0]*fEdge[25]-31.316701275974033*coeff[0]*fSkin[8]-12.2543613688594*coeff[0]*fEdge[8]-12.577882373436315*coeff[0]*fSkin[4]+12.577882373436315*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = -(34.09975027401227*coeff[0]*fSkin[36])-34.09975027401227*coeff[0]*fEdge[36]-19.6875*coeff[0]*fSkin[26]+19.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = -(34.09975027401227*coeff[0]*fSkin[39])-34.09975027401227*coeff[0]*fEdge[39]-19.6875*coeff[0]*fSkin[27]+19.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = -(61.5234375*coeff[0]*fSkin[28])-56.6015625*coeff[0]*fEdge[28]-34.09975027401227*coeff[0]*fSkin[14]+34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = -(34.09975027401227*coeff[0]*fSkin[41])-34.09975027401227*coeff[0]*fEdge[41]-19.6875*coeff[0]*fSkin[29]+19.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = -(34.09975027401227*coeff[0]*fSkin[42])-34.09975027401227*coeff[0]*fEdge[42]-19.6875*coeff[0]*fSkin[30]+19.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = -(70.53065765632414*coeff[0]*fSkin[44])+51.468317749209504*coeff[0]*fEdge[44]-61.5234375*coeff[0]*fSkin[31]-56.6015625*coeff[0]*fEdge[31]-34.09975027401226*coeff[0]*fSkin[18]+34.09975027401226*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = -(70.6640625*coeff[0]*fSkin[32])-3.1640625*coeff[0]*fEdge[32]-31.316701275974005*coeff[0]*fSkin[15]-12.2543613688594*coeff[0]*fEdge[15]-12.57788237343632*coeff[0]*fSkin[7]+12.57788237343632*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = -(61.5234375*coeff[0]*fSkin[33])-56.6015625*coeff[0]*fEdge[33]-34.09975027401227*coeff[0]*fSkin[22]+34.09975027401227*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[0]*fSkin[34])-56.6015625*coeff[0]*fEdge[34]-34.09975027401227*coeff[0]*fSkin[24]+34.09975027401227*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = -(70.6640625*coeff[0]*fSkin[35])-3.1640625*coeff[0]*fEdge[35]-31.316701275974005*coeff[0]*fSkin[16]-12.2543613688594*coeff[0]*fEdge[16]-12.57788237343632*coeff[0]*fSkin[9]+12.57788237343632*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = -(61.5234375*coeff[0]*fSkin[36])-56.6015625*coeff[0]*fEdge[36]-34.09975027401227*coeff[0]*fSkin[26]+34.09975027401227*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = -(70.6640625*coeff[0]*fSkin[37])-3.1640625*coeff[0]*fEdge[37]-31.316701275974005*coeff[0]*fSkin[17]-12.2543613688594*coeff[0]*fEdge[17]-12.57788237343632*coeff[0]*fSkin[10]+12.57788237343632*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = -(34.09975027401227*coeff[0]*fSkin[45])-34.09975027401227*coeff[0]*fEdge[45]-19.6875*coeff[0]*fSkin[38]+19.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = -(61.5234375*coeff[0]*fSkin[39])-56.6015625*coeff[0]*fEdge[39]-34.09975027401227*coeff[0]*fSkin[27]+34.09975027401227*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = -(34.09975027401227*coeff[0]*fSkin[46])-34.09975027401227*coeff[0]*fEdge[46]-19.6875*coeff[0]*fSkin[40]+19.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[0]*fSkin[41])-56.6015625*coeff[0]*fEdge[41]-34.09975027401227*coeff[0]*fSkin[29]+34.09975027401227*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = -(61.5234375*coeff[0]*fSkin[42])-56.6015625*coeff[0]*fEdge[42]-34.09975027401227*coeff[0]*fSkin[30]+34.09975027401227*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = -(34.09975027401227*coeff[0]*fSkin[47])-34.09975027401227*coeff[0]*fEdge[47]-19.6875*coeff[0]*fSkin[43]+19.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = -(70.6640625*coeff[0]*fSkin[44])-3.1640625*coeff[0]*fEdge[44]-31.316701275974033*coeff[0]*fSkin[31]-12.2543613688594*coeff[0]*fEdge[31]-12.577882373436315*coeff[0]*fSkin[18]+12.577882373436315*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = -(61.5234375*coeff[0]*fSkin[45])-56.6015625*coeff[0]*fEdge[45]-34.09975027401227*coeff[0]*fSkin[38]+34.09975027401227*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[0]*fSkin[46])-56.6015625*coeff[0]*fEdge[46]-34.09975027401227*coeff[0]*fSkin[40]+34.09975027401227*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[0]*fSkin[47])-56.6015625*coeff[0]*fEdge[47]-34.09975027401227*coeff[0]*fSkin[43]+34.09975027401227*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 19.062339907114627*coeff[0]*fSkin[11]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 19.062339907114634*coeff[0]*fSkin[19]-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 19.062339907114634*coeff[0]*fSkin[21]-4.921875*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 19.062339907114634*coeff[0]*fSkin[25]-4.921875*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 19.062339907114627*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[11]; 
  boundSurf_incr[15] = 19.062339907114627*coeff[0]*fSkin[32]-4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 19.062339907114627*coeff[0]*fSkin[35]-4.921875*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 19.062339907114627*coeff[0]*fSkin[37]-4.921875*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = 19.062339907114634*coeff[0]*fSkin[5]-73.828125*coeff[0]*fSkin[19]; 
  boundSurf_incr[20] = -(4.921875*coeff[0]*fSkin[20]); 
  boundSurf_incr[21] = 19.062339907114634*coeff[0]*fSkin[6]-73.828125*coeff[0]*fSkin[21]; 
  boundSurf_incr[23] = -(4.921875*coeff[0]*fSkin[23]); 
  boundSurf_incr[25] = 19.062339907114634*coeff[0]*fSkin[8]-73.828125*coeff[0]*fSkin[25]; 
  boundSurf_incr[28] = -(4.921875*coeff[0]*fSkin[28]); 
  boundSurf_incr[31] = 19.062339907114634*coeff[0]*fSkin[44]-4.921875*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 19.062339907114627*coeff[0]*fSkin[15]-73.828125*coeff[0]*fSkin[32]; 
  boundSurf_incr[33] = -(4.921875*coeff[0]*fSkin[33]); 
  boundSurf_incr[34] = -(4.921875*coeff[0]*fSkin[34]); 
  boundSurf_incr[35] = 19.062339907114627*coeff[0]*fSkin[16]-73.828125*coeff[0]*fSkin[35]; 
  boundSurf_incr[36] = -(4.921875*coeff[0]*fSkin[36]); 
  boundSurf_incr[37] = 19.062339907114627*coeff[0]*fSkin[17]-73.828125*coeff[0]*fSkin[37]; 
  boundSurf_incr[39] = -(4.921875*coeff[0]*fSkin[39]); 
  boundSurf_incr[41] = -(4.921875*coeff[0]*fSkin[41]); 
  boundSurf_incr[42] = -(4.921875*coeff[0]*fSkin[42]); 
  boundSurf_incr[44] = 19.062339907114634*coeff[0]*fSkin[31]-73.828125*coeff[0]*fSkin[44]; 
  boundSurf_incr[45] = -(4.921875*coeff[0]*fSkin[45]); 
  boundSurf_incr[46] = -(4.921875*coeff[0]*fSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[0]*fSkin[47]); 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[0]*fSkin[11])+35.21807064562169*coeff[0]*fEdge[11]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[11]-51.46831774920949*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(35.21807064562168*coeff[0]*fSkin[19])+35.21807064562168*coeff[0]*fEdge[19]+34.09975027401226*coeff[0]*fSkin[5]+34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[0]*fSkin[21])+35.21807064562168*coeff[0]*fEdge[21]+34.09975027401226*coeff[0]*fSkin[6]+34.09975027401226*coeff[0]*fEdge[6]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[0]*fSkin[25])+35.21807064562168*coeff[0]*fEdge[25]+34.09975027401226*coeff[0]*fSkin[8]+34.09975027401226*coeff[0]*fEdge[8]-19.6875*coeff[0]*fSkin[4]+19.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[0]*fSkin[19]-51.468317749209504*coeff[0]*fEdge[19]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 70.53065765632414*coeff[0]*fSkin[21]-51.468317749209504*coeff[0]*fEdge[21]-61.5234375*coeff[0]*fSkin[6]-56.6015625*coeff[0]*fEdge[6]+34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = -(35.21807064562169*coeff[0]*fSkin[32])+35.21807064562169*coeff[0]*fEdge[32]+34.09975027401226*coeff[0]*fSkin[15]+34.09975027401226*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[7]+19.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 70.53065765632414*coeff[0]*fSkin[25]-51.468317749209504*coeff[0]*fEdge[25]-61.5234375*coeff[0]*fSkin[8]-56.6015625*coeff[0]*fEdge[8]+34.09975027401226*coeff[0]*fSkin[4]-34.09975027401226*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = -(35.21807064562169*coeff[0]*fSkin[35])+35.21807064562169*coeff[0]*fEdge[35]+34.09975027401226*coeff[0]*fSkin[16]+34.09975027401226*coeff[0]*fEdge[16]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[0]*fSkin[37])+35.21807064562169*coeff[0]*fEdge[37]+34.09975027401226*coeff[0]*fSkin[17]+34.09975027401226*coeff[0]*fEdge[17]-19.6875*coeff[0]*fSkin[10]+19.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]+31.316701275974005*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 34.09975027401227*coeff[0]*fSkin[20]+34.09975027401227*coeff[0]*fEdge[20]-19.6875*coeff[0]*fSkin[12]+19.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 34.09975027401227*coeff[0]*fSkin[23]+34.09975027401227*coeff[0]*fEdge[23]-19.6875*coeff[0]*fSkin[13]+19.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[0]*fSkin[28]+34.09975027401227*coeff[0]*fEdge[28]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 70.53065765632411*coeff[0]*fSkin[32]-51.46831774920949*coeff[0]*fEdge[32]-61.5234375*coeff[0]*fSkin[15]-56.6015625*coeff[0]*fEdge[15]+34.09975027401226*coeff[0]*fSkin[7]-34.09975027401226*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 70.53065765632411*coeff[0]*fSkin[35]-51.46831774920949*coeff[0]*fEdge[35]-61.5234375*coeff[0]*fSkin[16]-56.6015625*coeff[0]*fEdge[16]+34.09975027401226*coeff[0]*fSkin[9]-34.09975027401226*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 70.53065765632411*coeff[0]*fSkin[37]-51.46831774920949*coeff[0]*fEdge[37]-61.5234375*coeff[0]*fSkin[17]-56.6015625*coeff[0]*fEdge[17]+34.09975027401226*coeff[0]*fSkin[10]-34.09975027401226*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = -(35.21807064562168*coeff[0]*fSkin[44])+35.21807064562168*coeff[0]*fEdge[44]+34.09975027401226*coeff[0]*fSkin[31]+34.09975027401226*coeff[0]*fEdge[31]-19.6875*coeff[0]*fSkin[18]+19.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = -(70.6640625*coeff[0]*fSkin[19])-3.1640625*coeff[0]*fEdge[19]+31.316701275974033*coeff[0]*fSkin[5]+12.2543613688594*coeff[0]*fEdge[5]-12.577882373436315*coeff[0]*fSkin[2]+12.577882373436315*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = -(61.5234375*coeff[0]*fSkin[20])-56.6015625*coeff[0]*fEdge[20]+34.09975027401227*coeff[0]*fSkin[12]-34.09975027401227*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = -(70.6640625*coeff[0]*fSkin[21])-3.1640625*coeff[0]*fEdge[21]+31.316701275974033*coeff[0]*fSkin[6]+12.2543613688594*coeff[0]*fEdge[6]-12.577882373436315*coeff[0]*fSkin[3]+12.577882373436315*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 34.09975027401227*coeff[0]*fSkin[33]+34.09975027401227*coeff[0]*fEdge[33]-19.6875*coeff[0]*fSkin[22]+19.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(61.5234375*coeff[0]*fSkin[23])-56.6015625*coeff[0]*fEdge[23]+34.09975027401227*coeff[0]*fSkin[13]-34.09975027401227*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 34.09975027401227*coeff[0]*fSkin[34]+34.09975027401227*coeff[0]*fEdge[34]-19.6875*coeff[0]*fSkin[24]+19.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = -(70.6640625*coeff[0]*fSkin[25])-3.1640625*coeff[0]*fEdge[25]+31.316701275974033*coeff[0]*fSkin[8]+12.2543613688594*coeff[0]*fEdge[8]-12.577882373436315*coeff[0]*fSkin[4]+12.577882373436315*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 34.09975027401227*coeff[0]*fSkin[36]+34.09975027401227*coeff[0]*fEdge[36]-19.6875*coeff[0]*fSkin[26]+19.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 34.09975027401227*coeff[0]*fSkin[39]+34.09975027401227*coeff[0]*fEdge[39]-19.6875*coeff[0]*fSkin[27]+19.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = -(61.5234375*coeff[0]*fSkin[28])-56.6015625*coeff[0]*fEdge[28]+34.09975027401227*coeff[0]*fSkin[14]-34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 34.09975027401227*coeff[0]*fSkin[41]+34.09975027401227*coeff[0]*fEdge[41]-19.6875*coeff[0]*fSkin[29]+19.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 34.09975027401227*coeff[0]*fSkin[42]+34.09975027401227*coeff[0]*fEdge[42]-19.6875*coeff[0]*fSkin[30]+19.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 70.53065765632414*coeff[0]*fSkin[44]-51.468317749209504*coeff[0]*fEdge[44]-61.5234375*coeff[0]*fSkin[31]-56.6015625*coeff[0]*fEdge[31]+34.09975027401226*coeff[0]*fSkin[18]-34.09975027401226*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = -(70.6640625*coeff[0]*fSkin[32])-3.1640625*coeff[0]*fEdge[32]+31.316701275974005*coeff[0]*fSkin[15]+12.2543613688594*coeff[0]*fEdge[15]-12.57788237343632*coeff[0]*fSkin[7]+12.57788237343632*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = -(61.5234375*coeff[0]*fSkin[33])-56.6015625*coeff[0]*fEdge[33]+34.09975027401227*coeff[0]*fSkin[22]-34.09975027401227*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[0]*fSkin[34])-56.6015625*coeff[0]*fEdge[34]+34.09975027401227*coeff[0]*fSkin[24]-34.09975027401227*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = -(70.6640625*coeff[0]*fSkin[35])-3.1640625*coeff[0]*fEdge[35]+31.316701275974005*coeff[0]*fSkin[16]+12.2543613688594*coeff[0]*fEdge[16]-12.57788237343632*coeff[0]*fSkin[9]+12.57788237343632*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = -(61.5234375*coeff[0]*fSkin[36])-56.6015625*coeff[0]*fEdge[36]+34.09975027401227*coeff[0]*fSkin[26]-34.09975027401227*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = -(70.6640625*coeff[0]*fSkin[37])-3.1640625*coeff[0]*fEdge[37]+31.316701275974005*coeff[0]*fSkin[17]+12.2543613688594*coeff[0]*fEdge[17]-12.57788237343632*coeff[0]*fSkin[10]+12.57788237343632*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 34.09975027401227*coeff[0]*fSkin[45]+34.09975027401227*coeff[0]*fEdge[45]-19.6875*coeff[0]*fSkin[38]+19.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = -(61.5234375*coeff[0]*fSkin[39])-56.6015625*coeff[0]*fEdge[39]+34.09975027401227*coeff[0]*fSkin[27]-34.09975027401227*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 34.09975027401227*coeff[0]*fSkin[46]+34.09975027401227*coeff[0]*fEdge[46]-19.6875*coeff[0]*fSkin[40]+19.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[0]*fSkin[41])-56.6015625*coeff[0]*fEdge[41]+34.09975027401227*coeff[0]*fSkin[29]-34.09975027401227*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = -(61.5234375*coeff[0]*fSkin[42])-56.6015625*coeff[0]*fEdge[42]+34.09975027401227*coeff[0]*fSkin[30]-34.09975027401227*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 34.09975027401227*coeff[0]*fSkin[47]+34.09975027401227*coeff[0]*fEdge[47]-19.6875*coeff[0]*fSkin[43]+19.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = -(70.6640625*coeff[0]*fSkin[44])-3.1640625*coeff[0]*fEdge[44]+31.316701275974033*coeff[0]*fSkin[31]+12.2543613688594*coeff[0]*fEdge[31]-12.577882373436315*coeff[0]*fSkin[18]+12.577882373436315*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = -(61.5234375*coeff[0]*fSkin[45])-56.6015625*coeff[0]*fEdge[45]+34.09975027401227*coeff[0]*fSkin[38]-34.09975027401227*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[0]*fSkin[46])-56.6015625*coeff[0]*fEdge[46]+34.09975027401227*coeff[0]*fSkin[40]-34.09975027401227*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[0]*fSkin[47])-56.6015625*coeff[0]*fEdge[47]+34.09975027401227*coeff[0]*fSkin[43]-34.09975027401227*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = -(19.062339907114627*coeff[0]*fSkin[11])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = -(19.062339907114634*coeff[0]*fSkin[19])-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = -(19.062339907114634*coeff[0]*fSkin[21])-4.921875*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = -(19.062339907114634*coeff[0]*fSkin[25])-4.921875*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = -(73.828125*coeff[0]*fSkin[11])-19.062339907114627*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = -(19.062339907114627*coeff[0]*fSkin[32])-4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = -(19.062339907114627*coeff[0]*fSkin[35])-4.921875*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = -(19.062339907114627*coeff[0]*fSkin[37])-4.921875*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = -(73.828125*coeff[0]*fSkin[19])-19.062339907114634*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = -(4.921875*coeff[0]*fSkin[20]); 
  boundSurf_incr[21] = -(73.828125*coeff[0]*fSkin[21])-19.062339907114634*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = -(4.921875*coeff[0]*fSkin[23]); 
  boundSurf_incr[25] = -(73.828125*coeff[0]*fSkin[25])-19.062339907114634*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = -(4.921875*coeff[0]*fSkin[28]); 
  boundSurf_incr[31] = -(19.062339907114634*coeff[0]*fSkin[44])-4.921875*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = -(73.828125*coeff[0]*fSkin[32])-19.062339907114627*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = -(4.921875*coeff[0]*fSkin[33]); 
  boundSurf_incr[34] = -(4.921875*coeff[0]*fSkin[34]); 
  boundSurf_incr[35] = -(73.828125*coeff[0]*fSkin[35])-19.062339907114627*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = -(4.921875*coeff[0]*fSkin[36]); 
  boundSurf_incr[37] = -(73.828125*coeff[0]*fSkin[37])-19.062339907114627*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = -(4.921875*coeff[0]*fSkin[39]); 
  boundSurf_incr[41] = -(4.921875*coeff[0]*fSkin[41]); 
  boundSurf_incr[42] = -(4.921875*coeff[0]*fSkin[42]); 
  boundSurf_incr[44] = -(73.828125*coeff[0]*fSkin[44])-19.062339907114634*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = -(4.921875*coeff[0]*fSkin[45]); 
  boundSurf_incr[46] = -(4.921875*coeff[0]*fSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[0]*fSkin[47]); 

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
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 

  return 0.;
}

