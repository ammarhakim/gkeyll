#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order6_boundary_surfy_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[12])+35.21807064562169*coeff[1]*fEdge[12]-34.09975027401226*coeff[1]*fSkin[2]-34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[20])+35.21807064562168*coeff[1]*fEdge[20]-34.09975027401226*coeff[1]*fSkin[5]-34.09975027401226*coeff[1]*fEdge[5]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(70.53065765632411*coeff[1]*fSkin[12])+51.46831774920949*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]-34.09975027401226*fSkin[0]*coeff[1]+34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[22])+35.21807064562168*coeff[1]*fEdge[22]-34.09975027401226*coeff[1]*fSkin[7]-34.09975027401226*coeff[1]*fEdge[7]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[1]*fSkin[26])+35.21807064562168*coeff[1]*fEdge[26]-34.09975027401226*coeff[1]*fSkin[9]-34.09975027401226*coeff[1]*fEdge[9]-19.6875*coeff[1]*fSkin[4]+19.6875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(70.53065765632414*coeff[1]*fSkin[20])+51.468317749209504*coeff[1]*fEdge[20]-61.5234375*coeff[1]*fSkin[5]-56.6015625*coeff[1]*fEdge[5]-34.09975027401226*coeff[1]*fSkin[1]+34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = -(35.21807064562169*coeff[1]*fSkin[33])+35.21807064562169*coeff[1]*fEdge[33]-34.09975027401226*coeff[1]*fSkin[15]-34.09975027401226*coeff[1]*fEdge[15]-19.6875*coeff[1]*fSkin[6]+19.6875*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = -(70.53065765632414*coeff[1]*fSkin[22])+51.468317749209504*coeff[1]*fEdge[22]-61.5234375*coeff[1]*fSkin[7]-56.6015625*coeff[1]*fEdge[7]-34.09975027401226*coeff[1]*fSkin[3]+34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = -(35.21807064562169*coeff[1]*fSkin[36])+35.21807064562169*coeff[1]*fEdge[36]-34.09975027401226*coeff[1]*fSkin[16]-34.09975027401226*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[8]+19.6875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(70.53065765632414*coeff[1]*fSkin[26])+51.468317749209504*coeff[1]*fEdge[26]-61.5234375*coeff[1]*fSkin[9]-56.6015625*coeff[1]*fEdge[9]-34.09975027401226*coeff[1]*fSkin[4]+34.09975027401226*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[1]*fSkin[38])+35.21807064562169*coeff[1]*fEdge[38]-34.09975027401226*coeff[1]*fSkin[18]-34.09975027401226*coeff[1]*fEdge[18]-19.6875*coeff[1]*fSkin[10]+19.6875*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = -(34.09975027401227*coeff[1]*fSkin[19])-34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[11]+19.6875*coeff[1]*fEdge[11]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]-31.316701275974005*coeff[1]*fSkin[2]-12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = -(34.09975027401227*coeff[1]*fSkin[24])-34.09975027401227*coeff[1]*fEdge[24]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = -(34.09975027401227*coeff[1]*fSkin[29])-34.09975027401227*coeff[1]*fEdge[29]-19.6875*coeff[1]*fSkin[14]+19.6875*coeff[1]*fEdge[14]; 
  edgeSurf_incr[15] = -(70.53065765632411*coeff[1]*fSkin[33])+51.46831774920949*coeff[1]*fEdge[33]-61.5234375*coeff[1]*fSkin[15]-56.6015625*coeff[1]*fEdge[15]-34.09975027401226*coeff[1]*fSkin[6]+34.09975027401226*coeff[1]*fEdge[6]; 
  edgeSurf_incr[16] = -(70.53065765632411*coeff[1]*fSkin[36])+51.46831774920949*coeff[1]*fEdge[36]-61.5234375*coeff[1]*fSkin[16]-56.6015625*coeff[1]*fEdge[16]-34.09975027401226*coeff[1]*fSkin[8]+34.09975027401226*coeff[1]*fEdge[8]; 
  edgeSurf_incr[17] = -(35.21807064562168*coeff[1]*fSkin[45])+35.21807064562168*coeff[1]*fEdge[45]-34.09975027401226*coeff[1]*fSkin[31]-34.09975027401226*coeff[1]*fEdge[31]-19.6875*coeff[1]*fSkin[17]+19.6875*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = -(70.53065765632411*coeff[1]*fSkin[38])+51.46831774920949*coeff[1]*fEdge[38]-61.5234375*coeff[1]*fSkin[18]-56.6015625*coeff[1]*fEdge[18]-34.09975027401226*coeff[1]*fSkin[10]+34.09975027401226*coeff[1]*fEdge[10]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*fSkin[19])-56.6015625*coeff[1]*fEdge[19]-34.09975027401227*coeff[1]*fSkin[11]+34.09975027401227*coeff[1]*fEdge[11]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*fSkin[20])-3.1640625*coeff[1]*fEdge[20]-31.316701275974033*coeff[1]*fSkin[5]-12.2543613688594*coeff[1]*fEdge[5]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[21] = -(34.09975027401227*coeff[1]*fSkin[32])-34.09975027401227*coeff[1]*fEdge[32]-19.6875*coeff[1]*fSkin[21]+19.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*fSkin[22])-3.1640625*coeff[1]*fEdge[22]-31.316701275974033*coeff[1]*fSkin[7]-12.2543613688594*coeff[1]*fEdge[7]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[23] = -(34.09975027401227*coeff[1]*fSkin[34])-34.09975027401227*coeff[1]*fEdge[34]-19.6875*coeff[1]*fSkin[23]+19.6875*coeff[1]*fEdge[23]; 
  edgeSurf_incr[24] = -(61.5234375*coeff[1]*fSkin[24])-56.6015625*coeff[1]*fEdge[24]-34.09975027401227*coeff[1]*fSkin[13]+34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[25] = -(34.09975027401227*coeff[1]*fSkin[35])-34.09975027401227*coeff[1]*fEdge[35]-19.6875*coeff[1]*fSkin[25]+19.6875*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*fSkin[26])-3.1640625*coeff[1]*fEdge[26]-31.316701275974033*coeff[1]*fSkin[9]-12.2543613688594*coeff[1]*fEdge[9]-12.577882373436315*coeff[1]*fSkin[4]+12.577882373436315*coeff[1]*fEdge[4]; 
  edgeSurf_incr[27] = -(34.09975027401227*coeff[1]*fSkin[40])-34.09975027401227*coeff[1]*fEdge[40]-19.6875*coeff[1]*fSkin[27]+19.6875*coeff[1]*fEdge[27]; 
  edgeSurf_incr[28] = -(34.09975027401227*coeff[1]*fSkin[41])-34.09975027401227*coeff[1]*fEdge[41]-19.6875*coeff[1]*fSkin[28]+19.6875*coeff[1]*fEdge[28]; 
  edgeSurf_incr[29] = -(61.5234375*coeff[1]*fSkin[29])-56.6015625*coeff[1]*fEdge[29]-34.09975027401227*coeff[1]*fSkin[14]+34.09975027401227*coeff[1]*fEdge[14]; 
  edgeSurf_incr[30] = -(34.09975027401227*coeff[1]*fSkin[43])-34.09975027401227*coeff[1]*fEdge[43]-19.6875*coeff[1]*fSkin[30]+19.6875*coeff[1]*fEdge[30]; 
  edgeSurf_incr[31] = -(70.53065765632414*coeff[1]*fSkin[45])+51.468317749209504*coeff[1]*fEdge[45]-61.5234375*coeff[1]*fSkin[31]-56.6015625*coeff[1]*fEdge[31]-34.09975027401226*coeff[1]*fSkin[17]+34.09975027401226*coeff[1]*fEdge[17]; 
  edgeSurf_incr[32] = -(61.5234375*coeff[1]*fSkin[32])-56.6015625*coeff[1]*fEdge[32]-34.09975027401227*coeff[1]*fSkin[21]+34.09975027401227*coeff[1]*fEdge[21]; 
  edgeSurf_incr[33] = -(70.6640625*coeff[1]*fSkin[33])-3.1640625*coeff[1]*fEdge[33]-31.316701275974005*coeff[1]*fSkin[15]-12.2543613688594*coeff[1]*fEdge[15]-12.57788237343632*coeff[1]*fSkin[6]+12.57788237343632*coeff[1]*fEdge[6]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[1]*fSkin[34])-56.6015625*coeff[1]*fEdge[34]-34.09975027401227*coeff[1]*fSkin[23]+34.09975027401227*coeff[1]*fEdge[23]; 
  edgeSurf_incr[35] = -(61.5234375*coeff[1]*fSkin[35])-56.6015625*coeff[1]*fEdge[35]-34.09975027401227*coeff[1]*fSkin[25]+34.09975027401227*coeff[1]*fEdge[25]; 
  edgeSurf_incr[36] = -(70.6640625*coeff[1]*fSkin[36])-3.1640625*coeff[1]*fEdge[36]-31.316701275974005*coeff[1]*fSkin[16]-12.2543613688594*coeff[1]*fEdge[16]-12.57788237343632*coeff[1]*fSkin[8]+12.57788237343632*coeff[1]*fEdge[8]; 
  edgeSurf_incr[37] = -(34.09975027401227*coeff[1]*fSkin[44])-34.09975027401227*coeff[1]*fEdge[44]-19.6875*coeff[1]*fSkin[37]+19.6875*coeff[1]*fEdge[37]; 
  edgeSurf_incr[38] = -(70.6640625*coeff[1]*fSkin[38])-3.1640625*coeff[1]*fEdge[38]-31.316701275974005*coeff[1]*fSkin[18]-12.2543613688594*coeff[1]*fEdge[18]-12.57788237343632*coeff[1]*fSkin[10]+12.57788237343632*coeff[1]*fEdge[10]; 
  edgeSurf_incr[39] = -(34.09975027401227*coeff[1]*fSkin[46])-34.09975027401227*coeff[1]*fEdge[46]-19.6875*coeff[1]*fSkin[39]+19.6875*coeff[1]*fEdge[39]; 
  edgeSurf_incr[40] = -(61.5234375*coeff[1]*fSkin[40])-56.6015625*coeff[1]*fEdge[40]-34.09975027401227*coeff[1]*fSkin[27]+34.09975027401227*coeff[1]*fEdge[27]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[1]*fSkin[41])-56.6015625*coeff[1]*fEdge[41]-34.09975027401227*coeff[1]*fSkin[28]+34.09975027401227*coeff[1]*fEdge[28]; 
  edgeSurf_incr[42] = -(34.09975027401227*coeff[1]*fSkin[47])-34.09975027401227*coeff[1]*fEdge[47]-19.6875*coeff[1]*fSkin[42]+19.6875*coeff[1]*fEdge[42]; 
  edgeSurf_incr[43] = -(61.5234375*coeff[1]*fSkin[43])-56.6015625*coeff[1]*fEdge[43]-34.09975027401227*coeff[1]*fSkin[30]+34.09975027401227*coeff[1]*fEdge[30]; 
  edgeSurf_incr[44] = -(61.5234375*coeff[1]*fSkin[44])-56.6015625*coeff[1]*fEdge[44]-34.09975027401227*coeff[1]*fSkin[37]+34.09975027401227*coeff[1]*fEdge[37]; 
  edgeSurf_incr[45] = -(70.6640625*coeff[1]*fSkin[45])-3.1640625*coeff[1]*fEdge[45]-31.316701275974033*coeff[1]*fSkin[31]-12.2543613688594*coeff[1]*fEdge[31]-12.577882373436315*coeff[1]*fSkin[17]+12.577882373436315*coeff[1]*fEdge[17]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[1]*fSkin[46])-56.6015625*coeff[1]*fEdge[46]-34.09975027401227*coeff[1]*fSkin[39]+34.09975027401227*coeff[1]*fEdge[39]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[1]*fSkin[47])-56.6015625*coeff[1]*fEdge[47]-34.09975027401227*coeff[1]*fSkin[42]+34.09975027401227*coeff[1]*fEdge[42]; 

  boundSurf_incr[2] = 19.062339907114627*coeff[1]*fSkin[12]-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = 19.062339907114634*coeff[1]*fSkin[20]-4.921875*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = 19.062339907114634*coeff[1]*fSkin[22]-4.921875*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = 19.062339907114634*coeff[1]*fSkin[26]-4.921875*coeff[1]*fSkin[9]; 
  boundSurf_incr[12] = 19.062339907114627*coeff[1]*fSkin[2]-73.828125*coeff[1]*fSkin[12]; 
  boundSurf_incr[15] = 19.062339907114627*coeff[1]*fSkin[33]-4.921875*coeff[1]*fSkin[15]; 
  boundSurf_incr[16] = 19.062339907114627*coeff[1]*fSkin[36]-4.921875*coeff[1]*fSkin[16]; 
  boundSurf_incr[18] = 19.062339907114627*coeff[1]*fSkin[38]-4.921875*coeff[1]*fSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*fSkin[19]); 
  boundSurf_incr[20] = 19.062339907114634*coeff[1]*fSkin[5]-73.828125*coeff[1]*fSkin[20]; 
  boundSurf_incr[22] = 19.062339907114634*coeff[1]*fSkin[7]-73.828125*coeff[1]*fSkin[22]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*fSkin[24]); 
  boundSurf_incr[26] = 19.062339907114634*coeff[1]*fSkin[9]-73.828125*coeff[1]*fSkin[26]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*fSkin[29]); 
  boundSurf_incr[31] = 19.062339907114634*coeff[1]*fSkin[45]-4.921875*coeff[1]*fSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*fSkin[32]); 
  boundSurf_incr[33] = 19.062339907114627*coeff[1]*fSkin[15]-73.828125*coeff[1]*fSkin[33]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*fSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*fSkin[35]); 
  boundSurf_incr[36] = 19.062339907114627*coeff[1]*fSkin[16]-73.828125*coeff[1]*fSkin[36]; 
  boundSurf_incr[38] = 19.062339907114627*coeff[1]*fSkin[18]-73.828125*coeff[1]*fSkin[38]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*fSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*fSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*fSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*fSkin[44]); 
  boundSurf_incr[45] = 19.062339907114634*coeff[1]*fSkin[31]-73.828125*coeff[1]*fSkin[45]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*fSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*fSkin[47]); 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*fSkin[12])+35.21807064562169*coeff[1]*fEdge[12]+34.09975027401226*coeff[1]*fSkin[2]+34.09975027401226*coeff[1]*fEdge[2]-19.6875*fSkin[0]*coeff[1]+19.6875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*fSkin[20])+35.21807064562168*coeff[1]*fEdge[20]+34.09975027401226*coeff[1]*fSkin[5]+34.09975027401226*coeff[1]*fEdge[5]-19.6875*coeff[1]*fSkin[1]+19.6875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 70.53065765632411*coeff[1]*fSkin[12]-51.46831774920949*coeff[1]*fEdge[12]-61.5234375*coeff[1]*fSkin[2]-56.6015625*coeff[1]*fEdge[2]+34.09975027401226*fSkin[0]*coeff[1]-34.09975027401226*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*fSkin[22])+35.21807064562168*coeff[1]*fEdge[22]+34.09975027401226*coeff[1]*fSkin[7]+34.09975027401226*coeff[1]*fEdge[7]-19.6875*coeff[1]*fSkin[3]+19.6875*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[1]*fSkin[26])+35.21807064562168*coeff[1]*fEdge[26]+34.09975027401226*coeff[1]*fSkin[9]+34.09975027401226*coeff[1]*fEdge[9]-19.6875*coeff[1]*fSkin[4]+19.6875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[1]*fSkin[20]-51.468317749209504*coeff[1]*fEdge[20]-61.5234375*coeff[1]*fSkin[5]-56.6015625*coeff[1]*fEdge[5]+34.09975027401226*coeff[1]*fSkin[1]-34.09975027401226*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = -(35.21807064562169*coeff[1]*fSkin[33])+35.21807064562169*coeff[1]*fEdge[33]+34.09975027401226*coeff[1]*fSkin[15]+34.09975027401226*coeff[1]*fEdge[15]-19.6875*coeff[1]*fSkin[6]+19.6875*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = 70.53065765632414*coeff[1]*fSkin[22]-51.468317749209504*coeff[1]*fEdge[22]-61.5234375*coeff[1]*fSkin[7]-56.6015625*coeff[1]*fEdge[7]+34.09975027401226*coeff[1]*fSkin[3]-34.09975027401226*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = -(35.21807064562169*coeff[1]*fSkin[36])+35.21807064562169*coeff[1]*fEdge[36]+34.09975027401226*coeff[1]*fSkin[16]+34.09975027401226*coeff[1]*fEdge[16]-19.6875*coeff[1]*fSkin[8]+19.6875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = 70.53065765632414*coeff[1]*fSkin[26]-51.468317749209504*coeff[1]*fEdge[26]-61.5234375*coeff[1]*fSkin[9]-56.6015625*coeff[1]*fEdge[9]+34.09975027401226*coeff[1]*fSkin[4]-34.09975027401226*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[1]*fSkin[38])+35.21807064562169*coeff[1]*fEdge[38]+34.09975027401226*coeff[1]*fSkin[18]+34.09975027401226*coeff[1]*fEdge[18]-19.6875*coeff[1]*fSkin[10]+19.6875*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = 34.09975027401227*coeff[1]*fSkin[19]+34.09975027401227*coeff[1]*fEdge[19]-19.6875*coeff[1]*fSkin[11]+19.6875*coeff[1]*fEdge[11]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*fSkin[12])-3.1640625*coeff[1]*fEdge[12]+31.316701275974005*coeff[1]*fSkin[2]+12.2543613688594*coeff[1]*fEdge[2]-12.57788237343632*fSkin[0]*coeff[1]+12.57788237343632*fEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = 34.09975027401227*coeff[1]*fSkin[24]+34.09975027401227*coeff[1]*fEdge[24]-19.6875*coeff[1]*fSkin[13]+19.6875*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[1]*fSkin[29]+34.09975027401227*coeff[1]*fEdge[29]-19.6875*coeff[1]*fSkin[14]+19.6875*coeff[1]*fEdge[14]; 
  edgeSurf_incr[15] = 70.53065765632411*coeff[1]*fSkin[33]-51.46831774920949*coeff[1]*fEdge[33]-61.5234375*coeff[1]*fSkin[15]-56.6015625*coeff[1]*fEdge[15]+34.09975027401226*coeff[1]*fSkin[6]-34.09975027401226*coeff[1]*fEdge[6]; 
  edgeSurf_incr[16] = 70.53065765632411*coeff[1]*fSkin[36]-51.46831774920949*coeff[1]*fEdge[36]-61.5234375*coeff[1]*fSkin[16]-56.6015625*coeff[1]*fEdge[16]+34.09975027401226*coeff[1]*fSkin[8]-34.09975027401226*coeff[1]*fEdge[8]; 
  edgeSurf_incr[17] = -(35.21807064562168*coeff[1]*fSkin[45])+35.21807064562168*coeff[1]*fEdge[45]+34.09975027401226*coeff[1]*fSkin[31]+34.09975027401226*coeff[1]*fEdge[31]-19.6875*coeff[1]*fSkin[17]+19.6875*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = 70.53065765632411*coeff[1]*fSkin[38]-51.46831774920949*coeff[1]*fEdge[38]-61.5234375*coeff[1]*fSkin[18]-56.6015625*coeff[1]*fEdge[18]+34.09975027401226*coeff[1]*fSkin[10]-34.09975027401226*coeff[1]*fEdge[10]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*fSkin[19])-56.6015625*coeff[1]*fEdge[19]+34.09975027401227*coeff[1]*fSkin[11]-34.09975027401227*coeff[1]*fEdge[11]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*fSkin[20])-3.1640625*coeff[1]*fEdge[20]+31.316701275974033*coeff[1]*fSkin[5]+12.2543613688594*coeff[1]*fEdge[5]-12.577882373436315*coeff[1]*fSkin[1]+12.577882373436315*coeff[1]*fEdge[1]; 
  edgeSurf_incr[21] = 34.09975027401227*coeff[1]*fSkin[32]+34.09975027401227*coeff[1]*fEdge[32]-19.6875*coeff[1]*fSkin[21]+19.6875*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*fSkin[22])-3.1640625*coeff[1]*fEdge[22]+31.316701275974033*coeff[1]*fSkin[7]+12.2543613688594*coeff[1]*fEdge[7]-12.577882373436315*coeff[1]*fSkin[3]+12.577882373436315*coeff[1]*fEdge[3]; 
  edgeSurf_incr[23] = 34.09975027401227*coeff[1]*fSkin[34]+34.09975027401227*coeff[1]*fEdge[34]-19.6875*coeff[1]*fSkin[23]+19.6875*coeff[1]*fEdge[23]; 
  edgeSurf_incr[24] = -(61.5234375*coeff[1]*fSkin[24])-56.6015625*coeff[1]*fEdge[24]+34.09975027401227*coeff[1]*fSkin[13]-34.09975027401227*coeff[1]*fEdge[13]; 
  edgeSurf_incr[25] = 34.09975027401227*coeff[1]*fSkin[35]+34.09975027401227*coeff[1]*fEdge[35]-19.6875*coeff[1]*fSkin[25]+19.6875*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*fSkin[26])-3.1640625*coeff[1]*fEdge[26]+31.316701275974033*coeff[1]*fSkin[9]+12.2543613688594*coeff[1]*fEdge[9]-12.577882373436315*coeff[1]*fSkin[4]+12.577882373436315*coeff[1]*fEdge[4]; 
  edgeSurf_incr[27] = 34.09975027401227*coeff[1]*fSkin[40]+34.09975027401227*coeff[1]*fEdge[40]-19.6875*coeff[1]*fSkin[27]+19.6875*coeff[1]*fEdge[27]; 
  edgeSurf_incr[28] = 34.09975027401227*coeff[1]*fSkin[41]+34.09975027401227*coeff[1]*fEdge[41]-19.6875*coeff[1]*fSkin[28]+19.6875*coeff[1]*fEdge[28]; 
  edgeSurf_incr[29] = -(61.5234375*coeff[1]*fSkin[29])-56.6015625*coeff[1]*fEdge[29]+34.09975027401227*coeff[1]*fSkin[14]-34.09975027401227*coeff[1]*fEdge[14]; 
  edgeSurf_incr[30] = 34.09975027401227*coeff[1]*fSkin[43]+34.09975027401227*coeff[1]*fEdge[43]-19.6875*coeff[1]*fSkin[30]+19.6875*coeff[1]*fEdge[30]; 
  edgeSurf_incr[31] = 70.53065765632414*coeff[1]*fSkin[45]-51.468317749209504*coeff[1]*fEdge[45]-61.5234375*coeff[1]*fSkin[31]-56.6015625*coeff[1]*fEdge[31]+34.09975027401226*coeff[1]*fSkin[17]-34.09975027401226*coeff[1]*fEdge[17]; 
  edgeSurf_incr[32] = -(61.5234375*coeff[1]*fSkin[32])-56.6015625*coeff[1]*fEdge[32]+34.09975027401227*coeff[1]*fSkin[21]-34.09975027401227*coeff[1]*fEdge[21]; 
  edgeSurf_incr[33] = -(70.6640625*coeff[1]*fSkin[33])-3.1640625*coeff[1]*fEdge[33]+31.316701275974005*coeff[1]*fSkin[15]+12.2543613688594*coeff[1]*fEdge[15]-12.57788237343632*coeff[1]*fSkin[6]+12.57788237343632*coeff[1]*fEdge[6]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[1]*fSkin[34])-56.6015625*coeff[1]*fEdge[34]+34.09975027401227*coeff[1]*fSkin[23]-34.09975027401227*coeff[1]*fEdge[23]; 
  edgeSurf_incr[35] = -(61.5234375*coeff[1]*fSkin[35])-56.6015625*coeff[1]*fEdge[35]+34.09975027401227*coeff[1]*fSkin[25]-34.09975027401227*coeff[1]*fEdge[25]; 
  edgeSurf_incr[36] = -(70.6640625*coeff[1]*fSkin[36])-3.1640625*coeff[1]*fEdge[36]+31.316701275974005*coeff[1]*fSkin[16]+12.2543613688594*coeff[1]*fEdge[16]-12.57788237343632*coeff[1]*fSkin[8]+12.57788237343632*coeff[1]*fEdge[8]; 
  edgeSurf_incr[37] = 34.09975027401227*coeff[1]*fSkin[44]+34.09975027401227*coeff[1]*fEdge[44]-19.6875*coeff[1]*fSkin[37]+19.6875*coeff[1]*fEdge[37]; 
  edgeSurf_incr[38] = -(70.6640625*coeff[1]*fSkin[38])-3.1640625*coeff[1]*fEdge[38]+31.316701275974005*coeff[1]*fSkin[18]+12.2543613688594*coeff[1]*fEdge[18]-12.57788237343632*coeff[1]*fSkin[10]+12.57788237343632*coeff[1]*fEdge[10]; 
  edgeSurf_incr[39] = 34.09975027401227*coeff[1]*fSkin[46]+34.09975027401227*coeff[1]*fEdge[46]-19.6875*coeff[1]*fSkin[39]+19.6875*coeff[1]*fEdge[39]; 
  edgeSurf_incr[40] = -(61.5234375*coeff[1]*fSkin[40])-56.6015625*coeff[1]*fEdge[40]+34.09975027401227*coeff[1]*fSkin[27]-34.09975027401227*coeff[1]*fEdge[27]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[1]*fSkin[41])-56.6015625*coeff[1]*fEdge[41]+34.09975027401227*coeff[1]*fSkin[28]-34.09975027401227*coeff[1]*fEdge[28]; 
  edgeSurf_incr[42] = 34.09975027401227*coeff[1]*fSkin[47]+34.09975027401227*coeff[1]*fEdge[47]-19.6875*coeff[1]*fSkin[42]+19.6875*coeff[1]*fEdge[42]; 
  edgeSurf_incr[43] = -(61.5234375*coeff[1]*fSkin[43])-56.6015625*coeff[1]*fEdge[43]+34.09975027401227*coeff[1]*fSkin[30]-34.09975027401227*coeff[1]*fEdge[30]; 
  edgeSurf_incr[44] = -(61.5234375*coeff[1]*fSkin[44])-56.6015625*coeff[1]*fEdge[44]+34.09975027401227*coeff[1]*fSkin[37]-34.09975027401227*coeff[1]*fEdge[37]; 
  edgeSurf_incr[45] = -(70.6640625*coeff[1]*fSkin[45])-3.1640625*coeff[1]*fEdge[45]+31.316701275974033*coeff[1]*fSkin[31]+12.2543613688594*coeff[1]*fEdge[31]-12.577882373436315*coeff[1]*fSkin[17]+12.577882373436315*coeff[1]*fEdge[17]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[1]*fSkin[46])-56.6015625*coeff[1]*fEdge[46]+34.09975027401227*coeff[1]*fSkin[39]-34.09975027401227*coeff[1]*fEdge[39]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[1]*fSkin[47])-56.6015625*coeff[1]*fEdge[47]+34.09975027401227*coeff[1]*fSkin[42]-34.09975027401227*coeff[1]*fEdge[42]; 

  boundSurf_incr[2] = -(19.062339907114627*coeff[1]*fSkin[12])-4.921875*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = -(19.062339907114634*coeff[1]*fSkin[20])-4.921875*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = -(19.062339907114634*coeff[1]*fSkin[22])-4.921875*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = -(19.062339907114634*coeff[1]*fSkin[26])-4.921875*coeff[1]*fSkin[9]; 
  boundSurf_incr[12] = -(73.828125*coeff[1]*fSkin[12])-19.062339907114627*coeff[1]*fSkin[2]; 
  boundSurf_incr[15] = -(19.062339907114627*coeff[1]*fSkin[33])-4.921875*coeff[1]*fSkin[15]; 
  boundSurf_incr[16] = -(19.062339907114627*coeff[1]*fSkin[36])-4.921875*coeff[1]*fSkin[16]; 
  boundSurf_incr[18] = -(19.062339907114627*coeff[1]*fSkin[38])-4.921875*coeff[1]*fSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*fSkin[19]); 
  boundSurf_incr[20] = -(73.828125*coeff[1]*fSkin[20])-19.062339907114634*coeff[1]*fSkin[5]; 
  boundSurf_incr[22] = -(73.828125*coeff[1]*fSkin[22])-19.062339907114634*coeff[1]*fSkin[7]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*fSkin[24]); 
  boundSurf_incr[26] = -(73.828125*coeff[1]*fSkin[26])-19.062339907114634*coeff[1]*fSkin[9]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*fSkin[29]); 
  boundSurf_incr[31] = -(19.062339907114634*coeff[1]*fSkin[45])-4.921875*coeff[1]*fSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*fSkin[32]); 
  boundSurf_incr[33] = -(73.828125*coeff[1]*fSkin[33])-19.062339907114627*coeff[1]*fSkin[15]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*fSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*fSkin[35]); 
  boundSurf_incr[36] = -(73.828125*coeff[1]*fSkin[36])-19.062339907114627*coeff[1]*fSkin[16]; 
  boundSurf_incr[38] = -(73.828125*coeff[1]*fSkin[38])-19.062339907114627*coeff[1]*fSkin[18]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*fSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*fSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*fSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*fSkin[44]); 
  boundSurf_incr[45] = -(73.828125*coeff[1]*fSkin[45])-19.062339907114634*coeff[1]*fSkin[31]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*fSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*fSkin[47]); 

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

