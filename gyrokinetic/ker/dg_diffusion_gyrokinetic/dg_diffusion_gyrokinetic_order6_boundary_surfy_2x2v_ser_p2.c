#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_surfy_2x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[1],6.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*qSkin[12])+35.21807064562169*coeff[1]*qEdge[12]-34.09975027401226*coeff[1]*qSkin[2]-34.09975027401226*coeff[1]*qEdge[2]-19.6875*qSkin[0]*coeff[1]+19.6875*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*qSkin[20])+35.21807064562168*coeff[1]*qEdge[20]-34.09975027401226*coeff[1]*qSkin[5]-34.09975027401226*coeff[1]*qEdge[5]-19.6875*coeff[1]*qSkin[1]+19.6875*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = -(70.53065765632411*coeff[1]*qSkin[12])+51.46831774920949*coeff[1]*qEdge[12]-61.5234375*coeff[1]*qSkin[2]-56.6015625*coeff[1]*qEdge[2]-34.09975027401226*qSkin[0]*coeff[1]+34.09975027401226*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*qSkin[22])+35.21807064562168*coeff[1]*qEdge[22]-34.09975027401226*coeff[1]*qSkin[7]-34.09975027401226*coeff[1]*qEdge[7]-19.6875*coeff[1]*qSkin[3]+19.6875*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[1]*qSkin[26])+35.21807064562168*coeff[1]*qEdge[26]-34.09975027401226*coeff[1]*qSkin[9]-34.09975027401226*coeff[1]*qEdge[9]-19.6875*coeff[1]*qSkin[4]+19.6875*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = -(70.53065765632414*coeff[1]*qSkin[20])+51.468317749209504*coeff[1]*qEdge[20]-61.5234375*coeff[1]*qSkin[5]-56.6015625*coeff[1]*qEdge[5]-34.09975027401226*coeff[1]*qSkin[1]+34.09975027401226*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = -(35.21807064562169*coeff[1]*qSkin[33])+35.21807064562169*coeff[1]*qEdge[33]-34.09975027401226*coeff[1]*qSkin[15]-34.09975027401226*coeff[1]*qEdge[15]-19.6875*coeff[1]*qSkin[6]+19.6875*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = -(70.53065765632414*coeff[1]*qSkin[22])+51.468317749209504*coeff[1]*qEdge[22]-61.5234375*coeff[1]*qSkin[7]-56.6015625*coeff[1]*qEdge[7]-34.09975027401226*coeff[1]*qSkin[3]+34.09975027401226*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = -(35.21807064562169*coeff[1]*qSkin[36])+35.21807064562169*coeff[1]*qEdge[36]-34.09975027401226*coeff[1]*qSkin[16]-34.09975027401226*coeff[1]*qEdge[16]-19.6875*coeff[1]*qSkin[8]+19.6875*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = -(70.53065765632414*coeff[1]*qSkin[26])+51.468317749209504*coeff[1]*qEdge[26]-61.5234375*coeff[1]*qSkin[9]-56.6015625*coeff[1]*qEdge[9]-34.09975027401226*coeff[1]*qSkin[4]+34.09975027401226*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[1]*qSkin[38])+35.21807064562169*coeff[1]*qEdge[38]-34.09975027401226*coeff[1]*qSkin[18]-34.09975027401226*coeff[1]*qEdge[18]-19.6875*coeff[1]*qSkin[10]+19.6875*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = -(34.09975027401227*coeff[1]*qSkin[19])-34.09975027401227*coeff[1]*qEdge[19]-19.6875*coeff[1]*qSkin[11]+19.6875*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*qSkin[12])-3.1640625*coeff[1]*qEdge[12]-31.316701275974005*coeff[1]*qSkin[2]-12.2543613688594*coeff[1]*qEdge[2]-12.57788237343632*qSkin[0]*coeff[1]+12.57788237343632*qEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = -(34.09975027401227*coeff[1]*qSkin[24])-34.09975027401227*coeff[1]*qEdge[24]-19.6875*coeff[1]*qSkin[13]+19.6875*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = -(34.09975027401227*coeff[1]*qSkin[29])-34.09975027401227*coeff[1]*qEdge[29]-19.6875*coeff[1]*qSkin[14]+19.6875*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = -(70.53065765632411*coeff[1]*qSkin[33])+51.46831774920949*coeff[1]*qEdge[33]-61.5234375*coeff[1]*qSkin[15]-56.6015625*coeff[1]*qEdge[15]-34.09975027401226*coeff[1]*qSkin[6]+34.09975027401226*coeff[1]*qEdge[6]; 
  edgeSurf_incr[16] = -(70.53065765632411*coeff[1]*qSkin[36])+51.46831774920949*coeff[1]*qEdge[36]-61.5234375*coeff[1]*qSkin[16]-56.6015625*coeff[1]*qEdge[16]-34.09975027401226*coeff[1]*qSkin[8]+34.09975027401226*coeff[1]*qEdge[8]; 
  edgeSurf_incr[17] = -(35.21807064562168*coeff[1]*qSkin[45])+35.21807064562168*coeff[1]*qEdge[45]-34.09975027401226*coeff[1]*qSkin[31]-34.09975027401226*coeff[1]*qEdge[31]-19.6875*coeff[1]*qSkin[17]+19.6875*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = -(70.53065765632411*coeff[1]*qSkin[38])+51.46831774920949*coeff[1]*qEdge[38]-61.5234375*coeff[1]*qSkin[18]-56.6015625*coeff[1]*qEdge[18]-34.09975027401226*coeff[1]*qSkin[10]+34.09975027401226*coeff[1]*qEdge[10]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*qSkin[19])-56.6015625*coeff[1]*qEdge[19]-34.09975027401227*coeff[1]*qSkin[11]+34.09975027401227*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*qSkin[20])-3.1640625*coeff[1]*qEdge[20]-31.316701275974033*coeff[1]*qSkin[5]-12.2543613688594*coeff[1]*qEdge[5]-12.577882373436315*coeff[1]*qSkin[1]+12.577882373436315*coeff[1]*qEdge[1]; 
  edgeSurf_incr[21] = -(34.09975027401227*coeff[1]*qSkin[32])-34.09975027401227*coeff[1]*qEdge[32]-19.6875*coeff[1]*qSkin[21]+19.6875*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*qSkin[22])-3.1640625*coeff[1]*qEdge[22]-31.316701275974033*coeff[1]*qSkin[7]-12.2543613688594*coeff[1]*qEdge[7]-12.577882373436315*coeff[1]*qSkin[3]+12.577882373436315*coeff[1]*qEdge[3]; 
  edgeSurf_incr[23] = -(34.09975027401227*coeff[1]*qSkin[34])-34.09975027401227*coeff[1]*qEdge[34]-19.6875*coeff[1]*qSkin[23]+19.6875*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = -(61.5234375*coeff[1]*qSkin[24])-56.6015625*coeff[1]*qEdge[24]-34.09975027401227*coeff[1]*qSkin[13]+34.09975027401227*coeff[1]*qEdge[13]; 
  edgeSurf_incr[25] = -(34.09975027401227*coeff[1]*qSkin[35])-34.09975027401227*coeff[1]*qEdge[35]-19.6875*coeff[1]*qSkin[25]+19.6875*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*qSkin[26])-3.1640625*coeff[1]*qEdge[26]-31.316701275974033*coeff[1]*qSkin[9]-12.2543613688594*coeff[1]*qEdge[9]-12.577882373436315*coeff[1]*qSkin[4]+12.577882373436315*coeff[1]*qEdge[4]; 
  edgeSurf_incr[27] = -(34.09975027401227*coeff[1]*qSkin[40])-34.09975027401227*coeff[1]*qEdge[40]-19.6875*coeff[1]*qSkin[27]+19.6875*coeff[1]*qEdge[27]; 
  edgeSurf_incr[28] = -(34.09975027401227*coeff[1]*qSkin[41])-34.09975027401227*coeff[1]*qEdge[41]-19.6875*coeff[1]*qSkin[28]+19.6875*coeff[1]*qEdge[28]; 
  edgeSurf_incr[29] = -(61.5234375*coeff[1]*qSkin[29])-56.6015625*coeff[1]*qEdge[29]-34.09975027401227*coeff[1]*qSkin[14]+34.09975027401227*coeff[1]*qEdge[14]; 
  edgeSurf_incr[30] = -(34.09975027401227*coeff[1]*qSkin[43])-34.09975027401227*coeff[1]*qEdge[43]-19.6875*coeff[1]*qSkin[30]+19.6875*coeff[1]*qEdge[30]; 
  edgeSurf_incr[31] = -(70.53065765632414*coeff[1]*qSkin[45])+51.468317749209504*coeff[1]*qEdge[45]-61.5234375*coeff[1]*qSkin[31]-56.6015625*coeff[1]*qEdge[31]-34.09975027401226*coeff[1]*qSkin[17]+34.09975027401226*coeff[1]*qEdge[17]; 
  edgeSurf_incr[32] = -(61.5234375*coeff[1]*qSkin[32])-56.6015625*coeff[1]*qEdge[32]-34.09975027401227*coeff[1]*qSkin[21]+34.09975027401227*coeff[1]*qEdge[21]; 
  edgeSurf_incr[33] = -(70.6640625*coeff[1]*qSkin[33])-3.1640625*coeff[1]*qEdge[33]-31.316701275974005*coeff[1]*qSkin[15]-12.2543613688594*coeff[1]*qEdge[15]-12.57788237343632*coeff[1]*qSkin[6]+12.57788237343632*coeff[1]*qEdge[6]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[1]*qSkin[34])-56.6015625*coeff[1]*qEdge[34]-34.09975027401227*coeff[1]*qSkin[23]+34.09975027401227*coeff[1]*qEdge[23]; 
  edgeSurf_incr[35] = -(61.5234375*coeff[1]*qSkin[35])-56.6015625*coeff[1]*qEdge[35]-34.09975027401227*coeff[1]*qSkin[25]+34.09975027401227*coeff[1]*qEdge[25]; 
  edgeSurf_incr[36] = -(70.6640625*coeff[1]*qSkin[36])-3.1640625*coeff[1]*qEdge[36]-31.316701275974005*coeff[1]*qSkin[16]-12.2543613688594*coeff[1]*qEdge[16]-12.57788237343632*coeff[1]*qSkin[8]+12.57788237343632*coeff[1]*qEdge[8]; 
  edgeSurf_incr[37] = -(34.09975027401227*coeff[1]*qSkin[44])-34.09975027401227*coeff[1]*qEdge[44]-19.6875*coeff[1]*qSkin[37]+19.6875*coeff[1]*qEdge[37]; 
  edgeSurf_incr[38] = -(70.6640625*coeff[1]*qSkin[38])-3.1640625*coeff[1]*qEdge[38]-31.316701275974005*coeff[1]*qSkin[18]-12.2543613688594*coeff[1]*qEdge[18]-12.57788237343632*coeff[1]*qSkin[10]+12.57788237343632*coeff[1]*qEdge[10]; 
  edgeSurf_incr[39] = -(34.09975027401227*coeff[1]*qSkin[46])-34.09975027401227*coeff[1]*qEdge[46]-19.6875*coeff[1]*qSkin[39]+19.6875*coeff[1]*qEdge[39]; 
  edgeSurf_incr[40] = -(61.5234375*coeff[1]*qSkin[40])-56.6015625*coeff[1]*qEdge[40]-34.09975027401227*coeff[1]*qSkin[27]+34.09975027401227*coeff[1]*qEdge[27]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[1]*qSkin[41])-56.6015625*coeff[1]*qEdge[41]-34.09975027401227*coeff[1]*qSkin[28]+34.09975027401227*coeff[1]*qEdge[28]; 
  edgeSurf_incr[42] = -(34.09975027401227*coeff[1]*qSkin[47])-34.09975027401227*coeff[1]*qEdge[47]-19.6875*coeff[1]*qSkin[42]+19.6875*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = -(61.5234375*coeff[1]*qSkin[43])-56.6015625*coeff[1]*qEdge[43]-34.09975027401227*coeff[1]*qSkin[30]+34.09975027401227*coeff[1]*qEdge[30]; 
  edgeSurf_incr[44] = -(61.5234375*coeff[1]*qSkin[44])-56.6015625*coeff[1]*qEdge[44]-34.09975027401227*coeff[1]*qSkin[37]+34.09975027401227*coeff[1]*qEdge[37]; 
  edgeSurf_incr[45] = -(70.6640625*coeff[1]*qSkin[45])-3.1640625*coeff[1]*qEdge[45]-31.316701275974033*coeff[1]*qSkin[31]-12.2543613688594*coeff[1]*qEdge[31]-12.577882373436315*coeff[1]*qSkin[17]+12.577882373436315*coeff[1]*qEdge[17]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[1]*qSkin[46])-56.6015625*coeff[1]*qEdge[46]-34.09975027401227*coeff[1]*qSkin[39]+34.09975027401227*coeff[1]*qEdge[39]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[1]*qSkin[47])-56.6015625*coeff[1]*qEdge[47]-34.09975027401227*coeff[1]*qSkin[42]+34.09975027401227*coeff[1]*qEdge[42]; 

  boundSurf_incr[2] = 19.062339907114627*coeff[1]*qSkin[12]-4.921875*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = 19.062339907114634*coeff[1]*qSkin[20]-4.921875*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = 19.062339907114634*coeff[1]*qSkin[22]-4.921875*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = 19.062339907114634*coeff[1]*qSkin[26]-4.921875*coeff[1]*qSkin[9]; 
  boundSurf_incr[12] = 19.062339907114627*coeff[1]*qSkin[2]-73.828125*coeff[1]*qSkin[12]; 
  boundSurf_incr[15] = 19.062339907114627*coeff[1]*qSkin[33]-4.921875*coeff[1]*qSkin[15]; 
  boundSurf_incr[16] = 19.062339907114627*coeff[1]*qSkin[36]-4.921875*coeff[1]*qSkin[16]; 
  boundSurf_incr[18] = 19.062339907114627*coeff[1]*qSkin[38]-4.921875*coeff[1]*qSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*qSkin[19]); 
  boundSurf_incr[20] = 19.062339907114634*coeff[1]*qSkin[5]-73.828125*coeff[1]*qSkin[20]; 
  boundSurf_incr[22] = 19.062339907114634*coeff[1]*qSkin[7]-73.828125*coeff[1]*qSkin[22]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*qSkin[24]); 
  boundSurf_incr[26] = 19.062339907114634*coeff[1]*qSkin[9]-73.828125*coeff[1]*qSkin[26]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*qSkin[29]); 
  boundSurf_incr[31] = 19.062339907114634*coeff[1]*qSkin[45]-4.921875*coeff[1]*qSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*qSkin[32]); 
  boundSurf_incr[33] = 19.062339907114627*coeff[1]*qSkin[15]-73.828125*coeff[1]*qSkin[33]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*qSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*qSkin[35]); 
  boundSurf_incr[36] = 19.062339907114627*coeff[1]*qSkin[16]-73.828125*coeff[1]*qSkin[36]; 
  boundSurf_incr[38] = 19.062339907114627*coeff[1]*qSkin[18]-73.828125*coeff[1]*qSkin[38]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*qSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*qSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*qSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*qSkin[44]); 
  boundSurf_incr[45] = 19.062339907114634*coeff[1]*qSkin[31]-73.828125*coeff[1]*qSkin[45]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*qSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*qSkin[47]); 

  } else { 

  edgeSurf_incr[0] = -(35.21807064562169*coeff[1]*qSkin[12])+35.21807064562169*coeff[1]*qEdge[12]+34.09975027401226*coeff[1]*qSkin[2]+34.09975027401226*coeff[1]*qEdge[2]-19.6875*qSkin[0]*coeff[1]+19.6875*qEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(35.21807064562168*coeff[1]*qSkin[20])+35.21807064562168*coeff[1]*qEdge[20]+34.09975027401226*coeff[1]*qSkin[5]+34.09975027401226*coeff[1]*qEdge[5]-19.6875*coeff[1]*qSkin[1]+19.6875*coeff[1]*qEdge[1]; 
  edgeSurf_incr[2] = 70.53065765632411*coeff[1]*qSkin[12]-51.46831774920949*coeff[1]*qEdge[12]-61.5234375*coeff[1]*qSkin[2]-56.6015625*coeff[1]*qEdge[2]+34.09975027401226*qSkin[0]*coeff[1]-34.09975027401226*qEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(35.21807064562168*coeff[1]*qSkin[22])+35.21807064562168*coeff[1]*qEdge[22]+34.09975027401226*coeff[1]*qSkin[7]+34.09975027401226*coeff[1]*qEdge[7]-19.6875*coeff[1]*qSkin[3]+19.6875*coeff[1]*qEdge[3]; 
  edgeSurf_incr[4] = -(35.21807064562168*coeff[1]*qSkin[26])+35.21807064562168*coeff[1]*qEdge[26]+34.09975027401226*coeff[1]*qSkin[9]+34.09975027401226*coeff[1]*qEdge[9]-19.6875*coeff[1]*qSkin[4]+19.6875*coeff[1]*qEdge[4]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[1]*qSkin[20]-51.468317749209504*coeff[1]*qEdge[20]-61.5234375*coeff[1]*qSkin[5]-56.6015625*coeff[1]*qEdge[5]+34.09975027401226*coeff[1]*qSkin[1]-34.09975027401226*coeff[1]*qEdge[1]; 
  edgeSurf_incr[6] = -(35.21807064562169*coeff[1]*qSkin[33])+35.21807064562169*coeff[1]*qEdge[33]+34.09975027401226*coeff[1]*qSkin[15]+34.09975027401226*coeff[1]*qEdge[15]-19.6875*coeff[1]*qSkin[6]+19.6875*coeff[1]*qEdge[6]; 
  edgeSurf_incr[7] = 70.53065765632414*coeff[1]*qSkin[22]-51.468317749209504*coeff[1]*qEdge[22]-61.5234375*coeff[1]*qSkin[7]-56.6015625*coeff[1]*qEdge[7]+34.09975027401226*coeff[1]*qSkin[3]-34.09975027401226*coeff[1]*qEdge[3]; 
  edgeSurf_incr[8] = -(35.21807064562169*coeff[1]*qSkin[36])+35.21807064562169*coeff[1]*qEdge[36]+34.09975027401226*coeff[1]*qSkin[16]+34.09975027401226*coeff[1]*qEdge[16]-19.6875*coeff[1]*qSkin[8]+19.6875*coeff[1]*qEdge[8]; 
  edgeSurf_incr[9] = 70.53065765632414*coeff[1]*qSkin[26]-51.468317749209504*coeff[1]*qEdge[26]-61.5234375*coeff[1]*qSkin[9]-56.6015625*coeff[1]*qEdge[9]+34.09975027401226*coeff[1]*qSkin[4]-34.09975027401226*coeff[1]*qEdge[4]; 
  edgeSurf_incr[10] = -(35.21807064562169*coeff[1]*qSkin[38])+35.21807064562169*coeff[1]*qEdge[38]+34.09975027401226*coeff[1]*qSkin[18]+34.09975027401226*coeff[1]*qEdge[18]-19.6875*coeff[1]*qSkin[10]+19.6875*coeff[1]*qEdge[10]; 
  edgeSurf_incr[11] = 34.09975027401227*coeff[1]*qSkin[19]+34.09975027401227*coeff[1]*qEdge[19]-19.6875*coeff[1]*qSkin[11]+19.6875*coeff[1]*qEdge[11]; 
  edgeSurf_incr[12] = -(70.6640625*coeff[1]*qSkin[12])-3.1640625*coeff[1]*qEdge[12]+31.316701275974005*coeff[1]*qSkin[2]+12.2543613688594*coeff[1]*qEdge[2]-12.57788237343632*qSkin[0]*coeff[1]+12.57788237343632*qEdge[0]*coeff[1]; 
  edgeSurf_incr[13] = 34.09975027401227*coeff[1]*qSkin[24]+34.09975027401227*coeff[1]*qEdge[24]-19.6875*coeff[1]*qSkin[13]+19.6875*coeff[1]*qEdge[13]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[1]*qSkin[29]+34.09975027401227*coeff[1]*qEdge[29]-19.6875*coeff[1]*qSkin[14]+19.6875*coeff[1]*qEdge[14]; 
  edgeSurf_incr[15] = 70.53065765632411*coeff[1]*qSkin[33]-51.46831774920949*coeff[1]*qEdge[33]-61.5234375*coeff[1]*qSkin[15]-56.6015625*coeff[1]*qEdge[15]+34.09975027401226*coeff[1]*qSkin[6]-34.09975027401226*coeff[1]*qEdge[6]; 
  edgeSurf_incr[16] = 70.53065765632411*coeff[1]*qSkin[36]-51.46831774920949*coeff[1]*qEdge[36]-61.5234375*coeff[1]*qSkin[16]-56.6015625*coeff[1]*qEdge[16]+34.09975027401226*coeff[1]*qSkin[8]-34.09975027401226*coeff[1]*qEdge[8]; 
  edgeSurf_incr[17] = -(35.21807064562168*coeff[1]*qSkin[45])+35.21807064562168*coeff[1]*qEdge[45]+34.09975027401226*coeff[1]*qSkin[31]+34.09975027401226*coeff[1]*qEdge[31]-19.6875*coeff[1]*qSkin[17]+19.6875*coeff[1]*qEdge[17]; 
  edgeSurf_incr[18] = 70.53065765632411*coeff[1]*qSkin[38]-51.46831774920949*coeff[1]*qEdge[38]-61.5234375*coeff[1]*qSkin[18]-56.6015625*coeff[1]*qEdge[18]+34.09975027401226*coeff[1]*qSkin[10]-34.09975027401226*coeff[1]*qEdge[10]; 
  edgeSurf_incr[19] = -(61.5234375*coeff[1]*qSkin[19])-56.6015625*coeff[1]*qEdge[19]+34.09975027401227*coeff[1]*qSkin[11]-34.09975027401227*coeff[1]*qEdge[11]; 
  edgeSurf_incr[20] = -(70.6640625*coeff[1]*qSkin[20])-3.1640625*coeff[1]*qEdge[20]+31.316701275974033*coeff[1]*qSkin[5]+12.2543613688594*coeff[1]*qEdge[5]-12.577882373436315*coeff[1]*qSkin[1]+12.577882373436315*coeff[1]*qEdge[1]; 
  edgeSurf_incr[21] = 34.09975027401227*coeff[1]*qSkin[32]+34.09975027401227*coeff[1]*qEdge[32]-19.6875*coeff[1]*qSkin[21]+19.6875*coeff[1]*qEdge[21]; 
  edgeSurf_incr[22] = -(70.6640625*coeff[1]*qSkin[22])-3.1640625*coeff[1]*qEdge[22]+31.316701275974033*coeff[1]*qSkin[7]+12.2543613688594*coeff[1]*qEdge[7]-12.577882373436315*coeff[1]*qSkin[3]+12.577882373436315*coeff[1]*qEdge[3]; 
  edgeSurf_incr[23] = 34.09975027401227*coeff[1]*qSkin[34]+34.09975027401227*coeff[1]*qEdge[34]-19.6875*coeff[1]*qSkin[23]+19.6875*coeff[1]*qEdge[23]; 
  edgeSurf_incr[24] = -(61.5234375*coeff[1]*qSkin[24])-56.6015625*coeff[1]*qEdge[24]+34.09975027401227*coeff[1]*qSkin[13]-34.09975027401227*coeff[1]*qEdge[13]; 
  edgeSurf_incr[25] = 34.09975027401227*coeff[1]*qSkin[35]+34.09975027401227*coeff[1]*qEdge[35]-19.6875*coeff[1]*qSkin[25]+19.6875*coeff[1]*qEdge[25]; 
  edgeSurf_incr[26] = -(70.6640625*coeff[1]*qSkin[26])-3.1640625*coeff[1]*qEdge[26]+31.316701275974033*coeff[1]*qSkin[9]+12.2543613688594*coeff[1]*qEdge[9]-12.577882373436315*coeff[1]*qSkin[4]+12.577882373436315*coeff[1]*qEdge[4]; 
  edgeSurf_incr[27] = 34.09975027401227*coeff[1]*qSkin[40]+34.09975027401227*coeff[1]*qEdge[40]-19.6875*coeff[1]*qSkin[27]+19.6875*coeff[1]*qEdge[27]; 
  edgeSurf_incr[28] = 34.09975027401227*coeff[1]*qSkin[41]+34.09975027401227*coeff[1]*qEdge[41]-19.6875*coeff[1]*qSkin[28]+19.6875*coeff[1]*qEdge[28]; 
  edgeSurf_incr[29] = -(61.5234375*coeff[1]*qSkin[29])-56.6015625*coeff[1]*qEdge[29]+34.09975027401227*coeff[1]*qSkin[14]-34.09975027401227*coeff[1]*qEdge[14]; 
  edgeSurf_incr[30] = 34.09975027401227*coeff[1]*qSkin[43]+34.09975027401227*coeff[1]*qEdge[43]-19.6875*coeff[1]*qSkin[30]+19.6875*coeff[1]*qEdge[30]; 
  edgeSurf_incr[31] = 70.53065765632414*coeff[1]*qSkin[45]-51.468317749209504*coeff[1]*qEdge[45]-61.5234375*coeff[1]*qSkin[31]-56.6015625*coeff[1]*qEdge[31]+34.09975027401226*coeff[1]*qSkin[17]-34.09975027401226*coeff[1]*qEdge[17]; 
  edgeSurf_incr[32] = -(61.5234375*coeff[1]*qSkin[32])-56.6015625*coeff[1]*qEdge[32]+34.09975027401227*coeff[1]*qSkin[21]-34.09975027401227*coeff[1]*qEdge[21]; 
  edgeSurf_incr[33] = -(70.6640625*coeff[1]*qSkin[33])-3.1640625*coeff[1]*qEdge[33]+31.316701275974005*coeff[1]*qSkin[15]+12.2543613688594*coeff[1]*qEdge[15]-12.57788237343632*coeff[1]*qSkin[6]+12.57788237343632*coeff[1]*qEdge[6]; 
  edgeSurf_incr[34] = -(61.5234375*coeff[1]*qSkin[34])-56.6015625*coeff[1]*qEdge[34]+34.09975027401227*coeff[1]*qSkin[23]-34.09975027401227*coeff[1]*qEdge[23]; 
  edgeSurf_incr[35] = -(61.5234375*coeff[1]*qSkin[35])-56.6015625*coeff[1]*qEdge[35]+34.09975027401227*coeff[1]*qSkin[25]-34.09975027401227*coeff[1]*qEdge[25]; 
  edgeSurf_incr[36] = -(70.6640625*coeff[1]*qSkin[36])-3.1640625*coeff[1]*qEdge[36]+31.316701275974005*coeff[1]*qSkin[16]+12.2543613688594*coeff[1]*qEdge[16]-12.57788237343632*coeff[1]*qSkin[8]+12.57788237343632*coeff[1]*qEdge[8]; 
  edgeSurf_incr[37] = 34.09975027401227*coeff[1]*qSkin[44]+34.09975027401227*coeff[1]*qEdge[44]-19.6875*coeff[1]*qSkin[37]+19.6875*coeff[1]*qEdge[37]; 
  edgeSurf_incr[38] = -(70.6640625*coeff[1]*qSkin[38])-3.1640625*coeff[1]*qEdge[38]+31.316701275974005*coeff[1]*qSkin[18]+12.2543613688594*coeff[1]*qEdge[18]-12.57788237343632*coeff[1]*qSkin[10]+12.57788237343632*coeff[1]*qEdge[10]; 
  edgeSurf_incr[39] = 34.09975027401227*coeff[1]*qSkin[46]+34.09975027401227*coeff[1]*qEdge[46]-19.6875*coeff[1]*qSkin[39]+19.6875*coeff[1]*qEdge[39]; 
  edgeSurf_incr[40] = -(61.5234375*coeff[1]*qSkin[40])-56.6015625*coeff[1]*qEdge[40]+34.09975027401227*coeff[1]*qSkin[27]-34.09975027401227*coeff[1]*qEdge[27]; 
  edgeSurf_incr[41] = -(61.5234375*coeff[1]*qSkin[41])-56.6015625*coeff[1]*qEdge[41]+34.09975027401227*coeff[1]*qSkin[28]-34.09975027401227*coeff[1]*qEdge[28]; 
  edgeSurf_incr[42] = 34.09975027401227*coeff[1]*qSkin[47]+34.09975027401227*coeff[1]*qEdge[47]-19.6875*coeff[1]*qSkin[42]+19.6875*coeff[1]*qEdge[42]; 
  edgeSurf_incr[43] = -(61.5234375*coeff[1]*qSkin[43])-56.6015625*coeff[1]*qEdge[43]+34.09975027401227*coeff[1]*qSkin[30]-34.09975027401227*coeff[1]*qEdge[30]; 
  edgeSurf_incr[44] = -(61.5234375*coeff[1]*qSkin[44])-56.6015625*coeff[1]*qEdge[44]+34.09975027401227*coeff[1]*qSkin[37]-34.09975027401227*coeff[1]*qEdge[37]; 
  edgeSurf_incr[45] = -(70.6640625*coeff[1]*qSkin[45])-3.1640625*coeff[1]*qEdge[45]+31.316701275974033*coeff[1]*qSkin[31]+12.2543613688594*coeff[1]*qEdge[31]-12.577882373436315*coeff[1]*qSkin[17]+12.577882373436315*coeff[1]*qEdge[17]; 
  edgeSurf_incr[46] = -(61.5234375*coeff[1]*qSkin[46])-56.6015625*coeff[1]*qEdge[46]+34.09975027401227*coeff[1]*qSkin[39]-34.09975027401227*coeff[1]*qEdge[39]; 
  edgeSurf_incr[47] = -(61.5234375*coeff[1]*qSkin[47])-56.6015625*coeff[1]*qEdge[47]+34.09975027401227*coeff[1]*qSkin[42]-34.09975027401227*coeff[1]*qEdge[42]; 

  boundSurf_incr[2] = -(19.062339907114627*coeff[1]*qSkin[12])-4.921875*coeff[1]*qSkin[2]; 
  boundSurf_incr[5] = -(19.062339907114634*coeff[1]*qSkin[20])-4.921875*coeff[1]*qSkin[5]; 
  boundSurf_incr[7] = -(19.062339907114634*coeff[1]*qSkin[22])-4.921875*coeff[1]*qSkin[7]; 
  boundSurf_incr[9] = -(19.062339907114634*coeff[1]*qSkin[26])-4.921875*coeff[1]*qSkin[9]; 
  boundSurf_incr[12] = -(73.828125*coeff[1]*qSkin[12])-19.062339907114627*coeff[1]*qSkin[2]; 
  boundSurf_incr[15] = -(19.062339907114627*coeff[1]*qSkin[33])-4.921875*coeff[1]*qSkin[15]; 
  boundSurf_incr[16] = -(19.062339907114627*coeff[1]*qSkin[36])-4.921875*coeff[1]*qSkin[16]; 
  boundSurf_incr[18] = -(19.062339907114627*coeff[1]*qSkin[38])-4.921875*coeff[1]*qSkin[18]; 
  boundSurf_incr[19] = -(4.921875*coeff[1]*qSkin[19]); 
  boundSurf_incr[20] = -(73.828125*coeff[1]*qSkin[20])-19.062339907114634*coeff[1]*qSkin[5]; 
  boundSurf_incr[22] = -(73.828125*coeff[1]*qSkin[22])-19.062339907114634*coeff[1]*qSkin[7]; 
  boundSurf_incr[24] = -(4.921875*coeff[1]*qSkin[24]); 
  boundSurf_incr[26] = -(73.828125*coeff[1]*qSkin[26])-19.062339907114634*coeff[1]*qSkin[9]; 
  boundSurf_incr[29] = -(4.921875*coeff[1]*qSkin[29]); 
  boundSurf_incr[31] = -(19.062339907114634*coeff[1]*qSkin[45])-4.921875*coeff[1]*qSkin[31]; 
  boundSurf_incr[32] = -(4.921875*coeff[1]*qSkin[32]); 
  boundSurf_incr[33] = -(73.828125*coeff[1]*qSkin[33])-19.062339907114627*coeff[1]*qSkin[15]; 
  boundSurf_incr[34] = -(4.921875*coeff[1]*qSkin[34]); 
  boundSurf_incr[35] = -(4.921875*coeff[1]*qSkin[35]); 
  boundSurf_incr[36] = -(73.828125*coeff[1]*qSkin[36])-19.062339907114627*coeff[1]*qSkin[16]; 
  boundSurf_incr[38] = -(73.828125*coeff[1]*qSkin[38])-19.062339907114627*coeff[1]*qSkin[18]; 
  boundSurf_incr[40] = -(4.921875*coeff[1]*qSkin[40]); 
  boundSurf_incr[41] = -(4.921875*coeff[1]*qSkin[41]); 
  boundSurf_incr[43] = -(4.921875*coeff[1]*qSkin[43]); 
  boundSurf_incr[44] = -(4.921875*coeff[1]*qSkin[44]); 
  boundSurf_incr[45] = -(73.828125*coeff[1]*qSkin[45])-19.062339907114634*coeff[1]*qSkin[31]; 
  boundSurf_incr[46] = -(4.921875*coeff[1]*qSkin[46]); 
  boundSurf_incr[47] = -(4.921875*coeff[1]*qSkin[47]); 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdx2Sq; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdx2Sq; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdx2Sq; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdx2Sq; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdx2Sq; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdx2Sq; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdx2Sq; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*rdx2Sq; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*rdx2Sq; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*rdx2Sq; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*rdx2Sq; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*rdx2Sq; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*rdx2Sq; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*rdx2Sq; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*rdx2Sq; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*rdx2Sq; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*rdx2Sq; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*rdx2Sq; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*rdx2Sq; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*rdx2Sq; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*rdx2Sq; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*rdx2Sq; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*rdx2Sq; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*rdx2Sq; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*rdx2Sq; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*rdx2Sq; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*rdx2Sq; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*rdx2Sq; 

  return 0.;
}

