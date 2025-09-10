#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *alpha_surf_edge, const double *alpha_surf_skin, 
    const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
    const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points.
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points.
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1.
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdz2 = 2.0/dxv[2];

  const double *alphaL = &alpha_surf_skin[48];
  const double *alphaR = &alpha_surf_edge[48];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[48];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[48];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[2];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[2];

  if (edge == -1) { 

  double fUpR[24] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fskin[3]+0.7071067811865475*fskin[0]; 
  fUpR[1] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[1]; 
  fUpR[2] = 1.224744871391589*fskin[8]+0.7071067811865475*fskin[2]; 
  fUpR[3] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[4]; 
  fUpR[4] = 1.224744871391589*fskin[14]+0.7071067811865475*fskin[5]; 
  fUpR[5] = 1.224744871391589*fskin[16]+0.7071067811865475*fskin[6]; 
  fUpR[6] = 1.224744871391589*fskin[18]+0.7071067811865475*fskin[9]; 
  fUpR[7] = 1.224744871391589*fskin[19]+0.7071067811865475*fskin[10]; 
  fUpR[8] = 1.224744871391589*fskin[21]+0.7071067811865475*fskin[12]; 
  fUpR[9] = 1.224744871391589*fskin[22]+0.7071067811865475*fskin[13]; 
  fUpR[10] = 1.224744871391589*fskin[25]+0.7071067811865475*fskin[15]; 
  fUpR[11] = 1.224744871391589*fskin[26]+0.7071067811865475*fskin[17]; 
  fUpR[12] = 1.224744871391589*fskin[27]+0.7071067811865475*fskin[20]; 
  fUpR[13] = 1.224744871391589*fskin[29]+0.7071067811865475*fskin[23]; 
  fUpR[14] = 1.224744871391589*fskin[30]+0.7071067811865475*fskin[24]; 
  fUpR[15] = 1.224744871391589*fskin[31]+0.7071067811865475*fskin[28]; 
  fUpR[16] = 1.224744871391589*fskin[35]+0.7071067811865475*fskin[32]; 
  fUpR[17] = 1.224744871391589*fskin[38]+0.7071067811865475*fskin[33]; 
  fUpR[18] = 1.224744871391589*fskin[39]+0.7071067811865475*fskin[34]; 
  fUpR[19] = 1.224744871391589*fskin[42]+0.7071067811865475*fskin[36]; 
  fUpR[20] = 1.224744871391589*fskin[43]+0.7071067811865475*fskin[37]; 
  fUpR[21] = 1.224744871391589*fskin[45]+0.7071067811865475*fskin[40]; 
  fUpR[22] = 1.224744871391589*fskin[46]+0.7071067811865475*fskin[41]; 
  fUpR[23] = 1.224744871391589*fskin[47]+0.7071067811865475*fskin[44]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[3]; 
  fUpR[1] = 0.7071067811865475*fedge[1]-1.224744871391589*fedge[7]; 
  fUpR[2] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[8]; 
  fUpR[3] = 0.7071067811865475*fedge[4]-1.224744871391589*fedge[11]; 
  fUpR[4] = 0.7071067811865475*fedge[5]-1.224744871391589*fedge[14]; 
  fUpR[5] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[16]; 
  fUpR[6] = 0.7071067811865475*fedge[9]-1.224744871391589*fedge[18]; 
  fUpR[7] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[19]; 
  fUpR[8] = 0.7071067811865475*fedge[12]-1.224744871391589*fedge[21]; 
  fUpR[9] = 0.7071067811865475*fedge[13]-1.224744871391589*fedge[22]; 
  fUpR[10] = 0.7071067811865475*fedge[15]-1.224744871391589*fedge[25]; 
  fUpR[11] = 0.7071067811865475*fedge[17]-1.224744871391589*fedge[26]; 
  fUpR[12] = 0.7071067811865475*fedge[20]-1.224744871391589*fedge[27]; 
  fUpR[13] = 0.7071067811865475*fedge[23]-1.224744871391589*fedge[29]; 
  fUpR[14] = 0.7071067811865475*fedge[24]-1.224744871391589*fedge[30]; 
  fUpR[15] = 0.7071067811865475*fedge[28]-1.224744871391589*fedge[31]; 
  fUpR[16] = 0.7071067811865475*fedge[32]-1.224744871391589*fedge[35]; 
  fUpR[17] = 0.7071067811865475*fedge[33]-1.224744871391589*fedge[38]; 
  fUpR[18] = 0.7071067811865475*fedge[34]-1.224744871391589*fedge[39]; 
  fUpR[19] = 0.7071067811865475*fedge[36]-1.224744871391589*fedge[42]; 
  fUpR[20] = 0.7071067811865475*fedge[37]-1.224744871391589*fedge[43]; 
  fUpR[21] = 0.7071067811865475*fedge[40]-1.224744871391589*fedge[45]; 
  fUpR[22] = 0.7071067811865475*fedge[41]-1.224744871391589*fedge[46]; 
  fUpR[23] = 0.7071067811865475*fedge[44]-1.224744871391589*fedge[47]; 
    } 
  } else { 
  double f_cr[24] = {0.};
  double f_rl[24] = {0.};
  double sgn_alphaUpR[6] = {0.};
  sgn_alphaUpR[0] = 0.2777777777777778*sgn_alpha_surfR[5]+0.4444444444444444*sgn_alpha_surfR[4]+0.2777777777777778*sgn_alpha_surfR[3]+0.2777777777777778*sgn_alpha_surfR[2]+0.4444444444444444*sgn_alpha_surfR[1]+0.2777777777777778*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[1] = 0.2777777777777778*sgn_alpha_surfR[5]+0.4444444444444444*sgn_alpha_surfR[4]+0.2777777777777778*sgn_alpha_surfR[3]-0.2777777777777778*sgn_alpha_surfR[2]-0.4444444444444444*sgn_alpha_surfR[1]-0.2777777777777778*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[2] = 0.3726779962499649*sgn_alpha_surfR[5]-0.3726779962499649*sgn_alpha_surfR[3]+0.3726779962499649*sgn_alpha_surfR[2]-0.3726779962499649*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[3] = 0.3726779962499649*sgn_alpha_surfR[5]-0.3726779962499649*sgn_alpha_surfR[3]-0.3726779962499649*sgn_alpha_surfR[2]+0.3726779962499649*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[4] = 0.2484519974999766*sgn_alpha_surfR[5]-0.4969039949999532*sgn_alpha_surfR[4]+0.2484519974999766*sgn_alpha_surfR[3]+0.2484519974999766*sgn_alpha_surfR[2]-0.4969039949999532*sgn_alpha_surfR[1]+0.2484519974999766*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[5] = 0.2484519974999767*sgn_alpha_surfR[5]-0.4969039949999535*sgn_alpha_surfR[4]+0.2484519974999767*sgn_alpha_surfR[3]-0.2484519974999767*sgn_alpha_surfR[2]+0.4969039949999535*sgn_alpha_surfR[1]-0.2484519974999767*sgn_alpha_surfR[0]; 

  f_cr[0] = 1.224744871391589*fskin[3]+0.7071067811865475*fskin[0]; 
  f_cr[1] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[1]; 
  f_cr[2] = 1.224744871391589*fskin[8]+0.7071067811865475*fskin[2]; 
  f_cr[3] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[4]; 
  f_cr[4] = 1.224744871391589*fskin[14]+0.7071067811865475*fskin[5]; 
  f_cr[5] = 1.224744871391589*fskin[16]+0.7071067811865475*fskin[6]; 
  f_cr[6] = 1.224744871391589*fskin[18]+0.7071067811865475*fskin[9]; 
  f_cr[7] = 1.224744871391589*fskin[19]+0.7071067811865475*fskin[10]; 
  f_cr[8] = 1.224744871391589*fskin[21]+0.7071067811865475*fskin[12]; 
  f_cr[9] = 1.224744871391589*fskin[22]+0.7071067811865475*fskin[13]; 
  f_cr[10] = 1.224744871391589*fskin[25]+0.7071067811865475*fskin[15]; 
  f_cr[11] = 1.224744871391589*fskin[26]+0.7071067811865475*fskin[17]; 
  f_cr[12] = 1.224744871391589*fskin[27]+0.7071067811865475*fskin[20]; 
  f_cr[13] = 1.224744871391589*fskin[29]+0.7071067811865475*fskin[23]; 
  f_cr[14] = 1.224744871391589*fskin[30]+0.7071067811865475*fskin[24]; 
  f_cr[15] = 1.224744871391589*fskin[31]+0.7071067811865475*fskin[28]; 
  f_cr[16] = 1.224744871391589*fskin[35]+0.7071067811865475*fskin[32]; 
  f_cr[17] = 1.224744871391589*fskin[38]+0.7071067811865475*fskin[33]; 
  f_cr[18] = 1.224744871391589*fskin[39]+0.7071067811865475*fskin[34]; 
  f_cr[19] = 1.224744871391589*fskin[42]+0.7071067811865475*fskin[36]; 
  f_cr[20] = 1.224744871391589*fskin[43]+0.7071067811865475*fskin[37]; 
  f_cr[21] = 1.224744871391589*fskin[45]+0.7071067811865475*fskin[40]; 
  f_cr[22] = 1.224744871391589*fskin[46]+0.7071067811865475*fskin[41]; 
  f_cr[23] = 1.224744871391589*fskin[47]+0.7071067811865475*fskin[44]; 

  f_rl[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[3]; 
  f_rl[1] = 0.7071067811865475*fedge[1]-1.224744871391589*fedge[7]; 
  f_rl[2] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[8]; 
  f_rl[3] = 0.7071067811865475*fedge[4]-1.224744871391589*fedge[11]; 
  f_rl[4] = 0.7071067811865475*fedge[5]-1.224744871391589*fedge[14]; 
  f_rl[5] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[16]; 
  f_rl[6] = 0.7071067811865475*fedge[9]-1.224744871391589*fedge[18]; 
  f_rl[7] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[19]; 
  f_rl[8] = 0.7071067811865475*fedge[12]-1.224744871391589*fedge[21]; 
  f_rl[9] = 0.7071067811865475*fedge[13]-1.224744871391589*fedge[22]; 
  f_rl[10] = 0.7071067811865475*fedge[15]-1.224744871391589*fedge[25]; 
  f_rl[11] = 0.7071067811865475*fedge[17]-1.224744871391589*fedge[26]; 
  f_rl[12] = 0.7071067811865475*fedge[20]-1.224744871391589*fedge[27]; 
  f_rl[13] = 0.7071067811865475*fedge[23]-1.224744871391589*fedge[29]; 
  f_rl[14] = 0.7071067811865475*fedge[24]-1.224744871391589*fedge[30]; 
  f_rl[15] = 0.7071067811865475*fedge[28]-1.224744871391589*fedge[31]; 
  f_rl[16] = 0.7071067811865475*fedge[32]-1.224744871391589*fedge[35]; 
  f_rl[17] = 0.7071067811865475*fedge[33]-1.224744871391589*fedge[38]; 
  f_rl[18] = 0.7071067811865475*fedge[34]-1.224744871391589*fedge[39]; 
  f_rl[19] = 0.7071067811865475*fedge[36]-1.224744871391589*fedge[42]; 
  f_rl[20] = 0.7071067811865475*fedge[37]-1.224744871391589*fedge[43]; 
  f_rl[21] = 0.7071067811865475*fedge[40]-1.224744871391589*fedge[45]; 
  f_rl[22] = 0.7071067811865475*fedge[41]-1.224744871391589*fedge[46]; 
  f_rl[23] = 0.7071067811865475*fedge[44]-1.224744871391589*fedge[47]; 

  fUpR[0] = sgn_alphaUpR[5]*(0.25*f_cr[17]-0.25*f_rl[17])+sgn_alphaUpR[4]*(0.25*f_cr[16]-0.25*f_rl[16])+sgn_alphaUpR[3]*(0.25*f_cr[6]-0.25*f_rl[6])+sgn_alphaUpR[2]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[1]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = sgn_alphaUpR[4]*(0.2500000000000001*f_cr[17]-0.2500000000000001*f_rl[17])+sgn_alphaUpR[5]*(0.2500000000000001*f_cr[16]-0.2500000000000001*f_rl[16])+sgn_alphaUpR[2]*(0.25*f_cr[6]-0.25*f_rl[6])+(0.25*f_cr[3]-0.25*f_rl[3])*sgn_alphaUpR[3]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[1]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = sgn_alphaUpR[5]*(0.2500000000000001*f_cr[20]-0.2500000000000001*f_rl[20])+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[18]-0.2500000000000001*f_rl[18])+sgn_alphaUpR[3]*(0.25*f_cr[11]-0.25*f_rl[11])+sgn_alphaUpR[2]*(0.25*f_cr[7]-0.25*f_rl[7])+sgn_alphaUpR[1]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[2]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[17]-0.223606797749979*f_rl[17])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[16]-0.223606797749979*f_rl[16])+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[6]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[6]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[4]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[3]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[3]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[3]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[2]; 
  fUpR[4] = sgn_alphaUpR[5]*(0.2500000000000001*f_cr[21]-0.2500000000000001*f_rl[21])+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[19]-0.2500000000000001*f_rl[19])+sgn_alphaUpR[3]*(0.25*f_cr[13]-0.25*f_rl[13])+sgn_alphaUpR[2]*(0.25*f_cr[10]-0.25*f_rl[10])+sgn_alphaUpR[1]*(0.25*f_cr[8]-0.25*f_rl[8])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[4]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[4]; 
  fUpR[5] = sgn_alphaUpR[4]*(0.25*f_cr[20]-0.25*f_rl[20])+sgn_alphaUpR[5]*(0.25*f_cr[18]-0.25*f_rl[18])+sgn_alphaUpR[2]*(0.25*f_cr[11]-0.25*f_rl[11])+sgn_alphaUpR[3]*(0.25*f_cr[7]-0.25*f_rl[7])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[5]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[5]+sgn_alphaUpR[1]*(0.25*f_cr[2]-0.25*f_rl[2]); 
  fUpR[6] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[17]-0.223606797749979*f_rl[17])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[16]-0.223606797749979*f_rl[16])+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[6]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[5]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[2]; 
  fUpR[7] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[20]-0.223606797749979*f_rl[20])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[18]-0.223606797749979*f_rl[18])+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[11]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[11]+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[7]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[7]+sgn_alphaUpR[3]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[2]; 
  fUpR[8] = sgn_alphaUpR[4]*(0.25*f_cr[21]-0.25*f_rl[21])+sgn_alphaUpR[5]*(0.25*f_cr[19]-0.25*f_rl[19])+sgn_alphaUpR[2]*(0.25*f_cr[13]-0.25*f_rl[13])+sgn_alphaUpR[3]*(0.25*f_cr[10]-0.25*f_rl[10])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[8]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[8]+sgn_alphaUpR[1]*(0.25*f_cr[4]-0.25*f_rl[4]); 
  fUpR[9] = sgn_alphaUpR[5]*(0.25*f_cr[23]-0.25*f_rl[23])+sgn_alphaUpR[4]*(0.25*f_cr[22]-0.25*f_rl[22])+sgn_alphaUpR[3]*(0.25*f_cr[15]-0.25*f_rl[15])+sgn_alphaUpR[2]*(0.25*f_cr[14]-0.25*f_rl[14])+sgn_alphaUpR[1]*(0.25*f_cr[12]-0.25*f_rl[12])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[9]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[9]; 
  fUpR[10] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[21]-0.223606797749979*f_rl[21])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[19]-0.223606797749979*f_rl[19])+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[13]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[13]+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[10]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[10]+sgn_alphaUpR[3]*(0.25*f_cr[8]-0.25*f_rl[8])+sgn_alphaUpR[2]*(0.25*f_cr[4]-0.25*f_rl[4]); 
  fUpR[11] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[20]-0.223606797749979*f_rl[20])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[18]-0.223606797749979*f_rl[18])+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[11]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[11]+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[7]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[7]+sgn_alphaUpR[2]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[3]; 
  fUpR[12] = sgn_alphaUpR[4]*(0.2500000000000001*f_cr[23]-0.2500000000000001*f_rl[23])+sgn_alphaUpR[5]*(0.2500000000000001*f_cr[22]-0.2500000000000001*f_rl[22])+sgn_alphaUpR[2]*(0.25*f_cr[15]-0.25*f_rl[15])+sgn_alphaUpR[3]*(0.25*f_cr[14]-0.25*f_rl[14])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[12]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[12]+sgn_alphaUpR[1]*(0.25*f_cr[9]-0.25*f_rl[9]); 
  fUpR[13] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[21]-0.223606797749979*f_rl[21])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[19]-0.223606797749979*f_rl[19])+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[13]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[13]+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[10]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[10]+sgn_alphaUpR[2]*(0.25*f_cr[8]-0.25*f_rl[8])+sgn_alphaUpR[3]*(0.25*f_cr[4]-0.25*f_rl[4]); 
  fUpR[14] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[23]-0.223606797749979*f_rl[23])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[22]-0.223606797749979*f_rl[22])+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[15]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[15]+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[14]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[14]+sgn_alphaUpR[3]*(0.25*f_cr[12]-0.25*f_rl[12])+sgn_alphaUpR[2]*(0.25*f_cr[9]-0.25*f_rl[9]); 
  fUpR[15] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[23]-0.223606797749979*f_rl[23])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[22]-0.223606797749979*f_rl[22])+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[15]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[15]+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[14]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[14]+sgn_alphaUpR[2]*(0.25*f_cr[12]-0.25*f_rl[12])+sgn_alphaUpR[3]*(0.25*f_cr[9]-0.25*f_rl[9]); 
  fUpR[16] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[17]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[17]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[16]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[16]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.2500000000000001*f_cr[1]-0.2500000000000001*f_rl[1])*sgn_alphaUpR[5]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3]); 
  fUpR[17] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[17]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[17]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[16]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[16]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[5]+(0.2500000000000001*f_cr[1]-0.2500000000000001*f_rl[1])*sgn_alphaUpR[4]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[3]; 
  fUpR[18] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[20]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[20]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[18]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[18]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[11]-0.223606797749979*f_rl[11])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+(0.25*f_cr[5]-0.25*f_rl[5])*sgn_alphaUpR[5]+(0.2500000000000001*f_cr[2]-0.2500000000000001*f_rl[2])*sgn_alphaUpR[4]; 
  fUpR[19] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[21]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[21]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[19]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[19]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[13]-0.223606797749979*f_rl[13])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[10]-0.223606797749979*f_rl[10])+sgn_alphaUpR[5]*(0.25*f_cr[8]-0.25*f_rl[8])+(0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])*sgn_alphaUpR[4]; 
  fUpR[20] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[20]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[20]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[18]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[18]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[11]-0.223606797749979*f_rl[11])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+(0.2500000000000001*f_cr[2]-0.2500000000000001*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.25*f_cr[5]-0.25*f_rl[5]); 
  fUpR[21] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[21]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[21]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[19]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[19]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[13]-0.223606797749979*f_rl[13])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[10]-0.223606797749979*f_rl[10])+sgn_alphaUpR[4]*(0.25*f_cr[8]-0.25*f_rl[8])+(0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])*sgn_alphaUpR[5]; 
  fUpR[22] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[23]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[23]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[22]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[22]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[15]-0.223606797749979*f_rl[15])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[14]-0.223606797749979*f_rl[14])+sgn_alphaUpR[5]*(0.2500000000000001*f_cr[12]-0.2500000000000001*f_rl[12])+sgn_alphaUpR[4]*(0.25*f_cr[9]-0.25*f_rl[9]); 
  fUpR[23] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[23]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[23]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[22]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[22]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[15]-0.223606797749979*f_rl[15])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[14]-0.223606797749979*f_rl[14])+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[12]-0.2500000000000001*f_rl[12])+sgn_alphaUpR[5]*(0.25*f_cr[9]-0.25*f_rl[9]); 

  } 
  double GhatR[24] = {0.};
  GhatR[0] = 0.25*(alphaR[6]*fUpR[6]+alphaR[3]*fUpR[3]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.25*(alphaR[6]*fUpR[11]+alphaR[3]*fUpR[7]+alphaR[1]*fUpR[5]+alphaR[0]*fUpR[2]); 
  GhatR[3] = 0.223606797749979*alphaR[6]*fUpR[17]+0.223606797749979*alphaR[3]*fUpR[16]+0.25*(alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.25*(alphaR[6]*fUpR[13]+alphaR[3]*fUpR[10]+alphaR[1]*fUpR[8]+alphaR[0]*fUpR[4]); 
  GhatR[5] = 0.25*(alphaR[3]*fUpR[11]+alphaR[6]*fUpR[7]+alphaR[0]*fUpR[5]+alphaR[1]*fUpR[2]); 
  GhatR[6] = 0.223606797749979*alphaR[3]*fUpR[17]+0.223606797749979*alphaR[6]*fUpR[16]+0.25*(alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.223606797749979*alphaR[6]*fUpR[20]+0.223606797749979*alphaR[3]*fUpR[18]+0.25*(alphaR[1]*fUpR[11]+alphaR[0]*fUpR[7]+fUpR[5]*alphaR[6]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[3]*fUpR[13]+alphaR[6]*fUpR[10]+alphaR[0]*fUpR[8]+alphaR[1]*fUpR[4]); 
  GhatR[9] = 0.25*(alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]+alphaR[1]*fUpR[12]+alphaR[0]*fUpR[9]); 
  GhatR[10] = 0.223606797749979*alphaR[6]*fUpR[21]+0.223606797749979*alphaR[3]*fUpR[19]+0.25*(alphaR[1]*fUpR[13]+alphaR[0]*fUpR[10]+alphaR[6]*fUpR[8]+alphaR[3]*fUpR[4]); 
  GhatR[11] = 0.223606797749979*alphaR[3]*fUpR[20]+0.223606797749979*alphaR[6]*fUpR[18]+0.25*(alphaR[0]*fUpR[11]+alphaR[1]*fUpR[7]+fUpR[2]*alphaR[6]+alphaR[3]*fUpR[5]); 
  GhatR[12] = 0.25*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]+alphaR[0]*fUpR[12]+alphaR[1]*fUpR[9]); 
  GhatR[13] = 0.223606797749979*alphaR[3]*fUpR[21]+0.223606797749979*alphaR[6]*fUpR[19]+0.25*(alphaR[0]*fUpR[13]+alphaR[1]*fUpR[10]+alphaR[3]*fUpR[8]+fUpR[4]*alphaR[6]); 
  GhatR[14] = 0.223606797749979*alphaR[6]*fUpR[23]+0.223606797749979*alphaR[3]*fUpR[22]+0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[6]*fUpR[12]+alphaR[3]*fUpR[9]); 
  GhatR[15] = 0.223606797749979*alphaR[3]*fUpR[23]+0.223606797749979*alphaR[6]*fUpR[22]+0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[3]*fUpR[12]+alphaR[6]*fUpR[9]); 
  GhatR[16] = 0.2500000000000001*alphaR[1]*fUpR[17]+0.25*alphaR[0]*fUpR[16]+0.223606797749979*(alphaR[6]*fUpR[6]+alphaR[3]*fUpR[3]); 
  GhatR[17] = 0.25*alphaR[0]*fUpR[17]+0.2500000000000001*alphaR[1]*fUpR[16]+0.223606797749979*(alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]); 
  GhatR[18] = 0.2500000000000001*alphaR[1]*fUpR[20]+0.25*alphaR[0]*fUpR[18]+0.223606797749979*(alphaR[6]*fUpR[11]+alphaR[3]*fUpR[7]); 
  GhatR[19] = 0.2500000000000001*alphaR[1]*fUpR[21]+0.25*alphaR[0]*fUpR[19]+0.223606797749979*(alphaR[6]*fUpR[13]+alphaR[3]*fUpR[10]); 
  GhatR[20] = 0.25*alphaR[0]*fUpR[20]+0.2500000000000001*alphaR[1]*fUpR[18]+0.223606797749979*(alphaR[3]*fUpR[11]+alphaR[6]*fUpR[7]); 
  GhatR[21] = 0.25*alphaR[0]*fUpR[21]+0.2500000000000001*alphaR[1]*fUpR[19]+0.223606797749979*(alphaR[3]*fUpR[13]+alphaR[6]*fUpR[10]); 
  GhatR[22] = 0.2500000000000001*alphaR[1]*fUpR[23]+0.25*alphaR[0]*fUpR[22]+0.223606797749979*(alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]); 
  GhatR[23] = 0.25*alphaR[0]*fUpR[23]+0.2500000000000001*alphaR[1]*fUpR[22]+0.223606797749979*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdz2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdz2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdz2; 
  out[3] += -1.224744871391589*GhatR[0]*rdz2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdz2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdz2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdz2; 
  out[7] += -1.224744871391589*GhatR[1]*rdz2; 
  out[8] += -1.224744871391589*GhatR[2]*rdz2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdz2; 
  out[10] += -0.7071067811865475*GhatR[7]*rdz2; 
  out[11] += -1.224744871391589*GhatR[3]*rdz2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdz2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdz2; 
  out[14] += -1.224744871391589*GhatR[4]*rdz2; 
  out[15] += -0.7071067811865475*GhatR[10]*rdz2; 
  out[16] += -1.224744871391589*GhatR[5]*rdz2; 
  out[17] += -0.7071067811865475*GhatR[11]*rdz2; 
  out[18] += -1.224744871391589*GhatR[6]*rdz2; 
  out[19] += -1.224744871391589*GhatR[7]*rdz2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdz2; 
  out[21] += -1.224744871391589*GhatR[8]*rdz2; 
  out[22] += -1.224744871391589*GhatR[9]*rdz2; 
  out[23] += -0.7071067811865475*GhatR[13]*rdz2; 
  out[24] += -0.7071067811865475*GhatR[14]*rdz2; 
  out[25] += -1.224744871391589*GhatR[10]*rdz2; 
  out[26] += -1.224744871391589*GhatR[11]*rdz2; 
  out[27] += -1.224744871391589*GhatR[12]*rdz2; 
  out[28] += -0.7071067811865475*GhatR[15]*rdz2; 
  out[29] += -1.224744871391589*GhatR[13]*rdz2; 
  out[30] += -1.224744871391589*GhatR[14]*rdz2; 
  out[31] += -1.224744871391589*GhatR[15]*rdz2; 
  out[32] += -0.7071067811865475*GhatR[16]*rdz2; 
  out[33] += -0.7071067811865475*GhatR[17]*rdz2; 
  out[34] += -0.7071067811865475*GhatR[18]*rdz2; 
  out[35] += -1.224744871391589*GhatR[16]*rdz2; 
  out[36] += -0.7071067811865475*GhatR[19]*rdz2; 
  out[37] += -0.7071067811865475*GhatR[20]*rdz2; 
  out[38] += -1.224744871391589*GhatR[17]*rdz2; 
  out[39] += -1.224744871391589*GhatR[18]*rdz2; 
  out[40] += -0.7071067811865475*GhatR[21]*rdz2; 
  out[41] += -0.7071067811865475*GhatR[22]*rdz2; 
  out[42] += -1.224744871391589*GhatR[19]*rdz2; 
  out[43] += -1.224744871391589*GhatR[20]*rdz2; 
  out[44] += -0.7071067811865475*GhatR[23]*rdz2; 
  out[45] += -1.224744871391589*GhatR[21]*rdz2; 
  out[46] += -1.224744871391589*GhatR[22]*rdz2; 
  out[47] += -1.224744871391589*GhatR[23]*rdz2; 

  } else { 

  double fUpL[24] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fedge[3]+0.7071067811865475*fedge[0]; 
  fUpL[1] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[1]; 
  fUpL[2] = 1.224744871391589*fedge[8]+0.7071067811865475*fedge[2]; 
  fUpL[3] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[4]; 
  fUpL[4] = 1.224744871391589*fedge[14]+0.7071067811865475*fedge[5]; 
  fUpL[5] = 1.224744871391589*fedge[16]+0.7071067811865475*fedge[6]; 
  fUpL[6] = 1.224744871391589*fedge[18]+0.7071067811865475*fedge[9]; 
  fUpL[7] = 1.224744871391589*fedge[19]+0.7071067811865475*fedge[10]; 
  fUpL[8] = 1.224744871391589*fedge[21]+0.7071067811865475*fedge[12]; 
  fUpL[9] = 1.224744871391589*fedge[22]+0.7071067811865475*fedge[13]; 
  fUpL[10] = 1.224744871391589*fedge[25]+0.7071067811865475*fedge[15]; 
  fUpL[11] = 1.224744871391589*fedge[26]+0.7071067811865475*fedge[17]; 
  fUpL[12] = 1.224744871391589*fedge[27]+0.7071067811865475*fedge[20]; 
  fUpL[13] = 1.224744871391589*fedge[29]+0.7071067811865475*fedge[23]; 
  fUpL[14] = 1.224744871391589*fedge[30]+0.7071067811865475*fedge[24]; 
  fUpL[15] = 1.224744871391589*fedge[31]+0.7071067811865475*fedge[28]; 
  fUpL[16] = 1.224744871391589*fedge[35]+0.7071067811865475*fedge[32]; 
  fUpL[17] = 1.224744871391589*fedge[38]+0.7071067811865475*fedge[33]; 
  fUpL[18] = 1.224744871391589*fedge[39]+0.7071067811865475*fedge[34]; 
  fUpL[19] = 1.224744871391589*fedge[42]+0.7071067811865475*fedge[36]; 
  fUpL[20] = 1.224744871391589*fedge[43]+0.7071067811865475*fedge[37]; 
  fUpL[21] = 1.224744871391589*fedge[45]+0.7071067811865475*fedge[40]; 
  fUpL[22] = 1.224744871391589*fedge[46]+0.7071067811865475*fedge[41]; 
  fUpL[23] = 1.224744871391589*fedge[47]+0.7071067811865475*fedge[44]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[3]; 
  fUpL[1] = 0.7071067811865475*fskin[1]-1.224744871391589*fskin[7]; 
  fUpL[2] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[8]; 
  fUpL[3] = 0.7071067811865475*fskin[4]-1.224744871391589*fskin[11]; 
  fUpL[4] = 0.7071067811865475*fskin[5]-1.224744871391589*fskin[14]; 
  fUpL[5] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[16]; 
  fUpL[6] = 0.7071067811865475*fskin[9]-1.224744871391589*fskin[18]; 
  fUpL[7] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[19]; 
  fUpL[8] = 0.7071067811865475*fskin[12]-1.224744871391589*fskin[21]; 
  fUpL[9] = 0.7071067811865475*fskin[13]-1.224744871391589*fskin[22]; 
  fUpL[10] = 0.7071067811865475*fskin[15]-1.224744871391589*fskin[25]; 
  fUpL[11] = 0.7071067811865475*fskin[17]-1.224744871391589*fskin[26]; 
  fUpL[12] = 0.7071067811865475*fskin[20]-1.224744871391589*fskin[27]; 
  fUpL[13] = 0.7071067811865475*fskin[23]-1.224744871391589*fskin[29]; 
  fUpL[14] = 0.7071067811865475*fskin[24]-1.224744871391589*fskin[30]; 
  fUpL[15] = 0.7071067811865475*fskin[28]-1.224744871391589*fskin[31]; 
  fUpL[16] = 0.7071067811865475*fskin[32]-1.224744871391589*fskin[35]; 
  fUpL[17] = 0.7071067811865475*fskin[33]-1.224744871391589*fskin[38]; 
  fUpL[18] = 0.7071067811865475*fskin[34]-1.224744871391589*fskin[39]; 
  fUpL[19] = 0.7071067811865475*fskin[36]-1.224744871391589*fskin[42]; 
  fUpL[20] = 0.7071067811865475*fskin[37]-1.224744871391589*fskin[43]; 
  fUpL[21] = 0.7071067811865475*fskin[40]-1.224744871391589*fskin[45]; 
  fUpL[22] = 0.7071067811865475*fskin[41]-1.224744871391589*fskin[46]; 
  fUpL[23] = 0.7071067811865475*fskin[44]-1.224744871391589*fskin[47]; 
    } 
  } else { 
  double f_lr[24] = {0.};
  double f_cl[24] = {0.};
  double sgn_alphaUpL[6] = {0.};
  sgn_alphaUpL[0] = 0.2777777777777778*sgn_alpha_surfL[5]+0.4444444444444444*sgn_alpha_surfL[4]+0.2777777777777778*sgn_alpha_surfL[3]+0.2777777777777778*sgn_alpha_surfL[2]+0.4444444444444444*sgn_alpha_surfL[1]+0.2777777777777778*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[1] = 0.2777777777777778*sgn_alpha_surfL[5]+0.4444444444444444*sgn_alpha_surfL[4]+0.2777777777777778*sgn_alpha_surfL[3]-0.2777777777777778*sgn_alpha_surfL[2]-0.4444444444444444*sgn_alpha_surfL[1]-0.2777777777777778*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[2] = 0.3726779962499649*sgn_alpha_surfL[5]-0.3726779962499649*sgn_alpha_surfL[3]+0.3726779962499649*sgn_alpha_surfL[2]-0.3726779962499649*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[3] = 0.3726779962499649*sgn_alpha_surfL[5]-0.3726779962499649*sgn_alpha_surfL[3]-0.3726779962499649*sgn_alpha_surfL[2]+0.3726779962499649*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[4] = 0.2484519974999766*sgn_alpha_surfL[5]-0.4969039949999532*sgn_alpha_surfL[4]+0.2484519974999766*sgn_alpha_surfL[3]+0.2484519974999766*sgn_alpha_surfL[2]-0.4969039949999532*sgn_alpha_surfL[1]+0.2484519974999766*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[5] = 0.2484519974999767*sgn_alpha_surfL[5]-0.4969039949999535*sgn_alpha_surfL[4]+0.2484519974999767*sgn_alpha_surfL[3]-0.2484519974999767*sgn_alpha_surfL[2]+0.4969039949999535*sgn_alpha_surfL[1]-0.2484519974999767*sgn_alpha_surfL[0]; 

  f_lr[0] = 1.224744871391589*fedge[3]+0.7071067811865475*fedge[0]; 
  f_lr[1] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[1]; 
  f_lr[2] = 1.224744871391589*fedge[8]+0.7071067811865475*fedge[2]; 
  f_lr[3] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[4]; 
  f_lr[4] = 1.224744871391589*fedge[14]+0.7071067811865475*fedge[5]; 
  f_lr[5] = 1.224744871391589*fedge[16]+0.7071067811865475*fedge[6]; 
  f_lr[6] = 1.224744871391589*fedge[18]+0.7071067811865475*fedge[9]; 
  f_lr[7] = 1.224744871391589*fedge[19]+0.7071067811865475*fedge[10]; 
  f_lr[8] = 1.224744871391589*fedge[21]+0.7071067811865475*fedge[12]; 
  f_lr[9] = 1.224744871391589*fedge[22]+0.7071067811865475*fedge[13]; 
  f_lr[10] = 1.224744871391589*fedge[25]+0.7071067811865475*fedge[15]; 
  f_lr[11] = 1.224744871391589*fedge[26]+0.7071067811865475*fedge[17]; 
  f_lr[12] = 1.224744871391589*fedge[27]+0.7071067811865475*fedge[20]; 
  f_lr[13] = 1.224744871391589*fedge[29]+0.7071067811865475*fedge[23]; 
  f_lr[14] = 1.224744871391589*fedge[30]+0.7071067811865475*fedge[24]; 
  f_lr[15] = 1.224744871391589*fedge[31]+0.7071067811865475*fedge[28]; 
  f_lr[16] = 1.224744871391589*fedge[35]+0.7071067811865475*fedge[32]; 
  f_lr[17] = 1.224744871391589*fedge[38]+0.7071067811865475*fedge[33]; 
  f_lr[18] = 1.224744871391589*fedge[39]+0.7071067811865475*fedge[34]; 
  f_lr[19] = 1.224744871391589*fedge[42]+0.7071067811865475*fedge[36]; 
  f_lr[20] = 1.224744871391589*fedge[43]+0.7071067811865475*fedge[37]; 
  f_lr[21] = 1.224744871391589*fedge[45]+0.7071067811865475*fedge[40]; 
  f_lr[22] = 1.224744871391589*fedge[46]+0.7071067811865475*fedge[41]; 
  f_lr[23] = 1.224744871391589*fedge[47]+0.7071067811865475*fedge[44]; 

  f_cl[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[3]; 
  f_cl[1] = 0.7071067811865475*fskin[1]-1.224744871391589*fskin[7]; 
  f_cl[2] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[8]; 
  f_cl[3] = 0.7071067811865475*fskin[4]-1.224744871391589*fskin[11]; 
  f_cl[4] = 0.7071067811865475*fskin[5]-1.224744871391589*fskin[14]; 
  f_cl[5] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[16]; 
  f_cl[6] = 0.7071067811865475*fskin[9]-1.224744871391589*fskin[18]; 
  f_cl[7] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[19]; 
  f_cl[8] = 0.7071067811865475*fskin[12]-1.224744871391589*fskin[21]; 
  f_cl[9] = 0.7071067811865475*fskin[13]-1.224744871391589*fskin[22]; 
  f_cl[10] = 0.7071067811865475*fskin[15]-1.224744871391589*fskin[25]; 
  f_cl[11] = 0.7071067811865475*fskin[17]-1.224744871391589*fskin[26]; 
  f_cl[12] = 0.7071067811865475*fskin[20]-1.224744871391589*fskin[27]; 
  f_cl[13] = 0.7071067811865475*fskin[23]-1.224744871391589*fskin[29]; 
  f_cl[14] = 0.7071067811865475*fskin[24]-1.224744871391589*fskin[30]; 
  f_cl[15] = 0.7071067811865475*fskin[28]-1.224744871391589*fskin[31]; 
  f_cl[16] = 0.7071067811865475*fskin[32]-1.224744871391589*fskin[35]; 
  f_cl[17] = 0.7071067811865475*fskin[33]-1.224744871391589*fskin[38]; 
  f_cl[18] = 0.7071067811865475*fskin[34]-1.224744871391589*fskin[39]; 
  f_cl[19] = 0.7071067811865475*fskin[36]-1.224744871391589*fskin[42]; 
  f_cl[20] = 0.7071067811865475*fskin[37]-1.224744871391589*fskin[43]; 
  f_cl[21] = 0.7071067811865475*fskin[40]-1.224744871391589*fskin[45]; 
  f_cl[22] = 0.7071067811865475*fskin[41]-1.224744871391589*fskin[46]; 
  f_cl[23] = 0.7071067811865475*fskin[44]-1.224744871391589*fskin[47]; 

  fUpL[0] = sgn_alphaUpL[5]*(0.25*f_lr[17]-0.25*f_cl[17])+sgn_alphaUpL[4]*(0.25*f_lr[16]-0.25*f_cl[16])+sgn_alphaUpL[3]*(0.25*f_lr[6]-0.25*f_cl[6])+sgn_alphaUpL[2]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[1]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = sgn_alphaUpL[4]*(0.2500000000000001*f_lr[17]-0.2500000000000001*f_cl[17])+sgn_alphaUpL[5]*(0.2500000000000001*f_lr[16]-0.2500000000000001*f_cl[16])+sgn_alphaUpL[2]*(0.25*f_lr[6]-0.25*f_cl[6])+(0.25*f_lr[3]-0.25*f_cl[3])*sgn_alphaUpL[3]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[1]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = sgn_alphaUpL[5]*(0.2500000000000001*f_lr[20]-0.2500000000000001*f_cl[20])+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[18]-0.2500000000000001*f_cl[18])+sgn_alphaUpL[3]*(0.25*f_lr[11]-0.25*f_cl[11])+sgn_alphaUpL[2]*(0.25*f_lr[7]-0.25*f_cl[7])+sgn_alphaUpL[1]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[17]-0.223606797749979*f_cl[17])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[16]-0.223606797749979*f_cl[16])+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[6]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[6]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[4]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[3]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[3]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[2]; 
  fUpL[4] = sgn_alphaUpL[5]*(0.2500000000000001*f_lr[21]-0.2500000000000001*f_cl[21])+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[19]-0.2500000000000001*f_cl[19])+sgn_alphaUpL[3]*(0.25*f_lr[13]-0.25*f_cl[13])+sgn_alphaUpL[2]*(0.25*f_lr[10]-0.25*f_cl[10])+sgn_alphaUpL[1]*(0.25*f_lr[8]-0.25*f_cl[8])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[4]; 
  fUpL[5] = sgn_alphaUpL[4]*(0.25*f_lr[20]-0.25*f_cl[20])+sgn_alphaUpL[5]*(0.25*f_lr[18]-0.25*f_cl[18])+sgn_alphaUpL[2]*(0.25*f_lr[11]-0.25*f_cl[11])+sgn_alphaUpL[3]*(0.25*f_lr[7]-0.25*f_cl[7])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[5]+sgn_alphaUpL[1]*(0.25*f_lr[2]-0.25*f_cl[2]); 
  fUpL[6] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[17]-0.223606797749979*f_cl[17])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[16]-0.223606797749979*f_cl[16])+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[6]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[6]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[5]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[2]; 
  fUpL[7] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[20]-0.223606797749979*f_cl[20])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[18]-0.223606797749979*f_cl[18])+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[11]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[11]+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[7]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[7]+sgn_alphaUpL[3]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[2]; 
  fUpL[8] = sgn_alphaUpL[4]*(0.25*f_lr[21]-0.25*f_cl[21])+sgn_alphaUpL[5]*(0.25*f_lr[19]-0.25*f_cl[19])+sgn_alphaUpL[2]*(0.25*f_lr[13]-0.25*f_cl[13])+sgn_alphaUpL[3]*(0.25*f_lr[10]-0.25*f_cl[10])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[8]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[8]+sgn_alphaUpL[1]*(0.25*f_lr[4]-0.25*f_cl[4]); 
  fUpL[9] = sgn_alphaUpL[5]*(0.25*f_lr[23]-0.25*f_cl[23])+sgn_alphaUpL[4]*(0.25*f_lr[22]-0.25*f_cl[22])+sgn_alphaUpL[3]*(0.25*f_lr[15]-0.25*f_cl[15])+sgn_alphaUpL[2]*(0.25*f_lr[14]-0.25*f_cl[14])+sgn_alphaUpL[1]*(0.25*f_lr[12]-0.25*f_cl[12])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[9]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[9]; 
  fUpL[10] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[21]-0.223606797749979*f_cl[21])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[19]-0.223606797749979*f_cl[19])+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[13]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[13]+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[10]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[10]+sgn_alphaUpL[3]*(0.25*f_lr[8]-0.25*f_cl[8])+sgn_alphaUpL[2]*(0.25*f_lr[4]-0.25*f_cl[4]); 
  fUpL[11] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[20]-0.223606797749979*f_cl[20])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[18]-0.223606797749979*f_cl[18])+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[11]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[11]+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[7]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[7]+sgn_alphaUpL[2]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[3]; 
  fUpL[12] = sgn_alphaUpL[4]*(0.2500000000000001*f_lr[23]-0.2500000000000001*f_cl[23])+sgn_alphaUpL[5]*(0.2500000000000001*f_lr[22]-0.2500000000000001*f_cl[22])+sgn_alphaUpL[2]*(0.25*f_lr[15]-0.25*f_cl[15])+sgn_alphaUpL[3]*(0.25*f_lr[14]-0.25*f_cl[14])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[12]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[12]+sgn_alphaUpL[1]*(0.25*f_lr[9]-0.25*f_cl[9]); 
  fUpL[13] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[21]-0.223606797749979*f_cl[21])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[19]-0.223606797749979*f_cl[19])+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[13]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[13]+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[10]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[10]+sgn_alphaUpL[2]*(0.25*f_lr[8]-0.25*f_cl[8])+sgn_alphaUpL[3]*(0.25*f_lr[4]-0.25*f_cl[4]); 
  fUpL[14] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[23]-0.223606797749979*f_cl[23])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[22]-0.223606797749979*f_cl[22])+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[15]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[15]+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[14]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[14]+sgn_alphaUpL[3]*(0.25*f_lr[12]-0.25*f_cl[12])+sgn_alphaUpL[2]*(0.25*f_lr[9]-0.25*f_cl[9]); 
  fUpL[15] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[23]-0.223606797749979*f_cl[23])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[22]-0.223606797749979*f_cl[22])+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[15]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[15]+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[14]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[14]+sgn_alphaUpL[2]*(0.25*f_lr[12]-0.25*f_cl[12])+sgn_alphaUpL[3]*(0.25*f_lr[9]-0.25*f_cl[9]); 
  fUpL[16] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[17]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[17]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[16]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[16]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[5]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3]); 
  fUpL[17] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[17]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[17]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[16]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[16]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[5]+(0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[4]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[3]; 
  fUpL[18] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[20]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[20]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[18]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[18]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[11]-0.223606797749979*f_cl[11])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+(0.25*f_lr[5]-0.25*f_cl[5])*sgn_alphaUpL[5]+(0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[4]; 
  fUpL[19] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[21]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[21]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[19]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[19]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[13]-0.223606797749979*f_cl[13])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[10]-0.223606797749979*f_cl[10])+sgn_alphaUpL[5]*(0.25*f_lr[8]-0.25*f_cl[8])+(0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])*sgn_alphaUpL[4]; 
  fUpL[20] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[20]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[20]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[18]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[18]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[11]-0.223606797749979*f_cl[11])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+(0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.25*f_lr[5]-0.25*f_cl[5]); 
  fUpL[21] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[21]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[21]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[19]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[19]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[13]-0.223606797749979*f_cl[13])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[10]-0.223606797749979*f_cl[10])+sgn_alphaUpL[4]*(0.25*f_lr[8]-0.25*f_cl[8])+(0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])*sgn_alphaUpL[5]; 
  fUpL[22] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[23]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[23]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[22]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[22]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[15]-0.223606797749979*f_cl[15])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[14]-0.223606797749979*f_cl[14])+sgn_alphaUpL[5]*(0.2500000000000001*f_lr[12]-0.2500000000000001*f_cl[12])+sgn_alphaUpL[4]*(0.25*f_lr[9]-0.25*f_cl[9]); 
  fUpL[23] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[23]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[23]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[22]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[22]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[15]-0.223606797749979*f_cl[15])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[14]-0.223606797749979*f_cl[14])+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[12]-0.2500000000000001*f_cl[12])+sgn_alphaUpL[5]*(0.25*f_lr[9]-0.25*f_cl[9]); 

  } 
  double GhatL[24] = {0.};
  GhatL[0] = 0.25*(alphaL[6]*fUpL[6]+alphaL[3]*fUpL[3]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.25*(alphaL[6]*fUpL[11]+alphaL[3]*fUpL[7]+alphaL[1]*fUpL[5]+alphaL[0]*fUpL[2]); 
  GhatL[3] = 0.223606797749979*alphaL[6]*fUpL[17]+0.223606797749979*alphaL[3]*fUpL[16]+0.25*(alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.25*(alphaL[6]*fUpL[13]+alphaL[3]*fUpL[10]+alphaL[1]*fUpL[8]+alphaL[0]*fUpL[4]); 
  GhatL[5] = 0.25*(alphaL[3]*fUpL[11]+alphaL[6]*fUpL[7]+alphaL[0]*fUpL[5]+alphaL[1]*fUpL[2]); 
  GhatL[6] = 0.223606797749979*alphaL[3]*fUpL[17]+0.223606797749979*alphaL[6]*fUpL[16]+0.25*(alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.223606797749979*alphaL[6]*fUpL[20]+0.223606797749979*alphaL[3]*fUpL[18]+0.25*(alphaL[1]*fUpL[11]+alphaL[0]*fUpL[7]+fUpL[5]*alphaL[6]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[3]*fUpL[13]+alphaL[6]*fUpL[10]+alphaL[0]*fUpL[8]+alphaL[1]*fUpL[4]); 
  GhatL[9] = 0.25*(alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]+alphaL[1]*fUpL[12]+alphaL[0]*fUpL[9]); 
  GhatL[10] = 0.223606797749979*alphaL[6]*fUpL[21]+0.223606797749979*alphaL[3]*fUpL[19]+0.25*(alphaL[1]*fUpL[13]+alphaL[0]*fUpL[10]+alphaL[6]*fUpL[8]+alphaL[3]*fUpL[4]); 
  GhatL[11] = 0.223606797749979*alphaL[3]*fUpL[20]+0.223606797749979*alphaL[6]*fUpL[18]+0.25*(alphaL[0]*fUpL[11]+alphaL[1]*fUpL[7]+fUpL[2]*alphaL[6]+alphaL[3]*fUpL[5]); 
  GhatL[12] = 0.25*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]+alphaL[0]*fUpL[12]+alphaL[1]*fUpL[9]); 
  GhatL[13] = 0.223606797749979*alphaL[3]*fUpL[21]+0.223606797749979*alphaL[6]*fUpL[19]+0.25*(alphaL[0]*fUpL[13]+alphaL[1]*fUpL[10]+alphaL[3]*fUpL[8]+fUpL[4]*alphaL[6]); 
  GhatL[14] = 0.223606797749979*alphaL[6]*fUpL[23]+0.223606797749979*alphaL[3]*fUpL[22]+0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[6]*fUpL[12]+alphaL[3]*fUpL[9]); 
  GhatL[15] = 0.223606797749979*alphaL[3]*fUpL[23]+0.223606797749979*alphaL[6]*fUpL[22]+0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[3]*fUpL[12]+alphaL[6]*fUpL[9]); 
  GhatL[16] = 0.2500000000000001*alphaL[1]*fUpL[17]+0.25*alphaL[0]*fUpL[16]+0.223606797749979*(alphaL[6]*fUpL[6]+alphaL[3]*fUpL[3]); 
  GhatL[17] = 0.25*alphaL[0]*fUpL[17]+0.2500000000000001*alphaL[1]*fUpL[16]+0.223606797749979*(alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]); 
  GhatL[18] = 0.2500000000000001*alphaL[1]*fUpL[20]+0.25*alphaL[0]*fUpL[18]+0.223606797749979*(alphaL[6]*fUpL[11]+alphaL[3]*fUpL[7]); 
  GhatL[19] = 0.2500000000000001*alphaL[1]*fUpL[21]+0.25*alphaL[0]*fUpL[19]+0.223606797749979*(alphaL[6]*fUpL[13]+alphaL[3]*fUpL[10]); 
  GhatL[20] = 0.25*alphaL[0]*fUpL[20]+0.2500000000000001*alphaL[1]*fUpL[18]+0.223606797749979*(alphaL[3]*fUpL[11]+alphaL[6]*fUpL[7]); 
  GhatL[21] = 0.25*alphaL[0]*fUpL[21]+0.2500000000000001*alphaL[1]*fUpL[19]+0.223606797749979*(alphaL[3]*fUpL[13]+alphaL[6]*fUpL[10]); 
  GhatL[22] = 0.2500000000000001*alphaL[1]*fUpL[23]+0.25*alphaL[0]*fUpL[22]+0.223606797749979*(alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]); 
  GhatL[23] = 0.25*alphaL[0]*fUpL[23]+0.2500000000000001*alphaL[1]*fUpL[22]+0.223606797749979*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdz2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdz2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdz2; 
  out[3] += -1.224744871391589*GhatL[0]*rdz2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdz2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdz2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdz2; 
  out[7] += -1.224744871391589*GhatL[1]*rdz2; 
  out[8] += -1.224744871391589*GhatL[2]*rdz2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdz2; 
  out[10] += 0.7071067811865475*GhatL[7]*rdz2; 
  out[11] += -1.224744871391589*GhatL[3]*rdz2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdz2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdz2; 
  out[14] += -1.224744871391589*GhatL[4]*rdz2; 
  out[15] += 0.7071067811865475*GhatL[10]*rdz2; 
  out[16] += -1.224744871391589*GhatL[5]*rdz2; 
  out[17] += 0.7071067811865475*GhatL[11]*rdz2; 
  out[18] += -1.224744871391589*GhatL[6]*rdz2; 
  out[19] += -1.224744871391589*GhatL[7]*rdz2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdz2; 
  out[21] += -1.224744871391589*GhatL[8]*rdz2; 
  out[22] += -1.224744871391589*GhatL[9]*rdz2; 
  out[23] += 0.7071067811865475*GhatL[13]*rdz2; 
  out[24] += 0.7071067811865475*GhatL[14]*rdz2; 
  out[25] += -1.224744871391589*GhatL[10]*rdz2; 
  out[26] += -1.224744871391589*GhatL[11]*rdz2; 
  out[27] += -1.224744871391589*GhatL[12]*rdz2; 
  out[28] += 0.7071067811865475*GhatL[15]*rdz2; 
  out[29] += -1.224744871391589*GhatL[13]*rdz2; 
  out[30] += -1.224744871391589*GhatL[14]*rdz2; 
  out[31] += -1.224744871391589*GhatL[15]*rdz2; 
  out[32] += 0.7071067811865475*GhatL[16]*rdz2; 
  out[33] += 0.7071067811865475*GhatL[17]*rdz2; 
  out[34] += 0.7071067811865475*GhatL[18]*rdz2; 
  out[35] += -1.224744871391589*GhatL[16]*rdz2; 
  out[36] += 0.7071067811865475*GhatL[19]*rdz2; 
  out[37] += 0.7071067811865475*GhatL[20]*rdz2; 
  out[38] += -1.224744871391589*GhatL[17]*rdz2; 
  out[39] += -1.224744871391589*GhatL[18]*rdz2; 
  out[40] += 0.7071067811865475*GhatL[21]*rdz2; 
  out[41] += 0.7071067811865475*GhatL[22]*rdz2; 
  out[42] += -1.224744871391589*GhatL[19]*rdz2; 
  out[43] += -1.224744871391589*GhatL[20]*rdz2; 
  out[44] += 0.7071067811865475*GhatL[23]*rdz2; 
  out[45] += -1.224744871391589*GhatL[21]*rdz2; 
  out[46] += -1.224744871391589*GhatL[22]*rdz2; 
  out[47] += -1.224744871391589*GhatL[23]*rdz2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.375*rdz2*cflFreq; 

} 
