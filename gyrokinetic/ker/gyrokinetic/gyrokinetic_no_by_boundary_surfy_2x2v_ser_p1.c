#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
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

  double rdz2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_skin[12];
  const double *alphaR = &alpha_surf_edge[12];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[12];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[12];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[1];

  if (edge == -1) { 

  double fUpR[12] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  fUpR[1] = 1.224744871391589*fskin[5]+0.7071067811865475*fskin[1]; 
  fUpR[2] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[3]; 
  fUpR[3] = 1.224744871391589*fskin[9]+0.7071067811865475*fskin[4]; 
  fUpR[4] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[6]; 
  fUpR[5] = 1.224744871391589*fskin[12]+0.7071067811865475*fskin[8]; 
  fUpR[6] = 1.224744871391589*fskin[14]+0.7071067811865475*fskin[10]; 
  fUpR[7] = 1.224744871391589*fskin[15]+0.7071067811865475*fskin[13]; 
  fUpR[8] = 1.224744871391589*fskin[18]+0.7071067811865475*fskin[16]; 
  fUpR[9] = 1.224744871391589*fskin[20]+0.7071067811865475*fskin[17]; 
  fUpR[10] = 1.224744871391589*fskin[22]+0.7071067811865475*fskin[19]; 
  fUpR[11] = 1.224744871391589*fskin[23]+0.7071067811865475*fskin[21]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[2]; 
  fUpR[1] = 0.7071067811865475*fedge[1]-1.224744871391589*fedge[5]; 
  fUpR[2] = 0.7071067811865475*fedge[3]-1.224744871391589*fedge[7]; 
  fUpR[3] = 0.7071067811865475*fedge[4]-1.224744871391589*fedge[9]; 
  fUpR[4] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[11]; 
  fUpR[5] = 0.7071067811865475*fedge[8]-1.224744871391589*fedge[12]; 
  fUpR[6] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[14]; 
  fUpR[7] = 0.7071067811865475*fedge[13]-1.224744871391589*fedge[15]; 
  fUpR[8] = 0.7071067811865475*fedge[16]-1.224744871391589*fedge[18]; 
  fUpR[9] = 0.7071067811865475*fedge[17]-1.224744871391589*fedge[20]; 
  fUpR[10] = 0.7071067811865475*fedge[19]-1.224744871391589*fedge[22]; 
  fUpR[11] = 0.7071067811865475*fedge[21]-1.224744871391589*fedge[23]; 
    } 
  } else { 
  double f_cr[12] = {0.};
  double f_rl[12] = {0.};
  double sgn_alphaUpR[6] = {0.};
  sgn_alphaUpR[0] = 0.2777777777777778*sgn_alpha_surfR[5]+0.4444444444444444*sgn_alpha_surfR[4]+0.2777777777777778*sgn_alpha_surfR[3]+0.2777777777777778*sgn_alpha_surfR[2]+0.4444444444444444*sgn_alpha_surfR[1]+0.2777777777777778*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[1] = 0.2777777777777778*sgn_alpha_surfR[5]+0.4444444444444444*sgn_alpha_surfR[4]+0.2777777777777778*sgn_alpha_surfR[3]-0.2777777777777778*sgn_alpha_surfR[2]-0.4444444444444444*sgn_alpha_surfR[1]-0.2777777777777778*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[2] = 0.3726779962499649*sgn_alpha_surfR[5]-0.3726779962499649*sgn_alpha_surfR[3]+0.3726779962499649*sgn_alpha_surfR[2]-0.3726779962499649*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[3] = 0.3726779962499649*sgn_alpha_surfR[5]-0.3726779962499649*sgn_alpha_surfR[3]-0.3726779962499649*sgn_alpha_surfR[2]+0.3726779962499649*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[4] = 0.2484519974999766*sgn_alpha_surfR[5]-0.4969039949999532*sgn_alpha_surfR[4]+0.2484519974999766*sgn_alpha_surfR[3]+0.2484519974999766*sgn_alpha_surfR[2]-0.4969039949999532*sgn_alpha_surfR[1]+0.2484519974999766*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[5] = 0.2484519974999767*sgn_alpha_surfR[5]-0.4969039949999535*sgn_alpha_surfR[4]+0.2484519974999767*sgn_alpha_surfR[3]-0.2484519974999767*sgn_alpha_surfR[2]+0.4969039949999535*sgn_alpha_surfR[1]-0.2484519974999767*sgn_alpha_surfR[0]; 

  f_cr[0] = 1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  f_cr[1] = 1.224744871391589*fskin[5]+0.7071067811865475*fskin[1]; 
  f_cr[2] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[3]; 
  f_cr[3] = 1.224744871391589*fskin[9]+0.7071067811865475*fskin[4]; 
  f_cr[4] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[6]; 
  f_cr[5] = 1.224744871391589*fskin[12]+0.7071067811865475*fskin[8]; 
  f_cr[6] = 1.224744871391589*fskin[14]+0.7071067811865475*fskin[10]; 
  f_cr[7] = 1.224744871391589*fskin[15]+0.7071067811865475*fskin[13]; 
  f_cr[8] = 1.224744871391589*fskin[18]+0.7071067811865475*fskin[16]; 
  f_cr[9] = 1.224744871391589*fskin[20]+0.7071067811865475*fskin[17]; 
  f_cr[10] = 1.224744871391589*fskin[22]+0.7071067811865475*fskin[19]; 
  f_cr[11] = 1.224744871391589*fskin[23]+0.7071067811865475*fskin[21]; 

  f_rl[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[2]; 
  f_rl[1] = 0.7071067811865475*fedge[1]-1.224744871391589*fedge[5]; 
  f_rl[2] = 0.7071067811865475*fedge[3]-1.224744871391589*fedge[7]; 
  f_rl[3] = 0.7071067811865475*fedge[4]-1.224744871391589*fedge[9]; 
  f_rl[4] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[11]; 
  f_rl[5] = 0.7071067811865475*fedge[8]-1.224744871391589*fedge[12]; 
  f_rl[6] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[14]; 
  f_rl[7] = 0.7071067811865475*fedge[13]-1.224744871391589*fedge[15]; 
  f_rl[8] = 0.7071067811865475*fedge[16]-1.224744871391589*fedge[18]; 
  f_rl[9] = 0.7071067811865475*fedge[17]-1.224744871391589*fedge[20]; 
  f_rl[10] = 0.7071067811865475*fedge[19]-1.224744871391589*fedge[22]; 
  f_rl[11] = 0.7071067811865475*fedge[21]-1.224744871391589*fedge[23]; 

  fUpR[0] = sgn_alphaUpR[5]*(0.25*f_cr[9]-0.25*f_rl[9])+sgn_alphaUpR[4]*(0.25*f_cr[8]-0.25*f_rl[8])+sgn_alphaUpR[3]*(0.25*f_cr[4]-0.25*f_rl[4])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[2]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[1]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = sgn_alphaUpR[4]*(0.2500000000000001*f_cr[9]-0.2500000000000001*f_rl[9])+sgn_alphaUpR[5]*(0.2500000000000001*f_cr[8]-0.2500000000000001*f_rl[8])+sgn_alphaUpR[2]*(0.25*f_cr[4]-0.25*f_rl[4])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[3]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[1]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[9]-0.223606797749979*f_rl[9])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[8]-0.223606797749979*f_rl[8])+(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])*sgn_alphaUpR[5]+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.25*f_cr[4]-0.25*f_rl[4])+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[3]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[2]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = sgn_alphaUpR[5]*(0.2500000000000001*f_cr[11]-0.2500000000000001*f_rl[11])+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[10]-0.2500000000000001*f_rl[10])+sgn_alphaUpR[3]*(0.25*f_cr[7]-0.25*f_rl[7])+sgn_alphaUpR[2]*(0.25*f_cr[6]-0.25*f_rl[6])+sgn_alphaUpR[1]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[3]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[9]-0.223606797749979*f_rl[9])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[8]-0.223606797749979*f_rl[8])+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[5]+(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])*sgn_alphaUpR[4]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[4]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[4]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[3]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.25*f_cr[2]-0.25*f_rl[2]); 
  fUpR[5] = sgn_alphaUpR[4]*(0.25*f_cr[11]-0.25*f_rl[11])+sgn_alphaUpR[5]*(0.25*f_cr[10]-0.25*f_rl[10])+sgn_alphaUpR[2]*(0.25*f_cr[7]-0.25*f_rl[7])+sgn_alphaUpR[3]*(0.25*f_cr[6]-0.25*f_rl[6])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[5]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[5]+sgn_alphaUpR[1]*(0.25*f_cr[3]-0.25*f_rl[3]); 
  fUpR[6] = sgn_alphaUpR[3]*(0.223606797749979*f_cr[11]-0.223606797749979*f_rl[11])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[10]-0.223606797749979*f_rl[10])+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[7]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[7]+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[6]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[6]+sgn_alphaUpR[3]*(0.25*f_cr[5]-0.25*f_rl[5])+sgn_alphaUpR[2]*(0.25*f_cr[3]-0.25*f_rl[3]); 
  fUpR[7] = sgn_alphaUpR[2]*(0.223606797749979*f_cr[11]-0.223606797749979*f_rl[11])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[10]-0.223606797749979*f_rl[10])+((-0.223606797749979*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[7]+(0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[7]+((-0.223606797749979*sgn_alphaUpR[5])-0.25*sgn_alphaUpR[1])*f_rl[6]+(0.223606797749979*sgn_alphaUpR[5]+0.25*sgn_alphaUpR[1])*f_cr[6]+sgn_alphaUpR[2]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.25*f_cr[3]-0.25*f_rl[3])*sgn_alphaUpR[3]; 
  fUpR[8] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[9]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[9]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[8]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[8]+(0.2500000000000001*f_cr[1]-0.2500000000000001*f_rl[1])*sgn_alphaUpR[5]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[2]; 
  fUpR[9] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[9]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[9]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[8]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[8]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[5]+(0.2500000000000001*f_cr[1]-0.2500000000000001*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[3]; 
  fUpR[10] = ((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[11]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[11]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[10]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[10]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+sgn_alphaUpR[2]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.25*f_cr[5]-0.25*f_rl[5])*sgn_alphaUpR[5]+(0.2500000000000001*f_cr[3]-0.2500000000000001*f_rl[3])*sgn_alphaUpR[4]; 
  fUpR[11] = ((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[11]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[11]+((-0.159719141249985*sgn_alphaUpR[5])-0.2500000000000001*sgn_alphaUpR[1])*f_rl[10]+(0.159719141249985*sgn_alphaUpR[5]+0.2500000000000001*sgn_alphaUpR[1])*f_cr[10]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+sgn_alphaUpR[3]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.2500000000000001*f_cr[3]-0.2500000000000001*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.25*f_cr[5]-0.25*f_rl[5]); 

  } 
  double GhatR[12] = {0.};
  GhatR[0] = 0.3535533905932737*(alphaR[4]*fUpR[4]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.3535533905932737*(alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.3162277660168379*alphaR[4]*fUpR[9]+0.3162277660168379*alphaR[2]*fUpR[8]+0.3535533905932737*(alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.3535533905932737*(alphaR[4]*fUpR[7]+alphaR[2]*fUpR[6]+alphaR[1]*fUpR[5]+alphaR[0]*fUpR[3]); 
  GhatR[4] = 0.3162277660168379*alphaR[2]*fUpR[9]+0.3162277660168379*alphaR[4]*fUpR[8]+0.3535533905932737*(alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[5] = 0.3535533905932737*(alphaR[2]*fUpR[7]+alphaR[4]*fUpR[6]+alphaR[0]*fUpR[5]+alphaR[1]*fUpR[3]); 
  GhatR[6] = 0.3162277660168379*alphaR[4]*fUpR[11]+0.3162277660168379*alphaR[2]*fUpR[10]+0.3535533905932737*(alphaR[1]*fUpR[7]+alphaR[0]*fUpR[6]+alphaR[4]*fUpR[5]+alphaR[2]*fUpR[3]); 
  GhatR[7] = 0.3162277660168379*alphaR[2]*fUpR[11]+0.3162277660168379*alphaR[4]*fUpR[10]+0.3535533905932737*(alphaR[0]*fUpR[7]+alphaR[1]*fUpR[6]+alphaR[2]*fUpR[5]+fUpR[3]*alphaR[4]); 
  GhatR[8] = 0.3535533905932737*(alphaR[1]*fUpR[9]+alphaR[0]*fUpR[8])+0.3162277660168379*(alphaR[4]*fUpR[4]+alphaR[2]*fUpR[2]); 
  GhatR[9] = 0.3535533905932737*(alphaR[0]*fUpR[9]+alphaR[1]*fUpR[8])+0.3162277660168379*(alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.3535533905932737*(alphaR[1]*fUpR[11]+alphaR[0]*fUpR[10])+0.3162277660168379*(alphaR[4]*fUpR[7]+alphaR[2]*fUpR[6]); 
  GhatR[11] = 0.3535533905932737*(alphaR[0]*fUpR[11]+alphaR[1]*fUpR[10])+0.3162277660168379*(alphaR[2]*fUpR[7]+alphaR[4]*fUpR[6]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdz2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdz2; 
  out[2] += -1.224744871391589*GhatR[0]*rdz2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdz2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdz2; 
  out[5] += -1.224744871391589*GhatR[1]*rdz2; 
  out[6] += -0.7071067811865475*GhatR[4]*rdz2; 
  out[7] += -1.224744871391589*GhatR[2]*rdz2; 
  out[8] += -0.7071067811865475*GhatR[5]*rdz2; 
  out[9] += -1.224744871391589*GhatR[3]*rdz2; 
  out[10] += -0.7071067811865475*GhatR[6]*rdz2; 
  out[11] += -1.224744871391589*GhatR[4]*rdz2; 
  out[12] += -1.224744871391589*GhatR[5]*rdz2; 
  out[13] += -0.7071067811865475*GhatR[7]*rdz2; 
  out[14] += -1.224744871391589*GhatR[6]*rdz2; 
  out[15] += -1.224744871391589*GhatR[7]*rdz2; 
  out[16] += -0.7071067811865475*GhatR[8]*rdz2; 
  out[17] += -0.7071067811865475*GhatR[9]*rdz2; 
  out[18] += -1.224744871391589*GhatR[8]*rdz2; 
  out[19] += -0.7071067811865475*GhatR[10]*rdz2; 
  out[20] += -1.224744871391589*GhatR[9]*rdz2; 
  out[21] += -0.7071067811865475*GhatR[11]*rdz2; 
  out[22] += -1.224744871391589*GhatR[10]*rdz2; 
  out[23] += -1.224744871391589*GhatR[11]*rdz2; 

  } else { 

  double fUpL[12] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  fUpL[1] = 1.224744871391589*fedge[5]+0.7071067811865475*fedge[1]; 
  fUpL[2] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[3]; 
  fUpL[3] = 1.224744871391589*fedge[9]+0.7071067811865475*fedge[4]; 
  fUpL[4] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[6]; 
  fUpL[5] = 1.224744871391589*fedge[12]+0.7071067811865475*fedge[8]; 
  fUpL[6] = 1.224744871391589*fedge[14]+0.7071067811865475*fedge[10]; 
  fUpL[7] = 1.224744871391589*fedge[15]+0.7071067811865475*fedge[13]; 
  fUpL[8] = 1.224744871391589*fedge[18]+0.7071067811865475*fedge[16]; 
  fUpL[9] = 1.224744871391589*fedge[20]+0.7071067811865475*fedge[17]; 
  fUpL[10] = 1.224744871391589*fedge[22]+0.7071067811865475*fedge[19]; 
  fUpL[11] = 1.224744871391589*fedge[23]+0.7071067811865475*fedge[21]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[2]; 
  fUpL[1] = 0.7071067811865475*fskin[1]-1.224744871391589*fskin[5]; 
  fUpL[2] = 0.7071067811865475*fskin[3]-1.224744871391589*fskin[7]; 
  fUpL[3] = 0.7071067811865475*fskin[4]-1.224744871391589*fskin[9]; 
  fUpL[4] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[11]; 
  fUpL[5] = 0.7071067811865475*fskin[8]-1.224744871391589*fskin[12]; 
  fUpL[6] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[14]; 
  fUpL[7] = 0.7071067811865475*fskin[13]-1.224744871391589*fskin[15]; 
  fUpL[8] = 0.7071067811865475*fskin[16]-1.224744871391589*fskin[18]; 
  fUpL[9] = 0.7071067811865475*fskin[17]-1.224744871391589*fskin[20]; 
  fUpL[10] = 0.7071067811865475*fskin[19]-1.224744871391589*fskin[22]; 
  fUpL[11] = 0.7071067811865475*fskin[21]-1.224744871391589*fskin[23]; 
    } 
  } else { 
  double f_lr[12] = {0.};
  double f_cl[12] = {0.};
  double sgn_alphaUpL[6] = {0.};
  sgn_alphaUpL[0] = 0.2777777777777778*sgn_alpha_surfL[5]+0.4444444444444444*sgn_alpha_surfL[4]+0.2777777777777778*sgn_alpha_surfL[3]+0.2777777777777778*sgn_alpha_surfL[2]+0.4444444444444444*sgn_alpha_surfL[1]+0.2777777777777778*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[1] = 0.2777777777777778*sgn_alpha_surfL[5]+0.4444444444444444*sgn_alpha_surfL[4]+0.2777777777777778*sgn_alpha_surfL[3]-0.2777777777777778*sgn_alpha_surfL[2]-0.4444444444444444*sgn_alpha_surfL[1]-0.2777777777777778*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[2] = 0.3726779962499649*sgn_alpha_surfL[5]-0.3726779962499649*sgn_alpha_surfL[3]+0.3726779962499649*sgn_alpha_surfL[2]-0.3726779962499649*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[3] = 0.3726779962499649*sgn_alpha_surfL[5]-0.3726779962499649*sgn_alpha_surfL[3]-0.3726779962499649*sgn_alpha_surfL[2]+0.3726779962499649*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[4] = 0.2484519974999766*sgn_alpha_surfL[5]-0.4969039949999532*sgn_alpha_surfL[4]+0.2484519974999766*sgn_alpha_surfL[3]+0.2484519974999766*sgn_alpha_surfL[2]-0.4969039949999532*sgn_alpha_surfL[1]+0.2484519974999766*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[5] = 0.2484519974999767*sgn_alpha_surfL[5]-0.4969039949999535*sgn_alpha_surfL[4]+0.2484519974999767*sgn_alpha_surfL[3]-0.2484519974999767*sgn_alpha_surfL[2]+0.4969039949999535*sgn_alpha_surfL[1]-0.2484519974999767*sgn_alpha_surfL[0]; 

  f_lr[0] = 1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  f_lr[1] = 1.224744871391589*fedge[5]+0.7071067811865475*fedge[1]; 
  f_lr[2] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[3]; 
  f_lr[3] = 1.224744871391589*fedge[9]+0.7071067811865475*fedge[4]; 
  f_lr[4] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[6]; 
  f_lr[5] = 1.224744871391589*fedge[12]+0.7071067811865475*fedge[8]; 
  f_lr[6] = 1.224744871391589*fedge[14]+0.7071067811865475*fedge[10]; 
  f_lr[7] = 1.224744871391589*fedge[15]+0.7071067811865475*fedge[13]; 
  f_lr[8] = 1.224744871391589*fedge[18]+0.7071067811865475*fedge[16]; 
  f_lr[9] = 1.224744871391589*fedge[20]+0.7071067811865475*fedge[17]; 
  f_lr[10] = 1.224744871391589*fedge[22]+0.7071067811865475*fedge[19]; 
  f_lr[11] = 1.224744871391589*fedge[23]+0.7071067811865475*fedge[21]; 

  f_cl[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[2]; 
  f_cl[1] = 0.7071067811865475*fskin[1]-1.224744871391589*fskin[5]; 
  f_cl[2] = 0.7071067811865475*fskin[3]-1.224744871391589*fskin[7]; 
  f_cl[3] = 0.7071067811865475*fskin[4]-1.224744871391589*fskin[9]; 
  f_cl[4] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[11]; 
  f_cl[5] = 0.7071067811865475*fskin[8]-1.224744871391589*fskin[12]; 
  f_cl[6] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[14]; 
  f_cl[7] = 0.7071067811865475*fskin[13]-1.224744871391589*fskin[15]; 
  f_cl[8] = 0.7071067811865475*fskin[16]-1.224744871391589*fskin[18]; 
  f_cl[9] = 0.7071067811865475*fskin[17]-1.224744871391589*fskin[20]; 
  f_cl[10] = 0.7071067811865475*fskin[19]-1.224744871391589*fskin[22]; 
  f_cl[11] = 0.7071067811865475*fskin[21]-1.224744871391589*fskin[23]; 

  fUpL[0] = sgn_alphaUpL[5]*(0.25*f_lr[9]-0.25*f_cl[9])+sgn_alphaUpL[4]*(0.25*f_lr[8]-0.25*f_cl[8])+sgn_alphaUpL[3]*(0.25*f_lr[4]-0.25*f_cl[4])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[2]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[1]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = sgn_alphaUpL[4]*(0.2500000000000001*f_lr[9]-0.2500000000000001*f_cl[9])+sgn_alphaUpL[5]*(0.2500000000000001*f_lr[8]-0.2500000000000001*f_cl[8])+sgn_alphaUpL[2]*(0.25*f_lr[4]-0.25*f_cl[4])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[3]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[1]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[9]-0.223606797749979*f_cl[9])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[8]-0.223606797749979*f_cl[8])+(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])*sgn_alphaUpL[5]+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.25*f_lr[4]-0.25*f_cl[4])+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[3]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[2]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = sgn_alphaUpL[5]*(0.2500000000000001*f_lr[11]-0.2500000000000001*f_cl[11])+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[10]-0.2500000000000001*f_cl[10])+sgn_alphaUpL[3]*(0.25*f_lr[7]-0.25*f_cl[7])+sgn_alphaUpL[2]*(0.25*f_lr[6]-0.25*f_cl[6])+sgn_alphaUpL[1]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[9]-0.223606797749979*f_cl[9])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[8]-0.223606797749979*f_cl[8])+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[5]+(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])*sgn_alphaUpL[4]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[4]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[3]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.25*f_lr[2]-0.25*f_cl[2]); 
  fUpL[5] = sgn_alphaUpL[4]*(0.25*f_lr[11]-0.25*f_cl[11])+sgn_alphaUpL[5]*(0.25*f_lr[10]-0.25*f_cl[10])+sgn_alphaUpL[2]*(0.25*f_lr[7]-0.25*f_cl[7])+sgn_alphaUpL[3]*(0.25*f_lr[6]-0.25*f_cl[6])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[5]+sgn_alphaUpL[1]*(0.25*f_lr[3]-0.25*f_cl[3]); 
  fUpL[6] = sgn_alphaUpL[3]*(0.223606797749979*f_lr[11]-0.223606797749979*f_cl[11])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[10]-0.223606797749979*f_cl[10])+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[7]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[7]+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[6]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[6]+sgn_alphaUpL[3]*(0.25*f_lr[5]-0.25*f_cl[5])+sgn_alphaUpL[2]*(0.25*f_lr[3]-0.25*f_cl[3]); 
  fUpL[7] = sgn_alphaUpL[2]*(0.223606797749979*f_lr[11]-0.223606797749979*f_cl[11])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[10]-0.223606797749979*f_cl[10])+(0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[7]+((-0.223606797749979*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[7]+(0.223606797749979*sgn_alphaUpL[5]+0.25*sgn_alphaUpL[1])*f_lr[6]+((-0.223606797749979*sgn_alphaUpL[5])-0.25*sgn_alphaUpL[1])*f_cl[6]+sgn_alphaUpL[2]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.25*f_lr[3]-0.25*f_cl[3])*sgn_alphaUpL[3]; 
  fUpL[8] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[9]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[9]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[8]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[8]+(0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[5]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[2]; 
  fUpL[9] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[9]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[9]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[8]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[8]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[5]+(0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[3]; 
  fUpL[10] = (0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[11]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[11]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[10]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[10]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+sgn_alphaUpL[2]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.25*f_lr[5]-0.25*f_cl[5])*sgn_alphaUpL[5]+(0.2500000000000001*f_lr[3]-0.2500000000000001*f_cl[3])*sgn_alphaUpL[4]; 
  fUpL[11] = (0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[11]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[11]+(0.159719141249985*sgn_alphaUpL[5]+0.2500000000000001*sgn_alphaUpL[1])*f_lr[10]+((-0.159719141249985*sgn_alphaUpL[5])-0.2500000000000001*sgn_alphaUpL[1])*f_cl[10]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+sgn_alphaUpL[3]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.2500000000000001*f_lr[3]-0.2500000000000001*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.25*f_lr[5]-0.25*f_cl[5]); 

  } 
  double GhatL[12] = {0.};
  GhatL[0] = 0.3535533905932737*(alphaL[4]*fUpL[4]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.3535533905932737*(alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.3162277660168379*alphaL[4]*fUpL[9]+0.3162277660168379*alphaL[2]*fUpL[8]+0.3535533905932737*(alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.3535533905932737*(alphaL[4]*fUpL[7]+alphaL[2]*fUpL[6]+alphaL[1]*fUpL[5]+alphaL[0]*fUpL[3]); 
  GhatL[4] = 0.3162277660168379*alphaL[2]*fUpL[9]+0.3162277660168379*alphaL[4]*fUpL[8]+0.3535533905932737*(alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[5] = 0.3535533905932737*(alphaL[2]*fUpL[7]+alphaL[4]*fUpL[6]+alphaL[0]*fUpL[5]+alphaL[1]*fUpL[3]); 
  GhatL[6] = 0.3162277660168379*alphaL[4]*fUpL[11]+0.3162277660168379*alphaL[2]*fUpL[10]+0.3535533905932737*(alphaL[1]*fUpL[7]+alphaL[0]*fUpL[6]+alphaL[4]*fUpL[5]+alphaL[2]*fUpL[3]); 
  GhatL[7] = 0.3162277660168379*alphaL[2]*fUpL[11]+0.3162277660168379*alphaL[4]*fUpL[10]+0.3535533905932737*(alphaL[0]*fUpL[7]+alphaL[1]*fUpL[6]+alphaL[2]*fUpL[5]+fUpL[3]*alphaL[4]); 
  GhatL[8] = 0.3535533905932737*(alphaL[1]*fUpL[9]+alphaL[0]*fUpL[8])+0.3162277660168379*(alphaL[4]*fUpL[4]+alphaL[2]*fUpL[2]); 
  GhatL[9] = 0.3535533905932737*(alphaL[0]*fUpL[9]+alphaL[1]*fUpL[8])+0.3162277660168379*(alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.3535533905932737*(alphaL[1]*fUpL[11]+alphaL[0]*fUpL[10])+0.3162277660168379*(alphaL[4]*fUpL[7]+alphaL[2]*fUpL[6]); 
  GhatL[11] = 0.3535533905932737*(alphaL[0]*fUpL[11]+alphaL[1]*fUpL[10])+0.3162277660168379*(alphaL[2]*fUpL[7]+alphaL[4]*fUpL[6]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdz2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdz2; 
  out[2] += -1.224744871391589*GhatL[0]*rdz2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdz2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdz2; 
  out[5] += -1.224744871391589*GhatL[1]*rdz2; 
  out[6] += 0.7071067811865475*GhatL[4]*rdz2; 
  out[7] += -1.224744871391589*GhatL[2]*rdz2; 
  out[8] += 0.7071067811865475*GhatL[5]*rdz2; 
  out[9] += -1.224744871391589*GhatL[3]*rdz2; 
  out[10] += 0.7071067811865475*GhatL[6]*rdz2; 
  out[11] += -1.224744871391589*GhatL[4]*rdz2; 
  out[12] += -1.224744871391589*GhatL[5]*rdz2; 
  out[13] += 0.7071067811865475*GhatL[7]*rdz2; 
  out[14] += -1.224744871391589*GhatL[6]*rdz2; 
  out[15] += -1.224744871391589*GhatL[7]*rdz2; 
  out[16] += 0.7071067811865475*GhatL[8]*rdz2; 
  out[17] += 0.7071067811865475*GhatL[9]*rdz2; 
  out[18] += -1.224744871391589*GhatL[8]*rdz2; 
  out[19] += 0.7071067811865475*GhatL[10]*rdz2; 
  out[20] += -1.224744871391589*GhatL[9]*rdz2; 
  out[21] += 0.7071067811865475*GhatL[11]*rdz2; 
  out[22] += -1.224744871391589*GhatL[10]*rdz2; 
  out[23] += -1.224744871391589*GhatL[11]*rdz2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.5303300858899105*rdz2*cflFreq; 

} 
