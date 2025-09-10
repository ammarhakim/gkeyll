#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
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

  double rdx2 = 2.0/dxv[0];

  const double *alphaL = &alpha_surf_skin[0];
  const double *alphaR = &alpha_surf_edge[0];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[0];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[0];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[0];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[0];

  if (edge == -1) { 

  double fUpR[6] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  fUpR[1] = 1.224744871391589*fskin[4]+0.7071067811865475*fskin[2]; 
  fUpR[2] = 1.224744871391589*fskin[5]+0.7071067811865475*fskin[3]; 
  fUpR[3] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[6]; 
  fUpR[4] = 1.224744871391589*fskin[9]+0.7071067811865475*fskin[8]; 
  fUpR[5] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[10]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[1]; 
  fUpR[1] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[4]; 
  fUpR[2] = 0.7071067811865475*fedge[3]-1.224744871391589*fedge[5]; 
  fUpR[3] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[7]; 
  fUpR[4] = 0.7071067811865475*fedge[8]-1.224744871391589*fedge[9]; 
  fUpR[5] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[11]; 
    } 
  } else { 
  double f_cr[6] = {0.};
  double f_rl[6] = {0.};
  double sgn_alphaUpR[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  f_cr[1] = 1.224744871391589*fskin[4]+0.7071067811865475*fskin[2]; 
  f_cr[2] = 1.224744871391589*fskin[5]+0.7071067811865475*fskin[3]; 
  f_cr[3] = 1.224744871391589*fskin[7]+0.7071067811865475*fskin[6]; 
  f_cr[4] = 1.224744871391589*fskin[9]+0.7071067811865475*fskin[8]; 
  f_cr[5] = 1.224744871391589*fskin[11]+0.7071067811865475*fskin[10]; 

  f_rl[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[1]; 
  f_rl[1] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[4]; 
  f_rl[2] = 0.7071067811865475*fedge[3]-1.224744871391589*fedge[5]; 
  f_rl[3] = 0.7071067811865475*fedge[6]-1.224744871391589*fedge[7]; 
  f_rl[4] = 0.7071067811865475*fedge[8]-1.224744871391589*fedge[9]; 
  f_rl[5] = 0.7071067811865475*fedge[10]-1.224744871391589*fedge[11]; 

  fUpR[0] = (0.25*f_cr[5]-0.25*f_rl[5])*sgn_alphaUpR[5]+(0.25*f_cr[4]-0.25*f_rl[4])*sgn_alphaUpR[4]+(0.25*f_cr[3]-0.25*f_rl[3])*sgn_alphaUpR[3]+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[2]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[1]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[5]-0.223606797749979*f_rl[5])+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[1]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[5]-0.2500000000000001*f_rl[5])+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[2]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[5]-0.223606797749979*f_rl[5])+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*((-0.223606797749979*f_rl[4])+0.223606797749979*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[3]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[3]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.25*f_cr[2]-0.25*f_rl[2]); 
  fUpR[4] = ((-0.159719141249985*f_rl[5])+0.159719141249985*f_cr[5]-0.2500000000000001*f_rl[2]+0.2500000000000001*f_cr[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.2500000000000001*f_cr[5]-0.2500000000000001*f_rl[5])+((-0.159719141249985*f_rl[4])+0.159719141249985*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[4]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[4]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[4]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[3]+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[1]; 
  fUpR[5] = ((-0.159719141249985*f_rl[4])+0.159719141249985*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[5]+((-0.159719141249985*sgn_alphaUpR[4])-0.25*sgn_alphaUpR[0]+0.5)*f_rl[5]+(0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.2500000000000001*f_cr[2]-0.2500000000000001*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3]); 

  } 
  double GhatR[6] = {0.};
  GhatR[0] = 0.5*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.5*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.5*(alphaR[1]*fUpR[3]+alphaR[0]*fUpR[2]); 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[5]+0.5*(alphaR[0]*fUpR[3]+alphaR[1]*fUpR[2]); 
  GhatR[4] = 0.5*alphaR[0]*fUpR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.5*alphaR[0]*fUpR[5]+0.447213595499958*alphaR[1]*fUpR[3]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -1.224744871391589*GhatR[1]*rdx2; 
  out[5] += -1.224744871391589*GhatR[2]*rdx2; 
  out[6] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[7] += -1.224744871391589*GhatR[3]*rdx2; 
  out[8] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[9] += -1.224744871391589*GhatR[4]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[11] += -1.224744871391589*GhatR[5]*rdx2; 

  } else { 

  double fUpL[6] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  fUpL[1] = 1.224744871391589*fedge[4]+0.7071067811865475*fedge[2]; 
  fUpL[2] = 1.224744871391589*fedge[5]+0.7071067811865475*fedge[3]; 
  fUpL[3] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[6]; 
  fUpL[4] = 1.224744871391589*fedge[9]+0.7071067811865475*fedge[8]; 
  fUpL[5] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[10]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[1]; 
  fUpL[1] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[4]; 
  fUpL[2] = 0.7071067811865475*fskin[3]-1.224744871391589*fskin[5]; 
  fUpL[3] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[7]; 
  fUpL[4] = 0.7071067811865475*fskin[8]-1.224744871391589*fskin[9]; 
  fUpL[5] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[11]; 
    } 
  } else { 
  double f_lr[6] = {0.};
  double f_cl[6] = {0.};
  double sgn_alphaUpL[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  f_lr[1] = 1.224744871391589*fedge[4]+0.7071067811865475*fedge[2]; 
  f_lr[2] = 1.224744871391589*fedge[5]+0.7071067811865475*fedge[3]; 
  f_lr[3] = 1.224744871391589*fedge[7]+0.7071067811865475*fedge[6]; 
  f_lr[4] = 1.224744871391589*fedge[9]+0.7071067811865475*fedge[8]; 
  f_lr[5] = 1.224744871391589*fedge[11]+0.7071067811865475*fedge[10]; 

  f_cl[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[1]; 
  f_cl[1] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[4]; 
  f_cl[2] = 0.7071067811865475*fskin[3]-1.224744871391589*fskin[5]; 
  f_cl[3] = 0.7071067811865475*fskin[6]-1.224744871391589*fskin[7]; 
  f_cl[4] = 0.7071067811865475*fskin[8]-1.224744871391589*fskin[9]; 
  f_cl[5] = 0.7071067811865475*fskin[10]-1.224744871391589*fskin[11]; 

  fUpL[0] = (0.25*f_lr[5]-0.25*f_cl[5])*sgn_alphaUpL[5]+(0.25*f_lr[4]-0.25*f_cl[4])*sgn_alphaUpL[4]+(0.25*f_lr[3]-0.25*f_cl[3])*sgn_alphaUpL[3]+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[2]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[1]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[5]-0.223606797749979*f_cl[5])+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[1]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[5]-0.2500000000000001*f_cl[5])+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[2]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[5]-0.223606797749979*f_cl[5])+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[3]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.25*f_lr[2]-0.25*f_cl[2]); 
  fUpL[4] = (0.159719141249985*f_lr[5]-0.159719141249985*f_cl[5]+0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.2500000000000001*f_lr[5]-0.2500000000000001*f_cl[5])+(0.159719141249985*f_lr[4]-0.159719141249985*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[4]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[4]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[3]+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[1]; 
  fUpL[5] = (0.159719141249985*f_lr[4]-0.159719141249985*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[5]+(0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[5]+((-0.159719141249985*sgn_alphaUpL[4])-0.25*sgn_alphaUpL[0]+0.5)*f_cl[5]+(0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3]); 

  } 
  double GhatL[6] = {0.};
  GhatL[0] = 0.5*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.5*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.5*(alphaL[1]*fUpL[3]+alphaL[0]*fUpL[2]); 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[5]+0.5*(alphaL[0]*fUpL[3]+alphaL[1]*fUpL[2]); 
  GhatL[4] = 0.5*alphaL[0]*fUpL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.5*alphaL[0]*fUpL[5]+0.447213595499958*alphaL[1]*fUpL[3]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += -1.224744871391589*GhatL[1]*rdx2; 
  out[5] += -1.224744871391589*GhatL[2]*rdx2; 
  out[6] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[7] += -1.224744871391589*GhatL[3]*rdx2; 
  out[8] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[9] += -1.224744871391589*GhatL[4]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[11] += -1.224744871391589*GhatL[5]*rdx2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.75*rdx2*cflFreq; 

} 
