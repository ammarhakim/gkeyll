#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_basis_tensor_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points.
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points.
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1.
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdvx2 = 2.0/dxv[2];

  const double *alphaL = &alpha_surf_skin[32];
  const double *alphaR = &alpha_surf_edge[32];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[32];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[32];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[2];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[2];

  if (edge == -1) { 

  double fUpR[16] = {0.};
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
    } 
  } else { 
  double f_cr[16] = {0.};
  double f_rl[16] = {0.};
  double sgn_alphaUpR[16] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

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

  fUpR[0] = (0.125*f_cr[15]-0.125*f_rl[15])*sgn_alphaUpR[15]+(0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[14]+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[13]+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[12]+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[11]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[10]+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[9]+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[8]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[7]+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[6]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[5]+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[4]+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[3]+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[2]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[1]+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[15]+sgn_alphaUpR[14]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[13]+sgn_alphaUpR[10]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[12]+sgn_alphaUpR[9]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[11]+sgn_alphaUpR[7]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[8]+sgn_alphaUpR[4]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[1]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[15]+sgn_alphaUpR[13]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[14]+sgn_alphaUpR[10]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[12]+sgn_alphaUpR[8]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[11]+sgn_alphaUpR[6]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[9]+sgn_alphaUpR[4]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[2]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[15]+sgn_alphaUpR[12]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[14]+sgn_alphaUpR[9]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[13]+sgn_alphaUpR[8]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[11]+sgn_alphaUpR[5]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[10]+sgn_alphaUpR[4]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[7]+sgn_alphaUpR[2]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[6]+sgn_alphaUpR[1]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[3]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[3]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = (0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[15]+sgn_alphaUpR[11]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[14]+sgn_alphaUpR[7]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[13]+sgn_alphaUpR[6]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[12]+sgn_alphaUpR[5]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[10]+sgn_alphaUpR[3]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[9]+sgn_alphaUpR[2]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[8]+sgn_alphaUpR[1]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[4]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[4]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[4]; 
  fUpR[5] = (0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[15]+sgn_alphaUpR[10]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[14]+sgn_alphaUpR[13]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[12]+sgn_alphaUpR[4]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[11]+sgn_alphaUpR[3]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[9]+sgn_alphaUpR[8]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[7]+sgn_alphaUpR[6]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[5]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[5]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.125*f_cr[2]-0.125*f_rl[2]); 
  fUpR[6] = (0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[15]+sgn_alphaUpR[9]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[14]+sgn_alphaUpR[12]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[13]+sgn_alphaUpR[4]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[11]+sgn_alphaUpR[2]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[10]+sgn_alphaUpR[8]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[7]+sgn_alphaUpR[5]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[6]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[6]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[7] = (0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[15]+sgn_alphaUpR[8]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[14]+sgn_alphaUpR[4]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[13]+sgn_alphaUpR[12]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[11]+sgn_alphaUpR[1]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[10]+sgn_alphaUpR[9]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[7]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[7]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[7]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[6]+sgn_alphaUpR[5]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[8] = (0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[15]+sgn_alphaUpR[7]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[14]+sgn_alphaUpR[11]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[13]+sgn_alphaUpR[3]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[12]+sgn_alphaUpR[2]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[10]+sgn_alphaUpR[6]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[9]+sgn_alphaUpR[5]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[8]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[8]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[8]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[9] = (0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[15]+sgn_alphaUpR[6]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[14]+sgn_alphaUpR[3]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[13]+sgn_alphaUpR[11]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[12]+sgn_alphaUpR[1]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[10]+sgn_alphaUpR[7]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[9]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[9]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[9]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[8]+sgn_alphaUpR[5]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[10] = (0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[15]+sgn_alphaUpR[5]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[14]+sgn_alphaUpR[2]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[13]+sgn_alphaUpR[1]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[12]+sgn_alphaUpR[11]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[10]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[10]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[10]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[9]+sgn_alphaUpR[7]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[8]+sgn_alphaUpR[6]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[11] = (0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[15]+sgn_alphaUpR[4]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[14]+sgn_alphaUpR[8]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[13]+sgn_alphaUpR[9]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[12]+sgn_alphaUpR[10]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[11]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[11]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[11]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[12] = (0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[15]+sgn_alphaUpR[3]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[14]+sgn_alphaUpR[6]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[13]+sgn_alphaUpR[7]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[12]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[12]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[12]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[11]+sgn_alphaUpR[10]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[9]+sgn_alphaUpR[1]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[8]+sgn_alphaUpR[2]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[13] = (0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[15]+sgn_alphaUpR[2]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[14]+sgn_alphaUpR[5]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[13]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[13]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[13]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[12]+sgn_alphaUpR[7]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[11]+sgn_alphaUpR[9]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[10]+sgn_alphaUpR[1]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[8]+sgn_alphaUpR[3]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[6]+sgn_alphaUpR[4]*(0.125*f_cr[6]-0.125*f_rl[6]); 
  fUpR[14] = (0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[15]+sgn_alphaUpR[1]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[14]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[14]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[14]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[13]+sgn_alphaUpR[5]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[12]+sgn_alphaUpR[6]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[11]+sgn_alphaUpR[8]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[10]+sgn_alphaUpR[2]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[9]+sgn_alphaUpR[3]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[7]+sgn_alphaUpR[4]*(0.125*f_cr[7]-0.125*f_rl[7]); 
  fUpR[15] = (0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[15]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[15]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[15]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[14]+sgn_alphaUpR[1]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[13]+sgn_alphaUpR[2]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[12]+sgn_alphaUpR[3]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[11]+sgn_alphaUpR[4]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[10]+sgn_alphaUpR[5]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[9]+sgn_alphaUpR[6]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[8]+sgn_alphaUpR[7]*(0.125*f_cr[8]-0.125*f_rl[8]); 

  } 
  double GhatR[16] = {0.};
  GhatR[0] = 0.25*(alphaR[14]*fUpR[14]+alphaR[10]*fUpR[10]+alphaR[9]*fUpR[9]+alphaR[7]*fUpR[7]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[14]*fUpR[15]+alphaR[10]*fUpR[13]+alphaR[9]*fUpR[12]+alphaR[7]*fUpR[11]+alphaR[4]*fUpR[8]+alphaR[3]*fUpR[6]+alphaR[2]*fUpR[5]+alphaR[0]*fUpR[1]); 
  GhatR[2] = 0.25*(alphaR[10]*fUpR[14]+fUpR[10]*alphaR[14]+alphaR[4]*fUpR[9]+fUpR[4]*alphaR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.25*(alphaR[9]*fUpR[14]+fUpR[9]*alphaR[14]+alphaR[4]*fUpR[10]+fUpR[4]*alphaR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.25*(alphaR[7]*fUpR[14]+fUpR[7]*alphaR[14]+alphaR[3]*fUpR[10]+fUpR[3]*alphaR[10]+alphaR[2]*fUpR[9]+fUpR[2]*alphaR[9]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[10]*fUpR[15]+fUpR[13]*alphaR[14]+alphaR[4]*fUpR[12]+alphaR[3]*fUpR[11]+fUpR[8]*alphaR[9]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.25*(alphaR[9]*fUpR[15]+fUpR[12]*alphaR[14]+alphaR[4]*fUpR[13]+alphaR[2]*fUpR[11]+fUpR[8]*alphaR[10]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.25*(alphaR[4]*fUpR[14]+fUpR[4]*alphaR[14]+alphaR[9]*fUpR[10]+fUpR[9]*alphaR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[7]*fUpR[15]+fUpR[11]*alphaR[14]+alphaR[3]*fUpR[13]+alphaR[2]*fUpR[12]+fUpR[6]*alphaR[10]+fUpR[5]*alphaR[9]+alphaR[0]*fUpR[8]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[3]*fUpR[14]+fUpR[3]*alphaR[14]+alphaR[7]*fUpR[10]+fUpR[7]*alphaR[10]+alphaR[0]*fUpR[9]+fUpR[0]*alphaR[9]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.25*(alphaR[2]*fUpR[14]+fUpR[2]*alphaR[14]+alphaR[0]*fUpR[10]+fUpR[0]*alphaR[10]+alphaR[7]*fUpR[9]+fUpR[7]*alphaR[9]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.25*(alphaR[4]*fUpR[15]+fUpR[8]*alphaR[14]+alphaR[9]*fUpR[13]+alphaR[10]*fUpR[12]+alphaR[0]*fUpR[11]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+alphaR[3]*fUpR[5]); 
  GhatR[12] = 0.25*(alphaR[3]*fUpR[15]+fUpR[6]*alphaR[14]+alphaR[7]*fUpR[13]+alphaR[0]*fUpR[12]+alphaR[10]*fUpR[11]+fUpR[1]*alphaR[9]+alphaR[2]*fUpR[8]+alphaR[4]*fUpR[5]); 
  GhatR[13] = 0.25*(alphaR[2]*fUpR[15]+fUpR[5]*alphaR[14]+alphaR[0]*fUpR[13]+alphaR[7]*fUpR[12]+alphaR[9]*fUpR[11]+fUpR[1]*alphaR[10]+alphaR[3]*fUpR[8]+alphaR[4]*fUpR[6]); 
  GhatR[14] = 0.25*(alphaR[0]*fUpR[14]+fUpR[0]*alphaR[14]+alphaR[2]*fUpR[10]+fUpR[2]*alphaR[10]+alphaR[3]*fUpR[9]+fUpR[3]*alphaR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[15] = 0.25*(alphaR[0]*fUpR[15]+fUpR[1]*alphaR[14]+alphaR[2]*fUpR[13]+alphaR[3]*fUpR[12]+alphaR[4]*fUpR[11]+fUpR[5]*alphaR[10]+fUpR[6]*alphaR[9]+alphaR[7]*fUpR[8]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdvx2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvx2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvx2; 
  out[3] += -1.224744871391589*GhatR[0]*rdvx2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdvx2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvx2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdvx2; 
  out[7] += -1.224744871391589*GhatR[1]*rdvx2; 
  out[8] += -1.224744871391589*GhatR[2]*rdvx2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdvx2; 
  out[10] += -0.7071067811865475*GhatR[7]*rdvx2; 
  out[11] += -1.224744871391589*GhatR[3]*rdvx2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdvx2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdvx2; 
  out[14] += -1.224744871391589*GhatR[4]*rdvx2; 
  out[15] += -0.7071067811865475*GhatR[10]*rdvx2; 
  out[16] += -1.224744871391589*GhatR[5]*rdvx2; 
  out[17] += -0.7071067811865475*GhatR[11]*rdvx2; 
  out[18] += -1.224744871391589*GhatR[6]*rdvx2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvx2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdvx2; 
  out[21] += -1.224744871391589*GhatR[8]*rdvx2; 
  out[22] += -1.224744871391589*GhatR[9]*rdvx2; 
  out[23] += -0.7071067811865475*GhatR[13]*rdvx2; 
  out[24] += -0.7071067811865475*GhatR[14]*rdvx2; 
  out[25] += -1.224744871391589*GhatR[10]*rdvx2; 
  out[26] += -1.224744871391589*GhatR[11]*rdvx2; 
  out[27] += -1.224744871391589*GhatR[12]*rdvx2; 
  out[28] += -0.7071067811865475*GhatR[15]*rdvx2; 
  out[29] += -1.224744871391589*GhatR[13]*rdvx2; 
  out[30] += -1.224744871391589*GhatR[14]*rdvx2; 
  out[31] += -1.224744871391589*GhatR[15]*rdvx2; 

  } else { 

  double fUpL[16] = {0.};
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
    } 
  } else { 
  double f_lr[16] = {0.};
  double f_cl[16] = {0.};
  double sgn_alphaUpL[16] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

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

  fUpL[0] = (0.125*f_lr[15]-0.125*f_cl[15])*sgn_alphaUpL[15]+(0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[14]+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[13]+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[12]+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[11]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[10]+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[9]+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[8]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[7]+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[6]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[5]+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[4]+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[3]+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[2]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[1]+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[15]+sgn_alphaUpL[14]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[13]+sgn_alphaUpL[10]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[12]+sgn_alphaUpL[9]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[11]+sgn_alphaUpL[7]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[8]+sgn_alphaUpL[4]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[1]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[15]+sgn_alphaUpL[13]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[14]+sgn_alphaUpL[10]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[12]+sgn_alphaUpL[8]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[11]+sgn_alphaUpL[6]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[9]+sgn_alphaUpL[4]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[2]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[15]+sgn_alphaUpL[12]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[14]+sgn_alphaUpL[9]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[13]+sgn_alphaUpL[8]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[11]+sgn_alphaUpL[5]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[10]+sgn_alphaUpL[4]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[7]+sgn_alphaUpL[2]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[6]+sgn_alphaUpL[1]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[3]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = (0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[15]+sgn_alphaUpL[11]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[14]+sgn_alphaUpL[7]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[13]+sgn_alphaUpL[6]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[12]+sgn_alphaUpL[5]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[10]+sgn_alphaUpL[3]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[9]+sgn_alphaUpL[2]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[8]+sgn_alphaUpL[1]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[4]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[4]; 
  fUpL[5] = (0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[15]+sgn_alphaUpL[10]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[14]+sgn_alphaUpL[13]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[12]+sgn_alphaUpL[4]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[11]+sgn_alphaUpL[3]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[9]+sgn_alphaUpL[8]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[7]+sgn_alphaUpL[6]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[5]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[5]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.125*f_lr[2]-0.125*f_cl[2]); 
  fUpL[6] = (0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[15]+sgn_alphaUpL[9]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[14]+sgn_alphaUpL[12]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[13]+sgn_alphaUpL[4]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[11]+sgn_alphaUpL[2]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[10]+sgn_alphaUpL[8]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[7]+sgn_alphaUpL[5]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[6]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[6]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[6]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[7] = (0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[15]+sgn_alphaUpL[8]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[14]+sgn_alphaUpL[4]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[13]+sgn_alphaUpL[12]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[11]+sgn_alphaUpL[1]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[10]+sgn_alphaUpL[9]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[7]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[7]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[7]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[6]+sgn_alphaUpL[5]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[8] = (0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[15]+sgn_alphaUpL[7]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[14]+sgn_alphaUpL[11]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[13]+sgn_alphaUpL[3]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[12]+sgn_alphaUpL[2]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[10]+sgn_alphaUpL[6]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[9]+sgn_alphaUpL[5]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[8]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[8]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[8]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[9] = (0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[15]+sgn_alphaUpL[6]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[14]+sgn_alphaUpL[3]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[13]+sgn_alphaUpL[11]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[12]+sgn_alphaUpL[1]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[10]+sgn_alphaUpL[7]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[9]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[9]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[9]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[8]+sgn_alphaUpL[5]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[10] = (0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[15]+sgn_alphaUpL[5]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[14]+sgn_alphaUpL[2]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[13]+sgn_alphaUpL[1]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[12]+sgn_alphaUpL[11]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[10]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[10]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[10]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[9]+sgn_alphaUpL[7]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[8]+sgn_alphaUpL[6]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[11] = (0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[15]+sgn_alphaUpL[4]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[14]+sgn_alphaUpL[8]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[13]+sgn_alphaUpL[9]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[12]+sgn_alphaUpL[10]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[11]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[11]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[11]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[12] = (0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[15]+sgn_alphaUpL[3]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[14]+sgn_alphaUpL[6]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[13]+sgn_alphaUpL[7]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[12]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[12]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[12]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[11]+sgn_alphaUpL[10]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[9]+sgn_alphaUpL[1]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[8]+sgn_alphaUpL[2]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[13] = (0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[15]+sgn_alphaUpL[2]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[14]+sgn_alphaUpL[5]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[13]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[13]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[13]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[12]+sgn_alphaUpL[7]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[11]+sgn_alphaUpL[9]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[10]+sgn_alphaUpL[1]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[8]+sgn_alphaUpL[3]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[6]+sgn_alphaUpL[4]*(0.125*f_lr[6]-0.125*f_cl[6]); 
  fUpL[14] = (0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[15]+sgn_alphaUpL[1]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[14]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[14]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[14]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[13]+sgn_alphaUpL[5]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[12]+sgn_alphaUpL[6]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[11]+sgn_alphaUpL[8]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[10]+sgn_alphaUpL[2]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[9]+sgn_alphaUpL[3]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[7]+sgn_alphaUpL[4]*(0.125*f_lr[7]-0.125*f_cl[7]); 
  fUpL[15] = (0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[15]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[15]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[15]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[14]+sgn_alphaUpL[1]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[13]+sgn_alphaUpL[2]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[12]+sgn_alphaUpL[3]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[11]+sgn_alphaUpL[4]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[10]+sgn_alphaUpL[5]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[9]+sgn_alphaUpL[6]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[8]+sgn_alphaUpL[7]*(0.125*f_lr[8]-0.125*f_cl[8]); 

  } 
  double GhatL[16] = {0.};
  GhatL[0] = 0.25*(alphaL[14]*fUpL[14]+alphaL[10]*fUpL[10]+alphaL[9]*fUpL[9]+alphaL[7]*fUpL[7]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[14]*fUpL[15]+alphaL[10]*fUpL[13]+alphaL[9]*fUpL[12]+alphaL[7]*fUpL[11]+alphaL[4]*fUpL[8]+alphaL[3]*fUpL[6]+alphaL[2]*fUpL[5]+alphaL[0]*fUpL[1]); 
  GhatL[2] = 0.25*(alphaL[10]*fUpL[14]+fUpL[10]*alphaL[14]+alphaL[4]*fUpL[9]+fUpL[4]*alphaL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.25*(alphaL[9]*fUpL[14]+fUpL[9]*alphaL[14]+alphaL[4]*fUpL[10]+fUpL[4]*alphaL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.25*(alphaL[7]*fUpL[14]+fUpL[7]*alphaL[14]+alphaL[3]*fUpL[10]+fUpL[3]*alphaL[10]+alphaL[2]*fUpL[9]+fUpL[2]*alphaL[9]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[10]*fUpL[15]+fUpL[13]*alphaL[14]+alphaL[4]*fUpL[12]+alphaL[3]*fUpL[11]+fUpL[8]*alphaL[9]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.25*(alphaL[9]*fUpL[15]+fUpL[12]*alphaL[14]+alphaL[4]*fUpL[13]+alphaL[2]*fUpL[11]+fUpL[8]*alphaL[10]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.25*(alphaL[4]*fUpL[14]+fUpL[4]*alphaL[14]+alphaL[9]*fUpL[10]+fUpL[9]*alphaL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[7]*fUpL[15]+fUpL[11]*alphaL[14]+alphaL[3]*fUpL[13]+alphaL[2]*fUpL[12]+fUpL[6]*alphaL[10]+fUpL[5]*alphaL[9]+alphaL[0]*fUpL[8]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[3]*fUpL[14]+fUpL[3]*alphaL[14]+alphaL[7]*fUpL[10]+fUpL[7]*alphaL[10]+alphaL[0]*fUpL[9]+fUpL[0]*alphaL[9]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.25*(alphaL[2]*fUpL[14]+fUpL[2]*alphaL[14]+alphaL[0]*fUpL[10]+fUpL[0]*alphaL[10]+alphaL[7]*fUpL[9]+fUpL[7]*alphaL[9]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.25*(alphaL[4]*fUpL[15]+fUpL[8]*alphaL[14]+alphaL[9]*fUpL[13]+alphaL[10]*fUpL[12]+alphaL[0]*fUpL[11]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+alphaL[3]*fUpL[5]); 
  GhatL[12] = 0.25*(alphaL[3]*fUpL[15]+fUpL[6]*alphaL[14]+alphaL[7]*fUpL[13]+alphaL[0]*fUpL[12]+alphaL[10]*fUpL[11]+fUpL[1]*alphaL[9]+alphaL[2]*fUpL[8]+alphaL[4]*fUpL[5]); 
  GhatL[13] = 0.25*(alphaL[2]*fUpL[15]+fUpL[5]*alphaL[14]+alphaL[0]*fUpL[13]+alphaL[7]*fUpL[12]+alphaL[9]*fUpL[11]+fUpL[1]*alphaL[10]+alphaL[3]*fUpL[8]+alphaL[4]*fUpL[6]); 
  GhatL[14] = 0.25*(alphaL[0]*fUpL[14]+fUpL[0]*alphaL[14]+alphaL[2]*fUpL[10]+fUpL[2]*alphaL[10]+alphaL[3]*fUpL[9]+fUpL[3]*alphaL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[15] = 0.25*(alphaL[0]*fUpL[15]+fUpL[1]*alphaL[14]+alphaL[2]*fUpL[13]+alphaL[3]*fUpL[12]+alphaL[4]*fUpL[11]+fUpL[5]*alphaL[10]+fUpL[6]*alphaL[9]+alphaL[7]*fUpL[8]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdvx2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvx2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvx2; 
  out[3] += -1.224744871391589*GhatL[0]*rdvx2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdvx2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvx2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdvx2; 
  out[7] += -1.224744871391589*GhatL[1]*rdvx2; 
  out[8] += -1.224744871391589*GhatL[2]*rdvx2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdvx2; 
  out[10] += 0.7071067811865475*GhatL[7]*rdvx2; 
  out[11] += -1.224744871391589*GhatL[3]*rdvx2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdvx2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdvx2; 
  out[14] += -1.224744871391589*GhatL[4]*rdvx2; 
  out[15] += 0.7071067811865475*GhatL[10]*rdvx2; 
  out[16] += -1.224744871391589*GhatL[5]*rdvx2; 
  out[17] += 0.7071067811865475*GhatL[11]*rdvx2; 
  out[18] += -1.224744871391589*GhatL[6]*rdvx2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvx2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdvx2; 
  out[21] += -1.224744871391589*GhatL[8]*rdvx2; 
  out[22] += -1.224744871391589*GhatL[9]*rdvx2; 
  out[23] += 0.7071067811865475*GhatL[13]*rdvx2; 
  out[24] += 0.7071067811865475*GhatL[14]*rdvx2; 
  out[25] += -1.224744871391589*GhatL[10]*rdvx2; 
  out[26] += -1.224744871391589*GhatL[11]*rdvx2; 
  out[27] += -1.224744871391589*GhatL[12]*rdvx2; 
  out[28] += 0.7071067811865475*GhatL[15]*rdvx2; 
  out[29] += -1.224744871391589*GhatL[13]*rdvx2; 
  out[30] += -1.224744871391589*GhatL[14]*rdvx2; 
  out[31] += -1.224744871391589*GhatL[15]*rdvx2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.375*rdvx2*cflFreq; 

} 
