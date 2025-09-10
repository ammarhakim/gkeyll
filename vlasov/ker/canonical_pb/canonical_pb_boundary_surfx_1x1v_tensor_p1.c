#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_basis_tensor_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
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

  double rdx2 = 2.0/dxv[0];

  const double *alphaL = &alpha_surf_skin[0];
  const double *alphaR = &alpha_surf_edge[0];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[0];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[0];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[0];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[0];

  if (edge == -1) { 

  double fUpR[2] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  fUpR[1] = 1.224744871391589*fskin[3]+0.7071067811865475*fskin[2]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[1]; 
  fUpR[1] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[3]; 
    } 
  } else { 
  double f_cr[2] = {0.};
  double f_rl[2] = {0.};
  double sgn_alphaUpR[2] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p1_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.224744871391589*fskin[1]+0.7071067811865475*fskin[0]; 
  f_cr[1] = 1.224744871391589*fskin[3]+0.7071067811865475*fskin[2]; 

  f_rl[0] = 0.7071067811865475*fedge[0]-1.224744871391589*fedge[1]; 
  f_rl[1] = 0.7071067811865475*fedge[2]-1.224744871391589*fedge[3]; 

  fUpR[0] = (0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 

  } 
  double GhatR[2] = {0.};
  GhatR[0] = 0.7071067811865475*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.7071067811865475*alphaR[0]*fUpR[1]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -1.224744871391589*GhatR[1]*rdx2; 

  } else { 

  double fUpL[2] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  fUpL[1] = 1.224744871391589*fedge[3]+0.7071067811865475*fedge[2]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[1]; 
  fUpL[1] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[3]; 
    } 
  } else { 
  double f_lr[2] = {0.};
  double f_cl[2] = {0.};
  double sgn_alphaUpL[2] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p1_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.224744871391589*fedge[1]+0.7071067811865475*fedge[0]; 
  f_lr[1] = 1.224744871391589*fedge[3]+0.7071067811865475*fedge[2]; 

  f_cl[0] = 0.7071067811865475*fskin[0]-1.224744871391589*fskin[1]; 
  f_cl[1] = 0.7071067811865475*fskin[2]-1.224744871391589*fskin[3]; 

  fUpL[0] = (0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[1]+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 

  } 
  double GhatL[2] = {0.};
  GhatL[0] = 0.7071067811865475*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.7071067811865475*alphaL[0]*fUpL[1]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += -1.224744871391589*GhatL[1]*rdx2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 1.060660171779821*rdx2*cflFreq; 

} 
