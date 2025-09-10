#include <gkyl_advection_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double advection_surfy_2x_ser_p1(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // ul/uc/ur:  Advection velocity in left/center/right cells.
  // fl/fc/fr:  Input function in left/center/right cells.
  // out:       Incremented function in center cell.
  const double dx1 = 2.0/dxv[1]; 
  const double *ul_0 = &ul[4]; 
  const double *uc_0 = &uc[4]; 
  const double *ur_0 = &ur[4]; 
  double uQuad_l[2] = {0.0};
  double uQuad_r[2] = {0.0};
  double uMax_l[2] = {0.0};;
  double uMax_r[2] = {0.0};
  double Ghat_l[2] = {0.0}; 
  double Ghat_r[2] = {0.0}; 
  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(ul_0); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(uc_0); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(uc_0); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(ur_0); 
  uQuad_l[0] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[0] = fmax(fabs(u_c_r), fabs(u_r_l)); 
  u_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(ul_0); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(uc_0); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(uc_0); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(ur_0); 
  uQuad_l[1] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[1] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uQuad_l, uMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uQuad_r, uMax_r); 
  Ghat_l[0] = 0.5303300858899105*fl[3]*ul_0[3]+0.3061862178478971*fl[1]*ul_0[3]+0.5303300858899105*fc[3]*uc_0[3]-0.3061862178478971*fc[1]*uc_0[3]+0.3061862178478971*ul_0[1]*fl[3]+0.4330127018922193*uMax_l[1]*fl[3]-0.3061862178478971*uc_0[1]*fc[3]+0.4330127018922193*uMax_l[1]*fc[3]+0.5303300858899105*fl[2]*ul_0[2]+0.3061862178478971*fl[0]*ul_0[2]+0.5303300858899105*fc[2]*uc_0[2]-0.3061862178478971*fc[0]*uc_0[2]+0.3061862178478971*ul_0[0]*fl[2]+0.4330127018922193*uMax_l[0]*fl[2]-0.3061862178478971*uc_0[0]*fc[2]+0.4330127018922193*uMax_l[0]*fc[2]+0.1767766952966368*fl[1]*ul_0[1]+0.1767766952966368*fc[1]*uc_0[1]+0.25*fl[1]*uMax_l[1]-0.25*fc[1]*uMax_l[1]+0.1767766952966368*fl[0]*ul_0[0]+0.1767766952966368*fc[0]*uc_0[0]+0.25*fl[0]*uMax_l[0]-0.25*fc[0]*uMax_l[0]; 
  Ghat_l[1] = 0.5303300858899105*fl[2]*ul_0[3]+0.3061862178478971*fl[0]*ul_0[3]+0.5303300858899105*fc[2]*uc_0[3]-0.3061862178478971*fc[0]*uc_0[3]+0.5303300858899105*ul_0[2]*fl[3]+0.3061862178478971*ul_0[0]*fl[3]+0.4330127018922193*uMax_l[0]*fl[3]+0.5303300858899105*uc_0[2]*fc[3]-0.3061862178478971*uc_0[0]*fc[3]+0.4330127018922193*uMax_l[0]*fc[3]+0.3061862178478971*fl[1]*ul_0[2]-0.3061862178478971*fc[1]*uc_0[2]+0.3061862178478971*ul_0[1]*fl[2]+0.4330127018922193*uMax_l[1]*fl[2]-0.3061862178478971*uc_0[1]*fc[2]+0.4330127018922193*uMax_l[1]*fc[2]+0.1767766952966368*fl[0]*ul_0[1]+0.1767766952966368*fc[0]*uc_0[1]+0.25*fl[0]*uMax_l[1]-0.25*fc[0]*uMax_l[1]+0.1767766952966368*ul_0[0]*fl[1]+0.25*uMax_l[0]*fl[1]+0.1767766952966368*uc_0[0]*fc[1]-0.25*uMax_l[0]*fc[1]; 

  Ghat_r[0] = 0.5303300858899105*fr[3]*ur_0[3]-0.3061862178478971*fr[1]*ur_0[3]+0.5303300858899105*fc[3]*uc_0[3]+0.3061862178478971*fc[1]*uc_0[3]-0.3061862178478971*ur_0[1]*fr[3]+0.4330127018922193*uMax_r[1]*fr[3]+0.3061862178478971*uc_0[1]*fc[3]+0.4330127018922193*uMax_r[1]*fc[3]+0.5303300858899105*fr[2]*ur_0[2]-0.3061862178478971*fr[0]*ur_0[2]+0.5303300858899105*fc[2]*uc_0[2]+0.3061862178478971*fc[0]*uc_0[2]-0.3061862178478971*ur_0[0]*fr[2]+0.4330127018922193*uMax_r[0]*fr[2]+0.3061862178478971*uc_0[0]*fc[2]+0.4330127018922193*uMax_r[0]*fc[2]+0.1767766952966368*fr[1]*ur_0[1]+0.1767766952966368*fc[1]*uc_0[1]-0.25*fr[1]*uMax_r[1]+0.25*fc[1]*uMax_r[1]+0.1767766952966368*fr[0]*ur_0[0]+0.1767766952966368*fc[0]*uc_0[0]-0.25*fr[0]*uMax_r[0]+0.25*fc[0]*uMax_r[0]; 
  Ghat_r[1] = 0.5303300858899105*fr[2]*ur_0[3]-0.3061862178478971*fr[0]*ur_0[3]+0.5303300858899105*fc[2]*uc_0[3]+0.3061862178478971*fc[0]*uc_0[3]+0.5303300858899105*ur_0[2]*fr[3]-0.3061862178478971*ur_0[0]*fr[3]+0.4330127018922193*uMax_r[0]*fr[3]+0.5303300858899105*uc_0[2]*fc[3]+0.3061862178478971*uc_0[0]*fc[3]+0.4330127018922193*uMax_r[0]*fc[3]-0.3061862178478971*fr[1]*ur_0[2]+0.3061862178478971*fc[1]*uc_0[2]-0.3061862178478971*ur_0[1]*fr[2]+0.4330127018922193*uMax_r[1]*fr[2]+0.3061862178478971*uc_0[1]*fc[2]+0.4330127018922193*uMax_r[1]*fc[2]+0.1767766952966368*fr[0]*ur_0[1]+0.1767766952966368*fc[0]*uc_0[1]-0.25*fr[0]*uMax_r[1]+0.25*fc[0]*uMax_r[1]+0.1767766952966368*ur_0[0]*fr[1]-0.25*uMax_r[0]*fr[1]+0.1767766952966368*uc_0[0]*fc[1]+0.25*uMax_r[0]*fc[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 

  return 0.;

} 
