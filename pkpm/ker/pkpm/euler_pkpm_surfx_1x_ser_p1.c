#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r,
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r,
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates.
  // dxv[NDIM]:              Cell spacing.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // prim_surf_l/c/r:        Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // p_ij_l/c/r:             Input volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in left/center/right cells.
  // euler_pkpm_l/c/r:       Input [rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // pkpm_lax_l/r:           Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m) on left/right surface.
  // pkpm_penalization_l/r:  Surface expansion of total momentum penalization on left/right surface.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[2]; 
  const double *rhouz_l = &euler_pkpm_l[4]; 

  const double *rhoux_c = &euler_pkpm_c[0]; 
  const double *rhouy_c = &euler_pkpm_c[2]; 
  const double *rhouz_c = &euler_pkpm_c[4]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[2]; 
  const double *rhouz_r = &euler_pkpm_r[4]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[2]; 
  const double *Pxz_l = &p_ij_l[4]; 

  const double *Pxx_c = &p_ij_c[0]; 
  const double *Pxy_c = &p_ij_c[2]; 
  const double *Pxz_c = &p_ij_c[4]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[2]; 
  const double *Pxz_r = &p_ij_r[4]; 

  const double *ux_surf_lr = &prim_surf_l[1]; 
  const double *uy_surf_lr = &prim_surf_l[3]; 
  const double *uz_surf_lr = &prim_surf_l[5]; 

  const double *ux_surf_cl = &prim_surf_c[0]; 
  const double *uy_surf_cl = &prim_surf_c[2]; 
  const double *uz_surf_cl = &prim_surf_c[4]; 

  const double *ux_surf_cr = &prim_surf_c[1]; 
  const double *uy_surf_cr = &prim_surf_c[3]; 
  const double *uz_surf_cr = &prim_surf_c[5]; 

  const double *ux_surf_rl = &prim_surf_r[0]; 
  const double *uy_surf_rl = &prim_surf_r[2]; 
  const double *uz_surf_rl = &prim_surf_r[4]; 

  const double *pkpm_penalization_rhoux_l = &pkpm_penalization_l[0]; 
  const double *pkpm_penalization_rhouy_l = &pkpm_penalization_l[1]; 
  const double *pkpm_penalization_rhouz_l = &pkpm_penalization_l[2]; 

  const double *pkpm_penalization_rhoux_r = &pkpm_penalization_r[0]; 
  const double *pkpm_penalization_rhouy_r = &pkpm_penalization_r[1]; 
  const double *pkpm_penalization_rhouz_r = &pkpm_penalization_r[2]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double uxl_r = ux_surf_lr[0]; 
  double uxc_l = ux_surf_cl[0]; 
  double uxc_r = ux_surf_cr[0]; 
  double uxr_l = ux_surf_rl[0]; 

  double uxl_r_sq = uxl_r*uxl_r; 
  double uxc_l_sq = uxc_l*uxc_l; 
  double uxc_r_sq = uxc_r*uxc_r; 
  double uxr_l_sq = uxr_l*uxr_l; 

  double uyl_r = uy_surf_lr[0]; 
  double uyc_l = uy_surf_cl[0]; 
  double uyc_r = uy_surf_cr[0]; 
  double uyr_l = uy_surf_rl[0]; 

  double uzl_r = uz_surf_lr[0]; 
  double uzc_l = uz_surf_cl[0]; 
  double uzc_r = uz_surf_cr[0]; 
  double uzr_l = uz_surf_rl[0]; 

  double max_speed_l = pkpm_lax_l[0]; 
  double max_speed_r = pkpm_lax_r[0]; 

  Ghat_rhoux_l = 0.1530931089239486*rho_l[1]*uxl_r_sq-0.1530931089239486*rho_c[1]*uxl_r_sq+0.0883883476483184*rho_l[0]*uxl_r_sq+0.0883883476483184*rho_c[0]*uxl_r_sq+0.3061862178478971*rho_l[1]*uxc_l*uxl_r-0.3061862178478971*rho_c[1]*uxc_l*uxl_r+0.1767766952966368*rho_l[0]*uxc_l*uxl_r+0.1767766952966368*rho_c[0]*uxc_l*uxl_r+0.3061862178478971*rho_l[1]*max_speed_l*uxl_r+0.3061862178478971*rho_c[1]*max_speed_l*uxl_r+0.1767766952966368*rho_l[0]*max_speed_l*uxl_r-0.1767766952966368*rho_c[0]*max_speed_l*uxl_r+0.1530931089239486*rho_l[1]*uxc_l_sq-0.1530931089239486*rho_c[1]*uxc_l_sq+0.0883883476483184*rho_l[0]*uxc_l_sq+0.0883883476483184*rho_c[0]*uxc_l_sq+0.3061862178478971*rho_l[1]*max_speed_l*uxc_l+0.3061862178478971*rho_c[1]*max_speed_l*uxc_l+0.1767766952966368*rho_l[0]*max_speed_l*uxc_l-0.1767766952966368*rho_c[0]*max_speed_l*uxc_l+0.6123724356957944*Pxx_l[1]-0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0] - pkpm_penalization_rhoux_l[0]; 
  Ghat_rhouy_l = 0.1530931089239486*rho_l[1]*uxl_r*uyl_r-0.1530931089239486*rho_c[1]*uxl_r*uyl_r+0.08838834764831843*rho_l[0]*uxl_r*uyl_r+0.08838834764831843*rho_c[0]*uxl_r*uyl_r+0.1530931089239486*rho_l[1]*uxc_l*uyl_r-0.1530931089239486*rho_c[1]*uxc_l*uyl_r+0.08838834764831843*rho_l[0]*uxc_l*uyl_r+0.08838834764831843*rho_c[0]*uxc_l*uyl_r+0.3061862178478972*rho_l[1]*max_speed_l*uyl_r+0.3061862178478972*rho_c[1]*max_speed_l*uyl_r+0.1767766952966369*rho_l[0]*max_speed_l*uyl_r-0.1767766952966369*rho_c[0]*max_speed_l*uyl_r+0.1530931089239486*rho_l[1]*uxl_r*uyc_l-0.1530931089239486*rho_c[1]*uxl_r*uyc_l+0.08838834764831843*rho_l[0]*uxl_r*uyc_l+0.08838834764831843*rho_c[0]*uxl_r*uyc_l+0.1530931089239486*rho_l[1]*uxc_l*uyc_l-0.1530931089239486*rho_c[1]*uxc_l*uyc_l+0.08838834764831843*rho_l[0]*uxc_l*uyc_l+0.08838834764831843*rho_c[0]*uxc_l*uyc_l+0.3061862178478972*rho_l[1]*max_speed_l*uyc_l+0.3061862178478972*rho_c[1]*max_speed_l*uyc_l+0.1767766952966369*rho_l[0]*max_speed_l*uyc_l-0.1767766952966369*rho_c[0]*max_speed_l*uyc_l+0.6123724356957945*Pxy_l[1]-0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0] - pkpm_penalization_rhouy_l[0]; 
  Ghat_rhouz_l = 0.1530931089239486*rho_l[1]*uxl_r*uzl_r-0.1530931089239486*rho_c[1]*uxl_r*uzl_r+0.08838834764831843*rho_l[0]*uxl_r*uzl_r+0.08838834764831843*rho_c[0]*uxl_r*uzl_r+0.1530931089239486*rho_l[1]*uxc_l*uzl_r-0.1530931089239486*rho_c[1]*uxc_l*uzl_r+0.08838834764831843*rho_l[0]*uxc_l*uzl_r+0.08838834764831843*rho_c[0]*uxc_l*uzl_r+0.3061862178478972*rho_l[1]*max_speed_l*uzl_r+0.3061862178478972*rho_c[1]*max_speed_l*uzl_r+0.1767766952966369*rho_l[0]*max_speed_l*uzl_r-0.1767766952966369*rho_c[0]*max_speed_l*uzl_r+0.1530931089239486*rho_l[1]*uxl_r*uzc_l-0.1530931089239486*rho_c[1]*uxl_r*uzc_l+0.08838834764831843*rho_l[0]*uxl_r*uzc_l+0.08838834764831843*rho_c[0]*uxl_r*uzc_l+0.1530931089239486*rho_l[1]*uxc_l*uzc_l-0.1530931089239486*rho_c[1]*uxc_l*uzc_l+0.08838834764831843*rho_l[0]*uxc_l*uzc_l+0.08838834764831843*rho_c[0]*uxc_l*uzc_l+0.3061862178478972*rho_l[1]*max_speed_l*uzc_l+0.3061862178478972*rho_c[1]*max_speed_l*uzc_l+0.1767766952966369*rho_l[0]*max_speed_l*uzc_l-0.1767766952966369*rho_c[0]*max_speed_l*uzc_l+0.6123724356957945*Pxz_l[1]-0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0] - pkpm_penalization_rhouz_l[0]; 
  Ghat_rhoux_r = (-0.1530931089239486*rho_r[1]*uxr_l_sq)+0.1530931089239486*rho_c[1]*uxr_l_sq+0.0883883476483184*rho_r[0]*uxr_l_sq+0.0883883476483184*rho_c[0]*uxr_l_sq-0.3061862178478971*rho_r[1]*uxc_r*uxr_l+0.3061862178478971*rho_c[1]*uxc_r*uxr_l+0.1767766952966368*rho_r[0]*uxc_r*uxr_l+0.1767766952966368*rho_c[0]*uxc_r*uxr_l+0.3061862178478971*rho_r[1]*max_speed_r*uxr_l+0.3061862178478971*rho_c[1]*max_speed_r*uxr_l-0.1767766952966368*rho_r[0]*max_speed_r*uxr_l+0.1767766952966368*rho_c[0]*max_speed_r*uxr_l-0.1530931089239486*rho_r[1]*uxc_r_sq+0.1530931089239486*rho_c[1]*uxc_r_sq+0.0883883476483184*rho_r[0]*uxc_r_sq+0.0883883476483184*rho_c[0]*uxc_r_sq+0.3061862178478971*rho_r[1]*max_speed_r*uxc_r+0.3061862178478971*rho_c[1]*max_speed_r*uxc_r-0.1767766952966368*rho_r[0]*max_speed_r*uxc_r+0.1767766952966368*rho_c[0]*max_speed_r*uxc_r-0.6123724356957944*Pxx_r[1]+0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0] - pkpm_penalization_rhoux_r[0]; 
  Ghat_rhouy_r = (-0.1530931089239486*rho_r[1]*uxr_l*uyr_l)+0.1530931089239486*rho_c[1]*uxr_l*uyr_l+0.08838834764831843*rho_r[0]*uxr_l*uyr_l+0.08838834764831843*rho_c[0]*uxr_l*uyr_l-0.1530931089239486*rho_r[1]*uxc_r*uyr_l+0.1530931089239486*rho_c[1]*uxc_r*uyr_l+0.08838834764831843*rho_r[0]*uxc_r*uyr_l+0.08838834764831843*rho_c[0]*uxc_r*uyr_l+0.3061862178478972*rho_r[1]*max_speed_r*uyr_l+0.3061862178478972*rho_c[1]*max_speed_r*uyr_l-0.1767766952966369*rho_r[0]*max_speed_r*uyr_l+0.1767766952966369*rho_c[0]*max_speed_r*uyr_l-0.1530931089239486*rho_r[1]*uxr_l*uyc_r+0.1530931089239486*rho_c[1]*uxr_l*uyc_r+0.08838834764831843*rho_r[0]*uxr_l*uyc_r+0.08838834764831843*rho_c[0]*uxr_l*uyc_r-0.1530931089239486*rho_r[1]*uxc_r*uyc_r+0.1530931089239486*rho_c[1]*uxc_r*uyc_r+0.08838834764831843*rho_r[0]*uxc_r*uyc_r+0.08838834764831843*rho_c[0]*uxc_r*uyc_r+0.3061862178478972*rho_r[1]*max_speed_r*uyc_r+0.3061862178478972*rho_c[1]*max_speed_r*uyc_r-0.1767766952966369*rho_r[0]*max_speed_r*uyc_r+0.1767766952966369*rho_c[0]*max_speed_r*uyc_r-0.6123724356957945*Pxy_r[1]+0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0] - pkpm_penalization_rhouy_r[0]; 
  Ghat_rhouz_r = (-0.1530931089239486*rho_r[1]*uxr_l*uzr_l)+0.1530931089239486*rho_c[1]*uxr_l*uzr_l+0.08838834764831843*rho_r[0]*uxr_l*uzr_l+0.08838834764831843*rho_c[0]*uxr_l*uzr_l-0.1530931089239486*rho_r[1]*uxc_r*uzr_l+0.1530931089239486*rho_c[1]*uxc_r*uzr_l+0.08838834764831843*rho_r[0]*uxc_r*uzr_l+0.08838834764831843*rho_c[0]*uxc_r*uzr_l+0.3061862178478972*rho_r[1]*max_speed_r*uzr_l+0.3061862178478972*rho_c[1]*max_speed_r*uzr_l-0.1767766952966369*rho_r[0]*max_speed_r*uzr_l+0.1767766952966369*rho_c[0]*max_speed_r*uzr_l-0.1530931089239486*rho_r[1]*uxr_l*uzc_r+0.1530931089239486*rho_c[1]*uxr_l*uzc_r+0.08838834764831843*rho_r[0]*uxr_l*uzc_r+0.08838834764831843*rho_c[0]*uxr_l*uzc_r-0.1530931089239486*rho_r[1]*uxc_r*uzc_r+0.1530931089239486*rho_c[1]*uxc_r*uzc_r+0.08838834764831843*rho_r[0]*uxc_r*uzc_r+0.08838834764831843*rho_c[0]*uxc_r*uzc_r+0.3061862178478972*rho_r[1]*max_speed_r*uzc_r+0.3061862178478972*rho_c[1]*max_speed_r*uzc_r-0.1767766952966369*rho_r[0]*max_speed_r*uzc_r+0.1767766952966369*rho_c[0]*max_speed_r*uzc_r-0.6123724356957945*Pxz_r[1]+0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0] - pkpm_penalization_rhouz_r[0]; 

  outrhou0[0] += (0.7071067811865475*Ghat_rhoux_l-0.7071067811865475*Ghat_rhoux_r)*dx1; 
  outrhou0[1] += -1.224744871391589*(Ghat_rhoux_r+Ghat_rhoux_l)*dx1; 

  outrhou1[0] += (0.7071067811865475*Ghat_rhouy_l-0.7071067811865475*Ghat_rhouy_r)*dx1; 
  outrhou1[1] += -1.224744871391589*(Ghat_rhouy_r+Ghat_rhouy_l)*dx1; 

  outrhou2[0] += (0.7071067811865475*Ghat_rhouz_l-0.7071067811865475*Ghat_rhouz_r)*dx1; 
  outrhou2[1] += -1.224744871391589*(Ghat_rhouz_r+Ghat_rhouz_l)*dx1; 

  return 0.;

} 
