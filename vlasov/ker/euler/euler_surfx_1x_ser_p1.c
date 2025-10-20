#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_surfx_1x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *u_surf_l, const double *u_surf_c, const double *u_surf_r,
    const double *p_surf_l, const double *p_surf_c, const double *p_surf_r,
    const double *fluid_l, const double *fluid_c, const double *fluid_r, 
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:      Cell-center coordinates.
  // dxv[NDIM]:    Cell spacing.
  // wv_eqn:       Wave equation for computing fluctuations at the interface for upwinding.
  // geom_l:       Geometry for the left surface update.
  // geom_r:       Geometry for the right surface update.
  // u_surf_l/c/r: Input surface expansion of flow velocity in left/center/right cells in each direction.
  // p_surf_l/c/r: Input surface expansion of pressure in left/center/right cells in each direction.
  //               [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 
  // fluid_l/c/r:  [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out:          Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rho_l = &fluid_l[0]; 
  const double *rhoux_l = &fluid_l[2]; 
  const double *rhouy_l = &fluid_l[4]; 
  const double *rhouz_l = &fluid_l[6]; 
  const double *energy_l = &fluid_l[8]; 

  const double *rho_c = &fluid_c[0]; 
  const double *rhoux_c = &fluid_c[2]; 
  const double *rhouy_c = &fluid_c[4]; 
  const double *rhouz_c = &fluid_c[6]; 
  const double *energy_c = &fluid_c[8]; 

  const double *rho_r = &fluid_r[0]; 
  const double *rhoux_r = &fluid_r[2]; 
  const double *rhouy_r = &fluid_r[4]; 
  const double *rhouz_r = &fluid_r[6]; 
  const double *energy_r = &fluid_r[8]; 

  const double *ux_surf_lr = &u_surf_l[1]; 
  const double *uy_surf_lr = &u_surf_l[3]; 
  const double *uz_surf_lr = &u_surf_l[5]; 

  const double *ux_surf_cl = &u_surf_c[0]; 
  const double *uy_surf_cl = &u_surf_c[2]; 
  const double *uz_surf_cl = &u_surf_c[4]; 

  const double *ux_surf_cr = &u_surf_c[1]; 
  const double *uy_surf_cr = &u_surf_c[3]; 
  const double *uz_surf_cr = &u_surf_c[5]; 

  const double *ux_surf_rl = &u_surf_r[0]; 
  const double *uy_surf_rl = &u_surf_r[2]; 
  const double *uz_surf_rl = &u_surf_r[4]; 

  const double *p_surf_lr = &p_surf_l[1]; 
  const double *p_surf_cl = &p_surf_c[0]; 
  const double *p_surf_cr = &p_surf_c[1]; 
  const double *p_surf_rl = &p_surf_r[0]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[2]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[6]; 
  double *outenergy = &out[8]; 

  double Ghat_rho_l = 0.0; 
  double Ghat_rho_r = 0.0; 
  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double Ghat_energy_l = 0.0; 
  double Ghat_energy_r = 0.0; 
  double uxl_r = ux_surf_lr[0]; 
  double uxc_l = ux_surf_cl[0]; 
  double uxc_r = ux_surf_cr[0]; 
  double uxr_l = ux_surf_rl[0]; 

  double pl_r = p_surf_lr[0]; 
  double pc_l = p_surf_cl[0]; 
  double pc_r = p_surf_cr[0]; 
  double pr_l = p_surf_rl[0]; 

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

  double q_lr[5] = {0.0}; 
  double q_cl[5] = {0.0}; 
  double q_cr[5] = {0.0}; 
  double q_rl[5] = {0.0}; 
  q_lr[0] = 1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0]; 
  q_lr[1] = 1.224744871391589*rhoux_l[1]+0.7071067811865475*rhoux_l[0]; 
  q_lr[2] = 1.224744871391589*rhouy_l[1]+0.7071067811865475*rhouy_l[0]; 
  q_lr[3] = 1.224744871391589*rhouz_l[1]+0.7071067811865475*rhouz_l[0]; 
  q_lr[4] = 1.224744871391589*energy_l[1]+0.7071067811865475*energy_l[0]; 
  q_cl[0] = 0.7071067811865475*rho_c[0]-1.224744871391589*rho_c[1]; 
  q_cl[1] = 0.7071067811865475*rhoux_c[0]-1.224744871391589*rhoux_c[1]; 
  q_cl[2] = 0.7071067811865475*rhouy_c[0]-1.224744871391589*rhouy_c[1]; 
  q_cl[3] = 0.7071067811865475*rhouz_c[0]-1.224744871391589*rhouz_c[1]; 
  q_cl[4] = 0.7071067811865475*energy_c[0]-1.224744871391589*energy_c[1]; 
  q_cr[0] = 1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0]; 
  q_cr[1] = 1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0]; 
  q_cr[2] = 1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0]; 
  q_cr[3] = 1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0]; 
  q_cr[4] = 1.224744871391589*energy_c[1]+0.7071067811865475*energy_c[0]; 
  q_rl[0] = 0.7071067811865475*rho_r[0]-1.224744871391589*rho_r[1]; 
  q_rl[1] = 0.7071067811865475*rhoux_r[0]-1.224744871391589*rhoux_r[1]; 
  q_rl[2] = 0.7071067811865475*rhouy_r[0]-1.224744871391589*rhouy_r[1]; 
  q_rl[3] = 0.7071067811865475*rhouz_r[0]-1.224744871391589*rhouz_r[1]; 
  q_rl[4] = 0.7071067811865475*energy_r[0]-1.224744871391589*energy_r[1]; 

  double q_lr_local[5] = {0.0}; 
  double q_cl_local[5] = {0.0}; 
  double q_cr_local[5] = {0.0}; 
  double q_rl_local[5] = {0.0}; 
  rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  double delta_l[5] = {0.0}; 
  double delta_r[5] = {0.0}; 
  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 

  double waves_l[15] = {0.0}; 
  double waves_r[15] = {0.0}; 
  double speeds_l[3] = {0.0}; 
  double speeds_r[3] = {0.0}; 
  double my_max_speed_l = wave_roe(wv_eqn, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  double my_max_speed_r = wave_roe(wv_eqn, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 

  double lenr_l = geom_l->lenr[0]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  double lenr_r = geom_r->lenr[0]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 

  double amdq_l_local[5] = {0.0}; 
  double apdq_l_local[5] = {0.0}; 
  double amdq_r_local[5] = {0.0}; 
  double apdq_r_local[5] = {0.0}; 
  qfluct_roe(wv_eqn, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  qfluct_roe(wv_eqn, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  double amdq_l[5] = {0.0}; 
  double apdq_l[5] = {0.0}; 
  double amdq_r[5] = {0.0}; 
  double apdq_r[5] = {0.0}; 
  rot_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], amdq_l_local, amdq_l); 
  rot_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], apdq_l_local, apdq_l); 
  rot_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], amdq_r_local, amdq_r); 
  rot_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], apdq_r_local, apdq_r); 

  Ghat_rho_l = 0.3061862178478972*rho_l[1]*uxl_r-0.3061862178478972*rho_c[1]*uxl_r+0.1767766952966369*rho_l[0]*uxl_r+0.1767766952966369*rho_c[0]*uxl_r+0.3061862178478972*rho_l[1]*uxc_l-0.3061862178478972*rho_c[1]*uxc_l+0.1767766952966369*rho_l[0]*uxc_l+0.1767766952966369*rho_c[0]*uxc_l-0.5*apdq_l[0]+0.5*amdq_l[0]; 
  Ghat_rhoux_l = 0.5*Ghat_rho_l*uxl_r+0.5*Ghat_rho_l*uxc_l+0.5*pl_r+0.5*pc_l-0.5*apdq_l[1]+0.5*amdq_l[1]; 
  Ghat_rhouy_l = 0.5*Ghat_rho_l*uyl_r+0.5*Ghat_rho_l*uyc_l-0.5*apdq_l[2]+0.5*amdq_l[2]; 
  Ghat_rhouz_l = 0.5*Ghat_rho_l*uzl_r+0.5*Ghat_rho_l*uzc_l-0.5*apdq_l[3]+0.5*amdq_l[3]; 
  Ghat_energy_l = 0.25*pl_r*uxl_r+0.25*pc_l*uxl_r+0.3061862178478972*energy_l[1]*uxl_r-0.3061862178478972*energy_c[1]*uxl_r+0.1767766952966369*energy_l[0]*uxl_r+0.1767766952966369*energy_c[0]*uxl_r+0.25*pl_r*uxc_l+0.25*pc_l*uxc_l+0.3061862178478972*energy_l[1]*uxc_l-0.3061862178478972*energy_c[1]*uxc_l+0.1767766952966369*energy_l[0]*uxc_l+0.1767766952966369*energy_c[0]*uxc_l-0.5*apdq_l[4]+0.5*amdq_l[4]; 
  Ghat_rho_r = (-0.3061862178478972*rho_r[1]*uxr_l)+0.3061862178478972*rho_c[1]*uxr_l+0.1767766952966369*rho_r[0]*uxr_l+0.1767766952966369*rho_c[0]*uxr_l-0.3061862178478972*rho_r[1]*uxc_r+0.3061862178478972*rho_c[1]*uxc_r+0.1767766952966369*rho_r[0]*uxc_r+0.1767766952966369*rho_c[0]*uxc_r-0.5*apdq_r[0]+0.5*amdq_r[0]; 
  Ghat_rhoux_r = 0.5*Ghat_rho_r*uxr_l+0.5*Ghat_rho_r*uxc_r+0.5*pr_l+0.5*pc_r-0.5*apdq_r[1]+0.5*amdq_r[1]; 
  Ghat_rhouy_r = 0.5*Ghat_rho_r*uyr_l+0.5*Ghat_rho_r*uyc_r-0.5*apdq_r[2]+0.5*amdq_r[2]; 
  Ghat_rhouz_r = 0.5*Ghat_rho_r*uzr_l+0.5*Ghat_rho_r*uzc_r-0.5*apdq_r[3]+0.5*amdq_r[3]; 
  Ghat_energy_r = 0.25*pr_l*uxr_l+0.25*pc_r*uxr_l-0.3061862178478972*energy_r[1]*uxr_l+0.3061862178478972*energy_c[1]*uxr_l+0.1767766952966369*energy_r[0]*uxr_l+0.1767766952966369*energy_c[0]*uxr_l+0.25*pr_l*uxc_r+0.25*pc_r*uxc_r-0.3061862178478972*energy_r[1]*uxc_r+0.3061862178478972*energy_c[1]*uxc_r+0.1767766952966369*energy_r[0]*uxc_r+0.1767766952966369*energy_c[0]*uxc_r-0.5*apdq_r[4]+0.5*amdq_r[4]; 

  outrho[0] += (0.7071067811865475*Ghat_rho_l-0.7071067811865475*Ghat_rho_r)*dx1; 
  outrho[1] += -1.224744871391589*(Ghat_rho_r+Ghat_rho_l)*dx1; 

  outrhoux[0] += (0.7071067811865475*Ghat_rhoux_l-0.7071067811865475*Ghat_rhoux_r)*dx1; 
  outrhoux[1] += -1.224744871391589*(Ghat_rhoux_r+Ghat_rhoux_l)*dx1; 

  outrhouy[0] += (0.7071067811865475*Ghat_rhouy_l-0.7071067811865475*Ghat_rhouy_r)*dx1; 
  outrhouy[1] += -1.224744871391589*(Ghat_rhouy_r+Ghat_rhouy_l)*dx1; 

  outrhouz[0] += (0.7071067811865475*Ghat_rhouz_l-0.7071067811865475*Ghat_rhouz_r)*dx1; 
  outrhouz[1] += -1.224744871391589*(Ghat_rhouz_r+Ghat_rhouz_l)*dx1; 

  outenergy[0] += (0.7071067811865475*Ghat_energy_l-0.7071067811865475*Ghat_energy_r)*dx1; 
  outenergy[1] += -1.224744871391589*(Ghat_energy_r+Ghat_energy_l)*dx1; 

  return 0.;

} 
