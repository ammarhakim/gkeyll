#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv,
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

  const double dx1 = 2.0/dxv[1]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[4]; 
  const double *rhouz_l = &euler_pkpm_l[8]; 

  const double *rhoux_c = &euler_pkpm_c[0]; 
  const double *rhouy_c = &euler_pkpm_c[4]; 
  const double *rhouz_c = &euler_pkpm_c[8]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[4]; 
  const double *rhouz_r = &euler_pkpm_r[8]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *ux_surf_lr = &prim_surf_l[18]; 
  const double *uy_surf_lr = &prim_surf_l[22]; 
  const double *uz_surf_lr = &prim_surf_l[26]; 

  const double *ux_surf_cl = &prim_surf_c[16]; 
  const double *uy_surf_cl = &prim_surf_c[20]; 
  const double *uz_surf_cl = &prim_surf_c[24]; 

  const double *ux_surf_cr = &prim_surf_c[18]; 
  const double *uy_surf_cr = &prim_surf_c[22]; 
  const double *uz_surf_cr = &prim_surf_c[26]; 

  const double *ux_surf_rl = &prim_surf_r[16]; 
  const double *uy_surf_rl = &prim_surf_r[20]; 
  const double *uz_surf_rl = &prim_surf_r[24]; 

  const double *pkpm_lax_dir_l = &pkpm_lax_l[2]; 
  const double *pkpm_lax_dir_r = &pkpm_lax_r[2]; 

  const double *pkpm_penalization_rhoux_l = &pkpm_penalization_l[6]; 
  const double *pkpm_penalization_rhouy_l = &pkpm_penalization_l[8]; 
  const double *pkpm_penalization_rhouz_l = &pkpm_penalization_l[10]; 

  const double *pkpm_penalization_rhoux_r = &pkpm_penalization_r[6]; 
  const double *pkpm_penalization_rhouy_r = &pkpm_penalization_r[8]; 
  const double *pkpm_penalization_rhouz_r = &pkpm_penalization_r[10]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 

  double flux_rho_l[2] = {0.0}; 
  double flux_rho_r[2] = {0.0}; 

  double avg_p_ij_x_l[2] = {0.0}; 
  double avg_p_ij_x_r[2] = {0.0}; 
  double avg_p_ij_y_l[2] = {0.0}; 
  double avg_p_ij_y_r[2] = {0.0}; 
  double avg_p_ij_z_l[2] = {0.0}; 
  double avg_p_ij_z_r[2] = {0.0}; 

  const double *Pxy_l = &p_ij_l[4]; 
  const double *Pyy_l = &p_ij_l[12]; 
  const double *Pyz_l = &p_ij_l[16]; 

  const double *Pxy_c = &p_ij_c[4]; 
  const double *Pyy_c = &p_ij_c[12]; 
  const double *Pyz_c = &p_ij_c[16]; 

  const double *Pxy_r = &p_ij_r[4]; 
  const double *Pyy_r = &p_ij_r[12]; 
  const double *Pyz_r = &p_ij_r[16]; 

  flux_rho_l[0] = 0.2165063509461096*uy_surf_lr[1]*rho_l[3]+0.2165063509461096*uy_surf_cl[1]*rho_l[3]+0.4330127018922193*pkpm_lax_dir_l[1]*rho_l[3]-0.2165063509461096*uy_surf_lr[1]*rho_c[3]-0.2165063509461096*uy_surf_cl[1]*rho_c[3]+0.4330127018922193*pkpm_lax_dir_l[1]*rho_c[3]+0.2165063509461096*uy_surf_lr[0]*rho_l[2]+0.2165063509461096*uy_surf_cl[0]*rho_l[2]+0.4330127018922193*pkpm_lax_dir_l[0]*rho_l[2]-0.2165063509461096*uy_surf_lr[0]*rho_c[2]-0.2165063509461096*uy_surf_cl[0]*rho_c[2]+0.4330127018922193*pkpm_lax_dir_l[0]*rho_c[2]+0.125*rho_l[1]*uy_surf_lr[1]+0.125*rho_c[1]*uy_surf_lr[1]+0.125*rho_l[1]*uy_surf_cl[1]+0.125*rho_c[1]*uy_surf_cl[1]+0.25*pkpm_lax_dir_l[1]*rho_l[1]-0.25*pkpm_lax_dir_l[1]*rho_c[1]+0.125*rho_l[0]*uy_surf_lr[0]+0.125*rho_c[0]*uy_surf_lr[0]+0.125*rho_l[0]*uy_surf_cl[0]+0.125*rho_c[0]*uy_surf_cl[0]+0.25*pkpm_lax_dir_l[0]*rho_l[0]-0.25*pkpm_lax_dir_l[0]*rho_c[0]; 
  flux_rho_l[1] = 0.2165063509461096*uy_surf_lr[0]*rho_l[3]+0.2165063509461096*uy_surf_cl[0]*rho_l[3]+0.4330127018922193*pkpm_lax_dir_l[0]*rho_l[3]-0.2165063509461096*uy_surf_lr[0]*rho_c[3]-0.2165063509461096*uy_surf_cl[0]*rho_c[3]+0.4330127018922193*pkpm_lax_dir_l[0]*rho_c[3]+0.2165063509461096*uy_surf_lr[1]*rho_l[2]+0.2165063509461096*uy_surf_cl[1]*rho_l[2]+0.4330127018922193*pkpm_lax_dir_l[1]*rho_l[2]-0.2165063509461096*uy_surf_lr[1]*rho_c[2]-0.2165063509461096*uy_surf_cl[1]*rho_c[2]+0.4330127018922193*pkpm_lax_dir_l[1]*rho_c[2]+0.125*rho_l[0]*uy_surf_lr[1]+0.125*rho_c[0]*uy_surf_lr[1]+0.125*rho_l[0]*uy_surf_cl[1]+0.125*rho_c[0]*uy_surf_cl[1]+0.125*uy_surf_lr[0]*rho_l[1]+0.125*uy_surf_cl[0]*rho_l[1]+0.25*pkpm_lax_dir_l[0]*rho_l[1]+0.125*uy_surf_lr[0]*rho_c[1]+0.125*uy_surf_cl[0]*rho_c[1]-0.25*pkpm_lax_dir_l[0]*rho_c[1]+0.25*rho_l[0]*pkpm_lax_dir_l[1]-0.25*rho_c[0]*pkpm_lax_dir_l[1]; 

  flux_rho_r[0] = (-0.2165063509461096*uy_surf_rl[1]*rho_r[3])-0.2165063509461096*uy_surf_cr[1]*rho_r[3]+0.4330127018922193*pkpm_lax_dir_r[1]*rho_r[3]+0.2165063509461096*uy_surf_rl[1]*rho_c[3]+0.2165063509461096*uy_surf_cr[1]*rho_c[3]+0.4330127018922193*pkpm_lax_dir_r[1]*rho_c[3]-0.2165063509461096*uy_surf_rl[0]*rho_r[2]-0.2165063509461096*uy_surf_cr[0]*rho_r[2]+0.4330127018922193*pkpm_lax_dir_r[0]*rho_r[2]+0.2165063509461096*uy_surf_rl[0]*rho_c[2]+0.2165063509461096*uy_surf_cr[0]*rho_c[2]+0.4330127018922193*pkpm_lax_dir_r[0]*rho_c[2]+0.125*rho_r[1]*uy_surf_rl[1]+0.125*rho_c[1]*uy_surf_rl[1]+0.125*rho_r[1]*uy_surf_cr[1]+0.125*rho_c[1]*uy_surf_cr[1]-0.25*pkpm_lax_dir_r[1]*rho_r[1]+0.25*pkpm_lax_dir_r[1]*rho_c[1]+0.125*rho_r[0]*uy_surf_rl[0]+0.125*rho_c[0]*uy_surf_rl[0]+0.125*rho_r[0]*uy_surf_cr[0]+0.125*rho_c[0]*uy_surf_cr[0]-0.25*pkpm_lax_dir_r[0]*rho_r[0]+0.25*pkpm_lax_dir_r[0]*rho_c[0]; 
  flux_rho_r[1] = (-0.2165063509461096*uy_surf_rl[0]*rho_r[3])-0.2165063509461096*uy_surf_cr[0]*rho_r[3]+0.4330127018922193*pkpm_lax_dir_r[0]*rho_r[3]+0.2165063509461096*uy_surf_rl[0]*rho_c[3]+0.2165063509461096*uy_surf_cr[0]*rho_c[3]+0.4330127018922193*pkpm_lax_dir_r[0]*rho_c[3]-0.2165063509461096*uy_surf_rl[1]*rho_r[2]-0.2165063509461096*uy_surf_cr[1]*rho_r[2]+0.4330127018922193*pkpm_lax_dir_r[1]*rho_r[2]+0.2165063509461096*uy_surf_rl[1]*rho_c[2]+0.2165063509461096*uy_surf_cr[1]*rho_c[2]+0.4330127018922193*pkpm_lax_dir_r[1]*rho_c[2]+0.125*rho_r[0]*uy_surf_rl[1]+0.125*rho_c[0]*uy_surf_rl[1]+0.125*rho_r[0]*uy_surf_cr[1]+0.125*rho_c[0]*uy_surf_cr[1]+0.125*uy_surf_rl[0]*rho_r[1]+0.125*uy_surf_cr[0]*rho_r[1]-0.25*pkpm_lax_dir_r[0]*rho_r[1]+0.125*uy_surf_rl[0]*rho_c[1]+0.125*uy_surf_cr[0]*rho_c[1]+0.25*pkpm_lax_dir_r[0]*rho_c[1]-0.25*rho_r[0]*pkpm_lax_dir_r[1]+0.25*rho_c[0]*pkpm_lax_dir_r[1]; 

  avg_p_ij_x_l[0] = 0.6123724356957944*Pxy_l[2]-0.6123724356957944*Pxy_c[2]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_x_l[1] = 0.6123724356957944*Pxy_l[3]-0.6123724356957944*Pxy_c[3]+0.3535533905932737*Pxy_l[1]+0.3535533905932737*Pxy_c[1]; 

  avg_p_ij_x_r[0] = (-0.6123724356957944*Pxy_r[2])+0.6123724356957944*Pxy_c[2]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_x_r[1] = (-0.6123724356957944*Pxy_r[3])+0.6123724356957944*Pxy_c[3]+0.3535533905932737*Pxy_r[1]+0.3535533905932737*Pxy_c[1]; 

  avg_p_ij_y_l[0] = 0.6123724356957944*Pyy_l[2]-0.6123724356957944*Pyy_c[2]+0.3535533905932737*Pyy_l[0]+0.3535533905932737*Pyy_c[0]; 
  avg_p_ij_y_l[1] = 0.6123724356957944*Pyy_l[3]-0.6123724356957944*Pyy_c[3]+0.3535533905932737*Pyy_l[1]+0.3535533905932737*Pyy_c[1]; 

  avg_p_ij_y_r[0] = (-0.6123724356957944*Pyy_r[2])+0.6123724356957944*Pyy_c[2]+0.3535533905932737*Pyy_r[0]+0.3535533905932737*Pyy_c[0]; 
  avg_p_ij_y_r[1] = (-0.6123724356957944*Pyy_r[3])+0.6123724356957944*Pyy_c[3]+0.3535533905932737*Pyy_r[1]+0.3535533905932737*Pyy_c[1]; 

  avg_p_ij_z_l[0] = 0.6123724356957944*Pyz_l[2]-0.6123724356957944*Pyz_c[2]+0.3535533905932737*Pyz_l[0]+0.3535533905932737*Pyz_c[0]; 
  avg_p_ij_z_l[1] = 0.6123724356957944*Pyz_l[3]-0.6123724356957944*Pyz_c[3]+0.3535533905932737*Pyz_l[1]+0.3535533905932737*Pyz_c[1]; 

  avg_p_ij_z_r[0] = (-0.6123724356957944*Pyz_r[2])+0.6123724356957944*Pyz_c[2]+0.3535533905932737*Pyz_r[0]+0.3535533905932737*Pyz_c[0]; 
  avg_p_ij_z_r[1] = (-0.6123724356957944*Pyz_r[3])+0.6123724356957944*Pyz_c[3]+0.3535533905932737*Pyz_r[1]+0.3535533905932737*Pyz_c[1]; 

  outrhou0[0] += ((-0.25*flux_rho_r[1]*ux_surf_rl[1])+0.25*flux_rho_l[1]*ux_surf_lr[1]-0.25*flux_rho_r[1]*ux_surf_cr[1]+0.25*flux_rho_l[1]*ux_surf_cl[1]-0.25*flux_rho_r[0]*ux_surf_rl[0]+0.25*flux_rho_l[0]*ux_surf_lr[0]-0.25*flux_rho_r[0]*ux_surf_cr[0]+0.25*flux_rho_l[0]*ux_surf_cl[0]+0.7071067811865475*pkpm_penalization_rhoux_r[0]-0.7071067811865475*(pkpm_penalization_rhoux_l[0]+avg_p_ij_x_r[0])+0.7071067811865475*avg_p_ij_x_l[0])*dx1; 
  outrhou0[1] += ((-0.25*flux_rho_r[0]*ux_surf_rl[1])+0.25*flux_rho_l[0]*ux_surf_lr[1]-0.25*flux_rho_r[0]*ux_surf_cr[1]+0.25*flux_rho_l[0]*ux_surf_cl[1]+0.7071067811865475*pkpm_penalization_rhoux_r[1]-0.7071067811865475*pkpm_penalization_rhoux_l[1]-0.25*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+0.25*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_x_r[1]+0.7071067811865475*avg_p_ij_x_l[1])*dx1; 
  outrhou0[2] += ((-0.4330127018922193*(flux_rho_r[1]*ux_surf_rl[1]+flux_rho_l[1]*ux_surf_lr[1]+flux_rho_r[1]*ux_surf_cr[1]+flux_rho_l[1]*ux_surf_cl[1]+flux_rho_r[0]*ux_surf_rl[0]+flux_rho_l[0]*ux_surf_lr[0]+flux_rho_r[0]*ux_surf_cr[0]+flux_rho_l[0]*ux_surf_cl[0]))+1.224744871391589*(pkpm_penalization_rhoux_r[0]+pkpm_penalization_rhoux_l[0])-1.224744871391589*(avg_p_ij_x_r[0]+avg_p_ij_x_l[0]))*dx1; 
  outrhou0[3] += ((-0.4330127018922193*(flux_rho_r[0]*ux_surf_rl[1]+flux_rho_l[0]*ux_surf_lr[1]+flux_rho_r[0]*ux_surf_cr[1]+flux_rho_l[0]*ux_surf_cl[1]))+1.224744871391589*(pkpm_penalization_rhoux_r[1]+pkpm_penalization_rhoux_l[1])-0.4330127018922193*((ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_x_r[1]+avg_p_ij_x_l[1]))*dx1; 

  outrhou1[0] += ((-0.25*flux_rho_r[1]*uy_surf_rl[1])+0.25*flux_rho_l[1]*uy_surf_lr[1]-0.25*flux_rho_r[1]*uy_surf_cr[1]+0.25*flux_rho_l[1]*uy_surf_cl[1]-0.25*flux_rho_r[0]*uy_surf_rl[0]+0.25*flux_rho_l[0]*uy_surf_lr[0]-0.25*flux_rho_r[0]*uy_surf_cr[0]+0.25*flux_rho_l[0]*uy_surf_cl[0]+0.7071067811865475*pkpm_penalization_rhouy_r[0]-0.7071067811865475*(pkpm_penalization_rhouy_l[0]+avg_p_ij_y_r[0])+0.7071067811865475*avg_p_ij_y_l[0])*dx1; 
  outrhou1[1] += ((-0.25*flux_rho_r[0]*uy_surf_rl[1])+0.25*flux_rho_l[0]*uy_surf_lr[1]-0.25*flux_rho_r[0]*uy_surf_cr[1]+0.25*flux_rho_l[0]*uy_surf_cl[1]+0.7071067811865475*pkpm_penalization_rhouy_r[1]-0.7071067811865475*pkpm_penalization_rhouy_l[1]-0.25*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+0.25*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_y_r[1]+0.7071067811865475*avg_p_ij_y_l[1])*dx1; 
  outrhou1[2] += ((-0.4330127018922193*(flux_rho_r[1]*uy_surf_rl[1]+flux_rho_l[1]*uy_surf_lr[1]+flux_rho_r[1]*uy_surf_cr[1]+flux_rho_l[1]*uy_surf_cl[1]+flux_rho_r[0]*uy_surf_rl[0]+flux_rho_l[0]*uy_surf_lr[0]+flux_rho_r[0]*uy_surf_cr[0]+flux_rho_l[0]*uy_surf_cl[0]))+1.224744871391589*(pkpm_penalization_rhouy_r[0]+pkpm_penalization_rhouy_l[0])-1.224744871391589*(avg_p_ij_y_r[0]+avg_p_ij_y_l[0]))*dx1; 
  outrhou1[3] += ((-0.4330127018922193*(flux_rho_r[0]*uy_surf_rl[1]+flux_rho_l[0]*uy_surf_lr[1]+flux_rho_r[0]*uy_surf_cr[1]+flux_rho_l[0]*uy_surf_cl[1]))+1.224744871391589*(pkpm_penalization_rhouy_r[1]+pkpm_penalization_rhouy_l[1])-0.4330127018922193*((uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_y_r[1]+avg_p_ij_y_l[1]))*dx1; 

  outrhou2[0] += ((-0.25*flux_rho_r[1]*uz_surf_rl[1])+0.25*flux_rho_l[1]*uz_surf_lr[1]-0.25*flux_rho_r[1]*uz_surf_cr[1]+0.25*flux_rho_l[1]*uz_surf_cl[1]-0.25*flux_rho_r[0]*uz_surf_rl[0]+0.25*flux_rho_l[0]*uz_surf_lr[0]-0.25*flux_rho_r[0]*uz_surf_cr[0]+0.25*flux_rho_l[0]*uz_surf_cl[0]+0.7071067811865475*pkpm_penalization_rhouz_r[0]-0.7071067811865475*(pkpm_penalization_rhouz_l[0]+avg_p_ij_z_r[0])+0.7071067811865475*avg_p_ij_z_l[0])*dx1; 
  outrhou2[1] += ((-0.25*flux_rho_r[0]*uz_surf_rl[1])+0.25*flux_rho_l[0]*uz_surf_lr[1]-0.25*flux_rho_r[0]*uz_surf_cr[1]+0.25*flux_rho_l[0]*uz_surf_cl[1]+0.7071067811865475*pkpm_penalization_rhouz_r[1]-0.7071067811865475*pkpm_penalization_rhouz_l[1]-0.25*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+0.25*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_z_r[1]+0.7071067811865475*avg_p_ij_z_l[1])*dx1; 
  outrhou2[2] += ((-0.4330127018922193*(flux_rho_r[1]*uz_surf_rl[1]+flux_rho_l[1]*uz_surf_lr[1]+flux_rho_r[1]*uz_surf_cr[1]+flux_rho_l[1]*uz_surf_cl[1]+flux_rho_r[0]*uz_surf_rl[0]+flux_rho_l[0]*uz_surf_lr[0]+flux_rho_r[0]*uz_surf_cr[0]+flux_rho_l[0]*uz_surf_cl[0]))+1.224744871391589*(pkpm_penalization_rhouz_r[0]+pkpm_penalization_rhouz_l[0])-1.224744871391589*(avg_p_ij_z_r[0]+avg_p_ij_z_l[0]))*dx1; 
  outrhou2[3] += ((-0.4330127018922193*(flux_rho_r[0]*uz_surf_rl[1]+flux_rho_l[0]*uz_surf_lr[1]+flux_rho_r[0]*uz_surf_cr[1]+flux_rho_l[0]*uz_surf_cl[1]))+1.224744871391589*(pkpm_penalization_rhouz_r[1]+pkpm_penalization_rhouz_l[1])-0.4330127018922193*((uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_z_r[1]+avg_p_ij_z_l[1]))*dx1; 

  return 0.;

} 
