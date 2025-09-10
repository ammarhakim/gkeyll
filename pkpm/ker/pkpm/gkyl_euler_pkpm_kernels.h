#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_wv_ten_moment_priv.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_1x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_1x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_1x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_1x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_1x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p1(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p1(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p1(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p3(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_3x_ser_p3(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_3x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_3x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_3x_ser_p3(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_3x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_3x_ser_p3(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p3(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p3(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p3(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_z_3x_ser_p3(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_z_3x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p3(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_p_force_2x_tensor_p2(const double *prim_c, const double *div_b, double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH int pkpm_vars_u_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm); 
GKYL_CU_DH void pkpm_vars_u_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u); 
GKYL_CU_DH int pkpm_vars_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf); 
GKYL_CU_DH void pkpm_vars_integrated_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT pkpm_int_vars); 
GKYL_CU_DH void pkpm_vars_io_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* p_ij, const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 
GKYL_CU_DH void euler_pkpm_em_coupling_set_2x_tensor_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT em); 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_2x_tensor_p2(int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double* pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) ; 
GKYL_CU_DH void euler_pkpm_source_2x_tensor_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_pkpm_vol_2x_tensor_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_tensor_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_x_2x_tensor_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_x_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_tensor_p2(const double *dxv, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *prim_c, const double *bvar_c, const double *nu_c, 
    double* GKYL_RESTRICT pkpm_accel); 
GKYL_CU_DH void pkpm_vars_penalization_y_2x_tensor_p2(double tol, bool force_lax, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r, 
    const double *p_ij_l, const double *p_ij_r, 
    const double *prim_l, const double *prim_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_r, 
    double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 
GKYL_CU_DH void euler_pkpm_limiter_y_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r); 
GKYL_CU_DH double euler_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r, 
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax_l, const double *pkpm_lax_r, 
    const double *pkpm_penalization_l, const double *pkpm_penalization_r,  
    double* GKYL_RESTRICT out); 

EXTERN_C_END 
