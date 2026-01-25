#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p1(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf) 
{
  // drag_diff_coeff_corrs: Corrections to be added to coeffs, function of config space only. 
  // drag_coeff: FPO drag coefficient. 
  // diff_coeff: FPO diffusion coefficient. 
 
  // Index into drag and diffusion coefficients. 
  double *ax = &drag_coeff[0]; 
  double *ay = &drag_coeff[40]; 
  double *az = &drag_coeff[80]; 
  double *Dxx = &diff_coeff[0]; 
  double *Dyy = &diff_coeff[160]; 
  double *Dzz = &diff_coeff[320]; 
 
  // Index into surface expansions. 
  double *ax_surf = &drag_coeff_surf[0]; 
  double *ay_surf = &drag_coeff_surf[8]; 
  double *az_surf = &drag_coeff_surf[16]; 
  double *Dxx_surf = &diff_coeff_surf[0]; 
  double *Dyy_surf = &diff_coeff_surf[64]; 
  double *Dzz_surf = &diff_coeff_surf[128]; 
 
  // Index into correction array. 
  const double* ax_corr = &drag_diff_coeff_corrs[0]; 
  const double* ay_corr = &drag_diff_coeff_corrs[2]; 
  const double* az_corr = &drag_diff_coeff_corrs[4]; 
  const double* D_corr = &drag_diff_coeff_corrs[6]; 
 
 
 
 
 
 
 
 
 
 
 
 
} 

