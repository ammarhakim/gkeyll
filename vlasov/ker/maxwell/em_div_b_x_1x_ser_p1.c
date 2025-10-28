#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_div_b_x_1x_ser_p1(const double *dxv, 
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector in left/center/right cells in each direction. 
  // bvar_c:          Input volume expansion of magnetic field unit vector in center cell. 
  // max_b:           Output surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|. 
  // div_b:           Output volume expansion of div(b).

  const double dx1 = 2.0/dxv[0]; 
  const double *b_c = &bvar_c[0]; 
  const double *b_surf_lr = &bvar_surf_l[1]; 
  const double *b_surf_cl = &bvar_surf_c[0]; 
  const double *b_surf_cr = &bvar_surf_c[1]; 
  const double *b_surf_rl = &bvar_surf_r[0]; 

  double *max_b_l = &max_b[0]; 
  double *max_b_r = &max_b[1]; 

  max_b_l[0] = fmax(fabs(b_surf_lr[0]), fabs(b_surf_cl[0])); 
  max_b_r[0] = fmax(fabs(b_surf_cr[0]), fabs(b_surf_rl[0])); 

  div_b[0] += (0.3535533905932737*b_surf_rl[0]-0.3535533905932737*b_surf_lr[0]+0.3535533905932737*b_surf_cr[0]-0.3535533905932737*b_surf_cl[0])*dx1; 
  div_b[1] += (0.6123724356957944*(b_surf_rl[0]+b_surf_lr[0]+b_surf_cr[0]+b_surf_cl[0])-1.732050807568877*b_c[0])*dx1; 

} 
