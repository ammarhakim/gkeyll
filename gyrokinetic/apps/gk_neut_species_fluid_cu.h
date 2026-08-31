#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

#ifdef GKYL_HAVE_CUDA
void gk_neut_fluid_diffusion_dirichlet_ghost_cu(int cdim, int dir,
  int skin_idx_dir, double mass, double number_density,
  const struct gkyl_range *local_ext, const struct gkyl_range *ghost,
  const struct gkyl_array *number_density_field,
  struct gkyl_array *density, struct gkyl_array *tensor);

void gk_neut_fluid_diffusion_floor_cu(int cdim, double nu_min,
  const struct gkyl_range *range, const struct gkyl_array *vt_sq,
  const struct gkyl_array *nu, struct gkyl_array *D);

void gk_neut_fluid_diffusion_vec_set_cu(const struct gkyl_range *range,
  struct gkyl_array *out, double a, const struct gkyl_array *x,
  double b, const struct gkyl_array *y);

void gk_neut_fluid_diffusion_p_update_cu(const struct gkyl_range *range,
  struct gkyl_array *p, const struct gkyl_array *r,
  const struct gkyl_array *v, double beta, double omega);

double gk_neut_fluid_diffusion_dot_cu(const struct gkyl_range *range,
  const struct gkyl_array *a, const struct gkyl_array *b, double *dot_cu);

void gk_neut_fluid_diffusion_positivity_cu(int cdim,
  const struct gkyl_range *range, struct gkyl_array *rho);
#endif
