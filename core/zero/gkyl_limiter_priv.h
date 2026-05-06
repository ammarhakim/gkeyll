#pragma once

#include <assert.h>

#include <gkyl_limiter.h>

typedef void (*limit_cell_func_t)(const struct gkyl_limiter *up,
  const struct gkyl_range *range, struct gkyl_array *f, const int *idx, double *fc);

typedef void (*limit_timestep_func_t)(const struct gkyl_limiter *up, const double *fc,
  struct gkyl_array *dfdt, double *dt, long lidx, double *dt_bound);

typedef void (*limit_func_t)(const struct gkyl_limiter *up,
  const struct gkyl_range *range, struct gkyl_array *f, struct gkyl_array *dfdt, double *dt);

struct gkyl_limiter {
  enum gkyl_limiter_type type;
  struct gkyl_basis basis;
  double cellav_fac;
  double dt_factor;
  double *fquad;
  limit_func_t limiter_func;
  limit_cell_func_t limit_cell_func;
  limit_timestep_func_t limit_timestep_func;
};