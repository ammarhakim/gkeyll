#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

typedef struct gkyl_limiter gkyl_limiter;

enum gkyl_limiter_type {
  GKYL_LIMITER_ZS = 0, // Zhang & Shu limiter.
  GKYL_LIMITER_MRS, // Moe, Rossmanith, and Seal limiter (https://arxiv.org/pdf/1507.03024)
};

struct gkyl_limiter_inp {
  struct gkyl_basis basis;
  enum gkyl_limiter_type type;
  double dt_factor; // safety factor used when computing timestep restriction
};

struct gkyl_limiter*
gkyl_limiter_new(struct gkyl_limiter_inp inp);

void
gkyl_limiter_advance(gkyl_limiter* up, const struct gkyl_range *range, struct gkyl_array *f);

void
gkyl_limiter_advance_timestep(gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt);

void
gkyl_limiter_release(gkyl_limiter* up);
