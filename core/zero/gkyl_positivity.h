#pragma once

#include <stdbool.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

typedef struct gkyl_positivity gkyl_positivity;

enum gkyl_positivity_type {
  GKYL_POSITIVITY_NONE = 0,
  GKYL_POSITIVITY_ZS, // Zhang & Shu limiter.
  GKYL_POSITIVITY_MRS, // Moe, Rossmanith, and Seal limiter (https://arxiv.org/pdf/1507.03024)
  GKYL_POSITIVITY_TIMESTEP_AVG, // Only apply timestep limiter based on cell-average value.
  GKYL_POSITIVITY_TIMESTEP_QUAD, // Only apply timestep limiter based on quadrature point values.
};

struct gkyl_positivity_inp {
  struct gkyl_basis basis;
  enum gkyl_positivity_type type;
  double dt_factor; // safety factor used when computing timestep restriction
};

struct gkyl_positivity*
gkyl_positivity_new(struct gkyl_positivity_inp inp);

void
gkyl_positivity_advance(gkyl_positivity* up, const struct gkyl_range *range, struct gkyl_array *f);

void
gkyl_positivity_advance_timestep(gkyl_positivity* up, const struct gkyl_range *range,
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt);

void
gkyl_positivity_release(gkyl_positivity* up);
