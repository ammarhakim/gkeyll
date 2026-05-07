#pragma once

#include <gkyl_positivity_fdot_restrict.h>

struct gkyl_positivity_fdot_restrict {
  struct gkyl_basis basis;
  enum gkyl_positivity_fdot_restrict_mode mode;
  double safety_factor;
  double *fquad; // Temporary array for quadrature values.
  int num_basis; // Number of basis functions.
};
