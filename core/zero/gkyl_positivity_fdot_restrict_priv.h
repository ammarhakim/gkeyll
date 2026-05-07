#pragma once

#include <gkyl_positivity_fdot_restrict.h>

struct gkyl_positivity_fdot_restrict {
  struct gkyl_basis basis;
  enum gkyl_positivity_fdot_restrict_type type;
  double safety_factor;
};
