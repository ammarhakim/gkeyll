#pragma once

#include <assert.h>

#include <gkyl_positivity_shift_zhang_shu.h>

struct gkyl_positivity_shift_zhang_shu {
  struct gkyl_basis basis;
  double cellav_fac;
};