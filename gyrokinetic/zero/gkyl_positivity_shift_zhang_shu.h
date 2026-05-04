#pragma once

#include <stdbool.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

typedef struct gkyl_positivity_shift_zhang_shu gkyl_positivity_shift_zhang_shu;

struct gkyl_positivity_shift_zhang_shu*
gkyl_positivity_shift_zhang_shu_new(struct gkyl_basis basis);

bool
gkyl_positivity_shift_zhang_shu_advance(gkyl_positivity_shift_zhang_shu* up,
  const struct gkyl_range *range, struct gkyl_array *f);

void
gkyl_positivity_shift_zhang_shu_release(gkyl_positivity_shift_zhang_shu* up);