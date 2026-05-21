/* -*- c++ -*- */

#include <cassert>
#include <cstdlib>

extern "C" {
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
}

void
gkyl_dg_gr_maxwell_divide_Jc_cu(const struct gkyl_basis *conf_basis,
  const struct gkyl_range *conf_range, const struct gkyl_array *det_h,
  const struct gkyl_array *field_with_J_con, struct gkyl_array *field_no_J_con)
{
  assert(!"GR DG Maxwell GPU divide_Jc is incomplete");
  abort();
}

void
gkyl_dg_gr_maxwell_rescale_Jc_cu(const struct gkyl_basis *conf_basis,
  const struct gkyl_range *conf_range, const struct gkyl_array *det_h,
  const struct gkyl_array *field_no_J_con, struct gkyl_array *field_with_J_con)
{
  assert(!"GR DG Maxwell GPU rescale_Jc is incomplete");
  abort();
}
