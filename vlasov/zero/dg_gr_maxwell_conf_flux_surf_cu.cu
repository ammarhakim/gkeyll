/* -*- c++ -*- */

#include <cassert>
#include <cstdlib>

extern "C" {
#include <gkyl_dg_gr_maxwell_conf_flux_surf.h>
}

struct gkyl_dg_gr_maxwell_conf_flux_surf*
gkyl_dg_gr_maxwell_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_gr_maxwell_conf_flux_surf_inp *inp)
{
  assert(!"GR DG Maxwell GPU configuration-space flux updater is incomplete");
  abort();
  return 0;
}

void
gkyl_dg_gr_maxwell_conf_flux_surf_advance_cu(struct gkyl_dg_gr_maxwell_conf_flux_surf *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext,
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift,
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *det_h,
  const struct gkyl_array *field_con, const struct gkyl_array *field_no_J_con, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf)
{
  assert(!"GR DG Maxwell GPU configuration-space flux advance is incomplete");
  abort();
}
