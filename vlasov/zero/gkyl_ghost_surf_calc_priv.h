#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_ghost_surf_calc {
  struct gkyl_rect_grid grid; // grid object
  const struct gkyl_dg_eqn *equation; // equation object
  int cdim; // number of configuration space dimensions

  uint32_t flags;
  struct gkyl_ghost_surf_calc *on_dev; // pointer to itself or device data
};
