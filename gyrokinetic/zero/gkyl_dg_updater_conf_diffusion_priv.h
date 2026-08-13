#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_hyper_dg.h>

struct gkyl_dg_updater_conf_diffusion {
  struct gkyl_dg_eqn *dgeqn;
  struct gkyl_hyper_dg *hyperdg;
  double diffusion_tm;
};
