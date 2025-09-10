#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_rad_vlasov.h>

struct gkyl_dg_updater_rad_vlasov {
  struct gkyl_dg_eqn *rad_drag; // Radiation drag equation
  struct gkyl_hyper_dg *drag; // solver for radiation drag term
  bool use_gpu;

  double drag_tm; // total time spent in computing radiation drag term
};
