#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_calc.h>

struct gkyl_dg_updater_moment {
  enum gkyl_model_id model_id; // Identifier for model (e.g., SR, see gkyl_eqn_type.h)
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_mom_type *type; // Moment type
  struct gkyl_mom_calc *up_moment; // Updater for computing moment

  double moment_tm; // total time spent in computing moment
};

struct gkyl_dg_updater_moment_tm {
  double moment_tm; // time for moment updates
};