#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_gk_anomalous_diffusion {
  struct gkyl_dg_eqn *dgeqn; // Equation object.
  struct gkyl_hyper_dg *hyperdg; // Solvers for specific diffusion equation.
  bool use_gpu;

  double diffusion_tm; // Total time spent in computing diffusion equation.
};
