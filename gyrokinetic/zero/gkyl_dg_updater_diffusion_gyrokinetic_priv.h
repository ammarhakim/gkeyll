#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_diffusion_gyrokinetic {
  struct gkyl_dg_eqn *dgeqn; // Equation object.
  struct gkyl_hyper_dg *hyperdg; // solvers for specific diffusion equation.
  bool use_gpu;

  double diffusion_tm; // total time spent in computing diffusion equation.
};
