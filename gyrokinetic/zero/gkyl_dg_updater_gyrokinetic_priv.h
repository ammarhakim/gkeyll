#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_gyrokinetic.h>

struct gkyl_dg_updater_gyrokinetic {
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_dg_eqn *eqn_gyrokinetic; // Equation object
  struct gkyl_hyper_dg *up_gyrokinetic; // solvers for specific gyrokinetic equation

  double gyrokinetic_tm; // total time spent in computing gyrokinetic equation
};
