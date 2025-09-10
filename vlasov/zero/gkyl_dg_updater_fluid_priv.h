#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_fluid {
  enum gkyl_eqn_type eqn_id; // Identifier for equation type (e.g., Euler, see gkyl_eqn_type.h)
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_dg_eqn *eqn_fluid; // Equation object
  struct gkyl_hyper_dg *up_fluid; // solvers for specific fluid equation

  double fluid_tm; // total time spent in computing fluid equation
};
