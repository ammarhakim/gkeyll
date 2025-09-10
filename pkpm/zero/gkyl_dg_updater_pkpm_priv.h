#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_pkpm {
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_dg_eqn *eqn_vlasov; // Equation object for kinetic equations in PKPM
  struct gkyl_dg_eqn *eqn_fluid; // Equation object for fluid equations in PKPM
  struct gkyl_hyper_dg *up_vlasov; // solvers for kinetic equations in PKPM
  struct gkyl_hyper_dg *up_fluid; // solvers for fluid equations in PKPM

  double vlasov_tm; // total time spent in computing kinetic equations
  double fluid_tm; // total time spent in computing fluid equations
};
