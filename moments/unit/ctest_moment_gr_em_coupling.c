// Unit test for the explicit (special-)relativistic multi-fluid + Maxwell
// coupling, gkyl_moment_em_coupling_gr_em_explicit_advance.
//
// As in ctest_moment_em_coupling.c we disregard the spatial (homogeneous)
// update and integrate ONLY the EM source terms in a single cell. With a
// single mobile electron fluid, immobile background, and vanishing magnetic
// field, a small perturbation of E_x sets up a relativistic Langmuir
// oscillation at the (relativistic) plasma frequency
//
//   omega_p^2 = (q/m)^2 rho / (eps0 h),
//
// where h = 1 + (p/rho) gamma/(gamma-1) is the specific enthalpy. The
// explicit SSP-RK3 solver should reproduce E_x(t) = E0 cos(omega_p t) to the
// integrator's (small-dt) accuracy, and conserve D exactly. Separately we
// check the stability-control return: at a dt that under-resolves omega_p the
// advance must report success=false with a sensible dt_suggested.
//
// Units are normalized: eps0 = mu0 = 1, electron mass m = 1, |charge| = 1,
// rest-mass density rho = 1.
#include <acutest.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Modular relativistic-Euler hydro state f = [D, S_x, S_y, S_z, tau], and the
// 8-component Maxwell field em = [E_x, E_y, E_z, B_x, B_y, B_z, phi, psi].
enum { F_D, F_SX, F_SY, F_SZ, F_TAU };
enum { EM_EX, EM_EY, EM_EZ, EM_BX, EM_BY, EM_BZ };

// Build the modular GR-Euler conserved state for a fluid at rest (v = 0, W = 1)
// with rest density rho and pressure p (ideal EOS, index gas_gamma):
//   D = rho, S_i = 0, tau = rho(h-1) - p = p/(gas_gamma-1).
static void
set_rest_state(double *f, double rho, double p, double gas_gamma)
{
  f[F_D] = rho;
  f[F_SX] = 0.0; f[F_SY] = 0.0; f[F_SZ] = 0.0;
  f[F_TAU] = p / (gas_gamma - 1.0);
}

// Integrate the relativistic plasma oscillation forward by nsteps of size dt,
// returning the final E_x and S_x, the maximum |D - D0| drift, and (from the
// first step) the reported omega_max / success.
static void
run_plasma_oscillation(double rho, double p, double gas_gamma, double E0,
  double dt, int nsteps, double *Ex_final, double *Sx_final, double *D_drift,
  double *omega_max_first, bool *success_first)
{
  const double eps0 = 1.0, mu0 = 1.0;
  const double mass = 1.0, charge = -1.0; // electrons

  struct gkyl_rect_grid grid;
  double lower[1] = { 0.0 }, upper[1] = { 1.0 };
  int cells[1] = { 1 };
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int rlower[1] = { 1 }, rupper[1] = { 1 };
  struct gkyl_range range;
  gkyl_range_init(&range, 1, rlower, rupper);

  struct gkyl_moment_em_coupling_inp inp = { 0 };
  inp.grid = &grid;
  inp.nfluids = 1;
  inp.param[0].type = GKYL_EQN_GR_EULER_TETRAD;
  inp.param[0].charge = charge;
  inp.param[0].mass = mass;
  inp.epsilon0 = eps0;
  inp.mu0 = mu0;
  inp.has_gr_em_coupling = true;
  inp.gr_em_eos[0] = gkyl_gr_euler_eos_ideal(gas_gamma);
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(inp);

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5, range.volume);
  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(em, 0.0);

  long cidx = gkyl_range_idx(&range, rlower);

  double *f0 = gkyl_array_fetch(fluid, cidx);
  set_rest_state(f0, rho, p, gas_gamma);
  double D0 = f0[F_D];

  double *e0 = gkyl_array_fetch(em, cidx);
  e0[EM_EX] = E0;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES] = { fluid };

  *D_drift = 0.0;
  double t_curr = 0.0;
  for (int n = 0; n < nsteps; n++) {
    struct gkyl_moment_em_coupling_status status =
      gkyl_moment_em_coupling_gr_em_explicit_advance(mom_em, t_curr, dt, &range, fluids, em, NULL);
    if (n == 0) {
      *omega_max_first = status.omega_max;
      *success_first = status.success;
    }
    t_curr += dt;

    const double *f = gkyl_array_cfetch(fluid, cidx);
    double drift = fabs(f[F_D] - D0);
    if (drift > *D_drift) *D_drift = drift;
  }

  const double *f = gkyl_array_cfetch(fluid, cidx);
  const double *e = gkyl_array_cfetch(em, cidx);
  *Ex_final = e[EM_EX];
  *Sx_final = f[F_SX];

  gkyl_array_release(fluid);
  gkyl_array_release(em);
  gkyl_moment_em_coupling_release(mom_em);
}

// Small time step: E_x should track the analytic oscillation E0 cos(omega_p t)
// over half and full periods, D is exactly conserved, and the solver reports a
// stable (success=true) step.
static void
test_plasma_oscillation_small_dt()
{
  const double rho = 1.0, p = 1.0e-4, gas_gamma = 5.0 / 3.0;
  const double E0 = 1.0e-3; // small => linear oscillation

  const double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));
  const double omega_p = sqrt((1.0) * rho / (1.0 * h)); // (q/m)^2 = 1, eps0 = 1
  const double period = 2.0 * M_PI / omega_p;

  // Half period: E_x should reach -E0.
  {
    const int nsteps = 100;
    const double dt = 0.5 * period / nsteps;
    double Ex, Sx, D_drift, omega_max;
    bool success;
    run_plasma_oscillation(rho, p, gas_gamma, E0, dt, nsteps, &Ex, &Sx, &D_drift, &omega_max, &success);

    TEST_CHECK(success);
    TEST_CHECK(gkyl_compare(omega_max, omega_p, 1.0e-6));
    TEST_MSG("omega_max = %g, omega_p = %g", omega_max, omega_p);

    TEST_CHECK(fabs(Ex - (-E0)) < 1.0e-2 * E0);
    TEST_MSG("half-period Ex = %g, expected %g", Ex, -E0);

    TEST_CHECK(D_drift < 1.0e-14);
    TEST_MSG("D drift = %g", D_drift);
  }

  // Full period: E_x should return to +E0 and momentum to ~0.
  {
    const int nsteps = 200;
    const double dt = period / nsteps;
    double Ex, Sx, D_drift, omega_max;
    bool success;
    run_plasma_oscillation(rho, p, gas_gamma, E0, dt, nsteps, &Ex, &Sx, &D_drift, &omega_max, &success);

    TEST_CHECK(success);
    TEST_CHECK(fabs(Ex - E0) < 1.0e-2 * E0);
    TEST_MSG("full-period Ex = %g, expected %g", Ex, E0);

    TEST_CHECK(fabs(Sx) < 1.0e-2 * E0);
    TEST_MSG("full-period Sx = %g, expected ~0", Sx);
  }
}

// Large time step: the explicit solver is conditionally stable. With dt far
// beyond the plasma-frequency limit the advance must report success=false and a
// dt_suggested that brings omega_max * dt back to the stability bound.
static void
test_stability_control_large_dt()
{
  const double eps0 = 1.0, mu0 = 1.0;
  const double mass = 1.0, charge = -1.0;
  const double rho = 1.0, p = 1.0e-4, gas_gamma = 5.0 / 3.0;

  struct gkyl_rect_grid grid;
  double lower[1] = { 0.0 }, upper[1] = { 1.0 };
  int cells[1] = { 1 };
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int rlower[1] = { 1 }, rupper[1] = { 1 };
  struct gkyl_range range;
  gkyl_range_init(&range, 1, rlower, rupper);

  struct gkyl_moment_em_coupling_inp inp = { 0 };
  inp.grid = &grid;
  inp.nfluids = 1;
  inp.param[0].type = GKYL_EQN_GR_EULER_TETRAD;
  inp.param[0].charge = charge;
  inp.param[0].mass = mass;
  inp.epsilon0 = eps0;
  inp.mu0 = mu0;
  inp.has_gr_em_coupling = true;
  inp.gr_em_eos[0] = gkyl_gr_euler_eos_ideal(gas_gamma);
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(inp);

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5, range.volume);
  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(em, 0.0);

  long cidx = gkyl_range_idx(&range, rlower);
  double *f0 = gkyl_array_fetch(fluid, cidx);
  set_rest_state(f0, rho, p, gas_gamma);
  double *e0 = gkyl_array_fetch(em, cidx);
  e0[EM_EX] = 1.0e-3;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES] = { fluid };

  const double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));
  const double omega_p = sqrt(rho / h);

  // dt chosen so omega_p * dt >> stability bound.
  const double dt = 100.0;
  struct gkyl_moment_em_coupling_status status =
    gkyl_moment_em_coupling_gr_em_explicit_advance(mom_em, 0.0, dt, &range, fluids, em, NULL);

  TEST_CHECK(!status.success);
  TEST_MSG("expected failure at dt = %g (omega_p dt = %g)", dt, omega_p * dt);

  TEST_CHECK(gkyl_compare(status.omega_max, omega_p, 1.0e-6));

  // dt_suggested must satisfy omega_max * dt_suggested == stability bound and be
  // a genuine reduction from the attempted dt.
  TEST_CHECK(status.dt_suggested < dt);
  TEST_CHECK(status.omega_max * status.dt_suggested <= 1.5 + 1.0e-12);
  TEST_MSG("dt_suggested = %g, omega_max dt_suggested = %g",
    status.dt_suggested, status.omega_max * status.dt_suggested);

  gkyl_array_release(fluid);
  gkyl_array_release(em);
  gkyl_moment_em_coupling_release(mom_em);
}

TEST_LIST = {
  { "plasma_oscillation_small_dt", test_plasma_oscillation_small_dt },
  { "stability_control_large_dt", test_stability_control_large_dt },
  { NULL, NULL },
};
