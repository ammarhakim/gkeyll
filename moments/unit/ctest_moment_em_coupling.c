// Unit test for moment_em_coupling: the plasma-oscillation problem of
// Wang et al., J. Comput. Phys. 415 (2020) 109510, Section 4.1.
//
// We disregard the spatial (homogeneous) update entirely and integrate ONLY
// the EM source terms in a single cell of a 1D domain — there is no spatial
// variation, so the problem reduces to the constant-coefficient ODE system
// Eqn. (18):
//
//   dJ/dt = omega_p^2 eps0 E + J x Omega,   eps0 dE/dt = -sum_s J_s.
//
// With immobile ions (a single electron fluid), vanishing magnetic field, and
// a perturbation imposed on E_x, the initial state is an eigenvector of
// omega = omega_pe. The system should then oscillate at the plasma frequency
// and exactly conserve the discrete energy invariant Eqn. (19),
//
//   I = sum_s (1/2) J_s^2 / (eps0 omega_s^2) + (eps0/2) E^2,
//
// which the time-centered (implicit-midpoint) solver in
// gkyl_moment_em_coupling_implicit_advance preserves regardless of dt.
//
// Units are normalized so that omega_pe = 1: eps0 = mu0 = 1, electron mass
// m = 1, |charge| = 1, and background mass density rho = 1, giving
// omega_pe^2 = rho (q/m)^2 / eps0 = 1.
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

// Component layout reminders. Five-moment (Euler) fluid:
//   f = [rho, mom_x, mom_y, mom_z, energy].
// Maxwell EM field (8 components):
//   em = [Ex, Ey, Ez, Bx, By, Bz, phi, psi].
enum { F_RHO, F_MX, F_MY, F_MZ, F_E };
enum { EM_EX, EM_EY, EM_EZ, EM_BX, EM_BY, EM_BZ };

// Discrete energy invariant Eqn. (19) for the single electron species, in the
// normalized units where eps0 = omega_pe = 1: I = (1/2)(E_x^2 + J_x^2). The
// current is J_x = (q/m) mom_x.
static double
energy_invariant(double eps0, double omega_p_sq, double q_over_m,
  double Ex, double mom_x)
{
  double Jx = q_over_m * mom_x;
  return 0.5 * eps0 * (Ex * Ex) + 0.5 * (Jx * Jx) / (eps0 * omega_p_sq);
}

// Integrate the plasma oscillation forward by nsteps of size dt, returning the
// final E_x and mom_x and the maximum drift of the energy invariant over the
// whole run. The cosine reference frequency is the modified numerical
// frequency from the scheme's dispersion relation, Eqn. (27):
//   (2/dt) tan(omega_num dt/2) = omega_p  =>  omega_num = (2/dt) atan(omega_p dt/2).
static void
run_plasma_oscillation(double dt, int nsteps, double E0,
  double *Ex_final, double *mom_x_final, double *max_inv_drift,
  double *omega_num)
{
  const double eps0 = 1.0, mu0 = 1.0;
  const double mass = 1.0, charge = -1.0;       // electrons
  const double q_over_m = charge / mass;
  const double rho = 1.0;                        // background mass density
  const double gas_gamma = 5.0 / 3.0;
  const double pressure = 1.0;                   // arbitrary; internal energy is conserved

  const double omega_p_sq = rho * (q_over_m * q_over_m) / eps0;  // = 1

  // Single-cell 1D grid.
  struct gkyl_rect_grid grid;
  double lower[1] = { 0.0 };
  double upper[1] = { 1.0 };
  int cells[1] = { 1 };
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int rlower[1] = { 1 };
  int rupper[1] = { 1 };
  struct gkyl_range range;
  gkyl_range_init(&range, 1, rlower, rupper);

  // Coupling object: one charged Euler fluid, fields free to evolve.
  struct gkyl_moment_em_coupling_inp inp = { 0 };
  inp.grid = &grid;
  inp.nfluids = 1;
  inp.param[0].type = GKYL_EQN_EULER;
  inp.param[0].charge = charge;
  inp.param[0].mass = mass;
  inp.epsilon0 = eps0;
  inp.mu0 = mu0;
  inp.static_field = false;
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(inp);

  // Fluid, EM, and (zeroed) auxiliary source arrays. The implicit advance
  // unconditionally fetches app_accel, p_rhs and nT_sources per species and
  // dereferences ext_em / app_current, so all must be valid (zeroed) arrays.
  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5, range.volume);
  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *app_accel = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);
  struct gkyl_array *p_rhs = gkyl_array_new(GKYL_DOUBLE, 5, range.volume);
  struct gkyl_array *nT_source = gkyl_array_new(GKYL_DOUBLE, 5, range.volume);
  struct gkyl_array *ext_em = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *app_current = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(em, 0.0);
  gkyl_array_clear(app_accel, 0.0);
  gkyl_array_clear(p_rhs, 0.0);
  gkyl_array_clear(nT_source, 0.0);
  gkyl_array_clear(ext_em, 0.0);
  gkyl_array_clear(app_current, 0.0);

  long cidx = gkyl_range_idx(&range, rlower);

  // Initial condition: stationary, uniform plasma (mom = 0), B = 0, with a
  // perturbation E_x = E0. This is an eigenvector of omega = omega_pe.
  double *f0 = gkyl_array_fetch(fluid, cidx);
  f0[F_RHO] = rho;
  f0[F_MX] = 0.0; f0[F_MY] = 0.0; f0[F_MZ] = 0.0;
  f0[F_E] = pressure / (gas_gamma - 1.0);  // mom = 0 => no kinetic part

  double *e0 = gkyl_array_fetch(em, cidx);
  e0[EM_EX] = E0;

  double internal_energy_0 = f0[F_E];
  double inv0 = energy_invariant(eps0, omega_p_sq, q_over_m, e0[EM_EX], f0[F_MX]);
  *max_inv_drift = 0.0;

  // Array-of-pointers views expected by the advance API.
  struct gkyl_array *fluids[GKYL_MAX_SPECIES] = { fluid };
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES] = { app_accel };
  const struct gkyl_array *p_rhss[GKYL_MAX_SPECIES] = { p_rhs };
  const struct gkyl_array *nT_sources[GKYL_MAX_SPECIES] = { nT_source };

  double t_curr = 0.0;
  for (int n = 0; n < nsteps; n++) {
    gkyl_moment_em_coupling_implicit_advance(mom_em, t_curr, dt, &range,
      fluids, app_accels, p_rhss, em, app_current, ext_em, nT_sources);
    t_curr += dt;

    const double *f = gkyl_array_cfetch(fluid, cidx);
    const double *e = gkyl_array_cfetch(em, cidx);

    double inv = energy_invariant(eps0, omega_p_sq, q_over_m, e[EM_EX], f[F_MX]);
    double drift = fabs(inv - inv0) / inv0;
    if (drift > *max_inv_drift) *max_inv_drift = drift;
  }

  const double *f = gkyl_array_cfetch(fluid, cidx);
  const double *e = gkyl_array_cfetch(em, cidx);
  *Ex_final = e[EM_EX];
  *mom_x_final = f[F_MX];

  // The Euler energy update keeps the internal (thermal) energy fixed while
  // exchanging only kinetic and EM energy, so the thermal part must be
  // unchanged: internal = energy - (1/2) mom^2 / rho.
  double internal_energy = f[F_E] - 0.5 * (f[F_MX]*f[F_MX] + f[F_MY]*f[F_MY]
    + f[F_MZ]*f[F_MZ]) / f[F_RHO];
  TEST_CHECK(gkyl_compare(internal_energy, internal_energy_0, 1e-12));

  *omega_num = (2.0 / dt) * atan(sqrt(omega_p_sq) * dt / 2.0);

  gkyl_array_release(fluid);
  gkyl_array_release(em);
  gkyl_array_release(app_accel);
  gkyl_array_release(p_rhs);
  gkyl_array_release(nT_source);
  gkyl_array_release(ext_em);
  gkyl_array_release(app_current);
  gkyl_moment_em_coupling_release(mom_em);
}

// Small time step: the solver should reproduce the analytic oscillation
// E_x(t) = E0 cos(omega_num t) and conserve the energy invariant to roundoff.
static void
test_plasma_oscillation_small_dt()
{
  const double E0 = 1.0;
  const double dt = 0.1;   // = 0.1 / omega_pe, as in Fig. 1 of the paper
  const int nsteps = 200;  // ~3 plasma-oscillation periods

  double Ex, mom_x, max_inv_drift, omega_num;
  run_plasma_oscillation(dt, nsteps, E0, &Ex, &mom_x, &max_inv_drift, &omega_num);

  // Energy invariant is conserved to machine precision by the implicit
  // midpoint scheme (Section 3.5).
  TEST_CHECK(max_inv_drift < 1e-12);
  TEST_MSG("max energy-invariant drift = %g", max_inv_drift);

  // E_x tracks the analytic oscillation at the modified numerical frequency.
  double Ex_expected = E0 * cos(omega_num * dt * nsteps);
  TEST_CHECK(gkyl_compare(Ex, Ex_expected, 1e-10));
  TEST_MSG("Ex = %g, expected %g", Ex, Ex_expected);
}

// Large time step: the scheme is NOT restricted by the plasma frequency (Section 3.4).
// With dt = 10/omega_pe the  oscillation is aliased, but the run must stay finite 
// and still conserve the energy invariant exactly.
static void
test_plasma_oscillation_large_dt()
{
  const double E0 = 1.0;
  const double dt = 10.0;  // = 10 / omega_pe, far beyond an explicit stability limit
  const int nsteps = 100;

  double Ex, mom_x, max_inv_drift, omega_num;
  run_plasma_oscillation(dt, nsteps, E0, &Ex, &mom_x, &max_inv_drift, &omega_num);

  TEST_CHECK(max_inv_drift < 1e-12);
  TEST_MSG("max energy-invariant drift = %g", max_inv_drift);

  // No instability: the field amplitude is bounded by the (conserved) initial
  // energy, |Ex| <= E0 (up to roundoff).
  TEST_CHECK(fabs(Ex) <= E0 + 1e-10);
  TEST_MSG("Ex = %g", Ex);
}

TEST_LIST = {
  { "plasma_oscillation_small_dt", test_plasma_oscillation_small_dt },
  { "plasma_oscillation_large_dt", test_plasma_oscillation_large_dt },
  { NULL, NULL },
};
