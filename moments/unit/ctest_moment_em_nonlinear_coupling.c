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

enum {F_RHO, F_MX, F_MY, F_MZ, F_E};
enum {EM_EX, EM_EY, EM_EZ, EM_BX, EM_BY, EM_BZ};

static void
run_plasma_oscillation(double dt, int nsteps, double m0, 
                       double *Ex_final, double *mom_x_final, 
                       double *omega_num)
{
  TEST_MSG("Hello");
  const double eps0 = 1.0, mu = 1.0;
  const double mass = 1.0, charge = -1.0;
  const double q_over_m = charge / mass;
  const double rho = 1.0;
  const double gas_gamma = 5.0 / 3.0;
  const double pressure = 1.0;

  const double omega_p_sq = rho * (q_over_m * q_over_m) / eps0;

  // Single Cell 1D grid;
  
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
  inp.mu0 = mu;
  inp.use_rel = true;
  inp.static_field = false;
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(inp);
  
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

  double *f0 = gkyl_array_fetch(fluid, cidx);
  f0[F_RHO] = rho;
  f0[F_MX] = m0; f0[F_MY] = 0.0; f0[F_MZ] = 0.0;
  f0[F_E] = pressure / (gas_gamma - 1.0);

  double *e0 = gkyl_array_fetch(em, cidx);
  e0[EM_EX] = 0.0;

  double internal_energy_0 = f0[F_E];
  
  struct gkyl_array *fluids[GKYL_MAX_SPECIES] = { fluid };
  struct gkyl_array *app_accels[GKYL_MAX_SPECIES] = { app_accel };
  struct gkyl_array *p_rhss[GKYL_MAX_SPECIES] = { p_rhs };
  struct gkyl_array *nT_sources[GKYL_MAX_SPECIES] = { nT_source };

  double t_curr = 0.0;
  for (int n = 0; n < nsteps; n++) {
    gkyl_moment_em_coupling_implicit_advance(mom_em, t_curr, dt, &range,
                                             fluids, app_accels, p_rhss, em,
                                             app_current, ext_em, nT_sources);
    t_curr += dt;
    const double *f = gkyl_array_cfetch(fluid, cidx);
    const double *e = gkyl_array_cfetch(em, cidx);
  }

  const double *f = gkyl_array_cfetch(fluid, cidx);
  const double *e = gkyl_array_cfetch(em, cidx);
  *Ex_final = e[EM_EX];
  *mom_x_final = f[F_MX];  
}

static void
test_plasma_oscillation_small_dt()
{
  const double m0 = 5.0;
  const double dt = 0.1;
  const int nsteps = 200;
  
  double Ex, mom_x, omega_num;
  run_plasma_oscillation(dt, nsteps, m0, &Ex, &mom_x, &omega_num);

  double Ex_expected = 0.0;
  TEST_CHECK(gkyl_compare(Ex, Ex_expected, 1e-10));
}


static void
test_plasma_oscillation_large_dt()
{
  const double m0 = 5.0;
  const double dt = 10;
  const int nsteps = 200;
  
  double Ex, mom_x, omega_num;
  run_plasma_oscillation(dt, nsteps, m0, &Ex, &mom_x, &omega_num);

  double Ex_expected = 0.0;
  TEST_CHECK(fabs(Ex) <= m0 + 1e-10);
  TEST_MSG("Ex = %g", Ex);
}


TEST_LIST = {
  { "plasma_oscillation_small_dt", test_plasma_oscillation_small_dt },
  { "plasma_oscillation_large_dt", test_plasma_oscillation_large_dt },
  { NULL, NULL },
};
