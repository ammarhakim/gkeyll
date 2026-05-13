// Unit tests for moment_spacetime_coupling. Strategy:
//
// (1) derive_products: build a small Cartesian grid, run derive_products,
//     and for each cell verify the products row matches a direct per-cell
//     evaluation of fill_products_analytic at the same (x,y,z). This is
//     the contract between the grid-level driver and the per-cell helper.
//
// (2) Source-step equivalence: for a fixed (x,y,z) on each of the three
//     reference spacetimes, build a packed 71-component q and a matching
//     mod 5-component q + products row from identical spacetime callbacks.
//     Run packed explicit_gr_euler_source_update_euler and mod
//     gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler. Verify the
//     hydro slice matches to <1e-12.
//
// (3) SSP-RK3 wrapper equivalence: same setup as (2) but driven through
//     the full advance routine on a gkyl_array, compared to packed
//     explicit_gr_euler_source_update.
#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_spacetime_coupling.h>
#include <gkyl_moment_spacetime_coupling_priv.h>
#include <gkyl_moment_spacetime_products.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_util.h>

// Shared helper used by all the source-step tests. Returns a packed 71-comp
// q[] populated from the spacetime callbacks at (x,y,z) plus the supplied
// hydro state; also writes the matching mod products row.
static void
build_packed_and_products(struct gkyl_gr_spacetime *spacetime,
  double x, double y, double z,
  double rho, double u, double v, double w, double p,
  double gas_gamma,
  double q_packed[71], double *prods)
{
  // Use the same per-cell helper the coupling object uses internally.
  gkyl_moment_spacetime_coupling_fill_products_analytic(spacetime,
    0.0, x, y, z, prods);

  // Now pack q_packed[5..66] from the same products values (the mapping is
  // a direct subrange copy by design: prods[k] == packed-q[5+k] for k<62).
  q_packed[5]  = prods[GKYL_GR_SP_LAPSE];
  q_packed[6]  = prods[GKYL_GR_SP_SHIFT + 0];
  q_packed[7]  = prods[GKYL_GR_SP_SHIFT + 1];
  q_packed[8]  = prods[GKYL_GR_SP_SHIFT + 2];
  for (int i = 0; i < 9; i++) q_packed[9 + i]  = prods[GKYL_GR_SP_GIJ + i];
  for (int i = 0; i < 9; i++) q_packed[18 + i] = prods[GKYL_GR_SP_KIJ + i];
  q_packed[27] = prods[GKYL_GR_SP_EXCISION];
  for (int i = 0; i < 3; i++)  q_packed[28 + i] = prods[GKYL_GR_SP_DALPHA + i];
  for (int i = 0; i < 9; i++)  q_packed[31 + i] = prods[GKYL_GR_SP_DBETA + i];
  for (int i = 0; i < 27; i++) q_packed[40 + i] = prods[GKYL_GR_SP_DGIJ + i];
  q_packed[67] = 0.0;
  q_packed[68] = x; q_packed[69] = y; q_packed[70] = z;

  // Hydro components, computed from the recovered primitive form.
  double spatial_det = prods[GKYL_GR_SP_SPATIAL_DET];
  double vel[3] = { u, v, w };
  double v_sq = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_sq += prods[GKYL_GR_SP_GIJ + 3*i + j] * vel[i] * vel[j];
  double W = 1.0 / sqrt(1.0 - v_sq);
  double h = 1.0 + (p / rho) * (gas_gamma / (gas_gamma - 1.0));

  q_packed[0] = sqrt(spatial_det) * rho * W;
  q_packed[1] = sqrt(spatial_det) * rho * h * (W*W) * u;
  q_packed[2] = sqrt(spatial_det) * rho * h * (W*W) * v;
  q_packed[3] = sqrt(spatial_det) * rho * h * (W*W) * w;
  q_packed[4] = sqrt(spatial_det) * ((rho * h * (W*W)) - p - (rho * W));
}

// (1) derive_products equivalence test.
static void
run_derive_products_check(struct gkyl_gr_spacetime *spacetime)
{
  // Tiny 4x4 grid on [-0.2, 0.2]^2.
  struct gkyl_rect_grid grid;
  double lower[2] = { -0.2, -0.2 };
  double upper[2] = {  0.2,  0.2 };
  int cells[2]    = { 4, 4 };
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int rlower[2] = { 1, 1 };
  int rupper[2] = { 4, 4 };
  struct gkyl_range update_range;
  gkyl_range_init(&update_range, 2, rlower, rupper);

  struct gkyl_array *prods = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, update_range.volume);

  struct gkyl_moment_spacetime_coupling_inp inp = {
    .grid               = &grid,
    .nfluids            = 1,
    .fluid_param        = {{ .type = GKYL_EQN_GR_EULER_MOD, .gas_gamma = 5.0/3.0 }},
    .is_static          = true,
    .has_tetrad         = false,
    .analytic_spacetime = spacetime,
    .spacetime_gauge    = GKYL_STATIC_GAUGE,
    .reinit_freq        = 0,
    .einstein_eqn       = NULL,
  };
  gkyl_moment_spacetime_coupling *st =
    gkyl_moment_spacetime_coupling_new(inp);

  gkyl_moment_spacetime_coupling_derive_products(
    st, 0.0, &update_range, NULL, prods);

  // For each cell verify the row matches a direct fill call at the same
  // cell-center coordinates.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&update_range, iter.idx);
    const double *got = gkyl_array_cfetch(prods, cidx);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    double expected[GKYL_GR_SP_NCOMP_BASE];
    gkyl_moment_spacetime_coupling_fill_products_analytic(spacetime,
      0.0, xc[0], xc[1], 0.0, expected);

    for (int k = 0; k < GKYL_GR_SP_NCOMP_BASE; k++)
      TEST_CHECK( gkyl_compare(got[k], expected[k], 1e-13) );
  }

  gkyl_array_release(prods);
  gkyl_moment_spacetime_coupling_release(st);
}

static void
test_derive_products_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_derive_products_check(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_derive_products_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 1.0, 1.0, 0.0);
  run_derive_products_check(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_derive_products_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 1.0, 1.0, 0.0);
  run_derive_products_check(spacetime);
  gkyl_gr_spacetime_release(spacetime);
}

// (2) Source-step equivalence: per-cell euler step matches packed on the
// hydro slice.
static void
run_source_euler_equivalence(struct gkyl_gr_spacetime *spacetime,
  int half_extent, double step)
{
  double gas_gamma = 5.0 / 3.0;

  // Build a minimal gkyl_moment_em_coupling to satisfy the packed source-
  // step signature. The packed function only reads gas_gamma from this
  // struct (via the caller passing it through explicitly here), so the
  // remaining fields are zero/false and that's fine.
  struct gkyl_rect_grid grid;
  double lower[2] = { 0.0, 0.0 };
  double upper[2] = { 1.0, 1.0 };
  int cells[2]    = { 4, 4 };
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_moment_em_coupling_inp em_inp = { 0 };
  em_inp.grid               = &grid;
  em_inp.nfluids            = 1;
  em_inp.param[0].type      = GKYL_EQN_GR_EULER;
  em_inp.param[0].mass      = 1.0;
  em_inp.has_gr_euler_sources = true;
  em_inp.gr_euler_gas_gamma   = gas_gamma;
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(em_inp);

  for (int x_ind = -half_extent; x_ind < half_extent + 1; x_ind++) {
    for (int y_ind = -half_extent; y_ind < half_extent + 1; y_ind++) {
      double x = step * x_ind;
      double y = step * y_ind;

      double q_packed[71];
      double prods[GKYL_GR_SP_NCOMP_BASE];

      // Use a fixed but non-trivial hydro state. Same numbers as the
      // existing wv_gr_euler unit test so we know they are inside the
      // physical regime for all three spacetimes.
      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      build_packed_and_products(spacetime, x, y, 0.0,
        rho, u, v, w, p, gas_gamma, q_packed, prods);

      // Skip excised cells: both paths short-circuit, so any difference is
      // harmless and not what we are testing here.
      if (q_packed[27] < 0.0) continue;

      double q_mod[5] = { q_packed[0], q_packed[1], q_packed[2], q_packed[3], q_packed[4] };

      const double dt = 1.0e-3;

      double f_packed_new[71];
      explicit_gr_euler_source_update_euler(mom_em, gas_gamma, 0.0, dt, q_packed, f_packed_new);

      double f_mod_new[5];
      gkyl_moment_spacetime_coupling_gr_euler_mod_source_euler(
        gas_gamma, 0.0, dt, prods, q_mod, f_mod_new);

      for (int i = 0; i < 5; i++)
        TEST_CHECK( gkyl_compare(f_mod_new[i], f_packed_new[i], 1e-12) );
    }
  }

  gkyl_moment_em_coupling_release(mom_em);
}

static void
test_source_euler_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_source_euler_equivalence(spacetime, 5, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_source_euler_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 1.0, 1.0, 0.0);
  run_source_euler_equivalence(spacetime, 5, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_source_euler_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 1.0, 1.0, 0.0);
  run_source_euler_equivalence(spacetime, 5, 0.1);
  gkyl_gr_spacetime_release(spacetime);
}

// (3) SSP-RK3 wrapper equivalence: full advance matches packed
// explicit_gr_euler_source_update on the hydro slice.
static void
run_explicit_advance_equivalence(struct gkyl_gr_spacetime *spacetime,
  double dt)
{
  double gas_gamma = 5.0 / 3.0;

  // Single-cell domain for simplicity; the wrapper just iterates cells.
  struct gkyl_rect_grid grid;
  double lower[1] = { 0.5 };
  double upper[1] = { 1.5 };
  int cells[1]    = { 1 };
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int rlower[1] = { 1 };
  int rupper[1] = { 1 };
  struct gkyl_range update_range;
  gkyl_range_init(&update_range, 1, rlower, rupper);

  // Cell center is x = 1.0 (consistent with the single-cell grid).
  double x = 1.0, y = 0.0, z = 0.0;

  double q_packed[71];
  double prods_row[GKYL_GR_SP_NCOMP_BASE];
  build_packed_and_products(spacetime, x, y, z,
    1.0, 0.1, 0.2, 0.3, 1.5, gas_gamma, q_packed, prods_row);

  if (q_packed[27] < 0.0) return;  // Skip if excised; nothing to compare.

  // Packed path: drive explicit_gr_euler_source_update on a 71-comp buffer.
  struct gkyl_moment_em_coupling_inp em_inp = { 0 };
  em_inp.grid                 = &grid;
  em_inp.nfluids              = 1;
  em_inp.param[0].type        = GKYL_EQN_GR_EULER;
  em_inp.param[0].mass        = 1.0;
  em_inp.has_gr_euler_sources = true;
  em_inp.gr_euler_gas_gamma   = gas_gamma;
  gkyl_moment_em_coupling *mom_em = gkyl_moment_em_coupling_new(em_inp);

  double f_packed[71];
  for (int i = 0; i < 71; i++) f_packed[i] = q_packed[i];

  double *fluid_packed_ptrs[GKYL_MAX_SPECIES] = { f_packed };
  explicit_gr_euler_source_update(mom_em, 0.0, dt, fluid_packed_ptrs);

  // Mod path: build a 1-cell gkyl_array, fill with the 5-comp hydro state,
  // and drive the explicit_advance.
  struct gkyl_array *fluid_mod = gkyl_array_new(GKYL_DOUBLE, 5, update_range.volume);
  struct gkyl_array *prods_arr = gkyl_array_new(GKYL_DOUBLE,
    GKYL_GR_SP_NCOMP_BASE, update_range.volume);

  {
    double *f = gkyl_array_fetch(fluid_mod, 0);
    for (int i = 0; i < 5; i++) f[i] = q_packed[i];
    double *p = gkyl_array_fetch(prods_arr, 0);
    for (int i = 0; i < GKYL_GR_SP_NCOMP_BASE; i++) p[i] = prods_row[i];
  }

  struct gkyl_moment_spacetime_coupling_inp st_inp = {
    .grid               = &grid,
    .nfluids            = 1,
    .fluid_param        = {{ .type = GKYL_EQN_GR_EULER_MOD, .gas_gamma = gas_gamma }},
    .is_static          = true,
    .has_tetrad         = false,
    .analytic_spacetime = spacetime,
    .spacetime_gauge    = GKYL_STATIC_GAUGE,
    .reinit_freq        = 0,
    .einstein_eqn       = NULL,
  };
  gkyl_moment_spacetime_coupling *st =
    gkyl_moment_spacetime_coupling_new(st_inp);

  struct gkyl_array *fluid_arrs[GKYL_MAX_SPECIES] = { fluid_mod };
  gkyl_moment_spacetime_coupling_explicit_advance(
    st, 0.0, dt, &update_range, fluid_arrs, prods_arr);

  const double *f_mod_final = gkyl_array_cfetch(fluid_mod, 0);
  for (int i = 0; i < 5; i++)
    TEST_CHECK( gkyl_compare(f_mod_final[i], f_packed[i], 1e-12) );

  gkyl_array_release(prods_arr);
  gkyl_array_release(fluid_mod);
  gkyl_moment_spacetime_coupling_release(st);
  gkyl_moment_em_coupling_release(mom_em);
}

static void
test_explicit_advance_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_explicit_advance_equivalence(spacetime, 1.0e-3);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_explicit_advance_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 1.0, 1.0, 0.0);
  run_explicit_advance_equivalence(spacetime, 1.0e-3);
  gkyl_gr_spacetime_release(spacetime);
}

static void
test_explicit_advance_kerr()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 1.0, 1.0, 0.0);
  run_explicit_advance_equivalence(spacetime, 1.0e-3);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "derive_products_minkowski",      test_derive_products_minkowski },
  { "derive_products_schwarzschild",  test_derive_products_schwarzschild },
  { "derive_products_kerr",           test_derive_products_kerr },
  { "source_euler_minkowski",         test_source_euler_minkowski },
  { "source_euler_schwarzschild",     test_source_euler_schwarzschild },
  { "source_euler_kerr",              test_source_euler_kerr },
  { "explicit_advance_minkowski",     test_explicit_advance_minkowski },
  { "explicit_advance_schwarzschild", test_explicit_advance_schwarzschild },
  { "explicit_advance_kerr",          test_explicit_advance_kerr },
  { NULL, NULL },
};
