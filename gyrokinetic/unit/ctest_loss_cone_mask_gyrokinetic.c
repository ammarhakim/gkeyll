#include <acutest.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_const.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct loss_cone_mask_test_ctx {
  int cdim; // Configuration space dimensionality.
  double eV; // Elementary charge.
  double R_m; // Mirror ratio.
  double B_m; // Maximum magnetic field amplitude.
  double z_m; // Location of B_m.
  double mass, charge; // Species mass and charge.
  double n0, T0, B0; // Reference parameters.
  double phi_fac; // phi(z=0) = phi_fac*T0/e;
  double psi_max; // For 2x: upper limit of psi (radial coordinate).
  double z_max, vpar_max, mu_max; // Upper grid extents.
  int Npsi, Nz, Nvpar, Nmu; // Number of cells in each direction.
  enum gkyl_quad_type quad_type; // Type of quadrature/nodes.
  int num_quad; // Number of quadrature points to use in projection, 1 or p+1.
  bool cellwise_trap_loss; // Whether a whole cell is either trapped or lost.
  bool is_tandem; // Whether this is a tandem mirror configuration.
  double B_tandem; // Field at tandem mirror (for tandem case).
  double z_tandem; // z-coordinate of tandem mirror (for tandem case).
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                          : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
mapc2p_3x(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bfield_func_3x(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct loss_cone_mask_test_ctx *params = ctx;
  double R_m = params->R_m; // Mirror ratio.
  double B_m = params->B_m; // Maximum magnetic field amplitude.

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B_m * (1.0 - ((R_m - 1.0) / R_m) * pow(cos(z), 2.0));
// fout[0] = (B_m/R_m) * (1.0 + (R_m-1.0)*pow(sin(z), 2.0));
}

void
phi_func_1x(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0];

  struct loss_cone_mask_test_ctx *params = ctx;
  double phi_fac = params->phi_fac;
  double T0 = params->T0;
  double eV = params->eV;

  fout[0] = 0.0; // 0.5 * phi_fac*T0/eV * (1.0 + cos(z));
}

// Non-zero electrostatic potential: peaked at center, zero at wall.
void
phi_func_1x_nonzero(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0];

  struct loss_cone_mask_test_ctx *params = ctx;
  double phi_fac = params->phi_fac;
  double T0 = params->T0;
  double eV = params->eV;
  double z_max = params->z_max;

  // Parabolic potential profile: phi(z) = phi_fac*T0/eV * (1 - (z/z_max)^2)
  // This gives phi=phi_fac*T0/eV at z=0, and phi=0 at z=+/-z_max.
  fout[0] = phi_fac * T0 / eV * (1.0 - pow(z / z_max, 2.0));
}

// Reference mask for nonzero phi case.
void
mask_ref_1x2v_nonzero_phi(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0], vpar = xc[1], mu = xc[2];
  struct loss_cone_mask_test_ctx *params = ctx;

  double z_m = params->z_m;
  double mass = params->mass;
  double charge = params->charge;

  double phi, phi_m;
  phi_func_1x_nonzero(t, xc, &phi, ctx);
  phi_func_1x_nonzero(t, &z_m, &phi_m, ctx);

  double bfield[3], bmag;
  double zinfl[3] = { 0.0 }, z_minfl[3] = { 0.0 };
  zinfl[2] = z, z_minfl[2] = z_m;
  bfield_func_3x(t, zinfl, bfield, ctx);
  bmag = bfield[2];

  double bfield_m[3], bmag_m;
  bfield_func_3x(t, z_minfl, bfield_m, ctx);
  bmag_m = bfield_m[2];

  // mu_bound = (0.5*m*vpar^2+q*(phi-phi_m))/(B*(B_max/B-1))
  double mu_bound = (0.5 * mass * pow(vpar,
    2) + charge * (phi - phi_m)) / (bmag * (bmag_m / bmag - 1));
  if (mu_bound < mu && fabs(z) < z_m) {
    fout[0] = 1.0;
  }
  else {
    fout[0] = 0.0;
  }
}

void
mask_ref_1x2v(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0], vpar = xc[1], mu = xc[2];
  struct loss_cone_mask_test_ctx *params = ctx;

  double z_m = params->z_m;
  double mass = params->mass;
  double charge = params->charge;

  double phi, phi_m;
  phi_func_1x(t, xc, &phi, ctx);
  phi_func_1x(t, &z_m, &phi_m, ctx);

  double bfield[3], bmag;
  double zinfl[3] = { 0.0 }, z_minfl[3] = { 0.0 };
  zinfl[2] = z, z_minfl[2] = z_m;
  bfield_func_3x(t, zinfl, bfield, ctx);
  bmag = bfield[2];

  double bfield_m[3], bmag_m;
  bfield_func_3x(t, z_minfl, bfield_m, ctx);
  bmag_m = bfield_m[2];

  // mu_bound = (0.5*m*vpar^2+q*(phi-phi_m))/(B*(B_max/B-1))
  double mu_bound = (0.5 * mass * pow(vpar,
    2) + charge * (phi - phi_m)) / (bmag * (bmag_m / bmag - 1));
  if (mu_bound < mu && fabs(z) < z_m)
    fout[0] = 1.0;
  else
    fout[0] = 0;
}

void
test_1x2v_gk(int poly_order, bool use_gpu)
{

  double eV = GKYL_ELEMENTARY_CHARGE;
  double mass_proton = GKYL_PROTON_MASS;

  // Set reference parameters.
  struct loss_cone_mask_test_ctx ctx = {
    .cdim = 1,
    .eV = eV,
    .R_m = 8.0,
    .B_m = 4.0,
    .z_m = M_PI / 2.0,
    .mass = 2.014 * mass_proton,
    .charge = eV,
    .n0 = 1e18,
    .T0 = 100 * eV,
    .phi_fac = 3.0,
    .z_max = M_PI,
    .Nz = 8,
    .Nvpar = 8,
    .Nmu = 4,
    .quad_type = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 2,
    .cellwise_trap_loss = true,
  };
  ctx.B0 = ctx.B_m / 2.0;
  ctx.vpar_max = 6.0 * sqrt(ctx.T0 / ctx.mass);
  ctx.mu_max = 0.5 * ctx.mass * pow(ctx.vpar_max, 2) / ctx.B0;

  double mass = ctx.mass;
  double lower[] = { -ctx.z_max, -ctx.vpar_max, 0.0 },
    upper[] = { ctx.z_max, ctx.vpar_max, ctx.mu_max };
  int cells[] = { ctx.Nz, ctx.Nvpar, ctx.Nmu };
  const int ndim = sizeof(cells) / sizeof(cells[0]);
  const int cdim = ctx.cdim;
  const int vdim = ndim - ctx.cdim;

  // Grids.
  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; d++) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = cells[cdim + d];
  }
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis, basis_conf;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *basis_on_dev_conf;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    basis_on_dev_conf = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1)
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(basis_on_dev_conf, cdim, poly_order);
#endif
  }
  else {
    basis_on_dev = &basis;
    basis_on_dev_conf = &basis_conf;
  }

  // Ranges.
  int ghost_conf[] = { 1, 1, 1 }; // 3 elements because it's used by geo.
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d = 0; d < cdim; d++) {
    ghost[d] = ghost_conf[d];
  }
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = { 0.0, 0.0 },
    .mapc2p = mapc2p_3x, // mapping of computational to physical space
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x, // magnetic field magnitude
    .bfield_ctx = &ctx,
    .position_map = pmap,
    .grid = grid_conf,
    .local = local_conf,
    .local_ext = local_ext_conf,
    .global = local_conf,
    .global_ext = local_ext_conf,
    .basis = basis_conf,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(grid_conf, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, ghost_conf, &geometry_input.geo_local_ext,
    &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry *gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // Deflate geometry if necessary.
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // Use array_dg_find_peaks to find bmag_max along the z direction.
  // Search along the parallel (z) direction, which is the last configuration space dimension.
  int search_dir = cdim - 1;
  struct gkyl_array_dg_find_peaks_inp peak_inp = {
    .basis = &basis_conf,
    .grid = &grid_conf,
    .range = &local_conf,
    .range_ext = &local_ext_conf,
    .search_dir = search_dir,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_dg_find_peaks *bmag_peak_finder =
    gkyl_array_dg_find_peaks_new(&peak_inp, gk_geom->geo_int.bmag);
  gkyl_array_dg_find_peaks_advance(bmag_peak_finder, gk_geom->geo_int.bmag);

  // Get the LOCAL_MAX peak (bmag maximum along z direction).
  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(bmag_peak_finder);
  int bmag_max_peak_idx = num_peaks - 2; // Edge is num_peaks-1, so maximum is one less
  const struct gkyl_array *bmag_max = gkyl_array_dg_find_peaks_acquire_vals(bmag_peak_finder,
    bmag_max_peak_idx);
  const struct gkyl_array *bmag_max_z_coord =
    gkyl_array_dg_find_peaks_acquire_coords(bmag_peak_finder, bmag_max_peak_idx);
  const struct gkyl_array *bmag_wall = gkyl_array_dg_find_peaks_acquire_vals(bmag_peak_finder,
    num_peaks - 1);                                                                                          // First peak is wall
  const struct gkyl_array *bmag_wall_z_coord =
    gkyl_array_dg_find_peaks_acquire_coords(bmag_peak_finder, num_peaks - 1);
  const struct gkyl_basis *bmag_max_basis = gkyl_array_dg_find_peaks_get_basis(bmag_peak_finder);
  const struct gkyl_range *bmag_max_range = gkyl_array_dg_find_peaks_get_range(bmag_peak_finder);
  const struct gkyl_range *bmag_max_range_ext =
    gkyl_array_dg_find_peaks_get_range_ext(bmag_peak_finder);

  // Allocate arrays for phi evaluated at all peak locations.
  struct gkyl_array **phi_at_peaks = gkyl_malloc(num_peaks * sizeof(struct gkyl_array *));
  for (int p = 0; p < num_peaks; p++) {
    phi_at_peaks[p] = mkarr(use_gpu, bmag_max_basis->num_basis, bmag_max_range_ext->volume);
  }

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry *gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);

  // Project the electostatic potential.
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size)
                               : gkyl_array_acquire(phi);

  gkyl_eval_on_nodes *evphi = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf, 1, phi_func_1x, &ctx);
  gkyl_eval_on_nodes_advance(evphi, 0.0, &local_conf, phi_ho);
  gkyl_eval_on_nodes_release(evphi);
  gkyl_array_copy(phi, phi_ho);

  // Project phi onto peak locations to get phi_m at the mirror throat.
  gkyl_array_dg_find_peaks_project_on_peaks(bmag_peak_finder, phi, phi_at_peaks);

  // Get phi at the mirror throat (bmag_max peak location).
  const struct gkyl_array *phi_m = phi_at_peaks[bmag_max_peak_idx];

  // Basis used to project the mask.
  struct gkyl_basis basis_mask;
  if (ctx.num_quad == 1 || ctx.cellwise_trap_loss) {
    gkyl_cart_modal_serendip(&basis_mask, ndim, 0);
  }
  else {
    if (poly_order == 1) {
      gkyl_cart_modal_gkhybrid(&basis_mask, cdim, vdim);
    }
    else {
      gkyl_cart_modal_serendip(&basis_mask, ndim, poly_order);
    }
  }

  // Create mask array.
  struct gkyl_array *mask = mkarr(use_gpu, basis_mask.num_basis, local_ext.volume);
  struct gkyl_array *mask_ho = use_gpu? mkarr(false, mask->ncomp, mask->size)
                                : gkyl_array_acquire(mask);

  // Project the loss cone mask.
  // Use bmag_max and bmag_max_z_coord arrays from find_peaks.
  struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &basis_conf,
    .phase_basis = &basis,
    .conf_range = &local_conf,
    .conf_range_ext = &local_ext_conf,
    .vel_range = &local_vel,
    .vel_map = gvm,
    .bmag = gk_geom->geo_int.bmag,
    .bmag_max_z_coord = bmag_max_z_coord,
    .bmag_max = bmag_max,
    .bmag_wall = bmag_wall,
    .bmag_wall_z_coord = bmag_wall_z_coord,
    .bmag_max_basis = bmag_max_basis,
    .bmag_max_range = bmag_max_range,
    .mass = ctx.mass,
    .charge = ctx.charge,
    .qtype = ctx.quad_type,
    .num_quad = ctx.num_quad,
    .cellwise_trap_loss = ctx.cellwise_trap_loss,
    .use_gpu = use_gpu,
  };
  struct gkyl_loss_cone_mask_gyrokinetic *proj_mask =
    gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

  gkyl_loss_cone_mask_gyrokinetic_advance(proj_mask, &local, &local_conf, phi, phi_m, phi_m, mask);

  gkyl_array_copy(mask_ho, mask);

  // Project expected mask.
  struct gkyl_array *mask_ref_ho = mkarr(false, basis_mask.num_basis, local_ext.volume);
  gkyl_proj_on_basis *evmask_ref = gkyl_proj_on_basis_new(&grid, &basis_mask,
    basis_mask.poly_order + 1, 1, mask_ref_1x2v, &ctx);
  gkyl_proj_on_basis_advance(evmask_ref, 0.0, &local, mask_ref_ho);
  gkyl_proj_on_basis_release(evmask_ref);

//// values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
// double p1_vals[] = {
// 7.2307139183122714e-03, 0.0000000000000000e+00, 1.9198293226362615e-04, -7.7970439910196674e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
// -2.0701958137127286e-05, 0.0000000000000000e+00, -1.4953406100022537e-04, 0.0000000000000000e+00, 1.6124599381836546e-05, 0.0000000000000000e+00,
// -8.2719200283232917e-19, 0.0000000000000000e+00, -3.4806248503322844e-20, 0.0000000000000000e+00, };
// double p2_vals[] = {
// 7.2307468609012666e-03, 0.0000000000000000e+00, 1.9198380692343289e-04, -7.8092230706225602e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
// -2.0734294852987710e-05, 3.6591823321385775e-18, -1.4953474226616330e-04, 3.7739922227981074e-05, 0.0000000000000000e+00, 7.0473141211557788e-19,
// 0.0000000000000000e+00, -4.8789097761847700e-19, 1.6149786206441256e-05, 0.0000000000000000e+00, 1.0020339643610290e-06, 5.4210108624275222e-20,
// 0.0000000000000000e+00, 0.0000000000000000e+00 };
//
// const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
// if (poly_order == 1) {
// for (int i=0; i<basis.num_basis; ++i) {
// TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-2) );
// }
// }
//
// if (poly_order == 2) {
// for (int i=0; i<basis.num_basis; ++i)
// TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-2) );
// }

  // Write mask to file.
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_loss_cone_mask_gyrokinetic_1x2v_p%d_dev.gkyl", poly_order);
  }
  else {
    sprintf(fname, "ctest_loss_cone_mask_gyrokinetic_1x2v_p%d_ho.gkyl", poly_order);
  }
  gkyl_grid_sub_array_write(&grid, &local, 0, mask_ho, fname);

  sprintf(fname, "ctest_loss_cone_mask_gyrokinetic_1x2v_p%d_ref.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, mask_ref_ho, fname);

  // Free phi_m and phi_at_peaks arrays.
  for (int p = 0; p < num_peaks; p++) {
    gkyl_array_release(phi_at_peaks[p]);
  }
  gkyl_free(phi_at_peaks);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(mask);
  gkyl_array_release(mask_ho);
  gkyl_array_release(mask_ref_ho);
  gkyl_loss_cone_mask_gyrokinetic_release(proj_mask);
  gkyl_velocity_map_release(gvm);
  // Release acquired peak arrays.
  gkyl_array_release(bmag_max);
  gkyl_array_release(bmag_max_z_coord);
  gkyl_array_release(bmag_wall);
  gkyl_array_release(bmag_wall_z_coord);
  gkyl_array_dg_find_peaks_release(bmag_peak_finder);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(basis_on_dev_conf);
  }
#endif
}

// Test with non-zero electrostatic potential.
// This tests that the trapped-passing boundary correctly accounts for
// the q*(phi-phi_m) term in the mu_bound calculation.
void
test_1x2v_nonzero_phi_gk(int poly_order, bool use_gpu)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mass_proton = GKYL_PROTON_MASS;

  // Set reference parameters.
  struct loss_cone_mask_test_ctx ctx = {
    .cdim = 1,
    .eV = eV,
    .R_m = 8.0,
    .B_m = 4.0,
    .z_m = M_PI / 2.0,
    .mass = 2.014 * mass_proton,
    .charge = eV,  // Positive ions.
    .n0 = 1e18,
    .T0 = 100 * eV,
    .phi_fac = 3.0,  // phi(z=0) = 3*T0/e = 300 V.
    .z_max = M_PI,
    .Nz = 8,
    .Nvpar = 8,
    .Nmu = 4,
    .quad_type = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 2,
    .cellwise_trap_loss = true,
  };
  ctx.B0 = ctx.B_m / 2.0;
  ctx.vpar_max = 6.0 * sqrt(ctx.T0 / ctx.mass);
  ctx.mu_max = 0.5 * ctx.mass * pow(ctx.vpar_max, 2) / ctx.B0;

  double lower[] = { -ctx.z_max, -ctx.vpar_max, 0.0 };
  double upper[] = { ctx.z_max, ctx.vpar_max, ctx.mu_max };
  int cells[] = { ctx.Nz, ctx.Nvpar, ctx.Nmu };
  const int ndim = sizeof(cells) / sizeof(cells[0]);
  const int cdim = ctx.cdim;
  const int vdim = ndim - ctx.cdim;

  // Grids.
  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; d++) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = cells[cdim + d];
  }
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis, basis_conf;
  if (poly_order == 1) {
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  else {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  }
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *basis_on_dev_conf;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    basis_on_dev_conf = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) {
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    }
    else {
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    }
    gkyl_cart_modal_serendip_cu_dev(basis_on_dev_conf, cdim, poly_order);
#endif
  }
  else {
    basis_on_dev = &basis;
    basis_on_dev_conf = &basis_conf;
  }

  // Ranges.
  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_ext_conf;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d = 0; d < cdim; d++) {
    ghost[d] = ghost_conf[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  // Initialize geometry.
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = { 0.0, 0.0 },
    .mapc2p = mapc2p_3x,
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x,
    .bfield_ctx = &ctx,
    .grid = grid_conf,
    .local = local_conf,
    .local_ext = local_ext_conf,
    .global = local_conf,
    .global_ext = local_ext_conf,
    .basis = basis_conf,
    .position_map = pmap,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(grid_conf, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, ghost_conf, &geometry_input.geo_local_ext,
    &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry *gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // Use array_dg_find_peaks to find bmag_max.
  int search_dir = cdim - 1;
  struct gkyl_array_dg_find_peaks_inp peak_inp = {
    .basis = &basis_conf,
    .grid = &grid_conf,
    .range = &local_conf,
    .range_ext = &local_ext_conf,
    .search_dir = search_dir,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_dg_find_peaks *bmag_peak_finder =
    gkyl_array_dg_find_peaks_new(&peak_inp, gk_geom->geo_int.bmag);
  gkyl_array_dg_find_peaks_advance(bmag_peak_finder, gk_geom->geo_int.bmag);

  int num_peaks = gkyl_array_dg_find_peaks_num_peaks(bmag_peak_finder);
  int bmag_max_peak_idx = num_peaks - 2;
  const struct gkyl_array *bmag_max = gkyl_array_dg_find_peaks_acquire_vals(bmag_peak_finder,
    bmag_max_peak_idx);
  const struct gkyl_array *bmag_max_z_coord =
    gkyl_array_dg_find_peaks_acquire_coords(bmag_peak_finder, bmag_max_peak_idx);
  const struct gkyl_array *bmag_wall = gkyl_array_dg_find_peaks_acquire_vals(bmag_peak_finder,
    num_peaks - 1);
  const struct gkyl_array *bmag_wall_z_coord =
    gkyl_array_dg_find_peaks_acquire_coords(bmag_peak_finder, num_peaks - 1);
  const struct gkyl_basis *bmag_max_basis = gkyl_array_dg_find_peaks_get_basis(bmag_peak_finder);
  const struct gkyl_range *bmag_max_range = gkyl_array_dg_find_peaks_get_range(bmag_peak_finder);
  const struct gkyl_range *bmag_max_range_ext =
    gkyl_array_dg_find_peaks_get_range_ext(bmag_peak_finder);

  // Allocate arrays for phi evaluated at peak locations.
  struct gkyl_array **phi_at_peaks = gkyl_malloc(num_peaks * sizeof(struct gkyl_array *));
  for (int p = 0; p < num_peaks; p++) {
    phi_at_peaks[p] = mkarr(use_gpu, bmag_max_basis->num_basis, bmag_max_range_ext->volume);
  }

  if (use_gpu) {
    struct gk_geometry *gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);

  // Project the electrostatic potential with NON-ZERO phi.
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_ho = use_gpu ? mkarr(false, phi->ncomp, phi->size)
                                      : gkyl_array_acquire(phi);

  gkyl_eval_on_nodes *evphi = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf, 1,
    phi_func_1x_nonzero, &ctx);
  gkyl_eval_on_nodes_advance(evphi, 0.0, &local_conf, phi_ho);
  gkyl_eval_on_nodes_release(evphi);
  gkyl_array_copy(phi, phi_ho);

  // Project phi onto peak locations to get phi_m.
  gkyl_array_dg_find_peaks_project_on_peaks(bmag_peak_finder, phi, phi_at_peaks);
  const struct gkyl_array *phi_m = phi_at_peaks[bmag_max_peak_idx];

  // Basis used to project the mask.
  struct gkyl_basis basis_mask;
  if (ctx.num_quad == 1 || ctx.cellwise_trap_loss) {
    gkyl_cart_modal_serendip(&basis_mask, ndim, 0);
  }
  else {
    if (poly_order == 1) {
      gkyl_cart_modal_gkhybrid(&basis_mask, cdim, vdim);
    }
    else {
      gkyl_cart_modal_serendip(&basis_mask, ndim, poly_order);
    }
  }

  // Create mask array.
  struct gkyl_array *mask = mkarr(use_gpu, basis_mask.num_basis, local_ext.volume);
  struct gkyl_array *mask_ho = use_gpu ? mkarr(false, mask->ncomp, mask->size)
                                       : gkyl_array_acquire(mask);

  // Project the loss cone mask.
  struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
    .phase_grid = &grid,
    .conf_basis = &basis_conf,
    .phase_basis = &basis,
    .conf_range = &local_conf,
    .conf_range_ext = &local_ext_conf,
    .vel_range = &local_vel,
    .vel_map = gvm,
    .bmag = gk_geom->geo_int.bmag,
    .bmag_max_z_coord = bmag_max_z_coord,
    .bmag_max = bmag_max,
    .bmag_wall = bmag_wall,
    .bmag_wall_z_coord = bmag_wall_z_coord,
    .bmag_max_basis = bmag_max_basis,
    .bmag_max_range = bmag_max_range,
    .mass = ctx.mass,
    .charge = ctx.charge,
    .qtype = ctx.quad_type,
    .num_quad = ctx.num_quad,
    .cellwise_trap_loss = ctx.cellwise_trap_loss,
    .use_gpu = use_gpu,
  };
  struct gkyl_loss_cone_mask_gyrokinetic *proj_mask =
    gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

  gkyl_loss_cone_mask_gyrokinetic_advance(proj_mask, &local, &local_conf, phi, phi_m, phi_m, mask);

  gkyl_array_copy(mask_ho, mask);

  // Verify physical properties of the mask:
  // 1. At the center (z≈0), high-mu particles should be trapped (mask=1)
  // 2. At the wall (|z| ≈ z_max), particles should not be in the trapped region
  // 3. Low-mu particles near center should be passing (mask=0)

  // Check specific cells to verify correct behavior.
  // Cell indices: [iz, ivpar, imu] where each starts at 1 in local range.
  // Grid: z in [-pi, pi], vpar in [-vpar_max, vpar_max], mu in [0, mu_max]
  // Central z cells are around iz=4,5 (8 cells, symmetric)
  // High mu cells are imu=3,4 (4 cells)
  // Low mu cells are imu=1

  int num_trapped_high_mu_center = 0;
  int num_passing_low_mu_center = 0;
  int total_high_mu_center = 0;
  int total_low_mu_center = 0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    int iz = iter.idx[0];
    int imu = iter.idx[2];

    // Determine if we're at center (iz = 4 or 5 for 8 cells in [-pi, pi])
    // and if we're at high mu (imu = 3 or 4) or low mu (imu = 1)
    bool is_center = (iz == 4 || iz == 5);
    bool is_high_mu = (imu == 3 || imu == 4);
    bool is_low_mu = (imu == 1);

    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *mask_val = gkyl_array_cfetch(mask_ho, linidx);

    if (is_center && is_high_mu) {
      total_high_mu_center++;
      if (mask_val[0] > 0.5) {
        num_trapped_high_mu_center++;
      }
    }
    if (is_center && is_low_mu) {
      total_low_mu_center++;
      if (mask_val[0] < 0.5) {
        num_passing_low_mu_center++;
      }
    }
  }

  // High mu particles at center should mostly be trapped.
  double trapped_frac = (double)num_trapped_high_mu_center / (double)total_high_mu_center;
  // printf("Trapped fraction for high-mu center particles: %g (%d / %d)\n",
  //   trapped_frac, num_trapped_high_mu_center, total_high_mu_center);
  TEST_CHECK(trapped_frac >= 0.5);
  if (trapped_frac < 0.5) {
    printf("High-mu center trapped fraction: %g (%d / %d)\n",
      trapped_frac, num_trapped_high_mu_center, total_high_mu_center);
  }

  // Low mu particles at center should mostly be passing.
  double passing_frac = (double)num_passing_low_mu_center / (double)total_low_mu_center;
  // printf("Passing fraction for low-mu center particles: %g (%d / %d)\n",
  //   passing_frac, num_passing_low_mu_center, total_low_mu_center);
  TEST_CHECK(passing_frac >= 0.5);
  if (passing_frac < 0.5) {
    printf("Low-mu center passing fraction: %g (%d / %d)\n",
      passing_frac, num_passing_low_mu_center, total_low_mu_center);
  }

  // Write output for debugging.
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_loss_cone_mask_gyrokinetic_1x2v_nonzero_phi_p%d_dev.gkyl", poly_order);
  }
  else {
    sprintf(fname, "ctest_loss_cone_mask_gyrokinetic_1x2v_nonzero_phi_p%d_ho.gkyl", poly_order);
  }
  gkyl_grid_sub_array_write(&grid, &local, 0, mask_ho, fname);

  // Cleanup.
  for (int p = 0; p < num_peaks; p++) {
    gkyl_array_release(phi_at_peaks[p]);
  }
  gkyl_free(phi_at_peaks);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(mask);
  gkyl_array_release(mask_ho);
  gkyl_loss_cone_mask_gyrokinetic_release(proj_mask);
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(bmag_max);
  gkyl_array_release(bmag_max_z_coord);
  gkyl_array_release(bmag_wall);
  gkyl_array_release(bmag_wall_z_coord);
  gkyl_array_dg_find_peaks_release(bmag_peak_finder);
  gkyl_position_map_release(pmap);
  gkyl_gk_geometry_release(gk_geom);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(basis_on_dev_conf);
  }
#endif
}

void test_1x2v_p1_gk_ho()
{
  test_1x2v_gk(1, false);
}

void test_1x2v_p1_nonzero_phi_gk_ho()
{
  test_1x2v_nonzero_phi_gk(1, false);
}

#ifdef GKYL_HAVE_CUDA
void test_1x2v_p1_gk_dev()
{
  test_1x2v_gk(1, true);
}

void test_1x2v_p1_nonzero_phi_gk_dev()
{
  test_1x2v_nonzero_phi_gk(1, true);
}

#endif

TEST_LIST = {
  { "test_1x2v_p1_gk_ho", test_1x2v_p1_gk_ho },
  { "test_1x2v_p1_nonzero_phi_gk_ho", test_1x2v_p1_nonzero_phi_gk_ho },

#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gk_dev", test_1x2v_p1_gk_dev },
  { "test_1x2v_p1_nonzero_phi_gk_dev", test_1x2v_p1_nonzero_phi_gk_dev },
#endif
  { NULL, NULL },
};
