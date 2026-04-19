#include <acutest.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_velocity_map.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct loss_cone_mask_test_ctx {
  int cdim;
  double eV;
  double R_m;
  double B_m;
  double z_m;
  double mass, charge;
  double n0, T0, B0;
  double phi_fac;
  double z_max, vpar_max, mu_max;
  int Nz, Nvpar, Nmu;
};

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
mapc2p_3x(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0];
  xp[1] = xc[1];
  xp[2] = xc[2];
}

void
bfield_func_3x(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[2];
  struct loss_cone_mask_test_ctx *params = ctx;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = params->B_m * (1.0 - ((params->R_m - 1.0) / params->R_m) * pow(cos(z), 2.0));
}

void
phi_func_1x_zero(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
phi_func_1x_nonzero(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0];
  struct loss_cone_mask_test_ctx *params = ctx;

  // fout[0] = params->phi_fac * params->T0 / params->eV * (1.0 - pow(z / params->z_max, 2.0));
  fout[0] = params->phi_fac * params->T0 / params->eV * (1.0 - pow(cos(4*z), 2.0));
}

static inline void
corner_coords_range(int ndim, int corner, double *eta)
{
  for (int d = 0; d < ndim; ++d) {
    int bit = ndim - 1 - d;
    eta[d] = ((corner >> bit) & 1) ? 1.0 : -1.0;
  }
}

static struct gkyl_array*
mk_basis_at_corners(const struct gkyl_basis *basis)
{
  int ndim = basis->ndim;
  int num_corners = 1 << ndim;
  struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, num_corners);

  double eta[GKYL_MAX_DIM] = { 0.0 };
  for (int corner = 0; corner < num_corners; ++corner) {
    corner_coords_range(ndim, corner, eta);
    basis->eval(eta, gkyl_array_fetch(out, corner));
  }

  return out;
}

static inline double
field_corner_val(const struct gkyl_array *arr, const struct gkyl_array *basis_at_corners,
  int num_basis, long linidx, int corner)
{
  const double *arr_d = gkyl_array_cfetch(arr, linidx);
  const double *basis_d = gkyl_array_cfetch(basis_at_corners, corner);

  double val = 0.0;
  for (int k = 0; k < num_basis; ++k) {
    val += arr_d[k] * basis_d[k];
  }
  return val;
}

static inline double
vel_corner_val(const struct gkyl_velocity_map *gvm, const int *vel_idx, int vd, int vel_corner,
  int vdim)
{
  long linidx_vel = gkyl_range_idx(&gvm->local_ext_vel, vel_idx);
  const double *vmap_d = gkyl_array_cfetch(gvm->vmap, linidx_vel);

  double vel_eta[GKYL_MAX_DIM] = { 0.0 };
  corner_coords_range(vdim, vel_corner, vel_eta);
  double xcomp[1] = { vel_eta[vd] };

  return gvm->vmap_basis->eval_expand(xcomp, vmap_d + vd * gvm->vmap_basis->num_basis);
}

static inline int
corner_z_endpoint_side_ref(int cdim, int corner, int zdim)
{
  double eta[GKYL_MAX_DIM] = { 0.0 };
  corner_coords_range(cdim, corner, eta);
  return eta[zdim] > 0.0 ? 1 : 0;
}

static inline int
corner_with_z_side_ref(int cdim, int anchor_corner, int zdim, int z_side)
{
  int num_conf_corners = 1 << cdim;

  double eta_anchor[GKYL_MAX_DIM] = { 0.0 };
  corner_coords_range(cdim, anchor_corner, eta_anchor);

  for (int cand = 0; cand < num_conf_corners; ++cand) {
    if (corner_z_endpoint_side_ref(cdim, cand, zdim) != z_side) {
      continue;
    }

    double eta_cand[GKYL_MAX_DIM] = { 0.0 };
    corner_coords_range(cdim, cand, eta_cand);

    bool same_transverse = true;
    for (int d = 0; d < cdim; ++d) {
      if (d == zdim) {
        continue;
      }
      if (eta_cand[d] != eta_anchor[d]) {
        same_transverse = false;
        break;
      }
    }

    if (same_transverse) {
      return cand;
    }
  }

  assert(false);
  return anchor_corner;
}

static void
escape_barriers_ref(const struct gkyl_array *phi, const struct gkyl_array *bmag,
  const struct gkyl_range *conf_range, const struct gkyl_array *basis_at_corners_conf,
  int cdim, int num_basis_conf, const int *base_idx, int target_z_cell,
  int anchor_corner, int anchor_z_side, double mu, double charge, double *barrier_left,
  double *barrier_right)
{
  int zdim = cdim - 1;
  int anchor_corner_node = corner_with_z_side_ref(cdim, anchor_corner, zdim, anchor_z_side);
  int z_upper_corner = corner_with_z_side_ref(cdim, anchor_corner, zdim, 1);
  int z_lower_corner = corner_with_z_side_ref(cdim, anchor_corner, zdim, 0);

  int scan_idx[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    scan_idx[d] = base_idx[d];
  }

  *barrier_left = -DBL_MAX;
  *barrier_right = -DBL_MAX;

  for (int iz = conf_range->lower[zdim]; iz <= conf_range->upper[zdim]; ++iz) {
    scan_idx[zdim] = iz;
    long linidx = gkyl_range_idx(conf_range, scan_idx);

    int left_corner = z_upper_corner;
    int right_corner = z_lower_corner;
    if (iz == target_z_cell) {
      left_corner = anchor_corner_node;
      right_corner = anchor_corner_node;
    }

    if (iz <= target_z_cell) {
      double phi_left = field_corner_val(phi, basis_at_corners_conf, num_basis_conf, linidx, left_corner);
      double bmag_left = field_corner_val(bmag, basis_at_corners_conf, num_basis_conf, linidx, left_corner);
      double u_left = mu * bmag_left + charge * phi_left;
      if (u_left > *barrier_left) {
        *barrier_left = u_left;
      }
    }
    if (iz >= target_z_cell) {
      double phi_right = field_corner_val(phi, basis_at_corners_conf, num_basis_conf, linidx, right_corner);
      double bmag_right = field_corner_val(bmag, basis_at_corners_conf, num_basis_conf, linidx, right_corner);
      double u_right = mu * bmag_right + charge * phi_right;
      if (u_right > *barrier_right) {
        *barrier_right = u_right;
      }
    }
  }
}

static void
build_reference_mask(const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_velocity_map *gvm,
  const struct gkyl_array *bmag, const struct gkyl_array *phi,
  double mass, double charge, struct gkyl_array *mask_ref)
{
  int cdim = conf_basis->ndim;
  int pdim = phase_range->ndim;
  int vdim = pdim - cdim;
  int num_basis_conf = conf_basis->num_basis;
  int num_phase_corners = 1 << pdim;

  struct gkyl_array *basis_at_corners_conf = mk_basis_at_corners(conf_basis);

  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    int conf_idx[GKYL_MAX_DIM] = { 0 };
    int vel_idx[GKYL_MAX_DIM] = { 0 };

    for (int d = 0; d < cdim; ++d) {
      conf_idx[d] = phase_iter.idx[d];
    }
    for (int d = 0; d < vdim; ++d) {
      vel_idx[d] = phase_iter.idx[cdim + d];
    }

    bool cell_trapped = true;
    for (int corner = 0; corner < num_phase_corners && cell_trapped; ++corner) {
      int conf_corner = corner / (1 << vdim);
      int vel_corner = corner % (1 << vdim);

      double vpar = vel_corner_val(gvm, vel_idx, 0, vel_corner, vdim);
      double mu = vel_corner_val(gvm, vel_idx, 1, vel_corner, vdim);

      long linidx_conf = gkyl_range_idx(conf_range, conf_idx);
      double bmag_curr = field_corner_val(bmag, basis_at_corners_conf, num_basis_conf,
        linidx_conf, conf_corner);
      double phi_curr = field_corner_val(phi, basis_at_corners_conf, num_basis_conf,
        linidx_conf, conf_corner);
      double h_curr = 0.5 * mass * vpar * vpar + mu * bmag_curr + charge * phi_curr;

      int zdim = cdim - 1;
      int target_z_cell = conf_idx[zdim];
      int anchor_z_side = corner_z_endpoint_side_ref(cdim, conf_corner, zdim);

      double barrier_left, barrier_right;
      escape_barriers_ref(phi, bmag, conf_range, basis_at_corners_conf, cdim,
        num_basis_conf, conf_idx, target_z_cell, conf_corner, anchor_z_side, mu, charge,
        &barrier_left, &barrier_right);

      cell_trapped = h_curr < GKYL_MIN2(barrier_left, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(phase_range, phase_iter.idx);
    double *m = gkyl_array_fetch(mask_ref, linidx_phase);
    m[0] = cell_trapped ? 1.0 : 0.0;
  }

  gkyl_array_release(basis_at_corners_conf);
}

static void
run_case_1x2v(int poly_order, bool use_gpu, bool use_nonzero_phi)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mass_proton = GKYL_PROTON_MASS;

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
    .phi_fac = 5.0,
    .z_max = M_PI-0.5,
    .Nz = 64,
    .Nvpar = 16,
    .Nmu = 16,
  };
  ctx.B0 = ctx.B_m / 2.0;
  ctx.vpar_max = 6.0 * sqrt(ctx.T0 / ctx.mass);
  ctx.mu_max = 0.5 * ctx.mass * pow(ctx.vpar_max, 2) / ctx.B0;

  double lower[] = { -ctx.z_max, -ctx.vpar_max, 0.0 };
  double upper[] = { ctx.z_max, ctx.vpar_max, ctx.mu_max };
  int cells[] = { ctx.Nz, ctx.Nvpar, ctx.Nmu };
  const int ndim = sizeof(cells) / sizeof(cells[0]);
  const int cdim = ctx.cdim;
  const int vdim = ndim - cdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; ++d) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; ++d) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = cells[cdim + d];
  }

  struct gkyl_rect_grid grid, grid_conf, grid_vel;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

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

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_ext_conf;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = ghost_conf[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_position_map *pmap = gkyl_position_map_null_new();

  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = { 0.0, 0.0 },
    .mapc2p = mapc2p_3x,
    .c2p_ctx = 0,
    .bfield_func = bfield_func_3x,
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

  struct gk_geometry *gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  if (use_gpu) {
    struct gk_geometry *gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);
  struct gkyl_velocity_map *gvm_ho = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, false);

  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);

  evalf_t phi_func = use_nonzero_phi ? phi_func_1x_nonzero : phi_func_1x_zero;
  gkyl_eval_on_nodes *evphi = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf, 1, phi_func, &ctx);
  gkyl_eval_on_nodes_advance(evphi, 0.0, &local_conf, phi_ho);
  gkyl_eval_on_nodes_release(evphi);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
    .conf_basis = &basis_conf,
    .conf_range = &local_conf,
    .vel_map = gvm,
    .mass = ctx.mass,
    .charge = ctx.charge,
  };
  struct gkyl_loss_cone_mask_gyrokinetic *proj_mask =
    gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

  struct gkyl_array *mask = mkarr(use_gpu, 1, local_ext.volume);
  struct gkyl_array *mask_ho = mkarr(false, 1, local_ext.volume);
  struct gkyl_array *mask_ref = mkarr(false, 1, local_ext.volume);

  gkyl_loss_cone_mask_gyrokinetic_advance(proj_mask, &local, &local_conf,
    gk_geom->geo_corn.bmag, phi, mask);
  gkyl_array_copy(mask_ho, mask);

  struct gkyl_array *bmag_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  gkyl_array_copy(bmag_ho, gk_geom->geo_corn.bmag);

  {
    char f_bmag[256], f_phi[256], f_mask[256];
    const char *suffix_dev = use_gpu ? "dev" : "ho";
    const char *suffix_phi = use_nonzero_phi ? "nonzero_phi" : "zero_phi";
    const struct gkyl_range arr_range_conf = local_conf;
    const struct gkyl_range arr_range = local;

    snprintf(f_bmag, sizeof(f_bmag), "ctest_loss_cone_mask_%s_%s_bmag.gkyl", suffix_dev,
      suffix_phi);
    snprintf(f_phi, sizeof(f_phi), "ctest_loss_cone_mask_%s_%s_phi.gkyl", suffix_dev,
      suffix_phi);
    snprintf(f_mask, sizeof(f_mask), "ctest_loss_cone_mask_%s_%s_mask.gkyl", suffix_dev,
      suffix_phi);

    gkyl_grid_sub_array_write(&grid_conf, &arr_range_conf, 0, bmag_ho, f_bmag);
    gkyl_grid_sub_array_write(&grid_conf, &arr_range_conf, 0, phi_ho, f_phi);
    gkyl_grid_sub_array_write(&grid, &arr_range, 0, mask_ho, f_mask);
  }

  gkyl_array_clear(mask_ref, 0.0);
  build_reference_mask(&local, &local_conf, &basis_conf, gvm_ho,
    bmag_ho, phi_ho, ctx.mass, ctx.charge, mask_ref);

  int mismatches = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *m = gkyl_array_cfetch(mask_ho, linidx);
    const double *r = gkyl_array_cfetch(mask_ref, linidx);

    bool ok = fabs(m[0] - r[0]) < 1e-12;
    TEST_CHECK(ok);
    if (!ok && mismatches < 8) {
      printf("mask mismatch idx=(%d,%d,%d): got=%g ref=%g\n",
        iter.idx[0], iter.idx[1], iter.idx[2], m[0], r[0]);
      mismatches++;
    }
  }

  gkyl_array_release(bmag_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(mask);
  gkyl_array_release(mask_ho);
  gkyl_array_release(mask_ref);

  gkyl_loss_cone_mask_gyrokinetic_release(proj_mask);
  gkyl_velocity_map_release(gvm);
  gkyl_velocity_map_release(gvm_ho);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_position_map_release(pmap);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(basis_on_dev_conf);
  }
#else
  (void) basis_on_dev;
  (void) basis_on_dev_conf;
#endif
}

void
test_1x2v_p1_gk_ho(void)
{
  run_case_1x2v(1, false, false);
}

void
test_1x2v_p1_nonzero_phi_gk_ho(void)
{
  run_case_1x2v(1, false, true);
}

#ifdef GKYL_HAVE_CUDA
void
test_1x2v_p1_gk_dev(void)
{
  run_case_1x2v(1, true, false);
}

void
test_1x2v_p1_nonzero_phi_gk_dev(void)
{
  run_case_1x2v(1, true, true);
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
