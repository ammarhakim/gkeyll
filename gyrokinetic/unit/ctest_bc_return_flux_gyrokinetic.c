// Test updater which applies excision BCs.
// 
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_bc_return_flux_gyrokinetic.h>
#include <assert.h>

// Allocate array (filled with zeros).
static struct gkyl_array *mkarr(bool use_gpu, long nc, long size) {
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost) {
  int ndim = parent->ndim;
  for (int d = 0; d < ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d], d,
      GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d], d,
      GKYL_UPPER_EDGE, parent, ghost);
  }
}

void evalFunc_2x2v(double t, const double *xn, double *restrict fout,
  void *ctx) {
  double x = xn[0], y = xn[1];
  double vx = xn[2], vy = xn[3];
  fout[0] = y * (vx - 1) * (vy - 2);
}

void evalFunc_3x2v(double t, const double *xn, double *restrict fout,
  void *ctx) {
  double x = xn[0], y = xn[1], z = xn[2];
  double vx = xn[3], vy = xn[4];
  fout[0] = y * (z + 1) * (vx - 1) * (vy - 2);
}

void test_bc(int cdim, int vdim, int poly_order, bool use_gpu) {
  int bc_dir = 0; // Direction perpendicular to boundary.
  enum gkyl_edge_loc bc_edge = GKYL_LOWER_EDGE; // Boundary edge.
  int disp_dir = 1; // Direction in which to displace the boundary flux.

  const int ndim = cdim + vdim;
  double lower[ndim], upper[ndim];
  int cells[ndim];
  for (int i = 0; i < ndim; i++) {
    lower[i] = -2.0;
    upper[i] =  2.0;
    cells[i] = i < cdim ? 8 : 2;
  }
  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int i = 0; i < cdim; i++) {
    lower_conf[i] = lower[i];
    upper_conf[i] = upper[i];
    cells_conf[i] = cells[i];
  }

  // Grids
  struct gkyl_rect_grid grid, grid_conf;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions
  struct gkyl_basis basis, basis_conf;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost[ndim];
  for (int d = 0; d < ndim; d++) ghost[d] = 0;
  for (int d = 0; d < cdim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create distribution function array
  struct gkyl_array *distf, *distf_ho;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);

  // Ghost and skin ranges at desired boundary.
  struct gkyl_range *ghost_r = bc_edge==GKYL_LOWER_EDGE? &skin_ghost.lower_ghost[bc_dir] : &skin_ghost.upper_ghost[bc_dir];
  struct gkyl_range *skin_r = bc_edge==GKYL_LOWER_EDGE? &skin_ghost.lower_skin[bc_dir] : &skin_ghost.upper_skin[bc_dir];
  struct gkyl_range ghost_nosub_r;
  gkyl_range_init(&ghost_nosub_r, ghost_r->ndim, ghost_r->lower, ghost_r->upper);

  // Allocate a boundary flux array.
  struct gkyl_array *bflux, *bflux_ho;
  bflux = mkarr(use_gpu, basis.num_basis, ghost_nosub_r.volume);
  bflux_ho = use_gpu? mkarr(false, bflux->ncomp, bflux->size) : gkyl_array_acquire(bflux);

  // Projection updater for boundary flux.dist-function
  gkyl_proj_on_basis *proj_bflux;
  if (cdim == 1 && vdim == 1) {
    assert(false); // NYI
  } else if (cdim == 1 && vdim == 2) {
    assert(false); // NYI
  } else if (cdim == 2 && vdim == 2) {
    proj_bflux = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1,
      evalFunc_2x2v, NULL);
  } else if (cdim == 3 && vdim == 2) {
    proj_bflux = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1,
      evalFunc_3x2v, NULL);
  }

  // Project boundary flux onto basis.
  gkyl_proj_on_basis_advance(proj_bflux, 0.0, &ghost_nosub_r, bflux_ho);
  gkyl_array_copy(bflux, bflux_ho);

  // Create and apply BC at the lower x-boundary with displacing the boundary flux in the y direction.
  struct gkyl_bc_return_flux_gyrokinetic *bclo = gkyl_bc_return_flux_gyrokinetic_new(bc_dir, bc_edge,
    disp_dir, &grid, &ghost_nosub_r, skin_r, use_gpu);

  // Copy ghost cells into buffer.
  gkyl_bc_return_flux_gyrokinetic_advance(bclo, bflux, distf);

  // Check the lower skin cells of distf after applying BC.
  gkyl_array_copy(distf_ho, distf);
  int gidxs[ndim]; // Shifted ghost index.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, skin_r);
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<ndim; d++) gidxs[d] = iter.idx[d];
    gidxs[bc_dir] = bc_edge==GKYL_LOWER_EDGE? iter.idx[bc_dir]-1 : iter.idx[bc_dir]+1;
    gidxs[disp_dir] = iter.idx[disp_dir] < cells[disp_dir]/2? iter.idx[disp_dir]+cells[disp_dir]/2
                                                            : iter.idx[disp_dir]-cells[disp_dir]/2;

    int linidx_skin  = gkyl_range_idx(skin_r, iter.idx);
    int linidx_ghost = gkyl_range_idx(&ghost_nosub_r, gidxs);

    const double *bflux_d = gkyl_array_cfetch(bflux_ho, linidx_ghost);
    const double *distf_d = gkyl_array_cfetch(distf_ho, linidx_skin);

    for (int i = 0; i < basis.num_basis; i++) {
      TEST_CHECK(gkyl_compare(bflux_d[i], distf_d[i], 1e-12));
    }
  }

  gkyl_bc_return_flux_gyrokinetic_release(bclo);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(bflux);
  gkyl_array_release(bflux_ho);
  gkyl_proj_on_basis_release(proj_bflux);
}

void test_bc_return_flux_gyrokinetic_2x2v_p1_ho() { test_bc(2, 2, 1, false); }
void test_bc_return_flux_gyrokinetic_3x2v_p1_ho() { test_bc(3, 2, 1, false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_return_flux_gyrokinetic_2x2v_p1_dev() { test_bc(2, 2, 1, true); }
void test_bc_return_flux_gyrokinetic_3x2v_p1_dev() { test_bc(3, 2, 1, true); }
#endif

TEST_LIST = {
  {"test_bc_return_flux_gyrokinetic_2x2v_p1_ho", test_bc_return_flux_gyrokinetic_2x2v_p1_ho},
  {"test_bc_return_flux_gyrokinetic_3x2v_p1_ho", test_bc_return_flux_gyrokinetic_3x2v_p1_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_bc_return_flux_gyrokinetic_2x2v_p1_dev", test_bc_return_flux_gyrokinetic_2x2v_p1_dev},
  {"test_bc_return_flux_gyrokinetic_3x2v_p1_dev", test_bc_return_flux_gyrokinetic_3x2v_p1_dev},
#endif
  {NULL, NULL},
};
