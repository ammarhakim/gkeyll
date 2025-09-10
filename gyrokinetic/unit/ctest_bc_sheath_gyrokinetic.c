// Test creation and deallocation of updater that applies the
// conducting sheath BC for gyrokinetics.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>

void
test_bc_sheath_gyrokinetic(bool use_gpu)
{

  double charge = -1., mass = 1.;

  int poly_order = 1;
  double lower[] = {-2.0, -2.0, 0.}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {6, 8, 4};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = 1;
  int dir = 0;
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  int vdim = ndim-cdim;
  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  // Basis functions.
  struct gkyl_basis *basis;
  if (poly_order == 1)  // Force gkhybrid for p=1.
    basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                   : gkyl_cart_modal_gkhybrid_new(cdim, vdim);
  else
    basis = use_gpu? gkyl_cart_modal_serendip_cu_dev_new(ndim, poly_order)
                   : gkyl_cart_modal_serendip_new(ndim, poly_order);

  int ghost[] = { 1, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Need the configuration space range to index into phi.
  struct gkyl_range conf_r;
  int rlo[cdim], rup[cdim];
  for (int d=0; d<cdim; d++) {
    rlo[d] = local_ext.lower[d];
    rup[d] = local_ext.upper[d];
  }
  gkyl_range_init(&conf_r, cdim, rlo, rup);

  double q2Dm = 2.*charge/mass;

  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, q2Dm, use_gpu);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);
}

void test_bc_gksheath(){ test_bc_sheath_gyrokinetic(false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_gksheath_cu(){ test_bc_sheath_gyrokinetic(true); }
#endif

TEST_LIST = {
  { "test_bc_gksheath", test_bc_gksheath },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_gksheath_cu", test_bc_gksheath_cu },
#endif
  { NULL, NULL },
};
