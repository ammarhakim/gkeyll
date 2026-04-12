// Test the projection onto an FEM basis that is continuous in the parallel
// direction.
//
#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_parproj_couplex.h>
#include <gkyl_array_integrate.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_bc_basic_gyrokinetic.h>

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  // Allocate array (filled with zeros)
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  // Create ghost and skin sub-ranges given a parent range
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  // Apply periodic BCs along parallel direction
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}

static void check_periodicity(struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  // Check continuity along last dim.
  if (basis.poly_order > 1) return;
  int ndim = basis.ndim;
  int pardir = ndim-1;
  const int num_nodes_perp_max = 4; // 3x p=1.
  int num_nodes_perp = 1;
  if (ndim == 2)
    num_nodes_perp = 2;
  else if (ndim == 3)
    num_nodes_perp = 4;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));

  int idx_lo[ndim];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[pardir] == range.lower[pardir]) {
      int *idx_up = iter.idx;
      for (int d=0; d<pardir; d++)
        idx_lo[d] = idx_up[d];
      idx_lo[pardir] = range.upper[pardir];

      long lidx_lo = gkyl_range_idx(&range, idx_lo);
      long lidx_up = gkyl_range_idx(&range, idx_up);

      double *arr_lo = gkyl_array_fetch(field, lidx_lo);
      double *arr_up = gkyl_array_fetch(field, lidx_up);

      double fn_lo[num_nodes_perp_max], fn_up[num_nodes_perp_max];
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node_lo = gkyl_array_cfetch(nodes, num_nodes_perp+i);
        fn_lo[i] = basis.eval_expand(node_lo, arr_lo);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node_up = gkyl_array_cfetch(nodes, i);
        fn_up[i] = basis.eval_expand(node_up, arr_up);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        TEST_CHECK( gkyl_compare(fn_lo[i], fn_up[i], 1e-9) );
        TEST_MSG( "idx_lo=%d, node %d: lower=%g upper=%g diff=%g\n", idx_lo[0], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
      }
    }
  }
  
  gkyl_array_release(nodes);
}

static void check_continuity_par(struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  // Check continuity along last dim.
  if (basis.poly_order > 1) return;
  int ndim = basis.ndim;
  int pardir = ndim-1;
  const int num_nodes_perp_max = 4; // 3x p=1.
  int num_nodes_perp = 1;
  if (ndim == 2)
    num_nodes_perp = 2;
  else if (ndim == 3)
    num_nodes_perp = 4;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));

  int idx_up[ndim];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[pardir] < range.upper[pardir]) {
      int *idx_lo = iter.idx;
      for (int d=0; d<pardir; d++)
        idx_up[d] = idx_lo[d];
      idx_up[pardir] = idx_lo[pardir] + 1;

      long lidx_lo = gkyl_range_idx(&range, idx_lo);
      long lidx_up = gkyl_range_idx(&range, idx_up);

      double *arr_lo = gkyl_array_fetch(field, lidx_lo);
      double *arr_up = gkyl_array_fetch(field, lidx_up);

      double fn_lo[num_nodes_perp_max], fn_up[num_nodes_perp_max];
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node_lo = gkyl_array_cfetch(nodes, num_nodes_perp+i);
        fn_lo[i] = basis.eval_expand(node_lo, arr_lo);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node_up = gkyl_array_cfetch(nodes, i);
        fn_up[i] = basis.eval_expand(node_up, arr_up);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        TEST_CHECK( gkyl_compare(fn_lo[i], fn_up[i], 1e-9) );
        TEST_MSG( "idx_lo=%d, node %d: lower=%g upper=%g diff=%g\n", idx_lo[0], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
      }
    }
  }
  
  gkyl_array_release(nodes);
}

static void check_continuity_perp(struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  // Check continuity along last dim.
  if (basis.poly_order > 1) return;
  int ndim = basis.ndim;
  int pardir = ndim-1;
  int perpdir = 0;
  const int num_nodes_yz_max = 4; // 3x p=1.
  int num_nodes_yz = 1;
  if (ndim == 2)
    num_nodes_yz = 2;
  else if (ndim == 3)
    num_nodes_yz = 4;

//  int node_idx_lo_2d[] = {1, 3};
//  int node_idx_up_2d[] = {0, 2};
//  int node_idx_lo_3d[] = {1, 3, 5, 7};
//  int node_idx_up_3d[] = {0, 2, 4, 6};

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));

  int idx_up[ndim];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[perpdir] < range.upper[perpdir]) {
      int *idx_lo = iter.idx;
      for (int d=0; d<ndim; d++)
        idx_up[d] = idx_lo[d];

      idx_up[perpdir] = idx_lo[perpdir] + 1;

      long lidx_lo = gkyl_range_idx(&range, idx_lo);
      long lidx_up = gkyl_range_idx(&range, idx_up);

      double *arr_lo = gkyl_array_fetch(field, lidx_lo);
      double *arr_up = gkyl_array_fetch(field, lidx_up);

      double fn_lo[num_nodes_yz_max], fn_up[num_nodes_yz_max];
      for (int i=0; i<num_nodes_yz; i++) {
        const double *node_lo = gkyl_array_cfetch(nodes, 2*i+1);
        fn_lo[i] = basis.eval_expand(node_lo, arr_lo);
      }
      for (int i=0; i<num_nodes_yz; i++) {
        const double *node_up = gkyl_array_cfetch(nodes, 2*i);
        fn_up[i] = basis.eval_expand(node_up, arr_up);
      }
      for (int i=0; i<num_nodes_yz; i++) {
        TEST_CHECK( gkyl_compare(fn_lo[i], fn_up[i], 1e-9) );
        if (ndim == 1)
          TEST_MSG( "idx_lo=%d, idx_up=%d, node %d: lower=%g upper=%g diff=%g\n", idx_lo[0], idx_up[0], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
        else if (ndim == 2)
          TEST_MSG( "idx_lo=%d,%d, idx_up=%d,%d, node %d: lower=%g upper=%g diff=%g\n", idx_lo[0], idx_lo[1], idx_up[0], idx_up[1], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
        else if (ndim == 3)
          TEST_MSG( "idx_lo=%d,%d,%d, idx_up=%d,%d,%d, node %d: lower=%g upper=%g diff=%g\n", idx_lo[0], idx_lo[1], idx_lo[2], idx_up[0], idx_up[1], idx_up[2], i, fn_lo[i], fn_up[i], fn_lo[i]-fn_up[i]);
      }
    }
  }
  
  gkyl_array_release(nodes);
}

void check_dirichlet_bc(struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_basis basis,
  enum gkyl_fem_parproj_bc_type bctype, struct gkyl_array *field_dg, struct gkyl_array *field_fem)
{
  // Check that two fields have the same boundary values in last dimension.
  if (basis.poly_order > 1) return;
  int ndim = basis.ndim;
  int pardir = ndim-1;
  const int num_nodes_perp_max = 4; // 3x p=1.
  int num_nodes_perp = 1;
  if (ndim == 2)
    num_nodes_perp = 2;
  else if (ndim == 3)
    num_nodes_perp = 4;

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));

  for (int e=0; e<2; e++) {

    struct gkyl_range perp_range;
    if (e == 0)
      gkyl_range_shorten_from_above(&perp_range, &local, pardir, 1);
    else
      gkyl_range_shorten_from_below(&perp_range, &local, pardir, 1);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);
    while (gkyl_range_iter_next(&iter)) {
      int diri_idx[ndim];
      for (int d=0; d<ndim; d++)
        diri_idx[d] = iter.idx[d];

      if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
        diri_idx[pardir] = e==0? iter.idx[pardir]-1 : iter.idx[pardir]+1;

      long lidx_diri = gkyl_range_idx(&local_ext, diri_idx);
      long lidx_skin = gkyl_range_idx(&local, iter.idx);

      double *arr_dg = gkyl_array_fetch(field_dg, lidx_diri);
      double *arr_fem = gkyl_array_fetch(field_fem, lidx_skin);

      double fn_dg[num_nodes_perp_max], fn_fem[num_nodes_perp_max];

      int off_diri;
      if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
        off_diri = e==0? num_nodes_perp : 0;
      else
        off_diri = e==0? 0 : num_nodes_perp;

      for (int i=0; i<num_nodes_perp; i++) {
        const double *node = gkyl_array_cfetch(nodes, off_diri+i);
        fn_dg[i] = basis.eval_expand(node, arr_dg);
      }

      int off_skin = e==0? 0 : num_nodes_perp;
      for (int i=0; i<num_nodes_perp; i++) {
        const double *node = gkyl_array_cfetch(nodes, off_skin+i);
        fn_fem[i] = basis.eval_expand(node, arr_fem);
      }
      for (int i=0; i<num_nodes_perp; i++) {
        TEST_CHECK( gkyl_compare(fn_dg[i], fn_fem[i], 1e-11) );
        if (ndim == 1)
          TEST_MSG( "idx=%d, node %d: dg=%g fem=%g diff=%g\n", iter.idx[0], i, fn_dg[i], fn_fem[i], fn_dg[i]-fn_fem[i]);
        else if (ndim == 2)
          TEST_MSG( "idx=%d,%d, node %d: dg=%g fem=%g diff=%g\n", iter.idx[0], iter.idx[1], i, fn_dg[i], fn_fem[i], fn_dg[i]-fn_fem[i]);
        else if (ndim == 3)
          TEST_MSG( "idx=%d,%d,%d, node %d: dg=%g fem=%g diff=%g\n", iter.idx[0], iter.idx[1], iter.idx[2], i, fn_dg[i], fn_fem[i], fn_dg[i]-fn_fem[i]);
      }
    }
  }
  gkyl_array_release(nodes);
}

void evalFunc1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin(2.*M_PI*x);
}

void evalFunc1x_dirichlet(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  // Test Dirichlet BCs with something that's not 0 at the boundary.
  fout[0] = cos(2.*M_PI*x);
}

void ghost_from_skin_surf(bool use_gpu, int dim, struct skin_ghost_ranges *sgr,
  struct gkyl_basis *basis, struct gkyl_array *rho)
{
  // The ghost range with the value of the skin at the boundary.
  struct gkyl_array *bc_buffer = mkarr(use_gpu, rho->ncomp, sgr->lower_ghost[dim-1].volume);
  
  struct gkyl_bc_basic_gyrokinetic* bc_op_lo = gkyl_bc_basic_gyrokinetic_new(dim-1, GKYL_LOWER_EDGE,
    GKYL_BC_GK_FIELD_BOUNDARY_VALUE, basis, &sgr->lower_skin[dim-1], &sgr->lower_ghost[dim-1],
    basis->num_basis, dim, use_gpu);
  gkyl_bc_basic_gyrokinetic_advance(bc_op_lo, bc_buffer, rho);
  gkyl_bc_basic_gyrokinetic_release(bc_op_lo);
  
  struct gkyl_bc_basic_gyrokinetic* bc_op_up = gkyl_bc_basic_gyrokinetic_new(dim-1, GKYL_UPPER_EDGE,
    GKYL_BC_GK_FIELD_BOUNDARY_VALUE, basis, &sgr->upper_skin[dim-1], &sgr->upper_ghost[dim-1],
    basis->num_basis, dim, use_gpu);
  gkyl_bc_basic_gyrokinetic_advance(bc_op_up, bc_buffer, rho);
  gkyl_bc_basic_gyrokinetic_release(bc_op_up);
  
  gkyl_array_release(bc_buffer);
}

void
test_1x(int poly_order, enum gkyl_fem_parproj_bc_type bctype, bool use_gpu)
{
  double lower[] = {-0.5}, upper[] = {0.5};
  int cells[] = {4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1,
    bctype==GKYL_FEM_PARPROJ_DIRICHLET_SKIN || bctype==GKYL_FEM_PARPROJ_DIRICHLET_SKIN? evalFunc1x_dirichlet : evalFunc1x,
    NULL);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  struct gkyl_array *rho_ho = use_gpu? mkarr(false, rho->ncomp, rho->size) : gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    // Fill the ghost cell so we can apply Dirichlet BCs.
    ghost_from_skin_surf(false, dim, &skin_ghost, &basis, rho_ho);

  gkyl_array_copy(rho, rho_ho);
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_parproj_couplex_1x_p2_rho_1.gkyl");

  // parallel FEM projection method.
  struct gkyl_fem_parproj_couplex *parproj = gkyl_fem_parproj_couplex_new(&localRange, &basis,
    bctype, 0, 0, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_couplex_set_rhs(parproj, rho, rho);

  // Solve the problem.
  gkyl_fem_parproj_couplex_solve(parproj, phi);
  gkyl_array_copy(phi_ho, phi);

  if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    struct gkyl_array *parbuff = mkarr(false, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
    apply_periodic_bc(parbuff, phi_ho, dim-1, skin_ghost);
    gkyl_array_release(parbuff);
  }
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_parproj_couplex_1x_p2_phi_1.gkyl");

  // Check continuity at cell boundaries.
  check_continuity_par(localRange, basis, phi_ho);

  if (poly_order == 1) {
    if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
      check_dirichlet_bc(localRange, localRange_ext, basis, bctype, rho_ho, phi_ho);
    } else if (bctype == GKYL_FEM_PARPROJ_NONE) {
      // Solution (checked visually, also checked that phi is actually continuous,
      // and checked that visually looks like results in g2):
      const double sol[8] = {-0.9089542445638024, -0.4554124667453318,
                             -0.8488758876834943,  0.4900987222626481,
                              0.8488758876834943,  0.490098722262648 ,
                              0.9089542445638024, -0.4554124667453318};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[0], 1, phi_p[0]);
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[1], 1, phi_p[1]);
      phi_p = gkyl_array_cfetch(phi_ho, 2);
      TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[2], 2, phi_p[0]);
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[3], 2, phi_p[1]);
      phi_p = gkyl_array_cfetch(phi_ho, 3);
      TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[4], 3, phi_p[0]);
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[5], 3, phi_p[1]);
      phi_p = gkyl_array_cfetch(phi_ho, 4);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[6], 4, phi_p[0]);
      TEST_MSG("Expected: %.13e in cell (%d) | Got: %.13e\n", sol[7], 4, phi_p[1]);
    } else if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
      check_periodicity(localRange, basis, phi_ho);

      // Solution (checked visually against g2):
      const double sol[8] = {-0.8638954769035714, -0.498770286141977,
                             -0.8638954769035713,  0.498770286141977,
                              0.8638954769035713,  0.498770286141977,
                              0.8638954769035713, -0.498770286141977};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_ho, 0);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 2);
      TEST_CHECK( gkyl_compare(sol[2], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[3], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 3);
      TEST_CHECK( gkyl_compare(sol[4], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 4);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 5);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
    }
  } if (poly_order == 2) {
    if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
      check_dirichlet_bc(localRange, localRange_ext, basis, bctype, rho_ho, phi_ho);
    } else if (bctype == GKYL_FEM_PARPROJ_NONE) {
      // Solution (checked visually against g2):
      const double sol[12] = {-0.9010465429057769, -0.4272439810948228,  0.0875367707148495,
                              -0.9039382020247494,  0.4172269800703625,  0.08107082435707  ,
                               0.9039382020247495,  0.4172269800703625, -0.0810708243570699,
                               0.9010465429057768, -0.4272439810948229, -0.0875367707148495};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 2);
      TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 3);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 4);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
    } else if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
      // Solution (checked visually against g2):
      const double sol[12] = {-0.9044201452112453, -0.418896480241106,   0.0799931666307734,
                              -0.9044201452112451,  0.418896480241106,   0.0799931666307734,
                               0.904420145211245 ,  0.418896480241106,  -0.0799931666307734,
                               0.9044201452112451, -0.418896480241106,  -0.0799931666307734};
      const double *phi_p;
      phi_p = gkyl_array_cfetch(phi_ho, 0);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 1);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 2);
      TEST_CHECK( gkyl_compare(sol[3], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[4], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[5], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 3);
      TEST_CHECK( gkyl_compare(sol[6], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[7], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[8], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 4);
      TEST_CHECK( gkyl_compare(sol[9], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[10], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[11], phi_p[2], 1e-14) );
      phi_p = gkyl_array_cfetch(phi_ho, 5);
      TEST_CHECK( gkyl_compare(sol[0], phi_p[0], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[1], phi_p[1], 1e-14) );
      TEST_CHECK( gkyl_compare(sol[2], phi_p[2], 1e-14) );
    }
  }

  gkyl_fem_parproj_couplex_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);

}

static void make_common_x_nodes_equal_at_z_boundaries(struct gkyl_range *range, struct gkyl_basis *basis, struct gkyl_array *field)
{
  // Assuming range only has 2 cells in x, this function makes the lower-x
  // nodes of the upper-x cell equal to the upper-x nodes of the lower-x cell,
  // at the upper and lower z cells.
  if (basis->poly_order > 1 || basis->ndim == 1) return;
  int ndim = basis->ndim;
  int pardir = ndim-1;
  int perpdir = 0;
  const int num_nodes_y_max = 2; // 3x p=1.
  int num_nodes_y;
  int node_idx_up[num_nodes_y_max];
  if (ndim == 2) {
    num_nodes_y = 1;
    node_idx_up[0] = 0;
  }
  else if (ndim == 3) {
    num_nodes_y = 2;
    node_idx_up[0] = 0;
    node_idx_up[1] = 2;
  }

  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, basis->num_basis);
  basis->node_list(gkyl_array_fetch(nodes, 0));

  int idx_up[ndim];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    if ( iter.idx[perpdir] == range->lower[perpdir] && 
        (iter.idx[pardir] == range->lower[pardir] || iter.idx[pardir] == range->upper[pardir])) {
      int *idx_lo = iter.idx;
      for (int d=0; d<ndim; d++)
        idx_up[d] = idx_lo[d];

      idx_up[perpdir] = idx_lo[perpdir] + 1;

      long lidx_lo = gkyl_range_idx(range, idx_lo);
      long lidx_up = gkyl_range_idx(range, idx_up);

      double *arr_lo = gkyl_array_fetch(field, lidx_lo);
      double *arr_up = gkyl_array_fetch(field, lidx_up);

      int off = iter.idx[pardir] == range->lower[pardir]? 0 : (ndim == 2 ? 2 : 4);
      double fn_lo[num_nodes_y_max], fn_up[num_nodes_y_max];
//      printf("  idx_lo=%d,%d,%d | idx_up=%d,%d,%d \n",idx_lo[0],idx_lo[1],idx_lo[2],idx_up[0],idx_up[1],idx_up[2]);
      for (int i=0; i<num_nodes_y; i++) {
        const double *node_lo = gkyl_array_cfetch(nodes, off+2*i+1);
        fn_lo[i] = basis->eval_expand(node_lo, arr_lo);
//        printf("    fn_lo[%d]=%.6e\n",i,fn_lo[i]);
      }
      for (int i=0; i<num_nodes_y; i++) {
        const double *node_up = gkyl_array_cfetch(nodes, off+2*i);
        fn_up[i] = basis->eval_expand(node_up, arr_up);
//        printf("    fn_up[%d]=%.6e\n",i,fn_up[i]);
      }
//      printf("Set nodes equal \n");
      for (int i=0; i<num_nodes_y; i++)
        fn_up[i] = fn_lo[i];

      double fn_up_all[basis->num_basis];
      for (int i=0; i<basis->num_basis; i++) {
        const double *node_curr = gkyl_array_cfetch(nodes, i);
        fn_up_all[i] = basis->eval_expand(node_curr, arr_up);
      }
      for (int i=0; i<num_nodes_y; i++) {
        fn_up_all[off+node_idx_up[i]] = fn_up[i];
      }

      // Convert back to modal.
      basis->nodal_to_modal(&fn_up_all[0], arr_up);
    }
  }
  
  gkyl_array_release(nodes);
}

void evalFunc2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double mu = .2;
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu,2))/(2.0*sig*sig))*sin(2.*M_PI*y);
}

void evalFunc2x_xcont(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double mu = .2;
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu,2))/(2.0*sig*sig));
}

void evalFunc2x_ydiscont(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 2.0+sin(2.*M_PI*y);
}

void
evalFunc2x_dirichlet(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  fout[0] = cos(x)*cos(5*z);
}

void
test_2x(int poly_order, enum gkyl_fem_parproj_bc_type bctype, bool use_gpu)
{
  double lower[] = {-2., -0.5}, upper[] = {2., 0.5};
  int cells[] = {3, 4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Create a range that's only 2 cells in x.
  struct gkyl_range solve_rng;
  int srng_lower[GKYL_MAX_CDIM], srng_upper[GKYL_MAX_CDIM];
  for (int d=0; d<dim; d++) {
    srng_lower[d] = localRange.lower[d];
    srng_upper[d] = localRange.upper[d];
  }
  srng_upper[0] = localRange.lower[0]+1;
  gkyl_sub_range_init(&solve_rng, &localRange, srng_lower, srng_upper);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, bctype==GKYL_FEM_PARPROJ_DIRICHLET_SKIN? evalFunc2x_dirichlet : evalFunc2x, NULL);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  struct gkyl_array *rho_ho = use_gpu? mkarr(false, rho->ncomp, rho->size) : gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);

  // Project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
    // Make the RHS source equal along x at the z boundary. 
    make_common_x_nodes_equal_at_z_boundaries(&solve_rng, &basis, rho_ho);
  }

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    // Fill the ghost cell so we can apply Dirichlet BCs.
    ghost_from_skin_surf(false, dim, &skin_ghost, &basis, rho_ho);

  gkyl_array_copy(rho, rho_ho);

//  // Project a function that is continuous in x but discontinuous in z.
//  gkyl_eval_on_nodes *evcont = gkyl_eval_on_nodes_new(&grid, &basis,
//    1, evalFunc2x_xcont, NULL);
//  gkyl_proj_on_basis *projdiscont = gkyl_proj_on_basis_new(&grid, &basis,
//    poly_order+1, 1, evalFunc2x_ydiscont, NULL);
//  gkyl_eval_on_nodes_advance(evcont, 0.0, &localRange, rho_ho);
//  gkyl_proj_on_basis_advance(projdiscont, 0.0, &localRange, phi_ho);
//  gkyl_proj_on_basis_release(projdiscont);
//  gkyl_eval_on_nodes_release(evcont);
//  gkyl_dg_mul_op(basis, 0, rho_ho, 0, phi_ho, 0, rho_ho);
//  gkyl_array_copy(rho, rho_ho);

//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_parproj_couplex_2x_p1_rho_1.gkyl");

  // Parallel FEM projection method.
  struct gkyl_fem_parproj_couplex *parproj = gkyl_fem_parproj_couplex_new(&solve_rng, &basis,
    bctype, 0, 0, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_couplex_set_rhs(parproj, rho, rho);

  // Solve the problem.
  gkyl_fem_parproj_couplex_solve(parproj, phi);
  gkyl_array_copy(phi_ho, phi);

  if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    struct gkyl_array *parbuff = mkarr(false, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
    apply_periodic_bc(parbuff, phi_ho, dim-1, skin_ghost);
    gkyl_array_release(parbuff);
  }
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_parproj_couplex_2x_p1_phi_1.gkyl");

  // Check continuity at cell boundaries.
  check_continuity_par(solve_rng, basis, phi_ho);
  check_continuity_perp(solve_rng, basis, phi_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
    check_dirichlet_bc(solve_rng, localRange_ext, basis, bctype, rho_ho, phi_ho);
  } else if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    // Check periodicity.
    check_periodicity(localRange, basis, phi_ho);
  }

  gkyl_fem_parproj_couplex_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);

}

void evalWeight2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double mu = 0.0;
  double sig = 0.3;
  double Lx = 4.0;

  fout[0] = cos((2.*M_PI/(2*Lx))*x);
  if (y < 0.0)
    fout[0] *= exp(-(pow(y-mu,2))/(2.0*pow(sig,2)));
  else
    fout[0] *= 3.0*exp(-(pow(y-mu,2))/(2.0*pow(sig,2)));
}

void
test_2x_weighted(int poly_order, enum gkyl_fem_parproj_bc_type bctype, bool use_gpu)
{
  double lower[] = {-2., -0.5}, upper[] = {2., 0.5};
  int cells[] = {3, 4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Create a range that's only 2 cells in x.
  struct gkyl_range solve_rng;
  int srng_lower[GKYL_MAX_CDIM], srng_upper[GKYL_MAX_CDIM];
  for (int d=0; d<dim; d++) {
    srng_lower[d] = localRange.lower[d];
    srng_upper[d] = localRange.upper[d];
  }
  srng_upper[0] = localRange.lower[0]+1;
  gkyl_sub_range_init(&solve_rng, &localRange, srng_lower, srng_upper);

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Weight in the projection operation.
  struct gkyl_array *jac = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  struct gkyl_array *rho_ho = use_gpu? mkarr(false, rho->ncomp, rho->size) : gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  struct gkyl_array *jac_ho = use_gpu? mkarr(false, jac->ncomp, jac->size) : gkyl_array_acquire(jac);

  // Project distribution function on basis.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, bctype==GKYL_FEM_PARPROJ_DIRICHLET_SKIN? evalFunc2x_dirichlet : evalFunc2x, NULL);
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
    // Make the RHS source equal along x at the z boundary. 
    make_common_x_nodes_equal_at_z_boundaries(&solve_rng, &basis, rho_ho);
  }

  gkyl_array_copy(rho, rho_ho);
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_parproj_couplex_2x_p1_rho_1.gkyl");

  // Project the weight onto the basis.
  gkyl_eval_on_nodes *proj_weight = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalWeight2x, NULL);
  gkyl_eval_on_nodes_advance(proj_weight, 0.0, &localRange, jac_ho);
  gkyl_array_copy(jac, jac_ho);
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, jac_ho,  "ctest_fem_parproj_couplex_2x_p1_jac_1.gkyl");

  // Parallel FEM projection method.
  struct gkyl_fem_parproj_couplex *parproj = gkyl_fem_parproj_couplex_new(&solve_rng, &basis,
    bctype, jac, jac, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_couplex_set_rhs(parproj, rho, rho);

  // Solve the problem.
  gkyl_fem_parproj_couplex_solve(parproj, phi);
  gkyl_array_copy(phi_ho, phi);

  if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    struct gkyl_array *parbuff = mkarr(false, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
    apply_periodic_bc(parbuff, phi_ho, dim-1, skin_ghost);
    gkyl_array_release(parbuff);
  }
//  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_parproj_couplex_2x_p1_phi_1.gkyl");

  // Check that the field is continuous.
  check_continuity_par(solve_rng, basis, phi_ho);
  check_continuity_perp(solve_rng, basis, phi_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN)
    check_dirichlet_bc(solve_rng, localRange_ext, basis, bctype, rho_ho, phi_ho);

  gkyl_fem_parproj_couplex_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_eval_on_nodes_release(proj_weight);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(jac);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(jac_ho);

}

void evalFunc2x_selfadjoint(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double mu = .2;
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu,2))/(2.0*sig*sig))*(2.0+cos(2.*M_PI*y));
}
void evalGunc2x_selfadjoint(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double mu = .1;
  double sig = 0.4;
  fout[0] = exp(-(pow(x-mu,2))/(2.0*sig*sig))*(2.0+y*y);
}

void
test_2x_selfadjoint(int poly_order, enum gkyl_fem_parproj_bc_type bctype, bool use_gpu)
{
  // Check that the operator is self-adjoint.
  double lower[] = {-2., -0.5}, upper[] = {2., 0.5};
  int cells[] = {3, 4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Create DG fields.
  struct gkyl_array *rho_dg = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_dg = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *prod = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // Create FEM fields.
  struct gkyl_array *rho_fem = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *phi_fem = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  struct gkyl_array *rho_ho = use_gpu? mkarr(false, rho_dg->ncomp, rho_dg->size) : gkyl_array_acquire(rho_dg);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi_dg->ncomp, phi_dg->size) : gkyl_array_acquire(phi_dg);

  // Project fields onto basis.
  gkyl_proj_on_basis *projob_rho = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc2x_selfadjoint, NULL);
  gkyl_proj_on_basis_advance(projob_rho, 0.0, &localRange, rho_ho);
  gkyl_proj_on_basis_release(projob_rho);
  gkyl_array_copy(rho_dg, rho_ho);

  gkyl_proj_on_basis *projob_phi = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalGunc2x_selfadjoint, NULL);
  gkyl_proj_on_basis_advance(projob_phi, 0.0, &localRange, phi_ho);
  gkyl_proj_on_basis_release(projob_phi);
  gkyl_array_copy(phi_dg, phi_ho);

  // Parallel FEM projection method.
  struct gkyl_fem_parproj_couplex *parproj = gkyl_fem_parproj_couplex_new(&localRange, &basis,
    bctype, 0, 0, use_gpu);

  struct gkyl_array_integrate* arr_int_op = gkyl_array_integrate_new(&grid, &basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  // Smooth rho_dg and integrate phi_dg*rho_fem.
  gkyl_fem_parproj_couplex_set_rhs(parproj, rho_dg, rho_dg);
  gkyl_fem_parproj_couplex_solve(parproj, rho_fem);
  gkyl_dg_mul_op(basis, 0, prod, 0, phi_dg, 0, rho_fem);
  double *int_prodA = use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  gkyl_array_integrate_advance(arr_int_op, prod, 1.0, 0, &localRange, 0, int_prodA);
  double int_prodA_ho[1];
  if (use_gpu)
    gkyl_cu_memcpy(int_prodA_ho, int_prodA, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(int_prodA_ho, int_prodA, sizeof(double));

  // Smooth phi_dg and integrate phi_fem*rho_dg.
  gkyl_fem_parproj_couplex_set_rhs(parproj, phi_dg, phi_dg);
  gkyl_fem_parproj_couplex_solve(parproj, phi_fem);
  gkyl_dg_mul_op(basis, 0, prod, 0, phi_fem, 0, rho_dg);
  double *int_prodB = use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  gkyl_array_integrate_advance(arr_int_op, prod, 1.0, 0, &localRange, 0, int_prodB);
  double int_prodB_ho[1];
  if (use_gpu)
    gkyl_cu_memcpy(int_prodB_ho, int_prodB, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(int_prodB_ho, int_prodB, sizeof(double));

  TEST_CHECK( gkyl_compare(int_prodA_ho[0],int_prodB_ho[0], 1e-14) );
  TEST_MSG("int phi_dg*rho_fem = %.13e | int phi_fem*rho_dg = %.13e", int_prodA_ho[0],int_prodB_ho[0]);
//  printf("\nint phi_dg*rho_fem = %.13e | int phi_fem*rho_dg = %.13e\n", int_prodA_ho[0],int_prodB_ho[0]);

  if (use_gpu) {
    gkyl_cu_free(int_prodA);
    gkyl_cu_free(int_prodB);
  }
  else {
    gkyl_free(int_prodA);
    gkyl_free(int_prodB);
  }
  gkyl_array_integrate_release(arr_int_op);
  gkyl_fem_parproj_couplex_release(parproj);
  gkyl_array_release(rho_dg);
  gkyl_array_release(phi_dg);
  gkyl_array_release(prod);
  gkyl_array_release(rho_fem);
  gkyl_array_release(phi_fem);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);

}

void evalFunc3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double mu[2] = {.2, 0.2};
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu[0],2)+pow(y-mu[1],2))/(2.0*sig*sig))*sin(2.*M_PI*z);
}

void evalFunc3x_dirichlet(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double mu[2] = {.2, 0.2};
  double sig = 0.3;
  fout[0] = exp(-(pow(x-mu[0],2)+pow(y-mu[1],2))/(2.0*sig*sig))*cos(2.*M_PI*z);
}

void
test_3x(const int poly_order, enum gkyl_fem_parproj_bc_type bctype, bool use_gpu)
{
  double lower[] = {-2., -2., -0.5}, upper[] = {2., 2., 0.5};
  int cells[] = {3, 3, 4};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Create a range that's only 2 cells in x.
  struct gkyl_range solve_rng;
  int srng_lower[GKYL_MAX_CDIM], srng_upper[GKYL_MAX_CDIM];
  for (int d=0; d<dim; d++) {
    srng_lower[d] = localRange.lower[d];
    srng_upper[d] = localRange.upper[d];
  }
  srng_upper[0] = localRange.lower[0]+1;
  gkyl_sub_range_init(&solve_rng, &localRange, srng_lower, srng_upper);

  // projection updater for DG field.
  gkyl_proj_on_basis *projob = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, bctype==GKYL_FEM_PARPROJ_DIRICHLET_SKIN? evalFunc3x_dirichlet : evalFunc3x, NULL);

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(use_gpu, basis.num_basis, localRange_ext.volume);

  struct gkyl_array *rho_ho = use_gpu? mkarr(false, rho->ncomp, rho->size) : gkyl_array_acquire(rho);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);

  // project distribution function on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
    // Make the RHS source equal along x at the z boundary. 
    make_common_x_nodes_equal_at_z_boundaries(&solve_rng, &basis, rho_ho);
  }

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST)
    // Fill the ghost cell so we can apply Dirichlet BCs.
    ghost_from_skin_surf(false, dim, &skin_ghost, &basis, rho_ho);

  gkyl_array_copy(rho, rho_ho);
  gkyl_grid_sub_array_write(&grid, &localRange, 0, rho_ho, "ctest_fem_parproj_couplex_3x_p1_rho_1.gkyl");

  // parallel FEM projection method.
  struct gkyl_fem_parproj_couplex *parproj = gkyl_fem_parproj_couplex_new(&solve_rng, &basis,
    bctype, 0, 0, use_gpu);

  // Set the RHS source.
  gkyl_fem_parproj_couplex_set_rhs(parproj, rho, rho);

  // Solve the problem.
  gkyl_fem_parproj_couplex_solve(parproj, phi);
  gkyl_array_copy(phi_ho, phi);

  if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    struct gkyl_array *parbuff = mkarr(false, basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
    apply_periodic_bc(parbuff, phi_ho, dim-1, skin_ghost);
    gkyl_array_release(parbuff);
  }
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_ho, "ctest_fem_parproj_couplex_3x_p1_phi_1.gkyl");

  // Check continuity at cell boundaries.
  check_continuity_par(solve_rng, basis, phi_ho);
  check_continuity_perp(solve_rng, basis, phi_ho);

  if (bctype == GKYL_FEM_PARPROJ_DIRICHLET_GHOST || bctype == GKYL_FEM_PARPROJ_DIRICHLET_SKIN) {
    check_dirichlet_bc(solve_rng, localRange_ext, basis, bctype, rho_ho, phi_ho);
  } else if (bctype == GKYL_FEM_PARPROJ_PERIODIC) {
    check_periodicity(solve_rng, basis, phi_ho);
  }

  gkyl_fem_parproj_couplex_release(parproj);
  gkyl_proj_on_basis_release(projob);
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(rho_ho);
  gkyl_array_release(phi_ho);

}

void test_1x_p1_bcnone_ho() {test_1x(1, GKYL_FEM_PARPROJ_NONE, false);}
void test_1x_p1_bcdirichlet_ho() {
  test_1x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_1x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_1x_p1_bcperiodic_ho() {test_1x(1, GKYL_FEM_PARPROJ_PERIODIC, false);}

void test_1x_p2_bcnone_ho() {test_1x(2, GKYL_FEM_PARPROJ_NONE, false);}
void test_1x_p2_bcdirichlet_ho() {
  test_1x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_1x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_1x_p2_bcperiodic_ho() {test_1x(2, GKYL_FEM_PARPROJ_PERIODIC, false);}

void test_2x_p1_bcnone_ho() {test_2x(1, GKYL_FEM_PARPROJ_NONE, false);}
void test_2x_p1_bcdirichlet_ho() {
  test_2x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_2x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_2x_p1_bcperiodic_ho() {test_2x(1, GKYL_FEM_PARPROJ_PERIODIC, false);}
void test_2x_p1_weighted_ho() {test_2x_weighted(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);}
void test_2x_p1_selfadjoint_ho() {test_2x_selfadjoint(1, GKYL_FEM_PARPROJ_NONE, false);}

void test_2x_p2_bcnone_ho() {test_2x(2, GKYL_FEM_PARPROJ_NONE, false);}
void test_2x_p2_bcdirichlet_ho() {
  test_2x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_2x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_2x_p2_bcperiodic_ho() {test_2x(2, GKYL_FEM_PARPROJ_PERIODIC, false);}

void test_3x_p1_bcnone_ho() {test_3x(1, GKYL_FEM_PARPROJ_NONE, false);}
void test_3x_p1_bcdirichlet_ho() {
  test_3x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_3x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_3x_p1_bcperiodic_ho() {test_3x(1, GKYL_FEM_PARPROJ_PERIODIC, false);}

void test_3x_p2_bcnone_ho() {test_3x(2, GKYL_FEM_PARPROJ_NONE, false);}
void test_3x_p2_bcdirichlet_ho() {
  test_3x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, false);
  test_3x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, false);
}
void test_3x_p2_bcperiodic_ho() {test_3x(2, GKYL_FEM_PARPROJ_PERIODIC, false);}

#ifdef GKYL_HAVE_CUDA
// ......... GPU tests ............ //
void test_1x_p1_bcnone_dev() {test_1x(1, GKYL_FEM_PARPROJ_NONE, true);}
void test_1x_p1_bcdirichlet_dev() {
  test_1x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_1x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_1x_p1_bcperiodic_dev() {test_1x(1, GKYL_FEM_PARPROJ_PERIODIC, true);}

void test_1x_p2_bcnone_dev() {test_1x(2, GKYL_FEM_PARPROJ_NONE, true);}
void test_1x_p2_bcdirichlet_dev() {
  test_1x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_1x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_1x_p2_bcperiodic_dev() {test_1x(2, GKYL_FEM_PARPROJ_PERIODIC, true);}

void test_2x_p1_bcnone_dev() {test_2x(1, GKYL_FEM_PARPROJ_NONE, true);}
void test_2x_p1_bcdirichlet_dev() {
  test_2x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_2x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_2x_p1_bcperiodic_dev() {test_2x(1, GKYL_FEM_PARPROJ_PERIODIC, true);}
void test_2x_p1_weighted_dev() {test_2x_weighted(1, GKYL_FEM_PARPROJ_NONE, true);}
void test_2x_p1_selfadjoint_dev() {test_2x_selfadjoint(1, GKYL_FEM_PARPROJ_NONE, true);}

void test_2x_p2_bcnone_dev() {test_2x(2, GKYL_FEM_PARPROJ_NONE, true);}
void test_2x_p2_bcdirichlet_dev() {
  test_2x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_2x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_2x_p2_bcperiodic_dev() {test_2x(2, GKYL_FEM_PARPROJ_PERIODIC, true);}

void test_3x_p1_bcnone_dev() {test_3x(1, GKYL_FEM_PARPROJ_NONE, true);}
void test_3x_p1_bcdirichlet_dev() {
  test_3x(1, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_3x(1, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_3x_p1_bcperiodic_dev() {test_3x(1, GKYL_FEM_PARPROJ_PERIODIC, true);}

void test_3x_p2_bcnone_dev() {test_3x(2, GKYL_FEM_PARPROJ_NONE, true);}
void test_3x_p2_bcdirichlet_dev() {
  test_3x(2, GKYL_FEM_PARPROJ_DIRICHLET_GHOST, true);
  test_3x(2, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, true);
}
void test_3x_p2_bcperiodic_dev() {test_3x(2, GKYL_FEM_PARPROJ_PERIODIC, true);}
#endif

TEST_LIST = {
  { "test_1x_p1_bcnone_ho", test_1x_p1_bcnone_ho },
  { "test_1x_p1_bcdirichlet_ho", test_1x_p1_bcdirichlet_ho },
  { "test_1x_p1_bcperiodic_ho", test_1x_p1_bcperiodic_ho },
//  { "test_1x_p2_bcnone_ho", test_1x_p2_bcnone_ho },
//  { "test_1x_p2_bcdirichlet_ho", test_1x_p2_bcdirichlet_ho },
//  { "test_1x_p2_bcperiodic_ho", test_1x_p2_bcperiodic_ho },
  { "test_2x_p1_bcnone_ho", test_2x_p1_bcnone_ho },
  { "test_2x_p1_bcdirichlet_ho", test_2x_p1_bcdirichlet_ho },
  { "test_2x_p1_bcperiodic_ho", test_2x_p1_bcperiodic_ho },
//  { "test_2x_p2_bcnone_ho", test_2x_p2_bcnone_ho },
//  { "test_2x_p2_bcdirichlet_ho", test_2x_p2_bcdirichlet_ho },
//  { "test_2x_p2_bcperiodic_ho", test_2x_p2_bcperiodic_ho },
  { "test_2x_p1_weighted_ho", test_2x_p1_weighted_ho},
  { "test_2x_p1_selfadjoint_ho", test_2x_p1_selfadjoint_ho},
  { "test_3x_p1_bcnone_ho", test_3x_p1_bcnone_ho },
  { "test_3x_p1_bcdirichlet_ho", test_3x_p1_bcdirichlet_ho },
  { "test_3x_p1_bcperiodic_ho", test_3x_p1_bcperiodic_ho },
//  { "test_3x_p2_bcnone_ho", test_3x_p2_bcnone_ho },
//  { "test_3x_p2_bcdirichlet_ho", test_3x_p2_bcdirichlet_ho },
//  { "test_3x_p2_bcperiodic_ho", test_3x_p2_bcperiodic_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_p1_bcnone_dev", test_1x_p1_bcnone_dev },
  { "test_1x_p1_bcdirichlet_dev", test_1x_p1_bcdirichlet_dev },
  { "test_1x_p1_bcperiodic_dev", test_1x_p1_bcperiodic_dev },
//  { "test_1x_p2_bcnone_dev", test_1x_p2_bcnone_dev },
//  { "test_1x_p2_bcdirichlet_dev", test_1x_p2_bcdirichlet_dev },
//  { "test_1x_p2_bcperiodic_dev", test_1x_p2_bcperiodic_dev },
  { "test_2x_p1_bcnone_dev", test_2x_p1_bcnone_dev },
  { "test_2x_p1_bcdirichlet_dev", test_2x_p1_bcdirichlet_dev },
  { "test_2x_p1_bcperiodic_dev", test_2x_p1_bcperiodic_dev },
//  { "test_2x_p2_bcnone_dev", test_2x_p2_bcnone_dev },
//  { "test_2x_p2_bcdirichlet_dev", test_2x_p2_bcdirichlet_dev },
//  { "test_2x_p2_bcperiodic_dev", test_2x_p2_bcperiodic_dev },
  { "test_2x_p1_weighted_dev_dev", test_2x_p1_weighted_dev},
  { "test_2x_p1_selfadjoint_dev", test_2x_p1_selfadjoint_dev},
  { "test_3x_p1_bcnone_dev", test_3x_p1_bcnone_dev },
  { "test_3x_p1_bcdirichlet_dev", test_3x_p1_bcdirichlet_dev },
  { "test_3x_p1_bcperiodic_dev", test_3x_p1_bcperiodic_dev },
//  { "test_3x_p2_bcnone_dev", test_3x_p2_bcnone_dev },
//  { "test_3x_p2_bcdirichlet_dev", test_3x_p2_bcdirichlet_dev },
//  { "test_3x_p2_bcperiodic_dev", test_3x_p2_bcperiodic_dev },
#endif
  { NULL, NULL },
};

