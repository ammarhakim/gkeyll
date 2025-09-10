#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_rect_decomp.h>

static void
test_wham(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[] = { 2.0e-6, 0.0, -2.0 };
  double cupper[] = { 3.0e-3, 2*M_PI, 2.0 };
  int cells[] = { 10, 16, 32 };

  const char *fname = "core/data/unit/wham_hires.geqdsk_psi.gkyl";
  
  // computational grid
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, 3, clower, cupper, cells);

  if (!gkyl_check_file_exists(fname)) {
    fprintf(stderr, "Unable to find file %s!\n", fname);
    goto cleanup;
  }
  
  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, fname);


  // create mirror geometry
  struct gkyl_mirror_grid_gen *mirror_grid =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = fl_coord,
        .include_axis = include_axis,
        .write_psi_cubic = false,
      }
    );

  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&comp_grid, nghost, &ext_range, &range);

  struct gkyl_mirror_geo_gen *mirror_geo = 
    gkyl_mirror_geo_gen_inew(&(struct gkyl_mirror_geo_gen_inp) {
        .comp_grid = &comp_grid,
        .mirror_grid = mirror_grid,
        .range = range,
        .basis = basis,
      }
    );

  // Create nodal ranges. Must extend to 3x for cartesian projection
  struct gkyl_range nodal_range_2x, nodal_range_3x;
  int nodes_2x[2];
  nodes_2x[0] = gkyl_range_shape(&range, 0) + 1;
  nodes_2x[1] = gkyl_range_shape(&range, 2) + 1;
  gkyl_range_init_from_shape(&nodal_range_2x, 2, nodes_2x);

  int nodes_3x[3];
  for (int d=0; d<3; ++d)
    nodes_3x[d] = gkyl_range_shape(&range, d) + 1;
  gkyl_range_init_from_shape(&nodal_range_3x, 3, nodes_3x);

  enum { PSI_IDX, AL_IDX, Z_IDX }; // arrangement of computational coordinates

  double dalpha = comp_grid.dx[AL_IDX];
  double alpha_lo = comp_grid.lower[AL_IDX] + (range.lower[AL_IDX] - range.lower[AL_IDX]) * comp_grid.dx[AL_IDX];

  // Loop over all node locations
  for (int ip=nodal_range_3x.lower[PSI_IDX]; ip<nodal_range_3x.upper[PSI_IDX]; ++ip) {
    for (int ia=nodal_range_3x.lower[AL_IDX]; ia<nodal_range_3x.upper[AL_IDX]; ++ia) {
      for (int iz=nodal_range_3x.lower[Z_IDX]; iz<nodal_range_3x.upper[Z_IDX]; ++iz) {
        // Find our location in the array
        int idx_2x[2] = { ip, iz };
        int idx_3x[3] = { ip, ia, iz };
        long loc_2x = gkyl_range_idx(&nodal_range_2x, idx_2x);
        long loc_3x = gkyl_range_idx(&nodal_range_3x, idx_3x);
        double alpha_curr = alpha_lo + ia*dalpha;

        const struct gkyl_mirror_geo_gen_geom *g =
          gkyl_array_cfetch(mirror_geo->nodes_geom, loc_3x);
        
        const struct gkyl_mirror_grid_gen_geom *grid = 
          gkyl_array_cfetch(mirror_grid->nodes_geom, loc_2x);

        const double *rz = g->rza_coord;
        // check R,Z coordinates
        const double *rz_grid = gkyl_array_cfetch(mirror_grid->nodes_rz, loc_2x);
        TEST_CHECK( gkyl_compare_double(rz[0], rz_grid[0], 1e-14) );
        TEST_CHECK( gkyl_compare_double(rz[1], rz_grid[1], 1e-14) );

        // check tangent vectors
        for (int i=0; i<3; ++i) {
          struct gkyl_vec3 tang_cart = gkyl_vec3_polar_con_to_cart(rz[0], alpha_curr, grid->tang[i]);
          TEST_CHECK( gkyl_compare_double(tang_cart.x[0], g->tang[i].x[0], 1e-14) );
          TEST_CHECK( gkyl_compare_double(tang_cart.x[1], g->tang[i].x[1], 1e-14) );
          TEST_CHECK( gkyl_compare_double(tang_cart.x[2], g->tang[i].x[2], 1e-14) );
        }

        // check dual vectors
        for (int i=0; i<3; ++i) {
          struct gkyl_vec3 dual_cart = gkyl_vec3_polar_con_to_cart(rz[0], alpha_curr, grid->dual[i]);
          TEST_CHECK( gkyl_compare_double(dual_cart.x[0], g->dual[i].x[0], 1e-14) );
          TEST_CHECK( gkyl_compare_double(dual_cart.x[1], g->dual[i].x[1], 1e-14) );
          TEST_CHECK( gkyl_compare_double(dual_cart.x[2], g->dual[i].x[2], 1e-14) );
        }

        // check normal vectors
        for (int i=0; i<3; ++i) {
          TEST_CHECK( gkyl_compare_double(g->normal[i].x[0], g->dual[i].x[0] / g->dualmag[i], 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->normal[i].x[1], g->dual[i].x[1] / g->dualmag[i], 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->normal[i].x[2], g->dual[i].x[2] / g->dualmag[i], 1e-14) );
        }

        // check dual vector magnitudes
        for (int i=0; i<3; ++i) {
          TEST_CHECK( gkyl_compare_double(g->dualmag[i], gkyl_vec3_len(g->dual[i]), 1e-14) );
        }

        // Check metric tensors
        // ignore if ip=0 and include_axis is true

        int count = 0;
        for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
            if (i > j)
              continue;
            double g_ij = gkyl_vec3_dot(g->tang[i], g->tang[j]);
            TEST_CHECK( gkyl_compare_double(g->metric_covar[count], g_ij, 1e-14) );
            TEST_CHECK( gkyl_compare_double(g->metric_covar_neut[count], g_ij, 1e-14) );
            double gij = gkyl_vec3_dot(g->dual[i], g->dual[j]);
            TEST_CHECK( gkyl_compare_double(g->metric_contr[count], gij, 1e-12) );
            TEST_CHECK( gkyl_compare_double(g->metric_contr_neut[count], gij, 1e-12) );
            TEST_MSG("ip = %d, ia = %d, it = %d\n", ip, ia, iz);
            TEST_MSG("g_ij = %e, g^ij = %e\n", g->metric_covar[count], g->metric_contr[count]);
            TEST_MSG("g_ij_geom = %e, g^ij_geom = %e\n", g_ij, gij);
            count++;
          }
        }

        // check magnetic field
        double Bmag = gkyl_vec3_len( grid->B);
        TEST_CHECK( gkyl_compare_double(g->Bmag, Bmag, 1e-14));
        TEST_CHECK( gkyl_compare_double(g->Bmag_inv, 1.0/Bmag, 1e-14) );
        TEST_CHECK( gkyl_compare_double(g->Bmag_inv_sq, 1.0/(Bmag*Bmag), 1e-14) );

        // check B vector
        struct gkyl_vec3 B_covar = gkyl_vec3_polar_con_to_cov(rz[0], grid->B);
        struct gkyl_vec3 b_covar = gkyl_vec3_norm(B_covar);
        TEST_CHECK( gkyl_compare_double(g->b_covar.x[0], b_covar.x[0], 1e-14) );
        TEST_CHECK( gkyl_compare_double(g->b_covar.x[1], b_covar.x[1] * rz[0]*rz[0], 1e-14) ); // B_1 = B^1 * r^2
        TEST_CHECK( gkyl_compare_double(g->b_covar.x[2], b_covar.x[2], 1e-14) );

        struct gkyl_vec3 B_cart = gkyl_vec3_polar_con_to_cart(rz[0], alpha_curr, grid->B);
        struct gkyl_vec3 b_cart = gkyl_vec3_norm(B_cart);
        TEST_CHECK( gkyl_compare_double(g->b_cart.x[0], b_cart.x[0], 1e-14) );
        TEST_CHECK( gkyl_compare_double(g->b_cart.x[1], b_cart.x[1], 1e-14) );
        TEST_CHECK( gkyl_compare_double(g->b_cart.x[2], b_cart.x[2], 1e-14) );


        // check Jacobian
        double Jac = gkyl_vec3_triple(g->tang[0], g->tang[1], g->tang[2]);

        if (g->rza_coord[0] > 0){
          TEST_CHECK( gkyl_compare_double(Jac, g->Jc, 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->Jc_inv, 1.0/Jac, 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->JB, g->Jc*Bmag, 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->JB_inv, 1.0/(g->Jc*Bmag), 1e-14) );
        }

        // check C = Jc*Bmag/sqrt(g33)
        double g33 = gkyl_vec3_dot( g->tang[2], g->tang[2]);

        if (g->rza_coord[0] > 0) {
          TEST_CHECK( gkyl_compare_double(g->C, g->JB/sqrt(g33), 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->C, g->Jc*Bmag/sqrt(g33), 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->C, g->Jc*g->Bmag/sqrt(g33), 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->C, g->Jc*g->Bmag/sqrt(g->metric_covar[5]), 1e-14) );
          TEST_CHECK( gkyl_compare_double(g->eps2, g->Jc*(g->metric_contr[5] - 1/g->metric_covar[5]), 1e-14) );
        }
      }
    }
  }

  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_mirror_geo_gen_release(mirror_geo);
  gkyl_array_release(psi);

  cleanup:

  return;
}

static void
test_wham_no_axis_psi(void)
{
  test_wham(false, GKYL_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_wham_with_axis_psi(void)
{
  test_wham(true, GKYL_MIRROR_GRID_GEN_PSI_CART_Z);
}

static void
test_wham_no_axis_sqrt_psi(void)
{
  test_wham(false, GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}

static void
test_wham_with_axis_sqrt_psi(void)
{
  test_wham(true, GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z);
}


TEST_LIST = {
  { "wham_no_axis_psi", test_wham_no_axis_psi },
  { "wham_with_axis_psi", test_wham_with_axis_psi },
  
  { "wham_no_axis_sqrt_psi", test_wham_no_axis_sqrt_psi },
  { "wham_with_axis_sqrt_psi", test_wham_with_axis_sqrt_psi },
  
  { NULL, NULL },
};
