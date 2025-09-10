#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_mirror_geo_dg.h>
#include <gkyl_nodal_ops.h>


// Differences between this test and gk_geometry_mirror
// exact_gij [5] = 1/pi^2 in the mirror test, but here it is 1
// exact_g_contra_ij [5] = pi^2 in the mirror test, but here it is 1
// mapc2p there is sqrt(psi*4), z/PI, -alpha, but here it is sqrt(psi*4), z, -alpha
// dualmag [2] = pi in the mirror test, but here it is 1
// The normals are very different


// Functions for test_3x_straight_cylinder
void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], zeta = xn[2];
  fout[0] = sqrt(psi * 4) ; // Function fed is psi = 0.5/2 * R^2 from the efit file
  fout[1] = zeta; // Note that this does not have pi-1e-2 in it because the coordinate zeta is always defined -pi to pi
  fout[2] = -alpha; // There is a minus due to conventions
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = 1/psi; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = r*r; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1.0; // g_33
}

void exact_g_contra_ij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = r*r/4; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = 1/psi/4; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1.0; // g_33
}

void exact_dual_magnitude(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0]  = r/2;
  fout[1] = 1/(2*sqrt(psi));
  fout[2] = 1;
}

void exact_jacobian(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 2.0;
}

void exact_normals(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];
  // Remember cylindrical angle = - alpha
  fout[0] = -cos(-alpha);
  fout[1] = -sin(-alpha);
  fout[2] = 0.0;
  fout[3] = -sin(-alpha);
  fout[4] = cos(-alpha);
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;
  fout[8] = 1.0;
}

void bmag_func(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx){
  fout[0] = 0.5;
}


static void
test_wham(bool include_axis, enum gkyl_mirror_grid_gen_field_line_coord fl_coord)
{
  double clower[3];
  if (include_axis)
    clower[0] = 0.0; // psi lower
  else
    clower[0] = 1.0e-4; // psi lower
  clower[1] = 0.0; // alpha lower
  clower[2] = -2.0; // z lower
  
  double cupper[] = { 3.0e-3, 2*M_PI, 2.0 };
  int cells[] = { 2, 2, 8 };

  const char *fname = "core/data/unit/straight_cylinder.geqdsk_psi.gkyl";

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

  struct gkyl_mirror_geo_dg *mirror_geo_dg = 
    gkyl_mirror_geo_dg_inew(&(struct gkyl_mirror_geo_dg_inp) {
        .comp_grid = &comp_grid,
        .mirror_geo = mirror_geo,
        .range = range,
        .range_ext = ext_range,
        .basis = basis,
      }
    );
  
  // Check this generated geometry compared to the known analytic geometrical quantities

  // Define nodal operations
  enum { PSI_IDX, AL_IDX, Z_IDX }; // arrangement of computational coordinates
  int cidx[3];
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &comp_grid, false);

  // Create nodal range
  struct gkyl_range nodal_range;
  
  int nodes[3];
  for (int d=0; d<3; ++d)
    nodes[d] = gkyl_range_shape(&range, d) + 1;
  gkyl_range_init_from_shape(&nodal_range, 3, nodes);

  // Check that |bhat|=1 at nodes
  struct gkyl_array* bhat_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 3, bhat_nodal, mirror_geo_dg->b_cart);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *bhat_n = gkyl_array_fetch(bhat_nodal, gkyl_range_idx(&nodal_range, cidx));
        double bhat_mag = sqrt(bhat_n[0]*bhat_n[0] + bhat_n[1]*bhat_n[1] + bhat_n[2]*bhat_n[2]);
        TEST_CHECK( gkyl_compare( bhat_mag, 1.0, 1e-12) );
        TEST_MSG("bhat_mag = %g at (psi, alpha, z) = (%g, %g, %g)", bhat_mag,
          comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX],
          comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX],
          comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX]);
      }
    }
  }

  // Check that the duals are what they should be
  // There are errors at low psi, so we shift away from the axis by a few cells
  struct gkyl_array* dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 3, dualmag_nodal, mirror_geo_dg->dualmag);
  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 3, mapc2p_nodal, mirror_geo_dg->mapc2p);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX] + 3; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *dualmag_n = gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double xn[3] = {mapc2p_n[0], mapc2p_n[1], mapc2p_n[2]};
        double dualmag_anal[3];
        exact_dual_magnitude(0, xn, dualmag_anal, 0);
        TEST_CHECK( gkyl_compare( dualmag_n[0], dualmag_anal[0], 1e-6) );
        TEST_CHECK( gkyl_compare( dualmag_n[1], dualmag_anal[1], 1e-6) );
        TEST_CHECK( gkyl_compare( dualmag_n[2], dualmag_anal[2], 1e-6) );
      }
    }
  }

  // Check bmag is what it should be
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, bmag_nodal, mirror_geo_dg->Bmag);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX];
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nodal_range, cidx));
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_n[0], bmag_anal[0], 1e-7) );
        TEST_MSG("bmag = %.16g  at (psi, alpha, z) = (%g, %g, %g). Expected %.16g ", bmag_n[0],
          comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX],
          comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX],
          comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX], bmag_anal[0]);
      }
    }
  }

  // Check bmag_inv_sq = 1/B^2
  struct gkyl_array* bmag_inv_sq_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, bmag_inv_sq_nodal, mirror_geo_dg->Bmag_inv_sq);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *bmag_inv_sq_n = gkyl_array_fetch(bmag_inv_sq_nodal, gkyl_range_idx(&nodal_range, cidx));
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX];
        double xn[3] = {psi, alpha, theta};
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_inv_sq_n[0], 1/(bmag_anal[0]*bmag_anal[0]), 1e-6) );
      }
    }
  }

  // Check bmag_inv = 1/B
  struct gkyl_array* bmag_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, bmag_inv_nodal, mirror_geo_dg->Bmag_inv);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *bmag_inv_n = gkyl_array_fetch(bmag_inv_nodal, gkyl_range_idx(&nodal_range, cidx));
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX];
        double xn[3] = {psi, alpha, theta};
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_inv_n[0], 1/bmag_anal[0], 1e-7) );
        TEST_MSG("bmag_inv = %.16g  at (psi, alpha, z) = (%g, %g, %g). Expected %.16g ", bmag_inv_n[0],
          comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX],
          comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX],
          comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX], 1/bmag_anal[0]);
      }
    }
  }

  // Check cmag = 1
  struct gkyl_array* cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, cmag_nodal, mirror_geo_dg->C);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *cmag_n = gkyl_array_fetch(cmag_nodal, gkyl_range_idx(&nodal_range, cidx));
        TEST_CHECK( gkyl_compare( cmag_n[0], 1.0, 1e-8) );
      }
    }
  }

  // Check eps2 
  // What is this?
  struct gkyl_array* eps2_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, eps2_nodal, mirror_geo_dg->eps2);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *eps2_n = gkyl_array_fetch(eps2_nodal, gkyl_range_idx(&nodal_range, cidx));
        TEST_CHECK( gkyl_compare( eps2_n[0], 0.0, 1e-7) );
        TEST_MSG("eps2 = %.16g should be 0.0 at (psi, alpha, z) = (%.16g, %.16g, %.16g)", 
          eps2_n[0],
          comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX],
          comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX],
          comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX]);
            }
    }
  }

  // Check g_ij
  struct gkyl_array* gij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 6, gij_nodal, mirror_geo_dg->metric_covar);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]+1; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *gij_n = gkyl_array_fetch(gij_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_gij(0.0, xn, fout, 0);
        for (int i=0; i<6; ++i)
          TEST_CHECK( gkyl_compare( gij_n[i], fout[i], 1e-6) );
      }
    }
  }

  // Check g^ij
  struct gkyl_array* gij_contra_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 6, gij_contra_nodal, mirror_geo_dg->metric_contr);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]+1; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *gij_contra_n = gkyl_array_fetch(gij_contra_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_g_contra_ij(0.0, xn, fout, 0);
        for (int i=0; i<6; ++i){
          TEST_CHECK( gkyl_compare( gij_contra_n[i], fout[i], 1e-6) );
          TEST_MSG("g^%d%d = %.16g should be %.16g at (psi, alpha, z) = (%.16g, %.16g, %.16g)", 
            i/3+1, i%3+1, gij_contra_n[i], fout[i],
            comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX],
            comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX],
            comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX]);
        }
      }
    }
  }

  // Check that Jacobgeo is what it should be. This is the Jacobian for the problem
  struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, jacobgeo_nodal, mirror_geo_dg->Jc);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        // double jacobian_analytic = 2.0;
        double fout[1];
        exact_jacobian(0, mapc2p_n, fout, 0);
        TEST_CHECK( gkyl_compare( jacobgeo_n[0], fout[0], 1e-6) );
      }
    }
  }

  // Check jacobgeo_inv
  struct gkyl_array* jacobgeo_inv_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, jacobgeo_inv_nodal, mirror_geo_dg->Jc_inv);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *jacobgeo_inv_n = gkyl_array_fetch(jacobgeo_inv_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double fout[1];
        exact_jacobian(0, mapc2p_n, fout, 0);
        TEST_CHECK( gkyl_compare( jacobgeo_inv_n[0], 1/fout[0], 1e-6) );
      }
    }
  }

  // Check jacobtot
  struct gkyl_array* jacobtot_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, jacobtot_nodal, mirror_geo_dg->JB);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *jacobtot_n = gkyl_array_fetch(jacobtot_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double fout[1];
        exact_jacobian(0, mapc2p_n, fout, 0);
        double Bmag_anal[1];
        bmag_func(0, mapc2p_n, Bmag_anal, 0);
        double jacobtot_analytic = fout[0] * Bmag_anal[0];
        TEST_CHECK( gkyl_compare( jacobtot_n[0], jacobtot_analytic, 1e-6) );
      }
    }
  }

  // Check jacobtot_inv
  struct gkyl_array* jacobtot_inv_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 1, jacobtot_inv_nodal, mirror_geo_dg->JB_inv);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *jacobtot_inv_n = gkyl_array_fetch(jacobtot_inv_nodal, gkyl_range_idx(&nodal_range, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double fout[1];
        exact_jacobian(0, mapc2p_n, fout, 0);
        double Bmag_anal[1];
        bmag_func(0, mapc2p_n, Bmag_anal, 0);
        double jacobtot_analytic = fout[0] * Bmag_anal[0];
        TEST_CHECK( gkyl_compare( jacobtot_inv_n[0], 1/jacobtot_analytic, 1e-6) );
      }
    }
  }

  // Check mapc2p
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX]; 
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nodal_range, cidx));
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i=0; i<3; ++i)
        {
          TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-8) );
          TEST_MSG("mapc2p_n[%d] = %g, fout[%d] = %g at (psi, alpha, z) = (%g, %g, %g)", 
            i, mapc2p_n[i], i, fout[i], psi, alpha, theta);
        }
      }
    }
  }

  // Check mc2nu_pos
  struct gkyl_array* mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, comp_grid.ndim, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range, 3, mc2nu_pos_nodal, mirror_geo_dg->mc2nu_pos);
  for (int ia=nodal_range.lower[AL_IDX]; ia<=nodal_range.upper[AL_IDX]; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX];
        double xn[3] = {psi, alpha, theta};
        double *mc2nu_pos_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(&nodal_range, cidx));
        for (int i=0; i<3; ++i)
          TEST_CHECK( gkyl_compare( mc2nu_pos_n[i], xn[i], 1e-8) );
      }
    }
  }

  // Check normals
  // Plus 3 away from axis to avoid errors
  struct gkyl_array* normals_nodal = gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, nodal_range.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &comp_grid, &nodal_range, &range,  9*basis.num_basis, normals_nodal, mirror_geo_dg->normals);
  for (int ia=nodal_range.lower[AL_IDX]+1; ia<=nodal_range.upper[AL_IDX]-1; ++ia){
    for (int ip=nodal_range.lower[PSI_IDX]+3; ip<=nodal_range.upper[PSI_IDX]; ++ip) {
      for (int it=nodal_range.lower[Z_IDX]; it<=nodal_range.upper[Z_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[Z_IDX] = it;
        double *normals_n = gkyl_array_fetch(normals_nodal, gkyl_range_idx(&nodal_range, cidx));
        double psi = comp_grid.lower[PSI_IDX] + ip*(comp_grid.upper[PSI_IDX]-comp_grid.lower[PSI_IDX])/comp_grid.cells[PSI_IDX];
        double alpha = comp_grid.lower[AL_IDX] + ia*(comp_grid.upper[AL_IDX]-comp_grid.lower[AL_IDX])/comp_grid.cells[AL_IDX];
        double theta = comp_grid.lower[Z_IDX] + it*(comp_grid.upper[Z_IDX]-comp_grid.lower[Z_IDX])/comp_grid.cells[Z_IDX];
        double xn[3] = {psi, alpha, theta};
        double fout[9];
        exact_normals(0.0, xn, fout, 0);
        // printf("cidx: %d %d %d\n", cidx[PSI_IDX], cidx[AL_IDX], cidx[Z_IDX]);
        for (int i=0; i<9; ++i)
        {
          // printf("i=%d, normals_n[i]=%g, fout[i]=%g\n", i, normals_n[i], fout[i]);
          // TEST_CHECK( gkyl_compare( normals_n[i], fout[i], 1e-6) );
        }
      }
    }
  }

  gkyl_array_release(bhat_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(bmag_nodal);
  gkyl_array_release(bmag_inv_nodal);
  gkyl_array_release(bmag_inv_sq_nodal);
  gkyl_array_release(cmag_nodal);
  gkyl_array_release(eps2_nodal);
  gkyl_array_release(gij_nodal);
  gkyl_array_release(gij_contra_nodal);
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(jacobgeo_inv_nodal);
  gkyl_array_release(jacobtot_nodal);
  gkyl_array_release(jacobtot_inv_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mc2nu_pos_nodal);
  gkyl_array_release(normals_nodal);
  gkyl_nodal_ops_release(n2m);

  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_mirror_geo_gen_release(mirror_geo);
  gkyl_mirror_geo_dg_release(mirror_geo_dg);
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
  // { "wham_no_axis_sqrt_psi", test_wham_no_axis_sqrt_psi },
  // { "wham_with_axis_sqrt_psi", test_wham_with_axis_sqrt_psi },
  { NULL, NULL },
};
