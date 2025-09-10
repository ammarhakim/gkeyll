#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_efit.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_tok_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_nodal_ops.h>

void
write_geometry(gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_basis basis, struct gkyl_range local, const char *name)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->mc2p, fileNm);
  sprintf(fileNm, fmt, name, "bmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag, fileNm);
  sprintf(fileNm, fmt, name, "g_ij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->g_ij, fileNm);
  sprintf(fileNm, fmt, name, "dxdz");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->dxdz, fileNm);
  sprintf(fileNm, fmt, name, "dzdx");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->dzdx, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo_inv, fileNm);
  sprintf(fileNm, fmt, name, "gij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gij, fileNm);
  sprintf(fileNm, fmt, name, "b_i");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->b_i, fileNm);
  sprintf(fileNm, fmt, name, "bcart");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bcart, fileNm);
  sprintf(fileNm, fmt, name, "cmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->cmag, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobtot, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobtot_inv, fileNm);
  sprintf(fileNm, fmt, name, "bmag_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag_inv, fileNm);
  sprintf(fileNm, fmt, name, "bmag_inv_sq");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag_inv_sq, fileNm);
  sprintf(fileNm, fmt, name, "gxxj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxxj, fileNm);
  sprintf(fileNm, fmt, name, "gxyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxyj, fileNm);
  sprintf(fileNm, fmt, name, "gyyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gyyj, fileNm);
  sprintf(fileNm, fmt, name, "gxzj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxzj, fileNm);
  sprintf(fileNm, fmt, name, "eps2");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->eps2, fileNm);

  // Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &local, 1);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &local, 3, mc2p_nodal, up->mc2p);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
  sprintf(fileNm, fmt, name, "nodes");
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, fileNm);
  gkyl_array_release(mc2p_nodal);
}






// Z is constant at 8.429
// R goes from 4.9 to 5.9
void horizontal_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = 8.429;
}

// Z is constant at -8.429
// R goes from 4.9 to 5.3
void horizontal_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = -8.429;
}

// R is constant at 4.9
// R goes from 8.0 to 8.5
void vertical_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = 8.0 + s/2;
}

void vertical_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = -8.0 - s/2;
}

// Actual Plate info:
// p1 = [5.7413,8.5471]
// p2 = [5.8544, 8.5229]
// p3 = [5.8549, 8.4258]
// Try a Different (slanted) plate instead
// p1 [5.151,8.516]
// p2 [5.852, 8.434]
// Reasonable Diagonal plate:
void shaped_pfunc_upper(double s, double* RZ){
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = 8.516 + (8.434 - 8.516)*s;
}

void shaped_pfunc_lower(double s, double* RZ){
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = -(8.516 + (8.434 - 8.516)*s);
}

//X-pt info
//upper xpt R, Z = 2.50498, 6.14226
//lower xpt R, Z = 2.50504, -6.14214





void
test_fixed_z()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
      .reflect = true,
    };


  double psisep = 1.5098198350000001;
  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  int cpoly_order = 1;
  struct gkyl_basis cbasis;

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rright = 6.2,
    .rleft = 1.1,
    .zmin = -8.3,
    .zmax = 8.3,
    .rmin = 1.1,
    .rmax = 6.2,
  }; 

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  //write_geometry(up, cgrid, cbasis, clocal, "step_outboard_fixed_z");

  // Check that |bhat|=1 at nodes
  int nodes[] = { 1, 1, 1 };
  for (int d=0; d<cgrid.ndim; ++d)
    nodes[d] = cgrid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, cgrid.ndim, nodes);
  struct gkyl_array* bhat_nodal_fld = gkyl_array_new(GKYL_DOUBLE, cgrid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&cbasis, &cgrid, false);
  gkyl_nodal_ops_m2n(n2m, &cbasis, &cgrid, &nrange, &clocal, 3, bhat_nodal_fld, up->bcart);
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  for(int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
      for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
          for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
              cidx[PSI_IDX] = ip;
              cidx[AL_IDX] = ia;
              cidx[TH_IDX] = it;
              double *bhat_n = gkyl_array_fetch(bhat_nodal_fld, gkyl_range_idx(&nrange, cidx));
              double bhat_mag = sqrt(bhat_n[0]*bhat_n[0] + bhat_n[1]*bhat_n[1] + bhat_n[2]*bhat_n[2]);
              TEST_CHECK( gkyl_compare( bhat_mag, 1.0, 1e-12) );
          }
      }
  }

  gkyl_array_release(bhat_nodal_fld);
  gkyl_nodal_ops_release(n2m);
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_horizontal_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {1.3, 0.01, M_PI-1e-14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  int cpoly_order = 1;
  struct gkyl_basis cbasis;

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .plate_spec = true,
    // can set plate func if you dont want a fixed zmin and zmax
    .plate_func_lower = horizontal_pfunc_lower,
    .plate_func_upper = horizontal_pfunc_upper,
  }; 

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 



  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_vertical_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {1.4688, 0.01, M_PI-1e-14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  int cpoly_order = 1;
  struct gkyl_basis cbasis;

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .plate_spec = true,
    // can set plate func if you dont want a fixed zmin and zmax
    .plate_func_lower = vertical_pfunc_lower,
    .plate_func_upper = vertical_pfunc_upper,
  }; 

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_shaped_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };


  double psisep = 1.5098198350000001;
  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  int cpoly_order = 1;
  struct gkyl_basis cbasis;

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .plate_spec = true,
    // can set plate func if you dont want a fixed zmin and zmax
    .plate_func_lower = shaped_pfunc_lower,
    .plate_func_upper = shaped_pfunc_upper,
  }; 

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 

  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

TEST_LIST = {
  { "test_fixed_z", test_fixed_z},
  //{ "test_horizontal_plate", test_horizontal_plate},
  //{ "test_vertical_plate", test_vertical_plate},
  //{ "test_shaped_plate", test_shaped_plate},
  { NULL, NULL },
};
