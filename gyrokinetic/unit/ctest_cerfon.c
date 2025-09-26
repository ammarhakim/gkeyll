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
#include <gkyl_tok_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>


  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/cerfon.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
      .reflect = true,
    };


void
test_11()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_CORE_R,
    .rclose = 6.0,
    .rright= 6.0,
    .rleft = 0.25,

  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_12()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_CORE_L,
    .rclose = 6.0,
    .rleft= 0.25,
    .rright= 6.0,
  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_1()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_PF_LO_R,
    .rright = 6.0,
    .rleft = 0.25,
    .zmin_left = -5.8,
    .zmin_right = -5.8,
  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_2()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_OUT_LO,
    .rright = 6.0,
    .rleft = 0.25,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_3()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_OUT_MID,
    .rright = 6.0,
    .rleft = 0.25,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_4()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_OUT_UP,
    .rright = 6.0,
    .rleft = 0.25,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_7()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_IN_UP,
    .rleft = 0.25,
    .rright= 6.0,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_8()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_IN_MID,
    .rleft = 0.25,
    .rright= 6.0,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_9()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_DN_SOL_IN_LO,
    .rleft = 0.25,
    .rright= 6.0,
    .zmin = -5.8,
    .zmax = 5.8,
  }; 

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_10()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_PF_LO_L,
    .rleft= 0.25,
    .zmin_left = -5.8,
    .zmin_right = -5.8,
    .rright= 6.0,
  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_5()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_PF_UP_R,
    .rright = 6.0,
    .rleft = 0.25,
    .zmax_left = 5.8,
    .zmax_right = 5.8,
  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
test_6()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 10.0,
    .ftype = GKYL_PF_UP_L,
    .rleft = 0.25,
    .rright= 6.0,
    .zmax_left = 5.8,
    .zmax_right = 5.8,
  }; 

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };
  gkyl_gk_geometry_tok_set_grid_extents(efit_inp, ginp, &clower[2], &cupper[2]);

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
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
  { "test_1", test_1}, //cmag not so great
  { "test_2", test_2},
  { "test_3", test_3},
  { "test_11", test_11},
  { "test_12", test_12},
  { "test_4", test_4},
  {"test_7", test_7},
  {"test_8", test_8},
  {"test_9", test_9},
  {"test_10", test_10}, //cmag not so great
  {"test_5", test_5}, //cmag not so great
  {"test_6", test_6}, //cmag pretty good but not great
  { NULL, NULL },
};
