#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

//void shaped_pfunc_lower_outer(double s, double* RZ){
//  RZ[0] = 1.0+0.7*s;
//  RZ[1] = -2.07993;
//}
//
//void shaped_pfunc_upper_outer(double s, double* RZ){
//  RZ[0] = 1.0+0.7*s;
//  RZ[1] = 2.07993;
//}

void shaped_pfunc_lower_outer(double s, double* RZ){
  if (s < 0.5) {
    RZ[0] = 1.15+0.6*s;
    RZ[1] = -2.282494;
  }
  else {
    RZ[0] = 1.3712897;
    RZ[1] = -2.282494 + 0.6*(s-0.5);
  }
}
void shaped_pfunc_upper_outer(double s, double* RZ){
  if (s < 0.5) {
    RZ[0] = 1.15+0.6*s;
    RZ[1] = 2.282494;
  }
  else {
    RZ[0] = 1.3712897;
    RZ[1] = 2.282494 - 0.6*(s-0.5);
  }
}

void shaped_pfunc_upper_inner(double s, double* RZ){
    RZ[0] = 0.5;
    RZ[1] = 1.0 + s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
    RZ[0] = 0.5;
    RZ[1] = -(1.0 + s);
}

struct gkyl_gk_block_geom*
create_gk_block_geom(void)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 12);

  /* Block layout and coordinates

   x  
   ^  
   |
   4  +------------------+------------------+------------------+
   |  |b1                |b2                |b3                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   3  +------------------+------------------+------------------+
   |  |b0               x|o b10            %|$ b4              |
   |  |lower outer PF   x|o outer core     %|$ upper outer PF  |
   |  |                 x|o                %|$                 |
   |  +------------------+------------------+------------------+
   2  +------------------+------------------+------------------+
   |  |b9               x|o b11            %|$ b5              |
   |  |lower inner PF   x|o inner core     %|$ upper inner PF  |
   |  |                 x|o                %|$                 |
   1  +------------------+------------------+------------------+
   |  |b8                |b7                |b6                |
   |  |lower inner SOL   |middle inner SOL  |upper inner SOL   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+

      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  


  struct gkyl_efit_inp efit_inp = {
    .filepath = "./hl2a.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  //double psisep = -0.1257861755585853;
  double dsep = 0.000;

  //double psi_up_outer_sol = -0.08;
  double psi_up_outer_sol = -0.10289308777929265;
  double psi_lo_outer_sol = psisep + dsep;

  //double psi_lo_core = -0.17;
  // For core equal cells as outboard : -0.17157235111717056
  // For core 1/2 cells as outboard: -0.14867926333787793
  // For core 1/3 cells as outboard: -0.14104823407811373
  double psi_lo_core = -0.14867926333787793;
  double psi_up_core = psisep - dsep;

  //double psi_lo_pf = -0.14;
  double psi_lo_pf = -0.1334172048183495;
  double psi_up_pf = psisep - dsep;

  //double psi_up_inner_sol = -0.11;
  double psi_up_inner_sol = -0.11052411703905686;
  double psi_lo_inner_sol = psisep + dsep;

  int npsi_outer_sol = 6;
  int npsi_core = 6;
  int npsi_pf = 2;
  int npsi_inner_sol = 4;

  double ntheta_lower_inner  = 2;
  double ntheta_middle_inner = 8;
  double ntheta_upper_inner  = 2;

  double ntheta_lower_outer = 4;
  double ntheta_middle_outer = 10;
  double ntheta_upper_outer = 4;

  double zinner = 1.492;
  double zouter = 2.282494;
  double rright_out = 1.27;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

  // block 0. Lower outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_outer },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_R,
          .compression_factor=1.5,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_LO,
          .compression_factor=1.5,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 2. Middle outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psi_up_outer_sol, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_MID,
          .compression_factor=1.5,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 10, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 3. Upper outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol, theta_up},
      .cells = { npsi_outer_sol, ntheta_upper_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_UP,
          .compression_factor=1.5,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
      }
    }
  );

  // block 4. Upper outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_R,
          .compression_factor=1.5,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );
  
  // block 5. Upper inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_L,
          .compression_factor=1.5,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 6, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 6. Upper inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_upper_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_UP,
          .compression_factor=1.5,
          .rleft = 0.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner, 
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 7. Middle inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_middle_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_MID,
          .compression_factor=1.5,
          .rleft = 0.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 11, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 8. Lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 8, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_lower_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_LO,
          .compression_factor=1.5,
          .rleft = 0.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 9, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 9, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_inner },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_L,
          .compression_factor=1.5,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 8, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );


  // block 10. outer core.
  gkyl_gk_block_geom_set_block(bgeom, 10, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_core, theta_lo},
      .upper = { psi_up_core,  theta_up},
      .cells = { npsi_core, ntheta_middle_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_R,
          .compression_factor=1.5,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.0,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE}, // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 11. Inner Core.
  gkyl_gk_block_geom_set_block(bgeom, 11, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_core, theta_lo },
      .upper = { psi_up_core, theta_up },
      .cells = { npsi_core, ntheta_middle_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_L,
          .compression_factor=1.5,
          .rclose = 0.0,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.0,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE}, // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
}

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massH0; // Hydrogen mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TH0; // neutral hydrogen temperature
  double vtIon;
  double vtElc;
  double vtH0;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double n0H0; // neutral hydrogen reference density
  // Source parameters
  double nsource;
  double Tsource;
  double cx;
  double cz;
  double xcenter;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_H0; // Velocity space extents in vparallel for H0
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};



struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = GKYL_PROTON_MASS; // ion mass
  double mH0 = GKYL_PROTON_MASS; // H0 mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge


  double temp_fac = 3.0;

  double Te = 1500.0/temp_fac*eV;
  double Ti = 1500.0/temp_fac*eV;
  double TH0 = 100.0*eV; 
  double B0 = 2.18; // Magnetic field magnitude in Tesla
  double n0 = 2.0e19; // Particle density in 1/m^3
  double n0H0 = n0*1.0e-2; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtH0 = sqrt(TH0/mH0);


  double nsource = 1.675e22*2.5*temp_fac;
  double Tsource = 3000.0*eV/temp_fac;
  double cx = 0.0065612;
  double cz = 0.4916200*1.4;
  double xcenter = -0.10;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 8.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_H0 = 6.0*vtH0;

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 12;

  double t_end = 5.0e-7; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massH0 = mH0,
    .Te = Te, 
    .Ti = Ti, 
    .TH0 = TH0, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtH0 = vtH0,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0H0 = n0H0,
    .nsource = nsource,
    .Tsource = Tsource,
    .cx = cx,
    .xcenter = xcenter,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .vpar_max_H0 = vpar_max_H0, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}


void
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void
init_density_core(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //1.5e20 at inner core and 2e19 at sep
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double slope = -5.678570e+21;
  double intercept = -6.942856e+20;
  double n = slope*x + intercept;
  fout[0] = 0.333*n;
}

void
init_density_outboard_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double cz = 5.640389838180728;
  double n = app->n0H0;
  if (z < 0.0)
    n = n*exp(-cz*(z+M_PI));
  else 
    n = n*exp(-cz*(M_PI-z));
  fout[0] = n;
}

void
init_density_inboard_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double cz = 10.166078433144026;
  double n = app->n0H0;
  if (z < 0.0)
    n = n*exp(-cz*(z+M_PI));
  else 
    n = n*exp(-cz*(M_PI-z));
  fout[0] = n;
}

void
init_density_pfup_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double n = app->n0H0;
  if (z < -1.467030) {
    double cz = 4.125110078248697;
    n = n*exp(-cz*(z+M_PI));
  }
  else {
    double cz = 1.4988763017083815;
    n = n*exp(-cz*(M_PI - z));
  }
  fout[0] = n;
}

void
init_density_pflo_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double n = app->n0H0;
  if (z < 1.467030) {
    double cz = 1.4988763017083815;
    n = n*exp(-cz*(z+M_PI));
  }
  else {
    double cz = 4.125110078248697;
    n = n*exp(-cz*(M_PI - z));
  }
  fout[0] = n;
}




void
init_density_empty_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n = app->n0H0*1e-5;
  fout[0] = n;
}


void
init_density_outer(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  // This is for -0.08 psi_out
  //double slope = -4.324449e+20;
  //double intercept = -3.439560e+19;
  // This is half width psi_out
  double slope = -8.648899e+20;
  double intercept = -8.879119e+19;

  double n = slope*x + intercept;
  fout[0] = 0.333*n;
}

void
init_density_inner(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double slope = -1.297335e+21;
  double intercept = -1.431868e+20;
  double n = slope*x + intercept;
  fout[0] = 0.333*n;
}

void
init_density_pf(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double slope = 1.297335e+21;
  double intercept = 1.831868e+20;
  double n = slope*x + intercept;
  fout[0] = 0.333*n;
}


void
source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double nsource = app->nsource;
  if(x <= -0.145)
    fout[0] = nsource;
  else 
    fout[0] = nsource*1.0e-5;
}


void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_udrift_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
init_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Te;
  fout[0] = T;
}

void
init_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Ti;
  fout[0] = T;
}
void
source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Tsource;
  fout[0] = T;
}


void
init_temp_H0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->TH0;
  fout[0] = T;
}

void
init_nu_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuElc;
}

void
init_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuIon;
}

static inline void
mapc2p_vel_elc(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_elc = app->mu_max_elc;
  double vpar_max_elc = app->vpar_max_elc;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_elc*2.0*pow(cvpar,2);

  mu = mu_max_elc * (cmu * cmu);

  // Set rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
}

static inline void
mapc2p_vel_ion(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_ion = app->mu_max_ion;
  double vpar_max_ion = app->vpar_max_ion;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_ion*2.0*pow(cvpar,2);

  mu = mu_max_ion * (cmu * cmu);

  // Set rescaled ion velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar ; vp[1] = mu;
}


void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_step_ctx *app = ctx;

  fout[0] = 0.6; // Diffusivity [m^2/s].
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, double dt, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_multib_app_save_dt(app, t_curr, dt);
  }
}

static void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_multib_app_write(app, t_curr, frame);
    gkyl_gyrokinetic_multib_app_write_field_energy(app);
    gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
    gkyl_gyrokinetic_multib_app_write_dt(app);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Construct block geometry.
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[12];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[8] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 8,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[9] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 9,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[10] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 10,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

  elc_blocks[11] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 11,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    //block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    //block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },


    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_elc,
      .num_cross_collisions = 1,
      .bmag_mid = 2.18,
      .collide_with = { "ion" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
    }, 

    //.react_neut = {
    //  .num_react = 2,
    //  .react_type = {
    //    { .react_id = GKYL_REACT_IZ,
    //      .type_self = GKYL_SELF_ELC,
    //      .ion_id = GKYL_ION_H,
    //      .elc_nm = "elc",
    //      .ion_nm = "ion", // ion is always the higher charge state
    //      .donor_nm = "H0", // interacts with elc to give up charge
    //      .charge_state = 0, // corresponds to lower charge state (donor)
    //      .ion_mass = ctx.massIon,
    //      .elc_mass = ctx.massElc,
    //    },
    //    { .react_id = GKYL_REACT_RECOMB,
    //      .type_self = GKYL_SELF_ELC,
    //      .ion_id = GKYL_ION_H,
    //      .elc_nm = "elc",
    //      .ion_nm = "ion",
    //      .recvr_nm = "H0",
    //      .charge_state = 0,
    //      .ion_mass = ctx.massIon,
    //      .elc_mass = ctx.massElc,
    //    },
    //  },
    //}, 


    .duplicate_across_blocks = false,
    .blocks = elc_blocks,
    .num_physical_bcs = 20,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[12];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[8] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 8,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[9] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 9,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

 ion_blocks[10] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 10,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };

 ion_blocks[11] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 11,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },

  };


  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},

    //block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    //block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },

    
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_ion,
      .num_cross_collisions = 1,
      .bmag_mid = 2.18,
      .collide_with = { "elc" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
    }, 


    .duplicate_across_blocks = false,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = 2.18,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},

    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},

    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},

    // block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0}},
    // block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0}},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 12, 
    .bcs = field_phys_bcs,
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "gk_multib_hl2a_nonuniform",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,
    .cfl_frac_omegaH = 1.7/2.0,

    .gk_block_geom = bgeom,
    
    .num_species = 2,
    .species = { elc, ion},

    .num_neut_species = 0,
    .neut_species = {  },

    .field = field,

    .comm = comm
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, -1.0, false);
  write_data(&trig_write, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_multib_update(app, dt);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, status.dt_actual, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_multib_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, status.dt_actual, t_curr, true);
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_multib_app_stat_write(app);

  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_multib_app_stat(app);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_multib_app_print_timings(app, stdout);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);
  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
