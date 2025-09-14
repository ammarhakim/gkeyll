#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

void shaped_pfunc_lower_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = -8.29;
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = 8.29;
}

void shaped_pfunc_upper_inner(double s, double* RZ){
    RZ[0] = 1.651 + (1.8 - 1.651)*s;
    RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
    RZ[0] = 1.65 + (1.8 - 1.65)*s;
    RZ[1] = -(6.33 + (6.777 - 6.33)*s);
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
    // psiRZ and related inputs
    .filepath = "gyrokinetic/data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  double psi_up_core = 1.8;
  double psi_up_pf = 1.8;
  double psi_lo_outer_sol = 0.934;
  double psi_lo_inner_sol = 1.45;

  int npsi_outer_sol = 4;
  int npsi_core = 4;
  int npsi_inner_sol = 4;
  int npsi_lower_pf = 4;
  int npsi_upper_pf = 4;

  double ntheta_lower  = 8;
  double ntheta_middle = 8;
  double ntheta_upper  = 8;

  // Note that for tokamak multi-block simulations, 
  // these theta limits are just placeholders and will be 
  // reset in the multi-block app.
  double theta_lo = -M_PI + 1e-14, theta_up = M_PI - 1e-14;

  // block 0. Lower outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_lower_pf, ntheta_lower },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_R,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zmin_right = -8.29,
          .zmin_left = -6.34,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 1. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_lower},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_LO,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zmin = -8.29,
          .zmax = 8.29,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 2. Middle outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_MID,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zmin = -8.29,
          .zmax = 8.29,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 10, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 3. Upper outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_UP,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zmin = -8.29,
          .zmax = 8.29,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // block 4. Upper outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_upper_pf, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_R,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zmax_right = 8.29,
          .zmax_left = 6.34,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );
  
  // block 5. Upper inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_upper_pf, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_L,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zmax_right = 8.29,
          .zmax_left = 6.34,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
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
      .upper = { psisep, theta_up },
      .cells = { npsi_inner_sol, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_UP,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zmin = -6.34,  
          .zmax = 6.34,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
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
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_MID,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zmin = -6.34,  
          .zmax = 6.34,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 11, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
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
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_lower},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_LO,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zmin = -6.34,  
          .zmax = 6.34,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 9, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_lower_pf, ntheta_lower },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_L,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zmin_right = -8.29,
          .zmin_left = -6.34,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 8, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // block 10. Outer core.
  gkyl_gk_block_geom_set_block(bgeom, 10, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_core, theta_up},
      .cells = { npsi_core, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_R,
          .rclose = 6.2,
          .rleft= 1.1,
          .rright= 6.2,
          .rmin=1.5,
          .rmax=6.2,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 11. Outer core.
  gkyl_gk_block_geom_set_block(bgeom, 11, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_core, theta_up},
      .cells = { npsi_core, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_L,
          .rclose = 1.1,
          .rleft= 1.1,
          .rright= 6.2,
          .rmin=1.5,
          .rmax=6.2,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
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
  double massAr; // Argon mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double vtIon;
  double vtElc;
  double vtAr;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // Argon reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
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
  double vpar_max_Ar; // Velocity space extents in vparallel for Ar
  double mu_max_Ar; // Velocity space extents in mu for Ar

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
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
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 100*2.8*eV;
  double Ti = 150*2.8*eV;
  double TAr = 40.0*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3
  double n0Ar = n0*0.0001/3.0; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(TAr/mAr);

  // Source parameters.
  double nsource = 3.9e23/2.8; // peak source rate in particles/m^3/s 
  double T_source = 285*eV*2.8;
  double cz = 1.5;
  double cx = 0.0065612*9*5;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 18.*mAr*vtAr*vtAr/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 8;

  double t_end = 5.0e-7; 
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
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
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti, 
    .TAr = TAr, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtAr = vtAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar,
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

static inline double sq(double x) { return x*x; }

void
initDensity(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->n0;
  fout[0] = n;
}

void
initDensity_0(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z - M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_1(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z - M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_2(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_3(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z + M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_4(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z + M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_5(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z - M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_6(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z - M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_7(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_8(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z + M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_9(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*0.33399718598613176;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.5093)/sq(2*cx)) * exp(-sq(z + M_PI)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_10(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*1.2645802513482558;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.8)/sq(2*cx)) * exp(-sq(z)/sq(2*cz));
  fout[0] = n;
}

void
initDensity_11(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double cx = input->cx;
  double cz = input->cz;
  double n0 = input->n0*1.2645802513482558;
  double x = xn[0];
  double z = xn[1];
  double n = n0*exp(-sq(x - 1.8)/sq(2*cx)) * exp(-sq(z)/sq(2*cz));
  fout[0] = n;
}

void
initDensityImpurity(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = 1.0e5;
}


void
initDensityNeutral(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->n0Ar;
  fout[0] = n;
}

void
initTempElc(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Te;
  fout[0] = T;
}

void
initTempIon(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Ti;
  fout[0] = T;
}

void
initTempAr(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->TAr;
  fout[0] = T;
}

void
initUpar(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
initUDriftNeutral(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
sourceDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->nsource;
  fout[0] = n;
}

void
sourceUpar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
sourceTemp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->T_source;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuIon;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_multib_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_multib_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_multib_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_multib_app_write_field_energy(app);
      gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
      gkyl_gyrokinetic_multib_app_write_dt(app);
    }
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_multib_app_write_phase(app, t_curr, frame);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  // construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempElc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };


  struct gkyl_gyrokinetic_block_physical_bcs elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 2,
      .collide_with = { "ion", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .radiation = {
      .radiation_id = GKYL_GK_RADIATION, 
      .num_cross_collisions = 1, 
      .collide_with = { "Ar1" },
      .z = 18,
      .charge_state = 1,
      .num_of_densities = 1, // Must be 1 for now
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1", // ion is always the higher charge state
          .donor_nm = "Ar0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    }, 

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 20,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[12];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_0,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_1,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_2,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_3,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_4,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_5,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_6,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_7,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[8] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 8,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_8,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[9] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 9,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_9,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[10] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 10,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_10,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };
  ion_blocks[11] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 11,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity_11,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };


  struct gkyl_gyrokinetic_block_physical_bcs ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
  
    .duplicate_across_blocks = false,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Ar+1 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb Ar1_blocks[1];
  Ar1_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0Ar,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensityImpurity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempAr,
    },

  };

  struct gkyl_gyrokinetic_block_physical_bcs Ar1_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = {  ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0Ar, // Density used to calculate coulomb logarithm
      .T_ref = ctx.TAr, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "ion" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .donor_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },

  
    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Neutral Ar0 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_neut_species_pb Ar0_blocks[1];
  Ar0_blocks[0] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensityNeutral,
      .ctx_udrift = &ctx,
      .udrift = initUDriftNeutral,
      .ctx_temp = &ctx,
      .temp = initTempAr,
    },

  };

  struct gkyl_gyrokinetic_block_physical_bcs Ar0_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_neut_species Ar0 = {
    .name = "Ar0",
    .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = {  ctx.vpar_max_Ar,  ctx.vpar_max_Ar,  ctx.vpar_max_Ar },
    .cells = { cells_v[0], cells_v[0], cells_v[0] },
    .is_static = true,
    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2},

    .duplicate_across_blocks = true,
    .blocks = Ar0_blocks,
    .num_physical_bcs = 20,
    .bcs = Ar0_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
  };

  struct gkyl_gyrokinetic_block_physical_bcs field_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    //{ .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    //{ .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    //{ .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    //{ .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    //{ .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET },
    //{ .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET },
    //{ .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET },
    //{ .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NONE},
    // block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NEUMANN},
    // block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NEUMANN},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 12,
    .bcs = field_phys_bcs,
  };



  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "gk_multib_step_2x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,

    .gk_block_geom = bgeom,
    .cfl_frac = 0.9,
    
    .enforce_positivity = false,

    .num_species = 3,
    .species = { elc, ion, Ar1},

    .num_neut_species = 1,
    .neut_species = { Ar0 },

    .field = field,

    .comm = comm
  };

  // Create app object.
  gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
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
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, app_args.is_restart, false, -1.0);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, app_args.is_restart, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
//  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
//    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
//    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
//    gkyl_gyrokinetic_multib_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
//
//    if (!status.success) {
//      gkyl_gyrokinetic_multib_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
//      break;
//    }
//
//    t_curr += status.dt_actual;
//    dt = status.dt_suggested;
//
//    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
//    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);
//
//    if (dt_init < 0.0) {
//      dt_init = status.dt_actual;
//    }
//    else if (status.dt_actual < dt_failure_tol * dt_init) {
//      num_failures += 1;
//
//      gkyl_gyrokinetic_multib_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
//      gkyl_gyrokinetic_multib_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
//      gkyl_gyrokinetic_multib_app_cout(app, stdout, " num_failures = %d\n", num_failures);
//      if (num_failures >= num_failures_max) {
//        gkyl_gyrokinetic_multib_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
//        gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
//        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
//        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
//        break;
//      }
//    }
//    else {
//      num_failures = 0;
//    }
//
//    step += 1;
//  }

  gkyl_gyrokinetic_multib_app_stat_write(app);

  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_multib_app_stat(app);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-3 failures %ld.\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of write calls %ld.\n", stat.n_io);
  gkyl_gyrokinetic_multib_app_print_timings(app, stdout);

freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);
  gkyl_gk_block_geom_release(bgeom);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
