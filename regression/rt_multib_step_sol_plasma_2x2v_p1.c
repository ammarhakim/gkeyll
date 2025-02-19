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

struct gkyl_block_geom*
create_block_geom(void)
{
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);

  /* Block layout and coordinates

   x  
   ^  
   |
   1  +------------------+------------------+------------------+
   |  |b0                |b1                |b2                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
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
    .filepath = "./data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  psisep = 1.4688;
  gkyl_efit_release(efit);
  double psi_up_core = 1.8;
  double psi_up_pf = 1.8;
  double psi_lo_outer_sol = 0.934;
  double psi_lo_inner_sol = 1.45;

  int npsi_outer_sol = 10;

  double ntheta_lower  = 4;
  double ntheta_middle = 8;
  double ntheta_upper  = 4;

  double theta_lo = -M_PI + 1e-14, theta_up = M_PI - 1e-14;

  // block 0. Lower outer SOL.
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_lower},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_SOL_DN_OUT_LO,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zmin = -8.29,
          .zmax = 8.29,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Middle outer SOL.
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_SOL_DN_OUT_MID,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
          .zmin = -8.29,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 2. Upper outer SOL.
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_SOL_DN_OUT_UP,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .plate_spec = false,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
          .zmin = -8.29,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
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
  double Te; // electron temperature
  double Ti; // ion temperature
  double vtIon;
  double vtElc;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
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
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 364*eV;
  double Ti = 534*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);

  // Source parameters.
  double nsource = 3.9e23; // peak source rate in particles/m^3/s 
  double T_source = 1037.0*eV;
  double cz = 0.8;
  double cx = 0.0065612*4;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);


  // Number of cells.
  int Nx = 4;
  int Nz = 4;
  int Nvpar = 12;
  int Nmu = 8;

  double t_end = 1.0e-4; 
  double num_frames = 10;
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
    .Te = Te, 
    .Ti = Ti, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
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
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
  }
}

static void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app,
  double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }
    gkyl_gyrokinetic_multib_app_write(app, t_curr, frame);
    gkyl_gyrokinetic_multib_app_write_mom(app, t_curr, frame);
  }
}

void
initDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-1)
    n = n0*1e-1;
  fout[0] = n;
}

void
sourceDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/nsource < 1e-3)
    n = nsource*1e-3;
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
initUpar(double t, const double *xn, double* restrict fout, void *ctx)
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


int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

struct gkyl_comm *comm = 0;
  if (app_args.use_mpi) {
#ifdef GKYL_HAVE_MPI
    MPI_Init(&argc, &argv);
        comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#endif
  }
  if (comm == 0)
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  // construct block geometry
  struct gkyl_block_geom *bgeom = create_block_geom();
  int nblocks = gkyl_block_geom_num_blocks(bgeom);


  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.
  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

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
        .upar = initUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };


  struct gkyl_block_physical_bcs elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .bmag_mid =  2.51,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 8,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[1];

  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity,
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
        .upar = initUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };

  struct gkyl_block_physical_bcs ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .bmag_mid =  2.51,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
  
    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 8,
    .bcs = ion_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .polarization_bmag = 2.51,
  };

  struct gkyl_block_physical_bcs field_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 6,
    .bcs = field_phys_bcs,
  };



  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "multib_step_sol_plasma_2x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,

    .block_geom = bgeom,
    .cfl_frac = 0.9,
    
    .enforce_positivity = false,

    .num_species = 2,
    .species = { elc, ion},

    .field = field,

    .comm = comm
  };

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
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
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
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;

}
