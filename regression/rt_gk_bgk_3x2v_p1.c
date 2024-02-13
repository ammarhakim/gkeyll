#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

struct gk_lapd_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double B0; // reference magnetic field
  double n0; // reference density
  // Source parameters
  double Te_source; // Source electron temperature
  double r_s; // Radial extent of source
  double L_s; // Length of source
  // Simulation parameters
  double Lx; // Box size in x
  double Ly; // Box size in y
  double Lperp; // Perpendicular length
  double Lz; // Box size in z
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double finalTime; // end time
  int numFrames; // number of output frames
};

static inline double
profileA(double r, double c_edge, double Lperp)
{
  if (r < Lperp/2.0)
    return (1.0-c_edge)*pow(1.0-(r/(Lperp/2.))*(r/(Lperp/2.)), 3.0) + c_edge;
  else
    return c_edge;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double r = sqrt(x*x + y*y);
  double n0 = app->n0;
  double Lperp = app->Lperp;
  pcg32_random_t rng; // RNG for use in IC
  double perturb = 2e-3*(1.0 - 0.5*gkyl_pcg32_rand_double(&rng));
  double n = n0*profileA(r, 1.0/20.0, Lperp)*(1.0 + perturb);
  fout[0] = n;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double r = sqrt(x*x + y*y);
  double massElc = app->massElc;
  double Lperp = app->Lperp;
  double T = app->Te*profileA(r, 1.0/5.0, Lperp);
  fout[0] = T;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double r = sqrt(x*x + y*y);
  double massIon = app->massIon;
  double Lperp = app->Lperp;
  double T = app->Ti;
  fout[0] = T;

}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double r = sqrt(x*x + y*y);
  double n0 = app->n0;
  double Lz = app->Lz;
  double c_s = app->c_s;
  double r_s = app->r_s;
  double L_s = app->L_s;
  double sourceFloor = 0.01;
  double S0 = 1.08*n0*c_s/Lz;
  double n = S0*(sourceFloor + (1.0-sourceFloor)*0.5*(1.0 - tanh((r-r_s)/L_s)));
  fout[0] = n;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double r = sqrt(x*x + y*y);
  double n0 = app->n0;
  double Te_source = app->Te_source;
  double Lperp = app->Lperp;
  double T = Te_source*profileA(r, 1.0/2.5, Lperp);
  fout[0] = T;
}

void
eval_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double T = app->Ti;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *gc = ctx;
  fout[0] = gc->B0;
}

struct gk_lapd_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 3.973*GKYL_PROTON_MASS; // ion mass
  double mRat = 400.0;
  double me = mi/mRat;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 6.0*eV;
  double Ti = 1.0*eV;
  double B0 = 0.0398; // Magnetic field magnitude in Tesla
  double n0 = 2.0e18; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Bias parameters.
  double V_bias = 0;  // Biasing voltage (V).
  double r_bias = 0.25; // Radius of the biased plate (m).

  // Source parameters.
  double Te_source = 6.8*eV;
  double L_s = 0.5*rho_s;
  double r_s = 20*rho_s;

  // Collision parameters.
  double nuFrac = 0.1;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lx = 100.0*rho_s;
  double Ly = 100.0*rho_s;
  double Lperp = Lx;
  double R  = 40.0*rho_s;
  double Lz = 36.0*R;       // approx 18.0492 m.

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 0.75*me*(4.0*vtElc)*(4.0*vtElc)/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 0.75*mi*(4.0*vtIon)*(4.0*vtIon)/(2.0*B0);

  double finalTime = 1.0e-6; 
  double numFrames = 100;

  struct gk_lapd_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .B0 = B0, 
    .n0 = n0, 
    .Te_source = Te_source, 
    .r_s = r_s, 
    .L_s = L_s,
    .Lx = Lx, 
    .Ly = Ly,  
    .Lperp = Lperp, 
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .finalTime = finalTime, 
    .numFrames = numFrames, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr-1);
    gkyl_gyrokinetic_app_calc_mom(app); gkyl_gyrokinetic_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_lapd_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 6); // 18
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 2); // 18
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 2); // 10
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 10); // 10
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 5); // 5

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
    },
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
    },
    
    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },

    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET}, 
                    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET}, 
                    .lo_value = {0.0, 0.0}, .up_value = {0.0, 0.0}}, 
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_bgk_3x2v_p1",

    .cdim = 3, .vdim = 2,
    .lower = { -ctx.Lx/2.0, -ctx.Ly/2.0, -ctx.Lz/2.0 },
    .upper = { ctx.Lx/2.0, ctx.Ly/2.0, ctx.Lz/2.0 },
    .cells = { NX, NY, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // mapping of computational to physical space
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.finalTime;
  double dt = tend-tcurr;
  int nframe = ctx.numFrames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
    }
    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
  gkyl_gyrokinetic_app_stat_write(app);
  
  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  
  return 0;
}
