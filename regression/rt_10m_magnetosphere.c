#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct magnetosphere_ctx
{
  double epsilon0;
  double c;
  double mu0;

  double R0;
  double D;
  double rc;

  double rho_in;
  double vx_in;
  double vy_in;
  double vz_in;
  double p_in;
  double Bx_in;
  double By_in;
  double Bz_in;

  double mirdip_r_ramp1;
  double mirdip_r_ramp2;
  double mirdip_x_mirror;

  double mi;
  double mi_me;
  double pi_pe;
  double me;
  double di;
  double qi;
  double qe;
  double de;

  int Nx; // Cell count (x-direction);
  int Ny; // Cell count (y-direction);
  double Lx_l; // Domain size (x-direction).
  double Lx_u; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double k0_elc; // Closure parameter for electrons.
  double k0_ion; // Closure parameter for ions.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct magnetosphere_ctx
create_ctx(void)
{
  double mu0 = GKYL_MU0;
  double c = 3.0e6;
  double epsilon0 = 1/mu0/(c*c);

  double R0 = 1737.4e3;
  double D = -100.0e-9*(R0*R0);
  double rc = 0.5*R0;

  double rho_in = 10.0e6*GKYL_PROTON_MASS;
  double vx_in = 150.0e3;
  double vy_in = 0.0;
  double vz_in = 0.0;
  double p_in = 0.1e-9;
  double Bx_in = 0.0;
  double By_in = -20.0e-9;
  double Bz_in = 0.0;

  double mirdip_r_ramp1 = 2.0*R0;
  double mirdip_r_ramp2 = 3.0*R0;
  double mirdip_x_mirror = -3.0*R0;
  
  double mi = GKYL_PROTON_MASS;
  double mi_me = 25.0;
  double pi_pe = 5.0;
  double me = mi/mi_me;  
  double qi = GKYL_ELEMENTARY_CHARGE;
  double qe = -qi;

  double di = 0.2*R0;
  double de = di*sqrt(me/mi);

  // Simulation parameters.
  int Nx = 300; // Cell count (x-direction).
  int Ny = 400; // Cell count (y-direction).
  double Lx_l = -30.0*R0;
  double Lx_u = 30.0*R0;
  double Ly = 40.0*R0; // Domain size (y-direction).
  double k0_elc = 1.0/de; // Closure parameter for electrons.
  double k0_ion = 1.0/di; // Closure parameter for ions.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 250.0; // Final simulation time.
  int num_frames = 40; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct magnetosphere_ctx ctx = {
    .mu0 = mu0,
    .c = c,
    .epsilon0 = epsilon0,
    .R0 = R0,
    .D = D,
    .rc = rc,
    .rho_in = rho_in,
    .vx_in = vx_in,
    .vy_in = vy_in,
    .vz_in = vz_in,
    .p_in = p_in,
    .Bx_in = Bx_in,
    .By_in = By_in,
    .Bz_in = Bz_in,
    .mirdip_r_ramp1 = mirdip_r_ramp1,
    .mirdip_r_ramp2 = mirdip_r_ramp2,
    .mirdip_x_mirror = mirdip_x_mirror,
    .mi = mi,
    .mi_me = mi_me,
    .pi_pe = pi_pe,
    .me = me,  
    .qi = qi,
    .qe = qe,
    .di = di,
    .de = de,
    .Nx = Nx,
    .Ny = Ny,
    .Lx_l = Lx_l,
    .Lx_u = Lx_u,
    .Ly = Ly,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalPhiInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT phi, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double R0 = app->R0;
  double xc = 0.0;
  double yc = 0.0;

  if (((x-xc)*(x-xc) + (y-yc)*(y-yc)) < R0*R0)
    phi[0] = -1.0;
  else
    phi[0] = 1.0;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double r = sqrt(x*x + y*y);
  double mi_me = app->mi_me;

  double r_ramp1 = app->mirdip_r_ramp1;
  double r_ramp2 = app->mirdip_r_ramp2;

  double rhoe = app->rho_in/(1 + mi_me);
  double vx_in = app->vx_in, vy_in = app->vy_in, vz_in = app->vz_in;
  double me = app->me;
  double R0 = app->R0;

  double pi_pe = app->pi_pe;

  double pe = app->p_in/(1 + pi_pe);

  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;
  double vz = vz_in*s;

  /* if ((x*x + y*y) < (R0*R0)) { */
  /*   rhoe = 1.0e-10*rhoe; */
  /*   vx = 0.0; */
  /*   vy = 0.0; */
  /*   vz = 0.0; */
  /*   pe = 1.0e-10*pe; */
  /* } */

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = rhoe*vx; fout[2] = rhoe*vy; fout[3] = rhoe*vz;
  // Set electron pressure tensor.
  fout[4] = pe + rhoe*vx*vx; fout[5] = rhoe*vx*vy; fout[6] = rhoe*vx*vz;
  fout[7] = pe + rhoe*vy*vy; fout[8] = rhoe*vy*vz; fout[9] = pe + rhoe*vz*vz;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double mi_me = app->mi_me;
  double R0 = app->R0;
  double r = sqrt(x*x + y*y);

  double r_ramp1 = app->mirdip_r_ramp1;
  double r_ramp2 = app->mirdip_r_ramp2;

  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);
  double vx_in = app->vx_in, vy_in = app->vy_in, vz_in = app->vz_in;
  double mi = app->mi;

  double pi_pe = app->pi_pe;

  double pi = app->p_in - app->p_in/(1 + pi_pe);

  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;
  double vz = vz_in*s;

  /* if ((x*x + y*y) < (R0*R0)) { */
  /*   rhoi = 1.0e-10*rhoi; */
  /*   vx = 0.0; */
  /*   vy = 0.0; */
  /*   vz = 0.0; */
  /*   pi = 1.0e-10*pi; */
  /* } */

  // Set electron mass density.
  fout[0] = rhoi;
  // Set electron momentum density.
  fout[1] = rhoi*vx; fout[2] = rhoi*vy; fout[3] = rhoi*vz;
  // Set electron pressure tensor.
  fout[4] = pi + rhoi*vx*vx; fout[5] = rhoi*vx*vy; fout[6] = rhoi*vx*vz;
  fout[7] = pi + rhoi*vy*vy; fout[8] = rhoi*vy*vz; fout[9] = pi + rhoi*vz*vz;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double D = app->D;
  double x_mirror = app->mirdip_x_mirror;
  double rc = app->rc;
  double r_ramp1 = app->mirdip_r_ramp1;
  double r_ramp2 = app->mirdip_r_ramp2;

  double Bx_in = app->Bx_in;
  double By_in = app->By_in;
  
  double vx_in = app->vx_in;
  double vy_in = app->vy_in;

  double x0 = 2.0*x_mirror;
  
  double r = sqrt(x*x + y*y);
  double rm = sqrt((x - x0)*(x - x0) + y*y);

  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;

  double Bxt = 0.0;
  double Byt = 0.0;
  double Bxt_m = 0.0;
  double Byt_m = 0.0;

  double Bxt_d = 0.0;
  double Byt_d = 0.0;
  if ((x >= x_mirror) && (x*x + y*y >= rc*rc)) {
    Bxt = (2.0*x*D*x + 2.0*x*D*y - D*r*r)/(r*r*r*r);
    Byt = (2.0*y*D*x + 2.0*y*D*y - D*r*r)/(r*r*r*r);
  }

  if ((x >= x_mirror) && ((x - x0)*(x - x0) + y*y >= rc*rc)) {
    Bxt_m = (-2.0*(x - x0)*D*(x - x0) + 2.0*(x - x0)*D*y + D*rm*rm)/(rm*rm*rm*rm);
    Byt_m = (-2.0*y*D*(x - x0) + 2.0*y*D*y - D*rm*rm)/(rm*rm*rm*rm);
  }

  if (x*x + y*y >= rc*rc) {
    Bxt_d = (2.0*x*D*x + 2.0*x*D*y - D*r*r)/(r*r*r*r);
    Byt_d = (2.0*y*D*x + 2.0*y*D*y - D*r*r)/(r*r*r*r);
    Bx = Bxt + Bxt_m + Bx_in - Bxt_d;
    By = Byt + Byt_m + By_in - Byt_d;
  }

  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = -vx*By_in + vy*Bx_in; // Total electric field (z-direction).

  // Set electric field.
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExternalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double D = app->D;
  double rc = app->rc;
  double r = sqrt(x*x + y*y);

  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0; // External magnetic field (z-direction).

  if ((x*x + y*y) >= (rc*rc)) {
    Bx = (2.0*x*D*x + 2.0*x*D*y - D*r*r)/(r*r*r*r);
    By = (2.0*y*D*x + 2.0*y*D*y - D*r*r)/(r*r*r*r);
  }
  // Set external electric field.
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // Set external magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
}

void
evalElcLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double mi_me = app->mi_me;

  double rhoe = app->rho_in/(1 + mi_me);
  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;
  double me = app->me;

  double pi_pe = app->pi_pe;

  double pe = app->p_in/(1 + pi_pe);

  // Set electron mass density.
  ghost[0] = rhoe;
  // Set electron momentum density.
  ghost[1] = rhoe*vx; ghost[2] = rhoe*vy; ghost[3] = rhoe*vz;
  // Set electron pressure tensor.
  ghost[4] = pe + rhoe*vx*vx; ghost[5] = rhoe*vx*vy; ghost[6] = rhoe*vx*vz;
  ghost[7] = pe + rhoe*vy*vy; ghost[8] = rhoe*vy*vz; ghost[9] = pe + rhoe*vz*vz;
}

void
evalIonLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double mi_me = app->mi_me;

  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);
  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;
  double mi = app->mi;

  double pi_pe = app->pi_pe;

  double pi = app->p_in - app->p_in/(1 + pi_pe);

  // Set electron mass density.
  ghost[0] = rhoi;
  // Set electron momentum density.
  ghost[1] = rhoi*vx; ghost[2] = rhoi*vy; ghost[3] = rhoi*vz;
  // Set electron pressure tensor.
  ghost[4] = pi + rhoi*vx*vx; ghost[5] = rhoi*vx*vy; ghost[6] = rhoi*vx*vz;
  ghost[7] = pi + rhoi*vy*vy; ghost[8] = rhoi*vy*vz; ghost[9] = pi + rhoi*vz*vz;
}

void
evalFieldLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;

  double Bx_in = app->Bx_in;
  double By_in = app->By_in;

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = -vx*By_in + vy*Bx_in; // Total electric field (z-direction).

  // Set electric field.
  ghost[0] = Ex, ghost[1] = Ey; ghost[2] = Ez;
  // Set magnetic field.
  ghost[3] = Bx_in, ghost[4] = By_in; ghost[5] = 0.0;
  // Set correction potentials.
  ghost[6] = skin[6]; ghost[7] = skin[7];
}

static void
evalInnerElc(const double *q, double *qphi, double *delta, void *ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double mi_me = app->mi_me;
  double rhoe = app->rho_in/(1 + mi_me);

  double pi_pe = app->pi_pe;
  double pe = app->p_in/(1 + pi_pe);

  qphi[0] = rhoe;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = 0.0;
  qphi[4] = pe;
  qphi[5] = 0.0;
  qphi[6] = 0.0;
  qphi[7] = pe;
  qphi[8] = 0.0;
  qphi[9] = pe;
}

static void
evalInnerIon(const double *q, double *qphi, double *delta, void *ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double mi_me = app->mi_me;
  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);

  double pi_pe = app->pi_pe;
  double pi = app->p_in - app->p_in/(1 + pi_pe);

  qphi[0] = rhoi;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = 0.0;
  qphi[4] = pi;
  qphi[5] = 0.0;
  qphi[6] = 0.0;
  qphi[7] = pi;
  qphi[8] = 0.0;
  qphi[9] = pi;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
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

  struct magnetosphere_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  struct gkyl_wv_embed_geo *embed_geo_elc = gkyl_wv_embed_geo_new(GKYL_EMBED_FUNC,
    evalPhiInit, evalInnerElc, &ctx);

  struct gkyl_wv_embed_geo *embed_geo_ion = gkyl_wv_embed_geo_new(GKYL_EMBED_FUNC,
    evalPhiInit, evalInnerIon, &ctx); 

  struct gkyl_wv_embed_geo *embed_geo_fld = gkyl_wv_embed_geo_new(GKYL_EMBED_COPY_B,
    evalPhiInit, NULL, &ctx);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_elc,
    false, embed_geo_elc, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_ion,
    false, embed_geo_ion, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .equation = elc_ten_moment,
    
    .init = evalElcInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },  
    .bcx_func = { evalElcLowerBC, NULL},
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .equation = ion_ten_moment,
    
    .init = evalIonInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },  
    .bcx_func = { evalIonLowerBC, NULL},
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .init = evalFieldInit,
    .ctx = &ctx,

    .embed_geo = embed_geo_fld,
    
    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,
    .ext_em_evolve = false,

    .bcx = { GKYL_FIELD_FUNC, GKYL_FIELD_COPY },  
    .bcx_func = { evalFieldLowerBC, NULL },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "10m_magnetosphere",

    .ndim = 2,
    .lower = { ctx.Lx_l, -ctx.Ly },
    .upper = { ctx.Lx_u, ctx.Ly }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_wv_embed_geo_release(embed_geo_elc);
  gkyl_wv_embed_geo_release(embed_geo_ion);
  gkyl_wv_embed_geo_release(embed_geo_fld);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
