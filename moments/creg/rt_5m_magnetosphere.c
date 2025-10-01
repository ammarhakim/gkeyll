#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct magnetosphere_ctx
{
  double epsilon0; // permittivity
  double c; // speed of light
  double mu0; // permeability
  double gas_gamma; // 

  double R0; // radius of planet
  double Dx; // dipole field strength (x-direction)
  double Dy; // dipole field strength (y-direction)
  double rc; // radius of inner boundary of field

  double rho_in; // solar wind density
  double vx_in; // solar wind speed (x-direction)
  double vy_in; // solar wind speed (x-direction)
  double vz_in; // solar wind speed (x-direction)
  double p_in; // solar wind pressure
  double Bx_in; // solar wind magnetic field (x-direction)
  double By_in; // solar wind magnetic field (y-direction)
  double Bz_in; // solar wind magnetic field (z-direction)

  double r_ramp1; // scaling factor for velocity near the planet
  double r_ramp2; // scaling factor for velocity near the planet
  double x_mirror; // location of mirror dipole

  double mi; // ion mass
  double mi_me; // ion/electron mass ratio
  double pi_pe; // ion/electron pressure ratio
  double me; // electron mass
  double di; // ion skin depth
  double qi; // ion charge
  double qe; // electron charge
  double de; // electron skin depth

  int Nx; // cell count (x-direction);
  int Ny; // cell count (y-direction);
  double Lx; // domain size (x-direction).
  double Ly; // domain size (y-direction).
  double k0_elc; // closure parameter for electrons.
  double k0_ion; // closure parameter for ions.
  double cfl_frac; // CFL coefficient.

  double t_end; // final simulation time.
  int num_frames; // number of output frames.
  int field_energy_calcs; // number of times to calculate field energy.
  int integrated_mom_calcs; // number of times to calculate integrated moments.
  double dt_failure_tol; // minimum allowable fraction of initial time-step.
  int num_failures_max; // maximum allowable number of consecutive small time-steps.
};

struct magnetosphere_ctx
create_ctx(void)
{
  // physical constants
  double mu0 = GKYL_MU0; // physical permeability (SI)
  double c = GKYL_SPEED_OF_LIGHT/100.0; // reduced speed of light
  double epsilon0 = 1/mu0/(c*c); // reduced permittivity
  double gas_gamma = 5.0/3.0;

  double R0 = 1737.4e3; // radius of the Moon
  double Dx = -100.0e-9*(R0*R0);
  double Dy = -100.0e-9*(R0*R0);
  double rc = 0.5*R0;

  // solar wind inflow parameters
  double rho_in = 10.0e6*GKYL_PROTON_MASS;
  double vx_in = 150.0e3;
  double vy_in = 0.0;
  double vz_in = 0.0;
  double p_in = 0.1e-9;
  double Bx_in = 0.0;
  double By_in = -20.0e-9;
  double Bz_in = 0.0;

  double r_ramp1 = 2.0*R0;
  double r_ramp2 = 3.0*R0;
  double x_mirror = -3.0*R0;

  // electron mass is artificially increased
  double mi = 16.0*GKYL_PROTON_MASS;
  double mi_me = 25.0;
  double pi_pe = 5.0;
  double me = mi/mi_me;  
  
  double di = 0.2*R0;
  double de = di*sqrt(me/mi);

  double qi = mi/di/sqrt(mu0*rho_in);
  double qe = -qi;

  // simulation parameters.
  int Nx = 300; // cell count (x-direction).
  int Ny = 400; // cell count (y-direction).
  double Lx = 30.0*R0; // domain size (x-direction).
  double Ly = 40.0*R0; // domain size (y-direction).
  double k0_elc = 1.0/de; // closure parameter for electrons.
  double k0_ion = 1.0/de; // closure parameter for ions.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 250.0; // final simulation time.
  int num_frames = 40; // number of output frames.
  int field_energy_calcs = INT_MAX; // number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // maximum allowable number of consecutive small time-steps.

  struct magnetosphere_ctx ctx = {
    .mu0 = mu0,
    .c = c,
    .epsilon0 = epsilon0,
    .gas_gamma = gas_gamma,
    .R0 = R0,
    .Dx = Dx,
    .Dy = Dy,
    .rc = rc,
    .rho_in = rho_in,
    .vx_in = vx_in,
    .vy_in = vy_in,
    .vz_in = vz_in,
    .p_in = p_in,
    .Bx_in = Bx_in,
    .By_in = By_in,
    .Bz_in = Bz_in,
    .r_ramp1 = r_ramp1,
    .r_ramp2 = r_ramp2,
    .x_mirror = x_mirror,
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
    .Lx = Lx,
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

  // masks the region inside the planet
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
  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;

  double r_ramp1 = app->r_ramp1;
  double r_ramp2 = app->r_ramp2;

  double rhoe = app->rho_in/(1 + mi_me);
  double vx_in = app->vx_in, vy_in = app->vy_in, vz_in = app->vz_in;
  double me = app->me;
  double R0 = app->R0;

  double pi_pe = app->pi_pe;

  double pe = app->p_in/(1 + pi_pe);

  // reduce solar wind speed near planet
  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;
  double vz = vz_in*s;

  double Ee_tot = pe/(gas_gamma - 1.0) + 0.5*rhoe*(vx*vx + vy*vy + vz*vz);

  // set electron mass density.
  fout[0] = rhoe;
  // set electron momentum density.
  fout[1] = rhoe*vx; fout[2] = rhoe*vy; fout[3] = rhoe*vz;
  // set electron pressure tensor.
  fout[4] = Ee_tot;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;
  double R0 = app->R0;
  double r = sqrt(x*x + y*y);

  double r_ramp1 = app->r_ramp1;
  double r_ramp2 = app->r_ramp2;

  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);
  double vx_in = app->vx_in, vy_in = app->vy_in, vz_in = app->vz_in;
  double mi = app->mi;

  double pi_pe = app->pi_pe;

  double pi = app->p_in - app->p_in/(1 + pi_pe);

  // reduce solar wind speed near planet
  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;
  double vz = vz_in*s;

  double Ei_tot = pi/(gas_gamma - 1.0) + 0.5*rhoi*(vx*vx + vy*vy + vz*vz);

  // set ion mass density.
  fout[0] = rhoi;
  // set ion momentum density.
  fout[1] = rhoi*vx; fout[2] = rhoi*vy; fout[3] = rhoi*vz;
  // set ion total energy density.
  fout[4] = Ei_tot;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double Dx = app->Dx;
  double Dy = app->Dy;
  double x_mirror = app->x_mirror;
  double rc = app->rc;
  double r_ramp1 = app->r_ramp1;
  double r_ramp2 = app->r_ramp2;
  
  double vx_in = app->vx_in;
  double vy_in = app->vy_in;
  double vz_in = app->vz_in;

  double x0 = 2.0*x_mirror;
  
  double r = sqrt(x*x + y*y);
  double rm = sqrt((x - x0)*(x - x0) + y*y);

  // total field without static dipole field
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;

  // total magnetic field
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;

  // dipole field in x > x_mirror
  double Bxt_p = 0.0;
  double Byt_p = 0.0;
  double Bzt_p = 0.0;

  // mirror dipole in x > x_mirror
  double Bxt_m = 0.0;
  double Byt_m = 0.0;
  double Bzt_m = 0.0;

  // solar wind field
  double Bx_in = app->Bx_in;
  double By_in = app->By_in;
  double Bz_in = app->Bz_in;

  // static dipole field
  double Bxt_d = 0.0;
  double Byt_d = 0.0;
  double Bzt_d = 0.0;

  if ((x >= x_mirror) && (x*x + y*y >= rc*rc)) {
    Bxt_p = (2.0*x*Dx*x + 2.0*x*Dy*y - Dx*r*r)/(r*r*r*r);
    Byt_p = (2.0*y*Dx*x + 2.0*y*Dy*y - Dy*r*r)/(r*r*r*r);
  }

  if ((x >= x_mirror) && ((x - x0)*(x - x0) + y*y >= rc*rc)) {
    Bxt_m = (-2.0*(x - x0)*Dx*(x - x0) + 2.0*(x - x0)*Dy*y + Dx*rm*rm)/(rm*rm*rm*rm);
    Byt_m = (-2.0*y*Dx*(x - x0) + 2.0*y*Dy*y - Dy*rm*rm)/(rm*rm*rm*rm);
  }

  if (x*x + y*y >= rc*rc) {
    Bxt_d = (2.0*x*Dx*x + 2.0*x*Dy*y - Dx*r*r)/(r*r*r*r);
    Byt_d = (2.0*y*Dx*x + 2.0*y*Dy*y - Dy*r*r)/(r*r*r*r);
    Bxt = Bxt_p + Bxt_m + Bx_in;
    Byt = Byt_p + Byt_m + By_in;
    Bx = Bxt - Bxt_d;
    By = Byt - Byt_d;
  }

  // reduce solar wind speed near planet
  double s = (r - r_ramp1)/(r_ramp2 - r_ramp1);
  s = fmax(s, 0.0);
  s = fmin(s, 1.0);

  double vx = vx_in*s;
  double vy = vy_in*s;
  double vz = vz_in*s;

  // electric fields induced by total magnetic field and solar wind speed
  double Ex = -vy*Bzt + vz*Byt; // total electric field (x-direction).
  double Ey = -vx*Bzt + vz*Bxt; // total electric field (y-direction).
  double Ez = -vx*Byt + vy*Bxt; // total electric field (z-direction).

  // set electric field.
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExternalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct magnetosphere_ctx *app = ctx;

  double Dx = app->Dx;
  double Dy = app->Dy;
  double rc = app->rc;
  double r = sqrt(x*x + y*y);

  // dipole field
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;

  if ((x*x + y*y) >= (rc*rc)) {
    Bx = (2.0*x*Dx*x + 2.0*x*Dy*y - Dx*r*r)/(r*r*r*r);
    By = (2.0*y*Dx*x + 2.0*y*Dy*y - Dy*r*r)/(r*r*r*r);
  }
  // set external electric field.
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // set external magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
}

// electron solar wind inflow boundary condition
void
evalElcLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;

  double rhoe = app->rho_in/(1 + mi_me);
  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;
  double me = app->me;

  double pi_pe = app->pi_pe;

  double pe = app->p_in/(1 + pi_pe);

  double Ee_tot = pe/(gas_gamma - 1.0) + 0.5*rhoe*(vx*vx + vy*vy + vz*vz);

  // set electron mass density.
  ghost[0] = rhoe;
  // set electron momentum density.
  ghost[1] = rhoe*vx; ghost[2] = rhoe*vy; ghost[3] = rhoe*vz;
  // set electron total energy density.
  ghost[4] = Ee_tot;
}

// ion solar wind inflow boundary condition
void
evalIonLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;

  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);
  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;
  double mi = app->mi;

  double pi_pe = app->pi_pe;

  double pi = app->p_in - app->p_in/(1 + pi_pe);

  double Ei_tot = pi/(gas_gamma - 1.0) + 0.5*rhoi*(vx*vx + vy*vy + vz*vz);

  // set ion mass density.
  ghost[0] = rhoi;
  // set ion momentum density.
  ghost[1] = rhoi*vx; ghost[2] = rhoi*vy; ghost[3] = rhoi*vz;
  // set ion total energy density.
  ghost[4] = Ei_tot;
}

// field solar wind inflow boundary condition
void
evalFieldLowerBC(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double vx = app->vx_in, vy = app->vy_in, vz = app->vz_in;

  double Bx_in = app->Bx_in;
  double By_in = app->By_in;
  double Bz_in = app->Bz_in;

  // electric fields induced by solar magnetic field and speed
  double Ex_in = -vy*Bz_in + vz*By_in; // total electric field (x-direction).
  double Ey_in = -vx*Bz_in + vz*Bx_in; // total electric field (y-direction).
  double Ez_in = -vx*By_in + vy*Bx_in; // total electric field (z-direction).

  // set electric field.
  ghost[0] = Ex_in, ghost[1] = Ey_in; ghost[2] = Ez_in;
  // set magnetic field.
  ghost[3] = Bx_in, ghost[4] = By_in; ghost[5] = Bz_in;
  // set correction potentials.
  ghost[6] = skin[6]; ghost[7] = skin[7];
}

// electron boundary condition at planet surface
static void
evalInnerElc(const double *q, double *qphi, double *delta, void *ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;
  double rhoe = app->rho_in/(1 + mi_me);

  double pi_pe = app->pi_pe;
  double pe = app->p_in/(1 + pi_pe);

  double Ee_tot = pe/(gas_gamma - 1.0);

  // set electron mass density.
  qphi[0] = rhoe;
  // set electron momentum density.
  qphi[1] = 0.0; qphi[2] = 0.0; qphi[3] = 0.0;
  // set electron total energy density.
  qphi[4] = Ee_tot;
}

// ion boundary condition at planet surface
static void
evalInnerIon(const double *q, double *qphi, double *delta, void *ctx)
{
  struct magnetosphere_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  double mi_me = app->mi_me;
  double rhoi = app->rho_in - app->rho_in/(1 + mi_me);

  double pi_pe = app->pi_pe;
  double pi = app->p_in - app->p_in/(1 + pi_pe);

  double Ei_tot = pi/(gas_gamma - 1.0);

  // set ion mass density.
  qphi[0] = rhoi;
  // set ion momentum density.
  qphi[1] = 0.0; qphi[2] = 0.0; qphi[3] = 0.0;
  // set ion total energy density.
  qphi[4] = Ei_tot;
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

  // electron embedded surface
  struct gkyl_wv_embed_geo *embed_geo_elc = gkyl_wv_embed_geo_new(GKYL_EMBED_FUNC,
    evalPhiInit, evalInnerElc, &ctx);

  // ion embedded surface
  struct gkyl_wv_embed_geo *embed_geo_ion = gkyl_wv_embed_geo_new(GKYL_EMBED_FUNC,
    evalPhiInit, evalInnerIon, &ctx); 

  // field embedded surface
  struct gkyl_wv_embed_geo *embed_geo_fld = gkyl_wv_embed_geo_new(GKYL_EMBED_COPY_B,
    evalPhiInit, NULL, &ctx);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma,
    embed_geo_elc, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma,
    embed_geo_ion, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .equation = elc_euler,
    
    .init = evalElcInit,
    .ctx = &ctx,

    .force_low_order_flux = true, // Use Lax fluxes.

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },  
    .bcx_func = { evalElcLowerBC, NULL},
    // bcy defaults to GKYL_SPECIES_COPY
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .equation = ion_euler,
    
    .init = evalIonInit,
    .ctx = &ctx,

    .force_low_order_flux = true, // Use Lax fluxes.

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_COPY },  
    .bcx_func = { evalIonLowerBC, NULL},
    // bcy defaults to GKYL_SPECIES_COPY
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
    // bcy defaults to GKYL_FIELD_COPY
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
    .name = "5m_magnetosphere",

    .ndim = 2,
    .lower = { -ctx.Lx, -ctx.Ly },
    .upper = { ctx.Lx, ctx.Ly }, 
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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  // Free all embedded geometry structures
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
