#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <stdbool.h>

// Lower-level inputs: in general this does not need to be set by the
// user. It is needed when the App is being created on a sub-range of
// the global range, and is meant for use in higher-level drivers that
// use MPI or other parallel mechanism.
struct gkyl_gk_low_inp {
  // local range over which App operates
  struct gkyl_range local_range;
  // communicator to used
  struct gkyl_comm *comm;
};

// Parameters for projection
struct gkyl_gyrokinetic_projection {
  enum gkyl_projection_id proj_id; // type of projection (see gkyl_eqn_type.h)

  // pointer and context to initialization function 
  void *ctx_func; 
  void (*func)(double t, const double *xn, double *fout, void *ctx); 

  // pointers and contexts to initialization functions for gk maxwellian projection
  void *ctx_density;
  void (*density)(double t, const double *xn, double *fout, void *ctx);
  void *ctx_upar;
  void (*upar)(double t, const double *xn, double *fout, void *ctx);
  // if projection is Maxwellian
  void *ctx_temp;
  void (*temp)(double t, const double *xn, double *fout, void *ctx);
  // if projection is bi-Maxwellian
  void *ctx_temppar;
  void (*temppar)(double t, const double *xn, double *fout, void *ctx);
  void *ctx_tempperp;
  void (*tempperp)(double t, const double *xn, double *fout, void *ctx);
};

// Parameters for species collisions
struct gkyl_gyrokinetic_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)

  void *ctx; // context for collision function
  // function for computing self-collision frequency
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);

  // inputs for Spitzer collisionality
  bool normNu; // Set to true if you want to rescale collision frequency
  double self_nu_fac; // Self collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double cross_nu_fac[GKYL_MAX_SPECIES]; // Cross collision freqs without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double bmag_mid; // bmag at the middle of the domain
  double nuFrac; // Parameter for rescaling collision frequency from SI values
  double hbar; // Planck's constant/2 pi 

  int num_cross_collisions; // number of species to cross-collide with
  char collide_with[GKYL_MAX_SPECIES][128]; // names of species to cross collide with
};

// Parameters for species diffusion
struct gkyl_gyrokinetic_diffusion {
  int num_diff_dir; // number of diffusion directions
  int diff_dirs[3]; // list of diffusion directions
  double D[3]; // constant diffusion coefficient in each direction
  int order; // integer for order of the diffusion (4 for grad^4, 6 for grad^6, default is grad^2)
};

// Parameters for species source
struct gkyl_gyrokinetic_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // optional parameter to write out source
  int num_sources;

  // sources using projection routine
  struct gkyl_gyrokinetic_projection projection[GKYL_MAX_SOURCES];
};

// Parameters for boundary conditions
struct gkyl_gyrokinetic_bc {
  enum gkyl_species_bc_type type;
  void *aux_ctx;
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx);  
  double aux_parameter;
};

struct gkyl_gyrokinetic_bcs {
  struct gkyl_gyrokinetic_bc lower, upper;
};

struct gkyl_gyrokinetic_geometry {
  enum gkyl_geometry_id geometry_id;

  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  void *bmag_ctx; // context for bmag function
  // pointer to bmag function
  void (*bmag_func)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_tok_geo_efit_inp *tok_efit_info; // context with RZ data such as efit file for a tokamak
  struct gkyl_tok_geo_grid_inp *tok_grid_info; // context for tokamak geometry with computational domain info

  struct gkyl_mirror_geo_efit_inp *mirror_efit_info; // context with RZ data such as efit file for a mirror
  struct gkyl_mirror_geo_grid_inp *mirror_grid_info; // context for mirror geometry with computational domain info

  double world[3]; // extra computational coordinates for cases with reduced dimensionality
};

// Parameters for species radiation
struct gkyl_gyrokinetic_radiation {
  enum gkyl_radiation_id radiation_id; // type of radiation

  int num_cross_collisions; // number of species to cross-collide with
  char collide_with[GKYL_MAX_SPECIES][128]; // names of species to cross collide with
  // fitting parameters associated with cross-collisions
  double a[GKYL_MAX_SPECIES];
  double alpha[GKYL_MAX_SPECIES];
  double beta[GKYL_MAX_SPECIES];
  double gamma[GKYL_MAX_SPECIES];
  double v0[GKYL_MAX_SPECIES];

  // Atomic z and charge state of species colliding with
  int z[GKYL_MAX_SPECIES];
  int charge_state[GKYL_MAX_SPECIES];
  int num_of_densities[GKYL_MAX_SPECIES]; // Max number of densities to use per charge state
};

struct gkyl_gyrokinetic_react_type {
  enum gkyl_react_id react_id; // what type of reaction (ionization, charge exchange, recombination)
  enum gkyl_react_self_type type_self; // what is the role of species in this reaction
  enum gkyl_ion_type ion_id; // what type of ion is reacting
  char elc_nm[128]; // names of electron species in reaction
  char ion_nm[128]; // name of ion species in reaction
  char donor_nm[128]; // name of donor species in reaction
  char recvr_nm[128]; // name of receiver species in reaction
  int charge_state; // charge state of species in reaction
  double ion_mass; // mass of ion in reaction
  double elc_mass; // mass of electron in reaction
};

// Parameters for species radiation
struct gkyl_gyrokinetic_react {
  int num_react; // number of reactions
  // 3 types of reactions supported currently
  // Ionization, Charge exchange, and Recombination
  // GKYL_MAX_SPECIES number of reactions supported per species (8 different reactions)
  struct gkyl_gyrokinetic_react_type react_type[GKYL_MAX_SPECIES];
};

// Parameters for gk species
struct gkyl_gyrokinetic_species {
  char name[128]; // species name

  enum gkyl_gkmodel_id gkmodel_id;
  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  // initial conditions using projection routine
  struct gkyl_gyrokinetic_projection projection;

  double polarization_density;

  int num_diag_moments; // number of diagnostic moments
  char diag_moments[16][16]; // list of diagnostic moments

  // collisions to include
  struct gkyl_gyrokinetic_collisions collisions;

  // diffusion coupling to include
  struct gkyl_gyrokinetic_diffusion diffusion;

  // source to include
  struct gkyl_gyrokinetic_source source;

  // radiation to include
  struct gkyl_gyrokinetic_radiation radiation;

  // reactions between plasma species to include
  struct gkyl_gyrokinetic_react react;
  // reactions with neutral species to include
  struct gkyl_gyrokinetic_react react_neut;

  // boundary conditions
  struct gkyl_gyrokinetic_bcs bcx, bcy, bcz;
};

// Parameters for neutral species
struct gkyl_gyrokinetic_neut_species {
  char name[128]; // species name

  double mass; // mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  // initial conditions using projection routine
  struct gkyl_gyrokinetic_projection projection;

  int num_diag_moments; // number of diagnostic moments
  char diag_moments[16][16]; // list of diagnostic moments

  // source to include
  struct gkyl_gyrokinetic_source source;

  // reactions with plasma species to include
  struct gkyl_gyrokinetic_react react_neut;

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameter for gk field
struct gkyl_gyrokinetic_field {
  enum gkyl_gkfield_id gkfield_id;
  double bmag_fac; 
  double kperpSq; // kperp^2 parameter for 1D field equations
  double xLCFS; // radial location of the LCFS.

  // parameters for adiabatic electrons simulations
  double electron_mass, electron_charge, electron_temp;

  enum gkyl_fem_parproj_bc_type fem_parbc;
  struct gkyl_poisson_bc poisson_bcs;

  void *phi_wall_lo_ctx; // context for biased wall potential on lower wall
  // pointer to biased wall potential on lower wall function
  void (*phi_wall_lo)(double t, const double *xn, double *phi_wall_lo_out, void *ctx);
  bool phi_wall_lo_evolve; // set to true if biased wall potential on lower wall function is time dependent  

  void *phi_wall_up_ctx; // context for biased wall potential on upper wall
  // pointer to biased wall potential on upper wall function
  void (*phi_wall_up)(double t, const double *xn, double *phi_wall_up_out, void *ctx);
  bool phi_wall_up_evolve; // set to true if biased wall potential on upper wall function is time dependent  
};

// Top-level app parameters
struct gkyl_gk {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  double lower[3], upper[3]; // lower, upper bounds of config-space
  int cells[3]; // config-space cells
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  struct gkyl_gyrokinetic_geometry geometry; // geometry input struct
                                             //
  double cfl_frac; // CFL fraction to use (default 1.0)

  bool use_gpu; // Flag to indicate if solver should use GPUs

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_gyrokinetic_species species[GKYL_MAX_SPECIES]; // species objects

  int num_neut_species; // number of species
  struct gkyl_gyrokinetic_neut_species neut_species[GKYL_MAX_SPECIES]; // species objects
  
  bool skip_field; // Skip field update -> phi = 0 for all time
  struct gkyl_gyrokinetic_field field; // field object

  // this should not be set by typical user-facing code but only by
  // higher-level drivers
  bool has_low_inp; // should one use low-level inputs?
  struct gkyl_gk_low_inp low_inp; // low-level inputs  
};

// Simulation statistics
struct gkyl_gyrokinetic_stat {
  bool use_gpu; // did this sim use GPU?
  
  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
    
  double total_tm; // time for simulation (not including ICs)
  double init_species_tm; // time to initialize all species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species collisionless RHS
  double field_rhs_tm; // time to compute field RHS

  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_bc_tm; // time to compute species BCs
  double field_bc_tm; // time to compute field

  long nspecies_omega_cfl; // number of times CFL-omega all-reduce is called
  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl

  long nfield_omega_cfl; // number of times CFL-omega for field all-reduce is called
  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for field

  long nmom; // calls to moment calculation
  double mom_tm; // time to compute moments

  long ndiag; // calls to diagnostics
  double diag_tm; // time to compute diagnostics

  long nio; // number of calls to IO
  double io_tm; // time to perform IO
};

// Object representing gk app
typedef struct gkyl_gyrokinetic_app gkyl_gyrokinetic_app;

/**
 * Construct a new gk app.
 *
 * @param gk App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New gk app object.
 */
gkyl_gyrokinetic_app* gkyl_gyrokinetic_app_new(struct gkyl_gk *gk);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic(gkyl_gyrokinetic_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_app_apply_ic_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/**
 * Initialize neutral species by projecting initial conditions on basis
 * functions. Neutral species index (sidx) is the same index used to specify
 * the neutral species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_app_apply_ic_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_app_calc_mom(gkyl_gyrokinetic_app *app);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_app_calc_integrated_mom(gkyl_gyrokinetic_app* app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_app_calc_field_energy(gkyl_gyrokinetic_app* app, double tm);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_field(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write source species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_source_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write collisional moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_coll_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write radiation drag coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_rad_drag(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_mom(gkyl_gyrokinetic_app *app, double tm, int frame);

/**
 * Write integrated diagnostic moments for species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_integrated_mom(gkyl_gyrokinetic_app *app);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_field_energy(gkyl_gyrokinetic_app* app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_app_stat_write(gkyl_gyrokinetic_app* app);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_gyrokinetic_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, ...);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_gyrokinetic_update(gkyl_gyrokinetic_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_app_stat(gkyl_gyrokinetic_app* app);

/**
 * Run the RHS for the species update. This is used to compute kernel
 * timers and is not otherwise a useful function for a full
 * simulation.
 *
 * @param app App object.
 * @param update_vol_term Set to 1 to update vol term also, 0 otherwise
 */
void gkyl_gyrokinetic_app_species_ktm_rhs(gkyl_gyrokinetic_app* app, int update_vol_term);

/**
 * Free gk app.
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_app_release(gkyl_gyrokinetic_app* app);
