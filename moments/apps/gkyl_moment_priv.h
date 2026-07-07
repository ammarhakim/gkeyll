// Private header for use in moment app: do not include in user-facing
// header files!
//
// Holds the app struct itself, app-level helpers (IO metadata, integrated
// quantities, wedge/corner BCs), and the time-stepper declarations. The
// per-component objects live in their own private headers, included below:
//   gkyl_moment_species_priv.h   -- struct moment_species  + moment_species_* API
//   gkyl_moment_field_priv.h     -- struct moment_field    + moment_field_* API
//   gkyl_moment_spacetime_priv.h -- struct moment_spacetime + moment_spacetime_* API
//   gkyl_moment_coupling_priv.h  -- struct moment_coupling / mhd_src + their APIs
// They are included here (rather than by each .c file individually) because
// the app struct holds every component by value.
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_elem_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_moment.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>

// The per-component objects and their APIs (see the header comment above).
#include <gkyl_moment_species_priv.h>
#include <gkyl_moment_field_priv.h>
#include <gkyl_moment_spacetime_priv.h>
#include <gkyl_moment_coupling_priv.h>

// number of components that various applied functions should return
enum {
  GKYL_MOM_APP_NUM_APPLIED_CURRENT = 3,
  GKYL_MOM_APP_NUM_EXT_EM = 6,
  GKYL_MOM_APP_NUM_APPLIED_ACCELERATION = 3,
  GKYL_MOM_APP_NUM_NT_SOURCE = 2
};

// Moment app object: used as opaque pointer in user code
struct gkyl_moment_app {
  char name[128]; // name of app
  int ndim; // space dimensions
  double tcurr; // current time
  double cfl; // CFL number

  enum gkyl_moment_scheme scheme_type; // scheme to use
  enum gkyl_wave_split_type split_type; // edge splitting to use
  enum gkyl_mp_recon mp_recon; // reconstruction scheme to use
  // should shock-hybrid scheme be used when using KEP?
  bool use_hybrid_flux_kep;

  bool has_braginskii; // has Braginskii transport
  double coll_fac; // multiplicative collisionality factor for Braginskii

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
  int nghost[3]; // number of ghost-cells in each direction

  int is_dir_skipped[3]; // flags to tell if update in direction are skipped

  struct gkyl_rect_grid grid; // grid
  struct gkyl_range local, local_ext; // local, local-ext ranges
  struct gkyl_range global, global_ext; // global, global-ext ranges

  struct gkyl_rect_decomp *decomp; // decomposition object
  struct gkyl_comm *comm;   // communicator object

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx;   // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers

  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  int has_field; // flag to indicate if we have a field
  struct moment_field field; // field data

  int has_spacetime; // flag to indicate if we have a spacetime component
  struct moment_spacetime spacetime; // spacetime data (modular GR fluids)

  // species data
  int num_species;
  struct moment_species *species; // species data

  // work arrays for use in the KEP and MP scheme: these are stored
  // here so they can be reused
  struct {
    struct gkyl_array *ql, *qr;     // expansions on left/right edge of cell
    struct gkyl_array *amdq, *apdq; // minus/plus fluctuations
  };

  int update_sources; // flag to indicate if sources are to be updated
  struct moment_coupling sources; // sources

  int update_mhd_source;
  struct mhd_src mhd_source;

  struct gkyl_moment_stat stat; // statistics

  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_moment_app *app, double dt0);

  bool has_collision; // has collisions
  // scaling factors for collision frequencies so that nu_sr=nu_base_sr/rho_s
  // nu_rs=nu_base_rs/rho_r, and nu_base_sr=nu_base_rs
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
};

// Meta-data for IO
struct moment_output_meta {
  int frame; // frame number
  double stime; // output time
};

/** Some common functions to species and fields */

// functions for use in integrated quantities calculation
static inline void
integ_unit(int nc, const double *qin, double *integ_out)
{
  for (int i = 0; i < nc; ++i)
    integ_out[i] = qin[i];
}
static inline void
integ_sq(int nc, const double *qin, double *integ_out)
{
  for (int i = 0; i < nc; ++i)
    integ_out[i] = qin[i] * qin[i];
}

// function for copy BC
static inline void
bc_copy(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin,
  double *GKYL_RESTRICT ghost, void *ctx)
{
  for (int c = 0; c < nc; ++c)
    ghost[c] = skin[c];
}

// function for skip BCs
static inline void
bc_skip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double *skin,
  double *GKYL_RESTRICT ghost, void *ctx)
{
}

/**
 * Read the grid/array header of a restart file and validate it against the
 * app's grid. Shared by the per-component restart readers.
 *
 * @param app App object.
 * @param fname File to read the header from.
 * @return Restart status: IO status plus the frame number and simulation
 *   time from the file's embedded metadata.
 */
struct gkyl_app_restart_status moment_app_header_from_file(gkyl_moment_app *app,
  const char *fname);

// Compute integrated quantities specified by i_func
void calc_integ_quant(const struct gkyl_wv_eqn *eqn, double vol,
  const struct gkyl_array *q,
  const struct gkyl_wave_geom *geom,
  struct gkyl_range update_rng, double *integ_q);

// Check array "q" for nans
bool check_for_nans(const struct gkyl_array *q, struct gkyl_range update_rng);

// Apply periodic BCs to corner cells of "f" (ONLY WORKS IN 2D)
void moment_apply_periodic_corner_sync_2d(const gkyl_moment_app *app,
  struct gkyl_array *f);

// Apply wedge-periodic BCs to array "f"
void moment_apply_wedge_bc(const gkyl_moment_app *app, double tcurr,
  const struct gkyl_range *update_rng,
  struct gkyl_array *bc_buffer, int dir,
  const struct gkyl_wv_apply_bc *lo,
  const struct gkyl_wv_apply_bc *up,
  struct gkyl_array *f);

/**
 * Return ghost cell layout for grid.
 *
 * @param app App object.
 * @param nghost On output, ghost-cells used for grid.
 *
 */
void gkyl_moment_app_nghost(gkyl_moment_app *app, int nghost[3]);

/** Top-level app API */

// Take a single time-step using a single-step time-stepper
struct gkyl_update_status moment_update_one_step(gkyl_moment_app *app,
  double dt0);

// Take a single time-step using a SSP-RK3 stepper
struct gkyl_update_status moment_update_ssp_rk3(gkyl_moment_app *app,
  double dt0);

/**
 * Create new array meta header from input struct. Free returned
 * object with moment_array_meta_release.
 *
 * @param meta Meta-data for output.
 * @return New meta object to pass to write method.
 */
struct gkyl_msgpack_data* moment_array_meta_new(struct moment_output_meta meta);

/**
 * Release meta struct
 *
 * @param mt Meta object to free
 */
void moment_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Read meta-data from mpack formated binary input
 *
 * @param mt Mpack encoded meta-data
 * @return Meta-data for simulation
 */
struct moment_output_meta moment_meta_from_mpack(struct gkyl_msgpack_data *mt);
