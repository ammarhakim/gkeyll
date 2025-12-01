#include <gkyl_bc_twistshift.h>
#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_comm_io.h>



/**
 * Project a DG field onto the parallel FEM basis to make it
 * continuous along z (or to solve a Poisson equation in 1x).
 *
 * @param app Application object
 * @param field Field object
 * @param arr_dg Input DG array (to be projected)
 * @param arr_fem Output FEM array (projected result)
 */
void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

/**
 * Initialize the FEM object for the Boltzmann field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 1x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 2x or 3x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Accumulate charge density for a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_boltzmann(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for an adiabatic response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for a Poisson solve.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s, struct gkyl_array **bflux);

/**
 * Solve for the field using a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_boltzmann_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 1x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_1x_poisson_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation for IWL sims.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_deflate_poisson_es_iwl_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 2x and 3x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_2x3x_poisson_perp_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Release resources specific to Boltzmann field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_boltzmann(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources specific to 1x field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release resources specific to 2x3x field solver.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f);


static void
gk_field_enforce_parallel_bc_enabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout)
{
  // Apply the periodicity along the field to fill ghost cells.
  int num_periodic_dir = 1;
  int par_dir = app->cdim - 1;
  int periodic_dirs[] = {par_dir};
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, periodic_dirs, finout); 
  
  // Update the lower z ghosts with twist-and-shift if we are in 3x2v
  if (app->cdim == 3) {
    gkyl_bc_twistshift_advance(field->bc_T_LU_lo, finout, finout);
  }

  // Sync ghost cells between MPI processes.
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, finout);

  // Force the lower skin surface value to match the ghost cell at the node position.
  gkyl_skin_surf_from_ghost_advance(field->ssfg_z_lo, finout);
}

static void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout)
{
  // Do nothing.
}

static void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0/dt, field->es_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
}


static void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  gkyl_deflated_fem_poisson_advance(field->flr_op, phi, phi, phi);
}

static void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
}

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}



static inline void gk_field_polarization_potential_new(struct gk_field *f, struct gkyl_gyrokinetic_app *app)
{
  // Project the initial potential onto a p+1 tensor basis and compute the polarization
  // density to use use by species in calculating the initial ion density.
  f->init_phi_pol = true;
  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order + 1);

  f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
  struct gkyl_array *phi_pol_ho = app->use_gpu ? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                               : gkyl_array_acquire(f->phi_pol);

  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = f->info.polarization_potential,
      .ctx = f->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
  });

  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
  gkyl_array_copy(f->phi_pol, phi_pol_ho);

  gkyl_eval_on_nodes_release(phi_pol_proj);
  gkyl_array_release(phi_pol_ho);
}

static inline void gk_field_polarization_potential_release(struct gk_field *f)
{
  gkyl_array_release(f->phi_pol);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 2 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }, 
        GKYL_GK_META_NONE,
        0
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}


static inline void gk_field_biassed_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->phi_wall_lo = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_lo = false;
  f->phi_wall_lo_evolve = false;
  if (f->info.phi_wall_lo)
  {
    f->has_phi_wall_lo = true;
    if (f->info.phi_wall_lo_evolve)
      f->phi_wall_lo_evolve = f->info.phi_wall_lo_evolve;

    f->phi_wall_lo_host = f->phi_wall_lo;
    if (app->use_gpu)
      f->phi_wall_lo_host = mkarr(false, f->phi_wall_lo->ncomp, f->phi_wall_lo->size);

    f->phi_wall_lo_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
                                                 1, f->info.phi_wall_lo, f->info.phi_wall_lo_ctx);

    // Compute phi_wall_lo at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_lo_proj, 0.0, &app->local_ext, f->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(f->phi_wall_lo, f->phi_wall_lo_host);
  }

  // Set up biased upper wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_up = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_up = false;
  f->phi_wall_up_evolve = false;
  if (f->info.phi_wall_up)
  {
    f->has_phi_wall_up = true;
    if (f->info.phi_wall_up_evolve)
      f->phi_wall_up_evolve = f->info.phi_wall_up_evolve;

    f->phi_wall_up_host = f->phi_wall_up;
    if (app->use_gpu)
      f->phi_wall_up_host = mkarr(false, f->phi_wall_up->ncomp, f->phi_wall_up->size);

    f->phi_wall_up_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
                                                 1, f->info.phi_wall_up, f->info.phi_wall_up_ctx);

    // Compute phi_wall_up at t = 0.
    gkyl_eval_on_nodes_advance(f->phi_wall_up_proj, 0.0, &app->local_ext, f->phi_wall_up_host);
    if (app->use_gpu) // Note: phi_wall_up_host is same as phi_wall_up when not on GPUs.
      gkyl_array_copy(f->phi_wall_up, f->phi_wall_up_host);
  }
}

static inline void gk_field_biassed_wall_release(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->phi_wall_lo);
  if (f->has_phi_wall_lo)
  {
    gkyl_eval_on_nodes_release(f->phi_wall_lo_proj);
    if (app->use_gpu)
      gkyl_array_release(f->phi_wall_lo_host);
  }

  gkyl_array_release(f->phi_wall_up);
  if (f->has_phi_wall_up)
  {
    gkyl_eval_on_nodes_release(f->phi_wall_up_proj);
    if (app->use_gpu)
      gkyl_array_release(f->phi_wall_up_host);
  }
}

static inline void gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  assert(app->cdim > 1);
  f->invert_flr = gk_field_invert_flr;

  double flr_weight = 0.0;
  for (int i = 0; i < app->num_species; ++i)
  {
    struct gk_species *s = &app->species[i];
    double gyroradius_bmag = s->info.flr.bmag ? s->info.flr.bmag : app->bmag_ref;
    flr_weight += s->info.flr.Tperp * s->info.mass / (pow(s->info.charge * gyroradius_bmag, 2.0));
  }
  // Initialize the weight in the Laplacian operator.
  f->flr_rhoSq_sum = mkarr(app->use_gpu, (2 * (app->cdim - 1) - 1) * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gxxj, 0 * app->basis.num_basis);
  if (app->cdim > 2)
  {
    gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gxyj, 1 * app->basis.num_basis);
    gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gyyj, 2 * app->basis.num_basis);
  }
  // Initialize the factor multiplying the field in the FLR operator.
  f->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(f->flr_kSq, -pow(sqrt(2.0), app->cdim), 0); // Sets kSq=-1.

  // If domain is not periodic use Dirichlet BCs.
  struct gkyl_poisson_bc flr_bc = {};

  bool bc_is_np[GKYL_MAX_CDIM]; // Is the BC in this direction non-periodic?
  for (int d = 0; d < app->cdim; ++d)
    bc_is_np[d] = true;
  for (int d = 0; d < app->num_periodic_dir; ++d)
    bc_is_np[app->periodic_dirs[d]] = false;

  for (int d = 0; d < app->cdim - 1; d++)
  {
    if (bc_is_np[d])
    {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0)
        flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0)
        flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
    }
    else
    {
      flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
    }
  }
  // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
  f->flr_op = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
                                            app->local, app->local, f->flr_rhoSq_sum, f->flr_kSq, flr_bc, NULL, app->use_gpu);
}

static inline void gk_field_flr_release(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->flr_rhoSq_sum);
  gkyl_array_release(f->flr_kSq);
  gkyl_deflated_fem_poisson_release(f->flr_op);
}

static inline void gk_field_time_rate_diags_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_active;
  if (app->use_gpu)
  {
    f->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_new, 0, sizeof(double[1]));
    gkyl_cu_memset(f->em_energy_red_old, 0, sizeof(double[1]));
  }
  else
  {
    f->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
    memset(f->em_energy_red_new, 0, sizeof(double[1]));
    memset(f->em_energy_red_old, 0, sizeof(double[1]));
  }
  f->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_dot_write_call = true;
}

static inline void gk_field_time_rate_diags_release(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  if (app->use_gpu)
  {
    gkyl_cu_free(f->em_energy_red_new);
    gkyl_cu_free(f->em_energy_red_old);
  }
  else
  {
    gkyl_free(f->em_energy_red_new);
    gkyl_free(f->em_energy_red_old);
  }
  gkyl_dynvec_release(f->integ_energy_dot);
}


static inline void gk_field_calc_energy_enabled(gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth,
                               1.0, field->es_energy_fac, &app->local, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = {0.0};
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->em_energy_red_global[0];

  if (app->cdim == 1)
    energy_global[0] *= field->es_energy_fac_1d;

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (field->info.time_rate_diagnostics)
  {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_old, field->em_energy_red_global);
    double energy_dot_global_old[1] = {0.0};
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_old[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_old[0] *= field->es_energy_fac_1d;

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_new, field->em_energy_red_global);
    double energy_dot_global_new[1] = {0.0};
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_new[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_new[0] *= field->es_energy_fac_1d;

    double energy_dot_global[1] = {0.0};
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }
}

static inline void gk_field_calc_energy_disabled(gkyl_gyrokinetic_app *app, const struct gk_field *field, double tm)
{
  // Do nothing.
}