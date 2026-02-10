#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_alloc.h>
#include <gkyl_dg_basis_ops.h>

void
gk_species_fdot_multiplier_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_fdot_multiplier_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for multiplier.
  struct gkyl_msgpack_map_elem mpe_mult[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity" },
  };
  int mpe_mult_len = sizeof(mpe_mult)/sizeof(mpe_mult[0]);
  // Update app basic metada with time/frame.
  gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
  // Package metadata.
  int io_meta_len[] = {app->io_meta_basic_len, mpe_mult_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, mpe_mult, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out the multiplicative function.
  const char *fmt = "%s-%s_fdot_multiplier_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gks->fdot_mult.multiplier_host, gks->fdot_mult.multiplier);
  }

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->fdot_mult.multiplier_host, fileNm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_fdot_multiplier_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_fdot_multiplier_write_enabled(app, gks, tm, frame);
  gks->fdot_mult.write_func = gk_species_fdot_multiplier_write_disabled;
}

void
gk_species_fdot_multiplier_advance_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  // Multiply cflrate by the multplier.
  gkyl_array_scale_by_cell(cflrate, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_loss_cone_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  // Find the potential at the mirror throat.
  gkyl_dg_basis_ops_eval_array_at_coord_comp(phi, fdmul->bmag_max_coord,
    app->basis_on_dev, &app->grid, &app->local, fdmul->phi_m);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, fdmul->phi_m, fdmul->phi_m_global);

  // Project the loss cone mask.
  gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_proj_op, &gks->local, &app->local,
    phi, fdmul->phi_m_global, fdmul->multiplier);

  // Multiply cflrate by the multplier.
  gkyl_array_scale_by_cell(cflrate, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
}

// Set the multiplier array based on cflrate and time dilation parameters.
// See https://arxiv.org/html/2510.09756
void
gk_species_fdot_multiplier_advance_time_dilation_cfl(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  // Compute omega_max - a ceiling on omega is a floor on dt.
  // WARNING: dt_omegaH is DBL_MAX for boltzmann and adiabatic fields!
  double omega_max = DBL_MAX;

  switch (fdmul->type) {
    case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH:
      // Use omega_H based CFL dt flooring.
      omega_max = (gks->dt_omegaH > 1e-30) ? 1.0 / gks->dt_omegaH : DBL_MAX;
      break;

    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD:
      // Use user-specified minimum dt value.
      omega_max = 1.0 / fdmul->cfl_dt_min_value;
      break;

    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL:
    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL: {
      // Use mask-based approach to find omega_max from masked cells.
      gkyl_dg_array_mask_advance(fdmul->cfl_mask, f);
      gkyl_dg_array_mask_scale_by_cell(fdmul->cfl_mask, cflrate);
      const struct gkyl_array *mask_array = gkyl_dg_array_mask_get_mask(fdmul->cfl_mask);

      double omega_max_local;
      if (app->use_gpu) {
      #ifdef GKYL_HAVE_CUDA
        gkyl_array_reduce(fdmul->omega_max_local_cu, mask_array, GKYL_MAX);
        gkyl_cu_memcpy(&omega_max_local, fdmul->omega_max_local_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
      #endif
      }
      else {
        gkyl_array_reduce(&omega_max_local, mask_array, GKYL_MAX);
      }
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &omega_max_local, &omega_max);
      break;
    }

    default:
      // No time dilation needed.
      break;
  }

  // Compute multiplier = min(1.0, omega_max / omega_cfl).
  gkyl_array_copy(fdmul->multiplier, cflrate);
  gkyl_array_invert_by_cell(fdmul->multiplier); // 1/omega_cfl
  gkyl_array_scale(fdmul->multiplier, omega_max); // omega_max / omega_cfl
  gkyl_array_min_by_cell(fdmul->multiplier, 1.0); // min(1.0, omega_max / omega_cfl)

  // Scale the CFL rate by the multiplier.
  gkyl_array_scale_by_cell(cflrate, fdmul->multiplier);
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

void
gk_species_fdot_multiplier_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul)
{
  fdmul->type = gks->info.time_rate_multiplier.type;
  fdmul->write_diagnostics = gks->info.time_rate_multiplier.write_diagnostics;
  fdmul->evolve = gks->info.time_rate_multiplier.evolve;

  // Default function pointers.
  fdmul->write_func = gk_species_fdot_multiplier_write_disabled;
  fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_disabled;
  fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_disabled;

  if (fdmul->type) {
    bool cellwise_const = gks->info.time_rate_multiplier.cellwise_const;
    assert(cellwise_const); // MF 2025/06/11: Limited to this for now.

    // Create a basis for the multiplier.
    struct gkyl_basis basis_mult;
    if (cellwise_const) {
      gkyl_cart_modal_serendip(&basis_mult, gks->basis.ndim, 0);
    } else {
      basis_mult = gks->basis;
    }

    // Allocate multiplier array.
    fdmul->multiplier = mkarr(app->use_gpu, basis_mult.num_basis, gks->local_ext.volume);
    fdmul->multiplier_host = app->use_gpu? mkarr(false, fdmul->multiplier->ncomp, fdmul->multiplier->size)
                                         : gkyl_array_acquire(fdmul->multiplier);

    // Context for c2p function passed to proj_on_basis.
    fdmul->proj_on_basis_c2p_ctx.cdim = app->cdim;
    fdmul->proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
    fdmul->proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
    fdmul->proj_on_basis_c2p_ctx.pos_map = app->position_map;

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {

      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &basis_mult,
          .num_quad = basis_mult.poly_order+1,
          .num_ret_vals = 1,
          .eval = gks->info.time_rate_multiplier.profile,
          .ctx = gks->info.time_rate_multiplier.profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &fdmul->proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, fdmul->multiplier_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(fdmul->multiplier, fdmul->multiplier_host);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics) {
        fdmul->write_func = gk_species_fdot_multiplier_write_init_only;
      } else {
        gkyl_array_release(fdmul->multiplier_host);
      }
    }
    else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
      // Available options:
      //   A) num_quad=1, qtype=GKYL_GAUSS_QUAD. Output: ncomp=1 array.
      //   B) num_quad>1, qtype=GKYL_GAUSS_QUAD or GKYL_GAUSS_LOBATTO_QUAD, cellwise_const=true. Output: ncomp=1 array.
      enum gkyl_quad_type qtype = GKYL_GAUSS_LOBATTO_QUAD;
      int num_quad = gks->basis.poly_order+1; // This can be p+1 or 1. Must be
                                              // at leat p+1 for Gauss-Lobatto.

      // Maximum bmag and its location.
      // NOTE: if the same max bmag occurs at multiple locations,
      // bmag_max_coord may have different values on different MPI processes.
      double bmag_max_coord_ho[GKYL_MAX_CDIM];
      double bmag_max_ho = gkyl_gk_geometry_reduce_arg_bmag(app->gk_geom, GKYL_MAX, bmag_max_coord_ho);
      double bmag_max_local = bmag_max_ho;
      double bmag_max_global;
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);
      double bmag_max_coord_local[app->cdim], bmag_max_coord_global[app->cdim];
      if (fabs(bmag_max_ho - bmag_max_global) < 1e-16) {
        for (int d=0; d<app->cdim; d++) {
          bmag_max_coord_local[d] = bmag_max_coord_ho[d];
        }
      }
      else {
        for (int d=0; d<app->cdim; d++) {
          bmag_max_coord_local[d] = -DBL_MAX;
        }
      }
      gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, app->cdim, bmag_max_coord_local, bmag_max_coord_global);

      if (app->use_gpu) {
        fdmul->bmag_max = gkyl_cu_malloc(sizeof(double));
        fdmul->bmag_max_coord = gkyl_cu_malloc(app->cdim*sizeof(double));
        gkyl_cu_memcpy(fdmul->bmag_max, &bmag_max_global, sizeof(double), GKYL_CU_MEMCPY_H2D);
        gkyl_cu_memcpy(fdmul->bmag_max_coord, bmag_max_coord_ho, app->cdim*sizeof(double), GKYL_CU_MEMCPY_H2D);
      }
      else {
        fdmul->bmag_max = gkyl_malloc(sizeof(double));
        fdmul->bmag_max_coord = gkyl_malloc(app->cdim*sizeof(double));
        memcpy(fdmul->bmag_max, &bmag_max_global, sizeof(double));
        memcpy(fdmul->bmag_max_coord, bmag_max_coord_ho, app->cdim*sizeof(double));
      }

      // Electrostatic potential at bmag_max_coord.
      if (app->use_gpu) {
        fdmul->phi_m = gkyl_cu_malloc(sizeof(double));
        fdmul->phi_m_global = gkyl_cu_malloc(sizeof(double));
      }
      else {
        fdmul->phi_m = gkyl_malloc(sizeof(double));
        fdmul->phi_m_global = gkyl_malloc(sizeof(double));
      }

      // Operator that projects the loss cone mask.
      struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
        .phase_grid = &gks->grid,
        .conf_basis = &app->basis,
        .phase_basis = &gks->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &gks->local_vel, 
        .vel_map = gks->vel_map,
        .bmag = app->gk_geom->geo_int.bmag,
        .bmag_max = fdmul->bmag_max,
        .bmag_max_loc = fdmul->bmag_max_coord,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .qtype = qtype,
        .num_quad = num_quad,
        .cellwise_trap_loss = cellwise_const,
        .c2p_pos_func = proj_on_basis_c2p_position_func,
        .c2p_pos_func_ctx = &fdmul->proj_on_basis_c2p_ctx,
        .use_gpu = app->use_gpu,
      };
      fdmul->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew( &inp_proj );

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_loss_cone_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics) {
        fdmul->write_func = gk_species_fdot_multiplier_write_enabled;
      } else {
        gkyl_array_release(fdmul->multiplier_host);
      }
    }
    else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
      // Copy input parameters to struct.
      fdmul->cfl_dt_min_value = gks->info.time_rate_multiplier.cfl_dt_min_value;
      fdmul->f_threshold = gks->info.time_rate_multiplier.f_threshold;


      // Create mask object if using mask-based time dilation.
      if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
          (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
          (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
        // Allocate GPU scratch space for reduce operation if using GPU.
        if (app->use_gpu) {
        #ifdef GKYL_HAVE_CUDA
          fdmul->omega_max_local_cu = (double *)gkyl_cu_malloc(sizeof(double));
        #endif
        }

        enum gkyl_dg_array_mask_types mask_type = GKYL_DG_ARRAY_MASK_NONE;
        switch (fdmul->type) {
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC_CONF;
            break;
          case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL:
            mask_type = GKYL_DG_ARRAY_MASK_C0_GREATER_FRAC;
            break;
          default:
            break;
        }

        struct gkyl_dg_array_mask_inp cfl_mask_inp = {
          .type = mask_type,
          .threshold = fdmul->f_threshold,
          .phase_rng = &gks->local,
          .phase_rng_ext = &gks->local_ext,
          .conf_rng = &app->local,
          .conf_rng_ext = &app->local_ext,
          .vel_rng = &gks->local_vel,
          .use_gpu = app->use_gpu,
        };
        fdmul->cfl_mask = gkyl_dg_array_mask_new(cfl_mask_inp);
      }

      // Initialize multiplier to 1.0.
      gkyl_array_clear(fdmul->multiplier, 1.0);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_time_dilation_cfl;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      if (fdmul->write_diagnostics) {
        fdmul->write_func = gk_species_fdot_multiplier_write_enabled;
      }
      else {
        gkyl_array_release(fdmul->multiplier_host);
      }
    }
  }
}

// Compute the initial static mask for time dilation (called after species IC and field are set).
// This function computes the CFL rate via a full RHS call (which internally computes the multiplier),
// and then switches to the static advance function for subsequent time steps.
void
gk_species_fdot_multiplier_apply_ic(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *f)
{
  // Only compute initial mask for time dilation types with evolve=false.
  bool is_time_dilation_type = (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
                               (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
                               (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
                               (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
                               (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL);

  if (!fdmul->evolve && is_time_dilation_type) {
    // Compute the full RHS to get an accurate CFL rate from all terms.
    // The RHS internally calls gk_species_fdot_multiplier_advance_times_cfl,
    // which will compute the multiplier since the function pointer is currently
    // set to the evolving version (gk_species_fdot_multiplier_advance_time_dilation_cfl).
    // Use f1 as scratch for rhs output, and bflux.f for boundary flux moments.
    gk_species_rhs(app, gks, f, gks->f1, gks->bflux.f);

    // Switch to the static advance function for subsequent time steps.
    // The multiplier was already computed inside gk_species_rhs.
    fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
  }
}

void
gk_species_fdot_multiplier_advance_times_cfl(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_cfl_func(app, gks, fdmul, phi, f, cflrate);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}
  
void
gk_species_fdot_multiplier_advance_times_rate(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, const struct gkyl_array *f, struct gkyl_array *cflrate)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_rate_func(app, gks, fdmul, phi, f, cflrate);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
  
}

void
gk_species_fdot_multiplier_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->fdot_mult.write_func(app, gks, tm, frame);
}

void
gk_species_fdot_multiplier_release(const struct gkyl_gyrokinetic_app *app, const struct gk_fdot_multiplier *fdmul)
{
  if (fdmul->type) {
    gkyl_array_release(fdmul->multiplier);
    if (fdmul->write_diagnostics) {
      gkyl_array_release(fdmul->multiplier_host);
    }

    if (fdmul->type == GKYL_GK_DAMPING_USER_INPUT) {
      // Nothing to release.
    }
    else if (fdmul->type == GKYL_GK_DAMPING_LOSS_CONE) {
      if (app->use_gpu) {
        gkyl_cu_free(fdmul->bmag_max);
        gkyl_cu_free(fdmul->bmag_max_coord);
        gkyl_cu_free(fdmul->phi_m);
        gkyl_cu_free(fdmul->phi_m_global);
      }
      else {
        gkyl_free(fdmul->bmag_max);
        gkyl_free(fdmul->bmag_max_coord);
        gkyl_free(fdmul->phi_m);
        gkyl_free(fdmul->phi_m_global);
      }
      gkyl_loss_cone_mask_gyrokinetic_release(fdmul->lcm_proj_op);
    }
    else if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
             (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
      if ((fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD) ||
          (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL) ||
          (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL)) {
        gkyl_dg_array_mask_release(fdmul->cfl_mask);
      }
      if (app->use_gpu) {
      #ifdef GKYL_HAVE_CUDA
        gkyl_cu_free(fdmul->omega_max_local_cu);
      #endif
      }
    }
  }
}

void
gk_species_fdot_multiplier_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp)
{
  gk_species_fdot_multiplier_release(app, fdmul);

  gks->info.time_rate_multiplier = fdot_mult_inp;
  gk_species_fdot_multiplier_init(app, gks, fdmul);
}
