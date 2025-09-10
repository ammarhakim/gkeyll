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
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = "tensor",
    }
  );

  // Write out the multiplicative function.
  const char *fmt = "%s-%s_fdot_multiplier_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->fdot_mult.multiplier_host, gks->fdot_mult.multiplier);

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
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  // Multiply out by the multplier.
  gkyl_array_scale_by_cell(out, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_loss_cone_mult(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  // Find the potential at the mirror throat.
  gkyl_dg_basis_ops_eval_array_at_coord_comp(phi, fdmul->bmag_max_coord,
    app->basis_on_dev, &app->grid, &app->local, fdmul->phi_m);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, fdmul->phi_m, fdmul->phi_m_global);

  // Project the loss cone mask.
  gkyl_loss_cone_mask_gyrokinetic_advance(fdmul->lcm_proj_op, &gks->local, &app->local,
    phi, fdmul->phi_m_global, fdmul->multiplier);

  // Multiply out by the multplier.
  gkyl_array_scale_by_cell(out, fdmul->multiplier);
}

void
gk_species_fdot_multiplier_advance_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

void
gk_species_fdot_multiplier_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul)
{
  fdmul->type = gks->info.time_rate_multiplier.type;

  int num_quad = gks->info.time_rate_multiplier.num_quad? gks->info.time_rate_multiplier.num_quad : 1; // Default is p=0.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  fdmul->write_func = gk_species_fdot_multiplier_write_disabled;
  fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_disabled;
  fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_disabled;

  if (fdmul->type) {
    // Allocate rate array.
    fdmul->multiplier = mkarr(app->use_gpu, num_quad==1? 1 : gks->basis.num_basis, gks->local_ext.volume);
    fdmul->multiplier_host = fdmul->multiplier;
    if (app->use_gpu)
      fdmul->multiplier_host = mkarr(false, fdmul->multiplier->ncomp, fdmul->multiplier->size); 

    if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_USER_INPUT) {
      struct gkyl_array *multiplier_high_order_host = num_quad==1? mkarr(false, gks->basis.num_basis, gks->local_ext.volume)
	                                                         : gkyl_array_acquire(fdmul->multiplier_host);

      struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
      proj_on_basis_c2p_ctx.cdim = app->cdim;
      proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
      proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &gks->basis,
          .num_quad = num_quad,
          .num_ret_vals = 1,
          .eval = gks->info.time_rate_multiplier.profile,
          .ctx = gks->info.time_rate_multiplier.profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, multiplier_high_order_host);
      gkyl_proj_on_basis_release(projup);

      if (num_quad == 1) {
        gkyl_array_set_offset(fdmul->multiplier_host, 1.0/pow(sqrt(2.0),gks->grid.ndim), multiplier_high_order_host, 0);
        gkyl_array_copy(fdmul->multiplier, fdmul->multiplier_host);
      }
      else
        gkyl_array_copy(fdmul->multiplier, multiplier_high_order_host);

      gkyl_array_release(multiplier_high_order_host);

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->write_func = gk_species_fdot_multiplier_write_init_only;

    }
    else if (fdmul->type == GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE) {
      // Maximum bmag.
      double *bmag_max_local;
      if (app->use_gpu) {
        bmag_max_local = gkyl_cu_malloc(sizeof(double));
        fdmul->bmag_max = gkyl_cu_malloc(sizeof(double));
      }
      else {
        bmag_max_local = gkyl_malloc(sizeof(double));
        fdmul->bmag_max = gkyl_malloc(sizeof(double));
      }
      gkyl_array_dg_reducec_range(bmag_max_local, app->gk_geom->bmag, 0, GKYL_MAX, app->basis_on_dev, &app->local);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, bmag_max_local, fdmul->bmag_max);
      if (app->use_gpu)
        gkyl_cu_free(bmag_max_local);
      else
        gkyl_free(bmag_max_local);

      // Get the location of maximum bmag.
      double bmag_max_coord_ho[GKYL_MAX_CDIM];
      if (app->cdim == 1) {
        bmag_max_coord_ho[0] = 0.98;
      }
      else if (app->cdim == 2) {
        bmag_max_coord_ho[0] = 1.0001*2.653090e-05;
        bmag_max_coord_ho[1] = 0.98;
      }
      if (app->use_gpu) {
        fdmul->bmag_max_coord = gkyl_cu_malloc(app->cdim*sizeof(double));
	gkyl_cu_memcpy(fdmul->bmag_max_coord, bmag_max_coord_ho, app->cdim*sizeof(double), GKYL_CU_MEMCPY_H2D);
      }
      else {
        fdmul->bmag_max_coord = gkyl_malloc(app->cdim*sizeof(double));
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
        .bmag = app->gk_geom->bmag,
        .bmag_max = fdmul->bmag_max,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .num_quad = num_quad,
        .use_gpu = app->use_gpu,
      };
      fdmul->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew( &inp_proj );

      fdmul->advance_times_cfl_func = gk_species_fdot_multiplier_advance_loss_cone_mult;
      fdmul->advance_times_rate_func = gk_species_fdot_multiplier_advance_mult;
      fdmul->write_func = gk_species_fdot_multiplier_write_enabled;
    }
  }
}

void
gk_species_fdot_multiplier_advance_times_cfl(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_cfl_func(app, gks, fdmul, phi, out);

  app->stat.species_fdot_mult_tm += gkyl_time_diff_now_sec(wst);
}
  
void
gk_species_fdot_multiplier_advance_times_rate(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_fdot_multiplier *fdmul, const struct gkyl_array *phi, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();

  fdmul->advance_times_rate_func(app, gks, fdmul, phi, out);

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
    if (app->use_gpu)
      gkyl_array_release(fdmul->multiplier_host);

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
  }
}
