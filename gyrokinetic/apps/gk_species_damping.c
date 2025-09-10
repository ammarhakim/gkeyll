#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_alloc.h>
#include <gkyl_dg_basis_ops.h>

void
gk_species_damping_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_damping_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = 0,
      .basis_type = "tensor",
    }
  );

  // Write out the damping rate.
  const char *fmt = "%s-%s_damping_rate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->damping.rate_host, gks->damping.rate);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->damping.rate_host, fileNm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_damping_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_damping_write_enabled(app, gks, tm, frame);
  gks->damping.write_func = gk_species_damping_write_disabled;
}

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

void
gk_species_damping_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_damping *damp)
{
  damp->type = gks->info.damping.type;
  damp->evolve = false; // Whether the rate is time dependent.

  int num_quad = gks->info.damping.num_quad? gks->info.damping.num_quad : 1; // Default is a p=0 mask.
  assert(num_quad == 1); // MF 2025/06/11: Limited to this for now.

  // Default function pointers.
  damp->write_func = gk_species_damping_write_disabled;

  if (damp->type) {
    // Allocate rate array.
    damp->rate = mkarr(app->use_gpu, num_quad==1? 1 : gks->basis.num_basis, gks->local_ext.volume);
    damp->rate_host = damp->rate;
    if (app->use_gpu)
      damp->rate_host = mkarr(false, damp->rate->ncomp, damp->rate->size); 

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
      proj_on_basis_c2p_ctx.cdim = app->cdim;
      proj_on_basis_c2p_ctx.vdim = gks->local_vel.ndim;
      proj_on_basis_c2p_ctx.vel_map = gks->vel_map;
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &gks->basis,
          .num_quad = num_quad,
          .num_ret_vals = 1,
          .eval = gks->info.damping.rate_profile,
          .ctx = gks->info.damping.rate_profile_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj_on_basis_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, damp->rate_host);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(damp->rate, damp->rate_host);

      if (num_quad == 1)
        gkyl_array_scale_range(damp->rate, 1.0/pow(sqrt(2.0),gks->grid.ndim), &gks->local);
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      damp->evolve = true; // Since the loss cone boundary is proportional to phi(t).

      // Maximum bmag.
      double *bmag_max_local;
      if (app->use_gpu) {
        bmag_max_local = gkyl_cu_malloc(sizeof(double));
        damp->bmag_max = gkyl_cu_malloc(sizeof(double));
      }
      else {
        bmag_max_local = gkyl_malloc(sizeof(double));
        damp->bmag_max = gkyl_malloc(sizeof(double));
      }
      gkyl_array_dg_reducec_range(bmag_max_local, app->gk_geom->bmag, 0, GKYL_MAX, app->basis_on_dev, &app->local);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, bmag_max_local, damp->bmag_max);
      if (app->use_gpu)
        gkyl_cu_free(bmag_max_local);
      else
        gkyl_free(bmag_max_local);

      // Get the location of maximum bmag.
      double bmag_max_coord_ho[] = {0.98};
      if (app->use_gpu) {
        damp->bmag_max_coord = gkyl_cu_malloc(app->cdim*sizeof(double));
	gkyl_cu_memcpy(damp->bmag_max_coord, bmag_max_coord_ho, app->cdim*sizeof(double), GKYL_CU_MEMCPY_H2D);
      }
      else {
        damp->bmag_max_coord = gkyl_malloc(app->cdim*sizeof(double));
	memcpy(damp->bmag_max_coord, bmag_max_coord_ho, app->cdim*sizeof(double));
      }

      // Electrostatic potential at bmag_max_coord.
      if (app->use_gpu) {
        damp->phi_m = gkyl_cu_malloc(sizeof(double));
        damp->phi_m_global = gkyl_cu_malloc(sizeof(double));
      }
      else {
        damp->phi_m = gkyl_malloc(sizeof(double));
        damp->phi_m_global = gkyl_malloc(sizeof(double));
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
        .bmag_max = damp->bmag_max,
        .mass = gks->info.mass,
        .charge = gks->info.charge,
        .num_quad = num_quad,
        .use_gpu = app->use_gpu,
      };
      damp->lcm_proj_op = gkyl_loss_cone_mask_gyrokinetic_inew( &inp_proj );

      // Project the conf-space rate profile provided.
      struct gkyl_array *scale_prof_high_order = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
      struct gkyl_array *scale_prof_high_order_ho = app->use_gpu? mkarr(false, scale_prof_high_order->ncomp, scale_prof_high_order->size)
                                                     : gkyl_array_acquire(scale_prof_high_order);
      
      gkyl_proj_on_basis *projup = gkyl_proj_on_basis_new(&gks->grid, &gks->basis, num_quad, 1, 
        gks->info.damping.rate_profile, gks->info.damping.rate_profile_ctx);
      gkyl_proj_on_basis_advance(projup, 0.0, &gks->local, scale_prof_high_order_ho);
      gkyl_proj_on_basis_release(projup);
      gkyl_array_copy(scale_prof_high_order, scale_prof_high_order_ho);

      damp->scale_prof = mkarr(app->use_gpu, num_quad == 1? 1 : gks->basis.num_basis, gks->local_ext.volume);
      gkyl_array_set_offset(damp->scale_prof, pow(sqrt(2.0),gks->grid.ndim), scale_prof_high_order, 0);

      gkyl_array_release(scale_prof_high_order_ho);
      gkyl_array_release(scale_prof_high_order);

      // Compute the initial damping rate (assuming phi=0 because phi hasn't been computed).
      // Find the potential at the mirror throat.
      gkyl_dg_basis_ops_eval_array_at_coord_comp(app->field->phi_smooth, damp->bmag_max_coord,
        app->basis_on_dev, &app->grid, &app->local, damp->phi_m);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, damp->phi_m, damp->phi_m_global);
      // Project the loss cone mask.
      gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
        app->field->phi_smooth, damp->phi_m_global, damp->rate);
      // Multiply by the user's scaling profile.
      gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);
    }

    // Set function pointers chosen at runtime.
    if (damp->evolve) {
      damp->write_func = gk_species_damping_write_enabled;
    }
    else {
      damp->write_func = gk_species_damping_write_init_only;
    }
  }
}

void
gk_species_damping_advance(gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_damping *damp, 
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *f_buffer,
  struct gkyl_array *rhs, struct gkyl_array *cflrate)
{
  if (damp->type) {
    struct timespec wst = gkyl_wall_clock();
    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      gkyl_array_set(f_buffer, 1.0, fin);
      gkyl_array_scale_by_cell(f_buffer, damp->rate);
      gkyl_array_accumulate(rhs, -1.0, f_buffer);
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      // Find the potential at the mirror throat.
      gkyl_dg_basis_ops_eval_array_at_coord_comp(phi, damp->bmag_max_coord,
        app->basis_on_dev, &app->grid, &app->local, damp->phi_m);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, damp->phi_m, damp->phi_m_global);

      // Project the loss cone mask.
      gkyl_loss_cone_mask_gyrokinetic_advance(damp->lcm_proj_op, &gks->local, &app->local,
        phi, damp->phi_m_global, damp->rate);

      // Assemble the damping term -scale_prof * mask * f.
      gkyl_array_set(f_buffer, 1.0, fin);
      gkyl_array_scale_by_cell(damp->rate, damp->scale_prof);
      gkyl_array_scale_by_cell(f_buffer, damp->rate);
      gkyl_array_accumulate(rhs, -1.0, f_buffer);

    }

    // Add the frequency to the CFL frequency.
    gkyl_array_accumulate(cflrate, 1.0, damp->rate);

    app->stat.species_damp_tm += gkyl_time_diff_now_sec(wst);
  }
}

void
gk_species_damping_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->damping.write_func(app, gks, tm, frame);
}

void
gk_species_damping_release(const struct gkyl_gyrokinetic_app *app, const struct gk_damping *damp)
{
  if (damp->type) {
    gkyl_array_release(damp->rate);
    if (app->use_gpu)
      gkyl_array_release(damp->rate_host);

    if (damp->type == GKYL_GK_DAMPING_USER_INPUT) {
      // Nothing to release.
    }
    else if (damp->type == GKYL_GK_DAMPING_LOSS_CONE) {
      if (app->use_gpu) {
        gkyl_cu_free(damp->bmag_max);
        gkyl_cu_free(damp->bmag_max_coord);
        gkyl_cu_free(damp->phi_m);
        gkyl_cu_free(damp->phi_m_global);
      }
      else {
        gkyl_free(damp->bmag_max);
        gkyl_free(damp->bmag_max_coord);
        gkyl_free(damp->phi_m);
        gkyl_free(damp->phi_m_global);
      }
      gkyl_loss_cone_mask_gyrokinetic_release(damp->lcm_proj_op);
      gkyl_array_release(damp->scale_prof);
    }
  }
}
