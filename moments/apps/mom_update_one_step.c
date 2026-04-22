#include <gkyl_moment_priv.h>

// internal function that takes a single time-step using a single-step
// Strang-split scheme
struct gkyl_update_status
moment_update_one_step(gkyl_moment_app* app, double dt0)
{
  int ns = app->num_species, ndim = app->ndim;
  bool have_nans_occured = false;
  
  double dt_suggested = DBL_MAX;
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    FIELD_UPDATE,
    SPECIES_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  double tcurr = app->tcurr, dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        state = FIRST_COUPLING_UPDATE; // next state
          
        // copy old solution in case we need to redo this step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].fdup, app->species[i].f[0]);
        if (app->has_field)
          gkyl_array_copy(app->field.fdup, app->field.f[0]);

        break;
          
      case FIRST_COUPLING_UPDATE:
        state = FIELD_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src1_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_coupling_update(app, &app->sources,
            0, tcurr, dt/2);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }
        if (app->update_mhd_source) {
          struct timespec src1_tm = gkyl_wall_clock();
          mhd_src_update(app, &app->mhd_source, 0, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }

        break;

      case FIELD_UPDATE:
        state = SPECIES_UPDATE; // next state

        if (app->has_field) {
          struct timespec fl_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_field_update(app, &app->field, tcurr, dt);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
            
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.field_tm += gkyl_time_diff_now_sec(fl_tm);
        }
          
        break;

      case SPECIES_UPDATE:
        state = SECOND_COUPLING_UPDATE; // next state

        struct timespec sp_tm = gkyl_wall_clock();
        for (int i=0; i<ns; ++i) {         
          struct gkyl_update_status s =
            moment_species_update(app, &app->species[i], tcurr, dt);

          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
        }
        app->stat.species_tm += gkyl_time_diff_now_sec(sp_tm);
         
        break;

      case SECOND_COUPLING_UPDATE:
        state = POST_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src2_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }
        if (app->update_mhd_source) {
          struct timespec src2_tm = gkyl_wall_clock();
          mhd_src_update(app, &app->mhd_source, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }

        break;

      case POST_UPDATE:
        state = UPDATE_DONE;

        // copy solution in prep for next time-step
        for (int i=0; i<ns; ++i) {
          // check for nans before copying
          if (check_for_nans(app->species[i].f[ndim], app->local))
            have_nans_occured = true;
          else // only copy in case no nans, so old solution can be written out
            gkyl_array_copy(app->species[i].f[0], app->species[i].f[ndim]);
        }

        
        // static int tov_step = 0;
        // // /* Einstein Toolkit 
        // // double tov_R_star = 11.14;
        // // double tov_pos_x  = 0.5 * (app->grid.upper[0] + app->grid.lower[0]);  // auto-compute center
        // // */
        // double tov_R_star = 13.07;
        // double tov_pos_x  = 100.0;
        // tov_step++;

        // static int eq_check = 0;
        // if (!eq_check && tov_step == 1) {
        //   double x_lo = app->grid.lower[0];
        //   double dx = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

        //   double x_pts[3] = {
        //       tov_pos_x,
        //       tov_pos_x + 0.5*tov_R_star,
        //       tov_pos_x + 0.9*tov_R_star
        //   };

        //   printf("\n=== EQUILIBRIUM CHECK after step 1 ===\n");
        //   printf("  tov_pos_x=%.4f tov_R_star=%.4f\n", tov_pos_x, tov_R_star);
        //   printf("%-10s %-12s\n", "x", "mom_x(t=dt)");

        //   int idx[GKYL_MAX_DIM] = {0};
        //   for (int p = 0; p < 3; p++) {
        //       int ci = (int)((x_pts[p] - x_lo) / dx) + app->local.lower[0];
        //       ci = GKYL_MAX2(app->local.lower[0], GKYL_MIN2(app->local.upper[0], ci));
        //       idx[0] = ci;
        //       const double *q = gkyl_array_cfetch(app->species[0].f[0],
        //           gkyl_range_idx(&app->local, idx));
        //       printf("%-10.3f %-12.4e\n", x_pts[p], q[1]);
        //   }
        //   eq_check = 1;
        // }

        // if (tov_step % 10 == 0 && ns > 0) {
        //   struct gkyl_array *fld = app->species[0].f[0];
        //   struct gkyl_range *rng = &app->local;

        //   double x_lo = app->grid.lower[0];
        //   double dx   = (app->grid.upper[0] - x_lo) / app->grid.cells[0];

        //   // total energy and max Etot
        //   double total_E = 0.0;
        //   double Etot_max = 0.0;
        //   double x_max = 0.0;

        //   int idx[GKYL_MAX_DIM] = {0};
        //   for (int ci = rng->lower[0]; ci <= rng->upper[0]; ci++) {
        //     idx[0] = ci;
        //     const double *q = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));
        //     double cell_E = q[0] * dx;
        //     if (cell_E > 0.0) total_E += cell_E;
        //     if (q[0] > Etot_max) {
        //         Etot_max = q[0];
        //         x_max = x_lo + (ci - rng->lower[0] + 0.5) * dx;
        //     }
        //   }

        //   printf("\n[TOV 1D step=%d t=%.4f Nx=%d]\n", tov_step, tcurr, app->grid.cells[0]);
        //   printf("  Total_E=%.6e  Etot_max=%.4e\n", total_E, Etot_max);
          
        //   // sample at fixed physical locations
        //   double x_sample[6] = {
        //     tov_pos_x,
        //     tov_pos_x + 0.5*tov_R_star,
        //     tov_pos_x + 0.9*tov_R_star,
        //     tov_pos_x + tov_R_star,
        //     tov_pos_x + 0.7*tov_R_star,  // mid-surface
        //     tov_pos_x + 1.5*tov_R_star,
        //   };
          
        //   printf("  %-10s %-12s %-12s %-12s %-6s\n", 
        //       "x", "Etot", "mom_x", "lapse", "reg");
        //   for (int p = 0; p < 6; p++) {
        //     int ci = (int)((x_sample[p] - x_lo) / dx) + rng->lower[0];
        //     ci = GKYL_MAX2(rng->lower[0], GKYL_MIN2(rng->upper[0], ci));
        //     idx[0] = ci;
        //     const double *q = gkyl_array_cfetch(fld, gkyl_range_idx(rng, idx));
        //     double r = fabs(x_sample[p] - tov_pos_x);
        //     printf("  %-10.4f %-12.4e %-12.4e %-12.4f %-6s\n",
        //         x_sample[p], q[0], q[1], q[4],
        //         (r < tov_R_star) ? "in" : "out");
        //   }
        // }
      
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.f[ndim]);
          
        break;

      case UPDATE_REDO:
        state = PRE_UPDATE; // start all-over again
          
        // restore solution and retake step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].f[0], app->species[i].fdup);
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.fdup);
          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return (struct gkyl_update_status) {
    .success = have_nans_occured ? false : true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}
