#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_2x3x_add_TSBC_and_SSFG_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Parallel direction index (handle 2x and 3x cases).
  int par_dir = app->cdim - 1;
  // Take the TS function from the parallel BC of the first species.
  struct gk_species *gks = &app->species[0];
  const struct gkyl_gyrokinetic_bc *par_lower_bc;
  for (int i = 0; i < 2*app->cdim; i++) {
    if (gks->info.bcs[i].dir == par_dir && gks->info.bcs[i].edge == GKYL_LOWER_EDGE) {
      par_lower_bc = (const struct gkyl_gyrokinetic_bc *) &gks->info.bcs[i];
      break;
    }
  }

  // TSBC updaters
  int ghost[] = {1, 1, 1};
  if (app->cdim == 3) {
    // TS BC updater for up to low TS for the lower edge
    // this sets ghost_L = T_LU(ghost_L)
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = app->local_par_ext_core,
      .num_ghost = ghost, // one ghost per config direction
      .basis = app->basis,
      .grid = app->grid,
      .shift_func = par_lower_bc->aux_profile,
      .shift_func_ctx = par_lower_bc->aux_ctx,
      .use_gpu = app->use_gpu,
    };
    // Add the forward TS updater to f
    f->bc_T_LU_lo = gkyl_bc_twistshift_new(&T_LU_lo);
  }

  // Add the SSFG updater for lower and upper application.
  f->ssfg_z_lo = gkyl_skin_surf_from_ghost_new(par_dir,  GKYL_LOWER_EDGE,
    app->basis, &app->lower_skin_par_core,  &app->lower_ghost_par_core, app->use_gpu);
}

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
gk_field_rhs_poisson_perp_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Smooth the charge density along z.
  gk_field_fem_projection_par(app, field, field->rho_c, field->rho_c);

  // Solve the Poisson equation.
  gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c);
  gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_smooth);

  // Smooth the potential along z.
  gk_field_fem_projection_par(app, field, field->phi_smooth, field->phi_smooth);

  // Finish the Poisson solve with FLR effects.
  field->invert_flr(app, field, field->phi_smooth);

  // Enforce a BC of the field in the parallel direction.
  field->enforce_parallel_bc_func(app, field, field->phi_smooth);
}

static void
gk_field_rhs_deflate_poisson_es_iwl(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Gather charge density into global array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);

  // Smooth the charge density. Input is rho_c_global_dg, globally smoothed in z,
  // and then output should be in *local* phi_smooth.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_core, field->rho_c_global_smooth);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_sol, field->rho_c_global_smooth);

  // Solve the Poisson equation.
  gkyl_deflated_fem_poisson_advance(field->fem_poisson_deflated, field->rho_c_global_smooth,
    field->phi_bc, field->phi_smooth);

  // Finish the Poisson solve with FLR effects.
  field->invert_flr(app, field, field->phi_smooth);

  // Enforce a BC of the field in the parallel direction.
  field->enforce_parallel_bc_func(app, field, field->phi_smooth);
}

static void 
gk_field_ohm_solve(struct gkyl_gyrokinetic_app *app, struct gk_field *field){
  struct timespec wst = gkyl_wall_clock();
  
  // The update of kSq is causing issues with leaks and valgrind.
  // gkyl_fem_poisson_perp_update_kSq(field->fem_apardot_solver, field->dApartdtSlvr_kSq);

  // Temporary solution to avoid superLU leak when updating kSq: we release and re-create the solver
  gkyl_fem_poisson_perp_release(field->fem_apardot_solver);
  field->fem_apardot_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &field->ampere_bcs, field->lapWeightAmpere, field->dApartdtSlvr_kSq, app->use_gpu);

  gkyl_fem_poisson_perp_set_rhs(field->fem_apardot_solver, field->currentDensdot);
  gkyl_fem_poisson_perp_solve(field->fem_apardot_solver, field->apardot);

  // gkyl_array_clear_range(field->apardot, 0.0, &app->local_ext);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field){
  struct timespec wst = gkyl_wall_clock();
  gkyl_fem_poisson_perp_set_rhs(field->fem_apar_solver, field->currentDens);
  gkyl_fem_poisson_perp_solve(field->fem_apar_solver, field->apar);
  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_none(gkyl_gyrokinetic_app *app, struct gk_field *field){
  // Do nothing.
}

static void
gk_field_em_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  gk_field_accumulate_current_dens_dot(app, field, rhs_in);
  gk_field_accumulate_ohms_kSq(app, field, f_in);
  gk_field_ohm_solve(app, field);
}

static void
gk_field_em_rhs_none(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  // Do nothing.
}

static void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);
  gkyl_array_release(f->apar);
  gkyl_array_release(f->apardot);

  if (f->is_em) {
    gkyl_array_release(f->apar_curr);
    gkyl_array_release(f->apar1);
    gkyl_array_release(f->aparnew);
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->currentDens);
    gkyl_array_release(f->currentDensdot);
    gkyl_array_release(f->lapWeightAmpere);
    gkyl_array_release(f->dApartdtSlvr_kSq);
    gkyl_fem_parproj_release(f->fem_apar_parproj);
    gkyl_fem_poisson_perp_release(f->fem_apar_solver);
    gkyl_fem_poisson_perp_release(f->fem_apardot_solver);
    if (app->use_gpu) {
      gkyl_array_release(f->apar_host);
      gkyl_array_release(f->apardot_host);
    }
  }

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_array_release(f->epsilon);

  gkyl_deflated_fem_poisson_release(f->fem_poisson_deflated);
  gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
  if (f->is_dirichletvar) {
    gkyl_array_release(f->phi_bc);
  }
  
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(f->fem_parproj_core);
    gkyl_fem_parproj_release(f->fem_parproj_sol);
  } else {
    gkyl_fem_parproj_release(f->fem_parproj);
  }

  gkyl_array_integrate_release(f->calc_em_energy);

  // Release TS BC and SSFG updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_T_LU_lo);
    }
    gkyl_skin_surf_from_ghost_release(f->ssfg_z_lo);
  }
  
  if (f->use_flr) {
    gk_field_flr_release(app, f);
  }
}

void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Setup electromagnetic variables if needed.
  if (f->is_em) {
    f->apar_curr = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apar1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->aparnew = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    
    f->apar_host = f->apar;
    f->apardot_host = f->apardot;
    if (app->use_gpu) {
      f->apar_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->apardot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    f->currentDens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDensdot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->lapWeightAmpere = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
    f->dApartdtSlvr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_adiabatic;
  } else {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  }

  double polarization_weight = 0.0;
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }
  // Allocate array for the polarization weight times geometric coefficients.
  f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  
  // Initialize the polarization weight.
  struct gkyl_array *Jgij[3] = {app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj};
  for (int i=0; i<app->cdim-2/app->cdim; i++) {
    gkyl_array_set_offset(f->epsilon, polarization_weight, Jgij[i], i*app->basis.num_basis);
  }

  if (f->is_em) {
    // Compute the weights in front of the laplacian operator (metric x 1/mu0)
    for (int i=0; i<app->cdim-2/app->cdim; i++) {
      gkyl_array_set_offset(f->lapWeightAmpere, 1.0/f->info.mu0, Jgij[i], i*app->basis.num_basis);
    }
    gkyl_array_set(f->apar_energy_fac, -0.5, f->lapWeightAmpere);
  }

  // Translate input file BCs into Poisson BCs.
  struct gkyl_poisson_bc poisson_bcs = { };

  bool bc_is_np[GKYL_MAX_CDIM]; // Is the BC in this direction non-periodic?
  for (int d=0; d<app->cdim; ++d) bc_is_np[d] = true;
  for (int d=0; d<app->num_periodic_dir; ++d) {
    bc_is_np[app->periodic_dirs[d]] = false;
  }

  for (int d=0; d<app->cdim-1; d++) {
    if (bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      }
    } else {
      poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
    }
  }

  // Detect if this process contains an edge in the z dimension.
  // for applying bias at the extremal z planes only.
  int ndim = app->grid.ndim;
  poisson_bcs.contains_lower_z_edge = f->global_sub_range.lower[ndim-1] == app->global.lower[ndim-1];
  poisson_bcs.contains_upper_z_edge = f->global_sub_range.upper[ndim-1] == app->global.upper[ndim-1];

  // Initialize the Poisson solver.
  f->fem_poisson_deflated = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
    app->local, f->global_sub_range, f->epsilon, 0, poisson_bcs, f->info.bias_plane_list, app->use_gpu);
  f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &poisson_bcs, f->epsilon, NULL, app->use_gpu);

  f->phi_bc = 0;
  f->is_dirichletvar = false;
  for (int i=0; i<2*app->cdim; i++) {
    f->is_dirichletvar = f->is_dirichletvar ||
                          (f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING ||
                          f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
  }

  if (f->is_dirichletvar) {
    // Project the spatially varying BC if the user specifies it.
    f->phi_bc = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    struct gkyl_array *phi_bc_ho = mkarr(false, f->phi_bc->ncomp, f->phi_bc->size);

    for (int d=0; d<app->cdim; d++) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_lo->aux_profile, bc_lo->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_up->aux_profile, bc_up->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
    }
    gkyl_array_copy(f->phi_bc, phi_bc_ho);
    gkyl_array_release(phi_bc_ho);
  }

  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    f->rhs_phi_func = gk_field_rhs_deflate_poisson_es_iwl;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core, fem_parproj_bc_sol;
    if (app->cdim == 2) {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_PERIODIC;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    } else {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_NONE;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    }

    f->fem_parproj_core = gkyl_fem_parproj_new(&app->global_core, &app->basis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&app->global_sol, &app->basis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
  } else {
    f->rhs_phi_func = gk_field_rhs_poisson_perp_2x3x;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d) {
      if (app->periodic_dirs[d] == app->cdim-1) {
        fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
      }
    }

    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      fem_parproj_bc, 0, 0, app->use_gpu);
  }

  f->ampere_solve = gk_field_ampere_solve_none;
  f->em_rhs_func = gk_field_em_rhs_none;
  // Setup EM solvers.
  if (f->is_em) {
    // Translate input file BCs into Ampere BCs.
    for (int d=0; d<app->cdim-1; d++) {
      if (bc_is_np[d]) {
        struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
        if (bc_lo != 0) {
          f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }

        struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
        if (bc_up != 0) {
          f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      } else {
        f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      }
    }
    // FEM solver for Ampere equation (solved at IC only).
    f->fem_apar_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
        &f->ampere_bcs, f->lapWeightAmpere, NULL, app->use_gpu);
    // FEM solver for Ohm's law (evolves d(Apart)/dt).
    f->fem_apardot_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
      &f->ampere_bcs, f->lapWeightAmpere, f->dApartdtSlvr_kSq, app->use_gpu);
    // FEM smoother for Aparallel.
    f->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_NONE;
    f->fem_apar_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      f->fem_parproj_ampere_bc, 0, 0, app->use_gpu);

    f->ampere_solve = gk_field_ampere_solve_enabled;
    f->em_rhs_func = gk_field_em_rhs_enabled;
  }

  gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);

  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type) {
      f->use_flr = f->use_flr || s->info.flr.type;
    }
  }
  if (f->use_flr) {
    gk_field_flr_new(app, f);
  }

  // Twist-and-shift boundary condition for phi and skin surface from ghost to impose phi periodicity at z=-pi.
  f->enforce_parallel_bc_func = gk_field_enforce_parallel_bc_disabled;
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gk_field_2x3x_add_TSBC_and_SSFG_updaters(app,f);
    f->enforce_parallel_bc_func = gk_field_enforce_parallel_bc_enabled;
  }

  f->solver_release_func = gk_field_fem_release_2x3x;
}