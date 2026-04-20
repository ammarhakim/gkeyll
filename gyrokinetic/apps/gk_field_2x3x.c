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
gk_field_2x3x_write_twistshift(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Write the discretized shift (for TS BCs) to file.
  int comm_rank, comm_size;
  gkyl_comm_get_rank(app->comm, &comm_rank);
  gkyl_comm_get_size(app->comm, &comm_size);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};
  const char *fmt = "%s-bc_%s%s_twistshift.gkyl";

  struct gk_species *gks = &app->species[0];
  for (int i = 0; i < 2*app->cdim; i++) {
    if (gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_IWL ||
        gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {

      int dir = gks->info.bcs[i].dir;
      int edi = gks->info.bcs[i].edge;
      if (comm_rank == 0 && edi == GKYL_LOWER_EDGE) {
        struct gkyl_bc_twistshift *bc_ts = f->bc_ts_lo;
        
        struct gkyl_rect_grid shear_grid;
        struct gkyl_range shear_r;
        struct gkyl_basis shift_b;
        struct gkyl_array *shift_dg = gkyl_bc_twistshift_get_shift_objects(bc_ts, &shear_grid, &shear_r, &shift_b);

        // Twistshift updater stores the shift on a restricted range (the core) but a full
        // grid. Create a restricted grid for I/O.
        struct gkyl_rect_grid shear_grid_core;
        double lower[1], upper[1];
        int cells[] = {shear_r.volume};
        if (app->gk_geom->geqdsk_sign_convention == 0) {
          // x increases towards SOL.
          lower[0] = shear_grid.lower[0];
          upper[0] = shear_grid.lower[0] + shear_grid.dx[0]*cells[0];
          gkyl_rect_grid_init(&shear_grid_core, shear_grid.ndim, lower, upper, cells);
        }
        else {
          // x increases towards SOL.
          lower[0] = shear_grid.upper[0] - shear_grid.dx[0]*cells[0];
          upper[0] = shear_grid.upper[0];
          gkyl_rect_grid_init(&shear_grid_core, shear_grid.ndim, lower, upper, cells);
        }

        // Package metadata for shift file.
        struct gkyl_msgpack_map_elem io_meta_shift_dg[] = {
          { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = shift_b.poly_order },
          { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = shift_b.id }
        };
        int io_meta_shift_dg_len = sizeof(io_meta_shift_dg)/sizeof(io_meta_shift_dg[0]);
        int io_meta_shift_len[] = {app->io_meta_basic_len, io_meta_shift_dg_len};
        const struct gkyl_msgpack_map_elem* io_meta_shift[] = {app->io_meta_basic, io_meta_shift_dg};
        struct gkyl_msgpack_data *mt_shift = gkyl_msgpack_create_union(sizeof(io_meta_shift_len)/sizeof(int),
          io_meta_shift_len, io_meta_shift);

        int sz = gkyl_calc_strlen(fmt, app->name, vars[dir], edge[edi]);
        char fileNm[sz+1]; // ensures no buffer overflow
        sprintf(fileNm, fmt, app->name, vars[dir], edge[edi]);

        gkyl_grid_sub_array_write(&shear_grid_core, &shear_r, mt_shift, shift_dg, fileNm);

        gkyl_array_release(shift_dg);
        gkyl_msgpack_data_release(mt_shift);
      }
    }
  }
}

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
    struct gkyl_bc_twistshift_inp ts_lo = {
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
    f->bc_ts_lo = gkyl_bc_twistshift_new(&ts_lo);
  }

  // Add the SSFG updater for lower and upper application.
  f->ssfg_z_lo = gkyl_skin_surf_from_ghost_new(par_dir,  GKYL_LOWER_EDGE,
    app->basis, &app->lower_skin_par_core,  &app->lower_ghost_par_core, app->use_gpu);

  // Write the discrete shift to file.
  if (app->cdim == 3)
    gk_field_2x3x_write_twistshift(app, f);
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
    gkyl_bc_twistshift_advance(field->bc_ts_lo, finout, finout);
  }

  // Sync ghost cells between MPI processes.
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, finout);

  // Force the lower skin surface value to match the ghost cell at the node position.
  gkyl_skin_surf_from_ghost_advance(field->ssfg_z_lo, finout);
}

static void
gk_field_rhs_poisson_perp_2x3x(struct gkyl_gyrokinetic_app *app,
  struct gk_field *gkf, struct gkyl_array *phi_out)
{
  // Smooth the charge density along z.
  gk_field_fem_projection_par(app, gkf, gkf->rho_c, gkf->rho_c);

  // Solve the Poisson equation.
  gkyl_fem_poisson_perp_set_rhs(gkf->fem_poisson_perp, gkf->rho_c);
  gkyl_fem_poisson_perp_solve(gkf->fem_poisson_perp, phi_out);

  // Smooth the potential along z.
  gk_field_fem_projection_par(app, gkf, phi_out, phi_out);

  // Finish the Poisson solve with FLR effects.
  gkf->invert_flr(app, gkf, phi_out);

  // Enforce a BC of the field in the parallel direction.
  gkf->enforce_parallel_bc_func(app, gkf, phi_out);
}

static void
gk_field_rhs_deflate_poisson_es_iwl(struct gkyl_gyrokinetic_app *app,
  struct gk_field *gkf, struct gkyl_array *phi_out)
{
  // Gather charge density into global array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, gkf->rho_c, gkf->rho_c_global_dg);

  // Smooth the charge density. Input is rho_c_global_dg, globally smoothed in z,
  // and then output should be in *local* phi.
  gkyl_fem_parproj_set_rhs(gkf->fem_parproj_core, gkf->rho_c_global_dg, gkf->rho_c_global_dg);
  gkyl_fem_parproj_solve(gkf->fem_parproj_core, gkf->rho_c_global_smooth);
  gkyl_fem_parproj_set_rhs(gkf->fem_parproj_sol, gkf->rho_c_global_dg, gkf->rho_c_global_dg);
  gkyl_fem_parproj_solve(gkf->fem_parproj_sol, gkf->rho_c_global_smooth);

  // Solve the Poisson equation.
  gkyl_deflated_fem_poisson_advance(gkf->fem_poisson_deflated, gkf->rho_c_global_smooth,
    gkf->phi_bc, phi_out);

  // Finish the Poisson solve with FLR effects.
  gkf->invert_flr(app, gkf, phi_out);

  // Enforce a BC of the field in the parallel direction.
  gkf->enforce_parallel_bc_func(app, gkf, phi_out);
}

static void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_release(gkf->rho_c);
  gkyl_array_release(gkf->rho_c_global_dg);
  gkyl_array_release(gkf->rho_c_global_smooth);
  gkyl_array_release(gkf->phi_fem);
  gkyl_array_release(gkf->phi);
  gkyl_array_release(gkf->phi1);
  gkyl_array_release(gkf->phinew);

  if (gkf->gkfield_id == GKYL_GK_FIELD_EM_APAR) {
    gkyl_array_release(gkf->apar_fem);
    gkyl_array_release(gkf->apardot_fem);
  }

  if (app->use_gpu) {
    gkyl_array_release(gkf->phi_host);
  }

  gkyl_array_release(gkf->epsilon);

  gkyl_deflated_fem_poisson_release(gkf->fem_poisson_deflated);
  gkyl_fem_poisson_perp_release(gkf->fem_poisson_perp);
  if (gkf->is_dirichletvar) {
    gkyl_array_release(gkf->phi_bc);
  }
  
  if (gkf->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(gkf->fem_parproj_core);
    gkyl_fem_parproj_release(gkf->fem_parproj_sol);
  } else {
    gkyl_fem_parproj_release(gkf->fem_parproj);
  }

  gkyl_array_integrate_release(gkf->calc_em_energy);

  // Release TS BC and SSFG updater
  if (gkf->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(gkf->bc_ts_lo);
    }
    gkyl_skin_surf_from_ghost_release(gkf->ssfg_z_lo);
  }
  
  if (gkf->use_flr) {
    gk_field_flr_release(app, gkf);
  }
}

void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&gkf->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  gkf->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  gkf->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  gkf->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  gkf->phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->phi1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->phinew = gkyl_array_acquire(gkf->phi1);

  // Package field arrays into array of array pointers.
  gkf->f[0] = gkf->phi;
  gkf->f1[0] = gkf->phi1;
  gkf->fnew[0] = gkf->phinew;

  // Allocate electromagnetic arrays if needed.
  if (gkf->gkfield_id == GKYL_GK_FIELD_EM_APAR) {
    gkf->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkf->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

//    // Package field arrays into array of array pointers.
//    gkf->f[1] = gkf->apar;
//    gkf->f1[1] = gkf->apar1;
//    gkf->fnew[1] = gkf->aparnew;
  }

  // Allocate phi_host for I/O.
  gkf->phi_host = gkf->phi;
  if (app->use_gpu) {
    gkf->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  if (gkf->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    gkf->accumulate_rhoc_func = gk_field_accumulate_rho_c_adiabatic;
  } else {
    gkf->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  }

  double polarization_weight = 0.0;
  double polarization_bmag = gkf->info.polarization_bmag ? gkf->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }
  // Allocate array for the polarization weight times geometric coefficients.
  gkf->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  
  // Initialize the polarization weight.
  struct gkyl_array *Jgij[3] = {app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj};
  for (int i=0; i<app->cdim-2/app->cdim; i++) {
    gkyl_array_set_offset(gkf->epsilon, polarization_weight, Jgij[i], i*app->basis.num_basis);
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
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(gkf->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(gkf->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
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
  poisson_bcs.contains_lower_z_edge = gkf->global_sub_range.lower[ndim-1] == app->global.lower[ndim-1];
  poisson_bcs.contains_upper_z_edge = gkf->global_sub_range.upper[ndim-1] == app->global.upper[ndim-1];

  // Initialize the Poisson solver.
  gkf->fem_poisson_deflated = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
    app->local, gkf->global_sub_range, gkf->epsilon, 0, poisson_bcs, gkf->info.bias_plane_list, app->use_gpu);
  gkf->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &poisson_bcs, gkf->epsilon, NULL, app->use_gpu);

  gkf->phi_bc = 0;
  gkf->is_dirichletvar = false;
  for (int i=0; i<2*app->cdim; i++) {
    gkf->is_dirichletvar = gkf->is_dirichletvar ||
                          (gkf->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING ||
                           gkf->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
  }

  if (gkf->is_dirichletvar) {
    // Project the spatially varying BC if the user specifies it.
    gkf->phi_bc = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    struct gkyl_array *phi_bc_ho = mkarr(false, gkf->phi_bc->ncomp, gkf->phi_bc->size);

    for (int d=0; d<app->cdim; d++) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(gkf->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_lo->aux_profile, bc_lo->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(gkf->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_up->aux_profile, bc_up->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
    }
    gkyl_array_copy(gkf->phi_bc, phi_bc_ho);
    gkyl_array_release(phi_bc_ho);
  }

  if (gkf->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkf->rhs_phi_func = gk_field_rhs_deflate_poisson_es_iwl;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core, fem_parproj_bc_sol;
    if (app->cdim == 2) {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_PERIODIC;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    } else {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_NONE;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    }

    gkf->fem_parproj_core = gkyl_fem_parproj_new(&app->global_core, &app->basis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    gkf->fem_parproj_sol = gkyl_fem_parproj_new(&app->global_sol, &app->basis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
  } else {
    gkf->rhs_phi_func = gk_field_rhs_poisson_perp_2x3x;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d) {
      if (app->periodic_dirs[d] == app->cdim-1) {
        fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
      }
    }

    gkf->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      fem_parproj_bc, 0, 0, app->use_gpu);
  }

  gkyl_array_set(gkf->es_energy_fac, 0.5, gkf->epsilon);

  gkf->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  gkf->use_flr = false;
  gkf->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type) {
      gkf->use_flr = gkf->use_flr || s->info.flr.type;
    }
  }
  if (gkf->use_flr) {
    gk_field_flr_new(app, gkf);
  }

  // Twist-and-shift boundary condition for phi and skin surface from ghost to impose phi periodicity at z=-pi.
  gkf->enforce_parallel_bc_func = gk_field_enforce_parallel_bc_disabled;
  if (gkf->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gk_field_2x3x_add_TSBC_and_SSFG_updaters(app, gkf);
    gkf->enforce_parallel_bc_func = gk_field_enforce_parallel_bc_enabled;
  }

  gkf->solver_release_func = gk_field_fem_release_2x3x;
}
