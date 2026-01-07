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
gk_field_2x3x_add_IWL_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Allocation ranges and updaters for IWL field solve.

  // Global skin and ghost ranges.
  int ghost[] = {1, 1, 1};
  int cdim = app->cdim;
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&f->global_lower_skin[dir], &f->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&f->global_upper_skin[dir], &f->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->global_ext, ghost);
  }
  // Create core and SOL parallel skin and ghost ranges.
  int par_dir = app->cdim-1; // Parallel direction index.
  int idx_LCFS_lo = app->gk_geom->idx_LCFS_lo;
  int len_core = idx_LCFS_lo;
  int len_sol = app->global.upper[0]-len_core;
  for (int e=0; e<2; e++) {
    gkyl_range_shorten_from_above(e==0? &f->global_lower_skin_par_core  : &f->global_upper_skin_par_core,
                                  e==0? &f->global_lower_skin[par_dir]  : &f->global_upper_skin[par_dir], 0, len_core);
    gkyl_range_shorten_from_above(e==0? &f->global_lower_ghost_par_core : &f->global_upper_ghost_par_core,
                                  e==0? &f->global_lower_ghost[par_dir] : &f->global_upper_ghost[par_dir], 0, len_core);
    gkyl_range_shorten_from_below(e==0? &f->global_lower_skin_par_sol   : &f->global_upper_skin_par_sol,
                                  e==0? &f->global_lower_skin[par_dir]  : &f->global_upper_skin[par_dir], 0, len_sol);
    gkyl_range_shorten_from_below(e==0? &f->global_lower_ghost_par_sol  : &f->global_upper_ghost_par_sol,
                                  e==0? &f->global_lower_ghost[par_dir] : &f->global_upper_ghost[par_dir], 0, len_sol);
  }
  
  if (app->cdim == 3) {
    // Take the TS function from the parallel BC of the first species.
    struct gk_species *gks = &app->species[0];
    const struct gkyl_gyrokinetic_bc *par_lower_bc;
    for (int i = 0; i < 2*app->cdim; i++) {
      if (gks->info.bcs[i].dir == par_dir && gks->info.bcs[i].edge == GKYL_LOWER_EDGE) {
        par_lower_bc = (const struct gkyl_gyrokinetic_bc *) &gks->info.bcs[i];
        break;
      }
    }

    // TS BC updater for up to low TS for the lower edge. This sets ghost_L = T_LU(ghost_L).
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
    f->bc_T_LU_lo = gkyl_bc_twistshift_new(&T_LU_lo);

    long buff_sz = GKYL_MAX2(f->global_lower_ghost_par_sol.volume, f->global_lower_ghost_par_core.volume);
    f->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    f->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      &app->basis, &f->global_upper_skin_par_core, &f->global_upper_ghost_par_core, 1, app->cdim, app->use_gpu);
  }

}

static void
gk_field_fem_projection_par_iwl_2x(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x),
  // using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  int par_dir = app->cdim - 1;
  // Apply periodicity in the core.
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &field->global_lower_ghost_par_core, &field->global_upper_skin_par_core);
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &field->global_upper_ghost_par_core, &field->global_lower_skin_par_core);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_iwl_3x(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x),
  // using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  int par_dir = app->cdim - 1;
  // Apply TS BC in the core lower parallel boundary, and
  // fill core upper parallel boundary ghost with skin boundary value.
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &field->global_lower_ghost_par_core, &field->global_upper_skin_par_core);
  gkyl_bc_twistshift_advance(field->bc_T_LU_lo, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_up, field->bc_buffer, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_rhs_poisson_perp_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Smooth the charge density along z.
  field->fem_projection_par_pre_func(app, field, field->rho_c, field->rho_c);

  // Solve the Poisson equation.
  gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c);
  gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_smooth);

  // Smooth the potential along z.
  field->fem_projection_par_post_func(app, field, field->phi_smooth, field->phi_smooth);

  // Finish the Poisson solve with FLR effects.
  field->invert_flr(app, field, field->phi_smooth);
}

static void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_array_release(f->epsilon);

  gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
  if (f->is_dirichletvar) {
    gkyl_array_release(f->phi_bc);
  }
  
  gkyl_fem_parproj_release(f->fem_parproj);
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(f->fem_parproj_core);
    gkyl_fem_parproj_release(f->fem_parproj_sol);
  }

  gkyl_array_integrate_release(f->calc_em_energy);

  // Release TS BC and SSFG updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_T_LU_lo);
      gkyl_bc_basic_gyrokinetic_release(f->gfss_bc_op_core_up);
      gkyl_array_release(f->bc_buffer);
    }
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

  // Allocate electromagnetic arrays if needed.
  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
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

  // Initialize the Poisson solver.
  f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &poisson_bcs, f->info.bias_line_list, f->epsilon, NULL, app->use_gpu);

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

  f->rhs_phi_func = gk_field_rhs_poisson_perp_2x3x;

  // Potential smoothing (in z) updater
  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d=0; d<app->num_periodic_dir; ++d)
    if (app->periodic_dirs[d] == app->cdim-1) fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;

  f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
    fem_parproj_bc, 0, 0, app->use_gpu);

  f->fem_projection_par_pre_func = gk_field_fem_projection_par;
  f->fem_projection_par_post_func = gk_field_fem_projection_par;

  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_sol = GKYL_FEM_PARPROJ_DIRICHLET_SKIN;

    f->fem_parproj_core = gkyl_fem_parproj_new(&app->global_core, &app->basis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&app->global_sol, &app->basis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
      
    if (app->cdim == 2)
      f->fem_projection_par_post_func = gk_field_fem_projection_par_iwl_2x;
    else if (app->cdim == 3)
      f->fem_projection_par_post_func = gk_field_fem_projection_par_iwl_3x;
  } 

  // Updater for field energy calculation.
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
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gk_field_2x3x_add_IWL_updaters(app,f);
  }

  f->solver_release_func = gk_field_fem_release_2x3x;
}
