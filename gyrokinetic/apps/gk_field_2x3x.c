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


void gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
}

void gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s)
{
  // Gyroaverage the density if needed.
  s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);
  gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
  // Add the background (electron) charge density.
  double n_s0 = field->info.electron_density;
  double q_s = field->info.electron_charge;
  double dg_norm = pow(sqrt(2), app->basis.ndim);
  gkyl_array_shiftc_range(field->rho_c, q_s * n_s0 * dg_norm, 0, &app->local);
}

void
gk_field_deflate_poisson_es_iwl_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
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
  field->enforce_zbc(app, field, field->phi_smooth);
}

void
gk_field_2x3x_poisson_perp_rhs(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
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
  field->enforce_zbc(app, field, field->phi_smooth);
}

void
gk_field_fem_init_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
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
  for (int d=0; d<app->cdim-1; d++) {
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
    f->rhs_phi_func = gk_field_deflate_poisson_es_iwl_rhs;
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
  } 
  else {
    f->rhs_phi_func = gk_field_2x3x_poisson_perp_rhs;
    enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d) {
      if (app->periodic_dirs[d] == app->cdim-1) {
        fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
      }
    }

    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      fem_parproj_bc, f->epsilon_global, 0, app->use_gpu);
  }

  gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);

  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);
}
