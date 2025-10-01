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


void
gk_field_poisson_solve_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve the Poisson equation in 1x with the parallel FEM projection.
  gk_field_fem_projection_par(app, field, field->rho_c, field->phi_smooth);
}

void
gk_field_fem_init_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  f->field_solve = gk_field_poisson_solve_1x;

  double polarization_weight = 0.0;
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }

  double es_energy_fac_1d_adiabatic = 0.0;
  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    // Add the contribution from adiabatic electrons
    double n_s0 = f->info.electron_density;
    double q_s = f->info.electron_charge;
    double T_s = f->info.electron_temp;
    double quasineut_contr = q_s*n_s0*q_s/T_s;
    es_energy_fac_1d_adiabatic = 0.5*quasineut_contr;
    
    struct gkyl_array *epsilon_adiab = mkarr(app->use_gpu, f->epsilon->ncomp, f->epsilon->size);
    gkyl_array_copy(epsilon_adiab, app->gk_geom->geo_int.jacobgeo);
    gkyl_array_scale(epsilon_adiab, quasineut_contr);
    gkyl_array_accumulate(f->epsilon, 1., epsilon_adiab);
    gkyl_array_release(epsilon_adiab);
  }

  // Allocate array for the polarization weight times geometric coefficients.
  f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);

  // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
  gkyl_array_copy(f->epsilon, app->gk_geom->geo_int.jacobgeo);
  gkyl_array_scale(f->epsilon, polarization_weight);
  gkyl_array_scale(f->epsilon, f->info.kperpSq);

  // Gather epsilon for (global) smoothing in z.
  f->epsilon_global = mkarr(app->use_gpu, f->epsilon->ncomp, app->global_ext.volume);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, f->epsilon, f->epsilon_global);

  // Potential smoothing (in z) updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core, fem_parproj_bc_sol;
    fem_parproj_bc_core = GKYL_FEM_PARPROJ_NONE;
    fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;

    f->fem_parproj_core = gkyl_fem_parproj_new(&app->global_core, &app->basis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&app->global_sol, &app->basis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
  } 
  else {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d)
      if (app->periodic_dirs[d] == app->cdim-1) fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;

    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      fem_parproj_bc, f->epsilon_global, 0, app->use_gpu);
  }

  f->es_energy_fac_1d = polarization_weight*f->info.kperpSq + es_energy_fac_1d_adiabatic;

  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);
}