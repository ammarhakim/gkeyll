#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_dg_updater_moment.h>

void 
vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_boundary_fluxes *bflux)
{ 
  // Allocate solver.
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&vms->grid, vms->eqn, app->cdim, app->use_gpu);
  int ndim = app->cdim + app->vdim;
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    cells[d] = vms->grid.cells[d];
    lower[d] = vms->grid.lower[d];
    upper[d] = vms->grid.upper[d];
    ghost[d] = 0;
  }

  // Initialize moment solver.
  for (int i=0; i<app->cdim; ++i) {
    cells[i] = 1;

    bflux->flux_arr[2*i] = mkarr(app->use_gpu, vms->basis.num_basis, vms->lower_ghost[i].volume);
    bflux->flux_arr[2*i+1] = mkarr(app->use_gpu, vms->basis.num_basis, vms->upper_ghost[i].volume);

    gkyl_range_init(&bflux->flux_r[2*i], ndim, vms->lower_ghost[i].lower, vms->lower_ghost[i].upper);
    gkyl_range_init(&bflux->flux_r[2*i+1], ndim, vms->upper_ghost[i].lower, vms->upper_ghost[i].upper);

    gkyl_range_init(&bflux->conf_r[2*i], app->cdim, vms->lower_ghost[i].lower,
      vms->lower_ghost[i].upper);
    gkyl_range_init(&bflux->conf_r[2*i+1], app->cdim, vms->upper_ghost[i].lower,
      vms->upper_ghost[i].upper);

    upper[i] = vms->grid.lower[i] + vms->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i], ndim, lower, upper, cells);

    upper[i] = vms->grid.upper[i];
    lower[i] = vms->grid.upper[i] - vms->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i+1], ndim, lower, upper, cells);

    bflux->integ_moms[2*i] = gkyl_dg_updater_moment_new(&bflux->boundary_grid[2*i],
      &app->basis, &vms->basis, &bflux->conf_r[2*i], &vms->local_vel, &vms->local, 
      vms->model_id, 0, GKYL_F_MOMENT_M0M1M2, true, app->use_gpu);
    bflux->integ_moms[2*i+1] = gkyl_dg_updater_moment_new(&bflux->boundary_grid[2*i+1],
      &app->basis, &vms->basis, &bflux->conf_r[2*i+1], &vms->local_vel, &vms->local, 
      vms->model_id, 0, GKYL_F_MOMENT_M0M1M2, true, app->use_gpu);

    cells[i] = vms->grid.cells[i];

    bflux->mom_arr[2*i] = mkarr(app->use_gpu, app->basis.num_basis, bflux->conf_r[2*i].volume);
    bflux->mom_arr[2*i+1] = mkarr(app->use_gpu, app->basis.num_basis, bflux->conf_r[2*i+1].volume);
  }
}

// computes rhs of the boundary flux
void
vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *vms, 
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // zero ghost cells before calculation to ensure there's no residual data
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_clear_range(rhs, 0.0, &(vms->lower_ghost[j]));
    gkyl_array_clear_range(rhs, 0.0, &(vms->upper_ghost[j]));
  }
  // ghost cells of the rhs array are filled with the bflux
  // This is overwritten by the boundary conditions and is not being stored,
  // it is only currently used to calculate moments for other applications
  if (app->use_gpu) {
    gkyl_ghost_surf_calc_advance_cu(bflux->flux_slvr, &vms->local_ext, fin, rhs);
  } else {
    gkyl_ghost_surf_calc_advance(bflux->flux_slvr, &vms->local_ext, fin, rhs);
  }

  // only calculating integrated moments for use in the bflux source for now,
  // others can be added if applications require
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j], rhs, &bflux->flux_r[2*j],
      &vms->lower_ghost[j]);
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j+1], rhs, &bflux->flux_r[2*j+1],
      &vms->upper_ghost[j]);
    
    gkyl_dg_updater_moment_advance(bflux->integ_moms[2*j], &bflux->flux_r[2*j],
      &bflux->conf_r[2*j], bflux->flux_arr[2*j], bflux->mom_arr[2*j]);
    gkyl_dg_updater_moment_advance(bflux->integ_moms[2*j+1], &bflux->flux_r[2*j+1],
      &bflux->conf_r[2*j+1], bflux->flux_arr[2*j+1], bflux->mom_arr[2*j+1]);
  }
}

void
vm_species_bflux_release(const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    gkyl_array_release(bflux->mom_arr[i]);
    gkyl_array_release(bflux->flux_arr[i]);
    gkyl_dg_updater_moment_release(bflux->integ_moms[i]);
  }
}
