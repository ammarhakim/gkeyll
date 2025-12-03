#include <gkyl_bc_twistshift.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

void
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

void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout)
{
  // Do nothing.
}
