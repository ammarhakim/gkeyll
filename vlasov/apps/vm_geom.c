#include <gkyl_app.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>


// Create the dg-gr-maxwell geometry
static void 
vm_dg_maxwell_geom_new(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_geom *vmg)
{

  // populate the geometry context struct with mass and spin
  struct gkyl_dg_gr_maxwell_geom_ctx ctx = {
    .mass_bh = vmg->mass_bh,
    .spin_bh = vmg->spin_bh
  };

  // Evaluation of geometry at surface and volume nodal points.
  // Lapse - \alpha in the ADM split
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* lapse_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &app->grid, &app->basis, 1, vm_app_inp->poly_order, gkyl_dg_gr_maxwell_preset_lapse(vmg->triad_preset_geom_type), &ctx);
  vmg->lapse = gkyl_surf_and_vol_node_arrays_new(lapse_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(lapse_proj, 0.0, &app->local, vmg->lapse);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(lapse_proj);

  // shift - \beta^i components in the ADM split (contravariant)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* shift_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &app->grid, &app->basis, 3, vm_app_inp->poly_order, gkyl_dg_gr_maxwell_preset_shift(vmg->triad_preset_geom_type), &ctx);
  vmg->shift = gkyl_surf_and_vol_node_arrays_new(shift_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(shift_proj, 0.0, &app->local, vmg->shift);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(shift_proj);

  // h_ij - Covariant components of the spatial metric (assumed to allways be a upper 
  // triangular matrix of 6 unique elements)
  // Allocate arrays for specified metric inverse
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &app->grid, &app->basis, 6, vm_app_inp->poly_order, gkyl_dg_gr_maxwell_preset_h_ij(vmg->triad_preset_geom_type), &ctx);
  vmg->h_ij = gkyl_surf_and_vol_node_arrays_new(h_ij_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_proj, 0.0, &app->local, vmg->h_ij);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_proj);

  // Allocate arrays for the metric determinant (computed from J = sqrt(det(h_ij)))
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* det_h_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &app->grid, &app->basis, 1, vm_app_inp->poly_order, gkyl_dg_gr_maxwell_preset_det_h(vmg->triad_preset_geom_type), &ctx);
  vmg->det_h = gkyl_surf_and_vol_node_arrays_new(det_h_proj, app->local_ext.volume, app->use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(det_h_proj, 0.0, &app->local, vmg->det_h);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(det_h_proj);

}


// Initialize geom object.
void
vm_geom_init(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_geom *vmg)
{

  *vmg = (struct vm_geom) { 0 };
  vmg->info = vm_app_inp->geom;

  // Set the black hole parameters
  vmg->mass_bh = vmg->info.mass_bh;
  vmg->spin_bh = vmg->info.spin_bh;

  // Set triad preset-geometry parameters (used by TRIAD/TRIAD_GR models).
  vmg->use_preset_geom = vmg->info.use_preset_geom;
  vmg->triad_preset_geom_type = vmg->info.triad_preset_geom_type;

  // Determine if gr-maxwell are usd
  vmg->has_gr_fields = false;
  if ( (vm_app_inp->skip_field == false) && ( vm_app_inp->field.field_id == GKYL_FIELD_GR_D_B ) ) {
    vmg->has_gr_fields = true;
  }

  // If the fields are included, and are GR, then build the gr maxwell geometry
  if ( vmg->has_gr_fields ) {
    vm_dg_maxwell_geom_new(vm_app_inp, app, vmg);
  }

}

void
vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *vmg)
{

  // Release all variables associated with vm geom for dg-gr-maxwells
  if ( vmg->has_gr_fields ) {
    gkyl_surf_and_vol_node_arrays_release(vmg->lapse);
    gkyl_surf_and_vol_node_arrays_release(vmg->shift);
    gkyl_surf_and_vol_node_arrays_release(vmg->h_ij);
    gkyl_surf_and_vol_node_arrays_release(vmg->det_h);
  }
}
