#include <gkyl_app.h>
#include <gkyl_vlasov_priv.h>


// Initialize species object.
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
}

void
vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *vmg)
{
    
}
