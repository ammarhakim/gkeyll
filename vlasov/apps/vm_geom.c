#include <gkyl_app.h>
#include <gkyl_vlasov_priv.h>


// Initialize species object.
void
vm_geom_init(struct gkyl_vm *vm_app_inp, struct gkyl_vlasov_app *app, struct vm_geom *vmg)
{

    vmg->info = vm_app_inp->geom;

    // Set the black hole parameters
    vmg->mass_bh = vmg->info.mass_bh;
    vmg->spin_bh = vmg->info.spin_bh;
}

void
vm_geom_release(const gkyl_vlasov_app* app, const struct vm_geom *vmg)
{
    
}
