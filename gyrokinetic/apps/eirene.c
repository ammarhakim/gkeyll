#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>

#include <assert.h>
#include <time.h>


struct gk_eirene*
gk_eirene_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gk *gk)
{

  struct gk_eirene *eirene = gkyl_malloc(sizeof(struct gk_eirene));

  eirene->info = gk->eirene;

  for (int i=0; i<eirene->info.num_coupling_species; ++i)
    eirene->coupling_species[i] = gk_find_species(app, eirene->info.coupling_species[i]);

  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    gks->bgk_src = (struct gk_source_bgk) { };
    gk_species_source_bgk_init(app, gks, &gks->bgk_src);
  }



  return eirene;
}

//gkyl_gyrokinetic_eirene_advance()
//{
//}
//
//gkyl_gyrokinetic_eirene_release()
//{
//}
//
//gkyl_gyrokinetic_eirene_read()
//{
//}
//
//gkyl_gyrokinetic_eirene_write()
//{
//  // Write Gkeyll data.
//
//
//  // Read new EIRENE data.
//  gkyl_gyrokinetic_eirene_read()
//
//  // Reset the BGK source.
//  gkyl_gk_species_source_bgk_reset(&gk_eirene.gk_src, ...)
//}
