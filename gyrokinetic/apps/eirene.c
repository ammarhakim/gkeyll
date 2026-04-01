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

void
gk_eirene_rhs(gkyl_gyrokinetic_app *app, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  struct gk_eirene *eirene = app->eirene;
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    int sidx = gk_find_species_idx(app, gks->info.name);
    gk_species_source_bgk_rhs(app, gks, bgk_src, fin[sidx], rhs[sidx]);
  }
}

void
gk_eirene_read(struct gkyl_gyrokinetic_app *app, struct gkyl_array *out, cstr fileNm)
{
  struct gkyl_array* nnodal = mkarr(false, 1, app->gk_geom->nrange_int.volume);

  FILE *ptr = fopen(fileNm.str,"r");
  size_t status;

  int nr = gkyl_range_shape(&app->gk_geom->nrange_int, 0);
  int nz = gkyl_range_shape(&app->gk_geom->nrange_int, 1);
  int idx[2];

  for(int ir = 0; ir < nr; ir++){
    idx[0] = ir;
    for(int iz = 0; iz < nz; iz++){
      idx[1] = iz;
      double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&app->gk_geom->nrange_int, idx));
      status = fscanf(ptr,"%lf", nnodal_n);
    }
  }

  fclose(ptr);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
  gkyl_nodal_ops_n2m(n2m, &app->basis, &app->grid, &app->gk_geom->nrange_int, &app->local, 1, nnodal, out, true);
  gkyl_array_release(nnodal);
}

void
gk_eirene_write(struct gkyl_gyrokinetic_app *app, double tm, int frame)
{
  struct gk_eirene *eirene = app->eirene;

  // Read new EIRENE data.
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    cstr fileNm = cstr_from_fmt("%s%s-%s_M0source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M0dot_host, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M1source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M1dot_host, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M2source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M2dot_host, fileNm);
    if(app->use_gpu) {
      gkyl_array_copy(bgk_src->M0dot, bgk_src->M0dot_host);
      gkyl_array_copy(bgk_src->M1dot, bgk_src->M1dot_host);
      gkyl_array_copy(bgk_src->M2dot, bgk_src->M2dot_host);
    }
  }

  // Write out the source moments
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    gk_species_source_bgk_write_diags(app, gks, bgk_src, tm, frame);
  }

}

void
gk_eirene_calc_integrated_diagnostics(struct gkyl_gyrokinetic_app *app, double tm)
{
  struct gk_eirene *eirene = app->eirene;
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    bgk_src->calc_integrated_diags_func(app, gks, bgk_src, tm);
  }
}

void
gk_eirene_write_integrated_diagnostics(struct gkyl_gyrokinetic_app *app)
{
  struct gk_eirene *eirene = app->eirene;
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    bgk_src->write_integrated_diags_func(app, gks, bgk_src);
  }
}

struct gk_eirene*
gk_eirene_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gk *gk)
{

  struct gk_eirene *eirene = gkyl_malloc(sizeof(struct gk_eirene));

  eirene->info = gk->eirene;

  for (int i=0; i<eirene->info.num_coupling_species; ++i)
    eirene->coupling_species[i] = gk_find_species(app, eirene->info.coupling_species[i]);

  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    eirene->bgk_src[i].coupling_time = eirene->info.coupling_time[i];
    eirene->bgk_src[i].damping_factor = eirene->info.damping_factor[i];
    eirene->bgk_src[i].core_coll_factor = eirene->info.core_coll_factor[i];
    eirene->bgk_src[i].source_bgk_id = GKYL_SOURCE_BGK_EXTERNAL;
    eirene->bgk_src[i].write_diagnostics = true;
    gk_species_source_bgk_init(app, gks, &eirene->bgk_src[i]);
  }

  // Read initial EIRENE data.
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    cstr fileNm = cstr_from_fmt("%s%s-%s_M0source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M0dot_host, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M1source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M1dot_host, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M2source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M2dot_host, fileNm);
    if(app->use_gpu) {
      gkyl_array_copy(bgk_src->M0dot, bgk_src->M0dot_host);
      gkyl_array_copy(bgk_src->M1dot, bgk_src->M1dot_host);
      gkyl_array_copy(bgk_src->M2dot, bgk_src->M2dot_host);
    }
  }


  return eirene;
}

void
gk_eirene_release(struct gkyl_gyrokinetic_app *app, struct gk_eirene *eirene)
{
  for (int i=0; i<eirene->info.num_coupling_species; ++i)
    gk_species_source_bgk_release(app, &eirene->bgk_src[i]);
  gkyl_free(eirene);
}
