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
    gk_species_source_bgk_rhs(app, gks, bgk_src, fin[i], rhs[i]);
  }
}

void
gk_eirene_read(struct gkyl_gyrokinetic_app *app, struct gkyl_array *out, cstr fileNm)
{
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &app->local, app->poly_order);
  struct gkyl_array* nnodal = mkarr(false, 1, nrange.volume);

  FILE *ptr = fopen(fileNm.str,"r");
  size_t status;

  int nr = gkyl_range_shape(&nrange, 0);
  int nz = gkyl_range_shape(&nrange, 1);
  int idx[2];

  for(int ir = 0; ir < nr; ir++){
    idx[0] = ir;
    for(int iz = 0; iz < nz; iz++){
      idx[1] = iz;
      double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", nnodal_n);
    }
  }

  fclose(ptr);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
  gkyl_nodal_ops_n2m(n2m, &app->basis, &app->grid, &nrange, &app->local, 1, nnodal, out, false);
  gkyl_array_release(nnodal);
}

void
gk_eirene_write(struct gkyl_gyrokinetic_app *app, int frame)
{
  struct gk_eirene *eirene = app->eirene;
  // Write Gkeyll frame.
  cstr fileNm = cstr_from_fmt("%snew_data_flag", eirene->info.output_data_path);
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    FILE *fp = fopen(fileNm.str, "w");
    if (fp == NULL)
        return;
    fprintf(fp, "%d\n", frame);
    fclose(fp);
  }

  // Read new EIRENE data.
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    cstr fileNm = cstr_from_fmt("%s%s-%s_M0source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M0dot, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M1source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M1dot, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M2source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M2dot, fileNm);
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
    eirene->bgk_src[i].coupling_time = eirene->info.coupling_time;
    eirene->bgk_src[i].source_bgk_id = GKYL_SOURCE_BGK_EXTERNAL;
    gk_species_source_bgk_init(app, gks, &eirene->bgk_src[i]);
  }

  // Read initial EIRENE data.
  for (int i=0; i<eirene->info.num_coupling_species; ++i) {
    struct gk_species *gks = eirene->coupling_species[i];
    struct gk_source_bgk *bgk_src = &eirene->bgk_src[i];
    cstr fileNm = cstr_from_fmt("%s%s-%s_M0source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M0dot, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M1source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M1dot, fileNm);
    fileNm = cstr_from_fmt("%s%s-%s_M2source.txt", eirene->info.input_data_path, app->name, gks->info.name);
    gk_eirene_read(app, bgk_src->M2dot, fileNm);
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
