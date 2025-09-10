#include <gkyl_array_rio.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkylt_mirrorgridgen.h>

#include <stc/cstr.h>

static void
split_mirror_grid_data(
  struct gkyl_range *range, const struct gkyl_mirror_grid_gen *mirror_grid,
  struct gkyl_array *e1, struct gkyl_array *e2, struct gkyl_array *e3,
  struct gkyl_array *de1, struct gkyl_array *de2, struct gkyl_array *de3,
  struct gkyl_array *B, struct gkyl_array *Bmag,
  struct gkyl_array *J)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    double *e1p = gkyl_array_fetch(e1, loc);
    double *e2p = gkyl_array_fetch(e2, loc);
    double *e3p = gkyl_array_fetch(e3, loc);

    double *de1p = gkyl_array_fetch(e1, loc);
    double *de2p = gkyl_array_fetch(e2, loc);
    double *de3p = gkyl_array_fetch(e3, loc);

    double *Bp = gkyl_array_fetch(B, loc);
    double *Bmagp = gkyl_array_fetch(Bmag, loc);
    
    double *Jp = gkyl_array_fetch(J, loc);

    const double *xc = gkyl_array_cfetch(mirror_grid->nodes_rz, loc);

    const struct gkyl_mirror_grid_gen_geom *g =
      gkyl_array_cfetch(mirror_grid->nodes_geom, loc);

    // just copy various fields over
    for (int d=0; d<3; ++d) {
      e1p[d] = g->tang[0].x[d];
      e2p[d] = g->tang[1].x[d];
      e3p[d] = g->tang[2].x[d];

      de1p[d] = g->dual[0].x[d];
      de2p[d] = g->dual[1].x[d];
      de3p[d] = g->dual[2].x[d];
      
      Bp[d] = g->B.x[d];
    }

    Bmagp[0] = gkyl_vec3_len(
      gkyl_vec3_polar_con_to_cart(xc[0], 0.0, g->B)
    );
    
    Jp[0] = g->Jc;
  }
}

bool
gkylt_mirrorgridgen(const struct gkylt_mirrorgridgen_inp *inp)
{
  double clower[3], cupper[3];
  int cells[3];
  for (int d=0; d<3; ++d) {
    clower[d] = inp->lower[d];
    cupper[d] = inp->upper[d];
    cells[d] = inp->cells[d];
  }
  
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, 3, clower, cupper, cells);

  bool status = true;  
  if (!gkyl_check_file_exists(inp->psiRZ_fname)) {
    fprintf(stderr, "Unable to find file %s!\n", inp->psiRZ_fname);
    status = false;
    goto cleanup;
  }

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, inp->psiRZ_fname);

    // write DG projection of mapc2p to file
  cstr cubic_fileNm = cstr_from_fmt("%s-psi-cubic.gkyl", inp->out_prefix);
  
  // create mirror geometry
  struct gkyl_mirror_grid_gen *mirror_grid =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = inp->fl_coord,
        .include_axis = inp->include_axis,
        .write_psi_cubic = inp->write_psi_cubic,
        .psi_cubic_fname = cubic_fileNm.str
      }
    );

  cstr_drop(&cubic_fileNm);

  // write out node coordinates
  struct gkyl_rect_grid nodal_grid;
  gkyl_rect_grid_init(&nodal_grid, 2,
    (double[]) { 0.0, 0.0 }, // arbitrary
    (double[]) { 1.0, 1.0 }, // arbitrary
    (int[]) { cells[0]+1, cells[2]+1 }
  );

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 2, (int[2]) { cells[0]+1, cells[2]+1 });

  // for outputing data
  struct gkyl_array *e1 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *e2 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *e3 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *de1 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *de2 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *de3 = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *B = gkyl_array_new(GKYL_DOUBLE, 3, node_range.volume);
  struct gkyl_array *Bmag = gkyl_array_new(GKYL_DOUBLE, 1, node_range.volume);
  struct gkyl_array *J = gkyl_array_new(GKYL_DOUBLE, 1, node_range.volume);

  split_mirror_grid_data(&node_range, mirror_grid,
    e1, e2, e3, de1, de2, de3, B, Bmag, J
  );

  struct gkyl_msgpack_data *mpd = gkyl_msgpack_create(2, (struct gkyl_msgpack_map_elem []) {
      { .key = "include_axis",
        .elem_type = GKYL_MP_INT,
        .ival = inp->include_axis },
      { .key = "field_line_coordinate",
        .elem_type = GKYL_MP_STRING,
        .cval = inp->fl_coord == GKYL_MIRROR_GRID_GEN_PSI_CART_Z ? "psi_cart_z" : "sqrt_psi_cart_z" },
    }
  );
  
  do {
    cstr fileNm = cstr_from_fmt("%s-nodal_coords.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, mirror_grid->nodes_rz, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);

  do {
    cstr fileNm = cstr_from_fmt("%s-e1.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, e1, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  do {
    cstr fileNm = cstr_from_fmt("%s-e2.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, e2, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  do {
    cstr fileNm = cstr_from_fmt("%s-e3.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, e3, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);

  do {
    cstr fileNm = cstr_from_fmt("%s-de1.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, de1, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  do {
    cstr fileNm = cstr_from_fmt("%s-de2.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, de2, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  do {
    cstr fileNm = cstr_from_fmt("%s-de3.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, de3, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);  

  do {
    cstr fileNm = cstr_from_fmt("%s-B.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, B, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  do {
    cstr fileNm = cstr_from_fmt("%s-Bmag.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, Bmag, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  
  do {
    cstr fileNm = cstr_from_fmt("%s-Jc.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, mpd, J, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);    
  
  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_array_release(psi);

  gkyl_array_release(e1);
  gkyl_array_release(e2);
  gkyl_array_release(e3);

  gkyl_array_release(de1);
  gkyl_array_release(de2);
  gkyl_array_release(de3);

  gkyl_array_release(B);
  gkyl_array_release(Bmag);  
  gkyl_array_release(J);

  gkyl_msgpack_data_release(mpd);

  cleanup:
  
  return status;
}
