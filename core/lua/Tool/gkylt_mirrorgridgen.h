#pragma once

#include <gkyl_mirror_grid_gen.h>

struct gkylt_mirrorgridgen_inp {
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)

  double lower[3], upper[3]; // lower and upper bounds of computational space
  int cells[3]; // number of cells in computational space
  
  bool write_psi_cubic; // set to true to write the cubic fit to file
  const char *psiRZ_fname; // name for file with psi(R,Z)
  const char *out_prefix; // output prefix
};

bool gkylt_mirrorgridgen(const struct gkylt_mirrorgridgen_inp *mginp);
