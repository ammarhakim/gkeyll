#include <assert.h>
#include <float.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_dg_reduce_dir.h>
#include <gkyl_array_dg_reduce_dir_priv.h>
#include <gkyl_nodal_ops.h>

/**
 * Find the extremum along the reduction direction for a given preserved-direction
 * node index, storing results in the nodal arrays.
 */
static void
find_extremum_for_preserved_node(struct gkyl_array_dg_reduce_dir *up, const struct gkyl_array *in_ho,
  int preserved_node_idx)
{
  int ndim = up->grid.ndim;
  int reduce_dir = up->reduce_dir;
  int poly_order = up->basis.poly_order;

  // Determine number of nodes along reduction direction.
  int num_cells_reduce = up->range.upper[reduce_dir] - up->range.lower[reduce_dir] + 1;
  int total_nodes_reduce = (poly_order == 1) ? num_cells_reduce + 1 : 2*num_cells_reduce + 1;

  // Allocate arrays to store values and coordinates along reduction direction.
  double *vals = gkyl_malloc(sizeof(double) * total_nodes_reduce);
  double *coords = gkyl_malloc(sizeof(double) * total_nodes_reduce);
  bool *visited = gkyl_malloc(sizeof(bool) * total_nodes_reduce);
  for (int i = 0; i < total_nodes_reduce; i++) {
    vals[i] = 0.0;
    coords[i] = 0.0;
    visited[i] = false;
  }

  // For 2D, determine the preserved direction.
  int preserved_dir = (ndim == 1) ? -1 : ((reduce_dir == 0) ? 1 : 0);

  // Iterate along cells in reduction direction and collect nodal values.
  for (int cell_idx = up->range.lower[reduce_dir]; cell_idx <= up->range.upper[reduce_dir]; cell_idx++) {
    // For 2D, determine which cells in the preserved direction contribute to this node.
    int pres_cell_start, pres_cell_end;
    if (ndim == 1) {
      pres_cell_start = 0;
      pres_cell_end = 0;
    }
    else {
      if (poly_order == 1) {
        if (preserved_node_idx == 0) {
          pres_cell_start = up->range.lower[preserved_dir];
          pres_cell_end = up->range.lower[preserved_dir];
        }
        else if (preserved_node_idx == up->out_nrange.upper[0]) {
          pres_cell_start = up->range.upper[preserved_dir];
          pres_cell_end = up->range.upper[preserved_dir];
        }
        else {
          pres_cell_start = up->range.lower[preserved_dir] + preserved_node_idx - 1;
          pres_cell_end = pres_cell_start + 1;
          if (pres_cell_end > up->range.upper[preserved_dir]) {
            pres_cell_end = up->range.upper[preserved_dir];
          }
        }
      }
      else {
        int cell_local = preserved_node_idx / 2;
        pres_cell_start = up->range.lower[preserved_dir] + cell_local;
        pres_cell_end = pres_cell_start;
        if (preserved_node_idx % 2 == 0 && preserved_node_idx > 0) {
          pres_cell_start--;
        }
        if (pres_cell_start < up->range.lower[preserved_dir]) {
          pres_cell_start = up->range.lower[preserved_dir];
        }
        if (pres_cell_end > up->range.upper[preserved_dir]) {
          pres_cell_end = up->range.upper[preserved_dir];
        }
      }
    }

    for (int pres_cell = pres_cell_start; pres_cell <= pres_cell_end; pres_cell++) {
      // Build index array for this cell.
      int idx[GKYL_MAX_DIM];
      if (ndim == 1) {
        idx[0] = cell_idx;
      }
      else {
        idx[preserved_dir] = pres_cell;
        idx[reduce_dir] = cell_idx;
      }

      long linidx = gkyl_range_idx(&up->range, idx);
      const double *f_d = gkyl_array_cfetch(in_ho, linidx);

      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->grid, idx, xc);

      // Evaluate at each node in this cell.
      for (int n = 0; n < up->basis.num_basis; n++) {
        const double *nod_log = gkyl_array_cfetch(up->nodes, n);

        // Check if this node corresponds to our preserved node index.
        if (ndim > 1) {
          int pres_node_offset;
          if (poly_order == 1) {
            pres_node_offset = (nod_log[preserved_dir] < 0) ? 0 : 1;
          }
          else {
            if (nod_log[preserved_dir] < -0.5) {
              pres_node_offset = 0;
            }
            else if (nod_log[preserved_dir] > 0.5) {
              pres_node_offset = 2;
            }
            else {
              pres_node_offset = 1;
            }
          }
          int pres_cell_local = pres_cell - up->range.lower[preserved_dir];
          int this_pres_node;
          if (poly_order == 1) {
            this_pres_node = pres_cell_local + pres_node_offset;
          }
          else {
            this_pres_node = 2*pres_cell_local + pres_node_offset;
          }

          if (this_pres_node != preserved_node_idx) {
            continue;
          }
        }

        // Determine node offset in reduction direction.
        int reduce_node_offset;
        if (poly_order == 1) {
          reduce_node_offset = (nod_log[reduce_dir] < 0) ? 0 : 1;
        }
        else {
          if (nod_log[reduce_dir] < -0.5) {
            reduce_node_offset = 0;
          }
          else if (nod_log[reduce_dir] > 0.5) {
            reduce_node_offset = 2;
          }
          else {
            reduce_node_offset = 1;
          }
        }

        int cell_local = cell_idx - up->range.lower[reduce_dir];
        int reduce_node_idx;
        if (poly_order == 1) {
          reduce_node_idx = cell_local + reduce_node_offset;
        }
        else {
          reduce_node_idx = 2*cell_local + reduce_node_offset;
        }

        if (!visited[reduce_node_idx]) {
          double val = up->basis.eval_expand(nod_log, f_d);
          double nod_phys[GKYL_MAX_DIM];
          dg_reduce_dir_log_to_comp(ndim, nod_log, up->grid.dx, xc, nod_phys);

          vals[reduce_node_idx] = val;
          coords[reduce_node_idx] = nod_phys[reduce_dir];
          visited[reduce_node_idx] = true;
        }
      }
    }
  }

  // Find the global extremum.
  int extremum_idx = 0;
  double extremum_val = vals[0];
  
  for (int i = 1; i < total_nodes_reduce; i++) {
    bool is_better;
    if (up->op == GKYL_REDUCE_OP_MAX) {
      is_better = (vals[i] > extremum_val);
    }
    else {
      is_better = (vals[i] < extremum_val);
    }
    
    if (is_better) {
      extremum_val = vals[i];
      extremum_idx = i;
    }
  }

  // Store the result.
  double *val_n = gkyl_array_fetch(up->out_vals_nodal, preserved_node_idx);
  double *coord_n = gkyl_array_fetch(up->out_coords_nodal, preserved_node_idx);
  val_n[0] = extremum_val;
  coord_n[0] = coords[extremum_idx];

  gkyl_free(vals);
  gkyl_free(coords);
  gkyl_free(visited);
}

/**
 * Evaluate an input array at the extremum coordinate for a given preserved-direction
 * node index, storing result in the nodal output array.
 */
static void
eval_at_extremum_for_preserved_node(struct gkyl_array_dg_reduce_dir *up,
  const struct gkyl_array *in_ho, int preserved_node_idx)
{
  int ndim = up->grid.ndim;
  int reduce_dir = up->reduce_dir;
  int poly_order = up->basis.poly_order;
  int preserved_dir = (ndim == 1) ? -1 : ((reduce_dir == 0) ? 1 : 0);

  // Get the extremum coordinate.
  const double *extremum_coord_n = gkyl_array_cfetch(up->out_coords_nodal, preserved_node_idx);
  double extremum_coord = extremum_coord_n[0];

  // Find the cell containing this coordinate in the reduction direction.
  double point[GKYL_MAX_DIM];
  int known_idx[GKYL_MAX_DIM];
  int cell_idx[GKYL_MAX_DIM];
  
  for (int d = 0; d < ndim; d++) {
    if (d == reduce_dir) {
      point[d] = extremum_coord;
      known_idx[d] = -1;
    }
    else {
      point[d] = 0.0;
      known_idx[d] = -1;
    }
  }
  
  // If 2D, determine preserved direction cell from preserved_node_idx.
  if (ndim > 1) {
    int pres_cell;
    if (poly_order == 1) {
      if (preserved_node_idx == 0) {
        pres_cell = up->range.lower[preserved_dir];
      }
      else {
        pres_cell = up->range.lower[preserved_dir] + preserved_node_idx - 1;
        if (pres_cell > up->range.upper[preserved_dir]) {
          pres_cell = up->range.upper[preserved_dir];
        }
      }
    }
    else {
      pres_cell = up->range.lower[preserved_dir] + preserved_node_idx / 2;
    }
    known_idx[preserved_dir] = pres_cell;
    
    int pres_cell_idx[GKYL_MAX_DIM];
    for (int d = 0; d < ndim; d++) {
      pres_cell_idx[d] = (d == preserved_dir) ? pres_cell : 1;
    }
    double xc_pres[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->grid, pres_cell_idx, xc_pres);
    point[preserved_dir] = xc_pres[preserved_dir];
  }
  
  gkyl_rect_grid_find_cell(&up->grid, point, true, known_idx, cell_idx);

  // Clamp cell_idx to interior range.
  for (int d = 0; d < up->grid.ndim; d++) {
    if (cell_idx[d] < up->range.lower[d]) {
      cell_idx[d] = up->range.lower[d];
    }
    if (cell_idx[d] > up->range.upper[d]) {
      cell_idx[d] = up->range.upper[d];
    }
  }

  // Get the DG coefficients at this cell.
  long linidx = gkyl_range_idx(&up->range_ext, cell_idx);
  const double *f_d = gkyl_array_cfetch(in_ho, linidx);

  // Get cell center.
  double xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&up->grid, cell_idx, xc);

  // Convert extremum coordinate to logical space.
  double nod_log[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; d++) {
    if (d == reduce_dir) {
      nod_log[d] = 2.0 * (extremum_coord - xc[d]) / up->grid.dx[d];
    }
    else if (ndim > 1) {
      if (poly_order == 1) {
        nod_log[d] = (preserved_node_idx == 0) ? -1.0 : 1.0;
      }
      else {
        int pres_node_offset = preserved_node_idx % 3;
        if (pres_node_offset == 0) {
          nod_log[d] = -1.0;
        }
        else if (pres_node_offset == 1) {
          nod_log[d] = 0.0;
        }
        else {
          nod_log[d] = 1.0;
        }
      }
    }
  }

  // Evaluate the DG expansion at this logical coordinate.
  double val = up->basis.eval_expand(nod_log, f_d);

  // Store the result.
  double *val_n = gkyl_array_fetch(up->out_eval_nodal, preserved_node_idx);
  val_n[0] = val;
}

struct gkyl_array_dg_reduce_dir*
gkyl_array_dg_reduce_dir_new(const struct gkyl_array_dg_reduce_dir_inp *inp)
{
  struct gkyl_array_dg_reduce_dir *up = gkyl_malloc(sizeof(*up));

  // Copy input parameters.
  up->grid = *inp->grid;
  up->basis = *inp->basis;
  up->range = *inp->range;
  up->range_ext = *inp->range_ext;
  up->reduce_dir = inp->reduce_dir;
  up->op = inp->op;
  up->use_gpu = inp->use_gpu;

  int ndim = inp->grid->ndim;
  int poly_order = inp->basis->poly_order;
  int out_dim = ndim - 1;

  assert(inp->reduce_dir >= 0 && inp->reduce_dir < ndim);

  // Set up output grid/basis/range.
  if (out_dim == 0) {
    // 1D -> 0D case.
    int cells_1d[1] = {1};
    double lower_1d[1] = {0.0};
    double upper_1d[1] = {1.0};
    gkyl_rect_grid_init(&up->out_grid, 1, lower_1d, upper_1d, cells_1d);
    gkyl_range_init(&up->out_range, 1, (int[]){1}, (int[]){1});
    gkyl_range_init(&up->out_range_ext, 1, (int[]){0}, (int[]){2});
    gkyl_cart_modal_serendip(&up->out_basis, 1, 0);

    int nodes_shape[1] = {1};
    gkyl_range_init_from_shape(&up->out_nrange, 1, nodes_shape);
  }
  else if (out_dim == 1) {
    // 2D -> 1D case.
    int preserved_dir = (inp->reduce_dir == 0) ? 1 : 0;

    int cells_out = inp->grid->cells[preserved_dir];
    double lower_out = inp->grid->lower[preserved_dir];
    double upper_out = inp->grid->upper[preserved_dir];

    gkyl_rect_grid_init(&up->out_grid, 1, &lower_out, &upper_out, &cells_out);

    int lower_idx[1] = {inp->range->lower[preserved_dir]};
    int upper_idx[1] = {inp->range->upper[preserved_dir]};
    gkyl_range_init(&up->out_range, 1, lower_idx, upper_idx);

    int lower_ext_idx[1] = {inp->range_ext->lower[preserved_dir]};
    int upper_ext_idx[1] = {inp->range_ext->upper[preserved_dir]};
    gkyl_range_init(&up->out_range_ext, 1, lower_ext_idx, upper_ext_idx);

    gkyl_cart_modal_serendip(&up->out_basis, 1, poly_order);

    int num_nodes = (poly_order == 1) ? gkyl_range_shape(&up->out_range, 0) + 1
                                      : 2*gkyl_range_shape(&up->out_range, 0) + 1;
    int nodes_shape[1] = {num_nodes};
    gkyl_range_init_from_shape(&up->out_nrange, 1, nodes_shape);
  }
  else {
    assert(false); // Unsupported dimension
  }

  // Store node locations for input basis.
  up->nodes = gkyl_array_new(GKYL_DOUBLE, ndim, inp->basis->num_basis);
  inp->basis->node_list(gkyl_array_fetch(up->nodes, 0));

  // Create nodal-to-modal converter.
  up->n2m = gkyl_nodal_ops_new(&up->out_basis, &up->out_grid, false);

  // Allocate output arrays.
  up->out_vals = gkyl_array_new(GKYL_DOUBLE, up->out_basis.num_basis, up->out_range_ext.volume);
  up->out_coords = gkyl_array_new(GKYL_DOUBLE, up->out_basis.num_basis, up->out_range_ext.volume);
  up->out_vals_nodal = gkyl_array_new(GKYL_DOUBLE, 1, up->out_nrange.volume);
  up->out_coords_nodal = gkyl_array_new(GKYL_DOUBLE, 1, up->out_nrange.volume);
  up->out_eval_nodal = gkyl_array_new(GKYL_DOUBLE, 1, up->out_nrange.volume);

  return up;
}

void
gkyl_array_dg_reduce_dir_advance(struct gkyl_array_dg_reduce_dir *up, const struct gkyl_array *in)
{
  // Needs a GPU implementation.

  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;

  // Find extremum for each preserved-direction node.
  int num_nodes_out = up->out_nrange.volume;
  for (int pres_node = 0; pres_node < num_nodes_out; pres_node++) {
    find_extremum_for_preserved_node(up, in, pres_node);
  }

  // Transform nodal to modal.
  if (out_dim == 0) {
    // 1D -> 0D case: modal = nodal (p=0 has no nodal_to_modal function).
    double *val_m = gkyl_array_fetch(up->out_vals, 0);
    double *coord_m = gkyl_array_fetch(up->out_coords, 0);
    const double *val_n = gkyl_array_cfetch(up->out_vals_nodal, 0);
    const double *coord_n = gkyl_array_cfetch(up->out_coords_nodal, 0);
    val_m[0] = val_n[0];
    coord_m[0] = coord_n[0];
  }
  else {
    // 2D -> 1D case: use nodal-to-modal transform.
    gkyl_nodal_ops_n2m(up->n2m, &up->out_basis, &up->out_grid,
      &up->out_nrange, &up->out_range, 1, up->out_vals_nodal, up->out_vals, false);
    gkyl_nodal_ops_n2m(up->n2m, &up->out_basis, &up->out_grid,
      &up->out_nrange, &up->out_range, 1, up->out_coords_nodal, up->out_coords, false);
  }
}

void
gkyl_array_dg_reduce_dir_eval_at_extremum(struct gkyl_array_dg_reduce_dir *up,
  const struct gkyl_array *in_array, struct gkyl_array *out_val)
{
  // Needs a GPU implementation.

  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;

  // Evaluate the input array at extremum locations for each preserved-direction node.
  int num_nodes_out = up->out_nrange.volume;
  for (int pres_node = 0; pres_node < num_nodes_out; pres_node++) {
    eval_at_extremum_for_preserved_node(up, in_array, pres_node);
  }

  // Transform nodal to modal.
  if (out_dim == 0) {
    // 1D -> 0D case: modal = nodal.
    double *val_m = gkyl_array_fetch(out_val, 0);
    const double *val_n = gkyl_array_cfetch(up->out_eval_nodal, 0);
    val_m[0] = val_n[0];
  }
  else {
    // 2D -> 1D case: use nodal-to-modal transform.
    gkyl_nodal_ops_n2m(up->n2m, &up->out_basis, &up->out_grid,
      &up->out_nrange, &up->out_range, 1, up->out_eval_nodal, out_val, false);
  }
}

const struct gkyl_basis*
gkyl_array_dg_reduce_dir_get_basis(const struct gkyl_array_dg_reduce_dir *up)
{
  return &up->out_basis;
}

const struct gkyl_rect_grid*
gkyl_array_dg_reduce_dir_get_grid(const struct gkyl_array_dg_reduce_dir *up)
{
  return &up->out_grid;
}

const struct gkyl_range*
gkyl_array_dg_reduce_dir_get_range(const struct gkyl_array_dg_reduce_dir *up)
{
  return &up->out_range;
}

const struct gkyl_range*
gkyl_array_dg_reduce_dir_get_range_ext(const struct gkyl_array_dg_reduce_dir *up)
{
  return &up->out_range_ext;
}

const struct gkyl_range*
gkyl_array_dg_reduce_dir_get_nodal_range(const struct gkyl_array_dg_reduce_dir *up)
{
  return &up->out_nrange;
}

const struct gkyl_array*
gkyl_array_dg_reduce_dir_get_vals(const struct gkyl_array_dg_reduce_dir *up)
{
  return up->out_vals;
}

const struct gkyl_array*
gkyl_array_dg_reduce_dir_get_vals_nodal(const struct gkyl_array_dg_reduce_dir *up)
{
  return up->out_vals_nodal;
}

const struct gkyl_array*
gkyl_array_dg_reduce_dir_get_coords(const struct gkyl_array_dg_reduce_dir *up)
{
  return up->out_coords;
}

const struct gkyl_array*
gkyl_array_dg_reduce_dir_get_coords_nodal(const struct gkyl_array_dg_reduce_dir *up)
{
  return up->out_coords_nodal;
}

void
gkyl_array_dg_reduce_dir_release(struct gkyl_array_dg_reduce_dir *up)
{
  gkyl_array_release(up->out_vals);
  gkyl_array_release(up->out_coords);
  gkyl_array_release(up->out_vals_nodal);
  gkyl_array_release(up->out_coords_nodal);
  gkyl_array_release(up->out_eval_nodal);
  gkyl_array_release(up->nodes);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_free(up);
}
