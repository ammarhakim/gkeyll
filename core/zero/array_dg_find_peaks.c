#include <assert.h>
#include <float.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_array_dg_find_peaks_priv.h>
#include <gkyl_nodal_ops.h>

/**
 * Scan along the search direction at a fixed preserved-direction coordinate
 * to count the number of peaks and determine their types.
 * 
 * @param up Updater (partially initialized - just grid/basis/range/search_dir)
 * @param in Input field
 * @param preserved_idx Index in the preserved direction (ignored for 1D)
 * @param num_peaks_out Output: number of peaks found
 * @param peak_types_out Output: array of peak types (must be size GKYL_DG_FIND_PEAKS_MAX)
 */
static void
count_peaks_along_dir(const struct gkyl_array_dg_find_peaks *up, const struct gkyl_array *in,
  int preserved_idx, int *num_peaks_out, enum gkyl_peak_type *peak_types_out)
{
  int ndim = up->grid.ndim;
  int search_dir = up->search_dir;
  int poly_order = up->basis.poly_order;

  // Get node locations.
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, up->basis.num_basis);
  up->basis.node_list(gkyl_array_fetch(nodes, 0));

  // Copy input to host if needed.
  struct gkyl_array *in_ho = gkyl_array_new(GKYL_DOUBLE, in->ncomp, in->size);
  gkyl_array_copy(in_ho, in);

  // Determine number of nodes along search direction.
  int num_cells_search = up->range.upper[search_dir] - up->range.lower[search_dir] + 1;
  int nodes_per_cell = (poly_order == 1) ? 2 : 3;
  int total_nodes_search = (poly_order == 1) ? num_cells_search + 1 : 2*num_cells_search + 1;

  // Allocate arrays to store values and coordinates along search direction.
  double *vals = gkyl_malloc(sizeof(double) * total_nodes_search);
  double *coords = gkyl_malloc(sizeof(double) * total_nodes_search);
  for (int i = 0; i < total_nodes_search; i++) {
    vals[i] = 0.0;
    coords[i] = 0.0;
  }

  // Iterate along cells in search direction and collect nodal values.
  for (int cell_idx = up->range.lower[search_dir]; cell_idx <= up->range.upper[search_dir]; cell_idx++) {
    // Build index array for this cell.
    int idx[GKYL_MAX_DIM];
    if (ndim == 1) {
      idx[0] = cell_idx;
    }
    else {
      int preserved_dir = (search_dir == 0) ? 1 : 0;
      idx[preserved_dir] = preserved_idx;
      idx[search_dir] = cell_idx;
    }

    long linidx = gkyl_range_idx(&up->range, idx);
    const double *f_d = gkyl_array_cfetch(in_ho, linidx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->grid, idx, xc);

    // Evaluate at each node in this cell.
    for (int n = 0; n < up->basis.num_basis; n++) {
      const double *nod_log = gkyl_array_cfetch(nodes, n);
      
      // Determine node offset in search direction.
      int node_offset;
      if (poly_order == 1) {
        node_offset = (nod_log[search_dir] < 0) ? 0 : 1;
      }
      else {
        if (nod_log[search_dir] < -0.5)
          node_offset = 0;
        else if (nod_log[search_dir] > 0.5)
          node_offset = 2;
        else
          node_offset = 1;
      }

      int cell_local = cell_idx - up->range.lower[search_dir];
      int search_node_idx;
      if (poly_order == 1)
        search_node_idx = cell_local + node_offset;
      else
        search_node_idx = 2*cell_local + node_offset;

      double val = up->basis.eval_expand(nod_log, f_d);
      double nod_phys[GKYL_MAX_DIM];
      dg_find_peaks_log_to_comp(ndim, nod_log, up->grid.dx, xc, nod_phys);

      // Only store if this is the first time we see this search node
      // (avoid duplicates at cell boundaries).
      if (vals[search_node_idx] == 0.0 && coords[search_node_idx] == 0.0) {
        vals[search_node_idx] = val;
        coords[search_node_idx] = nod_phys[search_dir];
      }
    }
  }

  // Now scan the values to find peaks.
  // A peak is: EDGE_LO at index 0, EDGE_HI at last index, LOCAL_MAX/MIN in between.
  int num_peaks = 0;
  
  // Always add lower edge.
  peak_types_out[num_peaks++] = GKYL_PEAK_EDGE_LO;

  // Scan for local maxima and minima (indices 1 to total_nodes_search-2).
  for (int i = 1; i < total_nodes_search - 1; i++) {
    double prev = vals[i-1];
    double curr = vals[i];
    double next = vals[i+1];

    if (curr > prev && curr > next) {
      // Local maximum.
      assert(num_peaks < GKYL_DG_FIND_PEAKS_MAX);
      peak_types_out[num_peaks++] = GKYL_PEAK_LOCAL_MAX;
    }
    else if (curr < prev && curr < next) {
      // Local minimum.
      assert(num_peaks < GKYL_DG_FIND_PEAKS_MAX);
      peak_types_out[num_peaks++] = GKYL_PEAK_LOCAL_MIN;
    }
  }

  // Always add upper edge.
  assert(num_peaks < GKYL_DG_FIND_PEAKS_MAX);
  peak_types_out[num_peaks++] = GKYL_PEAK_EDGE_HI;

  *num_peaks_out = num_peaks;

  gkyl_free(vals);
  gkyl_free(coords);
  gkyl_array_release(nodes);
  gkyl_array_release(in_ho);
}

/**
 * Find all peaks along the search direction for a given preserved-direction
 * node index, storing results in the nodal arrays.
 */
static void
find_peaks_for_preserved_node(struct gkyl_array_dg_find_peaks *up, const struct gkyl_array *in_ho,
  int preserved_node_idx)
{
  int ndim = up->grid.ndim;
  int search_dir = up->search_dir;
  int poly_order = up->basis.poly_order;

  // Determine number of nodes along search direction.
  int num_cells_search = up->range.upper[search_dir] - up->range.lower[search_dir] + 1;
  int total_nodes_search = (poly_order == 1) ? num_cells_search + 1 : 2*num_cells_search + 1;

  // Allocate arrays to store values and coordinates along search direction.
  double *vals = gkyl_malloc(sizeof(double) * total_nodes_search);
  double *coords = gkyl_malloc(sizeof(double) * total_nodes_search);
  bool *visited = gkyl_malloc(sizeof(bool) * total_nodes_search);
  for (int i = 0; i < total_nodes_search; i++) {
    vals[i] = 0.0;
    coords[i] = 0.0;
    visited[i] = false;
  }

  // For 2D, determine the preserved direction cell index from the node index.
  int preserved_dir = (ndim == 1) ? -1 : ((search_dir == 0) ? 1 : 0);

  // Iterate along cells in search direction and collect nodal values.
  for (int cell_idx = up->range.lower[search_dir]; cell_idx <= up->range.upper[search_dir]; cell_idx++) {
    // For 2D, we need to iterate over cells in the preserved direction that
    // contribute to this preserved node index.
    int pres_cell_start, pres_cell_end;
    if (ndim == 1) {
      pres_cell_start = 0;
      pres_cell_end = 0;
    }
    else {
      // Determine which cells contribute to this preserved node.
      if (poly_order == 1) {
        // Node i is shared by cells i and i+1 (0-indexed from lower).
        // preserved_node_idx 0 is only in cell lower[preserved_dir].
        // preserved_node_idx N is only in cell upper[preserved_dir].
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
          if (pres_cell_end > up->range.upper[preserved_dir])
            pres_cell_end = up->range.upper[preserved_dir];
        }
      }
      else { // poly_order == 2
        // Similar logic for p=2 nodes.
        int cell_local = preserved_node_idx / 2;
        pres_cell_start = up->range.lower[preserved_dir] + cell_local;
        pres_cell_end = pres_cell_start;
        if (preserved_node_idx % 2 == 0 && preserved_node_idx > 0) {
          pres_cell_start--;
        }
        if (pres_cell_start < up->range.lower[preserved_dir])
          pres_cell_start = up->range.lower[preserved_dir];
        if (pres_cell_end > up->range.upper[preserved_dir])
          pres_cell_end = up->range.upper[preserved_dir];
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
        idx[search_dir] = cell_idx;
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
            if (nod_log[preserved_dir] < -0.5)
              pres_node_offset = 0;
            else if (nod_log[preserved_dir] > 0.5)
              pres_node_offset = 2;
            else
              pres_node_offset = 1;
          }
          int pres_cell_local = pres_cell - up->range.lower[preserved_dir];
          int this_pres_node;
          if (poly_order == 1)
            this_pres_node = pres_cell_local + pres_node_offset;
          else
            this_pres_node = 2*pres_cell_local + pres_node_offset;

          if (this_pres_node != preserved_node_idx)
            continue;
        }

        // Determine node offset in search direction.
        int search_node_offset;
        if (poly_order == 1) {
          search_node_offset = (nod_log[search_dir] < 0) ? 0 : 1;
        }
        else {
          if (nod_log[search_dir] < -0.5)
            search_node_offset = 0;
          else if (nod_log[search_dir] > 0.5)
            search_node_offset = 2;
          else
            search_node_offset = 1;
        }

        int cell_local = cell_idx - up->range.lower[search_dir];
        int search_node_idx;
        if (poly_order == 1)
          search_node_idx = cell_local + search_node_offset;
        else
          search_node_idx = 2*cell_local + search_node_offset;

        if (!visited[search_node_idx]) {
          double val = up->basis.eval_expand(nod_log, f_d);
          double nod_phys[GKYL_MAX_DIM];
          dg_find_peaks_log_to_comp(ndim, nod_log, up->grid.dx, xc, nod_phys);

          vals[search_node_idx] = val;
          coords[search_node_idx] = nod_phys[search_dir];
          visited[search_node_idx] = true;
        }
      }
    }
  }

  // Now extract peaks based on peak_types.
  int peak_idx = 0;
  
  // EDGE_LO is always first peak at index 0.
  if (up->peak_types[peak_idx] == GKYL_PEAK_EDGE_LO) {
    double *val_n = gkyl_array_fetch(up->out_vals_nodal[peak_idx], preserved_node_idx);
    double *coord_n = gkyl_array_fetch(up->out_coords_nodal[peak_idx], preserved_node_idx);
    val_n[0] = vals[0];
    coord_n[0] = coords[0];
    peak_idx++;
  }

  // Find local maxima and minima.
  for (int i = 1; i < total_nodes_search - 1 && peak_idx < up->num_peaks - 1; i++) {
    double prev = vals[i-1];
    double curr = vals[i];
    double next = vals[i+1];

    bool is_max = (curr > prev && curr > next);
    bool is_min = (curr < prev && curr < next);

    if ((is_max && up->peak_types[peak_idx] == GKYL_PEAK_LOCAL_MAX) ||
        (is_min && up->peak_types[peak_idx] == GKYL_PEAK_LOCAL_MIN)) {
      double *val_n = gkyl_array_fetch(up->out_vals_nodal[peak_idx], preserved_node_idx);
      double *coord_n = gkyl_array_fetch(up->out_coords_nodal[peak_idx], preserved_node_idx);
      val_n[0] = curr;
      coord_n[0] = coords[i];
      peak_idx++;
    }
  }

  // EDGE_HI is always last peak.
  if (peak_idx < up->num_peaks && up->peak_types[peak_idx] == GKYL_PEAK_EDGE_HI) {
    double *val_n = gkyl_array_fetch(up->out_vals_nodal[peak_idx], preserved_node_idx);
    double *coord_n = gkyl_array_fetch(up->out_coords_nodal[peak_idx], preserved_node_idx);
    val_n[0] = vals[total_nodes_search - 1];
    coord_n[0] = coords[total_nodes_search - 1];
  }

  gkyl_free(vals);
  gkyl_free(coords);
  gkyl_free(visited);
}

struct gkyl_array_dg_find_peaks*
gkyl_array_dg_find_peaks_new(const struct gkyl_array_dg_find_peaks_inp *inp, const struct gkyl_array *field)
{
  struct gkyl_array_dg_find_peaks *up = gkyl_malloc(sizeof(*up));

  // Copy input parameters.
  up->grid = *inp->grid;
  up->basis = *inp->basis;
  up->range = *inp->range;
  up->range_ext = *inp->range_ext;
  up->search_dir = inp->search_dir;
  up->use_gpu = inp->use_gpu;

  int ndim = inp->grid->ndim;
  int poly_order = inp->basis->poly_order;
  int out_dim = ndim - 1;

  assert(inp->search_dir >= 0 && inp->search_dir < ndim);

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
    int preserved_dir = (inp->search_dir == 0) ? 1 : 0;

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
    assert(false && "dg_find_peaks: only 1D->0D and 2D->1D supported");
  }

  // Store node locations for input basis.
  up->nodes = gkyl_array_new(GKYL_DOUBLE, ndim, inp->basis->num_basis);
  inp->basis->node_list(gkyl_array_fetch(up->nodes, 0));

  // Create nodal-to-modal converter.
  up->n2m = gkyl_nodal_ops_new(&up->out_basis, &up->out_grid, false);

  // Count peaks at middle preserved coordinate.
  int mid_preserved_idx = 0;
  if (out_dim == 1) {
    int preserved_dir = (inp->search_dir == 0) ? 1 : 0;
    mid_preserved_idx = (inp->range->lower[preserved_dir] + inp->range->upper[preserved_dir]) / 2;
  }

  count_peaks_along_dir(up, field, mid_preserved_idx, &up->num_peaks, up->peak_types);

  // Allocate output arrays for each peak.
  for (int p = 0; p < up->num_peaks; p++) {
    up->out_vals[p] = gkyl_array_new(GKYL_DOUBLE, up->out_basis.num_basis, up->out_range_ext.volume);
    up->out_coords[p] = gkyl_array_new(GKYL_DOUBLE, up->out_basis.num_basis, up->out_range_ext.volume);
    up->out_vals_nodal[p] = gkyl_array_new(GKYL_DOUBLE, 1, up->out_nrange.volume);
    up->out_coords_nodal[p] = gkyl_array_new(GKYL_DOUBLE, 1, up->out_nrange.volume);
  }

  // Initialize unused peak arrays to NULL.
  for (int p = up->num_peaks; p < GKYL_DG_FIND_PEAKS_MAX; p++) {
    up->out_vals[p] = NULL;
    up->out_coords[p] = NULL;
    up->out_vals_nodal[p] = NULL;
    up->out_coords_nodal[p] = NULL;
  }

  return up;
}

void
gkyl_array_dg_find_peaks_advance(struct gkyl_array_dg_find_peaks *up, const struct gkyl_array *in)
{
  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;

  // Copy input to host if on GPU.
  struct gkyl_array *in_ho = gkyl_array_new(GKYL_DOUBLE, in->ncomp, in->size);
  gkyl_array_copy(in_ho, in);

  // Find peaks for each preserved-direction node.
  int num_nodes_out = up->out_nrange.volume;
  for (int pres_node = 0; pres_node < num_nodes_out; pres_node++) {
    find_peaks_for_preserved_node(up, in_ho, pres_node);
  }

  // Transform nodal to modal for each peak.
  if (out_dim == 0) {
    // 1D -> 0D case: modal = nodal (p=0 has no nodal_to_modal function).
    for (int p = 0; p < up->num_peaks; p++) {
      double *val_m = gkyl_array_fetch(up->out_vals[p], 0);
      double *coord_m = gkyl_array_fetch(up->out_coords[p], 0);
      const double *val_n = gkyl_array_cfetch(up->out_vals_nodal[p], 0);
      const double *coord_n = gkyl_array_cfetch(up->out_coords_nodal[p], 0);
      val_m[0] = val_n[0];
      coord_m[0] = coord_n[0];
    }
  }
  else {
    // 2D -> 1D case: use nodal-to-modal transform.
    for (int p = 0; p < up->num_peaks; p++) {
      gkyl_nodal_ops_n2m(up->n2m, &up->out_basis, &up->out_grid,
        &up->out_nrange, &up->out_range, 1, up->out_vals_nodal[p], up->out_vals[p], false);
      gkyl_nodal_ops_n2m(up->n2m, &up->out_basis, &up->out_grid,
        &up->out_nrange, &up->out_range, 1, up->out_coords_nodal[p], up->out_coords[p], false);
    }
  }

  gkyl_array_release(in_ho);
}

int
gkyl_array_dg_find_peaks_num_peaks(const struct gkyl_array_dg_find_peaks *up)
{
  return up->num_peaks;
}

enum gkyl_peak_type
gkyl_array_dg_find_peaks_get_type(const struct gkyl_array_dg_find_peaks *up, int peak_idx)
{
  assert(peak_idx >= 0 && peak_idx < up->num_peaks);
  return up->peak_types[peak_idx];
}

const struct gkyl_basis*
gkyl_array_dg_find_peaks_get_basis(const struct gkyl_array_dg_find_peaks *up)
{
  return &up->out_basis;
}

const struct gkyl_rect_grid*
gkyl_array_dg_find_peaks_get_grid(const struct gkyl_array_dg_find_peaks *up)
{
  return &up->out_grid;
}

const struct gkyl_range*
gkyl_array_dg_find_peaks_get_range(const struct gkyl_array_dg_find_peaks *up)
{
  return &up->out_range;
}

const struct gkyl_range*
gkyl_array_dg_find_peaks_get_range_ext(const struct gkyl_array_dg_find_peaks *up)
{
  return &up->out_range_ext;
}

const struct gkyl_array*
gkyl_array_dg_find_peaks_get_vals(const struct gkyl_array_dg_find_peaks *up, int peak_idx)
{
  assert(peak_idx >= 0 && peak_idx < up->num_peaks);
  return up->out_vals[peak_idx];
}

const struct gkyl_array*
gkyl_array_dg_find_peaks_get_coords(const struct gkyl_array_dg_find_peaks *up, int peak_idx)
{
  assert(peak_idx >= 0 && peak_idx < up->num_peaks);
  return up->out_coords[peak_idx];
}

void
gkyl_array_dg_find_peaks_release(struct gkyl_array_dg_find_peaks *up)
{
  for (int p = 0; p < up->num_peaks; p++) {
    gkyl_array_release(up->out_vals[p]);
    gkyl_array_release(up->out_coords[p]);
    gkyl_array_release(up->out_vals_nodal[p]);
    gkyl_array_release(up->out_coords_nodal[p]);
  }
  gkyl_array_release(up->nodes);
  gkyl_nodal_ops_release(up->n2m);
  gkyl_free(up);
}
