/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_dg_find_peaks.h>
#include <gkyl_array_dg_find_peaks_priv.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_nodal_ops.h>
#include <assert.h>
}

// Maximum number of nodes along the search direction.
// For p=2: total_nodes = 2*num_cells + 1. With up to 512 cells this is 1025.
#define MAX_SEARCH_NODES 1025

/**
 * CUDA kernel: find peaks along the search direction for each preserved-direction
 * node index. One thread per preserved_node_idx.
 *
 * Each thread:
 *  1. Scans all cells along the search direction, collecting nodal values/coords
 *     into thread-local arrays.
 *  2. Extracts peaks (EDGE_LO, LOCAL_MAX, LOCAL_MIN, EDGE_HI) and writes results
 *     into the nodal output arrays.
 */
__global__ void
gkyl_find_peaks_kernel(const struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in, int num_nodes_out)
{
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < num_nodes_out; tid += blockDim.x * gridDim.x)
  {
    int preserved_node_idx = (int)tid;

    int ndim = up->grid.ndim;
    int search_dir = up->search_dir;
    int poly_order = up->basis.poly_order;
    int num_basis = up->basis.num_basis;

    // Number of cells and nodes along the search direction.
    int num_cells_search = up->range.upper[search_dir] - up->range.lower[search_dir] + 1;
    int total_nodes_search = (poly_order == 1) ? num_cells_search + 1 : 2*num_cells_search + 1;

    // Thread-local storage for values, coordinates, and visited flags.
    double vals[MAX_SEARCH_NODES];
    double coords[MAX_SEARCH_NODES];
    bool visited[MAX_SEARCH_NODES];
    for (int i = 0; i < total_nodes_search; i++) {
      vals[i] = 0.0;
      coords[i] = 0.0;
      visited[i] = false;
    }

    // Preserved direction (only used for 2D).
    int preserved_dir = (ndim == 1) ? -1 : ((search_dir == 0) ? 1 : 0);

    // Iterate along cells in the search direction.
    for (int cell_idx = up->range.lower[search_dir];
         cell_idx <= up->range.upper[search_dir]; cell_idx++)
    {
      // For 2D, determine which cells in the preserved direction
      // contribute to this preserved_node_idx.
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
            if (pres_cell_end > up->range.upper[preserved_dir])
              pres_cell_end = up->range.upper[preserved_dir];
          }
        }
        else { // poly_order == 2
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
        // Build cell index.
        int idx[GKYL_MAX_DIM];
        if (ndim == 1) {
          idx[0] = cell_idx;
        }
        else {
          idx[preserved_dir] = pres_cell;
          idx[search_dir] = cell_idx;
        }

        long linidx = gkyl_range_idx(&up->range, idx);
        const double *f_d = (const double *)gkyl_array_cfetch(in, linidx);

        double xc[GKYL_MAX_DIM];
        gkyl_rect_grid_cell_center(&up->grid, idx, xc);

        // Evaluate at each node in this cell.
        for (int n = 0; n < num_basis; n++) {
          const double *nod_log = (const double *)gkyl_array_cfetch(up->nodes, n);

          // Check if this node belongs to our preserved_node_idx (2D only).
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

          // Determine node offset in the search direction.
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

    // Extract peaks based on peak_types and write to nodal output arrays.
    int peak_idx = 0;

    // EDGE_LO is always the first peak at index 0.
    if (up->peak_types[peak_idx] == GKYL_PEAK_EDGE_LO) {
      double *val_n = (double *)gkyl_array_fetch(up->out_vals_nodal[peak_idx],
        preserved_node_idx);
      double *coord_n = (double *)gkyl_array_fetch(up->out_coords_nodal[peak_idx],
        preserved_node_idx);
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
        double *val_n = (double *)gkyl_array_fetch(up->out_vals_nodal[peak_idx],
          preserved_node_idx);
        double *coord_n = (double *)gkyl_array_fetch(up->out_coords_nodal[peak_idx],
          preserved_node_idx);
        val_n[0] = curr;
        coord_n[0] = coords[i];
        peak_idx++;
      }
    }

    // EDGE_HI is always the last peak.
    if (peak_idx < up->num_peaks && up->peak_types[peak_idx] == GKYL_PEAK_EDGE_HI) {
      double *val_n = (double *)gkyl_array_fetch(up->out_vals_nodal[peak_idx],
        preserved_node_idx);
      double *coord_n = (double *)gkyl_array_fetch(up->out_coords_nodal[peak_idx],
        preserved_node_idx);
      val_n[0] = vals[total_nodes_search - 1];
      coord_n[0] = coords[total_nodes_search - 1];
    }
  }
}

/**
 * CUDA kernel: evaluate an input array at peak locations for given peak indices.
 * Writes results into out_eval_at_peaks_vals_nodal arrays on device.
 *
 * Thread mapping: one thread per (preserved_node_idx, peak_offset) pair.
 * total_threads = num_nodes_out * num_peaks_to_eval.
 *
 * @param up Device-side updater struct
 * @param in Device-side input array (DG field to evaluate)
 * @param num_nodes_out Number of preserved-direction nodes
 * @param peak_start First peak index to evaluate
 * @param num_peaks_to_eval Number of peaks to evaluate (starting from peak_start)
 */
__global__ void
gkyl_eval_at_peaks_kernel(const struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in, int num_nodes_out,
  int peak_start, int num_peaks_to_eval)
{
  unsigned long total_threads = (unsigned long)num_nodes_out * num_peaks_to_eval;

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x;
    tid < total_threads; tid += blockDim.x * gridDim.x)
  {
    int preserved_node_idx = (int)(tid / num_peaks_to_eval);
    int peak_offset = (int)(tid % num_peaks_to_eval);
    int peak_idx = peak_start + peak_offset;

    int ndim = up->grid.ndim;
    int search_dir = up->search_dir;
    int poly_order = up->basis.poly_order;
    int preserved_dir = (ndim == 1) ? -1 : ((search_dir == 0) ? 1 : 0);

    // Get the peak coordinate found during advance.
    const double *peak_coord_n = (const double *)gkyl_array_cfetch(
      up->out_coords_nodal[peak_idx], preserved_node_idx);
    double peak_coord_search = peak_coord_n[0];

    // Determine cell index containing the peak.
    // We compute the search-direction cell directly from the uniform grid
    // geometry (avoids calling gkyl_rect_grid_find_cell which is not available
    // as a device symbol).
    int cell_idx[GKYL_MAX_DIM];

    // Search direction: compute cell from coordinate on uniform grid.
    // cell = floor((x - lower) / dx) + 1  (1-based indexing).
    // Use pick_lower semantics: if exactly on a boundary, pick the lower cell.
    {
      double rel = (peak_coord_search - up->grid.lower[search_dir]) / up->grid.dx[search_dir];
      int c = (int)rel + 1; // 1-based
      // pick_lower: if exactly on upper boundary of cell c, rel is integer, pick c not c+1.
      // The (int) cast truncates toward zero which gives pick_lower behavior for positive rel.
      // Clamp to valid range.
      if (c < up->range.lower[search_dir])
        c = up->range.lower[search_dir];
      if (c > up->range.upper[search_dir])
        c = up->range.upper[search_dir];
      cell_idx[search_dir] = c;
    }

    // For 2D: determine preserved-direction cell from preserved_node_idx.
    if (ndim > 1) {
      int pres_cell;
      if (poly_order == 1) {
        if (preserved_node_idx == 0) {
          pres_cell = up->range.lower[preserved_dir];
        }
        else {
          pres_cell = up->range.lower[preserved_dir] + preserved_node_idx - 1;
          if (pres_cell > up->range.upper[preserved_dir])
            pres_cell = up->range.upper[preserved_dir];
        }
      }
      else {
        pres_cell = up->range.lower[preserved_dir] + preserved_node_idx / 2;
      }
      cell_idx[preserved_dir] = pres_cell;
    }

    // Fetch DG coefficients at this cell.
    long linidx = gkyl_range_idx(&up->range_ext, cell_idx);
    const double *f_d = (const double *)gkyl_array_cfetch(in, linidx);

    // Get cell center for logical coordinate conversion.
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->grid, cell_idx, xc);

    // Convert peak coordinate to logical space [-1, 1].
    double nod_log[GKYL_MAX_DIM];
    for (int d = 0; d < ndim; d++) {
      if (d == search_dir) {
        nod_log[d] = 2.0 * (peak_coord_search - xc[d]) / up->grid.dx[d];
      }
      else if (ndim > 1) {
        if (poly_order == 1) {
          nod_log[d] = (preserved_node_idx == 0) ? -1.0 : 1.0;
        }
        else {
          int pres_node_offset = preserved_node_idx % 3;
          if (pres_node_offset == 0)
            nod_log[d] = -1.0;
          else if (pres_node_offset == 1)
            nod_log[d] = 0.0;
          else
            nod_log[d] = 1.0;
        }
      }
    }

    // Evaluate the DG expansion and store result.
    double val = up->basis.eval_expand(nod_log, f_d);
    double *val_n = (double *)gkyl_array_fetch(
      up->out_eval_at_peaks_vals_nodal[peak_idx], preserved_node_idx);
    val_n[0] = val;
  }
}

// Host function to launch the project_on_peaks kernel and run nodal-to-modal transforms.
void
gkyl_array_dg_find_peaks_project_on_peaks_cu(struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in_array, struct gkyl_array **out_vals)
{
  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;
  int num_nodes_out = up->out_nrange.volume;
  int num_peaks = up->num_peaks;

  // Launch kernel: one thread per (preserved_node, peak) pair.
  long total_threads = (long)num_nodes_out * num_peaks;
  int nthreads = 256;
  int nblocks = (total_threads + nthreads - 1) / nthreads;

  gkyl_eval_at_peaks_kernel<<<nblocks, nthreads>>>(
    up->on_dev, in_array->on_dev, num_nodes_out, 0, num_peaks);

  // Transform nodal to modal for each peak.
  if (out_dim == 0) {
    for (int p = 0; p < num_peaks; p++) {
      gkyl_array_copy(out_vals[p], up->out_eval_at_peaks_vals_nodal[p]);
    }
  }
  else {
    for (int p = 0; p < num_peaks; p++) {
      gkyl_nodal_ops_n2m_cu(up->n2m, up->out_basis_on_dev, &up->out_grid,
        &up->out_nrange, &up->out_range, 1,
        up->out_eval_at_peaks_vals_nodal[p], out_vals[p]);
    }
  }
}

// Host function to launch the project_on_peak_idx kernel and run nodal-to-modal transform.
void
gkyl_array_dg_find_peaks_project_on_peak_idx_cu(struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in_array, int peak_idx, struct gkyl_array *out_val)
{
  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;
  int num_nodes_out = up->out_nrange.volume;

  // Launch kernel: one thread per preserved_node, single peak.
  int nthreads = 256;
  int nblocks = (num_nodes_out + nthreads - 1) / nthreads;

  gkyl_eval_at_peaks_kernel<<<nblocks, nthreads>>>(
    up->on_dev, in_array->on_dev, num_nodes_out, peak_idx, 1);

  // Transform nodal to modal.
  if (out_dim == 0) {
    gkyl_array_copy(out_val, up->out_eval_at_peaks_vals_nodal[peak_idx]);
  }
  else {
    gkyl_nodal_ops_n2m_cu(up->n2m, up->out_basis_on_dev, &up->out_grid,
      &up->out_nrange, &up->out_range, 1,
      up->out_eval_at_peaks_vals_nodal[peak_idx], out_val);
  }
}

// Host function to launch the advance kernel and run nodal-to-modal transforms.
void
gkyl_array_dg_find_peaks_advance_cu(struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in)
{
  int ndim = up->grid.ndim;
  int out_dim = ndim - 1;
  int num_nodes_out = up->out_nrange.volume;

  // Launch the kernel: one thread per preserved node.
  int nthreads = 256;
  int nblocks = (num_nodes_out + nthreads - 1) / nthreads;

  gkyl_find_peaks_kernel<<<nblocks, nthreads>>>(
    up->on_dev, in->on_dev, num_nodes_out);

  // Transform nodal to modal for each peak.
  if (out_dim == 0) {
    // 1D -> 0D case: modal = nodal (p=0, single value).
    // Copy from nodal to modal arrays on device.
    for (int p = 0; p < up->num_peaks; p++) {
      gkyl_array_copy(up->out_vals[p], up->out_vals_nodal[p]);
      gkyl_array_copy(up->out_coords[p], up->out_coords_nodal[p]);
    }
  }
  else {
    // 2D -> 1D case: use nodal-to-modal transform on GPU.
    for (int p = 0; p < up->num_peaks; p++) {
      gkyl_nodal_ops_n2m_cu(up->n2m, up->out_basis_on_dev, &up->out_grid,
        &up->out_nrange, &up->out_range, 1,
        up->out_vals_nodal[p], up->out_vals[p]);
      gkyl_nodal_ops_n2m_cu(up->n2m, up->out_basis_on_dev, &up->out_grid,
        &up->out_nrange, &up->out_range, 1,
        up->out_coords_nodal[p], up->out_coords[p]);
    }
  }
}

struct gkyl_array_dg_find_peaks*
gkyl_array_dg_find_peaks_new_cu(struct gkyl_array_dg_find_peaks *up_ho)
{
  struct gkyl_array_dg_find_peaks *up =
    (struct gkyl_array_dg_find_peaks *)gkyl_malloc(sizeof(*up));

  // Copy all scalar/struct fields from host object.
  up->grid = up_ho->grid;
  up->basis = up_ho->basis;
  up->range = up_ho->range;
  up->range_ext = up_ho->range_ext;
  up->search_dir = up_ho->search_dir;
  up->use_gpu = true;

  up->out_grid = up_ho->out_grid;
  up->out_basis = up_ho->out_basis;
  up->out_range = up_ho->out_range;
  up->out_range_ext = up_ho->out_range_ext;
  up->out_nrange = up_ho->out_nrange;

  up->num_peaks = up_ho->num_peaks;
  for (int p = 0; p < GKYL_DG_FIND_PEAKS_MAX; p++)
    up->peak_types[p] = up_ho->peak_types[p];

  int ndim = up_ho->basis.ndim;
  int poly_order = up_ho->basis.poly_order;
  int out_dim = ndim - 1;

  // Create a GPU copy of the nodes array so the kernel can access it.
  up->nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    up_ho->nodes->ncomp, up_ho->nodes->size);
  gkyl_array_copy(up->nodes, up_ho->nodes);

  // Create GPU-enabled nodal-to-modal converter.
  // Use the host basis (up->out_basis) here because gkyl_nodal_ops_new
  // calls cbasis->node_list on the host.
  up->n2m = gkyl_nodal_ops_new(&up->out_basis, &up->out_grid, true);

  // Create a device-resident basis with device-callable function pointers.
  // This is needed by gkyl_nodal_ops_n2m_cu which passes the basis pointer
  // directly to a CUDA kernel that dereferences cbasis->nodal_to_modal().
  if (out_dim > 0)
    up->out_basis_on_dev = gkyl_cart_modal_serendip_cu_dev_new(1, poly_order);
  else
    up->out_basis_on_dev = NULL;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_array_dg_find_peaks_free);

  // Allocate GPU output arrays for each peak and copy data from host arrays.
  for (int p = 0; p < up->num_peaks; p++) {
    up->out_vals[p] = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up_ho->out_vals[p]->ncomp, up_ho->out_vals[p]->size);
    gkyl_array_copy(up->out_vals[p], up_ho->out_vals[p]);

    up->out_coords[p] = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up_ho->out_coords[p]->ncomp, up_ho->out_coords[p]->size);
    gkyl_array_copy(up->out_coords[p], up_ho->out_coords[p]);

    up->out_vals_nodal[p] = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up_ho->out_vals_nodal[p]->ncomp, up_ho->out_vals_nodal[p]->size);
    gkyl_array_copy(up->out_vals_nodal[p], up_ho->out_vals_nodal[p]);

    up->out_coords_nodal[p] = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up_ho->out_coords_nodal[p]->ncomp, up_ho->out_coords_nodal[p]->size);
    gkyl_array_copy(up->out_coords_nodal[p], up_ho->out_coords_nodal[p]);

    up->out_eval_at_peaks_vals_nodal[p] = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up_ho->out_eval_at_peaks_vals_nodal[p]->ncomp,
      up_ho->out_eval_at_peaks_vals_nodal[p]->size);
    gkyl_array_copy(up->out_eval_at_peaks_vals_nodal[p],
      up_ho->out_eval_at_peaks_vals_nodal[p]);
  }

  // Initialize unused peak arrays to NULL.
  for (int p = up->num_peaks; p < GKYL_DG_FIND_PEAKS_MAX; p++) {
    up->out_vals[p] = NULL;
    up->out_coords[p] = NULL;
    up->out_vals_nodal[p] = NULL;
    up->out_coords_nodal[p] = NULL;
    up->out_eval_at_peaks_vals_nodal[p] = NULL;
  }

  // Copy struct to device, with on_dev array pointers and device-callable
  // basis function pointers swapped in.
  // Save host-side array pointers and basis structs.
  struct gkyl_array *ho_nodes = up->nodes;
  struct gkyl_basis ho_basis = up->basis;
  struct gkyl_basis ho_out_basis = up->out_basis;
  struct gkyl_array *ho_out_vals[GKYL_DG_FIND_PEAKS_MAX];
  struct gkyl_array *ho_out_coords[GKYL_DG_FIND_PEAKS_MAX];
  struct gkyl_array *ho_out_vals_nodal[GKYL_DG_FIND_PEAKS_MAX];
  struct gkyl_array *ho_out_coords_nodal[GKYL_DG_FIND_PEAKS_MAX];
  struct gkyl_array *ho_out_eval_at_peaks_vals_nodal[GKYL_DG_FIND_PEAKS_MAX];

  // Swap in device-callable basis function pointers for the H2D copy.
  gkyl_cart_modal_serendip_cu_dev(&up->basis, ndim, poly_order);
  if (out_dim == 0)
    gkyl_cart_modal_serendip_cu_dev(&up->out_basis, 1, 0);
  else
    gkyl_cart_modal_serendip_cu_dev(&up->out_basis, 1, poly_order);

  // Swap nodes to its device pointer.
  up->nodes = up->nodes->on_dev;

  for (int p = 0; p < up->num_peaks; p++) {
    ho_out_vals[p] = up->out_vals[p];
    ho_out_coords[p] = up->out_coords[p];
    ho_out_vals_nodal[p] = up->out_vals_nodal[p];
    ho_out_coords_nodal[p] = up->out_coords_nodal[p];
    ho_out_eval_at_peaks_vals_nodal[p] = up->out_eval_at_peaks_vals_nodal[p];

    // Swap in device pointers for the H2D copy.
    up->out_vals[p] = up->out_vals[p]->on_dev;
    up->out_coords[p] = up->out_coords[p]->on_dev;
    up->out_vals_nodal[p] = up->out_vals_nodal[p]->on_dev;
    up->out_coords_nodal[p] = up->out_coords_nodal[p]->on_dev;
    up->out_eval_at_peaks_vals_nodal[p] = up->out_eval_at_peaks_vals_nodal[p]->on_dev;
  }

  // Allocate device struct and copy host struct (with device pointers) to device.
  struct gkyl_array_dg_find_peaks *up_cu =
    (struct gkyl_array_dg_find_peaks *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_dg_find_peaks), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;

  // Restore host-side array pointers and basis so the returned object
  // has usable host handles and host-callable function pointers.
  up->nodes = ho_nodes;
  up->basis = ho_basis;
  up->out_basis = ho_out_basis;
  for (int p = 0; p < up->num_peaks; p++) {
    up->out_vals[p] = ho_out_vals[p];
    up->out_coords[p] = ho_out_coords[p];
    up->out_vals_nodal[p] = ho_out_vals_nodal[p];
    up->out_coords_nodal[p] = ho_out_coords_nodal[p];
    up->out_eval_at_peaks_vals_nodal[p] = ho_out_eval_at_peaks_vals_nodal[p];
  }

  return up;
}
