#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Find all peaks (local maxima, local minima, and boundary values) of a DG
 * field along one direction.
 * 
 * For a 2D input array f(psi, z), finding peaks along z (dir=1) gives arrays:
 *   out_val[k](psi) = value of k-th peak along z for each psi
 *   out_coord[k](psi) = z-coordinate of k-th peak for each psi
 * 
 * For a 1D input array f(z), finding peaks along z (dir=0) gives scalars:
 *   out_val[k] = value of k-th peak
 *   out_coord[k] = z-coordinate of k-th peak
 * 
 * Peaks are detected by sampling the field at nodal points along the search
 * direction and identifying:
 *   - EDGE_LO: Value at the lower boundary of the domain
 *   - LOCAL_MAX: Points where f increases then decreases
 *   - LOCAL_MIN: Points where f decreases then increases  
 *   - EDGE_HI: Value at the upper boundary of the domain
 * 
 * The number of peaks is determined by scanning along the search direction
 * at a middle preserved-direction coordinate.
 */
typedef struct gkyl_array_dg_find_peaks gkyl_array_dg_find_peaks;

/** Types of peaks that can be found. */
enum gkyl_peak_type {
  GKYL_PEAK_EDGE_LO,    // Value at lower boundary
  GKYL_PEAK_LOCAL_MAX,  // Local maximum
  GKYL_PEAK_LOCAL_MIN,  // Local minimum
  GKYL_PEAK_EDGE_HI,    // Value at upper boundary
};

/** Input parameters for dg_find_peaks updater. */
struct gkyl_array_dg_find_peaks_inp {
  const struct gkyl_basis *basis;       // Input basis (N-dimensional)
  const struct gkyl_rect_grid *grid;    // Input grid
  const struct gkyl_range *range;       // Input range (local)
  const struct gkyl_range *range_ext;   // Input extended range
  int search_dir;                       // Direction to search for peaks (0-indexed)
  bool use_gpu;                         // Whether to run on GPU
};

/**
 * Create a new peak finder updater. The number of peaks is determined by
 * scanning the input field along the search direction at a middle coordinate.
 * This must be called AFTER the input field is initialized, as it scans the
 * field to determine the number of peaks.
 * 
 * @param inp Input parameters
 * @param field Input field to scan for peak count determination
 * @return New updater pointer
 */
struct gkyl_array_dg_find_peaks* gkyl_array_dg_find_peaks_new(
  const struct gkyl_array_dg_find_peaks_inp *inp, const struct gkyl_array *field);

/**
 * Compute the peaks. For each point along the preserved dimensions,
 * find all peaks along the search direction.
 * 
 * @param up Updater object
 * @param in Input array (N-dimensional DG field)
 */
void gkyl_array_dg_find_peaks_advance(struct gkyl_array_dg_find_peaks *up, const struct gkyl_array *in);

/**
 * Get the number of peaks found.
 * 
 * @param up Updater object
 * @return Number of peaks
 */
int gkyl_array_dg_find_peaks_num_peaks(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the type of a specific peak (EDGE_LO, LOCAL_MAX, LOCAL_MIN, EDGE_HI).
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Type of the peak
 */
enum gkyl_peak_type gkyl_array_dg_find_peaks_get_type(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Get the output basis ((N-1)-dimensional, or p=0 1D for 1D->0D).
 * 
 * @param up Updater object
 * @return Pointer to output basis
 */
const struct gkyl_basis* gkyl_array_dg_find_peaks_get_basis(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the output grid.
 * 
 * @param up Updater object
 * @return Pointer to output grid
 */
const struct gkyl_rect_grid* gkyl_array_dg_find_peaks_get_grid(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the output range.
 * 
 * @param up Updater object
 * @return Pointer to output range
 */
const struct gkyl_range* gkyl_array_dg_find_peaks_get_range(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the output extended range.
 * 
 * @param up Updater object
 * @return Pointer to output extended range
 */
const struct gkyl_range* gkyl_array_dg_find_peaks_get_range_ext(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the output nodal range.
 * 
 * @param up Updater object
 * @return Pointer to output nodal range
 */
const struct gkyl_range*
gkyl_array_dg_find_peaks_get_nodal_range(const struct gkyl_array_dg_find_peaks *up);

/**
 * Get the output array containing peak values for a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output values array (modal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_vals(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Get the output array containing peak values in nodal basis for a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output values array (nodal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_vals_nodal(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Get the output array containing coordinates of a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output coordinates array (modal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_coords(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Get the output array containing coordinates in nodal basis of a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output coordinates array (nodal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_coords_nodal(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Project (evaluate) an arbitrary array onto the peak locations previously
 * found by gkyl_array_dg_find_peaks_advance.
 * 
 * For a 1D case with 5 peaks, this evaluates the input array at those 5 peak
 * locations and returns the values.
 * 
 * For a 2D case with peaks along lines (e.g., psi vs z with peaks in z),
 * this evaluates the input array along the contours defined by the peak
 * locations for each psi.
 * 
 * The peak locations must have been previously computed via
 * gkyl_array_dg_find_peaks_advance. This method evaluates the provided array
 * at those same locations.
 * 
 * Example usage:
 * @code
 * // 1. Find peaks in bmag along z direction
 * struct gkyl_array_dg_find_peaks *peak_finder = gkyl_array_dg_find_peaks_new(&inp, bmag);
 * gkyl_array_dg_find_peaks_advance(peak_finder, bmag);
 * 
 * // 2. Get bmag_max (LOCAL_MAX peak) location and value
 * int num_peaks = gkyl_array_dg_find_peaks_num_peaks(peak_finder);
 * int bmag_max_idx = -1;
 * for (int p = 0; p < num_peaks; p++) {
 *   if (gkyl_array_dg_find_peaks_get_type(peak_finder, p) == GKYL_PEAK_LOCAL_MAX) {
 *     bmag_max_idx = p;
 *     break;
 *   }
 * }
 * const struct gkyl_array *bmag_max = gkyl_array_dg_find_peaks_get_vals(peak_finder, bmag_max_idx);
 * const struct gkyl_array *z_max = gkyl_array_dg_find_peaks_get_coords(peak_finder, bmag_max_idx);
 * 
 * // 3. Evaluate phi at the same locations where bmag has peaks
 * struct gkyl_array *phi_at_peaks[num_peaks];
 * for (int p = 0; p < num_peaks; p++) {
 *   phi_at_peaks[p] = gkyl_array_new(GKYL_DOUBLE, out_basis.num_basis, out_range_ext.volume);
 * }
 * gkyl_array_dg_find_peaks_project_on_peaks(peak_finder, phi, phi_at_peaks);
 * 
 * // 4. Now phi_at_peaks[bmag_max_idx] contains phi evaluated at the mirror throat
 * @endcode
 * 
 * @param up Updater object (must have run advance first)
 * @param in_array Array to evaluate at peak locations (same grid/basis as original field)
 * @param out_vals Output: array of evaluated values for each peak
 *                 (must be pre-allocated with num_peaks elements, each matching out_range_ext)
 */
void gkyl_array_dg_find_peaks_project_on_peaks(struct gkyl_array_dg_find_peaks *up,
  const struct gkyl_array *in_array, struct gkyl_array **out_vals);

/**
 * Release the updater and all internal arrays.
 * 
 * @param up Updater to delete
 */
void gkyl_array_dg_find_peaks_release(struct gkyl_array_dg_find_peaks *up);

