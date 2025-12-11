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
 * Get the output array containing peak values for a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output values array (modal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_vals(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Get the output array containing coordinates of a specific peak.
 * 
 * @param up Updater object
 * @param peak_idx Index of the peak (0 to num_peaks-1)
 * @return Pointer to output coordinates array (modal DG expansion)
 */
const struct gkyl_array* gkyl_array_dg_find_peaks_get_coords(const struct gkyl_array_dg_find_peaks *up, int peak_idx);

/**
 * Release the updater and all internal arrays.
 * 
 * @param up Updater to delete
 */
void gkyl_array_dg_find_peaks_release(struct gkyl_array_dg_find_peaks *up);

