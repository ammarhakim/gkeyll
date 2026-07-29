#pragma once

// Private header for dg_lowpass_filter updater, not for direct use in user code.

#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_util.h>
#include <assert.h>
#include <math.h>

// Primary struct in this updater.
struct gkyl_dg_lowpass_filter {
  int ndim; // Dimensionality of the field.
  bool use_gpu; // Whether to use the GPU.
  int dir; // Direction along which to filter.
  int half_width; // Stencil half-width M (stencil spans 2M+1 cells).
  int num_basis; // Number of DG coefficients per cell.
  struct gkyl_rect_grid grid; // Grid the field is defined on.
  struct gkyl_range range; // Range to filter in.
  double *weights; // 2M+1 filter weights, normalized to sum to 1.
};

static inline int
dg_lpf_mirror_idx(int idx, int lo, int up)
{
  // Reflect an out-of-range index back in about the outer cell faces.
  while (idx < lo || idx > up) {
    if (idx < lo) idx = 2*lo - 1 - idx;
    if (idx > up) idx = 2*up + 1 - idx;
  }
  return idx;
}

static void
dg_lpf_calc_weights(int half_width, double fc, double *weights)
{
  // Sinc (cutoff fc in cycles/cell) times a Blackman window, normalized.
  int M = half_width;
  double wsum = 0.0;
  for (int k=-M; k<M+1; k++) {
    double hk = k == 0? 2.0*fc : sin(2.0*M_PI*fc*k)/(M_PI*k);
    double wk = 0.42 + 0.5*cos(M_PI*k/M) + 0.08*cos(2.0*M_PI*k/M);
    weights[k+M] = hk*wk;
    wsum += weights[k+M];
  }
  for (int k=0; k<2*M+1; k++)
    weights[k] /= wsum;
}
