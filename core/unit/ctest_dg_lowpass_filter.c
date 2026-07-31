// Test the dg_lowpass_filter updater: a low-pass Blackman-windowed sinc FIR
// filter applied along one direction of a DG field, meant e.g. to de-alias
// the twist-shift BC by removing content beyond the coarse-grid resolution
// before restriction.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lowpass_filter.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a = on_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static double
filter_gain(int M, double fc, double freq)
{
  // Frequency response of the normalized Blackman-windowed sinc kernel
  // at freq cycles/cell (real, since the kernel is symmetric).
  double wsum = 0.0, gain = 0.0;
  for (int k=-M; k<M+1; k++) {
    double hk = k == 0? 2.0*fc : sin(2.0*M_PI*fc*k)/(M_PI*k);
    double wk = 0.42 + 0.5*cos(M_PI*k/M) + 0.08*cos(2.0*M_PI*k/M);
    wsum += hk*wk;
    gain += hk*wk*cos(2.0*M_PI*freq*k);
  }
  return gain/wsum;
}

void eval_const_1x(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 3.0;
}

void eval_linear_1x(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.5 + 0.5*xn[0];
}

struct mode_ctx {
  double mode_num; // Mode number m; the domain is assumed to be [0,1].
};

void eval_mode_1x(double t, const double *xn, double *fout, void *ctx)
{
  struct mode_ctx *mctx = ctx;
  fout[0] = cos(2.0*M_PI*mctx->mode_num*xn[0]);
}

void eval_gauss_1x(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = exp(-pow((xn[0]-0.5)/0.04, 2));
}

void eval_gauss_edge_1x(double t, const double *xn, double *fout, void *ctx)
{
  // Peaks at the lower boundary, so the stencil there is heavily reflected.
  fout[0] = exp(-pow(xn[0]/0.04, 2));
}

void eval_mode_3x(double t, const double *xn, double *fout, void *ctx)
{
  // Separable function: an x-Nyquist mode (for a 32 cell grid on [0,1])
  // times a smooth (y,z) profile.
  struct mode_ctx *mctx = ctx;
  fout[0] = cos(2.0*M_PI*mctx->mode_num*xn[0])
    * (1.0 + 0.3*cos(2.0*M_PI*xn[1])) * (1.0 + 0.2*xn[2]);
}

static void
test_1x(bool use_gpu)
{
  int poly_order = 1;
  int cells[] = {32};
  double lower[] = {0.0}, upper[] = {1.0};
  int nghost[] = {1};
  int M = 8;
  double fc = 0.3;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_dg_lowpass_filter *lpf = gkyl_dg_lowpass_filter_new(0, M,
    grid.dx[0]/fc, &basis, &grid, &local, use_gpu);

  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fout = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size) : gkyl_array_acquire(fin);
  struct gkyl_array *fout_ho = use_gpu? mkarr(false, fout->ncomp, fout->size) : gkyl_array_acquire(fout);

  // a) A constant field is preserved exactly everywhere, including at
  // the boundaries where the stencil is reflected.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_const_1x, NULL);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fin, fin_ho);

  gkyl_array_clear(fout, 7.0);
  gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
  gkyl_array_copy(fout_ho, fout);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *in_c = gkyl_array_cfetch(fin_ho, linidx);
    const double *out_c = gkyl_array_cfetch(fout_ho, linidx);
    for (int c=0; c<basis.num_basis; c++)
      TEST_CHECK( gkyl_compare(in_c[c], out_c[c], 1e-12) );
  }
  // Ghost cells of fout are untouched.
  int idx_ghost[] = {local.lower[0]-1};
  const double *out_g = gkyl_array_cfetch(fout_ho, gkyl_range_idx(&local_ext, idx_ghost));
  for (int c=0; c<basis.num_basis; c++)
    TEST_CHECK( out_g[c] == 7.0 );

  // b) A linear field is preserved in the interior (where the symmetric
  // kernel's first moment vanishes).
  proj = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, eval_linear_1x, NULL);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fin, fin_ho);

  gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
  gkyl_array_copy(fout_ho, fout);

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[0] < local.lower[0]+M || iter.idx[0] > local.upper[0]-M)
      continue;
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *in_c = gkyl_array_cfetch(fin_ho, linidx);
    const double *out_c = gkyl_array_cfetch(fout_ho, linidx);
    for (int c=0; c<basis.num_basis; c++)
      TEST_CHECK( gkyl_compare(in_c[c], out_c[c], 1e-12) );
  }

  // c) A single Fourier mode is scaled by the kernel's frequency
  // response: its DG coefficients are discrete sinusoids in the cell
  // index, so in the interior fout = gain*fin to machine precision.
  // The Nyquist mode (m=16) checks the filter kills grid-scale content.
  double mode_nums[] = {2.0, 8.0, 16.0};
  for (int im=0; im<3; im++) {
    struct mode_ctx mctx = { .mode_num = mode_nums[im] };
    proj = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, eval_mode_1x, &mctx);
    gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
    gkyl_proj_on_basis_release(proj);
    gkyl_array_copy(fin, fin_ho);

    gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
    gkyl_array_copy(fout_ho, fout);

    double freq = mctx.mode_num/cells[0]; // Cycles/cell.
    double gain = filter_gain(M, fc, freq);

    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[0] < local.lower[0]+M || iter.idx[0] > local.upper[0]-M)
        continue;
      long linidx = gkyl_range_idx(&local, iter.idx);
      const double *in_c = gkyl_array_cfetch(fin_ho, linidx);
      const double *out_c = gkyl_array_cfetch(fout_ho, linidx);
      for (int c=0; c<basis.num_basis; c++)
        TEST_CHECK( fabs(out_c[c] - gain*in_c[c]) < 1e-12 );
    }

    if (im == 2) {
      // The Nyquist mode is strongly damped in absolute terms.
      TEST_CHECK( fabs(gain) < 1e-3 );
      gkyl_range_iter_init(&iter, &local);
      while (gkyl_range_iter_next(&iter)) {
        if (iter.idx[0] < local.lower[0]+M || iter.idx[0] > local.upper[0]-M)
          continue;
        const double *out_c = gkyl_array_cfetch(fout_ho, gkyl_range_idx(&local, iter.idx));
        for (int c=0; c<basis.num_basis; c++)
          TEST_CHECK( fabs(out_c[c]) < 1e-3 );
      }
    }
  }

  gkyl_dg_lowpass_filter_release(lpf);
  gkyl_array_release(fin);
  gkyl_array_release(fout);
  gkyl_array_release(fin_ho);
  gkyl_array_release(fout_ho);
}

static void
test_1x_conservation(bool use_gpu)
{
  // The total integral is conserved everywhere: the kernel weights sum to 1 (reflected stencil at the boundaries).
  int poly_order = 1;
  int cells[] = {64};
  double lower[] = {0.0}, upper[] = {1.0};
  int nghost[] = {1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fout = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size) : gkyl_array_acquire(fin);
  struct gkyl_array *fout_ho = use_gpu? mkarr(false, fout->ncomp, fout->size) : gkyl_array_acquire(fout);

  // A ramp peaks at one boundary, a boundary-hugging gaussian at the other.
  evalf_t evals[] = {eval_linear_1x, eval_gauss_edge_1x, eval_gauss_1x};
  int half_widths[] = {8, 8, 80}; // The last stencil is wider than the grid.

  for (int q=0; q<3; q++) {
    struct gkyl_dg_lowpass_filter *lpf = gkyl_dg_lowpass_filter_new(0, half_widths[q],
      grid.dx[0]/0.3, &basis, &grid, &local, use_gpu);

    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evals[q], NULL);
    gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
    gkyl_proj_on_basis_release(proj);
    gkyl_array_copy(fin, fin_ho);

    gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
    gkyl_array_copy(fout_ho, fout);

    double tot_in = 0.0, tot_out = 0.0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long linidx = gkyl_range_idx(&local, iter.idx);
      tot_in += ((const double *) gkyl_array_cfetch(fin_ho, linidx))[0];
      tot_out += ((const double *) gkyl_array_cfetch(fout_ho, linidx))[0];
    }
    TEST_CHECK( fabs(tot_out-tot_in) < 1e-12*fabs(tot_in) );
    TEST_MSG("case %d (M=%d): total in %.13e | out %.13e | rel change %.3e",
      q, half_widths[q], tot_in, tot_out, fabs(tot_out-tot_in)/fabs(tot_in));

    gkyl_dg_lowpass_filter_release(lpf);
  }

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  gkyl_array_release(fin_ho);
  gkyl_array_release(fout_ho);
}

static void
test_3x(bool use_gpu)
{
  // Filter along x of a 3D field. For a separable field g(x)*h(y,z) the
  // p=1 tensor coefficients factor too, so filtering in x scales the
  // coefficients by the kernel's response at g's frequency, leaving the
  // (y,z) dependence untouched.
  int poly_order = 1;
  int cells[] = {32, 8, 6};
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int nghost[] = {1, 1, 1};
  int M = 8;
  double fc = 0.3;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_dg_lowpass_filter *lpf = gkyl_dg_lowpass_filter_new(0, M,
    grid.dx[0]/fc, &basis, &grid, &local, use_gpu);

  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fout = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size) : gkyl_array_acquire(fin);
  struct gkyl_array *fout_ho = use_gpu? mkarr(false, fout->ncomp, fout->size) : gkyl_array_acquire(fout);

  struct mode_ctx mctx = { .mode_num = 16.0 }; // x-Nyquist mode for 32 cells.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_mode_3x, &mctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fin, fin_ho);

  gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
  gkyl_array_copy(fout_ho, fout);

  double freq = mctx.mode_num/cells[0]; // Cycles/cell.
  double gain = filter_gain(M, fc, freq);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[0] < local.lower[0]+M || iter.idx[0] > local.upper[0]-M)
      continue;
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *in_c = gkyl_array_cfetch(fin_ho, linidx);
    const double *out_c = gkyl_array_cfetch(fout_ho, linidx);
    for (int c=0; c<basis.num_basis; c++)
      TEST_CHECK( fabs(out_c[c] - gain*in_c[c]) < 1e-12 );
  }

  gkyl_dg_lowpass_filter_release(lpf);
  gkyl_array_release(fin);
  gkyl_array_release(fout);
  gkyl_array_release(fin_ho);
  gkyl_array_release(fout_ho);
}

static void
test_1x_reflection(bool use_gpu)
{
  // Verify that the reflexion at the boundary is done correctly.
  int poly_order = 1;
  int cells[] = {64};
  double lower[] = {0.0}, upper[] = {1.0};
  int nghost[] = {1};
  int M = 8;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_dg_lowpass_filter *lpf = gkyl_dg_lowpass_filter_new(0, M,
    grid.dx[0]/0.25, &basis, &grid, &local, use_gpu);

  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fout = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size) : gkyl_array_acquire(fin);
  struct gkyl_array *fout_ho = use_gpu? mkarr(false, fout->ncomp, fout->size) : gkyl_array_acquire(fout);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_linear_1x, NULL);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, fin_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fin, fin_ho);

  gkyl_dg_lowpass_filter_advance(lpf, fin, fout);
  gkyl_array_copy(fout_ho, fout);

  // Basis evaluated at the two cell edges, to take the jump across a face.
  double b_lo[8], b_up[8];
  basis.eval((double[]) {-1.0}, b_lo);
  basis.eval((double[]) { 1.0}, b_up);

  for (int i=local.lower[0]; i<local.upper[0]; i++) {
    const double *left = gkyl_array_cfetch(fout_ho, gkyl_range_idx(&local, (int[]) {i}));
    const double *right = gkyl_array_cfetch(fout_ho, gkyl_range_idx(&local, (int[]) {i+1}));
    double jump = 0.0;
    for (int c=0; c<basis.num_basis; c++)
      jump += right[c]*b_lo[c] - left[c]*b_up[c];
    TEST_CHECK( fabs(jump) < 1.0e-12 );
    TEST_MSG("interface %d: filtering a ramp opened a jump of %.3e", i, jump);
  }

  gkyl_dg_lowpass_filter_release(lpf);
  gkyl_array_release(fin);
  gkyl_array_release(fout);
  gkyl_array_release(fin_ho);
  gkyl_array_release(fout_ho);
}

void test_1x_reflection_ho() { test_1x_reflection(false); }
void test_1x_ho() { test_1x(false); }
void test_1x_conservation_ho() { test_1x_conservation(false); }
void test_3x_ho() { test_3x(false); }

#ifdef GKYL_HAVE_CUDA
void test_1x_reflection_cu() { test_1x_reflection(true); }
void test_1x_cu() { test_1x(true); }
void test_1x_conservation_cu() { test_1x_conservation(true); }
void test_3x_cu() { test_3x(true); }
#endif

TEST_LIST = {
  { "test_1x_reflection", test_1x_reflection_ho },
  { "test_1x", test_1x_ho },
  { "test_1x_conservation", test_1x_conservation_ho },
  { "test_3x", test_3x_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_reflection_cu", test_1x_reflection_cu },
  { "test_1x_cu", test_1x_cu },
  { "test_1x_conservation_cu", test_1x_conservation_cu },
  { "test_3x_cu", test_3x_cu },
#endif
  { NULL, NULL },
};
