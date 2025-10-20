#include <gkyl_array.h>
#include <gkyl_util.h>
#include <acutest.h>
#include <gkyl_array_rio.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include "math.h"

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_m0s_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx));
}
void eval_m0r_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx));
}
void eval_vtsqs_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  fout[0] = -x*x+1.5*pow(M_PI,2);
}
void eval_vtsqr_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  fout[0] = -x*x+3.5*pow(M_PI,2);
}
void eval_nu_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.*M_PI;
  double norm_nu = 1./3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_1x(t, xn, vtsqs, ctx);
  eval_m0r_1x(t, xn, m0r, ctx);
  eval_vtsqr_1x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0]/pow(vtsqs[0]+vtsqr[0],1.5);
}

void
test_1x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI;
  double lower[] = {-Lx/2.}, upper[] = {Lx/2.};
  int cells[] = {32};

  double norm_nu = 1./3.;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(basis.num_basis, local_ext.volume);
  vtsqs = mkarr(basis.num_basis, local_ext.volume);
  m0r = mkarr(basis.num_basis, local_ext.volume);
  vtsqr = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_cu, *vtsqs_cu, *m0r_cu, *vtsqr_cu;
  if (use_gpu) { // Create device copies
    m0s_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqs_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    m0r_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqr_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_m0s = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0s_1x, NULL);
  gkyl_proj_on_basis *proj_vtsqs = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqs_1x, NULL);
  gkyl_proj_on_basis *proj_m0r = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0r_1x, NULL);
  gkyl_proj_on_basis *proj_vtsqr = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqr_1x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(m0s_cu , m0s );
    gkyl_array_copy(vtsqs_cu, vtsqs);
    gkyl_array_copy(m0r_cu , m0r );
    gkyl_array_copy(vtsqr_cu, vtsqr);
  }

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_cu;
  nu = mkarr(basis.num_basis, local_ext.volume);
  if (use_gpu)  // create device copy.
    nu_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up = gkyl_spitzer_coll_freq_new(&basis, poly_order+1, nufrac, epsilon_0, hbar, use_gpu);

  if (use_gpu) {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs_cu, 0., m0r_cu, vtsqr_cu, 0., norm_nu, nu_cu);
    gkyl_array_copy(nu, nu_cu);
  } else {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs, 0., m0r, vtsqr, 0., norm_nu, nu);
  }

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA;
  nuA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_nu_1x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA);

  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&local, idx);
    const double *nu_p = gkyl_array_cfetch(nu, linidx);
    const double *nuA_p = gkyl_array_cfetch(nuA, linidx);
    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(nuA_p[m], nu_p[m], 1e-10) );
      TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
      TEST_MSG("Produced: %.13e", nu_p[m]);
    }
  }

  gkyl_array_release(m0s); gkyl_array_release(vtsqs); gkyl_array_release(m0r); gkyl_array_release(vtsqr);
  gkyl_array_release(nu); gkyl_array_release(nuA);
  if (use_gpu) {
    gkyl_array_release(m0s_cu); gkyl_array_release(vtsqs_cu); gkyl_array_release(m0r_cu); gkyl_array_release(vtsqr_cu);
    gkyl_array_release(nu_cu);
  }

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void eval_m0s_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly);
}
void eval_m0r_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly);
}
void eval_vtsqs_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  fout[0] = (-x*x+1.5*pow(M_PI,2))*0.5*(1.+cos(0.5*2.*M_PI*y/Ly));
}
void eval_vtsqr_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  fout[0] = (-x*x+3.5*pow(M_PI,2))*0.5*(1.+cos(0.5*2.*M_PI*y/Ly));
}
void eval_nu_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  double norm_nu = 1./3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_2x(t, xn, vtsqs, ctx);
  eval_m0r_2x(t, xn, m0r, ctx);
  eval_vtsqr_2x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0]/pow(vtsqs[0]+vtsqr[0],1.5);
}

void
test_2x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI, Ly = 8.*M_PI;
  double lower[] = {-Lx/2., -Ly/2.}, upper[] = {Lx/2., Ly/2.};
  int cells[] = {32, 128};

  double norm_nu = 1./3.;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(basis.num_basis, local_ext.volume);
  vtsqs = mkarr(basis.num_basis, local_ext.volume);
  m0r = mkarr(basis.num_basis, local_ext.volume);
  vtsqr = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_cu, *vtsqs_cu, *m0r_cu, *vtsqr_cu;
  if (use_gpu) { // Create device copies
    m0s_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqs_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    m0r_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqr_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_m0s = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0s_2x, NULL);
  gkyl_proj_on_basis *proj_vtsqs = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqs_2x, NULL);
  gkyl_proj_on_basis *proj_m0r = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0r_2x, NULL);
  gkyl_proj_on_basis *proj_vtsqr = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqr_2x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(m0s_cu , m0s );
    gkyl_array_copy(vtsqs_cu, vtsqs);
    gkyl_array_copy(m0r_cu , m0r );
    gkyl_array_copy(vtsqr_cu, vtsqr);
  }

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_cu;
  nu = mkarr(basis.num_basis, local_ext.volume);
  if (use_gpu)  // create device copy.
    nu_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up = gkyl_spitzer_coll_freq_new(&basis, poly_order+1, nufrac, epsilon_0, hbar, use_gpu);

  if (use_gpu) {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs_cu, 0., m0r_cu, vtsqr_cu, 0., norm_nu, nu_cu);
    gkyl_array_copy(nu, nu_cu);
  } else {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs, 0., m0r, vtsqr, 0., norm_nu, nu);
  }

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA;
  nuA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_nu_2x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA);

  for (int j=0; j<cells[0]; j++) {
    for (int k=0; k<cells[1]; k++) {
      int idx[] = {j+1,k+1};
      long linidx = gkyl_range_idx(&local, idx);
      const double *nu_p = gkyl_array_cfetch(nu, linidx);
      const double *nuA_p = gkyl_array_cfetch(nuA, linidx);
      for (int m=0; m<basis.num_basis; m++) {
        TEST_CHECK( gkyl_compare(nuA_p[m], nu_p[m], 1e-6) );
        TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
        TEST_MSG("Produced: %.13e", nu_p[m]);
      }
    }
  }

  gkyl_array_release(m0s); gkyl_array_release(vtsqs); gkyl_array_release(m0r); gkyl_array_release(vtsqr);
  gkyl_array_release(nu); gkyl_array_release(nuA);
  if (use_gpu) {
    gkyl_array_release(m0s_cu); gkyl_array_release(vtsqs_cu); gkyl_array_release(m0r_cu); gkyl_array_release(vtsqr_cu);
    gkyl_array_release(nu_cu);
  }

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void eval_m0s_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly)*(Lz/2.-0.5*fabs(z));
}
void eval_m0r_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  fout[0] = 0.5*(1.+cos(0.5*2.*M_PI*x/Lx))*(y+Ly)*(Lz/2.-0.5*fabs(z));;
}
void eval_vtsqs_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  fout[0] = (-x*x+1.5*pow(M_PI,2))*0.5*(1.+cos(0.5*2.*M_PI*y/Ly))*(1.5+tanh(z));
}
void eval_vtsqr_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  fout[0] = (-x*x+3.5*pow(M_PI,2))*0.5*(1.+cos(0.5*2.*M_PI*y/Ly))*(1.5+tanh(z));;
}
void eval_nu_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  double norm_nu = 1./3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_3x(t, xn, vtsqs, ctx);
  eval_m0r_3x(t, xn, m0r, ctx);
  eval_vtsqr_3x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0]/pow(vtsqs[0]+vtsqr[0],1.5);
}

void
test_3x(int poly_order, bool use_gpu)
{
  double Lx = 2.*M_PI, Ly = 8.*M_PI, Lz = 4.*M_PI;
  double lower[] = {-Lx/2., -Ly/2., -Lz/2.}, upper[] = {Lx/2., Ly/2., Lz/2.};
  int cells[] = {32, 128, 64};

  double norm_nu = 1./3.;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(basis.num_basis, local_ext.volume);
  vtsqs = mkarr(basis.num_basis, local_ext.volume);
  m0r = mkarr(basis.num_basis, local_ext.volume);
  vtsqr = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_cu, *vtsqs_cu, *m0r_cu, *vtsqr_cu;
  if (use_gpu) { // Create device copies
    m0s_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqs_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    m0r_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    vtsqr_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  }

  gkyl_proj_on_basis *proj_m0s = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0s_3x, NULL);
  gkyl_proj_on_basis *proj_vtsqs = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqs_3x, NULL);
  gkyl_proj_on_basis *proj_m0r = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_m0r_3x, NULL);
  gkyl_proj_on_basis *proj_vtsqr = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtsqr_3x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr);

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(m0s_cu , m0s );
    gkyl_array_copy(vtsqs_cu, vtsqs);
    gkyl_array_copy(m0r_cu , m0r );
    gkyl_array_copy(vtsqr_cu, vtsqr);
  }

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_cu;
  nu = mkarr(basis.num_basis, local_ext.volume);
  if (use_gpu)  // create device copy.
    nu_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up = gkyl_spitzer_coll_freq_new(&basis, poly_order+1, nufrac, epsilon_0, hbar, use_gpu);

  if (use_gpu) {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs_cu, 0., m0r_cu, vtsqr_cu, 0., norm_nu, nu_cu);
    gkyl_array_copy(nu, nu_cu);
  } else {
    gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, vtsqs, 0., m0r, vtsqr, 0., norm_nu, nu);
  }

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA;
  nuA = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_nu_3x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA);

  for (int i=0; i<cells[0]; i++) {
    for (int j=0; j<cells[1]; j++) {
      for (int k=0; k<cells[2]; k++) {
        int idx[] = {i+1,j+1,k+1};
        long linidx = gkyl_range_idx(&local, idx);
        const double *nu_p = gkyl_array_cfetch(nu, linidx);
        const double *nuA_p = gkyl_array_cfetch(nuA, linidx);
        for (int m=0; m<basis.num_basis; m++) {
          TEST_CHECK( gkyl_compare(nuA_p[m], nu_p[m], 1e-4) );
          TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
          TEST_MSG("Produced: %.13e", nu_p[m]);
        }
      }
    }
  }

  gkyl_array_release(m0s); gkyl_array_release(vtsqs); gkyl_array_release(m0r); gkyl_array_release(vtsqr);
  gkyl_array_release(nu); gkyl_array_release(nuA);
  if (use_gpu) {
    gkyl_array_release(m0s_cu); gkyl_array_release(vtsqs_cu); gkyl_array_release(m0r_cu); gkyl_array_release(vtsqr_cu);
    gkyl_array_release(nu_cu);
  }

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void test_1x_p1() { test_1x(1, false); }
void test_1x_p2() { test_1x(2, false); }

void test_2x_p1() { test_2x(1, false); }
void test_2x_p2() { test_2x(2, false); }

void test_3x_p1() { test_3x(1, false); }
void test_3x_p2() { test_3x(2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x_p1_gpu() { test_1x(1, true); }
void test_1x_p2_gpu() { test_1x(2, true); }

void test_2x_p1_gpu() { test_2x(1, true); }
void test_2x_p2_gpu() { test_2x(2, true); }

void test_3x_p1_gpu() { test_3x(1, true); }
void test_3x_p2_gpu() { test_3x(2, true); }
#endif

TEST_LIST = {
  { "test_1x_p1", test_1x_p1 },
  { "test_1x_p2", test_1x_p2 },

  { "test_2x_p1", test_2x_p1 },
  { "test_2x_p2", test_2x_p2 },

  { "test_3x_p1", test_3x_p1 },
  { "test_3x_p2", test_3x_p2 },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_p1_gpu", test_1x_p1_gpu },
  { "test_1x_p2_gpu", test_1x_p2_gpu },

  { "test_2x_p1_gpu", test_2x_p1_gpu },
  { "test_2x_p2_gpu", test_2x_p2_gpu },

  { "test_3x_p1_gpu", test_3x_p1_gpu },
  { "test_3x_p2_gpu", test_3x_p2_gpu },
#endif
  { NULL, NULL },
};
